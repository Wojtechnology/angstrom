"""Rigid holo redocking baseline with AutoDock Vina.

For every subset system the crystal (holo) receptor is kept rigid and the ligand, started from a
fresh random conformer, is redocked into a box around the crystal ligand. This is an *easier*
problem than co-folding (the pocket is already in its bound conformation) and serves as the
"physics-only" reference the co-folding poses are compared against.

Runs under the micromamba `dock` env (vina 1.2, meeko 0.8, rdkit):
  $HOME/.local/micromamba/envs/dock/bin/python dock_vina.py [--systems ID ...] [--workers N] [--exhaustiveness 16]

Outputs raw/docked/vina/<system_id>/pose.sdf (top pose, bond orders from the crystal ligand) and
result.json (all pose scores, box, timing) plus raw/docked/vina/log.jsonl.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import traceback
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")
HERE = Path(__file__).parent
RAW = HERE / "raw"
OUT = RAW / "docked" / "vina"

AMINO = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO",
         "SER", "THR", "TRP", "TYR", "VAL", "MSE", "SEC", "PYL", "HID", "HIE", "HIP", "CYX", "ASH", "GLH", "LYN"}


# --------------------------------------------------------------------------- receptor
def receptor_pdb_from_cif(cif_path: Path) -> str:
    """Protein heavy atoms of the crystal receptor as a PDB block (waters/ions/cofactors dropped)."""
    lines = []
    header = None
    with open(cif_path) as f:
        content = f.read()
    # minimal mmCIF _atom_site loop parser (gemmi is not in the docking env)
    block = content.split("loop_")
    cols, rows = None, []
    for b in block:
        s = b.strip().splitlines()
        if s and s[0].startswith("_atom_site."):
            cols = [l.strip() for l in s if l.startswith("_atom_site.")]
            rows = [l for l in s[len(cols):] if l and not l.startswith("#") and not l.startswith("_")]
            break
    if cols is None:
        raise RuntimeError("no _atom_site loop")
    idx = {c.split(".", 1)[1]: i for i, c in enumerate(cols)}
    import shlex
    serial = 0
    chain_map = {}
    for r in rows:
        t = shlex.split(r)
        if len(t) < len(cols):
            continue
        if t[idx["group_PDB"]] != "ATOM":
            continue
        resn = t[idx["label_comp_id"]]
        if resn not in AMINO:
            continue
        el = t[idx["type_symbol"]].capitalize()
        if el == "H" or el == "D":
            continue
        name = t[idx["label_atom_id"]].strip('"')
        chain = t[idx["auth_asym_id"]] if "auth_asym_id" in idx else t[idx["label_asym_id"]]
        if chain not in chain_map:
            chain_map[chain] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[len(chain_map) % 26]
        ch = chain_map[chain]
        resi = int(t[idx["auth_seq_id"]]) if "auth_seq_id" in idx and t[idx["auth_seq_id"]] not in (".", "?") else int(t[idx["label_seq_id"]])
        alt = t[idx["label_alt_id"]] if "label_alt_id" in idx else "."
        if alt not in (".", "?", "A"):
            continue
        x, y, z = (float(t[idx[k]]) for k in ("Cartn_x", "Cartn_y", "Cartn_z"))
        serial += 1
        nm = name if len(name) == 4 else f" {name:<3s}"
        lines.append(f"ATOM  {serial % 100000:5d} {nm:4s} {resn:>3s} {ch:1s}{resi % 10000:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2s}\n")
    return "".join(lines) + "END\n"


def parse_protein(pdb_text: str):
    """PDB -> RDKit mol; spurious proximity bonds (clashing atoms) removed until it sanitises."""
    from rdkit.Chem import rdMolTransforms
    mol = Chem.MolFromPDBBlock(pdb_text, removeHs=False, sanitize=False, proximityBonding=True)
    rw = Chem.RWMol(mol)
    for _ in range(50):
        probs = [p for p in Chem.DetectChemistryProblems(rw) if p.GetType() == "AtomValenceException"]
        if not probs:
            break
        conf = rw.GetConformer()
        for p in probs:
            a = rw.GetAtomWithIdx(p.GetAtomIdx())
            bonds = sorted(a.GetBonds(), key=lambda b: -rdMolTransforms.GetBondLength(conf, b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
            if bonds:
                rw.RemoveBond(bonds[0].GetBeginAtomIdx(), bonds[0].GetEndAtomIdx())
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    return mol


def ad4_type(atom: Chem.Atom) -> str | None:
    """AutoDock4 atom type used by the Vina scoring function."""
    el = atom.GetSymbol()
    if el == "H":
        nb = atom.GetNeighbors()
        return "HD" if nb and nb[0].GetSymbol() in ("N", "O", "S") else None  # non-polar H are merged
    if el == "C":
        return "A" if atom.GetIsAromatic() else "C"
    if el == "N":
        has_h = any(n.GetSymbol() == "H" for n in atom.GetNeighbors())
        # acceptor nitrogen: no H, not amide/aromatic-with-3-connections
        return "NA" if (not has_h and atom.GetDegree() < 3 and not atom.GetIsAromatic()) or (atom.GetIsAromatic() and atom.GetDegree() == 2 and not has_h) else "N"
    if el == "O":
        return "OA"
    if el == "S":
        return "SA"
    return {"P": "P", "F": "F", "Cl": "Cl", "Br": "Br", "I": "I", "Se": "S", "Mg": "Mg", "Zn": "Zn", "Ca": "Ca", "Fe": "Fe", "Mn": "Mn"}.get(el, el)


def receptor_pdbqt(pdb_text: str) -> str:
    """Rigid-receptor PDBQT: heavy atoms + polar hydrogens with AD4 types (charges unused by Vina)."""
    prot = parse_protein(pdb_text)
    protH = Chem.AddHs(prot, addCoords=True)
    conf = protH.GetConformer()
    out = []
    serial = 0
    for a in protH.GetAtoms():
        t = ad4_type(a)
        if t is None:
            continue
        ri = a.GetPDBResidueInfo()
        if ri is None:  # added hydrogen: borrow residue info from its heavy neighbour
            ri = a.GetNeighbors()[0].GetPDBResidueInfo()
            name = "H"
        else:
            name = ri.GetName().strip()
        p = conf.GetAtomPosition(a.GetIdx())
        serial += 1
        nm = name if len(name) == 4 else f" {name:<3s}"
        out.append(f"ATOM  {serial % 100000:5d} {nm:4s} {ri.GetResidueName():>3s} {ri.GetChainId() or 'A':1s}{ri.GetResidueNumber() % 10000:4d}    "
                   f"{p.x:8.3f}{p.y:8.3f}{p.z:8.3f}  1.00  0.00    {0.0:6.3f} {t:<2s}\n")
    return "".join(out)


# --------------------------------------------------------------------------- ligand
def ligand_pdbqt(gt_mol: Chem.Mol, seed: int = 0):
    """Random ETKDG conformer of the crystal ligand -> Meeko PDBQT string (crystal pose is not the input)."""
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    mol = Chem.AddHs(Chem.RemoveHs(gt_mol))
    cid = AllChem.EmbedMolecule(mol, randomSeed=seed, useRandomCoords=True)
    if cid < 0:
        mol = Chem.AddHs(Chem.RemoveHs(gt_mol), addCoords=True)  # fall back to crystal geometry + H
    else:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    prep = MoleculePreparation()
    setups = prep.prepare(mol)
    s, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise RuntimeError(f"meeko: {err}")
    return s, mol


def poses_to_mols(poses_pdbqt: str, template: Chem.Mol):
    from meeko import PDBQTMolecule, RDKitMolCreate
    pm = PDBQTMolecule(poses_pdbqt, is_dlg=False, skip_typing=True)
    mols = RDKitMolCreate.from_pdbqt_mol(pm)
    mol = mols[0]
    if mol is None:
        raise RuntimeError("meeko could not rebuild the docked ligand")
    mol = Chem.RemoveHs(mol)
    tmpl = Chem.RemoveHs(template)
    # make atom order / bond orders identical to the crystal ligand so downstream matching is trivial
    if mol.GetNumAtoms() == tmpl.GetNumAtoms():
        try:
            mol = AllChem.AssignBondOrdersFromTemplate(tmpl, mol)
        except Exception:  # noqa
            pass
    return mol


# --------------------------------------------------------------------------- per system
def dock_system(args):
    """Pool task: run one system in a child process so a Vina abort() cannot hang the pool.
    Retries with a larger box, which works around Vina's szv_grid assertion on extended ligands."""
    import subprocess
    sid, lig_chain, exhaustiveness, n_poses, overwrite = args
    d = OUT / sid
    res_path = d / "result.json"
    if res_path.exists() and not overwrite:
        return {"system_id": sid, "ok": True, "cached": True}
    d.mkdir(parents=True, exist_ok=True)
    last = None
    for padding in (8.0, 12.0, 16.0):
        cmd = [sys.executable, __file__, "--single", sid, lig_chain, str(exhaustiveness), str(n_poses), str(padding)]
        try:
            pr = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
        except subprocess.TimeoutExpired:
            last = {"system_id": sid, "ok": False, "error": f"timeout after 1800 s (padding {padding})"}
            continue
        if pr.returncode == 0 and res_path.exists():
            return json.loads(res_path.read_text())
        last = {"system_id": sid, "ok": False, "error": f"vina child exit {pr.returncode} (padding {padding}): {(pr.stderr or pr.stdout).strip()[-200:]}"}
    res_path.write_text(json.dumps(last, indent=1))
    return last


def dock_single(sid, lig_chain, exhaustiveness, n_poses, padding=8.0):
    t0 = time.time()
    d = OUT / sid
    res_path = d / "result.json"
    d.mkdir(parents=True, exist_ok=True)
    try:
        from vina import Vina
        gt_dir = RAW / "ground_truth" / sid
        gt = Chem.MolFromMolFile(str(gt_dir / "ligand_files" / f"{lig_chain}.sdf"), removeHs=False)
        pdb = receptor_pdb_from_cif(gt_dir / "receptor.cif")
        rec_pdbqt = receptor_pdbqt(pdb)
        rec_file = d / "receptor.pdbqt"
        rec_file.write_text(rec_pdbqt)
        lig_pdbqt, start_mol = ligand_pdbqt(gt, seed=0)
        xyz = Chem.RemoveHs(gt).GetConformer().GetPositions()
        center = xyz.mean(0)
        extent = xyz.max(0) - xyz.min(0)
        size = np.maximum(extent + 2 * padding, 20.0)  # ligand extent + padding each side, min 20 Å
        v = Vina(sf_name="vina", seed=0, cpu=1, verbosity=0)
        v.set_receptor(rigid_pdbqt_filename=str(rec_file))
        v.set_ligand_from_string(lig_pdbqt)
        v.compute_vina_maps(center=center.tolist(), box_size=size.tolist())
        v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
        energies = v.energies(n_poses=n_poses)
        poses = v.poses(n_poses=n_poses)
        (d / "poses.pdbqt").write_text(poses)
        mol = poses_to_mols(poses, gt)
        mol.SetProp("_Name", f"vina top pose {sid}")
        mol.SetProp("vina_affinity", str(float(energies[0][0])))
        w = Chem.SDWriter(str(d / "pose.sdf")); w.write(mol, confId=0); w.close()
        # all poses as a multi-record SDF too (cheap, useful for later)
        w = Chem.SDWriter(str(d / "poses.sdf"))
        for k in range(mol.GetNumConformers()):
            mol.SetProp("_Name", f"pose {k + 1}"); mol.SetProp("vina_affinity", str(float(energies[k][0]))); w.write(mol, confId=k)
        w.close()
        result = {"system_id": sid, "ok": True, "vina_score": float(energies[0][0]),
                  "pose_scores": [float(e[0]) for e in energies], "box_center": center.round(3).tolist(),
                  "box_size": size.round(1).tolist(), "box_padding": padding, "exhaustiveness": exhaustiveness, "n_receptor_atoms": rec_pdbqt.count("\n"),
                  "elapsed": round(time.time() - t0, 1)}
    except Exception as e:  # noqa
        result = {"system_id": sid, "ok": False, "error": f"{type(e).__name__}: {e}", "trace": traceback.format_exc()[-600:],
                  "elapsed": round(time.time() - t0, 1)}
    res_path.write_text(json.dumps(result, indent=1))
    return result


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--single":
        sid, chain, ex, npz, pad = sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), float(sys.argv[6])
        r = dock_single(sid, chain, ex, npz, pad)
        sys.exit(0 if r["ok"] else 1)
    ap = argparse.ArgumentParser()
    ap.add_argument("--systems", nargs="*")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) - 1))
    ap.add_argument("--exhaustiveness", type=int, default=16)
    ap.add_argument("--n-poses", type=int, default=9)
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()
    import pandas as pd
    subset = json.load(open(HERE / "subset.json"))
    systems = [x for a in (args.systems or []) for x in a.split()] or subset["systems"]
    ann = pd.read_csv(RAW / "annotations.csv", low_memory=False)
    ann = ann[ann["ligand_is_proper"] == True]  # noqa: E712
    chain = {sid: ann[ann.system_id == sid].iloc[0]["ligand_instance_chain"] for sid in systems}
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    jobs = [(sid, chain[sid], args.exhaustiveness, args.n_poses, args.overwrite) for sid in systems]
    n_ok = n_fail = 0
    with Pool(args.workers) as pool, open(OUT / "log.jsonl", "a") as log:
        for r in pool.imap_unordered(dock_system, jobs):
            log.write(json.dumps({k: v for k, v in r.items() if k != "trace"}) + "\n"); log.flush()
            if r["ok"]:
                n_ok += 1
                print(f"ok   {r['system_id']:34s} {r.get('vina_score', 0):7.2f} kcal/mol {r.get('elapsed', 0):6.1f}s{' (cached)' if r.get('cached') else ''}", flush=True)
            else:
                n_fail += 1
                print(f"FAIL {r['system_id']:34s} {r['error']}", flush=True)
    print(f"done: {n_ok} ok, {n_fail} failed, {(time.time() - t0) / 60:.1f} min wall")


if __name__ == "__main__":
    main()
