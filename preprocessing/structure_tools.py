"""Structure utilities for the Runs N' Poses memorisation visualiser.

Everything here works on RDKit molecules (ligands) and gemmi structures (receptors).
Predicted complexes are superposed onto the ground-truth binding site so that every
method's pose can be overlaid in one frame in the viewer.
"""
from __future__ import annotations

import io
import math
import re
from dataclasses import dataclass, field

import gemmi
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, rdMolTransforms
from rdkit.Chem import rdDetermineBonds  # noqa: F401  (import registers functionality)

VDW = {"H": 1.1, "C": 1.7, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.8, "S": 1.8, "Cl": 1.75,
       "Br": 1.85, "I": 1.98, "B": 1.92, "Se": 1.9, "Si": 2.1}
CLASH_RATIO = 0.75  # same threshold PoseBusters uses for minimum protein distance
POCKET_CUTOFF = 8.0  # Å from any ligand heavy atom -> pocket residue (for superposition)
CONTACT_CUTOFF = 4.5  # Å for the residue list shown in the viewer


# --------------------------------------------------------------------------- receptors
@dataclass
class Receptor:
    """Protein heavy atoms of one complex plus per-chain CA tables."""
    structure: gemmi.Structure
    chains: dict = field(default_factory=dict)  # chain name -> {"seq": str, "ca": (n,3), "resnums": [..]}

    @classmethod
    def from_file(cls, path: str, polymer_only: bool = True) -> "Receptor":
        st = gemmi.read_structure(str(path))
        st.setup_entities()
        st.remove_hydrogens()
        st.remove_waters()
        shorten_chain_names(st)
        rec = cls(structure=st)
        model = st[0]
        for ch in model:
            seq, ca, resnums = [], [], []
            for res in ch:
                if res.het_flag == "H" and polymer_only:
                    continue
                info = gemmi.find_tabulated_residue(res.name)
                if info is None or not info.is_amino_acid():
                    continue
                atom = res.find_atom("CA", "*")
                if atom is None:
                    continue
                seq.append(info.one_letter_code.upper() or "X")
                ca.append([atom.pos.x, atom.pos.y, atom.pos.z])
                resnums.append(res.seqid.num)
            if seq:
                rec.chains[ch.name] = {"seq": "".join(seq), "ca": np.array(ca), "resnums": resnums}
        return rec

    def heavy_atoms(self):
        """Yield (element, xyz, chain, resname, resnum, atom_name) for protein heavy atoms."""
        for ch in self.structure[0]:
            for res in ch:
                info = gemmi.find_tabulated_residue(res.name)
                if info is None or not info.is_amino_acid():
                    continue
                for at in res:
                    if at.element.name == "H":
                        continue
                    yield at.element.name, np.array([at.pos.x, at.pos.y, at.pos.z]), ch.name, res.name, res.seqid.num, at.name

    def transform(self, R: np.ndarray, t: np.ndarray):
        for ch in self.structure[0]:
            for res in ch:
                for at in res:
                    p = np.array([at.pos.x, at.pos.y, at.pos.z]) @ R.T + t
                    at.pos = gemmi.Position(*p)
        for c in self.chains.values():
            c["ca"] = c["ca"] @ R.T + t

    def to_pdb(self, protein_only: bool = True) -> str:
        st = self.structure.clone()
        if protein_only:
            for ch in st[0]:
                for i in range(len(ch) - 1, -1, -1):
                    info = gemmi.find_tabulated_residue(ch[i].name)
                    if info is None or not info.is_amino_acid():
                        del ch[i]
        st.remove_empty_chains()
        return st.make_pdb_string()


def shorten_chain_names(st: gemmi.Structure) -> dict:
    """Rename chains to single characters (PDB format) deterministically: '1.A' -> 'A' when free,
    otherwise the next unused letter/digit in file order. Applied once at load time so chain ids in
    the JSON diagnostics equal those in the served PDB files. Returns old -> new mapping."""
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    mapping, used = {}, set()
    for model in st:
        for ch in model:
            if ch.name in mapping:
                ch.name = mapping[ch.name]
                continue
            cand = ch.name.split(".")[-1]
            if not (len(cand) == 1 and cand not in used):
                cand = next(c for c in alphabet if c not in used)
            used.add(cand)
            mapping[ch.name] = cand
            ch.name = cand
    return mapping


def kabsch(P: np.ndarray, Q: np.ndarray):
    """Return R, t minimising |P@R.T + t - Q|."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    return R, qc - pc @ R.T


def align_seqs(a: str, b: str):
    """Global sequence alignment via gemmi; returns list of (i, j) matched index pairs."""
    res = gemmi.align_string_sequences(list(a), list(b), [])
    pairs = []
    i = j = 0
    for op in res.cigar_str().replace("M", "M ").replace("I", "I ").replace("D", "D ").split():
        n, kind = int(op[:-1]), op[-1]
        for _ in range(n):
            if kind == "M":
                pairs.append((i, j)); i += 1; j += 1
            elif kind == "I":
                i += 1
            else:
                j += 1
    return pairs, res.identity() if hasattr(res, "identity") else None


def pocket_residue_indices(rec: Receptor, lig_xyz: np.ndarray, cutoff: float):
    """Per chain, indices of CA table entries whose residue has any heavy atom within cutoff."""
    out = {}
    for ch in rec.structure[0]:
        if ch.name not in rec.chains:
            continue
        resnums = rec.chains[ch.name]["resnums"]
        idx_of = {r: k for k, r in enumerate(resnums)}
        hits = set()
        for res in ch:
            if res.seqid.num not in idx_of:
                continue
            for at in res:
                p = np.array([at.pos.x, at.pos.y, at.pos.z])
                if np.min(np.linalg.norm(lig_xyz - p, axis=1)) <= cutoff:
                    hits.add(idx_of[res.seqid.num]); break
        out[ch.name] = sorted(hits)
    return out


def superpose_on_pocket(pred: Receptor, gt: Receptor, gt_lig_xyz: np.ndarray):
    """Superpose predicted receptor onto the GT binding site (CA atoms within POCKET_CUTOFF).

    Chains are matched by sequence alignment. Returns dict with pocket/global CA RMSD and the
    chain mapping; the predicted receptor is transformed in place.
    """
    pocket = pocket_residue_indices(gt, gt_lig_xyz, POCKET_CUTOFF)
    P, Q, Pall, Qall = [], [], [], []
    mapping = {}
    used = set()
    for gname, g in gt.chains.items():
        best = None
        for pname, p in pred.chains.items():
            if pname in used:
                continue
            pairs, _ = align_seqs(g["seq"], p["seq"])
            ident = sum(g["seq"][i] == p["seq"][j] for i, j in pairs) / max(1, len(g["seq"]))
            if best is None or ident > best[0]:
                best = (ident, pname, pairs)
        if best is None or best[0] < 0.3:
            continue
        ident, pname, pairs = best
        used.add(pname)
        mapping[gname] = pname
        pk = set(pocket.get(gname, []))
        for i, j in pairs:
            Qall.append(g["ca"][i]); Pall.append(pred.chains[pname]["ca"][j])
            if i in pk:
                Q.append(g["ca"][i]); P.append(pred.chains[pname]["ca"][j])
    if len(P) < 3:
        P, Q = Pall, Qall
    if len(P) < 3:
        raise RuntimeError("could not match enough CA atoms for superposition")
    P, Q = np.array(P), np.array(Q)
    R, t = kabsch(P, Q)
    pred.transform(R, t)
    pocket_rmsd = float(np.sqrt(np.mean(np.sum((P @ R.T + t - Q) ** 2, axis=1))))
    Pall, Qall = np.array(Pall), np.array(Qall)
    global_rmsd = float(np.sqrt(np.mean(np.sum((Pall @ R.T + t - Qall) ** 2, axis=1))))
    return {"pocket_ca_rmsd": pocket_rmsd, "global_ca_rmsd_pocket_aligned": global_rmsd,
            "n_pocket_ca": int(len(P)), "chain_map": mapping, "R": R, "t": t}


# --------------------------------------------------------------------------- ligands
def ligand_from_structure(path: str, chain_name: str, template: Chem.Mol) -> Chem.Mol:
    """Extract ligand heavy atoms of one chain from a predicted CIF/PDB and assign bond orders
    from the ground-truth template molecule (same chemical species)."""
    st = gemmi.read_structure(str(path))
    st.remove_hydrogens()
    model = st[0]
    xyz, elems, names = [], [], []
    for ch in model:
        if ch.name != chain_name:
            continue
        for res in ch:
            for at in res:
                if at.element.name in ("H", "D"):
                    continue
                xyz.append([at.pos.x, at.pos.y, at.pos.z])
                elems.append(at.element.name)
                names.append(at.name)
    if not xyz:
        raise RuntimeError(f"no ligand atoms in chain {chain_name}")
    xyz = np.array(xyz)
    # Build a bare molecule with connectivity from distances, then take bond orders from template.
    raw = Chem.RWMol()
    conf = Chem.Conformer(len(xyz))
    for k, (e, p) in enumerate(zip(elems, xyz)):
        a = Chem.Atom(e.capitalize() if len(e) > 1 else e)
        a.SetNoImplicit(True)
        raw.AddAtom(a)
        conf.SetAtomPosition(k, p.tolist())
    raw.AddConformer(conf)
    mol = raw.GetMol()
    rdDetermineBonds.DetermineConnectivity(mol, useVdw=True)
    for a in mol.GetAtoms():
        a.SetNoImplicit(False)
    tmpl = Chem.RemoveHs(template)
    if tmpl.GetNumAtoms() != mol.GetNumAtoms():
        raise RuntimeError(f"atom count mismatch template={tmpl.GetNumAtoms()} model={mol.GetNumAtoms()}")
    mol = AllChem.AssignBondOrdersFromTemplate(tmpl, mol)
    Chem.SanitizeMol(mol)
    Chem.AssignStereochemistryFrom3D(mol)
    return mol


def transform_mol(mol: Chem.Mol, R: np.ndarray, t: np.ndarray):
    conf = mol.GetConformer()
    pos = conf.GetPositions() @ R.T + t
    for i, p in enumerate(pos):
        conf.SetAtomPosition(i, p.tolist())


def symmetry_rmsd(pred: Chem.Mol, ref: Chem.Mol) -> float:
    """Symmetry-corrected heavy-atom RMSD without re-alignment."""
    return float(rdMolAlign.CalcRMS(Chem.RemoveHs(pred), Chem.RemoveHs(ref)))


# --------------------------------------------------------------------------- diagnostics
def _bonded_distance_matrix(mol: Chem.Mol):
    return Chem.GetDistanceMatrix(mol)


def ligand_diagnostics(pred: Chem.Mol, ref: Chem.Mol, protein_atoms: list):
    """Atom-level localisation of geometry problems in a predicted ligand pose.

    Bond lengths/angles are compared with the ground-truth pose of the very same molecule,
    which is an experimentally refined geometry, so deviations are strain the model invented.
    """
    pred = Chem.RemoveHs(pred); ref = Chem.RemoveHs(ref)
    match = ref.GetSubstructMatch(pred, useChirality=False)
    if len(match) != pred.GetNumAtoms():
        match = tuple(range(pred.GetNumAtoms()))
    pc, rc = pred.GetConformer(), ref.GetConformer()
    ppos = pc.GetPositions()

    bonds = []
    for b in pred.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        lp = rdMolTransforms.GetBondLength(pc, i, j)
        lr = rdMolTransforms.GetBondLength(rc, match[i], match[j])
        ratio = lp / lr if lr > 0 else 1.0
        bonds.append({"atoms": [i, j], "pred": round(lp, 3), "ref": round(lr, 3), "ratio": round(ratio, 3),
                      "flag": bool(ratio < 0.85 or ratio > 1.15)})

    angles = []
    for a in pred.GetAtoms():
        nbrs = [n.GetIdx() for n in a.GetNeighbors()]
        for x in range(len(nbrs)):
            for y in range(x + 1, len(nbrs)):
                i, j, k = nbrs[x], a.GetIdx(), nbrs[y]
                ap = rdMolTransforms.GetAngleDeg(pc, i, j, k)
                ar = rdMolTransforms.GetAngleDeg(rc, match[i], match[j], match[k])
                dev = ap - ar
                angles.append({"atoms": [i, j, k], "pred": round(ap, 1), "ref": round(ar, 1),
                               "dev": round(dev, 1), "flag": bool(abs(dev) > 20)})

    # intra-ligand clashes: atoms >= 4 bonds apart closer than CLASH_RATIO * vdW sum
    topo = _bonded_distance_matrix(pred)
    intra = []
    n = pred.GetNumAtoms()
    for i in range(n):
        for j in range(i + 1, n):
            if topo[i, j] < 4:
                continue
            d = float(np.linalg.norm(ppos[i] - ppos[j]))
            lim = CLASH_RATIO * (VDW.get(pred.GetAtomWithIdx(i).GetSymbol(), 1.7) + VDW.get(pred.GetAtomWithIdx(j).GetSymbol(), 1.7))
            if d < lim:
                intra.append({"atoms": [i, j], "dist": round(d, 2), "limit": round(lim, 2)})

    # protein-ligand clashes and contacts
    clashes, contacts = [], {}
    if protein_atoms:
        P = np.array([p[1] for p in protein_atoms])
        for i in range(n):
            d = np.linalg.norm(P - ppos[i], axis=1)
            near = np.where(d < CONTACT_CUTOFF)[0]
            ri = VDW.get(pred.GetAtomWithIdx(i).GetSymbol(), 1.7)
            for k in near:
                el, _, ch, resn, resi, an = protein_atoms[k]
                key = f"{ch}:{resn}{resi}"
                contacts[key] = min(contacts.get(key, 99.0), float(d[k]))
                lim = CLASH_RATIO * (ri + VDW.get(el.capitalize() if len(el) > 1 else el, 1.7))
                if d[k] < lim:
                    clashes.append({"atom": i, "protein": {"chain": ch, "resname": resn, "resnum": int(resi), "atom": an},
                                    "dist": round(float(d[k]), 2), "limit": round(lim, 2)})

    # stereo: compare CIP labels of tetrahedral centres with the ground truth
    Chem.AssignStereochemistryFrom3D(pred); Chem.AssignStereochemistryFrom3D(ref)
    stereo = []
    for i in range(n):
        a, r = pred.GetAtomWithIdx(i), ref.GetAtomWithIdx(match[i])
        if a.GetAtomicNum() != 6:  # PoseBusters only judges carbon centres; P/S labels are ambiguous
            continue
        cp, cr = a.GetPropsAsDict().get("_CIPCode"), r.GetPropsAsDict().get("_CIPCode")
        if cp or cr:
            stereo.append({"atom": i, "pred": cp, "ref": cr, "flag": bool(cp != cr)})

    # aromatic ring flatness
    rings = []
    ri = pred.GetRingInfo()
    for ring in ri.AtomRings():
        if not all(pred.GetAtomWithIdx(k).GetIsAromatic() for k in ring):
            continue
        pts = ppos[list(ring)]
        c = pts.mean(0)
        _, _, vt = np.linalg.svd(pts - c)
        dev = float(np.max(np.abs((pts - c) @ vt[2])))
        rings.append({"atoms": list(ring), "max_dev": round(dev, 3), "flag": bool(dev > 0.25)})

    contact_list = [{"residue": k, "min_dist": round(v, 2)} for k, v in sorted(contacts.items(), key=lambda kv: kv[1])]
    flagged = set()
    for b in bonds:
        if b["flag"]:
            flagged.update(b["atoms"])
    for a in angles:
        if a["flag"]:
            flagged.add(a["atoms"][1])
    for c in intra:
        flagged.update(c["atoms"])
    for c in clashes:
        flagged.add(c["atom"])
    for s in stereo:
        if s["flag"]:
            flagged.add(s["atom"])
    for r in rings:
        if r["flag"]:
            flagged.update(r["atoms"])
    return {"bonds": bonds, "angles": [a for a in angles if a["flag"]], "n_angles": len(angles),
            "intra_clashes": intra, "protein_clashes": clashes, "stereo": stereo, "rings": rings,
            "contacts": contact_list, "flagged_atoms": sorted(flagged),
            "summary": {"bad_bonds": sum(b["flag"] for b in bonds), "bad_angles": sum(a["flag"] for a in angles),
                        "intra_clashes": len(intra), "protein_clashes": len(clashes),
                        "stereo_mismatches": sum(s["flag"] for s in stereo),
                        "nonplanar_rings": sum(r["flag"] for r in rings)}}


# --------------------------------------------------------------------------- minimisation
def _mmff(mol):
    props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
    if props is None:
        raise RuntimeError("MMFF parameters unavailable")
    return props


def _energy(mol, props, pos=None):
    ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=0)
    return float(ff.CalcEnergy(pos.ravel().tolist()) if pos is not None else ff.CalcEnergy())


def minimise_with_trajectory(mol: Chem.Mol, max_frames: int = 60, its_per_frame: int = 15, max_its: int = 2000, n_global: int = 30, seed: int = 42):
    """MMFF94s minimisation of a ligand pose recording a trajectory.

    Returns (summary dict, list of heavy-atom coordinate frames, mol with hydrogens at the local minimum).
    Strain is reported relative to (a) the pose after a restrained relaxation of hydrogens and tiny
    heavy-atom deviations and (b) the nearest local minimum and (c) the global minimum from a
    conformer search. Displacement per heavy atom localises where the pose had to move.
    """
    mh = Chem.AddHs(mol, addCoords=True)
    props = _mmff(mh)
    heavy = [a.GetIdx() for a in mh.GetAtoms() if a.GetAtomicNum() > 1]
    start = mh.GetConformer().GetPositions()

    # (a) relax hydrogens with heavy atoms restrained close to the predicted positions
    ff = AllChem.MMFFGetMoleculeForceField(mh, props)
    for i in heavy:
        ff.MMFFAddPositionConstraint(i, 0.05, 200.0)
    ff.Minimize(maxIts=500)
    e_pose = _energy(mh, props)
    pos_pose = mh.GetConformer().GetPositions()

    # (b) free local minimisation, one frame per `its_per_frame` iterations until convergence
    ff = AllChem.MMFFGetMoleculeForceField(mh, props)
    frames = [pos_pose[heavy].copy()]
    energies = [e_pose]
    converged = False
    total_its = 0
    while total_its < max_its:
        rc = ff.Minimize(maxIts=its_per_frame)
        total_its += its_per_frame
        pos = np.array(ff.Positions()).reshape(-1, 3)
        frames.append(pos[heavy].copy())
        energies.append(float(ff.CalcEnergy()))
        if rc == 0:
            converged = True
            break
    if len(frames) > max_frames:  # keep every k-th frame, always first and last
        keep = sorted(set(np.linspace(0, len(frames) - 1, max_frames).round().astype(int).tolist()))
        frames = [frames[i] for i in keep]
        energies = [energies[i] for i in keep]
    pos_local = np.array(ff.Positions()).reshape(-1, 3)
    conf = mh.GetConformer()
    for i, p in enumerate(pos_local):
        conf.SetAtomPosition(i, p.tolist())
    e_local = float(ff.CalcEnergy())

    # (c) global minimum from a conformer search
    mg = Chem.Mol(mh)
    mg.RemoveAllConformers()
    cids = AllChem.EmbedMultipleConfs(mg, numConfs=n_global, randomSeed=seed)
    e_global = e_local
    if len(cids):
        res = AllChem.MMFFOptimizeMoleculeConfs(mg, mmffVariant="MMFF94s", maxIters=2000)
        e_global = min(min(e for _, e in res), e_local)

    disp = np.linalg.norm(pos_local[heavy] - start[heavy], axis=1)
    rmsd_drift = float(np.sqrt(np.mean(disp ** 2)))
    summary = {
        "e_pose": round(e_pose, 2), "e_local": round(e_local, 2), "e_global": round(e_global, 2),
        "strain_local": round(e_pose - e_local, 2), "strain_global": round(e_pose - e_global, 2),
        "rmsd_drift": round(rmsd_drift, 3), "max_atom_displacement": round(float(disp.max()), 3),
        "atom_displacement": [round(float(d), 3) for d in disp],
        "n_frames": len(frames), "converged": converged, "total_iterations": total_its, "energies": [round(e, 2) for e in energies],
        "frame_rmsd": [round(float(np.sqrt(np.mean(np.sum((f - frames[0]) ** 2, axis=1)))), 3) for f in frames],
    }
    return summary, frames, mh


def frames_to_sdf(mol_heavy: Chem.Mol, frames: list, energies: list) -> str:
    """Multi-record SDF with one heavy-atom conformer per minimisation frame."""
    m = Chem.RemoveHs(mol_heavy)
    buf = io.StringIO()
    w = Chem.SDWriter(buf)
    for k, (f, e) in enumerate(zip(frames, energies)):
        conf = m.GetConformer()
        for i, p in enumerate(f):
            conf.SetAtomPosition(i, p.tolist())
        m.SetProp("_Name", f"frame {k}")
        m.SetProp("energy", str(e))
        w.write(m)
    w.close()
    return buf.getvalue()


def mol_to_sdf(mol: Chem.Mol, name: str = "") -> str:
    m = Chem.Mol(mol)
    if name:
        m.SetProp("_Name", name)
    return Chem.MolToMolBlock(m) + "$$$$\n"


# --------------------------------------------------------------------------- pocket minimisation
def _count_clashes(lig_xyz, lig_elems, pk_xyz, pk_elems) -> int:
    n = 0
    for i, p in enumerate(lig_xyz):
        d = np.linalg.norm(pk_xyz - p, axis=1)
        ri = VDW.get(lig_elems[i], 1.7)
        lim = CLASH_RATIO * (ri + np.array([VDW.get(e, 1.7) for e in pk_elems]))
        n += int(np.sum(d < lim))
    return n


def _pdb_atom_line(serial, name, resname, chain, resnum, icode, xyz, elem):
    nm = name if len(name) == 4 else f" {name:<3s}"
    return (f"ATOM  {serial:5d} {nm:4s} {resname:>3s} {chain:1s}{resnum:4d}{icode:1s}   "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {elem:>2s}\n")


def pocket_frames_to_pdb(pocket_info: list, frames: list) -> str:
    """Multi-MODEL PDB of pocket heavy atoms, one MODEL per minimisation frame.

    pocket_info: list of (name, resname, chain, resnum, icode, element) per pocket heavy atom.
    """
    out = []
    for k, f in enumerate(frames):
        out.append(f"MODEL     {k + 1:4d}\n")
        for s, (info, xyz) in enumerate(zip(pocket_info, f), start=1):
            out.append(_pdb_atom_line(s, *info[:5], xyz, info[5]))
        out.append("ENDMDL\n")
    out.append("END\n")
    return "".join(out)


def _parse_protein(pdb_text: str):
    """PDB -> RDKit mol with standard-residue bond orders. Proximity bonding occasionally adds a
    spurious bond between clashing atoms in predicted models (e.g. a backbone O gaining a third
    neighbour); such bonds are removed (longest first) until the molecule sanitises."""
    mol = Chem.MolFromPDBBlock(pdb_text, removeHs=False, sanitize=False, proximityBonding=True)
    if mol is None:
        return None
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
    try:
        mol = rw.GetMol()
        Chem.SanitizeMol(mol)
        return mol
    except Exception:  # noqa
        return None


def _extract_pocket(prot: Chem.Mol, keep: list):
    """Sub-molecule of the kept atoms that sanitises. Hydrogen counts are frozen from the sanitised
    full protein so aromatic HIS/TRP rings kekulise unambiguously; if sanitisation still fails, retry
    without kekulisation, then drop the offending residue(s). Returns (pocket, keep, fallback)."""
    frozen = Chem.RWMol(prot)
    for a in frozen.GetAtoms():
        a.SetNumExplicitHs(a.GetTotalNumHs()); a.SetNoImplicit(True)
    keep = list(keep)
    for attempt in range(6):
        em = Chem.RWMol(frozen)
        for idx in sorted(set(range(prot.GetNumAtoms())) - set(keep), reverse=True):
            em.RemoveAtom(idx)
        pocket = em.GetMol()
        pocket.UpdatePropertyCache(strict=False)
        try:
            Chem.SanitizeMol(pocket)
            return pocket, keep, None if attempt == 0 else f"dropped_residues:{attempt}"
        except Chem.KekulizeException as e:
            try:  # partial sanitisation without kekulisation
                p2 = Chem.Mol(pocket)
                Chem.SanitizeMol(p2, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE ^ Chem.SANITIZE_SETAROMATICITY)
                Chem.SetAromaticity(p2)
                if AllChem.MMFFHasAllMoleculeParams(p2):
                    return p2, keep, "no_kekulize"
            except Exception:  # noqa
                pass
            bad = [int(x) for x in re.findall(r"\d+", str(e).split("atoms:")[-1])] if "atoms:" in str(e) else []
            if not bad:
                raise
            # map pocket atom indices back to the full-protein indices and drop their residues
            res = set()
            for i in bad:
                ri = pocket.GetAtomWithIdx(i).GetPDBResidueInfo()
                res.add((ri.GetChainId(), ri.GetResidueNumber(), ri.GetInsertionCode()))
            keep = [k for k in keep if (prot.GetAtomWithIdx(k).GetPDBResidueInfo().GetChainId(),
                                        prot.GetAtomWithIdx(k).GetPDBResidueInfo().GetResidueNumber(),
                                        prot.GetAtomWithIdx(k).GetPDBResidueInfo().GetInsertionCode()) not in res]
    raise RuntimeError("pocket could not be sanitised")


def minimise_in_pocket(lig_mol: Chem.Mol, receptor_pdb_text: str, cutoff: float = 8.0, restrain_cut: float = 6.0,
                       its_per_frame: int = 40, max_its: int = 1000):
    """Restrained minimisation of the ligand inside its binding pocket (MMFF94s, UFF fallback).

    Whole residues with any heavy atom within `cutoff` Å of the ligand are kept. Backbone atoms and
    pocket heavy atoms further than `restrain_cut` Å from the ligand are fixed; other side-chain
    heavy atoms get a flat-bottom position restraint; the ligand and all hydrogens are free.
    Hydrogens are added with coordinates because the source structures carry none.

    Returns (summary, ligand_frames, pocket_frames, pocket_info, ligand_with_Hs).
    """
    lig = Chem.RemoveHs(lig_mol)
    prot = _parse_protein(receptor_pdb_text)
    if prot is None:
        return {"ok": False, "error": "protein PDB block could not be parsed"}, None, None, None, None
    lig_xyz = lig.GetConformer().GetPositions()
    pxyz = prot.GetConformer().GetPositions()
    d = np.min(np.linalg.norm(pxyz[:, None, :] - lig_xyz[None, :, :], axis=2), axis=1)
    keep_res = set()
    for a in prot.GetAtoms():
        if d[a.GetIdx()] <= cutoff:
            ri = a.GetPDBResidueInfo()
            keep_res.add((ri.GetChainId(), ri.GetResidueNumber(), ri.GetInsertionCode()))
    keep = [a.GetIdx() for a in prot.GetAtoms()
            if (a.GetPDBResidueInfo().GetChainId(), a.GetPDBResidueInfo().GetResidueNumber(), a.GetPDBResidueInfo().GetInsertionCode()) in keep_res]
    if len(keep) < 10:
        return {"ok": False, "error": "no pocket residues within cutoff"}, None, None, None, None
    pocket, keep, parse_fallback = _extract_pocket(prot, keep)
    d_keep = d[keep]
    res_labels = []
    seen = set()
    for a in pocket.GetAtoms():
        ri = a.GetPDBResidueInfo()
        key = (ri.GetChainId(), ri.GetResidueNumber(), ri.GetInsertionCode())
        if key not in seen:
            seen.add(key)
            res_labels.append(f"{ri.GetChainId()}:{ri.GetResidueName().strip()}{ri.GetResidueNumber()}")
    pocketH = Chem.AddHs(pocket, addCoords=True)
    ligH = Chem.AddHs(lig, addCoords=True)
    complex_ = Chem.CombineMols(pocketH, ligH)
    Chem.SanitizeMol(complex_)
    n_p = pocketH.GetNumAtoms()
    heavy_p = [a.GetIdx() for a in pocketH.GetAtoms() if a.GetAtomicNum() > 1]
    heavy_l = [a.GetIdx() for a in ligH.GetAtoms() if a.GetAtomicNum() > 1]
    pk_elems = [pocketH.GetAtomWithIdx(i).GetSymbol() for i in heavy_p]
    lg_elems = [ligH.GetAtomWithIdx(i).GetSymbol() for i in heavy_l]
    pocket_info = []
    for i in heavy_p:
        ri = pocketH.GetAtomWithIdx(i).GetPDBResidueInfo()
        pocket_info.append((ri.GetName().strip(), ri.GetResidueName().strip(), ri.GetChainId() or "A", ri.GetResidueNumber(),
                            ri.GetInsertionCode() or " ", pocketH.GetAtomWithIdx(i).GetSymbol()))

    # ---- force field factory (MMFF94s, distance-dependent dielectric 4; UFF fallback)
    ff_kind = None
    props = AllChem.MMFFGetMoleculeProperties(complex_, mmffVariant="MMFF94s")
    if props is not None:
        ff_kind = "MMFF94s"
        if hasattr(AllChem, "MMFFDielectricModel"):
            props.SetMMFFDielectricModel(AllChem.MMFFDielectricModel.MMFFDistDielectric)
        props.SetMMFFDielectricConstant(4.0)

        def make_ff(mol=complex_, p=props):
            return AllChem.MMFFGetMoleculeForceField(mol, p, nonBondedThresh=9.0, ignoreInterfragInteractions=False)

        def frag_energy(m):
            pp = AllChem.MMFFGetMoleculeProperties(m, mmffVariant="MMFF94s")
            if hasattr(AllChem, "MMFFDielectricModel"):
                pp.SetMMFFDielectricModel(AllChem.MMFFDielectricModel.MMFFDistDielectric)
            pp.SetMMFFDielectricConstant(4.0)
            return AllChem.MMFFGetMoleculeForceField(m, pp, nonBondedThresh=9.0, ignoreInterfragInteractions=False).CalcEnergy()

        def add_restraint(ff, i):
            ff.MMFFAddPositionConstraint(i, 0.3, 5.0)
    elif AllChem.UFFHasAllMoleculeParams(complex_):
        ff_kind = "UFF"

        def make_ff(mol=complex_):
            return AllChem.UFFGetMoleculeForceField(mol, vdwThresh=9.0, ignoreInterfragInteractions=False)

        def frag_energy(m):
            return AllChem.UFFGetMoleculeForceField(m, vdwThresh=9.0, ignoreInterfragInteractions=False).CalcEnergy()

        def add_restraint(ff, i):
            ff.UFFAddPositionConstraint(i, 0.3, 5.0)
    else:
        return {"ok": False, "error": "neither MMFF94s nor UFF could type the pocket–ligand complex",
                "n_pocket_residues": len(res_labels), "n_pocket_atoms": len(heavy_p), "pocket_residues": res_labels}, None, None, None, None

    def set_positions(pos):
        conf = complex_.GetConformer()
        for i, p in enumerate(pos):
            conf.SetAtomPosition(i, p.tolist())

    def frag_energies(pos):
        pH, lH = Chem.Mol(pocketH), Chem.Mol(ligH)
        cp, cl = pH.GetConformer(), lH.GetConformer()
        for i in range(n_p):
            cp.SetAtomPosition(i, pos[i].tolist())
        for i in range(ligH.GetNumAtoms()):
            cl.SetAtomPosition(i, pos[n_p + i].tolist())
        return frag_energy(pH), frag_energy(lH)

    try:
        # ---- stage 0: hydrogen-only relaxation (all heavy atoms fixed)
        ff = make_ff()
        for i in heavy_p:
            ff.AddFixedPoint(i)
        for i in heavy_l:
            ff.AddFixedPoint(n_p + i)
        ff.Initialize()
        ff.Minimize(maxIts=200)
        pos0 = np.array(ff.Positions()).reshape(-1, 3)
        set_positions(pos0)
        plain = make_ff()  # unrestrained, for reporting energies
        plain.Initialize()
        e_pose = float(plain.CalcEnergy())
        ep, el = frag_energies(pos0)
        e_int_pose = e_pose - ep - el

        # ---- stage 1: restrained pocket + free ligand minimisation with trajectory
        ff = make_ff()
        n_fixed = n_restr = 0
        movable_res = set()
        for k, i in enumerate(heavy_p):
            name = pocket_info[k][0]
            if name in ("N", "CA", "C", "O") or d_keep[i] > restrain_cut:
                ff.AddFixedPoint(i); n_fixed += 1
            else:
                add_restraint(ff, i); n_restr += 1
                movable_res.add(pocket_info[k][1:5])
        ff.Initialize()
        lig_frames = [pos0[n_p:][heavy_l].copy()]
        pk_frames = [pos0[heavy_p].copy()]
        energies = [e_pose]
        converged = False
        total_its = 0
        while total_its < max_its:  # uniform sampling: one frame per `its_per_frame` iterations
            rc = ff.Minimize(maxIts=its_per_frame)
            total_its += its_per_frame
            pos = np.array(ff.Positions()).reshape(-1, 3)
            lig_frames.append(pos[n_p:][heavy_l].copy()); pk_frames.append(pos[heavy_p].copy())
            energies.append(float(plain.CalcEnergy(pos.ravel().tolist())))
            if rc == 0:
                converged = True
                break
        pos1 = np.array(ff.Positions()).reshape(-1, 3)
        e_min = float(plain.CalcEnergy(pos1.ravel().tolist()))
        ep1, el1 = frag_energies(pos1)
        e_int_min = e_min - ep1 - el1
    except Exception as e:  # noqa
        return {"ok": False, "error": f"{ff_kind}: {type(e).__name__}: {e}", "n_pocket_residues": len(res_labels),
                "n_pocket_atoms": len(heavy_p), "pocket_residues": res_labels}, None, None, None, None

    # trim the pocket trajectory to residues that have at least one movable atom (the rest is static)
    keep_idx = [k for k, info in enumerate(pocket_info) if info[1:5] in movable_res]
    pocket_info = [pocket_info[k] for k in keep_idx]
    pk_frames_out = [f[keep_idx] for f in pk_frames]

    lig0, lig1 = lig_frames[0], lig_frames[-1]
    disp = np.linalg.norm(lig1 - lig0, axis=1)
    pk_disp = np.linalg.norm(pk_frames[-1] - pk_frames[0], axis=1)
    summary = {
        "ok": True, "force_field": ff_kind, "pocket_parse_fallback": parse_fallback,
        "n_pocket_residues": len(res_labels), "n_pocket_atoms": len(heavy_p), "pocket_residues": res_labels,
        "n_fixed": n_fixed, "n_restrained": n_restr,
        "e_complex_pose": round(e_pose, 2), "e_complex_min": round(e_min, 2),
        "e_interaction_pose": round(e_int_pose, 2), "e_interaction_min": round(e_int_min, 2),
        "ligand_rmsd_drift": round(float(np.sqrt(np.mean(disp ** 2))), 3),
        "ligand_atom_displacement": [round(float(x), 3) for x in disp],
        "pocket_heavy_rmsd": round(float(np.sqrt(np.mean(pk_disp ** 2))), 3),
        "max_pocket_atom_displacement": round(float(pk_disp.max()), 3),
        "clashes_pose": _count_clashes(lig0, lg_elems, pk_frames[0], pk_elems),
        "clashes_min": _count_clashes(lig1, lg_elems, pk_frames[-1], pk_elems),
        "n_frames": len(lig_frames), "total_iterations": total_its, "energies": [round(e, 2) for e in energies],
        "frame_ligand_rmsd": [round(float(np.sqrt(np.mean(np.sum((f - lig0) ** 2, axis=1)))), 3) for f in lig_frames],
        "converged": bool(converged),
    }
    return summary, lig_frames, pk_frames_out, pocket_info, ligH


# --------------------------------------------------------------------------- contacts & pocket hit
from rdkit.Chem import Lipinski, rdShapeHelpers  # noqa: E402

CONTACT_DIST = 4.0
HBOND_DIST = 3.5
PROTEIN_IONIC = {("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2"), ("LYS", "NZ"),
                 ("ARG", "NH1"), ("ARG", "NH2"), ("ARG", "NE"), ("HIS", "ND1"), ("HIS", "NE2")}
_IONIC_SMARTS = [Chem.MolFromSmarts(x) for x in (
    "[CX3](=O)[O-,OX2H1]",            # carboxylate / carboxylic acid
    "[CX3](=[NX2,NX3+])[NX3]",         # amidine / guanidine
    "[NX3;H2,H1,H0;!$(NC=O);!$(N-a);!$(N-[SX4]);!$(NC=N)]([CX4])",  # basic aliphatic amine
    "[PX4](=O)([O-,OX2H1])",           # phosphate / phosphonate
    "[SX4](=O)(=O)[O-,OX2H1]",         # sulfonate / sulfate
    "[+,-]",                           # explicit formal charge
)]


def best_atom_mapping(pred: Chem.Mol, ref: Chem.Mol) -> tuple:
    """pred atom index -> ref atom index; among symmetry-equivalent matches pick the one with the
    lowest RMSD on current coordinates (both heavy-atom molecules, same species)."""
    pred, ref = Chem.RemoveHs(pred), Chem.RemoveHs(ref)
    matches = ref.GetSubstructMatches(pred, uniquify=False, useChirality=False, maxMatches=2000)
    if not matches:
        return tuple(range(pred.GetNumAtoms()))
    P, R = pred.GetConformer().GetPositions(), ref.GetConformer().GetPositions()
    best, best_rmsd = matches[0], None
    for m in matches:
        r = float(np.sqrt(np.mean(np.sum((P - R[list(m)]) ** 2, axis=1))))
        if best_rmsd is None or r < best_rmsd:
            best, best_rmsd = m, r
    return best


def _ligand_atom_classes(mol: Chem.Mol):
    """Per heavy atom: polar (H-bond donor/acceptor N,O), hydrophobic C, ionic-group member."""
    mol = Chem.RemoveHs(mol)
    n = mol.GetNumAtoms()
    polar = set()
    for patt in (Lipinski.HDonorSmarts, Lipinski.HAcceptorSmarts):
        for m in mol.GetSubstructMatches(patt):
            polar.update(i for i in m if mol.GetAtomWithIdx(i).GetSymbol() in ("N", "O"))
    hydrophobic = {a.GetIdx() for a in mol.GetAtoms()
                   if a.GetSymbol() == "C" and not any(nb.GetSymbol() in ("N", "O") for nb in a.GetNeighbors())}
    ionic = set()
    for patt in _IONIC_SMARTS:
        for m in mol.GetSubstructMatches(patt):
            ionic.update(i for i in m if mol.GetAtomWithIdx(i).GetSymbol() in ("N", "O", "P", "S"))
    return {"polar": polar, "hydrophobic": hydrophobic, "ionic": ionic, "n": n}


def ligand_contacts(lig_xyz: np.ndarray, lig_elems: list, classes: dict, protein_atoms: list) -> dict:
    """{(ligand atom index, residue label): set(types)} for heavy-atom pairs within CONTACT_DIST."""
    if not protein_atoms:
        return {}
    P = np.array([p[1] for p in protein_atoms])
    P_el = [p[0] for p in protein_atoms]
    P_res = [f"{p[2]}:{p[3]}{p[4]}" for p in protein_atoms]
    P_ionic = [(p[3], p[5]) in PROTEIN_IONIC for p in protein_atoms]
    out = {}
    for i, x in enumerate(lig_xyz):
        d = np.linalg.norm(P - x, axis=1)
        for k in np.where(d <= CONTACT_DIST)[0]:
            key = (i, P_res[k])
            types = out.setdefault(key, {"any"})
            el = P_el[k]
            if i in classes["polar"] and el in ("N", "O") and d[k] <= HBOND_DIST:
                types.add("hbond")
            if i in classes["hydrophobic"] and el == "C":
                types.add("hydrophobic")
            if i in classes["ionic"] and P_ionic[k]:
                types.add("ionic")
    return out


def _specific_type(types: set) -> str:
    for t in ("ionic", "hbond", "hydrophobic"):
        if t in types:
            return t
    return "any"


def gt_contact_summary(gt_mol: Chem.Mol, protein_atoms: list) -> dict:
    gt = Chem.RemoveHs(gt_mol)
    classes = _ligand_atom_classes(gt)
    c = ligand_contacts(gt.GetConformer().GetPositions(), [a.GetSymbol() for a in gt.GetAtoms()], classes, protein_atoms)
    by_type = {t: sum(1 for v in c.values() if t in v) for t in ("any", "hbond", "hydrophobic", "ionic")}
    return {"cutoff": CONTACT_DIST, "total": len(c), "by_type": by_type, "residues": sorted({r for _, r in c})}


def contact_retention(pred_mol: Chem.Mol, gt_mol: Chem.Mol, protein_atoms: list) -> dict:
    """Which crystal-ligand contacts (ligand atom, residue) the predicted pose reproduces, both
    measured against the crystal receptor in the crystal frame."""
    gt, pred = Chem.RemoveHs(gt_mol), Chem.RemoveHs(pred_mol)
    classes = _ligand_atom_classes(gt)
    gt_c = ligand_contacts(gt.GetConformer().GetPositions(), [a.GetSymbol() for a in gt.GetAtoms()], classes, protein_atoms)
    mapping = best_atom_mapping(pred, gt)  # pred idx -> gt idx
    ppos = pred.GetConformer().GetPositions()
    xyz = np.zeros_like(gt.GetConformer().GetPositions())
    for i, j in enumerate(mapping):
        xyz[j] = ppos[i]
    pred_c = ligand_contacts(xyz, [a.GetSymbol() for a in gt.GetAtoms()], classes, protein_atoms)
    types = ("any", "hbond", "hydrophobic", "ionic")
    by_type = {t: {"gt": sum(1 for v in gt_c.values() if t in v),
                   "kept": sum(1 for k, v in gt_c.items() if t in v and k in pred_c and t in pred_c[k])} for t in types}
    residues_gt = {r for _, r in gt_c}
    residues_pred = {r for _, r in pred_c}
    return {
        "cutoff": CONTACT_DIST, "gt_total": len(gt_c), "kept": by_type["any"]["kept"],
        "retention": (round(by_type["any"]["kept"] / len(gt_c), 3) if gt_c else None),
        "by_type": by_type,
        "residues": [{"residue": r, "gt": r in residues_gt, "pred": r in residues_pred} for r in sorted(residues_gt | residues_pred)],
        "lost": [{"ligand_atom": int(k[0]), "residue": k[1], "type": _specific_type(v)} for k, v in sorted(gt_c.items()) if k not in pred_c],
        "new": [{"ligand_atom": int(k[0]), "residue": k[1], "type": _specific_type(v)} for k, v in sorted(pred_c.items()) if k not in gt_c],
    }


def pocket_hit(pred_mol: Chem.Mol, gt_mol: Chem.Mol) -> dict:
    """Did the pose land in the right pocket at all? Shape overlap on current coordinates + centroid distance."""
    pred, gt = Chem.RemoveHs(pred_mol), Chem.RemoveHs(gt_mol)
    try:
        overlap = 1.0 - float(rdShapeHelpers.ShapeTanimotoDist(pred, gt))
    except Exception:  # noqa
        overlap = 0.0
    cd = float(np.linalg.norm(pred.GetConformer().GetPositions().mean(0) - gt.GetConformer().GetPositions().mean(0)))
    return {"shape_overlap": round(overlap, 3), "centroid_distance": round(cd, 2), "hit": bool(overlap >= 0.05 or cd <= 4.0)}
