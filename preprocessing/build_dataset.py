"""Build the demo dataset served by the web app.

For every subset system and every method, the top-ranked prediction is:
  1. superposed onto the ground-truth binding site (CA atoms within 8 Å of the ligand),
  2. scored with a symmetry-corrected ligand RMSD (plus the benchmark's own BiSyRMSD / lDDT-PLI),
  3. checked with PoseBusters (redock config) and an atom-level diagnostic pass,
  4. minimised with MMFF94s to get strain energies and a short trajectory.

Outputs (api/data):
  index.json                       overview table (systems, per-method results, similarity per cutoff)
  systems/<system_id>.json         per-system detail (all methods)
  structures/<system_id>/*.gz      gzipped PDB/SDF files for the viewer

Usage:
  python build_dataset.py [--systems ID ...] [--methods M ...] [--workers N]
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

import structure_tools as st

RDLogger.DisableLog("rdApp.*")

HERE = Path(__file__).parent
RAW = HERE / "raw"
PRED_ROOT = RAW / "prediction_files"
GT_ROOT = RAW / "ground_truth"
OUT = HERE.parent / "api" / "data"
WORK = HERE / "work"

CUTOFFS = ["2018-01-01", "2019-01-01", "2020-01-01", "2021-01-01", "2021-09-30",
           "2022-06-01", "2023-01-01", "2023-06-01", "2024-01-01"]
DEFAULT_CUTOFF = "2021-09-30"

METHODS = [
    # id, display name, training cutoff, colour, dir name in prediction_files, predictions csv name
    {"id": "af3", "name": "AF3", "kind": "cofolding", "training_cutoff": "2021-09-30", "color": "#5e6ad2", "dir": "af3", "csv": "af3", "note": None},
    {"id": "boltz", "name": "Boltz-1", "kind": "cofolding", "training_cutoff": "2021-09-30", "color": "#d98c1c", "dir": "boltz", "csv": "boltz", "note": None},
    {"id": "boltz1x", "name": "Boltz-1x", "kind": "cofolding", "training_cutoff": "2021-09-30", "color": "#e8b45a", "dir": "boltz1x", "csv": "boltz1x", "note": None},
    {"id": "vina", "name": "Vina (rigid holo redocking)", "kind": "docking", "training_cutoff": None, "color": "#7a7a85", "dir": None, "csv": None,
     "note": "Rigid receptor, ligand redocked into the crystal (holo) pocket — an easier problem than co-folding."},
]
# Other Runs N' Poses methods are available in the raw data but are not part of the demo:
EXTRA_METHODS = [
    {"id": "af3_no_template", "name": "AF3 (no templates)", "kind": "cofolding", "training_cutoff": "2021-09-30", "color": "#9aa3e6", "dir": "af3_no_template", "csv": "af3_no_template", "note": None},
    {"id": "boltz2", "name": "Boltz-2", "kind": "cofolding", "training_cutoff": "2023-06-01", "color": "#b8651b", "dir": "boltz2", "csv": "boltz2", "note": None},
    {"id": "chai", "name": "Chai-1", "kind": "cofolding", "training_cutoff": "2021-09-30", "color": "#2f9e9e", "dir": "chai", "csv": "chai", "note": None},
    {"id": "protenix", "name": "Protenix", "kind": "cofolding", "training_cutoff": "2021-09-30", "color": "#c04a8a", "dir": "protenix", "csv": "protenix", "note": None},
    {"id": "rfaa", "name": "RFAA", "kind": "cofolding", "training_cutoff": "2021-09-30", "color": "#7a7a85", "dir": "rfaa", "csv": "rfaa", "note": None},
]
DOCKED = RAW / "docked"
PB_EXCLUDE = {"rmsd_≤_2å", "file", "molecule", "position"}


# --------------------------------------------------------------------------- inputs
def load_inputs():
    ann = pd.read_csv(RAW / "annotations.csv", low_memory=False)
    ann = ann[ann["ligand_is_proper"] == True]  # noqa: E712
    preds = {}
    for m in METHODS:
        if not m["csv"]:
            continue
        p = RAW / "predictions" / "predictions" / f"{m['csv']}.csv"
        if p.exists():
            df = pd.read_csv(p, low_memory=False)
            preds[m["id"]] = df[df["ligand_is_proper"] == True] if "ligand_is_proper" in df else df  # noqa: E712
    sim = pd.read_parquet(RAW / "subset_similarity_scores.parquet") if (RAW / "subset_similarity_scores.parquet").exists() else None
    return ann, preds, sim


def similarity_by_cutoff(sim: pd.DataFrame | None, system_id: str, column: str = "sucos_shape_pocket_qcov") -> dict:
    """Max similarity (per metric column) to any structure released before each cutoff date."""
    out = {}
    if sim is None or column not in sim:
        return out
    s = sim[sim["query_system"] == system_id]
    # a system's own PDB entry (and symmetry mates) is a trivial 100% hit once the cutoff passes its
    # release date; drop those so later cutoffs measure similarity to *other* structures
    s = s[s["target_system"].str[:4] != system_id[:4]]
    for c in CUTOFFS:
        t = s[s["target_release_date"] < pd.Timestamp(c)]
        v = float(t[column].max()) if len(t) else float("nan")
        out[c] = 0.0 if np.isnan(v) else round(v, 2)
    return out


def _sim_val(a, col, fallback=None):
    v = a.get(col) if hasattr(a, "get") else a[col]
    if v is None or pd.isna(v):
        if fallback is not None and fallback in a and not pd.isna(a[fallback]):
            return round(float(a[fallback]), 2)
        return 0.0
    return round(float(v), 2)


def system_meta(sid: str, ann: pd.DataFrame, inputs: dict, sim: pd.DataFrame | None) -> dict:
    a = ann[ann["system_id"] == sid].iloc[0]
    return {
        "system_id": sid, "pdb_id": a["entry_pdb_id"], "ligand_smiles": a["ligand_smiles"], "ccd": str(a["ligand_ccd_code"]),
        "ligand_chain": a["ligand_instance_chain"], "n_heavy": int(a["ligand_num_heavy_atoms"]),
        "seq_len": int(sum(len(x) for x in inputs[sid]["sequences"].values())), "n_protein_chains": int(a["num_protein_chains"]),
        "release_date": str(a["release_date"]), "cluster": str(a["cluster"]),
        "closest_training": {"system_id": None if pd.isna(a["target_system"]) else a["target_system"],
                             "pdb_id": None if pd.isna(a["target_system"]) else str(a["target_system"])[:4],
                             "release_date": None if pd.isna(a["target_release_date"]) else str(a["target_release_date"])[:10]},
        # SuCOS-pocket (ligand shape/colour x pocket coverage): the benchmark's headline similarity
        "similarity": _sim_val(a, "sucos_shape_pocket_qcov"),
        "similarity_by_cutoff": similarity_by_cutoff(sim, sid, "sucos_shape_pocket_qcov"),
        # ligand-only (Morgan Tanimoto) and pocket-only (pocket_qcov) views of the same question
        "ligand_similarity": _sim_val(a, "morgan_tanimoto"),
        "ligand_similarity_by_cutoff": similarity_by_cutoff(sim, sid, "morgan_tanimoto"),
        "pocket_similarity": _sim_val(a, "pocket_qcov", fallback="pocket_qcov_alone"),
        "pocket_similarity_by_cutoff": similarity_by_cutoff(sim, sid, "pocket_qcov"),
        "num_training_systems_with_similar_ccds": 0 if pd.isna(a["num_training_systems_with_similar_ccds"]) else int(a["num_training_systems_with_similar_ccds"]),
    }


META_KEYS = ("system_id", "pdb_id", "ligand_smiles", "ccd", "ligand_chain", "n_heavy", "seq_len", "n_protein_chains", "release_date",
             "cluster", "closest_training", "similarity", "similarity_by_cutoff", "ligand_similarity", "ligand_similarity_by_cutoff",
             "pocket_similarity", "pocket_similarity_by_cutoff", "num_training_systems_with_similar_ccds")


def top_ranked(df: pd.DataFrame, system_id: str, ligand_chain: str):
    rows = df[(df["target"] == system_id) & (df["ligand_instance_chain"] == ligand_chain)]
    if rows.empty:
        return None
    if "ranking_score" in rows and rows["ranking_score"].notna().any():
        rows = rows.sort_values("ranking_score", ascending=False)
    return rows.iloc[0]


def find_model_file(method: dict, system_id: str, seed, sample) -> Path | None:
    """Locate the structure file for a given seed/sample; layouts differ per method."""
    d = PRED_ROOT / method["dir"] / system_id
    if not d.exists():
        return None
    files = [p for p in d.rglob("*") if p.suffix in (".cif", ".pdb")]
    seed_s, sample_s = str(seed), str(sample)
    patterns = [
        rf"seed-{seed_s}_sample-{sample_s}\.cif$",                      # af3
        rf"seed[_-]{seed_s}.*model[_-]?(idx[_-])?{sample_s}\.(cif|pdb)$",  # chai / protenix / boltz variants
        rf"{seed_s}[/_].*model[_-]{sample_s}\.(cif|pdb)$",
        rf"seed[_-]{seed_s}.*sample[_-]{sample_s}\.(cif|pdb)$",
        rf"{seed_s}.*_{sample_s}\.(cif|pdb)$",
    ]
    for pat in patterns:
        for p in files:
            if re.search(pat, str(p.relative_to(d))):
                return p
    if len(files) == 1:  # single-model methods (rfaa)
        return files[0]
    return None


def ligand_chain_in_model(row, path: Path, gt_mol, n_prot_chains: int) -> str:
    """Chain id of the ligand in the model file: use the benchmark's mapping when given, otherwise
    pick the hetero chain whose heavy-atom count matches the template."""
    for col in ("model_ligand_chain_rmsd", "model_ligand_chain_lddt_pli", "model_ligand_chain"):
        if col in row and isinstance(row[col], str) and row[col]:
            return row[col]
    import gemmi
    s = gemmi.read_structure(str(path))
    n = Chem.RemoveHs(gt_mol).GetNumAtoms()
    for ch in s[0]:
        heavy = sum(1 for r in ch for a in r if a.element.name != "H" and gemmi.find_tabulated_residue(r.name) is None)
        if heavy == n:
            return ch.name
    raise RuntimeError("ligand chain not found")


# --------------------------------------------------------------------------- per-system work
def gz_write(path: Path, text: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", compresslevel=9) as f:
        f.write(text)


def process_system(system_id: str, meta: dict, method_rows: dict, method_list: list, reuse: bool = True) -> dict:
    """Runs in a worker process. Returns the per-system detail dict (also written to disk).

    With reuse=True, per-method results already present in systems/<id>.json are kept (so the
    build can run incrementally while prediction files are still being extracted)."""
    t0 = time.time()
    existing = {}
    prev_path = OUT / "systems" / f"{system_id}.json"
    if reuse and prev_path.exists():
        try:
            existing = json.loads(prev_path.read_text()).get("methods", {})
        except Exception:  # noqa
            existing = {}
    gt_dir = GT_ROOT / system_id
    lig_chain = meta["ligand_chain"]
    gt_mol = Chem.MolFromMolFile(str(gt_dir / "ligand_files" / f"{lig_chain}.sdf"), removeHs=False)
    gt_rec = st.Receptor.from_file(gt_dir / "receptor.cif")
    gt_xyz = Chem.RemoveHs(gt_mol).GetConformer().GetPositions()
    sdir = OUT / "structures" / system_id
    prev_gt = {}
    if reuse and prev_path.exists():
        try:
            prev_gt = json.loads(prev_path.read_text()).get("gt", {})
        except Exception:  # noqa
            prev_gt = {}
    gt_files = {"receptor": "gt_receptor.pdb", "ligand": "gt_ligand.sdf", "traj": "gt_traj.sdf",
                "pocket_traj": "gt_pocket_traj.sdf", "pocket_traj_pdb": "gt_pocket_traj.pdb"}
    if prev_gt.get("pocket_minimisation", {}).get("ok") and (sdir / "gt_pocket_traj.pdb.gz").exists():
        gt_block = prev_gt  # ground-truth analysis is deterministic; reuse it
    else:
        gz_write(sdir / "gt_receptor.pdb.gz", gt_rec.to_pdb())
        gz_write(sdir / "gt_ligand.sdf.gz", st.mol_to_sdf(Chem.RemoveHs(gt_mol), "ground truth"))
        gt_min = None
        try:
            gt_min, frames, mh = st.minimise_with_trajectory(Chem.RemoveHs(gt_mol))
            gz_write(sdir / "gt_traj.sdf.gz", st.frames_to_sdf(mh, frames, gt_min["energies"]))
        except Exception:  # noqa
            gt_min = None
        gt_pocket = _pocket_minimisation(Chem.RemoveHs(gt_mol), gt_rec.to_pdb(), sdir, "gt")
        gt_prot_atoms = list(gt_rec.heavy_atoms())
        pocket = sorted({f"{ch}:{resn}{resi}" for el, xyz, ch, resn, resi, an in gt_prot_atoms
                         if np.min(np.linalg.norm(gt_xyz - xyz, axis=1)) <= st.CONTACT_CUTOFF})
        gt_block = {"minimisation": gt_min, "pocket_minimisation": gt_pocket, "files": gt_files, "pocket_residues": pocket}

    detail = {**meta, "gt": gt_block,
              "methods": {k: v for k, v in existing.items() if k not in {m["id"] for m in method_list}}}

    from posebusters import PoseBusters
    pb = PoseBusters(config="redock")
    work = WORK / system_id
    work.mkdir(parents=True, exist_ok=True)

    for m in method_list:
        mid = m["id"]
        row = method_rows.get(mid)
        res = {"ok": False}
        prev = existing.get(mid)
        if prev and prev.get("ok") and (prev.get("pocket_minimisation") or {}).get("ok") and (sdir / f"{mid}_ligand.sdf.gz").exists():
            detail["methods"][mid] = prev
            continue
        try:
            if row is None:
                raise RuntimeError("no prediction rows")
            seed, sample = row.get("seed", ""), row.get("sample", "")
            if m["kind"] == "docking":
                # docked into the crystal receptor: pose is already in the ground-truth frame
                path = DOCKED / mid / system_id / "pose.sdf"
                if not path.exists():
                    raise RuntimeError("docked pose not found")
                lig = Chem.MolFromMolFile(str(path), removeHs=False)
                if lig is None:
                    raise RuntimeError("docked pose SDF unreadable")
                lig = Chem.RemoveHs(lig)
                tmpl = Chem.RemoveHs(gt_mol)
                if lig.GetNumAtoms() != tmpl.GetNumAtoms():
                    raise RuntimeError(f"atom count mismatch docked={lig.GetNumAtoms()} crystal={tmpl.GetNumAtoms()}")
                lig = AllChem.AssignBondOrdersFromTemplate(tmpl, lig)
                Chem.AssignStereochemistryFrom3D(lig)
                pred_rec = gt_rec
                sup = {"pocket_ca_rmsd": 0.0, "global_ca_rmsd_pocket_aligned": 0.0, "n_pocket_ca": 0, "chain_map": {}}
            else:
                path = find_model_file(m, system_id, seed, sample)
                if path is None:
                    raise RuntimeError(f"model file not found for seed={seed} sample={sample}")
                pred_rec = st.Receptor.from_file(path)
                sup = st.superpose_on_pocket(pred_rec, gt_rec, gt_xyz)
                chain = ligand_chain_in_model(row, path, gt_mol, meta["n_protein_chains"])
                lig = st.ligand_from_structure(path, chain, gt_mol)
                st.transform_mol(lig, sup["R"], sup["t"])
            rmsd = st.symmetry_rmsd(lig, gt_mol)
            prot_atoms = list(pred_rec.heavy_atoms())
            diag = st.ligand_diagnostics(lig, gt_mol, prot_atoms)
            pdb_text = pred_rec.to_pdb()
            rec_path = work / f"{mid}_receptor.pdb"
            rec_path.write_text(pdb_text)
            pb_df = pb.bust([lig], gt_mol, str(rec_path), full_report=False)
            pb_res = {}
            for k, v in pb_df.iloc[0].items():
                if k in PB_EXCLUDE:
                    continue
                pb_res[k] = None if pd.isna(v) else bool(v)
            pb_pass = all(v for v in pb_res.values() if v is not None)
            minim, frames, mh = None, None, None
            try:
                minim, frames, mh = st.minimise_with_trajectory(lig)
            except Exception as e:  # noqa
                minim = None
            gz_write(sdir / f"{mid}_receptor.pdb.gz", pdb_text)
            gz_write(sdir / f"{mid}_ligand.sdf.gz", st.mol_to_sdf(Chem.RemoveHs(lig), f"{mid} prediction"))
            if minim is not None:
                gz_write(sdir / f"{mid}_traj.sdf.gz", st.frames_to_sdf(mh, frames, minim["energies"]))
            pocket_min = _pocket_minimisation(lig, pdb_text, sdir, mid)
            res = {
                "ok": True, "seed": str(seed), "sample": str(sample),
                "ranking_score": _f(row.get("ranking_score")),
                "rmsd": round(rmsd, 3), "rmsd_ref": _f(row.get("rmsd")), "lddt_pli": _f(row.get("lddt_pli")),
                "lddt_lp": _f(row.get("lddt_lp")), "bb_rmsd": _f(row.get("bb_rmsd")),
                "superposition": {k: (round(v, 3) if isinstance(v, float) else v) for k, v in sup.items() if k not in ("R", "t")},
                "posebusters": pb_res, "pb_pass": pb_pass, "diagnostics": diag, "minimisation": minim,
                "pocket_minimisation": pocket_min,
                "files": {"receptor": f"{mid}_receptor.pdb", "ligand": f"{mid}_ligand.sdf", "traj": f"{mid}_traj.sdf",
                          "pocket_traj": f"{mid}_pocket_traj.sdf", "pocket_traj_pdb": f"{mid}_pocket_traj.pdb"},
                "model_file": str(path.relative_to(RAW)),
                "vina_score": _f(row.get("vina_score")), "pose_scores": row.get("pose_scores"),
            }
        except Exception as e:  # noqa
            res = {"ok": False, "error": f"{type(e).__name__}: {e}"}
        detail["methods"][mid] = res

    (OUT / "systems").mkdir(parents=True, exist_ok=True)
    (OUT / "systems" / f"{system_id}.json").write_text(json.dumps(detail, separators=(",", ":")))
    detail["_elapsed"] = round(time.time() - t0, 1)
    return detail


def _pocket_minimisation(lig, receptor_pdb_text, sdir, prefix):
    """Ligand + pocket restrained minimisation; writes the two trajectory files and returns the summary."""
    try:
        summ, lig_frames, pk_frames, pocket_info, ligH = st.minimise_in_pocket(lig, receptor_pdb_text)
    except Exception as e:  # noqa
        return {"ok": False, "error": f"{type(e).__name__}: {e}"}
    if summ.get("ok"):
        gz_write(sdir / f"{prefix}_pocket_traj.sdf.gz", st.frames_to_sdf(ligH, lig_frames, summ["energies"]))
        gz_write(sdir / f"{prefix}_pocket_traj.pdb.gz", st.pocket_frames_to_pdb(pocket_info, pk_frames))
    return summ


def _f(x):
    try:
        if x is None or pd.isna(x):
            return None
        return round(float(x), 4)
    except Exception:  # noqa
        return None


# --------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--systems", nargs="*")
    ap.add_argument("--methods", nargs="*")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--index-only", action="store_true", help="rebuild index.json from existing systems/*.json")
    ap.add_argument("--force", action="store_true", help="recompute methods even if a previous ok result exists")
    args = ap.parse_args()

    subset = json.load(open(HERE / "subset.json"))
    systems = [x for a in (args.systems or []) for x in a.split()] or subset["systems"]
    methods_arg = [x for a in (args.methods or []) for x in a.split()]
    method_list = [m for m in METHODS if not methods_arg or m["id"] in methods_arg]
    ann, preds, sim = load_inputs()
    inputs = json.load(open(RAW / "inputs.json"))

    metas, rows_by_system = {}, {}
    for sid in systems:
        a = ann[ann["system_id"] == sid].iloc[0]
        metas[sid] = system_meta(sid, ann, inputs, sim)
        rows_by_system[sid] = {m["id"]: (None if m["id"] not in preds else top_ranked(preds[m["id"]], sid, a["ligand_instance_chain"]))
                               for m in method_list}
        rows_by_system[sid] = {k: (None if v is None else v.to_dict()) for k, v in rows_by_system[sid].items()}
        for m in method_list:
            if m["kind"] == "docking":
                rp = DOCKED / m["id"] / sid / "result.json"
                r = json.loads(rp.read_text()) if rp.exists() else None
                rows_by_system[sid][m["id"]] = ({"seed": m["id"], "sample": "0", "ranking_score": r["vina_score"], "vina_score": r["vina_score"],
                                                "pose_scores": r["pose_scores"]} if r and r.get("ok") else None)

    OUT.mkdir(parents=True, exist_ok=True)
    if not args.index_only:
        with ProcessPoolExecutor(args.workers) as ex:
            futs = {ex.submit(process_system, sid, metas[sid], rows_by_system[sid], method_list, not args.force): sid for sid in systems}
            for fut in as_completed(futs):
                sid = futs[fut]
                try:
                    d = fut.result()
                    oks = [k for k, v in d["methods"].items() if v["ok"]]
                    fails = {k: v["error"] for k, v in d["methods"].items() if not v["ok"]}
                    print(f"[{d['_elapsed']:6.1f}s] {sid}: ok={oks} fail={fails}", flush=True)
                except Exception:
                    print(f"FAILED {sid}\n{traceback.format_exc()}", flush=True)

    # ---- index.json from whatever system files exist
    results, sys_rows = [], []
    for sid in subset["systems"]:
        p = OUT / "systems" / f"{sid}.json"
        if not p.exists():
            continue
        d = json.loads(p.read_text())
        fresh = system_meta(sid, ann, inputs, sim)  # keep per-system meta in sync with the index
        if any(d.get(k) != v for k, v in fresh.items()):
            d.update(fresh)
            p.write_text(json.dumps(d, separators=(",", ":")))
        sys_rows.append({k: d[k] for k in META_KEYS})
        for mid, r in d["methods"].items():
            if not r["ok"]:
                results.append({"system_id": sid, "method": mid, "ok": False, "error": r.get("error")})
                continue
            mn = r.get("minimisation") or {}
            pm = r.get("pocket_minimisation") or {}
            results.append({
                "system_id": sid, "method": mid, "ok": True, "seed": r["seed"], "sample": r["sample"],
                "ranking_score": r["ranking_score"], "rmsd": r["rmsd"], "rmsd_ref": r["rmsd_ref"], "lddt_pli": r["lddt_pli"],
                "pb_pass": r["pb_pass"], "pb_fail_count": sum(1 for v in r["posebusters"].values() if v is False),
                "strain_local": mn.get("strain_local"), "strain_global": mn.get("strain_global"), "rmsd_drift": mn.get("rmsd_drift"),
                "protein_clashes": r["diagnostics"]["summary"]["protein_clashes"], "flagged_atoms": len(r["diagnostics"]["flagged_atoms"]),
                "e_interaction_pose": pm.get("e_interaction_pose"), "e_interaction_min": pm.get("e_interaction_min"),
                "pocket_ligand_drift": pm.get("ligand_rmsd_drift"), "clashes_pose": pm.get("clashes_pose"), "clashes_min": pm.get("clashes_min"),
                "vina_score": r.get("vina_score"),
            })
    index = {
        "generated": time.strftime("%Y-%m-%d %H:%M"),
        "methods": [{k: m[k] for k in ("id", "name", "kind", "training_cutoff", "color", "note")} for m in METHODS
                    if any(r["method"] == m["id"] and r["ok"] for r in results)],
        "cutoffs": CUTOFFS, "default_cutoff": DEFAULT_CUTOFF, "buckets": subset["buckets"],
        "systems": sys_rows, "results": results,
    }
    (OUT / "index.json").write_text(json.dumps(index, separators=(",", ":")))
    print(f"index.json: {len(sys_rows)} systems, {sum(r['ok'] for r in results)} ok results, {sum(not r['ok'] for r in results)} failed")


if __name__ == "__main__":
    main()
