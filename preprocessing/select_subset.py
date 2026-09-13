"""Select a small, similarity-stratified demo subset of Runs N' Poses systems.

Writes subset.json with the chosen system_ids and metadata. Selection criteria keep the
demo dataset small (single proper ligand, <=2 protein chains, short sequences) while
covering the full range of similarity to the training set.
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path

RAW = Path(__file__).parent / "raw"
METHODS_REQUIRED = ["af3", "af3_no_template", "boltz", "boltz1x", "chai", "protenix"]
METHODS_OPTIONAL = ["boltz2", "rfaa"]
BUCKETS = [(0, 20), (20, 40), (40, 60), (60, 80), (80, 101)]
PER_BUCKET = 30
MAX_SEQ_LEN = 450
SEED = 0

ann = pd.read_csv(RAW / "annotations.csv")
inputs = json.load(open(RAW / "inputs.json"))
seq_len = {k: sum(len(s) for s in v["sequences"].values()) for k, v in inputs.items()}

present = {}
for m in METHODS_REQUIRED + METHODS_OPTIONAL:
    df = pd.read_csv(RAW / "predictions" / "predictions" / f"{m}.csv", usecols=["target"])
    present[m] = set(df["target"].unique())

cand = ann[
    (ann["ligand_is_proper"] == True)
    & (ann["num_proper_ligand_chains"] == 1)
    & (ann["num_protein_chains"] <= 2)
    & (ann["ligand_num_heavy_atoms"].between(12, 50))
].copy()
cand["seq_len"] = cand["system_id"].map(seq_len)
cand = cand[cand["seq_len"] <= MAX_SEQ_LEN]
for m in METHODS_REQUIRED:
    cand = cand[cand["system_id"].isin(present[m])]
cand["has_boltz2"] = cand["system_id"].isin(present["boltz2"])
cand["has_rfaa"] = cand["system_id"].isin(present["rfaa"])
# NaN similarity = no training system with detectable pocket+ligand similarity -> treat as 0
cand["similarity"] = cand["sucos_shape_pocket_qcov"].fillna(0.0)
cand = cand.drop_duplicates("system_id")
print("candidates:", len(cand))

rng = np.random.default_rng(SEED)
chosen = []
used_pdb, used_cluster = set(), set()
# keep systems from an earlier, smaller selection so already-processed data is reused
prior = Path(__file__).parent / "subset_50.json"
prior_ids = set(json.load(open(prior))["systems"]) if prior.exists() else set()
for _, r in cand[cand["system_id"].isin(prior_ids)].iterrows():
    chosen.append(r); used_pdb.add(r["entry_pdb_id"]); used_cluster.add(r["cluster"])
for lo, hi in BUCKETS:
    b = cand[(cand["similarity"] >= lo) & (cand["similarity"] < hi)].copy()
    b["rand"] = rng.random(len(b))
    # prefer systems covered by every method, then random
    b = b.sort_values(["has_boltz2", "has_rfaa", "rand"], ascending=[False, False, True])
    n = int(b["system_id"].isin(prior_ids).sum())
    for _, r in b.iterrows():
        if r["system_id"] in prior_ids:
            continue
        if r["entry_pdb_id"] in used_pdb or r["cluster"] in used_cluster:
            continue
        chosen.append(r)
        used_pdb.add(r["entry_pdb_id"]); used_cluster.add(r["cluster"])
        n += 1
        if n >= PER_BUCKET:
            break
    print(f"bucket {lo}-{hi}: {n} chosen of {len(b)} candidates")

sel = pd.DataFrame(chosen)
out = {
    "systems": sel["system_id"].tolist(),
    "methods": METHODS_REQUIRED + METHODS_OPTIONAL,
    "buckets": BUCKETS,
}
json.dump(out, open(Path(__file__).parent / "subset.json", "w"), indent=1)
print(sel[["system_id", "similarity", "seq_len", "ligand_num_heavy_atoms", "has_boltz2", "has_rfaa"]].to_string())
