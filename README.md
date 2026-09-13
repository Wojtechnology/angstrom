# Angstrom

**Which co-folding models learn physical priors?** An interactive explorer built on the
[Runs N' Poses](https://github.com/plinder-org/runs-n-poses) benchmark. A model that has learned physics should keep
producing physically plausible poses — no clashes, low ligand strain, a favourable interaction energy with the pocket —
even when a system looks nothing like its training set. A model that has mostly learned to recall should degrade on
both accuracy *and* plausibility as similarity to the training set drops. For every system in a demo subset, the app
shows per method whether the top-ranked pose is right (RMSD), whether it is plausible (PoseBusters, clashes, strain,
interaction energy before and after a restrained pocket relaxation), and *where* in the ligand things go wrong. A rigid
holo-redocking baseline (AutoDock Vina) marks what a physics-only method achieves on the same systems — an easier
problem, since the pocket is already in its bound conformation, and shown with a badge in the UI.

Two views:

| View | What it shows |
| --- | --- |
| **Overview** (`/`) | Method chips (selection persisted in `?m=`) and a **Metric** select: RMSD ≤ 2 Å success rate, PoseBusters-valid rate, clash-free-at-pose rate, median relaxation ΔE (interaction energy pose − min), median local ligand strain. Bar chart of the metric per SuCOS-pocket similarity bucket; scatter of a selectable y (RMSD on a log axis, relaxation ΔE, local / global strain, pocket clashes, PoseBusters failures; extreme energies are clamped and drawn as triangles) vs. similarity with click-through to the system page; compact table (PDB id, CCD code linking to the RCSB ligand page, similarity, RMSD per method with an amber dot for PoseBusters failures and a red dot for pocket clashes). Energies are labelled lower-is-better. |
| **System** (`/system/<id>?m=<method>`) | 3Dmol.js scene of the crystal complex with the selected method's pose superposed on the binding site. 2×2 checkboxes (ground truth / predicted × ligand / receptor), pocket residues, violation markers, and a trajectory select (off / ligand only / ligand + pocket) with playback. A **Physical violations** panel under the viewer lists bond, angle, clash, stereo and ring problems with two-way hover between entries and viewer atoms. Side cards: accuracy, pocket relaxation (interaction energy pose → min, clashes before / after, ligand drift, pocket RMSD, ground-truth baseline), ligand-only strain, PoseBusters. |

Similarity throughout is the benchmark's `sucos_shape_pocket_qcov` (training cutoff 30 Sep 2021): SuCOS shape/colour
overlap of the ligand with the closest training-set ligand, scaled by pocket coverage, 0–100. Systems with no detectable
training neighbour count as 0.

## Architecture

```
angstrom/
├── preprocessing/   heavy Python env: rdkit, posebusters, gemmi, pandas  (runs offline)
│   ├── select_subset.py               pick the demo systems -> subset.json (keeps subset_50.json's systems)
│   ├── stream_extract_predictions.py  stream the 40 GB prediction tarball, keep only af3/boltz/boltz1x × subset
│   ├── dock_vina.py                   Vina rigid holo-redocking baseline (separate micromamba env)
│   ├── structure_tools.py             superposition, RMSD, diagnostics, ligand + pocket minimisation
│   └── build_dataset.py               writes api/data/ (incremental; --force to recompute)
├── api/             light Python env: FastAPI only (one Vercel function)
│   ├── index.py                       /api/index, /api/systems/{id}, /api/structures/{id}/{file}
│   └── data/                          precomputed JSON + gzipped PDB/SDF (committed, no database)
├── web/             React 19 + Vite + Tailwind v4 + Radix primitives + 3Dmol.js + Plotly (basic bundle)
└── vercel.json      static frontend + /api rewrite to the Python function
```

Two Python environments are deliberate: preprocessing needs RDKit, PoseBusters and gemmi; the web function needs only
FastAPI, so the serverless bundle is the data plus a few MB of code. Docking runs in a third, micromamba-managed env
(`dock`: vina 1.2, meeko 0.8, rdkit) because Vina's Python bindings are not pip-installable on every platform.

## Methods

| Method | Kind | Training cutoff | Source |
| --- | --- | --- | --- |
| AF3 | co-folding | 2021-09-30 | Runs N' Poses `prediction_files/af3` |
| Boltz-1 | co-folding | 2021-09-30 | `prediction_files/boltz` |
| Boltz-1x | co-folding | 2021-09-30 | `prediction_files/boltz1x` |
| Vina (rigid holo redocking) | docking baseline | – | `dock_vina.py`, AutoDock Vina 1.2 + Meeko |

For the co-folding methods the top-ranked sample (highest `ranking_score` in the benchmark's `predictions/*.csv`) is
used. The other Runs N' Poses methods (AF3 without templates, Boltz-2, Chai-1, Protenix, RFAA) are listed in
`build_dataset.py` as `EXTRA_METHODS` but are not extracted or processed.

**Vina baseline.** The crystal (holo) receptor is kept rigid (protein heavy atoms + polar hydrogens, AD4 typing;
waters, ions and cofactors dropped). The ligand starts from a fresh random ETKDG conformer, never from the crystal
pose, and is docked into a box around the crystal ligand (ligand extent + 8 Å per side, min 20 Å) with
exhaustiveness 16; the top-scoring pose is taken and its bond orders assigned from the crystal ligand. It is already in
the ground-truth frame, so it goes through the same scoring, PoseBusters and relaxation as the co-folding poses.

## Demo subset

`select_subset.py` picks 150 systems from the 2 600 in Runs N' Poses, 30 per SuCOS-pocket similarity bucket
(0–20, 20–40, 40–60, 60–80, 80–100):

- one proper ligand (no ions/artifacts), ≤ 2 protein chains, ≤ 450 residues in total, 12–50 ligand heavy atoms;
- benchmark predictions exist for AF3, AF3 no-templates, Boltz-1, Boltz-1x, Chai-1 and Protenix;
- at most one system per PDB entry and per similarity cluster, preferring systems every method covers, otherwise
  random with seed 0. The original 50-system selection (`subset_50.json`) is kept verbatim so already-processed data is
  reused.

## What is precomputed per (system, method)

| Quantity | How |
| --- | --- |
| Superposition | Chains matched by sequence alignment (gemmi); Kabsch fit on Cα atoms within 8 Å of the crystal ligand. All predicted coordinates are stored in the ground-truth frame so methods overlay in the viewer. Chain ids are shortened to single letters at load time so JSON labels and PDB files agree. |
| Ligand RMSD | Symmetry-corrected heavy-atom RMSD (`rdMolAlign.CalcRMS`) without re-alignment. Bond orders for the predicted ligand come from the crystal ligand template. |
| Benchmark metrics | BiSyRMSD (`rmsd`), lDDT-PLI, lDDT-LP, backbone RMSD and ranking score copied from the benchmark CSVs (co-folding only; Vina reports its affinity score). |
| PoseBusters | `PoseBusters(config="redock")` on the pose vs. crystal ligand and the (superposed) receptor. `pb_pass` = every check true. |
| Atom-level diagnostics | Bond lengths (ratio outside 0.85–1.15) and angles (> 20°) vs. the crystal pose of the same ligand; intra-ligand and protein–ligand clashes at 0.75 × vdW-sum (the PoseBusters threshold); carbon stereocentre mismatches; aromatic ring planarity; residues within 4.5 Å. |
| Ligand-only strain | MMFF94s in vacuum. Hydrogens relaxed with heavy atoms restrained (`e_pose`); free minimisation recorded in frames of 15 iterations (≤ 30 frames, `e_local`); 30 ETKDG conformers minimised for the global reference (`e_global`). `strain_local = e_pose − e_local`, `strain_global = e_pose − e_global`, plus per-atom displacement and RMSD drift. |
| Pocket relaxation | Restrained minimisation of the ligand inside its pocket, see below. Reports complex and interaction energies at the pose and at the minimum, clash counts before/after, ligand drift, pocket heavy-atom RMSD and a trajectory. |

The crystal complex goes through the same strain and pocket-relaxation protocol and is shown as the baseline.

### Pocket relaxation protocol (`structure_tools.minimise_in_pocket`)

1. Whole residues with any heavy atom within 8 Å of the ligand form the pocket; the receptor PDB is parsed with
   proximity bonding and spurious bonds between clashing atoms are removed until it sanitises.
2. Hydrogens are added with coordinates to pocket and ligand (no source structure carries any; protonation is
   RDKit's neutral default). The complex is typed with MMFF94s using a distance-dependent dielectric, ε = 4,
   9 Å non-bonded threshold; UFF is the fallback when MMFF cannot type the complex.
3. Stage 0: hydrogen-only relaxation with every heavy atom fixed → `e_complex_pose`,
   `e_interaction_pose = E(complex) − E(pocket) − E(ligand)`.
4. Stage 1: backbone atoms (N, CA, C, O) and pocket atoms further than 6 Å from the ligand are fixed; remaining
   side-chain heavy atoms get a flat-bottom position restraint (0.3 Å, 5 kcal/mol/Å²); ligand and hydrogens are free.
   Minimised in frames of 40 iterations, up to 16 frames, then to convergence → `e_complex_min`, `e_interaction_min`.
5. Clashes are ligand–pocket heavy-atom pairs closer than 0.75 × vdW-sum, counted at the pose and at the minimum.

Output files under `api/data/structures/<system_id>/` per prefix (`gt`, `af3`, `boltz`, `boltz1x`, `vina`):
`<prefix>_receptor.pdb.gz`, `<prefix>_ligand.sdf.gz`, `<prefix>_traj.sdf.gz` (ligand-only frames),
`<prefix>_pocket_traj.sdf.gz` and `<prefix>_pocket_traj.pdb.gz` (ligand and movable pocket residues per frame). The API
serves them with `Content-Encoding: gzip`. Receptors are protein heavy atoms only; waters, ions and cofactors are
stripped everywhere.

`index.json` also carries, per system, similarity recomputed from `all_similarity_scores.parquet` for nine cutoff
dates (SuCOS-pocket, Morgan Tanimoto and pocket coverage, excluding the system's own PDB entry). The UI currently uses
only the benchmark's 2021-09-30 SuCOS-pocket value.

## Local development

Node ≥ 20 and [uv](https://docs.astral.sh/uv/).

```bash
# API (light env)
cd api && uv venv --python 3.12 .venv && uv pip install --python .venv/bin/python -r requirements.txt uvicorn
.venv/bin/uvicorn index:app --port 8000 --reload

# Frontend (proxies /api to :8000, see web/vite.config.ts)
cd web && npm ci && npm run dev        # http://localhost:5173
```

## Rebuilding the data

```bash
cd preprocessing
uv venv --python 3.12 .venv && uv pip install --python .venv/bin/python -e .
mkdir -p raw && cd raw
Z=https://zenodo.org/api/records/18366081/files
curl -L $Z/annotations.csv/content -o annotations.csv
curl -L $Z/inputs.json/content -o inputs.json
curl -L $Z/predictions.tar.gz/content -o predictions.tar.gz && mkdir predictions && tar xzf predictions.tar.gz -C predictions
curl -L $Z/ground_truth.tar.gz/content -o ground_truth.tar.gz && tar xzf ground_truth.tar.gz     # 414 MB
curl -L $Z/all_similarity_scores.parquet/content -o all_similarity_scores.parquet               # 376 MB
cd ..

.venv/bin/python select_subset.py                   # -> subset.json (150 systems, deterministic)
.venv/bin/python stream_extract_predictions.py      # streams the 40 GB tarball, keeps af3/boltz/boltz1x × subset
$HOME/.local/micromamba/envs/dock/bin/python dock_vina.py --workers 8     # -> raw/docked/vina/<id>/pose.sdf
.venv/bin/python build_dataset.py --workers 4       # -> ../api/data  (incremental: keeps ok results, --force recomputes,
                                                    #    --index-only rebuilds index.json, --systems/--methods to restrict)
```

`stream_extract_predictions.py` never stores the archive: it prefetches 32 MB HTTP ranges with 12 connections and feeds
them to a streaming tar reader, writing only members under the `af3`, `boltz` and `boltz1x` method directories whose
path contains a subset `system_id`. The tarball is grouped by method, so it stops as soon as those three directories
have been passed. gzip cannot be resumed mid-stream, so Zenodo outages (504s, short reads) are retried indefinitely
with capped backoff rather than abandoned. `raw/` is git-ignored.

The docking env:

```bash
micromamba create -n dock -c conda-forge python=3.11 vina=1.2 meeko=0.8 rdkit numpy pandas   # -> ~/.local/micromamba/envs/dock
```

`build_dataset.py` reads `raw/subset_similarity_scores.parquet`, a filter of `all_similarity_scores.parquet` on the
subset's systems, if present. Create it with:

```bash
.venv/bin/python -c "
import json, pandas as pd
ids = json.load(open('subset.json'))['systems']
cols = ['group_key','query_system','target_system','target_release_date','sucos_shape_pocket_qcov','pocket_qcov','pli_qcov','morgan_tanimoto']
df = pd.read_parquet('raw/all_similarity_scores.parquet', columns=cols)
df[df.query_system.isin(ids)].to_parquet('raw/subset_similarity_scores.parquet')"
```

## Deploying to Vercel

Import the repository with the **root directory left at the repo root**; `vercel.json` does the rest:

- `framework: "vite"` pins the frontend preset so Vercel does not switch to its FastAPI preset (which would take over
  every route). `buildCommand` runs `npm ci && npm run build` inside `web/`; `outputDirectory` is `web/dist`.
- No `installCommand` override. The Python builder disables its own `pip`/`uv` install when a custom install command is
  set, which would leave the function without FastAPI. The root `package.json` exists only so the default install step
  has something harmless to run.
- `api/index.py` is a file-based Python function. Its dependencies come from `api/requirements.txt` (the builder looks in
  the entrypoint directory first) and `api/.python-version` pins 3.12. The builder bundles the repo checkout minus
  `excludeFiles`; `web/**` and `preprocessing/**` are excluded so only `api/` (code + data) ships. The Python bundle
  limit is 500 MB uncompressed; the demo data is tens of MB.
- Rewrites send `/api/*` to the function (FastAPI routes are declared with the `/api` prefix) and every other path
  that is not a static file to `index.html` for client-side routing.

Deploy with `vercel` (CLI ≥ 48) or via the Git integration. Check `/api/health` after the first deploy.

## Data and attribution

- Runs N' Poses dataset: Zenodo [10.5281/zenodo.14794785](https://doi.org/10.5281/zenodo.14794785)
  (this build used record version 18366081). Ground truth is derived from the PDB; predictions are the outputs of the
  respective models and are covered by `OUTPUT_TERMS_OF_USE.md` in the Zenodo record (in particular the AlphaFold 3
  output terms). The demo subset in `api/data` redistributes a small fraction of those files for visualisation.
- Preprint: *Have protein-ligand co-folding methods moved beyond memorisation?* (PLINDER / Runs N' Poses team),
  bioRxiv 2025, [10.1101/2025.02.03.636309](https://doi.org/10.1101/2025.02.03.636309).
- Tools: [RDKit](https://www.rdkit.org/), [PoseBusters](https://github.com/maabuu/posebusters),
  [gemmi](https://gemmi.readthedocs.io/), [AutoDock Vina](https://vina.scripps.edu/), [Meeko](https://github.com/forlilab/Meeko),
  [3Dmol.js](https://3dmol.csb.pitt.edu/), [Plotly](https://plotly.com/javascript/).
- Ground-truth links point to `https://www.rcsb.org/structure/<PDB id>` and ligands to `https://www.rcsb.org/ligand/<CCD>`.

## Known limitations

- 150 systems (30 per bucket) show the trend; they do not reproduce the paper's numbers.
- Only the top-ranked sample per co-folding method and the top Vina pose are analysed.
- The Vina baseline redocks into the crystal pocket: it cannot fail on receptor conformation, so its numbers are an
  upper bound for a physics-only method, not a like-for-like competitor.
- Force-field energies are MMFF94s with a distance-dependent dielectric, neutral protonation, no waters or ions, and
  pocket residues restrained; read them relative to the crystal-complex baseline shown alongside, not as absolute
  binding energies. Interaction energy > 0 means the pose is repulsive with its own pocket.
- Bond/angle deviations are measured against the crystal pose, itself a refined model, not an ideal geometry.
- Stereo mismatches are evaluated on carbon centres only (matching PoseBusters).
- Predicted receptors are stored as protein heavy atoms; ions and cofactors present in the models are not scored.
