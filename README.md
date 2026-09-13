# Angstrom

An interactive explorer for how protein–ligand **co-folding models memorise rather than generalise**, built on the
[Runs N' Poses](https://github.com/plinder-org/runs-n-poses) benchmark. For every system in a demo subset it shows,
per method, whether the top-ranked pose is right (RMSD), whether it is physically plausible (PoseBusters, strain),
and *where* in the ligand things go wrong, all as a function of how similar the system is to the model's training set.

Two views:

| View | What it shows |
| --- | --- |
| **Overview** (`/`) | Success rate (ligand RMSD ≤ 2 Å) per method, with a slider on max similarity to the training set and a selector for the assumed training-cutoff date. Bar chart by similarity bucket, line chart of success vs. cutoff date, and a systems table linking to the 3D view and to RCSB. |
| **System** (`/system/<id>`) | 3Dmol.js scene of the ground-truth receptor + ligand with the selected method's pose superposed on the binding site. Toggles for the predicted receptor, pocket residues, violation markers and an MMFF minimisation trajectory. Side panel with accuracy, strain energies, PoseBusters checks and atom-level diagnostics that highlight offending atoms on hover. |

## Architecture

```
angstrom/
├── preprocessing/   heavy Python env: rdkit, posebusters, gemmi, pandas  (runs offline)
│   ├── select_subset.py               pick the demo systems -> subset.json
│   ├── stream_extract_predictions.py  stream the 40 GB prediction tarball, keep only the subset
│   ├── structure_tools.py             superposition, RMSD, diagnostics, minimisation
│   └── build_dataset.py               writes api/data/
├── api/             light Python env: FastAPI only (one Vercel function)
│   ├── index.py                       /api/index, /api/systems/{id}, /api/structures/{id}/{file}
│   └── data/                          precomputed JSON + gzipped PDB/SDF (committed, no database)
├── web/             React 19 + Vite + Tailwind v4 + Radix primitives + 3Dmol.js + Plotly (basic bundle)
└── vercel.json      static frontend + /api rewrite to the Python function
```

Two Python environments are deliberate: the preprocessing env pulls in RDKit and PoseBusters, the web function
only needs FastAPI so the serverless bundle stays small (data plus a few MB of code).

## Demo subset

`select_subset.py` picks 50 systems from the 2 600 in Runs N' Poses:

- one proper ligand (no ions/artifacts), ≤ 2 protein chains, ≤ 450 residues in total, 12–50 ligand heavy atoms;
- predictions exist for all core methods (AF3, AF3 no-templates, Boltz-1, Boltz-1x, Chai-1, Protenix); Boltz-2 and RFAA
  are included where available;
- 10 systems per SuCOS-pocket similarity bucket (0–20, 20–40, 40–60, 60–80, 80–100), at most one per PDB entry and per
  similarity cluster, preferring systems every method covers. Systems with no detectable training neighbour are
  treated as similarity 0.

Similarity is `sucos_shape_pocket_qcov` from `annotations.csv`: the SuCOS shape/colour overlap of the ligand with the
closest training-set ligand, scaled by pocket coverage (0–100).

## What is precomputed per (system, method)

For the top-ranked sample (highest `ranking_score` in the benchmark's `predictions/*.csv`):

| Quantity | How |
| --- | --- |
| Superposition | Chains matched by sequence alignment (gemmi); Kabsch fit on Cα atoms within 8 Å of the ligand. All predicted coordinates are stored in the ground-truth frame, so methods overlay in the viewer. |
| Ligand RMSD | Symmetry-corrected heavy-atom RMSD (`rdMolAlign.CalcRMS`) without re-alignment, after the pocket superposition. Bond orders are assigned to the predicted ligand from the crystal ligand template. |
| Benchmark metrics | BiSyRMSD (`rmsd`), lDDT-PLI, lDDT-LP, backbone RMSD and ranking score copied from the benchmark CSVs. |
| PoseBusters | `PoseBusters(config="redock")` on the predicted ligand vs. crystal ligand and superposed predicted receptor. `pb_pass` = all checks true. |
| Atom-level diagnostics | Bond lengths (ratio outside 0.85–1.15) and angles (> 20°) compared with the crystal pose of the same ligand; intra-ligand and protein–ligand clashes at 0.75 × vdW-sum (the PoseBusters threshold); carbon stereocentre mismatches; aromatic ring planarity; residues within 4.5 Å. |
| Minimisation | MMFF94s in RDKit. Hydrogens are relaxed with heavy atoms restrained (`e_pose`), then a free minimisation is recorded in frames of 15 iterations (≤ 30 frames, `e_local`), then 30 ETKDG conformers are minimised for the global reference (`e_global`). `strain_local = e_pose − e_local`, `strain_global = e_pose − e_global`, plus per-atom displacement and RMSD drift. The crystal ligand goes through the same protocol as a baseline. |

Output files under `api/data/structures/<system_id>/`: `gt_receptor.pdb.gz`, `gt_ligand.sdf.gz`, `gt_traj.sdf.gz`
and `<method>_receptor.pdb.gz`, `<method>_ligand.sdf.gz`, `<method>_traj.sdf.gz`. The API serves them with
`Content-Encoding: gzip` so the browser decompresses transparently.

### Similarity vs. training cutoff

The benchmark's headline similarity assumes a 30 Sep 2021 structural cutoff. `all_similarity_scores.parquet` holds every
query–target pair against the whole PDB, so for each subset system the pipeline recomputes the max
`sucos_shape_pocket_qcov` over targets released before each of nine cutoff dates (2018-01-01 … 2024-01-01, including
2021-09-30 and Boltz-2's 2023-06-01). The overview's cutoff selector and the "success vs. cutoff" plot use these
values; the recomputed 2021-09-30 and 2023-06-01 columns reproduce the benchmark's own annotations exactly.

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

.venv/bin/python select_subset.py                   # -> subset.json (deterministic, seed 0)
.venv/bin/python stream_extract_predictions.py      # streams the 40 GB tarball, keeps subset files only (~1 h)
.venv/bin/python build_dataset.py --workers 4       # -> ../api/data  (add --index-only to just rebuild index.json)
```

`stream_extract_predictions.py` never stores the archive: it prefetches HTTP ranges in parallel and feeds them to a
streaming tar reader, writing only members whose path contains a subset `system_id`. `raw/` is git-ignored.

`build_dataset.py` reads `raw/subset_similarity_scores.parquet`, a filter of `all_similarity_scores.parquet` on the
subset's systems, if present (without it the cutoff selector falls back to the benchmark similarity). Create it with:

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
  [gemmi](https://gemmi.readthedocs.io/), [3Dmol.js](https://3dmol.csb.pitt.edu/), [Plotly](https://plotly.com/javascript/).
- Ground-truth links point to `https://www.rcsb.org/structure/<PDB id>`.

## Known limitations

- 50 systems is enough to see the trend, not to reproduce the paper's numbers; bucket sizes are 10 systems.
- Only the top-ranked sample per method is analysed; the other 24 samples' benchmark metrics are not shown.
- Strain is ligand-only MMFF94s in vacuum (no pocket during minimisation); clashes with the receptor are reported by the
  static diagnostics and PoseBusters instead. Absolute strain values depend on protonation and are best read relative to
  the crystal-ligand baseline shown alongside.
- Bond/angle deviations are measured against the crystal pose, which is itself a refined model, not an ideal geometry.
- Stereo mismatches are only evaluated on carbon centres (matching PoseBusters); phosphorus/sulfur pseudo-centres are ignored.
- The overview aggregates RMSD only; PoseBusters validity and strain are available per system and in `index.json`
  (`pb_pass`, `strain_local`, …) for further plots.
- Predicted receptors are stored as protein heavy atoms only; waters, ions and cofactors from the models are dropped.
