// Types mirror the JSON written by preprocessing/build_dataset.py.

export interface MethodInfo {
  id: string
  name: string
  training_cutoff: string
  color: string
}

export interface SystemSummary {
  system_id: string
  pdb_id: string
  ligand_smiles: string
  ccd: string
  ligand_chain: string
  n_heavy: number
  seq_len: number
  n_protein_chains: number
  release_date: string
  cluster: string
  closest_training: { system_id: string | null; pdb_id: string | null; release_date: string | null }
  similarity: number
  similarity_by_cutoff: Record<string, number>
}

export interface ResultRow {
  system_id: string
  method: string
  ok: boolean
  seed?: string | number
  sample?: string | number
  ranking_score?: number | null
  rmsd?: number | null
  rmsd_ref?: number | null
  lddt_pli?: number | null
  pb_pass?: boolean | null
  pb_fail_count?: number | null
  strain_local?: number | null
  strain_global?: number | null
  rmsd_drift?: number | null
  protein_clashes?: number | null
  flagged_atoms?: number | null
  e_interaction_pose?: number | null
  e_interaction_min?: number | null
  pocket_ligand_drift?: number | null
  clashes_pose?: number | null
  clashes_min?: number | null
  error?: string
}

export interface IndexData {
  generated: string
  methods: MethodInfo[]
  cutoffs: string[]
  default_cutoff: string
  buckets: [number, number][]
  systems: SystemSummary[]
  results: ResultRow[]
}

export interface Diagnostics {
  bonds: { atoms: number[]; pred: number; ref: number; ratio: number; flag: boolean }[]
  angles: { atoms: number[]; pred: number; ref: number; dev: number; flag: boolean }[]
  n_angles: number
  intra_clashes: { atoms: number[]; dist: number; limit: number }[]
  protein_clashes: { atom: number; protein: { chain: string; resname: string; resnum: number; atom: string }; dist: number; limit: number }[]
  stereo: { atom: number; pred: string | null; ref: string | null; flag: boolean }[]
  rings: { atoms: number[]; max_dev: number; flag: boolean }[]
  contacts: { residue: string; min_dist: number }[]
  flagged_atoms: number[]
  summary: Record<string, number>
}

export interface Minimisation {
  e_pose: number
  e_local: number
  e_global: number
  strain_local: number
  strain_global: number
  rmsd_drift: number
  max_atom_displacement: number
  atom_displacement: number[]
  n_frames: number
  converged: boolean
  energies: number[]
  frame_rmsd: number[]
}

export interface PocketMinimisation {
  ok: boolean
  error?: string
  n_pocket_residues: number
  n_pocket_atoms: number
  pocket_residues: string[]
  e_complex_pose: number
  e_complex_min: number
  e_interaction_pose: number
  e_interaction_min: number
  ligand_rmsd_drift: number
  ligand_atom_displacement: number[]
  pocket_heavy_rmsd: number
  max_pocket_atom_displacement: number
  clashes_pose: number
  clashes_min: number
  n_frames: number
  energies: number[]
  frame_ligand_rmsd: number[]
  converged: boolean
}

export interface StructureFiles { receptor: string; ligand: string; traj: string; pocket_traj?: string; pocket_traj_pdb?: string }

export interface MethodDetail {
  ok: boolean
  error?: string
  seed?: string | number
  sample?: string | number
  ranking_score?: number | null
  rmsd?: number
  rmsd_ref?: number | null
  lddt_pli?: number | null
  lddt_lp?: number | null
  bb_rmsd?: number | null
  superposition?: { pocket_ca_rmsd: number; global_ca_rmsd_pocket_aligned: number; n_pocket_ca: number; chain_map: Record<string, string> }
  posebusters?: Record<string, boolean | null>
  pb_pass?: boolean
  diagnostics?: Diagnostics
  minimisation?: Minimisation
  pocket_minimisation?: PocketMinimisation | null
  files?: StructureFiles
}

export interface SystemDetail extends SystemSummary {
  gt: { minimisation: Minimisation | null; pocket_minimisation?: PocketMinimisation | null; files: StructureFiles; pocket_residues: string[] }
  methods: Record<string, MethodDetail>
}

const cache = new Map<string, Promise<unknown>>()

function cached<T>(key: string, fn: () => Promise<T>): Promise<T> {
  if (!cache.has(key)) cache.set(key, fn())
  return cache.get(key) as Promise<T>
}

export function fetchIndex(): Promise<IndexData> {
  return cached('index', async () => {
    const r = await fetch('/api/index')
    if (!r.ok) throw new Error(`index: ${r.status}`)
    return r.json()
  })
}

export function fetchSystem(id: string): Promise<SystemDetail> {
  return cached(`system:${id}`, async () => {
    const r = await fetch(`/api/systems/${encodeURIComponent(id)}`)
    if (!r.ok) throw new Error(`system ${id}: ${r.status}`)
    return r.json()
  })
}

export function fetchStructure(id: string, name: string): Promise<string> {
  return cached(`structure:${id}/${name}`, async () => {
    const r = await fetch(`/api/structures/${encodeURIComponent(id)}/${encodeURIComponent(name)}`)
    if (!r.ok) throw new Error(`structure ${name}: ${r.status}`)
    return r.text()
  })
}

export const rcsbUrl = (pdbId: string) => `https://www.rcsb.org/structure/${pdbId.toUpperCase()}`

export const PB_CHECK_LABELS: Record<string, string> = {
  mol_pred_loaded: 'Prediction loads',
  mol_true_loaded: 'Ground truth loads',
  mol_cond_loaded: 'Receptor loads',
  sanitization: 'Sanitises',
  inchi_convertible: 'InChI convertible',
  all_atoms_connected: 'All atoms connected',
  no_radicals: 'No radicals',
  molecular_formula: 'Formula matches',
  molecular_bonds: 'Bonds match',
  double_bond_stereochemistry: 'Double-bond stereo',
  tetrahedral_chirality: 'Tetrahedral chirality',
  bond_lengths: 'Bond lengths',
  bond_angles: 'Bond angles',
  internal_steric_clash: 'No internal clash',
  aromatic_ring_flatness: 'Aromatic rings flat',
  'non-aromatic_ring_non-flatness': 'Aliphatic rings puckered',
  double_bond_flatness: 'Double bonds flat',
  internal_energy: 'Internal energy',
  'protein-ligand_maximum_distance': 'Ligand near protein',
  minimum_distance_to_protein: 'No protein clash',
  minimum_distance_to_organic_cofactors: 'No cofactor clash',
  minimum_distance_to_inorganic_cofactors: 'No ion clash',
  minimum_distance_to_waters: 'No water clash',
  volume_overlap_with_protein: 'Volume overlap protein',
  volume_overlap_with_organic_cofactors: 'Volume overlap cofactor',
  volume_overlap_with_inorganic_cofactors: 'Volume overlap ion',
  volume_overlap_with_waters: 'Volume overlap water',
  'rmsd_≤_2å': 'RMSD ≤ 2 Å',
}

// The checks PoseBusters counts towards "PB-valid" (redock config, excluding the RMSD line).
export const PB_VALIDITY_CHECKS = Object.keys(PB_CHECK_LABELS).filter((k) => k !== 'rmsd_≤_2å')
