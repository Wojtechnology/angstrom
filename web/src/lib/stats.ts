import type { IndexData, ResultRow, SystemSummary } from './api'

export const RMSD_SUCCESS = 2.0

// ---------------------------------------------------------------- similarity
/** Benchmark SuCOS-pocket similarity (training cutoff 2021-09-30) of a system to its closest training structure, 0–100. */
export function similarityAt(sys: SystemSummary): number {
  return sys.similarity
}

export const SIMILARITY_AXIS_TITLE = 'SuCOS-pocket similarity to closest training structure'

// ---------------------------------------------------------------- aggregate metrics
export type MetricKey = 'success' | 'pb_valid' | 'clash_free' | 'pocket_hit' | 'contact_retention' | 'excess_relax_de' | 'excess_strain_local'

export interface MetricDef {
  value: MetricKey
  label: string
  short: string
  kind: 'rate' | 'median'
  unit: string
  /** one-sentence explanation shown under the tiles */
  description: string
  /** value per result row, or null when not applicable */
  of: (r: ResultRow) => number | null
  /** for rates: predicate */
  ok?: (r: ResultRow) => boolean
}

export const METRICS: MetricDef[] = [
  { value: 'success', label: `RMSD ≤ ${RMSD_SUCCESS} Å (success rate)`, short: `% RMSD ≤ ${RMSD_SUCCESS} Å`, kind: 'rate', unit: '%', description: `Fraction of systems whose top-ranked pose has ligand RMSD ≤ ${RMSD_SUCCESS} Å.`, of: (r) => r.rmsd ?? null, ok: (r) => r.rmsd != null && r.rmsd <= RMSD_SUCCESS },
  { value: 'pb_valid', label: 'PoseBusters valid (rate)', short: '% PoseBusters valid', kind: 'rate', unit: '%', description: 'Fraction of poses passing every PoseBusters redock check.', of: (r) => (r.pb_pass == null ? null : r.pb_pass ? 1 : 0), ok: (r) => r.pb_pass === true },
  { value: 'clash_free', label: 'Clash-free at pose (rate)', short: '% clash-free at pose', kind: 'rate', unit: '%', description: 'Fraction of poses with no ligand–protein heavy-atom clash before relaxation.', of: (r) => (r.clashes_pose == null ? null : r.clashes_pose === 0 ? 1 : 0), ok: (r) => r.clashes_pose === 0 },
  { value: 'pocket_hit', label: 'Pocket hit rate', short: '% in the crystal pocket', kind: 'rate', unit: '%', description: 'Fraction of poses that land in the crystal binding site, judged by shape overlap with the crystal ligand (≥ 0.05) or centroid distance ≤ 4 Å.', of: (r) => (r.pocket_hit == null ? null : r.pocket_hit ? 1 : 0), ok: (r) => r.pocket_hit === true },
  { value: 'contact_retention', label: 'Median contact retention', short: 'median contact retention', kind: 'median', unit: '', description: "Median fraction of the crystal ligand's protein contacts (heavy-atom pairs ≤ 4 Å, per ligand atom and residue) that the predicted pose keeps.", of: (r) => r.contact_retention ?? null },
  { value: 'excess_relax_de', label: 'Median excess relaxation ΔE', short: 'median excess relaxation ΔE vs crystal (kcal/mol)', kind: 'median', unit: 'kcal/mol', description: "Median of the relaxation ΔE minus the crystal pose's own ΔE under the same force field (MMFF94s, kcal/mol): the crystal complex is minimised the same way and its ΔE subtracted, so ≈ 0 means as physically reasonable as the experimental pose and large positive values mean the model added strain or clashes beyond what the crystal itself carries.", of: (r) => r.excess_relaxation_de ?? null },
  { value: 'excess_strain_local', label: 'Median excess ligand strain', short: 'median excess ligand strain vs crystal (kcal/mol)', kind: 'median', unit: 'kcal/mol', description: "Median of the local ligand strain minus the crystal ligand's own strain under the same force field (MMFF94s, kcal/mol): ≈ 0 means as physically reasonable as the experimental pose; large positive values mean the model added strain beyond what the crystal itself carries.", of: (r) => r.excess_strain_local ?? null },
]

export function metricDef(key: MetricKey): MetricDef {
  return METRICS.find((m) => m.value === key) ?? METRICS[0]
}

export function median(xs: number[]): number | null {
  if (!xs.length) return null
  const s = [...xs].sort((a, b) => a - b)
  const m = Math.floor(s.length / 2)
  return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2
}

export interface MethodStats {
  method: string
  n: number
  /** rate (0–1) or median, depending on the metric */
  value: number | null
  successRate: number | null
  medianRmsd: number | null
}

export function okRows(data: IndexData, systems: SystemSummary[], method: string): ResultRow[] {
  const ids = new Set(systems.map((s) => s.system_id))
  return data.results.filter((r) => r.method === method && r.ok && ids.has(r.system_id) && r.rmsd != null)
}

export function methodStats(data: IndexData, systems: SystemSummary[], method: string, metric: MetricKey = 'success'): MethodStats {
  const rows = okRows(data, systems, method)
  const rmsds = rows.map((r) => r.rmsd as number)
  const def = metricDef(metric)
  const applicable = rows.filter((r) => def.of(r) != null)
  const value = def.kind === 'rate'
    ? (applicable.length ? applicable.filter((r) => def.ok!(r)).length / applicable.length : null)
    : median(applicable.map((r) => def.of(r) as number))
  return {
    method,
    n: applicable.length,
    value,
    successRate: rows.length ? rmsds.filter((x) => x <= RMSD_SUCCESS).length / rows.length : null,
    medianRmsd: median(rmsds),
  }
}

export function formatMetric(def: MetricDef, v: number | null): string {
  if (v == null) return '–'
  if (def.kind === 'rate') return `${Math.round(v * 100)}%`
  return def.unit === '' ? v.toFixed(2) : v.toFixed(1)
}

// ---------------------------------------------------------------- scatter y axes
export type ScatterY = 'rmsd' | 'contact_retention' | 'shape_overlap' | 'excess_relax_de' | 'excess_strain_local' | 'clashes_pose' | 'pb_fails'

export interface ScatterYDef {
  value: ScatterY
  label: string
  axisTitle: string
  log: boolean
  unit: string
  /** one-sentence explanation shown under the scatter */
  description: string
  of: (r: ResultRow) => number | null
}

export const SCATTER_Y: ScatterYDef[] = [
  { value: 'rmsd', label: 'Ligand RMSD (Å, log)', axisTitle: 'ligand RMSD (Å)', log: true, unit: 'Å', description: 'Symmetry-corrected heavy-atom RMSD between the predicted and crystal ligand after superposing the predicted receptor on the crystal binding-site Cα atoms; ≤ 2 Å counts as success.', of: (r) => r.rmsd ?? null },
  { value: 'contact_retention', label: 'Contact retention', axisTitle: 'contact retention (fraction of crystal contacts kept)', log: false, unit: '', description: "Fraction of the crystal ligand's protein contacts (heavy-atom pairs ≤ 4 Å, per ligand atom and residue) that the predicted pose keeps; 1 means every crystal contact is reproduced.", of: (r) => r.contact_retention ?? null },
  { value: 'shape_overlap', label: 'Shape overlap with crystal ligand', axisTitle: 'shape overlap with crystal ligand (0–1)', log: false, unit: '', description: 'Volume overlap between the predicted and crystal ligand poses (0–1); values near 0 mean the pose missed the binding site.', of: (r) => r.shape_overlap ?? null },
  { value: 'excess_relax_de', label: 'Excess relaxation ΔE', axisTitle: 'excess relaxation ΔE vs crystal pose (kcal/mol)', log: false, unit: 'kcal/mol', description: "Relaxation ΔE of the prediction minus the crystal pose's own ΔE under the same force field; ≈ 0 means as physically reasonable as the experimental pose, large positive values mean strain or clashes the model added.", of: (r) => r.excess_relaxation_de ?? null },
  { value: 'excess_strain_local', label: 'Excess ligand strain', axisTitle: 'excess local ligand strain vs crystal ligand (kcal/mol)', log: false, unit: 'kcal/mol', description: "Local ligand strain of the prediction minus the crystal ligand's own strain under the same force field; ≈ 0 means as physically reasonable as the experimental pose.", of: (r) => r.excess_strain_local ?? null },
  { value: 'clashes_pose', label: 'Pocket clashes at pose', axisTitle: 'pocket clashes at pose', log: false, unit: '', description: 'Number of ligand–protein heavy-atom pairs closer than 0.75 × the sum of their van der Waals radii, before any relaxation.', of: (r) => r.clashes_pose ?? null },
  { value: 'pb_fails', label: 'PoseBusters failures', axisTitle: 'PoseBusters checks failed', log: false, unit: '', description: 'Number of failed PoseBusters redock checks (geometry, stereochemistry, intra- and intermolecular clashes) for the predicted pose.', of: (r) => r.pb_fail_count ?? null },
]

export function scatterYDef(key: ScatterY): ScatterYDef {
  return SCATTER_Y.find((y) => y.value === key) ?? SCATTER_Y[0]
}

/** Clamp extreme outliers for linear energy axes (robust IQR fence) so a few huge clashes don't squash the axis. */
export function clampOutliers(ys: number[], enable: boolean): { ys: number[]; cap: number | null; nClamped: number } {
  if (!enable || ys.length < 8) return { ys, cap: null, nClamped: 0 }
  const s = [...ys].sort((a, b) => a - b)
  const q = (f: number) => s[Math.min(s.length - 1, Math.floor(s.length * f))]
  const iqr = q(0.75) - q(0.25)
  const cap = Math.max(q(0.75) + 3 * iqr, q(0.9) * 1.5, 25)
  const n = ys.filter((y) => y > cap).length
  return n ? { ys: ys.map((y) => Math.min(y, cap)), cap, nClamped: n } : { ys, cap: null, nClamped: 0 }
}

// ---------------------------------------------------------------- misc
export function bucketLabel([lo, hi]: [number, number]): string {
  return `${lo}–${Math.min(hi, 100)}`
}

export function bucketOf(sim: number, buckets: [number, number][]): number {
  for (let i = 0; i < buckets.length; i++) {
    const [lo, hi] = buckets[i]
    if (sim >= lo && sim < hi) return i
  }
  return buckets.length - 1
}

export function resultFor(data: IndexData, systemId: string, method: string): ResultRow | undefined {
  return data.results.find((r) => r.system_id === systemId && r.method === method)
}

export const fmt = (x: number | null | undefined, d = 2) => (x == null || Number.isNaN(x) ? '–' : x.toFixed(d))
export const pct = (x: number | null | undefined) => (x == null ? '–' : `${Math.round(x * 100)}%`)
