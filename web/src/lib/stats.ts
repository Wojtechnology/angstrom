import type { IndexData, ResultRow, SystemSummary } from './api'

export const RMSD_SUCCESS = 2.0

// ---------------------------------------------------------------- similarity
/** Benchmark SuCOS-pocket similarity (training cutoff 2021-09-30) of a system to its closest training structure, 0–100. */
export function similarityAt(sys: SystemSummary): number {
  return sys.similarity
}

export const SIMILARITY_AXIS_TITLE = 'SuCOS-pocket similarity to closest training structure'

// ---------------------------------------------------------------- aggregate metrics
export type MetricKey = 'success' | 'pb_valid' | 'clash_free' | 'e_int_pose' | 'relax_de' | 'strain_local'

export interface MetricDef {
  value: MetricKey
  label: string
  short: string
  kind: 'rate' | 'median'
  unit: string
  /** value per result row, or null when not applicable */
  of: (r: ResultRow) => number | null
  /** for rates: predicate */
  ok?: (r: ResultRow) => boolean
}

export const METRICS: MetricDef[] = [
  { value: 'success', label: `RMSD ≤ ${RMSD_SUCCESS} Å (success rate)`, short: `% RMSD ≤ ${RMSD_SUCCESS} Å`, kind: 'rate', unit: '%', of: (r) => r.rmsd ?? null, ok: (r) => r.rmsd != null && r.rmsd <= RMSD_SUCCESS },
  { value: 'pb_valid', label: 'PoseBusters valid (rate)', short: '% PoseBusters valid', kind: 'rate', unit: '%', of: (r) => (r.pb_pass == null ? null : r.pb_pass ? 1 : 0), ok: (r) => r.pb_pass === true },
  { value: 'clash_free', label: 'Clash-free at pose (rate)', short: '% clash-free at pose', kind: 'rate', unit: '%', of: (r) => (r.clashes_pose == null ? null : r.clashes_pose === 0 ? 1 : 0), ok: (r) => r.clashes_pose === 0 },
  { value: 'e_int_pose', label: 'Median interaction energy at pose (kcal/mol, lower is better)', short: 'median E_int at pose (kcal/mol, lower is better)', kind: 'median', unit: 'kcal/mol', of: (r) => r.e_interaction_pose ?? null },
  { value: 'relax_de', label: 'Median relaxation ΔE, pose − min (lower is better)', short: 'median relaxation ΔE (kcal/mol, lower is better)', kind: 'median', unit: 'kcal/mol', of: (r) => (r.e_interaction_pose != null && r.e_interaction_min != null ? r.e_interaction_pose - r.e_interaction_min : null) },
  { value: 'strain_local', label: 'Median ligand strain, local (kcal/mol, lower is better)', short: 'median local strain (kcal/mol, lower is better)', kind: 'median', unit: 'kcal/mol', of: (r) => r.strain_local ?? null },
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
  return def.kind === 'rate' ? `${Math.round(v * 100)}%` : `${v.toFixed(1)}`
}

// ---------------------------------------------------------------- scatter y axes
export type ScatterY = 'rmsd' | 'e_int_pose' | 'e_int_min' | 'relax_de' | 'strain_local' | 'strain_global' | 'clashes_pose' | 'pb_fails'

export interface ScatterYDef {
  value: ScatterY
  label: string
  axisTitle: string
  log: boolean
  unit: string
  of: (r: ResultRow) => number | null
}

export const SCATTER_Y: ScatterYDef[] = [
  { value: 'rmsd', label: 'Ligand RMSD (Å, log)', axisTitle: 'ligand RMSD (Å)', log: true, unit: 'Å', of: (r) => r.rmsd ?? null },
  { value: 'e_int_pose', label: 'Interaction energy at pose', axisTitle: 'E_int at pose (kcal/mol, lower is better)', log: false, unit: 'kcal/mol', of: (r) => r.e_interaction_pose ?? null },
  { value: 'e_int_min', label: 'Interaction energy after relaxation', axisTitle: 'E_int after relaxation (kcal/mol, lower is better)', log: false, unit: 'kcal/mol', of: (r) => r.e_interaction_min ?? null },
  { value: 'relax_de', label: 'Relaxation ΔE (pose − min)', axisTitle: 'relaxation ΔE (kcal/mol, lower is better)', log: false, unit: 'kcal/mol', of: (r) => (r.e_interaction_pose != null && r.e_interaction_min != null ? r.e_interaction_pose - r.e_interaction_min : null) },
  { value: 'strain_local', label: 'Ligand strain, local', axisTitle: 'local ligand strain (kcal/mol, lower is better)', log: false, unit: 'kcal/mol', of: (r) => r.strain_local ?? null },
  { value: 'strain_global', label: 'Ligand strain, global', axisTitle: 'global ligand strain (kcal/mol, lower is better)', log: false, unit: 'kcal/mol', of: (r) => r.strain_global ?? null },
  { value: 'clashes_pose', label: 'Pocket clashes at pose', axisTitle: 'pocket clashes at pose', log: false, unit: '', of: (r) => r.clashes_pose ?? null },
  { value: 'pb_fails', label: 'PoseBusters failures', axisTitle: 'PoseBusters checks failed', log: false, unit: '', of: (r) => r.pb_fail_count ?? null },
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
