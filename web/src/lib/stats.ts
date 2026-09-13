import type { IndexData, ResultRow, SystemSummary } from './api'

export const RMSD_SUCCESS = 2.0

export function median(xs: number[]): number | null {
  if (!xs.length) return null
  const s = [...xs].sort((a, b) => a - b)
  const m = Math.floor(s.length / 2)
  return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2
}

export function similarityAt(sys: SystemSummary, cutoff: string): number {
  const v = sys.similarity_by_cutoff[cutoff]
  return v == null ? sys.similarity : v
}

/** Systems whose max similarity to any training structure released before `cutoff` is <= threshold. */
export function filterSystems(data: IndexData, cutoff: string, threshold: number): SystemSummary[] {
  return data.systems.filter((s) => similarityAt(s, cutoff) <= threshold)
}

export interface MethodStats {
  method: string
  n: number
  successRate: number | null
  medianRmsd: number | null
  pbValidRate: number | null
  medianStrain: number | null
}

export function methodStats(data: IndexData, systems: SystemSummary[], method: string): MethodStats {
  const ids = new Set(systems.map((s) => s.system_id))
  const rows = data.results.filter((r) => r.method === method && r.ok && ids.has(r.system_id) && r.rmsd != null)
  const rmsds = rows.map((r) => r.rmsd as number)
  const pb = rows.filter((r) => r.pb_pass != null)
  const strains = rows.map((r) => r.strain_local).filter((x): x is number => x != null)
  return {
    method,
    n: rows.length,
    successRate: rows.length ? rmsds.filter((x) => x <= RMSD_SUCCESS).length / rows.length : null,
    medianRmsd: median(rmsds),
    pbValidRate: pb.length ? pb.filter((r) => r.pb_pass).length / pb.length : null,
    medianStrain: median(strains),
  }
}

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
