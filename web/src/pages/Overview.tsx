import { useEffect, useMemo, useState } from 'react'
import { Link } from 'react-router-dom'
import { ArrowUpRight, ExternalLink } from 'lucide-react'
import type { Data } from 'plotly.js-basic-dist-min'
import { fetchIndex, rcsbUrl, type IndexData, type SystemSummary } from '../lib/api'
import { bucketLabel, bucketOf, filterSystems, fmt, methodStats, pct, resultFor, RMSD_SUCCESS, similarityAt } from '../lib/stats'
import { ErrorBox, Label, Select, Slider, Spinner, Stat, Tip } from '../components/ui'
import Plot from '../components/Plot'

type SortKey = 'similarity' | 'pdb' | string

export default function Overview() {
  const [data, setData] = useState<IndexData | null>(null)
  const [error, setError] = useState<unknown>(null)
  const [threshold, setThreshold] = useState(100)
  const [cutoff, setCutoff] = useState<string>('')
  const [sort, setSort] = useState<SortKey>('similarity')

  useEffect(() => {
    fetchIndex().then((d) => { setData(d); setCutoff(d.default_cutoff) }).catch(setError)
  }, [])

  const systems = useMemo(() => (data && cutoff ? filterSystems(data, cutoff, threshold) : []), [data, cutoff, threshold])
  const stats = useMemo(() => (data ? data.methods.map((m) => methodStats(data, systems, m.id)) : []), [data, systems])

  const bucketPlot = useMemo<Data[]>(() => {
    if (!data || !cutoff) return []
    return data.methods.map((m) => {
      const ys = data.buckets.map((_, bi) => {
        const inBucket = data.systems.filter((s) => bucketOf(similarityAt(s, cutoff), data.buckets) === bi)
        return methodStats(data, inBucket, m.id).successRate
      })
      const ns = data.buckets.map((_, bi) => data.systems.filter((s) => bucketOf(similarityAt(s, cutoff), data.buckets) === bi).length)
      return {
        type: 'bar', name: m.name, x: data.buckets.map(bucketLabel), y: ys.map((y) => (y == null ? null : y * 100)),
        marker: { color: m.color }, hovertemplate: `${m.name}<br>%{x} similarity: %{y:.0f}% ≤ ${RMSD_SUCCESS} Å<br>n=%{customdata}<extra></extra>`,
        customdata: ns,
      } as Data
    })
  }, [data, cutoff])

  const cutoffPlot = useMemo<Data[]>(() => {
    if (!data) return []
    return data.methods.map((m) => {
      const ys = data.cutoffs.map((c) => methodStats(data, filterSystems(data, c, threshold), m.id).successRate)
      const ns = data.cutoffs.map((c) => filterSystems(data, c, threshold).length)
      return {
        type: 'scatter', mode: 'lines+markers', name: m.name, x: data.cutoffs, y: ys.map((y) => (y == null ? null : y * 100)),
        line: { color: m.color, width: 2 }, marker: { size: 6 },
        hovertemplate: `${m.name}<br>cutoff %{x}: %{y:.0f}% ≤ ${RMSD_SUCCESS} Å<br>n=%{customdata}<extra></extra>`, customdata: ns,
      } as Data
    })
  }, [data, threshold])

  const rows = useMemo(() => {
    if (!data || !cutoff) return []
    const list = [...systems]
    if (sort === 'similarity') list.sort((a, b) => similarityAt(a, cutoff) - similarityAt(b, cutoff))
    else if (sort === 'pdb') list.sort((a, b) => a.pdb_id.localeCompare(b.pdb_id))
    else list.sort((a, b) => (resultFor(data, a.system_id, sort)?.rmsd ?? 99) - (resultFor(data, b.system_id, sort)?.rmsd ?? 99))
    return list
  }, [data, systems, sort, cutoff])

  if (error) return <ErrorBox error={error} />
  if (!data || !cutoff) return <Spinner label="Loading dataset" />

  const cutoffOptions = data.cutoffs.map((c) => ({ value: c, label: c === data.default_cutoff ? `${c} (benchmark)` : c }))

  return (
    <div className="flex flex-col gap-5">
      <div className="flex flex-wrap items-end gap-4 justify-between">
        <div>
          <h1 className="text-[18px] font-semibold">Do co-folding models generalise, or remember?</h1>
          <p className="text-fg-2 mt-1 max-w-[720px]">
            Top-ranked poses from {data.methods.length} methods on {data.systems.length} Runs N&apos; Poses systems, bucketed by the
            SuCOS-pocket similarity of each system to its closest structure in the training set. Success is ligand RMSD ≤ {RMSD_SUCCESS} Å after
            binding-site superposition.
          </p>
        </div>
      </div>

      <div className="card px-4 py-3 grid gap-4 md:grid-cols-[1fr_auto] items-center">
        <div className="flex flex-col gap-2">
          <div className="flex items-center justify-between">
            <Label>Max similarity to training set</Label>
            <span className="mono text-fg-2">≤ {threshold} · {systems.length} of {data.systems.length} systems</span>
          </div>
          <Slider value={threshold} onChange={setThreshold} min={0} max={100} step={1} />
          <div className="flex justify-between text-[11px] text-fg-3"><span>novel</span><span>seen before</span></div>
        </div>
        <div className="flex flex-col gap-2">
          <Label>Training cutoff date</Label>
          <Select value={cutoff} onChange={setCutoff} options={cutoffOptions} width={190} />
        </div>
      </div>

      <div className="grid gap-3" style={{ gridTemplateColumns: `repeat(auto-fit, minmax(150px, 1fr))` }}>
        {stats.map((s) => {
          const m = data.methods.find((x) => x.id === s.method)!
          return (
            <Stat
              key={s.method}
              label={<span className="flex items-center gap-1.5"><span className="w-2 h-2 rounded-full" style={{ background: m.color }} />{m.name}</span>}
              value={pct(s.successRate)}
              sub={`median ${fmt(s.medianRmsd)} Å · n=${s.n}`}
            />
          )
        })}
      </div>

      <div className="grid gap-4 lg:grid-cols-2">
        <div className="card p-4">
          <div className="flex items-baseline justify-between mb-1">
            <h2 className="font-medium">Success rate by similarity bucket</h2>
            <span className="text-[11px] text-fg-3">cutoff {cutoff} · all systems</span>
          </div>
          <Plot data={bucketPlot} layout={{ barmode: 'group', yaxis: { title: { text: `% RMSD ≤ ${RMSD_SUCCESS} Å` }, range: [0, 100] }, xaxis: { title: { text: 'SuCOS-pocket similarity to closest training structure' } } }} />
        </div>
        <div className="card p-4">
          <div className="flex items-baseline justify-between mb-1">
            <h2 className="font-medium">Success rate vs. assumed training cutoff</h2>
            <span className="text-[11px] text-fg-3">systems with similarity ≤ {threshold} under each cutoff</span>
          </div>
          <Plot data={cutoffPlot} layout={{ yaxis: { title: { text: `% RMSD ≤ ${RMSD_SUCCESS} Å` }, range: [0, 100] }, xaxis: { type: 'category', tickvals: data.cutoffs, ticktext: data.cutoffs.map(shortDate), tickfont: { size: 10 } } }} />
        </div>
      </div>

      <div className="card overflow-hidden">
        <div className="px-4 py-3 border-b hairline flex items-center gap-3">
          <h2 className="font-medium">Systems</h2>
          <span className="text-fg-3 text-[12px] hidden md:inline-flex items-center gap-1.5">click a row for the 3D breakdown · cells show ligand RMSD (Å) · <span className="w-1.5 h-1.5 rounded-full bg-warn inline-block" /> PoseBusters violation · <span className="w-1.5 h-1.5 rounded-full bg-bad inline-block" /> pocket clash at pose</span>
          <div className="ml-auto flex items-center gap-2">
            <Label>sort</Label>
            <Select value={sort} onChange={setSort} options={[{ value: 'similarity', label: 'Similarity' }, { value: 'pdb', label: 'PDB id' }, ...data.methods.map((m) => ({ value: m.id, label: `${m.name} RMSD` }))]} width={150} />
          </div>
        </div>
        <div className="overflow-auto max-h-[560px]">
          <table className="data">
            <thead>
              <tr>
                <th>System</th>
                <th>Ligand</th>
                <th className="text-right">Similarity</th>
                {data.methods.map((m) => <th key={m.id} className="text-right">{m.name}</th>)}
                <th></th>
              </tr>
            </thead>
            <tbody>
              {rows.map((s) => <Row key={s.system_id} s={s} data={data} cutoff={cutoff} />)}
            </tbody>
          </table>
        </div>
      </div>
    </div>
  )
}

function Row({ s, data, cutoff }: { s: SystemSummary; data: IndexData; cutoff: string }) {
  const sim = similarityAt(s, cutoff)
  return (
    <tr>
      <td>
        <Link to={`/system/${s.system_id}`} className="font-medium hover:text-accent">{s.pdb_id.toUpperCase()}</Link>
        <span className="mono text-fg-3 ml-2 hidden xl:inline">{s.system_id}</span>
      </td>
      <td>
        <Tip content={<span className="mono">{s.ligand_smiles}</span>}>
          <span className="mono text-fg-2">{s.ccd} · {s.n_heavy} heavy atoms</span>
        </Tip>
      </td>
      <td className="text-right tabular-nums">
        <span className="inline-flex items-center gap-2 justify-end">
          <span className="w-14 h-1.5 rounded-full bg-line overflow-hidden inline-block"><span className="block h-full bg-accent" style={{ width: `${sim}%` }} /></span>
          <span className="mono w-9 inline-block">{sim.toFixed(0)}</span>
        </span>
      </td>
      {data.methods.map((m) => {
        const r = resultFor(data, s.system_id, m.id)
        return <td key={m.id} className="text-right tabular-nums"><RmsdCell rmsd={r?.ok ? r.rmsd : null} pb={r?.pb_pass} clashes={r?.clashes_pose} /></td>
      })}
      <td className="text-right whitespace-nowrap">
        <a className="btn" style={{ border: 'none', height: 24 }} href={rcsbUrl(s.pdb_id)} target="_blank" rel="noreferrer" title="Ground truth on RCSB">
          <ExternalLink size={12} /> RCSB
        </a>
        <Link className="btn ml-1" style={{ border: 'none', height: 24 }} to={`/system/${s.system_id}`}>
          <ArrowUpRight size={12} /> 3D
        </Link>
      </td>
    </tr>
  )
}

export function RmsdCell({ rmsd, pb, clashes }: { rmsd: number | null | undefined; pb?: boolean | null; clashes?: number | null }) {
  if (rmsd == null) return <span className="text-fg-3">–</span>
  const good = rmsd <= RMSD_SUCCESS
  const notes = [pb === false ? 'PoseBusters violations' : null, clashes ? `${clashes} pocket clash${clashes > 1 ? 'es' : ''} at pose` : null].filter(Boolean).join(' · ')
  return (
    <span className={`chip ${good ? 'chip-ok' : 'chip-bad'}`} title={notes}>
      {rmsd.toFixed(2)}
      {pb === false && <span className="w-1.5 h-1.5 rounded-full bg-warn" />}
      {!!clashes && clashes > 0 && <span className="w-1.5 h-1.5 rounded-full bg-bad" />}
    </span>
  )
}

function shortDate(iso: string): string {
  const [y, m] = iso.split('-')
  const months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
  return `${months[Number(m) - 1]} ${y}`
}
