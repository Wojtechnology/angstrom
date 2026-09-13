import { useCallback, useEffect, useMemo, useState } from 'react'
import { Link, useNavigate } from 'react-router-dom'
import { ArrowUpRight, ExternalLink } from 'lucide-react'
import type { Data, PlotMouseEvent } from 'plotly.js-basic-dist-min'
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
  const navigate = useNavigate()

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

  const scatterPlot = useMemo<Data[]>(() => {
    if (!data || !cutoff) return []
    return data.methods.map((m) => {
      const pts = systems
        .map((s) => ({ s, r: resultFor(data, s.system_id, m.id) }))
        .filter((x): x is { s: SystemSummary; r: NonNullable<ReturnType<typeof resultFor>> } => !!x.r && x.r.ok && x.r.rmsd != null)
      const flagged = pts.map(({ r }) => (r.clashes_pose ?? 0) > 0 || r.pb_pass === false)
      return {
        type: 'scatter', mode: 'markers', name: m.name,
        x: pts.map(({ s }) => similarityAt(s, cutoff)),
        y: pts.map(({ r }) => r.rmsd as number),
        customdata: pts.map(({ s }) => [s.system_id, m.id, s.pdb_id.toUpperCase()]),
        marker: { size: 7, color: m.color, line: { width: flagged.map((f) => (f ? 1.5 : 1)), color: flagged.map((f) => (f ? '#d64545' : '#ffffff')) } },
        hovertemplate: `%{customdata[2]} · ${m.name} · %{y:.2f} Å · similarity %{x:.0f}<extra></extra>`,
      } as Data
    })
  }, [data, cutoff, systems])

  const onPointClick = useCallback((e: PlotMouseEvent) => {
    const cd = e.points?.[0]?.customdata as unknown as [string, string, string] | undefined
    if (cd) navigate(`/system/${cd[0]}?m=${cd[1]}`)
  }, [navigate])

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
          <Plot data={bucketPlot} layout={{ barmode: 'group', margin: { l: 44, r: 12, t: 8, b: 44 }, legend: { orientation: 'h', y: -0.38, x: 0 }, yaxis: { title: { text: `% RMSD ≤ ${RMSD_SUCCESS} Å` }, range: [0, 100] }, xaxis: { title: { text: 'SuCOS-pocket similarity to closest training structure', standoff: 6 } } }} />
        </div>
        <div className="card p-4">
          <div className="flex items-baseline justify-between mb-1">
            <h2 className="font-medium">Ligand RMSD vs. similarity to training set</h2>
            <span className="text-[11px] text-fg-3">cutoff {cutoff} · similarity ≤ {threshold} · click a point to open it</span>
          </div>
          <Plot
            data={scatterPlot}
            onClick={onPointClick}
            layout={{
              margin: { l: 44, r: 12, t: 8, b: 44 }, legend: { orientation: 'h', y: -0.38, x: 0 }, hovermode: 'closest',
              xaxis: { title: { text: 'SuCOS-pocket similarity at this cutoff', standoff: 6 }, range: [-3, 103] },
              yaxis: { title: { text: 'ligand RMSD (Å)' }, type: 'log', range: [Math.log10(0.2), Math.log10(50)], tickvals: [0.2, 0.5, 1, 2, 5, 10, 20, 50], ticktext: ['0.2', '0.5', '1', '2', '5', '10', '20', '50'] },
              shapes: [{ type: 'line', xref: 'paper', x0: 0, x1: 1, y0: RMSD_SUCCESS, y1: RMSD_SUCCESS, line: { color: '#8a8a95', width: 1, dash: 'dash' } }],
              annotations: [{ xref: 'paper', x: 1, y: Math.log10(RMSD_SUCCESS), xanchor: 'right', yanchor: 'bottom', text: `${RMSD_SUCCESS} Å`, showarrow: false, font: { size: 10, color: '#8a8a95' } }],
            }}
          />
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
