import { useCallback, useEffect, useMemo, useState } from 'react'
import { Link, useNavigate, useSearchParams } from 'react-router-dom'
import { ArrowUpRight } from 'lucide-react'
import type { Data, Layout, PlotMouseEvent } from 'plotly.js-basic-dist-min'
import { fetchIndex, isDocking, methodPlotName, methodShortName, methodTileLabel, type IndexData, type SystemSummary } from '../lib/api'
import {
  bucketLabel, bucketOf, clampOutliers, fmt, formatMetric, metricDef, methodStats, METRICS, resultFor,
  RMSD_SUCCESS, SCATTER_Y, scatterYDef, SIMILARITY_AXIS_TITLE, similarityAt, type MetricKey, type ScatterY,
} from '../lib/stats'
import { ErrorBox, Label, MethodBadge, Select, Spinner, Stat, Tip } from '../components/ui'
import Plot from '../components/Plot'

type SortKey = 'similarity' | 'pdb' | string

export default function Overview() {
  const [data, setData] = useState<IndexData | null>(null)
  const [error, setError] = useState<unknown>(null)
  const [params, setParams] = useSearchParams()
  const [metric, setMetric] = useState<MetricKey>('success')
  const [scatterY, setScatterY] = useState<ScatterY>('rmsd')
  const [sort, setSort] = useState<SortKey>('similarity')
  const navigate = useNavigate()

  // method multi-select, persisted in the URL as ?m=af3,boltz (absent = all)
  const methods = useMemo(() => {
    if (!data) return []
    const q = params.get('m')
    const wanted = q ? new Set(q.split(',').filter(Boolean)) : null
    const sel = data.methods.filter((m) => !wanted || wanted.has(m.id))
    return sel.length ? sel : data.methods
  }, [data, params])
  const toggleMethod = useCallback((id: string) => {
    if (!data) return
    const cur = new Set(methods.map((m) => m.id))
    if (cur.has(id)) { if (cur.size === 1) return; cur.delete(id) } else cur.add(id)
    const next = new URLSearchParams(params)
    if (cur.size === data.methods.length) next.delete('m')
    else next.set('m', data.methods.filter((m) => cur.has(m.id)).map((m) => m.id).join(','))
    setParams(next, { replace: true })
  }, [data, methods, params, setParams])

  useEffect(() => {
    fetchIndex().then(setData).catch(setError)
  }, [])

  const metricD = metricDef(metric)
  const yDef = scatterYDef(scatterY)

  const systems = useMemo(() => (data ? data.systems : []), [data])
  const stats = useMemo(() => (data ? methods.map((m) => methodStats(data, systems, m.id, metric)) : []), [data, methods, systems, metric])

  const bucketPlot = useMemo<Data[]>(() => {
    if (!data) return []
    const byBucket = data.buckets.map((_, bi) => systems.filter((s) => bucketOf(similarityAt(s), data.buckets) === bi))
    return methods.map((m) => {
      const st = byBucket.map((inBucket) => methodStats(data, inBucket, m.id, metric))
      const ys = st.map((x) => (x.value == null ? null : metricD.kind === 'rate' ? x.value * 100 : x.value))
      const hover = metricD.kind === 'rate'
        ? `${m.name}<br>%{x} similarity: %{y:.0f}% ${metricD.short.replace(/^% /, '')}<br>n=%{customdata}<extra></extra>`
        : `${m.name}<br>%{x} similarity: %{y:.1f} ${metricD.unit}<br>n=%{customdata}<extra></extra>`
      const dock = isDocking(m)
      const marker = dock
        ? { color: m.color, opacity: 0.35, line: { color: m.color, width: 2 }, pattern: { shape: '/', fgcolor: m.color, bgcolor: '#ffffff', size: 6, solidity: 0.4 } }
        : { color: m.color }
      return { type: 'bar', name: methodPlotName(m), x: data.buckets.map(bucketLabel), y: ys, marker, hovertemplate: hover, customdata: st.map((x) => x.n) } as Data
    })
  }, [data, methods, systems, metric, metricD])

  const scatter = useMemo(() => {
    if (!data) return { traces: [] as Data[], cap: null as number | null, nClamped: 0 }
    const all: { s: SystemSummary; m: string; y: number; flagged: boolean }[] = []
    for (const m of methods) {
      for (const s of systems) {
        const r = resultFor(data, s.system_id, m.id)
        if (!r || !r.ok) continue
        const y = yDef.of(r)
        if (y == null || (yDef.log && y <= 0)) continue
        all.push({ s, m: m.id, y, flagged: (r.clashes_pose ?? 0) > 0 || r.pb_pass === false })
      }
    }
    const clamped = clampOutliers(all.map((p) => p.y), !yDef.log)
    const traces = methods.map((m) => {
      const idx = all.map((p, i) => (p.m === m.id ? i : -1)).filter((i) => i >= 0)
      const pts = idx.map((i) => all[i])
      return {
        type: 'scatter', mode: 'markers', name: methodPlotName(m),
        x: pts.map(({ s }) => similarityAt(s)),
        y: idx.map((i) => clamped.ys[i]),
        customdata: pts.map(({ s }, k) => [s.system_id, m.id, s.pdb_id.toUpperCase(), all[idx[k]].y]),
        marker: { size: isDocking(m) ? 8 : 7, color: m.color, symbol: idx.map((i) => (clamped.cap != null && all[i].y > clamped.cap ? 'triangle-up' : isDocking(m) ? 'diamond' : 'circle')), line: { width: pts.map((p) => (p.flagged ? 1.5 : 1)), color: pts.map((p) => (p.flagged ? '#d64545' : '#ffffff')) } },
        hovertemplate: `%{customdata[2]} · ${m.name} · %{customdata[3]:.2f} ${yDef.unit} · similarity %{x:.0f}<extra></extra>`,
      } as Data
    })
    return { traces, cap: clamped.cap, nClamped: clamped.nClamped }
  }, [data, methods, systems, yDef])

  const scatterLayout = useMemo<Partial<Layout>>(() => {
    const base: Partial<Layout> = {
      margin: { l: 48, r: 12, t: 8, b: 44 }, legend: { orientation: 'h', y: -0.38, x: 0 }, hovermode: 'closest',
      xaxis: { title: { text: 'SuCOS-pocket similarity to closest training structure', standoff: 6 }, range: [-3, 103] },
    }
    if (yDef.log) {
      return {
        ...base,
        yaxis: { title: { text: yDef.axisTitle }, type: 'log', range: [Math.log10(0.2), Math.log10(50)], tickvals: [0.2, 0.5, 1, 2, 5, 10, 20, 50], ticktext: ['0.2', '0.5', '1', '2', '5', '10', '20', '50'] },
        shapes: [{ type: 'line', xref: 'paper', x0: 0, x1: 1, y0: RMSD_SUCCESS, y1: RMSD_SUCCESS, line: { color: '#8a8a95', width: 1, dash: 'dash' } }],
        annotations: [{ xref: 'paper', x: 1, y: Math.log10(RMSD_SUCCESS), xanchor: 'right', yanchor: 'bottom', text: `${RMSD_SUCCESS} Å`, showarrow: false, font: { size: 10, color: '#8a8a95' } }],
      }
    }
    const energy = yDef.unit === 'kcal/mol'
    return {
      ...base,
      yaxis: { title: { text: yDef.axisTitle }, autorange: true, zeroline: energy, zerolinecolor: '#8a8a95', zerolinewidth: 1, ...(energy ? {} : { dtick: 1, rangemode: 'tozero' as const }) },
      shapes: energy ? [{ type: 'line', xref: 'paper', x0: 0, x1: 1, y0: 0, y1: 0, line: { color: '#8a8a95', width: 1, dash: 'dash' } }] : [],
    }
  }, [yDef])

  const onPointClick = useCallback((e: PlotMouseEvent) => {
    const cd = e.points?.[0]?.customdata as unknown as [string, string, string, number] | undefined
    if (cd) navigate(`/system/${cd[0]}?m=${cd[1]}`)
  }, [navigate])

  const rows = useMemo(() => {
    if (!data) return []
    const list = [...systems]
    if (sort === 'similarity') list.sort((a, b) => similarityAt(a) - similarityAt(b))
    else if (sort === 'pdb') list.sort((a, b) => a.pdb_id.localeCompare(b.pdb_id))
    else list.sort((a, b) => (resultFor(data, a.system_id, sort)?.rmsd ?? 99) - (resultFor(data, b.system_id, sort)?.rmsd ?? 99))
    return list
  }, [data, systems, sort])

  if (error) return <ErrorBox error={error} />
  if (!data) return <Spinner label="Loading dataset" />

  const metricOptions = METRICS.map((m) => ({ value: m.value, label: m.label }))
  const yOptions = SCATTER_Y.map((y) => ({ value: y.value, label: y.label }))

  return (
    <div className="flex flex-col gap-5">
      <div className="flex flex-wrap items-end gap-4 justify-between">
        <div>
          <h1 className="text-[18px] font-semibold">Which co-folding models learn physical priors?</h1>
          <p className="text-fg-2 mt-1 max-w-[760px]">
            Top-ranked poses from {data.methods.length} methods on {data.systems.length} Runs N&apos; Poses systems, bucketed by SuCOS-pocket
            similarity to the closest training structure. A model that has learned physical priors keeps producing plausible poses (no clashes,
            low strain, favourable interaction energy) as similarity drops, while a model leaning on memorised complexes degrades; rigid holo
            redocking with Vina is shown as a physics-only baseline. RMSD success (≤ {RMSD_SUCCESS} Å) is shown by default and the Metric
            selector switches to the physics checks.
          </p>
        </div>
      </div>

      <div className="card px-4 py-3 grid gap-4 lg:grid-cols-[1fr_auto] items-end">
        <div className="flex flex-col gap-2 min-w-0">
          <Label>Methods</Label>
          <div className="flex flex-wrap gap-1.5">
            {data.methods.map((m) => {
              const on = methods.some((x) => x.id === m.id)
              return (
                <button key={m.id} type="button" className="btn" data-active={on} style={on ? undefined : { color: 'var(--color-fg-3)' }} onClick={() => toggleMethod(m.id)} aria-pressed={on}>
                  <span className="w-2 h-2 rounded-full" style={{ background: m.color, opacity: on ? 1 : 0.4 }} />
                  {m.name}
                  <MethodBadge method={m} />
                </button>
              )
            })}
          </div>
        </div>
        <div className="flex flex-col gap-2">
          <Label>Metric</Label>
          <Select value={metric} onChange={setMetric} options={metricOptions} width={250} />
        </div>
      </div>

      <div className="grid gap-3" style={{ gridTemplateColumns: `repeat(auto-fit, minmax(150px, 1fr))` }}>
        {stats.map((s) => {
          const m = data.methods.find((x) => x.id === s.method)!
          return (
            <Stat
              key={s.method}
              label={<span className="flex items-center gap-1.5 flex-wrap"><span className="w-2 h-2 rounded-full" style={{ background: m.color }} />{methodTileLabel(m)}<MethodBadge method={m} /></span>}
              value={formatMetric(metricD, s.value)}
              sub={metric === 'success' ? `median ${fmt(s.medianRmsd)} Å · n=${s.n}` : `${metricD.kind === 'rate' ? 'rate' : `median ${metricD.unit} · lower is better`} · ${fmt(s.successRate == null ? null : s.successRate * 100, 0)}% ≤ ${RMSD_SUCCESS} Å · n=${s.n}`}
            />
          )
        })}
      </div>

      <div className="text-[12px] text-fg-3 -mt-2">{metricD.description}</div>

      <div className="grid gap-4 lg:grid-cols-2">
        <div className="card p-4">
          <div className="flex items-baseline justify-between mb-1 gap-2">
            <h2 className="font-medium truncate">{metricD.kind === 'rate' ? metricD.short.replace(/^% /, '') : metricD.short} by similarity to the closest training structure</h2>
            <span className="text-[11px] text-fg-3 whitespace-nowrap">training cutoff {data.default_cutoff} · {systems.length} systems</span>
          </div>
          <Plot
            data={bucketPlot}
            layout={{
              barmode: 'group', margin: { l: 48, r: 12, t: 8, b: 44 }, legend: { orientation: 'h', y: -0.38, x: 0 },
              yaxis: metricD.kind === 'rate' ? { title: { text: metricD.short }, range: [0, 100] } : { title: { text: metricD.short }, zeroline: true, zerolinecolor: '#8a8a95' },
              xaxis: { title: { text: SIMILARITY_AXIS_TITLE, standoff: 6 } },
            }}
          />
        </div>
        <div className="card p-4">
          <div className="flex items-center justify-between mb-1 gap-2 flex-wrap">
            <h2 className="font-medium">{yDef.label.replace(/ \(.*\)$/, '')} vs. similarity to the closest training structure</h2>
            <div className="flex items-center gap-2">
              <span className="text-[11px] text-fg-3 hidden xl:inline">click a point to open it</span>
              <Select value={scatterY} onChange={setScatterY} options={yOptions} width={210} />
            </div>
          </div>
          <div className="text-[11px] text-fg-3 mb-1">{yDef.description}</div>
          <Plot data={scatter.traces} onClick={onPointClick} layout={scatterLayout} />
          {scatter.cap != null && (
            <div className="text-[11px] text-fg-3 mt-1">{scatter.nClamped} extreme value{scatter.nClamped > 1 ? 's' : ''} clamped to {scatter.cap.toFixed(0)} {yDef.unit} (shown as ▲); hover shows the true value.</div>
          )}
        </div>
      </div>

      <div className="card overflow-hidden">
        <div className="px-4 py-3 border-b hairline flex items-center gap-3">
          <h2 className="font-medium">Systems</h2>
          <span className="text-fg-3 text-[12px] hidden md:inline-flex items-center gap-1.5">ligand RMSD (Å) · <span className="w-1.5 h-1.5 rounded-full bg-warn inline-block" /> PoseBusters violation · <span className="w-1.5 h-1.5 rounded-full bg-bad inline-block" /> pocket clash at pose</span>
          <div className="ml-auto flex items-center gap-2">
            <Label>sort</Label>
            <Select value={sort} onChange={setSort} options={[{ value: 'similarity', label: 'Similarity' }, { value: 'pdb', label: 'PDB id' }, ...methods.map((m) => ({ value: m.id, label: `${m.name} RMSD` }))]} width={150} />
          </div>
        </div>
        <div className="overflow-auto max-h-[560px]">
          <table className="data">
            <thead>
              <tr>
                <th>System</th>
                <th>Ligand</th>
                <th className="text-right" title={SIMILARITY_AXIS_TITLE}>Similarity</th>
                {methods.map((m) => <th key={m.id} className="text-right"><span className="inline-flex items-center gap-1.5 justify-end">{isDocking(m) ? methodShortName(m) : m.name}<MethodBadge method={m} /></span></th>)}
                <th></th>
              </tr>
            </thead>
            <tbody>
              {rows.map((s) => <Row key={s.system_id} s={s} data={data} methods={methods} />)}
            </tbody>
          </table>
        </div>
      </div>
    </div>
  )
}

function Row({ s, data, methods }: { s: SystemSummary; data: IndexData; methods: IndexData['methods'] }) {
  const sim = similarityAt(s)
  return (
    <tr>
      <td>
        <Link to={`/system/${s.system_id}`} className="font-medium hover:text-accent">{s.pdb_id.toUpperCase()}</Link>
      </td>
      <td>
        <Tip content={<span className="mono">{s.ligand_smiles}</span>}>
          <a className="mono text-fg-2 hover:text-accent" href={`https://www.rcsb.org/ligand/${encodeURIComponent(s.ccd)}`} target="_blank" rel="noreferrer">{s.ccd}</a>
        </Tip>
      </td>
      <td className="text-right tabular-nums">
        <span className="inline-flex items-center gap-2 justify-end">
          <span className="w-14 h-1.5 rounded-full bg-line overflow-hidden inline-block"><span className="block h-full bg-accent" style={{ width: `${sim}%` }} /></span>
          <span className="mono w-9 inline-block">{sim.toFixed(0)}</span>
        </span>
      </td>
      {methods.map((m) => {
        const r = resultFor(data, s.system_id, m.id)
        return <td key={m.id} className="text-right tabular-nums"><RmsdCell rmsd={r?.ok ? r.rmsd : null} pb={r?.pb_pass} clashes={r?.clashes_pose} /></td>
      })}
      <td className="text-right whitespace-nowrap">
        <Link className="btn" style={{ border: 'none', height: 24 }} to={`/system/${s.system_id}`}>
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
