import { useCallback, useEffect, useMemo, useRef, useState } from 'react'
import { Link, useParams, useSearchParams } from 'react-router-dom'
import { ArrowLeft, ExternalLink, Pause, Play, Focus, Check, X, Minus } from 'lucide-react'
import type { Data } from 'plotly.js-basic-dist-min'
import { fetchIndex, fetchStructure, fetchSystem, PB_CHECK_LABELS, PB_VALIDITY_CHECKS, rcsbUrl, type IndexData, type MethodDetail, type SystemDetail } from '../lib/api'
import { fmt, RMSD_SUCCESS } from '../lib/stats'
import { ErrorBox, Label, Spinner, Switch, Tip } from '../components/ui'
import Viewer3D, { type ViewerHandle } from '../components/Viewer3D'
import Plot from '../components/Plot'

interface Structures { gtReceptor: string; gtLigand: string; predReceptor: string | null; predLigand: string | null; predTraj: string | null }

export default function SystemView() {
  const { id = '' } = useParams()
  const [params, setParams] = useSearchParams()
  const [index, setIndex] = useState<IndexData | null>(null)
  const [sys, setSys] = useState<SystemDetail | null>(null)
  const [error, setError] = useState<unknown>(null)
  const [structures, setStructures] = useState<Structures | null>(null)
  const [showPredReceptor, setShowPredReceptor] = useState(false)
  const [showGtLigand, setShowGtLigand] = useState(true)
  const [showPocket, setShowPocket] = useState(true)
  const [showViolations, setShowViolations] = useState(true)
  const [trajectoryMode, setTrajectoryMode] = useState(false)
  const [playing, setPlaying] = useState(false)
  const [frame, setFrame] = useState(0)
  const [highlight, setHighlight] = useState<number[] | null>(null)
  const viewer = useRef<ViewerHandle>(null)

  const method = params.get('m') || ''
  const setMethod = useCallback((m: string) => setParams({ m }, { replace: true }), [setParams])

  useEffect(() => {
    Promise.all([fetchIndex(), fetchSystem(id)]).then(([i, s]) => {
      setIndex(i); setSys(s)
      if (!params.get('m')) {
        const first = i.methods.find((m) => s.methods[m.id]?.ok)?.id ?? i.methods[0].id
        setParams({ m: first }, { replace: true })
      }
    }).catch(setError)
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [id])

  const detail: MethodDetail | undefined = sys?.methods[method]

  useEffect(() => {
    if (!sys || !method) return
    let cancelled = false
    const d = sys.methods[method]
    const files = d?.ok && d.files ? d.files : null
    Promise.all([
      fetchStructure(id, sys.gt.files.receptor),
      fetchStructure(id, sys.gt.files.ligand),
      files ? fetchStructure(id, files.receptor) : Promise.resolve(null),
      files ? fetchStructure(id, files.ligand) : Promise.resolve(null),
      files ? fetchStructure(id, files.traj).catch(() => null) : Promise.resolve(null),
    ]).then(([gtReceptor, gtLigand, predReceptor, predLigand, predTraj]) => {
      if (!cancelled) setStructures({ gtReceptor, gtLigand, predReceptor, predLigand, predTraj })
    }).catch(setError)
    return () => { cancelled = true }
  }, [sys, method, id])

  // trajectory playback
  const nFrames = detail?.minimisation?.n_frames ?? 0
  useEffect(() => {
    if (!playing || !trajectoryMode || nFrames < 2) return
    const t = setInterval(() => setFrame((f) => (f + 1) % nFrames), 120)
    return () => clearInterval(t)
  }, [playing, trajectoryMode, nFrames])
  useEffect(() => { setFrame(0); setPlaying(false); setHighlight(null) }, [method])

  const energyPlot = useMemo<Data[]>(() => {
    if (!detail?.minimisation) return []
    const m = detail.minimisation
    const gt = sys?.gt.minimisation
    const rel = (e: number[]) => e.map((x) => x - e[e.length - 1])
    const out: Data[] = [{ type: 'scatter', mode: 'lines', name: 'prediction', x: m.energies.map((_, i) => i), y: rel(m.energies), line: { color: methodColor(index, method), width: 2 }, hovertemplate: 'frame %{x}: +%{y:.1f} kcal/mol<extra></extra>' }]
    if (gt) out.push({ type: 'scatter', mode: 'lines', name: 'ground truth', x: gt.energies.map((_, i) => i), y: rel(gt.energies), line: { color: '#2f9e6b', width: 1.5, dash: 'dot' }, hovertemplate: 'frame %{x}: +%{y:.1f} kcal/mol<extra></extra>' })
    if (trajectoryMode) out.push({ type: 'scatter', mode: 'markers', name: 'current', x: [frame], y: [rel(m.energies)[frame]], marker: { color: '#1c1c22', size: 8 }, showlegend: false, hoverinfo: 'skip' })
    return out
  }, [detail, sys, frame, trajectoryMode, index, method])

  if (error) return <ErrorBox error={error} />
  if (!index || !sys || !method) return <Spinner label="Loading system" />

  const color = methodColor(index, method)
  const pbRows = detail?.posebusters ? PB_VALIDITY_CHECKS.filter((k) => k in detail.posebusters!) : []
  const pbFails = pbRows.filter((k) => detail!.posebusters![k] === false)

  return (
    <div className="flex flex-col gap-4">
      <div className="flex flex-wrap items-center gap-3">
        <Link to="/" className="btn" style={{ border: 'none' }}><ArrowLeft size={14} /> Overview</Link>
        <h1 className="text-[16px] font-semibold flex items-center gap-2">
          {sys.pdb_id.toUpperCase()} <span className="mono text-fg-3 font-normal">{sys.system_id}</span>
        </h1>
        <a className="btn" href={rcsbUrl(sys.pdb_id)} target="_blank" rel="noreferrer"><ExternalLink size={12} /> RCSB {sys.pdb_id.toUpperCase()}</a>
        {sys.closest_training.pdb_id && (
          <Tip content={`Closest training system ${sys.closest_training.system_id} (released ${sys.closest_training.release_date})`}>
            <a className="btn" href={rcsbUrl(sys.closest_training.pdb_id)} target="_blank" rel="noreferrer"><ExternalLink size={12} /> closest train {sys.closest_training.pdb_id.toUpperCase()}</a>
          </Tip>
        )}
        <span className="chip chip-muted">similarity {sys.similarity.toFixed(0)}</span>
        <span className="chip chip-muted">{sys.ccd} · {sys.n_heavy} heavy atoms</span>
        <span className="chip chip-muted">{sys.seq_len} residues · {sys.n_protein_chains} chain{sys.n_protein_chains > 1 ? 's' : ''}</span>
      </div>

      {/* method switcher */}
      <div className="flex flex-wrap gap-1.5">
        {index.methods.map((m) => {
          const d = sys.methods[m.id]
          const active = m.id === method
          return (
            <button key={m.id} className="btn" data-active={active} disabled={!d?.ok} onClick={() => setMethod(m.id)} title={d?.ok ? '' : d?.error ?? 'not available'}>
              <span className="w-2 h-2 rounded-full" style={{ background: m.color }} />
              {m.name}
              {d?.ok && <span className={`mono ${d.rmsd! <= RMSD_SUCCESS ? 'text-ok' : 'text-bad'}`}>{fmt(d.rmsd)} Å</span>}
              {d?.ok && d.pb_pass === false && <span className="w-1.5 h-1.5 rounded-full bg-warn" title="PoseBusters violations" />}
            </button>
          )
        })}
      </div>

      <div className="grid gap-4 lg:grid-cols-[minmax(0,1fr)_380px] lg:items-start">
        {/* viewer */}
        <div className="card overflow-hidden flex flex-col lg:sticky lg:top-16">
          <div className="flex flex-wrap items-center gap-x-4 gap-y-2 px-4 py-2.5 border-b hairline">
            <Switch checked={showGtLigand} onChange={setShowGtLigand} label={<span><span className="inline-block w-2 h-2 rounded-full bg-ok mr-1" />ground truth ligand</span>} />
            <Switch checked={showPredReceptor} onChange={setShowPredReceptor} label="predicted receptor" />
            <Switch checked={showPocket} onChange={setShowPocket} label="pocket residues" />
            <Switch checked={showViolations} onChange={setShowViolations} label="violations" />
            <Switch checked={trajectoryMode} onChange={(v) => { setTrajectoryMode(v); setPlaying(v) }} label="minimisation trajectory" />
            <button className="btn ml-auto" style={{ border: 'none' }} onClick={() => viewer.current?.zoomToLigand()}><Focus size={13} /> ligand</button>
          </div>
          <div className="relative" style={{ height: 520 }}>
            {structures ? (
              <Viewer3D
                ref={viewer}
                gtReceptor={structures.gtReceptor}
                gtLigand={structures.gtLigand}
                predReceptor={structures.predReceptor}
                predLigand={structures.predLigand}
                predTraj={structures.predTraj}
                predColor={color}
                showPredReceptor={showPredReceptor}
                showGtLigand={showGtLigand}
                showPocket={showPocket}
                showViolations={showViolations}
                trajectoryMode={trajectoryMode}
                frame={frame}
                diagnostics={detail?.diagnostics}
                atomDisplacement={detail?.minimisation?.atom_displacement}
                highlightAtoms={highlight}
              />
            ) : <Spinner label="Loading structures" />}
            <div className="absolute left-3 bottom-3 flex items-center gap-3 text-[11px] text-fg-2 bg-panel/85 backdrop-blur px-2.5 py-1.5 rounded-md border hairline">
              <span className="flex items-center gap-1"><span className="w-2 h-2 rounded-sm" style={{ background: color }} />{methodName(index, method)}</span>
              <span className="flex items-center gap-1"><span className="w-2 h-2 rounded-sm bg-ok" />ground truth</span>
              <span className="flex items-center gap-1"><span className="w-2 h-2 rounded-full bg-bad/50" />violation</span>
            </div>
          </div>
          {trajectoryMode && detail?.minimisation && (
            <div className="border-t hairline px-4 py-3 grid gap-3 md:grid-cols-[auto_1fr_auto] items-center">
              <button className="btn" onClick={() => setPlaying((p) => !p)}>{playing ? <Pause size={13} /> : <Play size={13} />}{playing ? 'pause' : 'play'}</button>
              <input type="range" min={0} max={Math.max(0, nFrames - 1)} value={frame} onChange={(e) => { setPlaying(false); setFrame(Number(e.target.value)) }} className="w-full accent-[#5e6ad2]" />
              <span className="mono text-fg-2 whitespace-nowrap">
                frame {frame + 1}/{nFrames} · ΔE {fmt(detail.minimisation.energies[frame] - detail.minimisation.e_local, 1)} kcal/mol · moved {fmt(detail.minimisation.frame_rmsd[frame])} Å
              </span>
            </div>
          )}
        </div>

        {/* side panel */}
        <div className="flex flex-col gap-4 min-w-0">
          {!detail?.ok ? (
            <div className="card px-4 py-3 text-fg-2">No usable prediction for this method{detail?.error ? `: ${detail.error}` : ''}.</div>
          ) : (
            <>
              <div className="card px-4 py-3">
                <Label>Accuracy</Label>
                <div className="grid grid-cols-3 gap-3 mt-2">
                  <Metric label="RMSD" value={`${fmt(detail.rmsd)} Å`} good={detail.rmsd! <= RMSD_SUCCESS} sub="pocket-aligned" />
                  <Metric label="lDDT-PLI" value={fmt(detail.lddt_pli)} sub="benchmark" />
                  <Metric label="pocket Cα" value={`${fmt(detail.superposition?.pocket_ca_rmsd)} Å`} sub={`${detail.superposition?.n_pocket_ca ?? 0} atoms`} />
                </div>
                <div className="text-[11px] text-fg-3 mt-2">
                  seed {String(detail.seed)} · sample {String(detail.sample)} · ranking score {fmt(detail.ranking_score, 3)} · benchmark RMSD {fmt(detail.rmsd_ref)} Å
                </div>
              </div>

              <div className="card px-4 py-3">
                <div className="flex items-center justify-between">
                  <Label>Minimisation strain (MMFF94s)</Label>
                  <span className={`chip ${detail.minimisation && detail.minimisation.strain_local > (sys.gt.minimisation?.strain_local ?? 30) * 1.5 ? 'chip-warn' : 'chip-muted'}`}>
                    GT {fmt(sys.gt.minimisation?.strain_local, 1)} kcal/mol
                  </span>
                </div>
                {detail.minimisation ? (
                  <>
                    <div className="grid grid-cols-3 gap-3 mt-2">
                      <Metric label="local strain" value={`${fmt(detail.minimisation.strain_local, 1)}`} sub="kcal/mol to nearest min" />
                      <Metric label="global strain" value={`${fmt(detail.minimisation.strain_global, 1)}`} sub="kcal/mol to best conf" />
                      <Metric label="drift" value={`${fmt(detail.minimisation.rmsd_drift)} Å`} sub={`max atom ${fmt(detail.minimisation.max_atom_displacement)} Å`} />
                    </div>
                    <div className="mt-2 -mx-2"><Plot data={energyPlot} height={150} layout={{ margin: { l: 40, r: 8, t: 24, b: 28 }, yaxis: { title: { text: 'kcal/mol above min' } }, xaxis: { title: { text: 'frame' } }, legend: { orientation: 'h', y: 1.3, x: 0 } }} /></div>
                  </>
                ) : <div className="text-fg-3 mt-2">not available</div>}
              </div>

              <div className="card px-4 py-3">
                <div className="flex items-center justify-between">
                  <Label>PoseBusters</Label>
                  {detail.pb_pass != null && <span className={`chip ${detail.pb_pass ? 'chip-ok' : 'chip-bad'}`}>{detail.pb_pass ? 'valid' : `${pbFails.length} failed`}</span>}
                </div>
                <div className="grid grid-cols-2 gap-x-3 mt-2">
                  {pbRows.map((k) => {
                    const v = detail.posebusters![k]
                    return (
                      <div key={k} className={`flex items-center gap-1.5 py-[3px] text-[12px] ${v === false ? 'text-bad' : v == null ? 'text-fg-3' : 'text-fg-2'}`}>
                        {v === true ? <Check size={12} className="text-ok" /> : v === false ? <X size={12} /> : <Minus size={12} />}
                        <span className="truncate">{PB_CHECK_LABELS[k] ?? k}</span>
                      </div>
                    )
                  })}
                </div>
              </div>

              {detail.diagnostics && <DiagnosticsPanel d={detail.diagnostics} onHover={setHighlight} />}
            </>
          )}
        </div>
      </div>
    </div>
  )
}

function DiagnosticsPanel({ d, onHover }: { d: NonNullable<MethodDetail['diagnostics']>; onHover: (atoms: number[] | null) => void }) {
  const s = d.summary
  const items: { label: string; count: number; atoms: number[][]; detail: (i: number) => string }[] = [
    { label: 'bond lengths off', count: s.bad_bonds, atoms: d.bonds.filter((b) => b.flag).map((b) => b.atoms), detail: (i) => { const b = d.bonds.filter((x) => x.flag)[i]; return `${b.pred} Å vs ${b.ref} Å in crystal (×${b.ratio})` } },
    { label: 'bond angles off', count: s.bad_angles, atoms: d.angles.map((a) => a.atoms), detail: (i) => `${d.angles[i].pred}° vs ${d.angles[i].ref}° (${d.angles[i].dev > 0 ? '+' : ''}${d.angles[i].dev}°)` },
    { label: 'internal clashes', count: s.intra_clashes, atoms: d.intra_clashes.map((c) => c.atoms), detail: (i) => `${d.intra_clashes[i].dist} Å < ${d.intra_clashes[i].limit} Å` },
    { label: 'protein clashes', count: s.protein_clashes, atoms: d.protein_clashes.map((c) => [c.atom]), detail: (i) => { const c = d.protein_clashes[i]; return `${c.protein.resname}${c.protein.resnum}:${c.protein.atom} at ${c.dist} Å (< ${c.limit})` } },
    { label: 'stereo mismatches', count: s.stereo_mismatches, atoms: d.stereo.filter((x) => x.flag).map((x) => [x.atom]), detail: (i) => { const x = d.stereo.filter((y) => y.flag)[i]; return `${x.pred ?? '?'} predicted, ${x.ref ?? '?'} in crystal` } },
    { label: 'non-planar aromatic rings', count: s.nonplanar_rings, atoms: d.rings.filter((r) => r.flag).map((r) => r.atoms), detail: (i) => `${d.rings.filter((r) => r.flag)[i].max_dev} Å out of plane` },
  ]
  const total = items.reduce((a, b) => a + b.count, 0)
  return (
    <div className="card px-4 py-3">
      <div className="flex items-center justify-between">
        <Label>Where it goes wrong</Label>
        <span className={`chip ${total ? 'chip-warn' : 'chip-ok'}`}>{total ? `${d.flagged_atoms.length} atoms flagged` : 'clean geometry'}</span>
      </div>
      <div className="text-[11px] text-fg-3 mt-1">Hover an item to highlight the atoms in the viewer. Bond/angle deviations are measured against the crystal pose of the same ligand.</div>
      <div className="mt-2 flex flex-col">
        {items.map((it) => (
          <div key={it.label}>
            <div className={`flex items-center justify-between py-1 text-[12px] ${it.count ? 'text-fg' : 'text-fg-3'}`} onMouseEnter={() => onHover(it.atoms.flat())} onMouseLeave={() => onHover(null)}>
              <span>{it.label}</span>
              <span className={`mono ${it.count ? 'text-bad' : ''}`}>{it.count}</span>
            </div>
            {it.count > 0 && (
              <div className="pl-3 pb-1 flex flex-col gap-0.5">
                {it.atoms.slice(0, 8).map((atoms, i) => (
                  <div key={i} className="text-[11px] text-fg-2 mono cursor-default hover:text-accent" onMouseEnter={() => onHover(atoms)} onMouseLeave={() => onHover(null)}>
                    atoms {atoms.map((a) => a + 1).join('–')} · {it.detail(i)}
                  </div>
                ))}
                {it.count > 8 && <div className="text-[11px] text-fg-3">+{it.count - 8} more</div>}
              </div>
            )}
          </div>
        ))}
      </div>
      {d.contacts.length > 0 && (
        <div className="mt-2 pt-2 border-t hairline">
          <div className="text-[11px] text-fg-3 mb-1">Residues within 4.5 Å of the predicted ligand</div>
          <div className="flex flex-wrap gap-1">
            {d.contacts.slice(0, 30).map((c) => <span key={c.residue} className="chip chip-muted mono" title={`${c.min_dist} Å`}>{c.residue}</span>)}
          </div>
        </div>
      )}
    </div>
  )
}

function Metric({ label, value, sub, good }: { label: string; value: string; sub?: string; good?: boolean }) {
  return (
    <div className="min-w-0">
      <div className="text-[11px] text-fg-3">{label}</div>
      <div className={`text-[15px] font-semibold tabular-nums ${good === true ? 'text-ok' : good === false ? 'text-bad' : ''}`}>{value}</div>
      {sub && <div className="text-[10px] text-fg-3 truncate">{sub}</div>}
    </div>
  )
}

function methodColor(index: IndexData | null, id: string) {
  return index?.methods.find((m) => m.id === id)?.color ?? '#5e6ad2'
}
function methodName(index: IndexData | null, id: string) {
  return index?.methods.find((m) => m.id === id)?.name ?? id
}
