import { useCallback, useEffect, useMemo, useRef, useState } from 'react'
import { Link, useParams, useSearchParams } from 'react-router-dom'
import { ArrowLeft, ExternalLink, Pause, Play, Focus, Check, X, Minus } from 'lucide-react'
import type { Data } from 'plotly.js-basic-dist-min'
import { fetchIndex, fetchStructure, fetchSystem, isDocking, PB_CHECK_LABELS, PB_VALIDITY_CHECKS, rcsbUrl, type IndexData, type MethodDetail, type PocketMinimisation, type SystemDetail } from '../lib/api'
import { fmt, RMSD_SUCCESS } from '../lib/stats'
import { Checkbox, ErrorBox, Label, MethodBadge, Select, Spinner, Tip } from '../components/ui'
import Viewer3D, { type HoveredAtom, type ViewerHandle } from '../components/Viewer3D'
import Plot from '../components/Plot'

interface Structures { gtReceptor: string; gtLigand: string; predReceptor: string | null; predLigand: string | null; predTraj: string | null; pocketTraj: string | null; pocketTrajPdb: string | null }
type TrajMode = 'off' | 'ligand' | 'pocket'

export default function SystemView() {
  const { id = '' } = useParams()
  const [params, setParams] = useSearchParams()
  const [index, setIndex] = useState<IndexData | null>(null)
  const [sys, setSys] = useState<SystemDetail | null>(null)
  const [error, setError] = useState<unknown>(null)
  const [structures, setStructures] = useState<Structures | null>(null)
  const [showGtLigand, setShowGtLigand] = useState(true)
  const [showGtReceptor, setShowGtReceptor] = useState(true)
  const [showPredLigand, setShowPredLigand] = useState(true)
  const [showPredReceptor, setShowPredReceptor] = useState(true)
  const [hoveredAtom, setHoveredAtom] = useState<HoveredAtom | null>(null)
  const [showPocket, setShowPocket] = useState(true)
  const [showViolations, setShowViolations] = useState(true)
  const [trajectoryMode, setTrajectoryMode] = useState<TrajMode>('off')
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
      files?.pocket_traj ? fetchStructure(id, files.pocket_traj).catch(() => null) : Promise.resolve(null),
      files?.pocket_traj_pdb ? fetchStructure(id, files.pocket_traj_pdb).catch(() => null) : Promise.resolve(null),
    ]).then(([gtReceptor, gtLigand, predReceptor, predLigand, predTraj, pocketTraj, pocketTrajPdb]) => {
      if (!cancelled) setStructures({ gtReceptor, gtLigand, predReceptor, predLigand, predTraj, pocketTraj, pocketTrajPdb })
    }).catch(setError)
    return () => { cancelled = true }
  }, [sys, method, id])

  // trajectory playback
  const pocket = detail?.pocket_minimisation?.ok ? detail.pocket_minimisation : null
  const hasPocketTraj = !!pocket && !!detail?.files?.pocket_traj
  const activeRun = trajectoryMode === 'pocket' ? pocket : trajectoryMode === 'ligand' ? detail?.minimisation ?? null : null
  const nFrames = activeRun?.n_frames ?? 0
  useEffect(() => {
    if (!playing || trajectoryMode === 'off' || nFrames < 2) return
    const t = setInterval(() => setFrame((f) => (f + 1) % nFrames), 120)
    return () => clearInterval(t)
  }, [playing, trajectoryMode, nFrames])
  useEffect(() => { setFrame(0); setPlaying(false); setHighlight(null); setHoveredAtom(null) }, [method])
  const onAtomHover = useCallback((a: HoveredAtom | null) => setHoveredAtom(a), [])
  useEffect(() => { setFrame(0); setPlaying(trajectoryMode !== 'off') }, [trajectoryMode])

  const energyPlot = useMemo<Data[]>(() => {
    if (!detail?.minimisation) return []
    const m = detail.minimisation
    const gt = sys?.gt.minimisation
    const rel = (e: number[]) => e.map((x) => x - e[e.length - 1])
    const out: Data[] = [{ type: 'scatter', mode: 'lines', name: 'prediction', x: m.energies.map((_, i) => i), y: rel(m.energies), line: { color: methodColor(index, method), width: 2 }, hovertemplate: 'frame %{x}: +%{y:.1f} kcal/mol<extra></extra>' }]
    if (gt) out.push({ type: 'scatter', mode: 'lines', name: 'ground truth', x: gt.energies.map((_, i) => i), y: rel(gt.energies), line: { color: '#2f9e6b', width: 1.5, dash: 'dot' }, hovertemplate: 'frame %{x}: +%{y:.1f} kcal/mol<extra></extra>' })
    if (trajectoryMode === 'ligand') out.push({ type: 'scatter', mode: 'markers', name: 'current', x: [frame], y: [rel(m.energies)[frame]], marker: { color: '#1c1c22', size: 8 }, showlegend: false, hoverinfo: 'skip' })
    return out
  }, [detail, sys, frame, trajectoryMode, index, method])

  const pocketPlot = useMemo<Data[]>(() => {
    if (!pocket) return []
    const gt = sys?.gt.pocket_minimisation?.ok ? sys.gt.pocket_minimisation : null
    const rel = (e: number[]) => e.map((x) => x - e[e.length - 1])
    const out: Data[] = [{ type: 'scatter', mode: 'lines', name: 'prediction', x: pocket.energies.map((_, i) => i), y: rel(pocket.energies), line: { color: methodColor(index, method), width: 2 }, hovertemplate: 'frame %{x}: +%{y:.1f} kcal/mol<extra></extra>' }]
    if (gt) out.push({ type: 'scatter', mode: 'lines', name: 'ground truth', x: gt.energies.map((_, i) => i), y: rel(gt.energies), line: { color: '#2f9e6b', width: 1.5, dash: 'dot' }, hovertemplate: 'frame %{x}: +%{y:.1f} kcal/mol<extra></extra>' })
    if (trajectoryMode === 'pocket') out.push({ type: 'scatter', mode: 'markers', name: 'current', x: [frame], y: [rel(pocket.energies)[frame]], marker: { color: '#1c1c22', size: 8 }, showlegend: false, hoverinfo: 'skip' })
    return out
  }, [pocket, sys, frame, trajectoryMode, index, method])

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
              <MethodBadge method={m} />
              {d?.ok && <span className={`mono ${d.rmsd! <= RMSD_SUCCESS ? 'text-ok' : 'text-bad'}`}>{fmt(d.rmsd)} Å</span>}
              {d?.ok && d.pb_pass === false && <span className="w-1.5 h-1.5 rounded-full bg-warn" title="PoseBusters violations" />}
            </button>
          )
        })}
      </div>

      {isDocking(index.methods.find((m) => m.id === method)) && (
        <div className="text-[12px] text-fg-2 bg-warn-2 border border-[#f3dfb5] rounded-md px-3 py-2">
          Rigid holo redocking baseline: the receptor is the experimental structure, so only the ligand pose is predicted.
        </div>
      )}

      <div className="grid gap-4 lg:grid-cols-[minmax(0,1fr)_380px] lg:items-start">
        {/* viewer */}
        <div className="flex flex-col gap-4 min-w-0">
        <div className="card overflow-hidden flex flex-col">
          <div className="relative" style={{ height: 520 }}>
            {structures ? (
              <Viewer3D
                ref={viewer}
                gtReceptor={structures.gtReceptor}
                gtLigand={structures.gtLigand}
                predReceptor={structures.predReceptor}
                predLigand={structures.predLigand}
                predTraj={structures.predTraj}
                pocketTraj={structures.pocketTraj}
                pocketTrajPdb={structures.pocketTrajPdb}
                predColor={color}
                showGtReceptor={showGtReceptor}
                showGtLigand={showGtLigand}
                showPredLigand={showPredLigand}
                showPredReceptor={showPredReceptor}
                showPocket={showPocket}
                showViolations={showViolations}
                trajectoryMode={trajectoryMode}
                frame={frame}
                diagnostics={detail?.diagnostics}
                atomDisplacement={trajectoryMode === 'pocket' ? pocket?.ligand_atom_displacement : detail?.minimisation?.atom_displacement}
                highlightAtoms={highlight}
                onAtomHover={onAtomHover}
              />
            ) : <Spinner label="Loading structures" />}
            <div className="absolute left-3 bottom-3 flex items-center gap-3 text-[11px] text-fg-2 bg-panel/85 backdrop-blur px-2.5 py-1.5 rounded-md border hairline">
              <span className="flex items-center gap-1"><span className="w-2 h-2 rounded-sm" style={{ background: color }} />{methodName(index, method)}</span>
              <span className="flex items-center gap-1"><span className="w-2 h-2 rounded-sm bg-ok" />ground truth</span>
              <span className="flex items-center gap-1"><span className="w-2 h-2 rounded-full bg-bad/50" />violation</span>
              {hoveredAtom && <span className="mono text-fg">{hoveredAtom.elem}{hoveredAtom.index + 1}</span>}
            </div>
          </div>
          {/* control strip */}
          <div className="border-t hairline px-4 py-2.5 flex flex-wrap items-center gap-x-6 gap-y-2">
            <div className="grid grid-cols-[auto_auto_auto] gap-x-3 gap-y-1 items-center">
              <span className="text-[11px] text-fg-3 uppercase tracking-wider">Ground truth</span>
              <Checkbox checked={showGtLigand} onChange={setShowGtLigand} label="ligand" dot="#2f9e6b" />
              <Checkbox checked={showGtReceptor} onChange={setShowGtReceptor} label="receptor" dot="#d9d9de" />
              <span className="text-[11px] text-fg-3 uppercase tracking-wider">Predicted</span>
              <Checkbox checked={showPredLigand} onChange={setShowPredLigand} label="ligand" dot={color} />
              <Checkbox checked={showPredReceptor} onChange={setShowPredReceptor} label="receptor" dot={color} />
            </div>
            <div className="flex flex-col gap-1">
              <Checkbox checked={showPocket} onChange={setShowPocket} label="pocket residues" />
              <Checkbox checked={showViolations} onChange={setShowViolations} label="violations" />
            </div>
            <span className="flex items-center gap-2 text-[12px] text-fg-2">
              trajectory
              <Select<TrajMode>
                value={trajectoryMode}
                onChange={setTrajectoryMode}
                width={150}
                options={[{ value: 'off', label: 'off' }, { value: 'ligand', label: 'ligand only' }, ...(hasPocketTraj ? [{ value: 'pocket' as TrajMode, label: 'ligand + pocket' }] : [])]}
              />
            </span>
            <button className="btn ml-auto" style={{ border: 'none' }} onClick={() => viewer.current?.zoomToLigand()}><Focus size={13} /> ligand</button>
          </div>
          {trajectoryMode !== 'off' && activeRun && (
            <div className="border-t hairline px-4 py-3 grid gap-3 md:grid-cols-[auto_1fr_auto] items-center">
              <button className="btn" onClick={() => setPlaying((p) => !p)}>{playing ? <Pause size={13} /> : <Play size={13} />}{playing ? 'pause' : 'play'}</button>
              <input type="range" min={0} max={Math.max(0, nFrames - 1)} value={frame} onChange={(e) => { setPlaying(false); setFrame(Number(e.target.value)) }} className="w-full accent-[#5e6ad2]" />
              <span className="mono text-fg-2 whitespace-nowrap">
                {trajectoryMode === 'pocket' ? 'ligand + pocket' : 'ligand only'} · frame {frame + 1}/{nFrames} · ΔE {fmt(activeRun.energies[frame] - activeRun.energies[activeRun.energies.length - 1], 1)} kcal/mol · ligand moved {fmt(('frame_ligand_rmsd' in activeRun ? activeRun.frame_ligand_rmsd : activeRun.frame_rmsd)[frame])} Å
              </span>
            </div>
          )}
        </div>
        {detail?.ok && detail.diagnostics && <DiagnosticsPanel d={detail.diagnostics} onHover={setHighlight} hoveredAtom={hoveredAtom?.index ?? null} />}
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

              <PocketCard pocket={detail.pocket_minimisation ?? null} gt={sys.gt.pocket_minimisation ?? null} plot={pocketPlot} showPlot={trajectoryMode === 'pocket'} />

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
                      <Metric label="local strain" value={`${fmt(detail.minimisation.strain_local, 1)}`} sub="kcal/mol to nearest min · lower is better" />
                      <Metric label="global strain" value={`${fmt(detail.minimisation.strain_global, 1)}`} sub="kcal/mol to best conf · lower is better" />
                      <Metric label="drift" value={`${fmt(detail.minimisation.rmsd_drift)} Å`} sub={`max atom ${fmt(detail.minimisation.max_atom_displacement)} Å`} />
                    </div>
                    <div className="text-[11px] text-fg-3 mt-2">Strain: energy the pose must release to reach a minimum; lower is better.</div>
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

            </>
          )}
        </div>
      </div>
    </div>
  )
}

function PocketCard({ pocket, gt, plot, showPlot }: { pocket: PocketMinimisation | null; gt: PocketMinimisation | null; plot: Data[]; showPlot: boolean }) {
  const ok = !!pocket?.ok
  const gtOk = !!gt?.ok
  const clashing = ok && pocket!.e_interaction_pose > 0
  return (
    <div className="card px-4 py-3">
      <div className="flex items-center justify-between gap-2">
        <Label>Pocket relaxation (MMFF94s, protein pocket restrained)</Label>
        {gtOk && <span className="chip chip-muted whitespace-nowrap" title="ground-truth complex, same protocol">GT {fmt(gt!.e_interaction_pose, 1)} → {fmt(gt!.e_interaction_min, 1)} kcal/mol</span>}
      </div>
      {!ok ? (
        <div className="text-fg-3 mt-2">not available{pocket?.error ? `: ${pocket.error}` : ''}</div>
      ) : (
        <>
          <div className="flex items-center gap-2 mt-2 flex-wrap">
            <span className={`chip ${clashing ? 'chip-bad' : 'chip-ok'}`}>{clashing ? 'clashing' : 'favourable'} interaction at pose</span>
            <span className={`chip ${pocket!.clashes_pose > 0 ? 'chip-bad' : 'chip-muted'}`}>{pocket!.clashes_pose} clash{pocket!.clashes_pose === 1 ? '' : 'es'} → {pocket!.clashes_min} after relaxation</span>
          </div>
          <div className="grid grid-cols-3 gap-3 mt-3">
            <Metric label="interaction E" value={`${fmt(pocket!.e_interaction_pose, 1)} → ${fmt(pocket!.e_interaction_min, 1)}`} sub="kcal/mol pose → min · lower is better" good={clashing ? false : undefined} />
            <Metric label="ligand drift" value={`${fmt(pocket!.ligand_rmsd_drift)} Å`} sub={`max atom ${fmt(Math.max(0, ...pocket!.ligand_atom_displacement))} Å`} />
            <Metric label="pocket RMSD" value={`${fmt(pocket!.pocket_heavy_rmsd)} Å`} sub={`${pocket!.n_pocket_residues} residues · max ${fmt(pocket!.max_pocket_atom_displacement)} Å`} />
          </div>
          <div className="text-[11px] text-fg-3 mt-2">Interaction energy: lower is better; &gt; 0 means the pose is repulsive (clashing) with the pocket.</div>
          {showPlot && (
            <div className="mt-2 -mx-2"><Plot data={plot} height={150} layout={{ margin: { l: 40, r: 8, t: 24, b: 28 }, yaxis: { title: { text: 'kcal/mol above min' } }, xaxis: { title: { text: 'frame' } }, legend: { orientation: 'h', y: 1.3, x: 0 } }} /></div>
          )}
        </>
      )}
    </div>
  )
}

function DiagnosticsPanel({ d, onHover, hoveredAtom }: { d: NonNullable<MethodDetail['diagnostics']>; onHover: (atoms: number[] | null) => void; hoveredAtom: number | null }) {
  const s = d.summary
  const panel = useRef<HTMLDivElement>(null)
  // reverse hover: scroll the first sidebar entry containing the hovered viewer atom into view
  useEffect(() => {
    if (hoveredAtom == null || !panel.current) return
    const first = panel.current.querySelector<HTMLElement>('[data-hit="true"]')
    first?.scrollIntoView({ block: 'nearest', behavior: 'smooth' })
  }, [hoveredAtom])
  const hit = (atoms: number[]) => hoveredAtom != null && atoms.includes(hoveredAtom)
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
    <div className="card px-4 py-3" ref={panel}>
      <div className="flex items-center justify-between">
        <Label>Physical violations</Label>
        <span className={`chip ${total ? 'chip-warn' : 'chip-ok'}`}>{total ? `${d.flagged_atoms.length} atom${d.flagged_atoms.length === 1 ? "" : "s"} flagged` : 'clean geometry'}</span>
      </div>
      <div className="text-[11px] text-fg-3 mt-1">Hover an item to highlight its atoms in the viewer, or hover a ligand atom in the viewer to find it here. Bond/angle deviations are measured against the crystal pose of the same ligand.</div>
      <div className="mt-2 grid gap-x-6 md:grid-cols-2">
        {items.map((it) => (
          <div key={it.label}>
            <div className={`flex items-center justify-between py-1 text-[12px] ${it.count ? 'text-fg' : 'text-fg-3'}`} onMouseEnter={() => onHover(it.atoms.flat())} onMouseLeave={() => onHover(null)}>
              <span>{it.label}</span>
              <span className={`mono ${it.count ? 'text-bad' : ''}`}>{it.count}</span>
            </div>
            {it.count > 0 && (
              <div className={`pl-3 pb-1 flex flex-col gap-0.5 ${it.count > 12 ? 'overflow-y-auto pr-1' : ''}`} style={it.count > 12 ? { maxHeight: 260 } : undefined}>
                {it.atoms.map((atoms, i) => (
                  <div key={i} data-hit={hit(atoms) ? 'true' : undefined} className={`text-[11px] mono cursor-default hover:text-accent rounded px-1 -mx-1 border-l-2 ${hit(atoms) ? 'bg-warn-2 text-fg border-warn' : 'text-fg-2 border-transparent'}`} onMouseEnter={() => onHover(atoms)} onMouseLeave={() => onHover(null)}>
                    atoms {atoms.map((a) => a + 1).join('–')} · {it.detail(i)}
                  </div>
                ))}
              </div>
            )}
          </div>
        ))}
      </div>
      {d.contacts.length > 0 && (
        <div className="mt-2 pt-2 border-t hairline">
          <div className="text-[11px] text-fg-3 mb-1">Residues within 4.5 Å of the predicted ligand</div>
          <div className="flex flex-wrap gap-1">
            {d.contacts.map((c) => <span key={c.residue} className="chip chip-muted mono" title={`${c.min_dist} Å`}>{c.residue}</span>)}
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
