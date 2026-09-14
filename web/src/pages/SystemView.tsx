import { useCallback, useEffect, useMemo, useRef, useState } from 'react'
import { Link, useParams, useSearchParams } from 'react-router-dom'
import { ArrowLeft, ExternalLink, Pause, Play, Focus, Check, X, Minus } from 'lucide-react'
import type { Data } from 'plotly.js-basic-dist-min'
import { fetchIndex, fetchStructure, fetchSystem, isDocking, PB_CHECK_LABELS, PB_VALIDITY_CHECKS, rcsbUrl, type Contacts, type IndexData, type MethodDetail, type PocketMinimisation, type SystemDetail } from '../lib/api'
import { fmt, RMSD_SUCCESS } from '../lib/stats'
import { Checkbox, ErrorBox, Label, MethodBadge, Select, Spinner, Tip } from '../components/ui'
import Viewer3D, { type ClashRef, type HoveredAtom, type ViewerHandle } from '../components/Viewer3D'
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
  const [hoveredClash, setHoveredClash] = useState<ClashRef | null>(null)   // viewer -> panel
  const [highlightClash, setHighlightClash] = useState<ClashRef | null>(null) // panel -> viewer
  const hoverClearTimer = useRef<number | null>(null)
  const [highlightResidue, setHighlightResidue] = useState<string | null>(null)
  const [highlightGtAtoms, setHighlightGtAtoms] = useState<number[] | null>(null)
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
  // viewer atom hover: set immediately, clear with a short debounce so moving between atoms does not flicker
  const onAtomHover = useCallback((a: HoveredAtom | null) => {
    if (hoverClearTimer.current) { window.clearTimeout(hoverClearTimer.current); hoverClearTimer.current = null }
    if (a) setHoveredAtom(a)
    else hoverClearTimer.current = window.setTimeout(() => { hoverClearTimer.current = null; setHoveredAtom(null) }, 150)
  }, [])
  // violations that involve the atom hovered in the viewer -> same highlight the panel hover produces
  const viewerHighlight = useMemo(() => {
    const d = detail?.ok ? detail.diagnostics : null
    if (!d || hoveredAtom == null) return null
    const i = hoveredAtom.index
    const atoms = new Set<number>([i])
    const clashes: ClashRef[] = []
    for (const b of d.bonds) if (b.flag && b.atoms.includes(i)) b.atoms.forEach((x) => atoms.add(x))
    for (const a of d.angles) if (a.atoms.includes(i)) a.atoms.forEach((x) => atoms.add(x))
    for (const r of d.rings) if (r.flag && r.atoms.includes(i)) r.atoms.forEach((x) => atoms.add(x))
    for (const s of d.stereo) if (s.flag && s.atom === i) atoms.add(i)
    d.intra_clashes.forEach((c, k) => { if (c.atoms.includes(i)) { c.atoms.forEach((x) => atoms.add(x)); clashes.push({ kind: 'intra', index: k }) } })
    d.protein_clashes.forEach((c, k) => { if (c.atom === i) clashes.push({ kind: 'protein', index: k }) })
    // a clean atom (not part of any violation) produces no highlight at all
    if (atoms.size === 1 && clashes.length === 0 && !d.flagged_atoms.includes(i)) return null
    return { atoms: [...atoms], clashes }
  }, [detail, hoveredAtom])
  const effHighlightAtoms = highlight ?? viewerHighlight?.atoms ?? null
  const effHighlightClash: ClashRef | ClashRef[] | null = highlight ? highlightClash : (viewerHighlight?.clashes.length ? viewerHighlight.clashes : null)
  const onClashHover = useCallback((c: ClashRef | null) => setHoveredClash(c), [])
  useEffect(() => { setHighlightClash(null); setHoveredClash(null); setHighlightResidue(null); setHighlightGtAtoms(null) }, [method])
  const onContactHover = useCallback((residue: string | null, gtAtoms: number[] | null) => { setHighlightResidue(residue); setHighlightGtAtoms(gtAtoms) }, [])
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
                highlightAtoms={effHighlightAtoms}
                highlightClash={effHighlightClash}
                highlightResidue={highlightResidue}
                highlightGtAtoms={highlightGtAtoms}
                onAtomHover={onAtomHover}
                onClashHover={onClashHover}
              />
            ) : <Spinner label="Loading structures" />}
            <div className="absolute left-3 bottom-3 flex items-center gap-3 text-[11px] text-fg-2 bg-panel/85 backdrop-blur px-2.5 py-1.5 rounded-md border hairline">
              <span className="flex items-center gap-1"><span className="w-2 h-2 rounded-sm" style={{ background: color }} />{methodName(index, method)}</span>
              <span className="flex items-center gap-1"><span className="w-2 h-2 rounded-sm bg-ok" />ground truth</span>
              <span className="flex items-center gap-1"><span className="w-2 h-2 rounded-full bg-bad/50" />violation</span>
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
                options={[{ value: 'off', label: 'off' }, ...(hasPocketTraj ? [{ value: 'pocket' as TrajMode, label: 'ligand + pocket' }] : []), { value: 'ligand', label: 'ligand only' }]}
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
        {detail?.ok && detail.diagnostics && <DiagnosticsPanel d={detail.diagnostics} contacts={detail.contacts ?? null} onHover={setHighlight} onHoverClash={setHighlightClash} onHoverContact={onContactHover} hoveredAtom={highlight || !viewerHighlight ? null : hoveredAtom?.index ?? null} hoveredClash={highlight ? null : hoveredClash} />}
        </div>

        {/* side panel */}
        <div className="flex flex-col gap-4 min-w-0">
          {!detail?.ok ? (
            <div className="card px-4 py-3 text-fg-2">No usable prediction for this method{detail?.error ? `: ${detail.error}` : ''}.</div>
          ) : (
            <>
              <div className="card px-4 py-3">
                <div className="flex items-center justify-between">
                  <Label>Accuracy</Label>
                  {detail.pocket_hit && !detail.pocket_hit.hit && <span className="chip chip-bad" title={`shape overlap ${fmt(detail.pocket_hit.shape_overlap)} · centroid ${fmt(detail.pocket_hit.centroid_distance, 1)} Å from the crystal ligand`}>missed the pocket</span>}
                </div>
                <div className="grid grid-cols-3 gap-3 mt-2">
                  <Metric label="RMSD" value={`${fmt(detail.rmsd)} Å`} good={detail.rmsd! <= RMSD_SUCCESS} sub="pocket-aligned" />
                  <Metric label="lDDT-PLI" value={fmt(detail.lddt_pli)} sub="benchmark" />
                  <Metric label="pocket Cα" value={`${fmt(detail.superposition?.pocket_ca_rmsd)} Å`} sub={`${detail.superposition?.n_pocket_ca ?? 0} atoms`} />
                  {detail.contacts && (
                    <Metric label="contacts kept" value={`${detail.contacts.kept} / ${detail.contacts.gt_total}`} sub={detail.contacts.retention == null ? 'no crystal contacts' : `${Math.round(detail.contacts.retention * 100)}% of crystal contacts`} good={detail.contacts.retention == null ? undefined : detail.contacts.retention >= 0.5} />
                  )}
                  {detail.pocket_hit && (
                    <Metric label="shape overlap" value={fmt(detail.pocket_hit.shape_overlap)} sub={`centroid ${fmt(detail.pocket_hit.centroid_distance, 1)} Å away`} good={detail.pocket_hit.hit ? undefined : false} />
                  )}
                </div>
                {detail.contacts && (
                  <div className="text-[11px] text-fg-3 mt-2">
                    H-bonds {detail.contacts.by_type.hbond?.kept ?? 0}/{detail.contacts.by_type.hbond?.gt ?? 0} · hydrophobic {detail.contacts.by_type.hydrophobic?.kept ?? 0}/{detail.contacts.by_type.hydrophobic?.gt ?? 0} · ionic {detail.contacts.by_type.ionic?.kept ?? 0}/{detail.contacts.by_type.ionic?.gt ?? 0}
                  </div>
                )}
                <div className="text-[11px] text-fg-3 mt-2">
                  seed {String(detail.seed)} · sample {String(detail.sample)} · ranking score {fmt(detail.ranking_score, 3)} · benchmark RMSD {fmt(detail.rmsd_ref)} Å
                </div>
              </div>

              <PocketCard pocket={detail.pocket_minimisation ?? null} gt={sys.gt.pocket_minimisation ?? null} plot={pocketPlot} showPlot={true} excess={excessFor(index, sys.system_id, method, 'excess_relaxation_de') ?? (detail.pocket_minimisation?.ok && sys.gt.pocket_minimisation?.ok ? (detail.pocket_minimisation.e_interaction_pose - detail.pocket_minimisation.e_interaction_min) - (sys.gt.pocket_minimisation.e_interaction_pose - sys.gt.pocket_minimisation.e_interaction_min) : null)} />

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
                    <div className="mt-2"><ExcessChip value={excessFor(index, sys.system_id, method, 'excess_strain_local') ?? (sys.gt.minimisation ? detail.minimisation.strain_local - sys.gt.minimisation.strain_local : null)} label="local strain, excess vs crystal" /></div>
                    <div className="text-[11px] text-fg-3 mt-2">Strain: energy the pose must release to reach a minimum; lower is better. If the crystal pose itself is strained under this force field, the model cannot be expected to do better; the excess is the part the model added.</div>
                    <div className="mt-2 -mx-2"><Plot data={energyPlot} height={150} layout={{ margin: { l: 44, r: 8, t: 24, b: 28 }, yaxis: { title: { text: 'ΔE vs. minimum (kcal/mol)' } }, xaxis: { title: { text: 'frame' } }, legend: { orientation: 'h', y: 1.3, x: 0 } }} /></div>
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

/** excess (prediction − crystal) value from index.json result rows when the preprocessing provides it */
function excessFor(index: IndexData | null, systemId: string, method: string, key: 'excess_relaxation_de' | 'excess_strain_local'): number | null {
  const r = index?.results.find((x) => x.system_id === systemId && x.method === method)
  const v = r?.[key]
  return v == null ? null : v
}

function ExcessChip({ value, label }: { value: number | null; label: string }) {
  if (value == null) return null
  const cls = value <= 5 ? 'chip-muted' : value <= 25 ? 'chip-warn' : 'chip-bad'  // small excess is neutral, never 'good'
  return <span className={`chip ${cls} self-start`} title="prediction − crystal pose, same force field">{label} {value >= 0 ? '+' : ''}{value.toFixed(1)} kcal/mol</span>
}

function PocketCard({ pocket, gt, plot, showPlot, excess }: { pocket: PocketMinimisation | null; gt: PocketMinimisation | null; plot: Data[]; showPlot: boolean; excess: number | null }) {
  const ok = !!pocket?.ok
  const gtOk = !!gt?.ok
  const dE = ok ? pocket!.e_interaction_pose - pocket!.e_interaction_min : null
  const gtDE = gtOk ? gt!.e_interaction_pose - gt!.e_interaction_min : null
  const clashing = ok && pocket!.clashes_pose > 0
  const strained = ok && !clashing && dE != null && dE > 50
  return (
    <div className="card px-4 py-3">
      <div className="flex items-center justify-between gap-2">
        <Label>Pocket relaxation (MMFF94s, protein pocket restrained)</Label>
        {gtOk && <span className="chip chip-muted whitespace-nowrap" title="ground-truth complex, same protocol">GT ΔE {fmt(gtDE, 1)} kcal/mol</span>}
      </div>
      {!ok ? (
        <div className="text-fg-3 mt-2">not available{pocket?.error ? `: ${pocket.error}` : ''}</div>
      ) : (
        <>
          <div className="flex items-center gap-2 mt-2 flex-wrap">
            <span className={`chip ${clashing ? 'chip-bad' : strained ? 'chip-warn' : 'chip-ok'}`}>{clashing ? 'clashing pose' : strained ? 'strained pose' : 'relaxed pose'}</span>
            <span className={`chip ${pocket!.clashes_pose > 0 ? 'chip-bad' : 'chip-muted'}`}>{pocket!.clashes_pose} clash{pocket!.clashes_pose === 1 ? '' : 'es'} → {pocket!.clashes_min} after relaxation</span>
          </div>
          <div className="grid grid-cols-3 gap-3 mt-3">
            <Metric label="relaxation ΔE" value={`${fmt(dE, 1)}`} sub="kcal/mol released · lower is better" good={clashing || strained ? false : undefined} />
            <Metric label="ligand drift" value={`${fmt(pocket!.ligand_rmsd_drift)} Å`} sub={`max atom ${fmt(Math.max(0, ...pocket!.ligand_atom_displacement))} Å`} />
            <Metric label="pocket RMSD" value={`${fmt(pocket!.pocket_heavy_rmsd)} Å`} sub={`${pocket!.n_pocket_residues} residues · max ${fmt(pocket!.max_pocket_atom_displacement)} Å`} />
          </div>
          <div className="mt-2"><ExcessChip value={excess} label="relaxation ΔE, excess vs crystal" /></div>
          <div className="text-[11px] text-fg-3 mt-2">Relaxation ΔE is the energy released when the pocket and ligand are allowed to relax; large values mean the pose was strained or clashing. If the crystal pose itself is strained under this force field, the model cannot be expected to do better; the excess is the part the model added.</div>
          {showPlot && (
            <div className="mt-2 -mx-2"><Plot data={plot} height={150} layout={{ margin: { l: 44, r: 8, t: 24, b: 28 }, yaxis: { title: { text: 'ΔE vs. minimum (kcal/mol)' } }, xaxis: { title: { text: 'frame' } }, legend: { orientation: 'h', y: 1.3, x: 0 } }} /></div>
          )}
        </>
      )}
    </div>
  )
}

function DiagnosticsPanel({ d, contacts, onHover, onHoverClash, onHoverContact, hoveredAtom, hoveredClash }: { d: NonNullable<MethodDetail['diagnostics']>; contacts: Contacts | null; onHover: (atoms: number[] | null) => void; onHoverClash: (c: ClashRef | null) => void; onHoverContact: (residue: string | null, gtAtoms: number[] | null) => void; hoveredAtom: number | null; hoveredClash: ClashRef | null }) {
  const s = d.summary
  const panel = useRef<HTMLDivElement>(null)
  // reverse hover: scroll the first sidebar entry containing the hovered viewer atom into view
  useEffect(() => {
    if ((hoveredAtom == null && hoveredClash == null) || !panel.current) return
    const first = panel.current.querySelector<HTMLElement>('[data-hit="true"]')
    first?.scrollIntoView({ block: 'nearest' })
  }, [hoveredAtom, hoveredClash])
  const hit = (atoms: number[], clash?: ClashRef) =>
    (hoveredClash != null && clash != null && hoveredClash.kind === clash.kind && hoveredClash.index === clash.index)
    || (hoveredClash == null && hoveredAtom != null && atoms.includes(hoveredAtom))
  const items: { label: string; count: number; atoms: number[][]; clashKind?: ClashRef['kind']; detail: (i: number) => string }[] = [
    { label: 'bond lengths off', count: s.bad_bonds, atoms: d.bonds.filter((b) => b.flag).map((b) => b.atoms), detail: (i) => { const b = d.bonds.filter((x) => x.flag)[i]; return `${b.pred} Å vs ${b.ref} Å in crystal (×${b.ratio})` } },
    { label: 'bond angles off', count: s.bad_angles, atoms: d.angles.map((a) => a.atoms), detail: (i) => `${d.angles[i].pred}° vs ${d.angles[i].ref}° (${d.angles[i].dev > 0 ? '+' : ''}${d.angles[i].dev}°)` },
    { label: 'internal clashes', count: s.intra_clashes, atoms: d.intra_clashes.map((c) => c.atoms), clashKind: 'intra', detail: (i) => `${d.intra_clashes[i].dist} Å < ${d.intra_clashes[i].limit} Å` },
    { label: 'protein clashes', count: s.protein_clashes, atoms: d.protein_clashes.map((c) => [c.atom]), clashKind: 'protein', detail: (i) => { const c = d.protein_clashes[i]; return `${c.protein.resname}${c.protein.resnum}:${c.protein.atom} at ${c.dist} Å (< ${c.limit})` } },
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
                {it.atoms.map((atoms, i) => {
                  const clash: ClashRef | undefined = it.clashKind ? { kind: it.clashKind, index: i } : undefined
                  const on = hit(atoms, clash)
                  return (
                  <div key={i} data-hit={on ? 'true' : undefined} className={`text-[11px] mono cursor-default hover:text-accent rounded px-1 -mx-1 border-l-2 ${on ? 'bg-warn-2 text-fg border-warn' : 'text-fg-2 border-transparent'}`}
                    onMouseEnter={() => { onHover(atoms); if (clash) onHoverClash(clash) }} onMouseLeave={() => { onHover(null); if (clash) onHoverClash(null) }}>
                    atoms {atoms.map((a) => a + 1).join('–')} · {it.detail(i)}
                  </div>
                  )
                })}
              </div>
            )}
          </div>
        ))}
      </div>
      {contacts && contacts.residues.length > 0 ? (
        <div className="mt-2 pt-2 border-t hairline">
          <div className="text-[11px] text-fg-3 mb-1 flex flex-wrap items-center gap-x-3 gap-y-1">
            <span>Contact residues (≤ {contacts.cutoff} Å) · hover to show the side chain</span>
            <span className="inline-flex items-center gap-1"><span className="chip chip-ok" style={{ height: 14, padding: '0 5px', fontSize: 10 }}>kept</span> in crystal and prediction</span>
            <span className="inline-flex items-center gap-1"><span className="chip chip-muted line-through" style={{ height: 14, padding: '0 5px', fontSize: 10 }}>lost</span> crystal only</span>
            <span className="inline-flex items-center gap-1"><span className="chip" style={{ height: 14, padding: '0 5px', fontSize: 10, border: '1px solid var(--color-line-2)', color: 'var(--color-fg-2)' }}>new</span> prediction only</span>
          </div>
          <div className="flex flex-wrap gap-1">
            {contacts.residues.map((c) => {
              const kind = c.gt && c.pred ? 'kept' : c.gt ? 'lost' : 'new'
              const lostAtoms = kind === 'lost' ? contacts.lost.filter((l) => l.residue === c.residue).map((l) => l.ligand_atom) : null
              const cls = kind === 'kept' ? 'chip chip-ok mono' : kind === 'lost' ? 'chip chip-muted mono line-through opacity-80' : 'chip mono'
              const style = kind === 'new' ? { border: '1px solid var(--color-line-2)', color: 'var(--color-fg-2)' } : undefined
              return (
                <span key={c.residue} className={`${cls} cursor-default hover:ring-2 hover:ring-accent-2`} style={style} title={kind}
                  onMouseEnter={() => onHoverContact(c.residue, lostAtoms)} onMouseLeave={() => onHoverContact(null, null)}>
                  {c.residue}
                </span>
              )
            })}
          </div>
        </div>
      ) : d.contacts.length > 0 && (
        <div className="mt-2 pt-2 border-t hairline">
          <div className="text-[11px] text-fg-3 mb-1">Residues within 4.5 Å of the predicted ligand · hover to show the side chain</div>
          <div className="flex flex-wrap gap-1">
            {d.contacts.map((c) => <span key={c.residue} className="chip chip-muted mono cursor-default" title={`${c.min_dist} Å`} onMouseEnter={() => onHoverContact(c.residue, null)} onMouseLeave={() => onHoverContact(null, null)}>{c.residue}</span>)}
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
