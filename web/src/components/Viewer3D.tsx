import { useEffect, useImperativeHandle, useRef, useState, forwardRef } from 'react'
import * as $3Dmol from '3dmol'
import type { Diagnostics } from '../lib/api'

export interface ViewerHandle {
  zoomToLigand: () => void
  setFrame: (i: number) => void
}

export interface HoveredAtom { index: number; elem: string }
/** identifies one clash entry of the diagnostics: protein_clashes[i] or intra_clashes[i] */
export interface ClashRef { kind: 'protein' | 'intra'; index: number }

export interface ViewerProps {
  gtReceptor: string
  gtLigand: string
  predReceptor?: string | null
  predLigand?: string | null
  predTraj?: string | null
  pocketTraj?: string | null
  pocketTrajPdb?: string | null
  predColor: string
  showGtReceptor: boolean
  showGtLigand: boolean
  showPredLigand: boolean
  showPredReceptor: boolean
  showPocket: boolean
  showViolations: boolean
  trajectoryMode: 'off' | 'ligand' | 'pocket'
  frame: number
  diagnostics?: Diagnostics | null
  atomDisplacement?: number[] | null
  highlightAtoms?: number[] | null
  /** clash entries to emphasise (panel hover or viewer-atom hover); everything else is dimmed */
  highlightClash?: ClashRef | ClashRef[] | null
  /** residue label like "A:PHE46" whose side chain is drawn in the accent colour (contact chip hover) */
  highlightResidue?: string | null
  /** crystal-ligand atom indices to emphasise (lost contacts) */
  highlightGtAtoms?: number[] | null
  /** fired (throttled) when the mouse enters / leaves a predicted-ligand atom */
  onAtomHover?: (atom: HoveredAtom | null) => void
  /** fired when the mouse enters / leaves a clash line or protein-atom marker in the viewer */
  onClashHover?: (clash: ClashRef | null) => void
}

const GT_LIGAND = '#2f9e6b'
const POCKET_CUTOFF = 4.5

/**
 * 3Dmol.js scene: ground-truth receptor (cartoon) + ligand (sticks) and the selected
 * prediction, all in the ground-truth frame. Violations are drawn as translucent spheres
 * on the offending ligand atoms and dashed lines for protein–ligand clashes.
 */
const Viewer3D = forwardRef<ViewerHandle, ViewerProps>(function Viewer3D(props, ref) {
  const el = useRef<HTMLDivElement>(null)
  const viewer = useRef<$3Dmol.GLViewer | null>(null)
  const models = useRef<Models>({ labels: [] })
  const loadedKey = useRef<string>('')
  const propsRef = useRef(props)
  propsRef.current = props
  const applyHoverable = useRef<() => void>(() => {})
  const lastHover = useRef<number | null>(null)
  const cancelHover = useRef<() => void>(() => {})
  const pointerInside = useRef(false)
  // timestamps used to tell a real unhover (pointer moved) from one caused by a geometry rebuild after a restyle
  const lastMoveTs = useRef(0)
  const lastRestyleTs = useRef(0)
  const [glError, setGlError] = useState<string | null>(null)

  useImperativeHandle(ref, () => ({
    zoomToLigand: () => {
      const v = viewer.current
      if (!v) return
      const m = models.current.predLig ?? models.current.gtLig
      if (m) v.zoomTo({ model: m }, 300)
    },
    setFrame: (i: number) => {
      const v = viewer.current
      if (v && (models.current.traj || models.current.pocketTraj)) { v.setFrame(i); v.render() }
    },
  }))

  // create viewer once
  useEffect(() => {
    if (!el.current) return
    el.current.replaceChildren()  // StrictMode mounts twice in dev; drop any canvas left by a torn-down viewer
    let v: $3Dmol.GLViewer
    try {
      v = $3Dmol.createViewer(el.current, { backgroundColor: '#fbfbfb', antialias: true })
    } catch (e) {
      setGlError('WebGL is not available in this browser, so the 3D view cannot be shown.')
      return
    }
    viewer.current = v
    loadedKey.current = ''
    const ro = new ResizeObserver(() => { if (el.current && el.current.clientHeight > 0) { v.resize(); v.render() } })
    ro.observe(el.current)
    const node = el.current
    const onEnter = () => { pointerInside.current = true }
    const onMove = () => { lastMoveTs.current = performance.now() }
    const onLeave = () => { pointerInside.current = false; cancelHover.current(); clearLabels(v, models.current); v.render(); propsRef.current.onClashHover?.(null) }
    // capture phase so the timestamp is updated before 3Dmol's own mousemove handling runs
    node.addEventListener('mouseenter', onEnter); node.addEventListener('mousemove', onMove, true); node.addEventListener('mouseleave', onLeave)
    return () => { ro.disconnect(); node.removeEventListener('mouseenter', onEnter); node.removeEventListener('mousemove', onMove, true); node.removeEventListener('mouseleave', onLeave); v.clear(); viewer.current = null }
  }, [])

  // (re)load models when the structure strings change
  useEffect(() => {
    const v = viewer.current
    if (!v) return
    const key = [props.gtReceptor.length, props.gtLigand.length, props.predReceptor?.length, props.predLigand?.length, props.predTraj?.length, props.pocketTraj?.length, props.pocketTrajPdb?.length, props.predColor].join('|')
    if (key === loadedKey.current) return
    loadedKey.current = key
    v.removeAllModels(); v.removeAllShapes(); v.removeAllLabels()
    const m: Models = { labels: [] }
    m.gtRec = v.addModel(props.gtReceptor, 'pdb')
    m.gtLig = v.addModel(props.gtLigand, 'sdf')
    if (props.predReceptor) m.predRec = v.addModel(props.predReceptor, 'pdb')
    if (props.predLigand) m.predLig = v.addModel(props.predLigand, 'sdf')
    if (props.predTraj) m.traj = v.addModelsAsFrames(props.predTraj, 'sdf')
    if (props.pocketTraj) m.pocketTraj = v.addModelsAsFrames(props.pocketTraj, 'sdf')
    if (props.pocketTrajPdb) m.pocketRec = v.addModelsAsFrames(props.pocketTrajPdb, 'pdb')
    // pocket residues are computed geometrically (chain names differ between CIF-derived lists and the PDB files)
    const ligRef = m.predLig ?? m.gtLig
    m.gtPocket = pocketSelection(m.gtRec, ligRef)
    m.predPocket = m.predRec ? pocketSelection(m.predRec, ligRef) : []
    models.current = m

    // reverse hover: ligand atom under the mouse -> label + callback (throttled by atom identity)
    applyHoverable.current = () => {}
    if (m.predLig) {
      const enter = (atom: $3Dmol.AtomSpec) => {
        if (atom.index == null || atom.index === lastHover.current) return
        // ignore hovers whose delayed timer fires after the pointer already left the canvas
        if (!pointerInside.current) return
        lastHover.current = atom.index
        propsRef.current.onAtomHover?.({ index: atom.index, elem: atom.elem ?? '' })
      }
      const leave = (force = false) => {
        if (lastHover.current == null) return
        // a restyle rebuilds the geometry and makes 3Dmol report an unhover although the pointer never moved: ignore those
        if (!force && pointerInside.current && lastRestyleTs.current > lastMoveTs.current) return
        lastHover.current = null
        propsRef.current.onAtomHover?.(null)
      }
      cancelHover.current = () => leave(true)
      v.setHoverDuration(60)
      // re-applied after every restyle: setStyle can rebuild atom records and drop the hoverable flag
      applyHoverable.current = () => v.setHoverable({ model: m.predLig }, true, enter, () => leave(false))
      applyHoverable.current()
    }
    if (import.meta.env.DEV) (window as unknown as { __angstromViewer?: $3Dmol.GLViewer; __angstromModels?: Models }).__angstromViewer = v
    if (import.meta.env.DEV) (window as unknown as { __angstromModels?: Models }).__angstromModels = m

    styleAll(v, m, propsRef.current)
    v.zoomTo({ model: m.predLig ?? m.gtLig })
    v.zoom(0.8)
    v.render()
  }, [props.gtReceptor, props.gtLigand, props.predReceptor, props.predLigand, props.predTraj, props.pocketTraj, props.pocketTrajPdb, props.predColor])

  // restyle on toggles
  useEffect(() => {
    const v = viewer.current
    if (!v || !models.current.gtRec) return
    styleAll(v, models.current, props)
    lastRestyleTs.current = performance.now()
    v.render()
  }, [props.showGtReceptor, props.showGtLigand, props.showPredLigand, props.showPredReceptor, props.showPocket, props.showViolations, props.trajectoryMode, props.diagnostics, props.highlightAtoms, props.highlightClash, props.highlightResidue, props.highlightGtAtoms, props.atomDisplacement, props.predColor])

  useEffect(() => {
    const v = viewer.current
    if (!v || props.trajectoryMode === 'off') return
    if (!models.current.traj && !models.current.pocketTraj) return
    v.setFrame(props.frame); v.render()
  }, [props.frame, props.trajectoryMode])

  return (
    <div className="w-full h-full relative">
      <div ref={el} className="viewer-canvas absolute inset-0" />
      {glError && <div className="absolute inset-0 flex items-center justify-center text-fg-3 text-[12px] px-6 text-center">{glError}</div>}
    </div>
  )
})

export default Viewer3D

type ResidueSel = { chain: string; resi: number[] }
type Models = {
  gtRec?: $3Dmol.GLModel; gtLig?: $3Dmol.GLModel; predRec?: $3Dmol.GLModel; predLig?: $3Dmol.GLModel
  traj?: $3Dmol.GLModel; pocketTraj?: $3Dmol.GLModel; pocketRec?: $3Dmol.GLModel
  gtPocket?: ResidueSel[]; predPocket?: ResidueSel[]
  labels: $3Dmol.Label[]
}

/** Residues of `rec` with any heavy atom within POCKET_CUTOFF of any atom of `lig`, grouped by chain. */
function pocketSelection(rec: $3Dmol.GLModel, lig: $3Dmol.GLModel | undefined): ResidueSel[] {
  if (!lig) return []
  const L = (lig.selectedAtoms({}) as $3Dmol.AtomSpec[]).filter((a) => a.x != null)
  const R = rec.selectedAtoms({}) as $3Dmol.AtomSpec[]
  const c2 = POCKET_CUTOFF * POCKET_CUTOFF
  const hits = new Map<string, Set<number>>()
  for (const a of R) {
    if (a.x == null || a.resi == null) continue
    for (const l of L) {
      const dx = a.x - l.x!, dy = a.y! - l.y!, dz = a.z! - l.z!
      if (dx * dx + dy * dy + dz * dz <= c2) {
        const ch = a.chain ?? ''
        if (!hits.has(ch)) hits.set(ch, new Set())
        hits.get(ch)!.add(Number(a.resi))
        break
      }
    }
  }
  return [...hits.entries()].map(([chain, s]) => ({ chain, resi: [...s] }))
}

function styleAll(v: $3Dmol.GLViewer, m: Models, p: ViewerProps) {
  v.removeAllShapes()
  clearLabels(v, m)
  const hide: $3Dmol.AtomStyleSpec = {}  // empty style = not drawn

  const pocketMode = p.trajectoryMode === 'pocket' && !!m.pocketTraj
  const ligandMode = p.trajectoryMode === 'ligand' && !!m.traj
  const showTraj = pocketMode || ligandMode
  const pocketStick = { stick: { radius: 0.12, colorscheme: 'whiteCarbon' } }

  // ground-truth receptor: light cartoon; pocket residues as thin sticks (static ones hidden while the pocket animates)
  if (m.gtRec) {
    v.setStyle({ model: m.gtRec }, p.showGtReceptor ? { cartoon: { color: '#d9d9de', opacity: 0.9 } } : hide)
    if (p.showGtReceptor && p.showPocket && !pocketMode) {
      for (const sel of m.gtPocket ?? []) v.addStyle({ model: m.gtRec, chain: sel.chain, resi: sel.resi }, pocketStick)
    }
  }
  if (m.predRec) {
    v.setStyle({ model: m.predRec }, p.showPredReceptor ? { cartoon: { color: p.predColor, opacity: 0.55 } } : hide)
    if (p.showPredReceptor && p.showPocket && !pocketMode) {
      for (const sel of m.predPocket ?? []) v.addStyle({ model: m.predRec, chain: sel.chain, resi: sel.resi }, { stick: { radius: 0.12, colorscheme: { prop: 'elem', map: elemMap(lighten(p.predColor, 0.45)) } } })
    }
  }
  if (m.gtLig) {
    v.setStyle({ model: m.gtLig }, p.showGtLigand ? { stick: { radius: 0.18, colorscheme: { prop: 'elem', map: elemMap(GT_LIGAND) } } } : hide)
  }
  const predVisible = p.showPredLigand
  if (m.predLig) {
    v.setStyle({ model: m.predLig }, showTraj || !predVisible ? hide : { stick: { radius: 0.22, colorscheme: { prop: 'elem', map: elemMap(p.predColor) } } })
  }
  if (m.traj) {
    v.setStyle({ model: m.traj }, ligandMode && predVisible ? { stick: { radius: 0.22, colorscheme: { prop: 'elem', map: elemMap(p.predColor) } } } : hide)
  }
  if (m.pocketTraj) {
    v.setStyle({ model: m.pocketTraj }, pocketMode && predVisible ? { stick: { radius: 0.22, colorscheme: { prop: 'elem', map: elemMap(p.predColor) } } } : hide)
  }
  if (m.pocketRec) {
    // animated pocket heavy atoms (protein restrained, so they move only slightly)
    v.setStyle({ model: m.pocketRec }, pocketMode && p.showPocket ? pocketStick : hide)
  }

  const lig = pocketMode ? m.pocketTraj : ligandMode ? m.traj : m.predLig
  if (!lig || !predVisible) return

  if (showTraj && p.atomDisplacement && p.atomDisplacement.length) {
    // colour trajectory atoms by how far they move during minimisation
    const max = Math.max(0.5, ...p.atomDisplacement)
    const disp = p.atomDisplacement
    v.setStyle({ model: lig }, {
      stick: {
        radius: 0.22,
        colorfunc: (atom: $3Dmol.AtomSpec) => {
          const i = atom.index ?? 0
          const t = Math.min(1, (disp[i] ?? 0) / max)
          return lerpColor(p.predColor, '#d64545', t)
        },
      },
    })
  }

  if (p.showViolations && p.diagnostics && !showTraj) {
    const d = p.diagnostics
    const hcs = p.highlightClash == null ? [] : Array.isArray(p.highlightClash) ? p.highlightClash : [p.highlightClash]
    const hc = hcs.length ? hcs : null
    const isHl = (kind: ClashRef['kind'], i: number) => hcs.some((c) => c.kind === kind && c.index === i)
    const hlAtoms = new Set(p.highlightAtoms ?? [])
    const dimming = hc != null || hlAtoms.size > 0  // something is emphasised: fade every other violation marker
    if (d.flagged_atoms.length) {
      v.addStyle({ model: lig, index: d.flagged_atoms }, { sphere: { radius: 0.55, color: '#d64545', opacity: dimming ? 0.15 : 0.5 } })
    }
    const ligAtoms = lig.selectedAtoms({}) as $3Dmol.AtomSpec[]
    const hoverCbs = (ref: ClashRef) => ({
      hoverable: true,
      hover_callback: () => propsOnClashHover(p, ref),
      unhover_callback: () => propsOnClashHover(p, null),
    })
    const xyz = (a: $3Dmol.AtomSpec) => ({ x: a.x!, y: a.y!, z: a.z! })
    const labelAt = (a: $3Dmol.AtomSpec, b: $3Dmol.AtomSpec, text: string) => {
      m.labels.push(v.addLabel(text, {
        position: { x: (a.x! + b.x!) / 2, y: (a.y! + b.y!) / 2, z: (a.z! + b.z!) / 2 }, backgroundColor: '#f2b01e', backgroundOpacity: 0.95,
        fontColor: '#1c1c22', fontSize: 11, borderThickness: 0, inFront: true,
      } as $3Dmol.LabelSpec))
    }
    d.protein_clashes.forEach((c, i) => {
      const a = ligAtoms[c.atom]
      // clash was measured against the predicted receptor; draw against it, falling back to the GT receptor
      const sel = { chain: c.protein.chain, resi: c.protein.resnum, atom: c.protein.atom }
      const ref = (m.predRec?.selectedAtoms(sel) as $3Dmol.AtomSpec[] | undefined) ?? []
      const b = ref.length ? ref[0] : ((m.gtRec?.selectedAtoms(sel) as $3Dmol.AtomSpec[] | undefined) ?? [])[0]
      if (!(a && b && a.x != null && b.x != null)) return
      const hl = hc ? isHl('protein', i) : hlAtoms.has(c.atom)
      const atomHl = false
      const cb = hoverCbs({ kind: 'protein', index: i })
      v.addCylinder({ start: xyz(a), end: xyz(b), radius: hl ? 0.14 : 0.06, dashed: true, color: hl ? '#f2b01e' : '#d64545', alpha: dimming && !hl ? 0.2 : 1, fromCap: 1, toCap: 1, ...cb })
      v.addSphere({ center: xyz(b), radius: hl || atomHl ? 0.55 : 0.35, color: hl || atomHl ? '#f2b01e' : '#d64545', alpha: hl || atomHl ? 0.9 : dimming ? 0.12 : 0.35, ...cb })
      if (hl) { v.addSphere({ center: xyz(a), radius: 0.7, color: '#f2b01e', alpha: 0.85 }); labelAt(a, b, `${c.dist.toFixed(2)} Å`) }
    })
    d.intra_clashes.forEach((c, i) => {
      const a = ligAtoms[c.atoms[0]], b = ligAtoms[c.atoms[1]]
      if (!(a && b && a.x != null && b.x != null)) return
      const hl = hc ? isHl('intra', i) : (hlAtoms.has(c.atoms[0]) || hlAtoms.has(c.atoms[1]))
      const cb = hoverCbs({ kind: 'intra', index: i })
      v.addCylinder({ start: xyz(a), end: xyz(b), radius: hl ? 0.14 : 0.06, dashed: true, color: hl ? '#f2b01e' : '#d98c1c', alpha: dimming && !hl ? 0.2 : 1, ...cb })
      if (hl) { for (const q of [a, b]) v.addSphere({ center: xyz(q), radius: 0.7, color: '#f2b01e', alpha: 0.85 }); labelAt(a, b, `${c.dist.toFixed(2)} Å`) }
    })
  }

  // hover highlight goes last so it wins over the violation sphere on the same atom
  // contact-chip hover: residue side chain in the accent colour on whichever receptors are shown, plus crystal-ligand atoms
  if (p.highlightResidue) {
    const sel = parseResidue(p.highlightResidue)
    if (sel) {
      for (const rec of [p.showPredReceptor ? m.predRec : undefined, p.showGtReceptor ? m.gtRec : undefined]) {
        if (!rec) continue
        const chains = new Set((rec.selectedAtoms({}) as $3Dmol.AtomSpec[]).map((a) => a.chain))
        const s: $3Dmol.AtomSelectionSpec = chains.has(sel.chain) ? { model: rec, chain: sel.chain, resi: sel.resi } : { model: rec, resi: sel.resi, resn: sel.resn }
        v.addStyle(s, { stick: { radius: 0.3, colorscheme: { prop: 'elem', map: elemMap('#f2b01e') } } })
      }
    }
  }
  if (p.highlightGtAtoms && p.highlightGtAtoms.length && m.gtLig && p.showGtLigand) {
    v.addStyle({ model: m.gtLig, index: p.highlightGtAtoms }, { sphere: { radius: 0.6, color: '#2f9e6b', opacity: 0.6 } })
  }

  if (p.highlightAtoms && p.highlightAtoms.length) {
    v.addStyle({ model: lig, index: p.highlightAtoms }, { sphere: { radius: 0.7, color: '#f2b01e', opacity: 0.8 } })
  }

}

/** remove every label this component added (distance labels); called on every restyle and on leave */
function clearLabels(v: $3Dmol.GLViewer, m: Models) {
  for (const l of m.labels) v.removeLabel(l)
  m.labels = []
  if (import.meta.env.DEV) (window as unknown as { __angstromLabels?: () => number }).__angstromLabels = () => m.labels.length
}

/** clash-shape hover callbacks run inside 3Dmol's event loop; route them to the latest props */
function propsOnClashHover(p: ViewerProps, ref: ClashRef | null) {
  p.onClashHover?.(ref)
}

/** "A:PHE46" or "1.A:PHE46" -> { chain, resn, resi } */
function parseResidue(label: string): { chain: string; resn: string; resi: number } | null {
  const m = /^(.*?):([A-Za-z]+)(-?\d+)$/.exec(label)
  return m ? { chain: m[1], resn: m[2].toUpperCase(), resi: Number(m[3]) } : null
}

function elemMap(carbon: string): Record<string, string> {
  return { C: carbon, N: '#3f6fd6', O: '#d64545', S: '#d9b21c', P: '#e8842a', F: '#4cbf9a', Cl: '#3fb26b', Br: '#b0532c', I: '#8a3fb2', H: '#cfcfd6' }
}

function lerpColor(a: string, b: string, t: number): string {
  const pa = parseInt(a.slice(1), 16), pb = parseInt(b.slice(1), 16)
  const ch = (s: number) => Math.round(((pa >> s) & 255) * (1 - t) + ((pb >> s) & 255) * t)
  return `#${[16, 8, 0].map((s) => ch(s).toString(16).padStart(2, '0')).join('')}`
}

const lighten = (c: string, t: number) => lerpColor(c, '#ffffff', t)
