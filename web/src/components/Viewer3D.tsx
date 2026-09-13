import { useEffect, useImperativeHandle, useRef, useState, forwardRef } from 'react'
import * as $3Dmol from '3dmol'
import type { Diagnostics } from '../lib/api'

export interface ViewerHandle {
  zoomToLigand: () => void
  setFrame: (i: number) => void
}

export interface HoveredAtom { index: number; elem: string }

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
  /** fired (throttled) when the mouse enters / leaves a predicted-ligand atom */
  onAtomHover?: (atom: HoveredAtom | null) => void
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
  const models = useRef<Models>({})
  const loadedKey = useRef<string>('')
  const propsRef = useRef(props)
  propsRef.current = props
  const hoverLabel = useRef<$3Dmol.Label | null>(null)
  const lastHover = useRef<number | null>(null)
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
    return () => { ro.disconnect(); v.clear(); viewer.current = null }
  }, [])

  // (re)load models when the structure strings change
  useEffect(() => {
    const v = viewer.current
    if (!v) return
    const key = [props.gtReceptor.length, props.gtLigand.length, props.predReceptor?.length, props.predLigand?.length, props.predTraj?.length, props.pocketTraj?.length, props.pocketTrajPdb?.length, props.predColor].join('|')
    if (key === loadedKey.current) return
    loadedKey.current = key
    v.removeAllModels(); v.removeAllShapes(); v.removeAllLabels()
    hoverLabel.current = null
    const m: Models = {}
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
    if (m.predLig) {
      const enter = (atom: $3Dmol.AtomSpec) => {
        if (atom.index == null || atom.index === lastHover.current) return
        lastHover.current = atom.index
        if (hoverLabel.current) { v.removeLabel(hoverLabel.current); hoverLabel.current = null }
        hoverLabel.current = v.addLabel(`${atom.elem ?? ''}${atom.index + 1}`, {
          position: { x: atom.x!, y: atom.y!, z: atom.z! }, backgroundColor: '#1c1c22', backgroundOpacity: 0.85,
          fontColor: '#ffffff', fontSize: 11, borderThickness: 0, inFront: true, alignment: 'bottomLeft', screenOffset: { x: 8, y: -8 },
        } as $3Dmol.LabelSpec)
        v.render()
        propsRef.current.onAtomHover?.({ index: atom.index, elem: atom.elem ?? '' })
      }
      const leave = () => {
        if (lastHover.current == null) return
        lastHover.current = null
        if (hoverLabel.current) { v.removeLabel(hoverLabel.current); hoverLabel.current = null; v.render() }
        propsRef.current.onAtomHover?.(null)
      }
      v.setHoverDuration(60)
      v.setHoverable({ model: m.predLig }, true, enter, leave)
    }

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
    v.render()
  }, [props.showGtReceptor, props.showGtLigand, props.showPredLigand, props.showPredReceptor, props.showPocket, props.showViolations, props.trajectoryMode, props.diagnostics, props.highlightAtoms, props.atomDisplacement, props.predColor])

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

  if (p.showViolations && p.diagnostics && !showTraj) {
    const d = p.diagnostics
    if (d.flagged_atoms.length) {
      v.addStyle({ model: lig, index: d.flagged_atoms }, { sphere: { radius: 0.55, color: '#d64545', opacity: 0.5 } })
    }
    const ligAtoms = lig.selectedAtoms({}) as $3Dmol.AtomSpec[]
    for (const c of d.protein_clashes) {
      const a = ligAtoms[c.atom]
      // clash was measured against the predicted receptor; draw against it, falling back to the GT receptor
      const sel = { chain: c.protein.chain, resi: c.protein.resnum, atom: c.protein.atom }
      const ref = (m.predRec?.selectedAtoms(sel) as $3Dmol.AtomSpec[] | undefined) ?? []
      const b = ref.length ? ref[0] : ((m.gtRec?.selectedAtoms(sel) as $3Dmol.AtomSpec[] | undefined) ?? [])[0]
      if (a && b && a.x != null && b.x != null) {
        v.addCylinder({ start: { x: a.x, y: a.y!, z: a.z! }, end: { x: b.x, y: b.y!, z: b.z! }, radius: 0.06, dashed: true, color: '#d64545', fromCap: 1, toCap: 1 })
        const hl = p.highlightAtoms?.includes(c.atom)
        v.addSphere({ center: { x: b.x, y: b.y!, z: b.z! }, radius: hl ? 0.5 : 0.35, color: hl ? '#f2b01e' : '#d64545', alpha: hl ? 0.85 : 0.35 })
      }
    }
    for (const c of d.intra_clashes) {
      const a = ligAtoms[c.atoms[0]], b = ligAtoms[c.atoms[1]]
      if (a && b && a.x != null && b.x != null) v.addCylinder({ start: { x: a.x, y: a.y!, z: a.z! }, end: { x: b.x, y: b.y!, z: b.z! }, radius: 0.06, dashed: true, color: '#d98c1c' })
    }
  }

  // hover highlight goes last so it wins over the violation sphere on the same atom
  if (p.highlightAtoms && p.highlightAtoms.length) {
    v.addStyle({ model: lig, index: p.highlightAtoms }, { sphere: { radius: 0.75, color: '#f2b01e', opacity: 0.85 } })
  }

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
