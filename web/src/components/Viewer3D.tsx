import { useEffect, useImperativeHandle, useRef, useState, forwardRef } from 'react'
import * as $3Dmol from '3dmol'
import type { Diagnostics } from '../lib/api'

export interface ViewerHandle {
  zoomToLigand: () => void
  setFrame: (i: number) => void
}

export interface ViewerProps {
  gtReceptor: string
  gtLigand: string
  predReceptor?: string | null
  predLigand?: string | null
  predTraj?: string | null
  pocketTraj?: string | null
  pocketTrajPdb?: string | null
  predColor: string
  showPredReceptor: boolean
  showGtLigand: boolean
  showPocket: boolean
  showViolations: boolean
  trajectoryMode: 'off' | 'ligand' | 'pocket'
  frame: number
  diagnostics?: Diagnostics | null
  atomDisplacement?: number[] | null
  highlightAtoms?: number[] | null
  onAtomHover?: (label: string | null) => void
}

const GT_LIGAND = '#2f9e6b'

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
    const m: typeof models.current = {}
    m.gtRec = v.addModel(props.gtReceptor, 'pdb')
    m.gtLig = v.addModel(props.gtLigand, 'sdf')
    if (props.predReceptor) m.predRec = v.addModel(props.predReceptor, 'pdb')
    if (props.predLigand) m.predLig = v.addModel(props.predLigand, 'sdf')
    if (props.predTraj) m.traj = v.addModelsAsFrames(props.predTraj, 'sdf')
    if (props.pocketTraj) m.pocketTraj = v.addModelsAsFrames(props.pocketTraj, 'sdf')
    if (props.pocketTrajPdb) m.pocketRec = v.addModelsAsFrames(props.pocketTrajPdb, 'pdb')
    models.current = m
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
  }, [props.showPredReceptor, props.showGtLigand, props.showPocket, props.showViolations, props.trajectoryMode, props.diagnostics, props.highlightAtoms, props.atomDisplacement, props.predColor])

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

type Models = { gtRec?: $3Dmol.GLModel; gtLig?: $3Dmol.GLModel; predRec?: $3Dmol.GLModel; predLig?: $3Dmol.GLModel; traj?: $3Dmol.GLModel; pocketTraj?: $3Dmol.GLModel; pocketRec?: $3Dmol.GLModel }

function styleAll(v: $3Dmol.GLViewer, m: Models, p: ViewerProps) {
  v.removeAllShapes()
  v.removeAllLabels()
  const hide: $3Dmol.AtomStyleSpec = {}  // empty style = not drawn

  const pocketMode = p.trajectoryMode === 'pocket' && !!m.pocketTraj
  const ligandMode = p.trajectoryMode === 'ligand' && !!m.traj
  const showTraj = pocketMode || ligandMode

  // ground-truth receptor: light cartoon; pocket residues as thin sticks (static ones hidden while the pocket animates)
  if (m.gtRec) {
    v.setStyle({ model: m.gtRec }, { cartoon: { color: '#d9d9de', opacity: 0.9 } })
    if (p.showPocket && m.gtLig && !pocketMode) {
      v.addStyle({ model: m.gtRec, within: { distance: 4.5, sel: { model: m.predLig ?? m.gtLig } } } as $3Dmol.AtomSelectionSpec, { stick: { radius: 0.12, colorscheme: 'whiteCarbon' } })
    }
  }
  if (m.predRec) {
    v.setStyle({ model: m.predRec }, p.showPredReceptor ? { cartoon: { color: p.predColor, opacity: 0.45 } } : hide)
  }
  if (m.gtLig) {
    v.setStyle({ model: m.gtLig }, p.showGtLigand ? { stick: { radius: 0.18, colorscheme: { prop: 'elem', map: elemMap(GT_LIGAND) } } } : hide)
  }
  if (m.predLig) {
    v.setStyle({ model: m.predLig }, showTraj ? hide : { stick: { radius: 0.22, colorscheme: { prop: 'elem', map: elemMap(p.predColor) } } })
  }
  if (m.traj) {
    v.setStyle({ model: m.traj }, ligandMode ? { stick: { radius: 0.22, colorscheme: { prop: 'elem', map: elemMap(p.predColor) } } } : hide)
  }
  if (m.pocketTraj) {
    v.setStyle({ model: m.pocketTraj }, pocketMode ? { stick: { radius: 0.22, colorscheme: { prop: 'elem', map: elemMap(p.predColor) } } } : hide)
  }
  if (m.pocketRec) {
    // animated pocket heavy atoms (protein restrained, so they move only slightly)
    v.setStyle({ model: m.pocketRec }, pocketMode ? { stick: { radius: 0.12, colorscheme: 'whiteCarbon' } } : hide)
  }

  const lig = pocketMode ? m.pocketTraj : ligandMode ? m.traj : m.predLig
  if (!lig) return

  if (p.showViolations && p.diagnostics && !showTraj) {
    const d = p.diagnostics
    if (d.flagged_atoms.length) {
      v.addStyle({ model: lig, index: d.flagged_atoms }, { sphere: { radius: 0.55, color: '#d64545', opacity: 0.5 } })
    }
    const ligAtoms = lig.selectedAtoms({}) as $3Dmol.AtomSpec[]
    for (const c of d.protein_clashes) {
      const a = ligAtoms[c.atom]
      const target = m.gtRec && !p.showPredReceptor ? m.gtRec : (m.predRec ?? m.gtRec)
      const pa = target?.selectedAtoms({ chain: c.protein.chain, resi: c.protein.resnum, atom: c.protein.atom }) as $3Dmol.AtomSpec[] | undefined
      // clash was measured against the predicted receptor; draw against it when available
      const ref = (m.predRec?.selectedAtoms({ chain: c.protein.chain, resi: c.protein.resnum, atom: c.protein.atom }) as $3Dmol.AtomSpec[] | undefined) ?? pa
      const b = ref && ref.length ? ref[0] : undefined
      if (a && b && a.x != null && b.x != null) {
        v.addCylinder({ start: { x: a.x, y: a.y!, z: a.z! }, end: { x: b.x, y: b.y!, z: b.z! }, radius: 0.06, dashed: true, color: '#d64545', fromCap: 1, toCap: 1 })
        v.addSphere({ center: { x: b.x, y: b.y!, z: b.z! }, radius: 0.35, color: '#d64545', alpha: 0.35 })
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
