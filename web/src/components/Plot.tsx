import { useEffect, useRef } from 'react'
import Plotly from 'plotly.js-basic-dist-min'
import type { Data, Layout, Config } from 'plotly.js-basic-dist-min'

const BASE_LAYOUT: Partial<Layout> = {
  font: { family: 'Inter, ui-sans-serif, system-ui, sans-serif', size: 11, color: '#55555f' },
  paper_bgcolor: 'rgba(0,0,0,0)',
  plot_bgcolor: 'rgba(0,0,0,0)',
  margin: { l: 44, r: 12, t: 8, b: 36 },
  legend: { orientation: 'h', y: -0.22, x: 0, font: { size: 11 } },
  xaxis: { gridcolor: '#eeeef1', zerolinecolor: '#e6e6e9', linecolor: '#e6e6e9', tickfont: { size: 11 } },
  yaxis: { gridcolor: '#eeeef1', zerolinecolor: '#e6e6e9', linecolor: '#e6e6e9', tickfont: { size: 11 } },
  hoverlabel: { bgcolor: '#1c1c22', bordercolor: '#1c1c22', font: { color: '#fff', size: 11, family: 'Inter, system-ui, sans-serif' } },
}

const CONFIG: Partial<Config> = { displayModeBar: false, responsive: true }

export default function Plot({ data, layout, height = 260, onClick }: { data: Data[]; layout?: Partial<Layout>; height?: number; onClick?: (e: Plotly.PlotMouseEvent) => void }) {
  const ref = useRef<HTMLDivElement>(null)
  useEffect(() => {
    const el = ref.current
    if (!el) return
    const merged: Partial<Layout> = {
      ...BASE_LAYOUT,
      ...layout,
      height,
      xaxis: { ...BASE_LAYOUT.xaxis, ...layout?.xaxis },
      yaxis: { ...BASE_LAYOUT.yaxis, ...layout?.yaxis },
    }
    // StrictMode mounts/unmounts twice in dev; a purge between the two makes the first
    // react() promise reject when it tries to emit on the purged element. Harmless.
    const gd = el as unknown as Plotly.PlotlyHTMLElement
    Plotly.react(el, data, merged, CONFIG).then(() => {
      if (onClick && typeof gd.on === 'function') { gd.removeAllListeners?.('plotly_click'); gd.on('plotly_click', onClick) }
    }).catch(() => {})
    return () => { gd.removeAllListeners?.('plotly_click') }
  }, [data, layout, height, onClick])
  useEffect(() => {
    const el = ref.current
    return () => { if (el) Plotly.purge(el) }
  }, [])
  return <div ref={ref} style={{ width: '100%', height }} />
}
