import * as RSlider from '@radix-ui/react-slider'
import * as RSelect from '@radix-ui/react-select'
import * as RSwitch from '@radix-ui/react-switch'
import * as RTooltip from '@radix-ui/react-tooltip'
import { Check, ChevronDown } from 'lucide-react'
import type { ReactNode } from 'react'
import { DOCKING_BADGE, DOCKING_NOTE_DEFAULT, isDocking, type MethodInfo } from '../lib/api'

export function Slider({ value, onChange, min = 0, max = 100, step = 1 }: { value: number; onChange: (v: number) => void; min?: number; max?: number; step?: number }) {
  return (
    <RSlider.Root className="slider-root" value={[value]} min={min} max={max} step={step} onValueChange={(v) => onChange(v[0])}>
      <RSlider.Track className="slider-track">
        <RSlider.Range className="slider-range" />
      </RSlider.Track>
      <RSlider.Thumb className="slider-thumb" aria-label="threshold" />
    </RSlider.Root>
  )
}

export function Select<T extends string>({ value, onChange, options, width }: { value: T; onChange: (v: T) => void; options: { value: T; label: string }[]; width?: number }) {
  return (
    <RSelect.Root value={value} onValueChange={(v) => onChange(v as T)}>
      <RSelect.Trigger className="select-trigger" style={width ? { minWidth: width } : undefined}>
        <RSelect.Value />
        <RSelect.Icon><ChevronDown size={14} className="text-fg-3" /></RSelect.Icon>
      </RSelect.Trigger>
      <RSelect.Portal>
        <RSelect.Content className="select-content" position="popper" sideOffset={4}>
          <RSelect.Viewport>
            {options.map((o) => (
              <RSelect.Item key={o.value} value={o.value} className="select-item">
                <RSelect.ItemText>{o.label}</RSelect.ItemText>
                <RSelect.ItemIndicator className="ml-auto"><Check size={12} /></RSelect.ItemIndicator>
              </RSelect.Item>
            ))}
          </RSelect.Viewport>
        </RSelect.Content>
      </RSelect.Portal>
    </RSelect.Root>
  )
}

export function Switch({ checked, onChange, label }: { checked: boolean; onChange: (v: boolean) => void; label: ReactNode }) {
  return (
    <label className="flex items-center gap-2 cursor-pointer select-none text-fg-2">
      <RSwitch.Root className="switch-root" checked={checked} onCheckedChange={onChange}>
        <RSwitch.Thumb className="switch-thumb" />
      </RSwitch.Root>
      <span className="text-[12px]">{label}</span>
    </label>
  )
}

export function Checkbox({ checked, onChange, label, dot }: { checked: boolean; onChange: (v: boolean) => void; label: ReactNode; dot?: string }) {
  return (
    <label className="inline-flex items-center gap-1.5 cursor-pointer select-none text-[12px] text-fg-2 whitespace-nowrap">
      <input type="checkbox" checked={checked} onChange={(e) => onChange(e.target.checked)} className="w-3.5 h-3.5 rounded-[3px] border border-line-2 accent-[#5e6ad2] cursor-pointer" />
      {dot && <span className="w-2 h-2 rounded-full inline-block" style={{ background: dot }} />}
      {label}
    </label>
  )
}

export function Tip({ children, content }: { children: ReactNode; content: ReactNode }) {
  return (
    <RTooltip.Root>
      <RTooltip.Trigger asChild>{children}</RTooltip.Trigger>
      <RTooltip.Portal>
        <RTooltip.Content className="tooltip" sideOffset={6}>{content}</RTooltip.Content>
      </RTooltip.Portal>
    </RTooltip.Root>
  )
}

/** Small warn-style badge marking a docking baseline; full note in the tooltip. Renders nothing for co-folding methods. */
export function MethodBadge({ method, className = '' }: { method: MethodInfo | undefined | null; className?: string }) {
  if (!isDocking(method)) return null
  return (
    <Tip content={method!.note ?? DOCKING_NOTE_DEFAULT}>
      <span className={`chip chip-warn ${className}`} style={{ height: 16, fontSize: 10, padding: '0 5px', textTransform: 'none', letterSpacing: 0, fontWeight: 500 }}>{DOCKING_BADGE}</span>
    </Tip>
  )
}

export function Label({ children }: { children: ReactNode }) {
  return <div className="text-[11px] uppercase tracking-wider text-fg-3 font-medium">{children}</div>
}

export function Stat({ label, value, sub }: { label: ReactNode; value: ReactNode; sub?: ReactNode }) {
  return (
    <div className="card px-4 py-3 min-w-0">
      <div className="text-[11px] text-fg-3 font-medium">{label}</div>
      <div className="text-[20px] font-semibold tracking-tight mt-0.5 tabular-nums">{value}</div>
      {sub && <div className="text-[11px] text-fg-3 mt-0.5 truncate">{sub}</div>}
    </div>
  )
}

export function Spinner({ label = 'Loading' }: { label?: string }) {
  return (
    <div className="flex items-center gap-2 text-fg-3 text-[12px] py-10 justify-center">
      <span className="w-3 h-3 rounded-full border-2 border-line-2 border-t-accent animate-spin" /> {label}
    </div>
  )
}

export function ErrorBox({ error }: { error: unknown }) {
  return <div className="card px-4 py-3 text-bad text-[12px]">{String(error)}</div>
}
