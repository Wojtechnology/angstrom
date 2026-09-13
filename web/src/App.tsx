import { Routes, Route, Link, useLocation } from 'react-router-dom'
import { Database } from 'lucide-react'
import Overview from './pages/Overview'
import SystemView from './pages/SystemView'

function Header() {
  const loc = useLocation()
  return (
    <header className="h-12 border-b hairline bg-panel/80 backdrop-blur sticky top-0 z-20">
      <div className="max-w-[1400px] mx-auto h-full px-5 flex items-center gap-6">
        <Link to="/" className="flex items-center gap-2 font-semibold text-[13px]">
          <span className="w-5 h-5 rounded-md bg-accent inline-flex items-center justify-center">
            <span className="w-2 h-2 rounded-full bg-white" />
          </span>
          Angstrom
        </Link>
        <nav className="flex items-center gap-1 text-fg-2">
          <Link to="/" className="btn" data-active={loc.pathname === '/'} style={{ border: 'none' }}>
            Overview
          </Link>
        </nav>
        <div className="ml-auto flex items-center gap-3 text-fg-3 text-[12px]">
          <span className="hidden sm:inline">Runs N&apos; Poses · demo subset</span>
          <a className="btn" style={{ border: 'none' }} href="https://github.com/plinder-org/runs-n-poses" target="_blank" rel="noreferrer">
            <Database size={14} /> dataset
          </a>
        </div>
      </div>
    </header>
  )
}

export default function App() {
  return (
    <div className="min-h-full flex flex-col">
      <Header />
      <main className="flex-1 max-w-[1400px] w-full mx-auto px-5 py-6">
        <Routes>
          <Route path="/" element={<Overview />} />
          <Route path="/system/:id" element={<SystemView />} />
        </Routes>
      </main>
    </div>
  )
}
