import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'
import tailwindcss from '@tailwindcss/vite'

// https://vite.dev/config/
export default defineConfig({
  plugins: [react(), tailwindcss()],
  server: {
    proxy: {
      // local dev: FastAPI via `uvicorn index:app --port 8000` from ./api
      '/api': `http://127.0.0.1:${process.env.API_PORT || 8000}`,
    },
  },
  build: {
    chunkSizeWarningLimit: 2000,
  },
})
