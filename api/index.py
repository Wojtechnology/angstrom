"""Minimal FastAPI backend serving the precomputed Runs N' Poses demo dataset.

Everything is static JSON / gzipped structure files under api/data, so the app has no
database and runs as a single Vercel serverless function.
"""
from __future__ import annotations

import json
import re
from functools import lru_cache
from pathlib import Path

from fastapi import FastAPI, HTTPException
from fastapi.responses import Response

DATA = Path(__file__).parent / "data"
SAFE = re.compile(r"^[A-Za-z0-9_.\-]+$")

app = FastAPI(title="Angstrom API", docs_url=None, redoc_url=None)


@lru_cache(maxsize=4)
def _index_cached(mtime: float) -> dict:
    return json.loads((DATA / "index.json").read_text())


def _index() -> dict:
    # keyed on mtime so a rebuilt index.json is picked up without restarting the server
    return _index_cached((DATA / "index.json").stat().st_mtime)


def _check(name: str) -> str:
    if not SAFE.match(name) or ".." in name:
        raise HTTPException(400, "bad name")
    return name


@app.get("/api/health")
def health():
    return {"ok": True, "systems": len(_index().get("systems", []))}


@app.get("/api/index")
def index():
    return _index()


@app.get("/api/systems/{system_id}")
def system(system_id: str):
    path = DATA / "systems" / f"{_check(system_id)}.json"
    if not path.exists():
        raise HTTPException(404, "unknown system")
    return Response(path.read_bytes(), media_type="application/json",
                    headers={"Cache-Control": "public, max-age=86400"})


@app.get("/api/structures/{system_id}/{name}")
def structure(system_id: str, name: str):
    """Serve pre-gzipped PDB / SDF files; the browser transparently decompresses them."""
    path = DATA / "structures" / _check(system_id) / (_check(name) + ".gz")
    if not path.exists():
        raise HTTPException(404, "unknown structure")
    media = "chemical/x-pdb" if name.endswith(".pdb") else "chemical/x-mdl-sdfile"
    return Response(path.read_bytes(), media_type=media,
                    headers={"Content-Encoding": "gzip", "Cache-Control": "public, max-age=86400"})
