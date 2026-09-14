"""Stream prediction_files.tar.gz (~40 GB) from Zenodo and extract only the demo subset.

gzip is not seekable, so the whole archive has to be streamed, but nothing except the
wanted members is written to disk. Chunks are prefetched in parallel HTTP range requests
and fed in order to a streaming tarfile reader.
"""
import io
import json
import sys
import time
import tarfile
import requests
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

URL = "https://zenodo.org/api/records/18366081/files/prediction_files.tar.gz/content"
CHUNK = 32 * 1024 * 1024
PREFETCH = 12
HERE = Path(__file__).parent
OUT = HERE / "raw" / "prediction_files"
LOG = HERE / "raw" / "stream_extract.log"
MANIFEST = HERE / "raw" / "prediction_files_manifest.txt"

subset = json.load(open(HERE / "subset.json"))
WANTED = set(subset["systems"])
WANTED_METHODS = {"af3", "boltz", "boltz1x"}  # only these method dirs are extracted


def fetch(start, end):
    """Retry forever with capped backoff: gzip cannot be resumed mid-stream, so a Zenodo
    outage must be waited out rather than abandoned."""
    attempt = 0
    while True:
        try:
            r = requests.get(URL, headers={"Range": f"bytes={start}-{end}"}, timeout=180, stream=True)
            r.raise_for_status()
            data = r.content
            if len(data) != end - start + 1:
                raise IOError(f"short read {len(data)} != {end-start+1}")
            return data
        except Exception as e:  # noqa
            attempt += 1
            log(f"retry {attempt} for {start}-{end}: {e}")
            time.sleep(min(120, 5 * attempt))


def _pos(stream):
    return stream.pos if hasattr(stream, "pos") else stream.tell()


def log(msg):
    with open(LOG, "a") as f:
        f.write(f"{time.strftime('%H:%M:%S')} {msg}\n")


class RangeStream(io.RawIOBase):
    def __init__(self, total):
        self.total = total
        self.pool = ThreadPoolExecutor(PREFETCH)
        self.futures = {}
        self.next_submit = 0
        self.pos = 0
        self.buf = b""
        for _ in range(PREFETCH):
            self._submit()

    def _submit(self):
        if self.next_submit >= self.total:
            return
        s = self.next_submit
        e = min(s + CHUNK, self.total) - 1
        self.futures[s] = self.pool.submit(fetch, s, e)
        self.next_submit = e + 1

    def readable(self):
        return True

    def readinto(self, b):
        n = len(b)
        while len(self.buf) < n and self.pos < self.total:
            fut = self.futures.pop(self.pos)
            data = fut.result()
            self.pos += len(data)
            self.buf += data
            self._submit()
        out = self.buf[:n]
        self.buf = self.buf[n:]
        b[: len(out)] = out
        return len(out)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--local", help="read a locally downloaded prediction_files.tar.gz instead of streaming from Zenodo")
    ap.add_argument("--methods", nargs="*", help="override WANTED_METHODS (e.g. --methods chai protenix)")
    args = ap.parse_args()
    global WANTED_METHODS
    if args.methods:
        WANTED_METHODS = set(args.methods)
    if args.local:
        run(open(args.local, "rb"), Path(args.local).stat().st_size, buffer_size=32 * 1024 * 1024)
        return
    total = 0
    while total < 1_000_000_000:  # Zenodo returns a small error page on 504; wait it out
        try:
            r = requests.head(URL, allow_redirects=True, timeout=120)
            r.raise_for_status()
            total = int(r.headers["Content-Length"])
        except Exception as e:  # noqa
            log(f"HEAD failed: {e}")
            time.sleep(60)
    log(f"total bytes {total}")
    run(RangeStream(total), total)


def run(raw_stream, total, buffer_size=8 * 1024 * 1024):
    stream = raw_stream
    buffered = io.BufferedReader(stream, buffer_size=buffer_size)
    tar = tarfile.open(fileobj=buffered, mode="r|gz")
    n_ok = 0
    seen_wanted, seen_other = set(), set()
    t0 = time.time()
    last_log = 0
    with open(MANIFEST, "w") as man:
        for member in tar:
            man.write(member.name + "\n")
            if time.time() - last_log > 60:
                last_log = time.time()
                log(f"{100 * _pos(stream) / total:.1f}% {_pos(stream) / (time.time() - t0) / 1e6:.1f} MB/s extracted={n_ok} at={member.name.split('prediction_files/')[-1][:60]}")
            parts = member.name.split("prediction_files/", 1)
            if len(parts) < 2 or not member.isfile():
                continue
            rel = parts[1]
            segs = rel.split("/")
            if segs and segs[0] not in WANTED_METHODS:
                seen_other.add(segs[0])
                # tar is grouped by method: once every wanted dir has been passed, stop streaming
                if WANTED_METHODS <= seen_wanted:
                    log(f"all wanted method dirs done; stopping at {segs[0]}")
                    break
                continue
            if segs:
                seen_wanted.add(segs[0])
            if len(segs) < 3 or segs[1] not in WANTED:
                continue
            dest = OUT / rel
            dest.parent.mkdir(parents=True, exist_ok=True)
            with tar.extractfile(member) as src, open(dest, "wb") as dst:
                dst.write(src.read())
            n_ok += 1
    log(f"done extracted={n_ok} in {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
