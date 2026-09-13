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


def fetch(start, end, retries=8):
    for attempt in range(retries):
        try:
            r = requests.get(URL, headers={"Range": f"bytes={start}-{end}"}, timeout=120, stream=True)
            r.raise_for_status()
            data = r.content
            if len(data) != end - start + 1:
                raise IOError(f"short read {len(data)} != {end-start+1}")
            return data
        except Exception as e:  # noqa
            log(f"retry {attempt} for {start}-{end}: {e}")
            time.sleep(5 * (attempt + 1))
    raise RuntimeError("giving up")


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
    total = int(requests.head(URL, allow_redirects=True).headers["Content-Length"])
    log(f"total bytes {total}")
    stream = RangeStream(total)
    buffered = io.BufferedReader(stream, buffer_size=8 * 1024 * 1024)
    tar = tarfile.open(fileobj=buffered, mode="r|gz")
    n_ok = 0
    t0 = time.time()
    last_log = 0
    with open(MANIFEST, "w") as man:
        for member in tar:
            man.write(member.name + "\n")
            if time.time() - last_log > 60:
                last_log = time.time()
                log(f"{100 * stream.pos / total:.1f}% {stream.pos / (time.time() - t0) / 1e6:.1f} MB/s extracted={n_ok} at={member.name.split('prediction_files/')[-1][:60]}")
            parts = member.name.split("prediction_files/", 1)
            if len(parts) < 2 or not member.isfile():
                continue
            rel = parts[1]
            segs = rel.split("/")
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
