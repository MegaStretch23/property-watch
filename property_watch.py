#!/usr/bin/env python3
# /// script
# requires-python = ">=3.9"
# dependencies = ["requests>=2.31"]
# ///
"""property-watch: look for large changes at one place over a 12-month window using free Sentinel-2 imagery.

  uv run property_watch.py "Hoover Dam" --end 2026-10 --out report.html
  python3 property_watch.py --lat 36.0157 --lon -114.7375 --months 12 --size 800 --json result.json

What it does: geocode an address (free Nominatim) or take lat/lon, query the free Element84 Earth Search STAC API for Sentinel-2 L2A
scenes over a small box, keep the clearest scene per month (lowest cloud cover, checked again on the clip itself), read only the
true-colour (TCI) pixels of that box out of the cloud-optimised GeoTIFF with HTTP range requests, compare the monthly clips with a
pixel-difference score and write ONE self-contained HTML report (before/after slider + change heatmap).

Sentinel-2 TCI is 10 m per pixel. It can show large construction, clearing, new roofs or paving and big ground changes. It cannot
see anything a few pixels wide (a pool, a shed, a car, a fence). Standard library plus requests only."""
from __future__ import annotations

import argparse
import base64
import datetime as dt
import html
import json
import math
import re
import struct
import sys
import zlib
from itertools import accumulate
from pathlib import Path

__version__ = "0.1.0"
UA = "property-watch/0.1 (+https://github.com/MegaStretch23/property-watch)"
STAC_URL = "https://earth-search.aws.element84.com/v1/search"
NOMINATIM = "https://nominatim.openstreetmap.org/search"
PIXEL_M = 10.0


class WatchError(Exception):
    """Raised for anything the run cannot continue past (network, bad input, unsupported file)."""


# ----------------------------------------------------------------------------- PNG (stdlib)
def png_encode(width: int, height: int, data: bytes, channels: int = 3) -> bytes:
    if channels not in (3, 4) or len(data) != width * height * channels:
        raise ValueError("png_encode: data size does not match width*height*channels")
    raw = b"".join(b"\x00" + bytes(data[y * width * channels:(y + 1) * width * channels]) for y in range(height))

    def chunk(tag: bytes, body: bytes) -> bytes:
        return struct.pack(">I", len(body)) + tag + body + struct.pack(">I", zlib.crc32(tag + body) & 0xFFFFFFFF)

    ihdr = struct.pack(">IIBBBBB", width, height, 8, 2 if channels == 3 else 6, 0, 0, 0)
    return b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", ihdr) + chunk(b"IDAT", zlib.compress(raw, 9)) + chunk(b"IEND", b"")


def png_decode(blob: bytes) -> tuple[int, int, int, bytes]:
    """Decode an 8-bit, non-interlaced RGB or RGBA PNG -> (width, height, channels, data)."""
    if blob[:8] != b"\x89PNG\r\n\x1a\n":
        raise WatchError("not a PNG")
    pos, idat, hdr = 8, b"", None
    while pos < len(blob):
        n, tag = struct.unpack(">I4s", blob[pos:pos + 8])
        body = blob[pos + 8:pos + 8 + n]
        pos += 12 + n
        if tag == b"IHDR":
            hdr = struct.unpack(">IIBBBBB", body)
        elif tag == b"IDAT":
            idat += body
    if not hdr or hdr[2] != 8 or hdr[3] not in (2, 6) or hdr[6] != 0:
        raise WatchError("unsupported PNG (need 8-bit RGB/RGBA, not interlaced)")
    w, h, ch = hdr[0], hdr[1], 3 if hdr[3] == 2 else 4
    raw, stride, out, prev = zlib.decompress(idat), w * ch, bytearray(), bytearray(w * ch)
    for y in range(h):
        f = raw[y * (stride + 1)]
        line = bytearray(raw[y * (stride + 1) + 1:(y + 1) * (stride + 1)])
        for i in range(stride):
            a = line[i - ch] if i >= ch else 0
            b = prev[i]
            c = prev[i - ch] if i >= ch else 0
            if f == 1:
                line[i] = (line[i] + a) & 255
            elif f == 2:
                line[i] = (line[i] + b) & 255
            elif f == 3:
                line[i] = (line[i] + (a + b) // 2) & 255
            elif f == 4:
                p = a + b - c
                pa, pb, pc = abs(p - a), abs(p - b), abs(p - c)
                pr = a if pa <= pb and pa <= pc else (b if pb <= pc else c)
                line[i] = (line[i] + pr) & 255
        out += line
        prev = line
    return w, h, ch, bytes(out)


class Image:
    """RGB image, row-major bytes."""

    def __init__(self, w: int, h: int, data: bytes):
        if len(data) != w * h * 3:
            raise ValueError("Image: bad data length")
        self.w, self.h, self.data = w, h, bytes(data)

    def png(self) -> bytes:
        return png_encode(self.w, self.h, self.data, 3)

    @staticmethod
    def from_png(blob: bytes) -> "Image":
        w, h, ch, d = png_decode(blob)
        if ch == 4:
            d = b"".join(d[i:i + 3] for i in range(0, len(d), 4))
        return Image(w, h, d)


# ----------------------------------------------------------------------------- projection
def utm_forward(lat: float, lon: float, epsg: int) -> tuple[float, float]:
    """WGS84 lat/lon -> UTM easting/northing for EPSG 326zz (north) or 327zz (south). Kruger series, sub-millimetre."""
    if 32601 <= epsg <= 32660:
        zone, south = epsg - 32600, False
    elif 32701 <= epsg <= 32760:
        zone, south = epsg - 32700, True
    else:
        raise WatchError(f"unsupported CRS EPSG:{epsg} (only WGS84 UTM zones)")
    a, f, k0 = 6378137.0, 1 / 298.257223563, 0.9996
    n = f / (2 - f)
    A = a / (1 + n) * (1 + n ** 2 / 4 + n ** 4 / 64 + n ** 6 / 256)
    alpha = [n / 2 - 2 * n ** 2 / 3 + 5 * n ** 3 / 16 + 41 * n ** 4 / 180,
             13 * n ** 2 / 48 - 3 * n ** 3 / 5 + 557 * n ** 4 / 1440,
             61 * n ** 3 / 240 - 103 * n ** 4 / 140,
             49561 * n ** 4 / 161280]
    phi, lam = math.radians(lat), math.radians(lon - ((zone - 1) * 6 - 180 + 3))
    e = math.sqrt(f * (2 - f))
    t = math.sinh(math.atanh(math.sin(phi)) - e * math.atanh(e * math.sin(phi)))
    xi = math.atan2(t, math.cos(lam))
    eta = math.atanh(math.sin(lam) / math.sqrt(1 + t * t))
    x = eta + sum(alpha[j] * math.cos(2 * (j + 1) * xi) * math.sinh(2 * (j + 1) * eta) for j in range(4))
    y = xi + sum(alpha[j] * math.sin(2 * (j + 1) * xi) * math.cosh(2 * (j + 1) * eta) for j in range(4))
    return 500000.0 + k0 * A * x, (10000000.0 if south else 0.0) + k0 * A * y


# ----------------------------------------------------------------------------- HTTP
def http_get(url: str, headers: dict | None = None, params: dict | None = None, tries: int = 3):
    import requests
    last = None
    for i in range(tries):
        try:
            r = requests.get(url, headers={"User-Agent": UA, **(headers or {})}, params=params, timeout=60)
            if r.status_code < 500 and r.status_code != 429:
                return r
            last = f"HTTP {r.status_code}"
        except requests.RequestException as e:   # network error: retry, then fail loudly
            last = f"{type(e).__name__}: {e}"
        import time
        time.sleep(1.5 * (i + 1))
    raise WatchError(f"GET {url} failed after {tries} tries: {last}")


def http_range(url: str, start: int, length: int) -> bytes:
    r = http_get(url, headers={"Range": f"bytes={start}-{start + length - 1}"})
    if r.status_code == 206:
        return r.content
    if r.status_code == 200:   # server ignored Range: slice locally
        return r.content[start:start + length]
    raise WatchError(f"range read {url} -> HTTP {r.status_code}")


def http_post_json(url: str, body: dict, tries: int = 3) -> dict:
    import requests
    last = None
    for i in range(tries):
        try:
            r = requests.post(url, json=body, headers={"User-Agent": UA}, timeout=60)
            if r.status_code == 200:
                return r.json()
            last = f"HTTP {r.status_code} {r.text[:160]}"
            if r.status_code < 500 and r.status_code != 429:
                break
        except requests.RequestException as e:
            last = f"{type(e).__name__}: {e}"
        import time
        time.sleep(1.5 * (i + 1))
    raise WatchError(f"POST {url} failed: {last}")


# ----------------------------------------------------------------------------- geocoding
LATLON_RX = re.compile(r"^\s*(-?\d{1,2}(?:\.\d+)?)\s*[,; ]\s*(-?\d{1,3}(?:\.\d+)?)\s*$")


def parse_latlon(text: str):
    m = LATLON_RX.match(text or "")
    if not m:
        return None
    lat, lon = float(m.group(1)), float(m.group(2))
    if not (-90 <= lat <= 90 and -180 <= lon <= 180):
        raise WatchError(f"lat/lon out of range: {lat}, {lon}")
    return lat, lon


def geocode(query: str) -> tuple[float, float, str]:
    r = http_get(NOMINATIM, params={"q": query, "format": "json", "limit": 1})
    if r.status_code != 200:
        raise WatchError(f"Nominatim HTTP {r.status_code}")
    rows = r.json()
    if not rows:
        raise WatchError(f"Nominatim found nothing for {query!r}; pass --lat and --lon instead")
    return float(rows[0]["lat"]), float(rows[0]["lon"]), rows[0].get("display_name", query)


# ----------------------------------------------------------------------------- months / STAC
def month_list(end: str, n: int) -> list[str]:
    y, m = int(end[:4]), int(end[5:7])
    out = []
    for _ in range(n):
        out.append(f"{y:04d}-{m:02d}")
        m -= 1
        if m == 0:
            y, m = y - 1, 12
    return out[::-1]


def default_end(today: dt.date | None = None) -> str:
    t = today or dt.date.today()
    y, m = (t.year, t.month - 1) if t.month > 1 else (t.year - 1, 12)
    return f"{y:04d}-{m:02d}"


def bbox_around(lat: float, lon: float, size_m: float) -> list[float]:
    dlat = size_m / 2 / 111320.0
    dlon = dlat / max(0.05, math.cos(math.radians(lat)))
    return [round(lon - dlon, 6), round(lat - dlat, 6), round(lon + dlon, 6), round(lat + dlat, 6)]


def stac_search(lat: float, lon: float, size_m: float, months: list[str], post=http_post_json) -> list[dict]:
    first, last = months[0], months[-1]
    ly, lm = int(last[:4]), int(last[5:7])
    last_day = (dt.date(ly + (lm == 12), lm % 12 + 1, 1) - dt.timedelta(days=1)).day
    body = {"collections": ["sentinel-2-l2a"], "bbox": bbox_around(lat, lon, size_m), "limit": 200,
            "datetime": f"{first}-01T00:00:00Z/{last}-{last_day:02d}T23:59:59Z"}
    feats: list[dict] = []
    for _ in range(8):   # follow 'next' pages
        page = post(STAC_URL, body)
        feats += page.get("features", [])
        nxt = next((l for l in page.get("links", []) if l.get("rel") == "next"), None)
        if not nxt or not page.get("features"):
            break
        body = nxt.get("body") or body
    return feats


def window_fits(item: dict, lat: float, lon: float, size_px: int) -> bool:
    """Metadata-only check that the clip box lies fully inside the scene raster (STAC proj:transform / proj:shape on the visual asset)."""
    try:
        a = item["assets"]["visual"]
        tr, shape = a["proj:transform"], a["proj:shape"]
        epsg = item["properties"].get("proj:epsg") or int(str(item["properties"].get("proj:code", "EPSG:0"))[5:])
        x, y = utm_forward(lat, lon, int(epsg))
    except (KeyError, ValueError, TypeError):
        return True   # no metadata: let the pixel reader decide
    c, r = (x - tr[2]) / tr[0], (y - tr[5]) / tr[4]
    half = size_px / 2 + 1
    return half <= c <= shape[1] - half and half <= r <= shape[0] - half


def pick_tile(items: list[dict], months: list[str]) -> str | None:
    """The MGRS tile that has scenes in the most months (ties: more scenes). One tile = one pixel grid for every month."""
    score: dict[str, tuple[set, int]] = {}
    for it in items:
        g = it["properties"].get("grid:code") or it["properties"].get("s2:mgrs_tile") or "?"
        s = score.setdefault(g, (set(), 0))
        s[0].add(it["properties"]["datetime"][:7])
        score[g] = (s[0], s[1] + 1)
    wanted = set(months)
    best = sorted(score.items(), key=lambda kv: (-len(kv[1][0] & wanted), -kv[1][1], kv[0]))
    return best[0][0] if best else None


def candidates_by_month(items: list[dict], months: list[str], tile: str | None) -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = {m: [] for m in months}
    for it in items:
        p = it["properties"]
        g = p.get("grid:code") or p.get("s2:mgrs_tile") or "?"
        if tile and g != tile:
            continue
        m = p["datetime"][:7]
        if m in out and "visual" in it.get("assets", {}):
            out[m].append(it)
    for m in out:
        out[m].sort(key=lambda it: (it["properties"].get("eo:cloud_cover", 100), it["properties"]["datetime"]))
    return out


# ----------------------------------------------------------------------------- Cloud-optimised GeoTIFF window reader
class Cog:
    """Reads a pixel window of a tiled, deflate-compressed, 8-bit RGB GeoTIFF (what Earth Search serves as TCI.tif) with range requests."""

    def __init__(self, url: str, rng=http_range):
        self.url, self.rng = url, rng
        head = rng(url, 0, 65536)
        if head[:2] == b"II":
            self.bo = "<"
        elif head[:2] == b"MM":
            self.bo = ">"
        else:
            raise WatchError("not a TIFF")
        magic, off = struct.unpack(self.bo + "HI", head[2:8])
        if magic != 42:
            raise WatchError("BigTIFF is not supported")
        self._head = head
        n = struct.unpack(self.bo + "H", self._bytes(off, 2))[0]
        ent = self._bytes(off + 2, n * 12)
        self.tags = {}
        for i in range(n):
            tag, typ, cnt = struct.unpack(self.bo + "HHI", ent[i * 12:i * 12 + 8])
            self.tags[tag] = (typ, cnt, ent[i * 12 + 8:i * 12 + 12])
        g = self.first
        self.width, self.height = g(256), g(257)
        self.tw, self.th = g(322), g(323)
        if not (self.tw and self.th):
            raise WatchError("TIFF is not tiled")
        if g(259) not in (8, 32946):
            raise WatchError(f"unsupported TIFF compression {g(259)} (need deflate)")
        if g(284, 1) != 1 or list(self.values(258)) != [8, 8, 8] or g(277) != 3:
            raise WatchError("unsupported TIFF layout (need chunky 8-bit RGB)")
        self.predictor = g(317, 1)
        if self.predictor not in (1, 2):
            raise WatchError(f"unsupported TIFF predictor {self.predictor}")
        self.offsets, self.counts = self.values(324), self.values(325)
        scale, tie = self.values(33550), self.values(33922)
        if len(scale) < 2 or len(tie) < 6:
            raise WatchError("TIFF has no GeoTIFF pixel scale / tiepoint")
        self.x0, self.y0, self.sx, self.sy = tie[3] - tie[0] * scale[0], tie[4] + tie[1] * scale[1], scale[0], scale[1]
        self.across = (self.width + self.tw - 1) // self.tw

    def _bytes(self, off: int, n: int) -> bytes:
        if off + n <= len(self._head):
            return self._head[off:off + n]
        return self.rng(self.url, off, n)

    SIZES = {1: 1, 2: 1, 3: 2, 4: 4, 5: 8, 12: 8, 16: 8}
    FMT = {1: "B", 2: "c", 3: "H", 4: "I", 12: "d", 16: "Q"}

    def values(self, tag: int) -> list:
        if tag not in self.tags:
            return []
        typ, cnt, raw = self.tags[tag]
        size = self.SIZES.get(typ)
        if size is None or typ not in self.FMT:
            raise WatchError(f"unsupported TIFF field type {typ} in tag {tag}")
        blob = raw[:size * cnt] if size * cnt <= 4 else self._bytes(struct.unpack(self.bo + "I", raw)[0], size * cnt)
        return list(struct.unpack(self.bo + f"{cnt}{self.FMT[typ]}", blob))

    def first(self, tag: int, default=None):
        v = self.values(tag)
        return v[0] if v else default

    def colrow(self, x: float, y: float) -> tuple[float, float]:
        return (x - self.x0) / self.sx, (self.y0 - y) / self.sy

    def window(self, col0: int, row0: int, w: int, h: int) -> bytes:
        """RGB bytes of the w x h pixel window whose top-left pixel is (col0,row0). Raises if it leaves the raster."""
        if col0 < 0 or row0 < 0 or col0 + w > self.width or row0 + h > self.height:
            raise WatchError("window falls outside this scene")
        out = bytearray(w * h * 3)
        for ty in range(row0 // self.th, (row0 + h - 1) // self.th + 1):
            for tx in range(col0 // self.tw, (col0 + w - 1) // self.tw + 1):
                idx = ty * self.across + tx
                raw = zlib.decompress(self._bytes(self.offsets[idx], self.counts[idx]))
                if len(raw) != self.tw * self.th * 3:
                    raise WatchError("tile size mismatch")
                gx0, gy0 = tx * self.tw, ty * self.th
                xa, xb = max(col0, gx0), min(col0 + w, gx0 + self.tw)
                ya, yb = max(row0, gy0), min(row0 + h, gy0 + self.th)
                for gy in range(ya, yb):
                    base = (gy - gy0) * self.tw * 3
                    row = raw[base:base + (xb - gx0) * 3]
                    if self.predictor == 2:
                        chans = [bytes(v & 255 for v in accumulate(row[c::3])) for c in range(3)]
                    else:
                        chans = [row[c::3] for c in range(3)]
                    o = ((gy - row0) * w + (xa - col0)) * 3
                    for c in range(3):
                        out[o + c:o + (xb - xa) * 3:3] = chans[c][xa - gx0:xb - gx0]
        return bytes(out)


def fetch_clip(item: dict, lat: float, lon: float, size_px: int, cog_factory=Cog) -> Image:
    props = item["properties"]
    epsg = props.get("proj:epsg")
    if not epsg and str(props.get("proj:code", "")).startswith("EPSG:"):
        epsg = int(props["proj:code"][5:])
    if not epsg:
        raise WatchError(f"scene {item.get('id')} has no proj:epsg")
    cog = cog_factory(item["assets"]["visual"]["href"])
    x, y = utm_forward(lat, lon, int(epsg))
    c, r = cog.colrow(x, y)
    col0, row0 = int(round(c)) - size_px // 2, int(round(r)) - size_px // 2
    return Image(size_px, size_px, cog.window(col0, row0, size_px, size_px))


# ----------------------------------------------------------------------------- change maths
def is_cloudlike(r: int, g: int, b: int) -> bool:
    """Very bright and nearly colourless: cloud, haze, snow (or a bright white roof)."""
    lo = min(r, g, b)
    return lo >= 190 and max(r, g, b) - lo <= 40


def clip_quality(img: Image) -> dict:
    d, n = img.data, img.w * img.h
    nod = cld = 0
    for i in range(0, len(d), 3):
        r, g, b = d[i], d[i + 1], d[i + 2]
        if r == 0 and g == 0 and b == 0:
            nod += 1
        elif is_cloudlike(r, g, b):
            cld += 1
    return {"nodata": nod / n, "cloudlike": cld / n}


def _median(v: list[int]) -> float:
    v = sorted(v)
    n = len(v)
    return 0.0 if not n else (v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2)


def compare(a: Image, b: Image, pixel_threshold: float = 30.0) -> dict:
    """Pixel-difference score of clip b against clip a.

    1. Mask pixels that are nodata or cloud-like in either clip.
    2. Remove the global brightness/colour shift (median per-channel difference over valid pixels), so a hazier day or a sun-angle
       change does not light up the whole box.
    3. Per pixel: largest absolute channel difference; smooth with a 3x3 mean over valid neighbours (kills single-pixel registration noise).
    4. changed = smoothed >= pixel_threshold.  score = changed / valid pixels."""
    if (a.w, a.h) != (b.w, b.h):
        raise ValueError("clips differ in size")
    w, h, n = a.w, a.h, a.w * a.h
    A, B = a.data, b.data
    valid = [False] * n
    for i in range(n):
        p = i * 3
        ra, ga, ba, rb, gb, bb = A[p], A[p + 1], A[p + 2], B[p], B[p + 1], B[p + 2]
        valid[i] = not ((ra == ga == ba == 0) or (rb == gb == bb == 0) or is_cloudlike(ra, ga, ba) or is_cloudlike(rb, gb, bb))
    nv = sum(valid)
    if nv == 0:
        return {"valid_frac": 0.0, "changed_frac": 0.0, "mean_delta": 0.0, "smooth": [0.0] * n, "valid": valid, "changed": [False] * n}
    shift = [_median([B[i * 3 + c] - A[i * 3 + c] for i in range(n) if valid[i]]) for c in range(3)]
    d = [0.0] * n
    for i in range(n):
        if valid[i]:
            p = i * 3
            d[i] = max(abs(B[p + c] - shift[c] - A[p + c]) for c in range(3))
    smooth = [0.0] * n
    for y in range(h):
        for x in range(w):
            s = k = 0
            for yy in range(max(0, y - 1), min(h, y + 2)):
                for xx in range(max(0, x - 1), min(w, x + 2)):
                    j = yy * w + xx
                    if valid[j]:
                        s += d[j]
                        k += 1
            smooth[y * w + x] = s / k if k else 0.0
    changed = [valid[i] and smooth[i] >= pixel_threshold for i in range(n)]
    return {"valid_frac": nv / n, "changed_frac": sum(changed) / nv, "mean_delta": sum(smooth[i] for i in range(n) if valid[i]) / nv,
            "smooth": smooth, "valid": valid, "changed": changed, "shift": shift}


def heatmap_png(cmp: dict, w: int, h: int, pixel_threshold: float) -> bytes:
    """RGBA overlay: transparent where unchanged or masked, yellow -> red as the difference grows past the threshold."""
    out = bytearray()
    for i in range(w * h):
        if cmp["changed"][i]:
            v = min(1.0, (cmp["smooth"][i] - pixel_threshold) / 90.0)
            out += bytes((255, int(230 * (1 - v)), 0, int(120 + 135 * v)))
        else:
            out += b"\x00\x00\x00\x00"
    return png_encode(w, h, bytes(out), 4)


# ----------------------------------------------------------------------------- the scan
def scan(lat: float, lon: float, months: list[str], size_m: float = 800, pixel_threshold: float = 30.0, flag_frac: float = 0.05,
         max_clip_cloud: float = 0.2, outlier: float = 1.5, tries_per_month: int = 3, search=stac_search, clip=fetch_clip, log=lambda s: None) -> dict:
    """Run the whole comparison. Returns {'result': json-able dict, 'clips': {month: Image}, 'heat': {pairkey: png bytes}}."""
    size_px = max(16, int(round(size_m / PIXEL_M)))
    items = [it for it in search(lat, lon, size_m, months) if window_fits(it, lat, lon, size_px)]
    if not items:
        raise WatchError("no Sentinel-2 scene fully contains the box (it may sit on a tile edge or the window has no scenes); try a smaller --size or a nearby centre")
    tile = pick_tile(items, months)
    cands = candidates_by_month(items, months, tile)
    rows, clips = [], {}
    for m in months:
        row = {"month": m, "status": "no-scene", "scene": None, "date": None, "scene_cloud": None, "clip_cloudlike": None, "tried": 0, "note": ""}
        if not cands[m]:
            row["note"] = "no Sentinel-2 L2A scene over the box in this month"
        for it in cands[m][:tries_per_month]:
            row["tried"] += 1
            try:
                img = clip(it, lat, lon, size_px)
            except WatchError as e:
                row.update(status="error", note=f"{it['id']}: {e}")
                log(f"  {m}: {it['id']}: {e}")
                continue
            q = clip_quality(img)
            if q["nodata"] > 0.2 or q["cloudlike"] > max_clip_cloud:
                row.update(status="cloudy", note=f"best tries were cloudy/nodata (last: {it['id']} cloudlike {q['cloudlike']:.0%}, nodata {q['nodata']:.0%})",
                           scene_cloud=it["properties"].get("eo:cloud_cover"))
                continue
            row.update(status="ok", scene=it["id"], date=it["properties"]["datetime"][:10], scene_cloud=it["properties"].get("eo:cloud_cover"),
                       clip_cloudlike=round(q["cloudlike"], 4), note="")
            clips[m] = img
            break
        log(f"  {m}: {row['status']} {row['scene'] or row['note']}")
        rows.append(row)
    ok = [r for r in rows if r["status"] == "ok"]
    heat: dict[str, bytes] = {}

    def pair(a_m: str, b_m: str) -> dict:
        c = compare(clips[a_m], clips[b_m], pixel_threshold)
        key = f"{a_m}..{b_m}"
        heat[key] = heatmap_png(c, size_px, size_px, pixel_threshold)
        usable = c["valid_frac"] >= 0.5
        return {"against": a_m, "changed_frac": round(c["changed_frac"], 4), "valid_frac": round(c["valid_frac"], 4), "mean_delta": round(c["mean_delta"], 2),
                "usable": usable, "flagged": False, "heat": key}

    for i, r in enumerate(ok):
        r["vs_prev"] = pair(ok[i - 1]["month"], r["month"]) if i else None
        r["vs_first"] = pair(ok[0]["month"], r["month"]) if i > 1 else (r["vs_prev"] if i == 1 else None)
    # A month is flagged when its change against the previous clip is above the absolute floor AND clearly above this place's own usual
    # month-to-month churn (crops, shadows and water move every month). The median is taken over all usable consecutive pairs.
    fr = [r["vs_prev"]["changed_frac"] for r in ok if r.get("vs_prev") and r["vs_prev"]["usable"]]
    typical = _median(fr) if fr else 0.0
    for r in ok:
        p = r.get("vs_prev")
        if p and p["usable"]:
            p["flagged"] = bool(p["changed_frac"] >= flag_frac and p["changed_frac"] >= outlier * typical)
    flagged = [r["month"] for r in ok if r.get("vs_prev") and r["vs_prev"]["flagged"]]
    summary = {"months_with_clip": len(ok), "months_requested": len(months), "flagged_months": flagged,
               "first_vs_last": ok[-1]["vs_first"] if len(ok) >= 2 else None, "errors": [r["month"] for r in rows if r["status"] == "error"]}
    result = {"tool": "property-watch", "version": __version__, "tile": tile, "size_m": size_m, "size_px": size_px, "pixel_m": PIXEL_M,
              "lat": lat, "lon": lon, "months": rows, "params": {"pixel_threshold": pixel_threshold, "flag_frac": flag_frac, "max_clip_cloud": max_clip_cloud},
              "summary": summary}
    result["params"]["outlier"] = outlier
    summary["typical_changed_frac"] = round(typical, 4)
    return {"result": result, "clips": clips, "heat": heat}


# ----------------------------------------------------------------------------- HTML report
def render_report(out: dict, title: str) -> str:
    res, clips, heat = out["result"], out["clips"], out["heat"]
    ok = [r for r in res["months"] if r["status"] == "ok"]
    uri = lambda b: "data:image/png;base64," + base64.b64encode(b).decode()
    data = {"clips": {m: uri(img.png()) for m, img in clips.items()}, "heat": {k: uri(v) for k, v in heat.items()},
            "pairs": []}
    pairs = []
    for r in ok:
        if r.get("vs_prev"):
            pairs.append((r["vs_prev"]["against"], r["month"], r["vs_prev"]))
    if len(ok) >= 3:
        pairs.append((ok[0]["month"], ok[-1]["month"], ok[-1]["vs_first"]))
    for a, b, p in pairs:
        data["pairs"].append({"a": a, "b": b, "key": p["heat"], "label": f"{a} -> {b}: {p['changed_frac']:.1%} changed" + (" (flagged)" if p["flagged"] else "") +
                              ("" if p["usable"] else " (unreliable: too much cloud/nodata)")})
    blob = json.dumps(data).replace("</", "<\\/")
    esc = html.escape
    trs = []
    for r in res["months"]:
        vp, vf = r.get("vs_prev"), r.get("vs_first")
        trs.append("<tr class='%s'><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>" % (
            "flag" if vp and vp["flagged"] else r["status"], esc(r["month"]), esc(r["status"]), esc(r["date"] or ""), esc(r["scene"] or esc(r["note"])),
            "" if r["scene_cloud"] is None else f"{r['scene_cloud']:.0f}%",
            "" if not vp else f"{vp['changed_frac']:.1%} vs {vp['against']}" + (" FLAG" if vp["flagged"] else "") + ("" if vp["usable"] else " (unreliable)"),
            "" if not vf else f"{vf['changed_frac']:.1%}"))
    fl = res["summary"]["flagged_months"]
    fvl = res["summary"]["first_vs_last"]
    head = ("Months flagged for significant change: <b>%s</b>." % esc(", ".join(fl)) if fl else "No month crossed the change threshold.") + (
        " First clip vs last clip: <b>%.1f%%</b> of valid pixels changed." % (fvl["changed_frac"] * 100) if fvl else "")
    p = res["params"]
    return f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>property-watch: {esc(title)}</title><style>
body{{font:15px/1.5 system-ui,sans-serif;max-width:860px;margin:1.5rem auto;padding:0 1rem;color:#1b1f23}}
h1{{font-size:1.35rem}} table{{border-collapse:collapse;width:100%;font-size:13px}} td,th{{border:1px solid #d0d7de;padding:3px 6px;text-align:left}}
tr.flag{{background:#ffe3e0}} tr.cloudy,tr.no-scene,tr.error{{color:#777}} .stage{{position:relative;width:min(100%,520px);aspect-ratio:1;margin:.6rem 0;border:1px solid #888;background:#222;overflow:hidden}}
.stage img{{position:absolute;inset:0;width:100%;height:100%;image-rendering:pixelated}} #top{{clip-path:inset(0 50% 0 0)}} #heat{{display:none;pointer-events:none}}
#line{{position:absolute;top:0;bottom:0;width:2px;background:#fff;left:50%;box-shadow:0 0 3px #000;pointer-events:none}}
.lab{{position:absolute;top:4px;font:12px monospace;background:#000a;color:#fff;padding:1px 5px}} .note{{background:#f6f8fa;border:1px solid #d0d7de;padding:.6rem .8rem;font-size:13px}}
input[type=range]{{width:min(100%,520px)}}</style></head><body>
<h1>property-watch: {esc(title)}</h1>
<p>{esc('%.5f, %.5f' % (res['lat'], res['lon']))} &middot; {res['size_m']:.0f} m box ({res['size_px']}&times;{res['size_px']} px at {res['pixel_m']:.0f} m) &middot; {esc(res['months'][0]['month'])} to {esc(res['months'][-1]['month'])} &middot; tile {esc(str(res['tile']))} &middot; Sentinel-2 L2A</p>
<p>{head}</p>
<div id="ui" style="display:none"><label>Compare <select id="pair"></select></label> <label><input type="checkbox" id="hm"> change heatmap</label>
<div class="stage" id="stage"><img id="bot" alt="before"><img id="top" alt="after"><img id="heat" alt="change heatmap"><div id="line"></div>
<span class="lab" id="la" style="left:4px"></span><span class="lab" id="lb" style="right:4px"></span></div>
<input type="range" id="sl" min="0" max="100" value="50" aria-label="before/after slider"></div>
<p id="none" style="display:none">Fewer than two usable monthly clips, so there is nothing to compare.</p>
<h2 style="font-size:1.05rem">Months</h2>
<table><tr><th>month</th><th>status</th><th>scene date</th><th>scene id / note</th><th>scene cloud</th><th>changed vs previous clip</th><th>changed vs first clip</th></tr>{''.join(trs)}</table>
<h2 style="font-size:1.05rem">How to read this</h2>
<div class="note">Left of the slider handle is the earlier month, right is the later one. A pixel counts as changed when its colour differs by at least {p['pixel_threshold']:.0f}/255 (largest RGB channel, 3&times;3 smoothed, after removing the box-wide brightness shift), ignoring cloud-like and no-data pixels.
A month is flagged when at least {p['flag_frac']:.0%} of the valid pixels changed against the previous usable clip and that is at least {p['outlier']:.1f}&times; this place's typical month-to-month change ({res['summary']['typical_changed_frac']:.1%} here). Sentinel-2 pixels are {res['pixel_m']:.0f} m: this finds large construction, clearing, new roofs or paving and big ground changes.
It cannot see objects a few pixels wide, so it says nothing about a pool, a shed or a fence. Crops, trees, water level, snow and shadows change the colour too: a flag is a prompt to look, not proof of a building change. Generated by property-watch {__version__}; imagery: contains modified Copernicus Sentinel data, via Element84 Earth Search.</div>
<script>var D=JSON.parse({json.dumps(blob)});</script>
<script>
(function(){{var sel=document.getElementById('pair'),bot=document.getElementById('bot'),top=document.getElementById('top'),heat=document.getElementById('heat'),
hm=document.getElementById('hm'),sl=document.getElementById('sl'),line=document.getElementById('line');
if(D.pairs.length<1){{document.getElementById('none').style.display='block';return}}
document.getElementById('ui').style.display='block';
D.pairs.forEach(function(p,i){{var o=document.createElement('option');o.value=i;o.textContent=p.label;sel.appendChild(o)}});
function show(){{var p=D.pairs[sel.value];bot.src=D.clips[p.a];top.src=D.clips[p.b];heat.src=D.heat[p.key];
document.getElementById('la').textContent=p.a;document.getElementById('lb').textContent=p.b;heat.style.display=hm.checked?'block':'none'}}
function slide(){{top.style.clipPath='inset(0 0 0 '+sl.value+'%)';line.style.left=sl.value+'%'}}
sel.onchange=show;hm.onchange=show;sl.oninput=slide;
var f=D.pairs.length-1;sel.value=f;show();slide();}})();
</script></body></html>
"""


# ----------------------------------------------------------------------------- CLI
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="property-watch", description="Find large changes at a place over a 12-month window with free Sentinel-2 imagery.")
    ap.add_argument("place", nargs="?", help="address / place name (geocoded with Nominatim) or 'lat,lon'")
    ap.add_argument("--lat", type=float)
    ap.add_argument("--lon", type=float)
    ap.add_argument("--end", default=None, help="last month of the window, YYYY-MM (default: last full month)")
    ap.add_argument("--months", type=int, default=12, help="window length in months (default 12)")
    ap.add_argument("--size", type=float, default=800, help="box edge in metres (default 800 = 80x80 px)")
    ap.add_argument("--threshold", type=float, default=30, help="pixel difference (0-255) that counts as changed (default 30)")
    ap.add_argument("--flag", type=float, default=0.05, help="fraction of valid pixels changed that flags a month (default 0.05)")
    ap.add_argument("--outlier", type=float, default=1.5, help="a flagged month must also be this many times the place's median month-to-month change (default 1.5)")
    ap.add_argument("--max-cloud", type=float, default=0.2, help="reject a clip with more than this cloud-like fraction (default 0.2)")
    ap.add_argument("--out", default="property-watch-report.html", help="HTML report path")
    ap.add_argument("--json", default=None, help="also write the numbers as JSON")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)
    log = (lambda s: None) if a.quiet else (lambda s: print(s, file=sys.stderr))
    try:
        if a.lat is not None and a.lon is not None:
            lat, lon, name = a.lat, a.lon, f"{a.lat:.5f}, {a.lon:.5f}"
        elif a.place and parse_latlon(a.place):
            lat, lon = parse_latlon(a.place)
            name = f"{lat:.5f}, {lon:.5f}"
        elif a.place:
            lat, lon, name = geocode(a.place)
            log(f"geocoded {a.place!r} -> {lat:.5f}, {lon:.5f} ({name})")
        else:
            ap.error("give a place or --lat and --lon")
            return 2
        if not (1 <= a.months <= 36) or a.size < 200 or a.size > 5000:
            raise WatchError("--months must be 1-36 and --size 200-5000 m")
        months = month_list(a.end or default_end(), a.months)
        log(f"window {months[0]}..{months[-1]}, box {a.size:.0f} m")
        out = scan(lat, lon, months, a.size, a.threshold, a.flag, a.max_cloud, a.outlier, log=log)
        res = out["result"]
        res["query"] = {"place": a.place, "display_name": name, "lat": lat, "lon": lon}
        Path(a.out).write_text(render_report(out, a.place or name), encoding="utf-8")
        if a.json:
            Path(a.json).write_text(json.dumps(res, indent=2), encoding="utf-8")
    except WatchError as e:
        print(f"FAILED: {e}", file=sys.stderr)
        return 1
    s = res["summary"]
    print(f"wrote {a.out}: {s['months_with_clip']}/{s['months_requested']} months with a clip, flagged: {', '.join(s['flagged_months']) or 'none'}")
    if s["months_with_clip"] < 2:
        print("FAILED: fewer than two usable monthly clips, nothing to compare", file=sys.stderr)
        return 1
    if s["errors"]:
        print(f"FAILED: {len(s['errors'])} month(s) could not be read ({', '.join(s['errors'])}); see the report table", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
