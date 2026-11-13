"""Offline tests: no network. Synthetic arrays for the diff logic, a mocked STAC response, a synthetic tiled GeoTIFF for the range reader."""
import json
import struct
import sys
import zlib
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import property_watch as pw  # noqa: E402


# ------------------------------------------------------------------ helpers
def solid(w, h, rgb):
    return pw.Image(w, h, bytes(rgb) * (w * h))


def paint(img, x0, y0, x1, y1, rgb):
    d = bytearray(img.data)
    for y in range(y0, y1):
        for x in range(x0, x1):
            d[(y * img.w + x) * 3:(y * img.w + x) * 3 + 3] = bytes(rgb)
    return pw.Image(img.w, img.h, bytes(d))


def textured(w, h):
    """Deterministic mid-tone texture (not cloud-like, not nodata)."""
    d = bytearray()
    for y in range(h):
        for x in range(w):
            d += bytes((60 + (x * 3 + y) % 40, 90 + (x + y * 2) % 40, 50 + (x * y) % 30))
    return pw.Image(w, h, bytes(d))


def build_tiff(width, height, tw, th, pixels, x0, y0, scale=10.0, predictor=2):
    """Tiled little-endian deflate RGB GeoTIFF. pixels: bytes row-major RGB, width*height*3."""
    across, down = (width + tw - 1) // tw, (height + th - 1) // th
    tiles = []
    for ty in range(down):
        for tx in range(across):
            rows = bytearray()
            for y in range(ty * th, ty * th + th):
                row = bytearray()
                for x in range(tx * tw, tx * tw + tw):
                    row += pixels[(y * width + x) * 3:(y * width + x) * 3 + 3] if x < width and y < height else b"\x00\x00\x00"
                if predictor == 2:
                    row = bytearray((row[i] - (row[i - 3] if i >= 3 else 0)) & 255 for i in range(len(row)))
                rows += row
            tiles.append(zlib.compress(bytes(rows)))
    buf = bytearray(b"II*\x00\x08\x00\x00\x00")
    ifd_pos = 8
    n_tags = 14
    buf += b"\x00" * (2 + n_tags * 12 + 4)
    extra = {}

    def put(blob):
        off = len(buf)
        buf.extend(blob)
        return off

    offsets, counts = [], []
    for t in tiles:
        offsets.append(put(t))
        counts.append(len(t))
    bits = put(struct.pack("<3H", 8, 8, 8))
    off_arr = put(struct.pack(f"<{len(offsets)}I", *offsets))
    cnt_arr = put(struct.pack(f"<{len(counts)}I", *counts))
    sc = put(struct.pack("<3d", scale, scale, 0))
    tp = put(struct.pack("<6d", 0, 0, 0, x0, y0, 0))
    short = lambda v: struct.pack("<H", v) + b"\x00\x00"
    long_ = lambda v: struct.pack("<I", v)
    ent = [(256, 4, 1, long_(width)), (257, 4, 1, long_(height)), (258, 3, 3, struct.pack("<I", bits)), (259, 3, 1, short(8)),
           (262, 3, 1, short(2)), (277, 3, 1, short(3)), (284, 3, 1, short(1)), (317, 3, 1, short(predictor)), (322, 3, 1, short(tw)),
           (323, 3, 1, short(th)), (324, 4, len(offsets), struct.pack("<I", off_arr)), (325, 4, len(counts), struct.pack("<I", cnt_arr)),
           (33550, 12, 3, struct.pack("<I", sc)), (33922, 12, 6, struct.pack("<I", tp))]
    assert len(ent) == n_tags
    ent.sort()
    body = struct.pack("<H", n_tags) + b"".join(struct.pack("<HHI", t, ty, c) + v for t, ty, c, v in ent) + b"\x00\x00\x00\x00"
    buf[ifd_pos:ifd_pos + len(body)] = body
    return bytes(buf)


def fake_range(blob, log=None):
    def rng(url, start, length):
        if log is not None:
            log.append((start, length))
        return blob[start:start + length]
    return rng


# ------------------------------------------------------------------ png / projection / inputs
def test_png_roundtrip_rgb_and_rgba():
    img = textured(13, 7)
    assert pw.Image.from_png(img.png()).data == img.data
    rgba = bytes(range(40)) * 3
    w, h, ch, d = pw.png_decode(pw.png_encode(5, 6, rgba[:120], 4))
    assert (w, h, ch, d) == (5, 6, 4, rgba[:120])


def test_png_bad_inputs():
    with pytest.raises(ValueError):
        pw.png_encode(2, 2, b"\x00" * 5, 3)
    with pytest.raises(pw.WatchError):
        pw.png_decode(b"nope")


def test_utm_known_values():
    assert pw.utm_forward(0, 3, 32631) == pytest.approx((500000.0, 0.0), abs=1e-6)
    e, n = pw.utm_forward(36.0157189, -114.7374858, 32611)   # reference values from pyproj
    assert (e, n) == pytest.approx((703891.0416, 3988059.7435), abs=0.01)
    e, n = pw.utm_forward(-33.8688, 151.2093, 32756)          # southern hemisphere false northing
    assert (e, n) == pytest.approx((334368.6336, 6250948.3454), abs=0.01)
    with pytest.raises(pw.WatchError):
        pw.utm_forward(0, 0, 3857)


def test_months_and_inputs():
    assert pw.month_list("2026-10", 12)[0] == "2025-11" and pw.month_list("2026-10", 12)[-1] == "2026-10"
    assert pw.month_list("2026-02", 3) == ["2025-12", "2026-01", "2026-02"]
    import datetime as dt
    assert pw.default_end(dt.date(2026, 1, 5)) == "2025-12" and pw.default_end(dt.date(2026, 11, 13)) == "2026-10"
    assert pw.parse_latlon("36.01, -114.73") == (36.01, -114.73)
    assert pw.parse_latlon("10 Downing Street London") is None
    with pytest.raises(pw.WatchError):
        pw.parse_latlon("95.0,10.0")
    b = pw.bbox_around(0.0, 10.0, 1000)
    assert b[2] - b[0] == pytest.approx(b[3] - b[1], rel=1e-3) and b[0] < 10 < b[2]


# ------------------------------------------------------------------ diff logic (synthetic arrays)
def test_identical_clips_have_no_change():
    a = textured(40, 40)
    c = pw.compare(a, a)
    assert c["changed_frac"] == 0 and c["valid_frac"] == 1


def test_global_brightness_shift_is_removed():
    a = textured(40, 40)
    b = pw.Image(40, 40, bytes(min(255, v + 20) for v in a.data))
    assert pw.compare(a, b)["changed_frac"] == 0


def test_block_change_is_found_with_right_size():
    a = textured(50, 50)
    b = paint(a, 10, 10, 30, 30, (200, 60, 40))   # 20x20 = 400 of 2500 px = 16 %
    c = pw.compare(a, b)
    assert 0.14 <= c["changed_frac"] <= 0.24   # 3x3 smoothing grows the block by about a pixel each side
    assert c["changed"][20 * 50 + 20] and not c["changed"][2 * 50 + 2]


def test_single_pixel_noise_is_smoothed_away():
    a = textured(40, 40)
    b = paint(a, 20, 20, 21, 21, (250, 0, 0))
    assert pw.compare(a, b, pixel_threshold=30)["changed_frac"] == 0   # 1 px of ~190 over a 3x3 mean = ~21 < 30


def test_cloud_and_nodata_pixels_are_masked():
    a = textured(40, 40)
    cloud = paint(a, 0, 0, 20, 40, (240, 240, 240))
    nod = paint(a, 20, 0, 40, 20, (0, 0, 0))
    c = pw.compare(a, cloud)
    assert c["valid_frac"] == pytest.approx(0.5) and c["changed_frac"] == 0
    assert pw.compare(a, nod)["valid_frac"] == pytest.approx(0.75)
    assert pw.compare(solid(10, 10, (240, 240, 240)), solid(10, 10, (240, 240, 240)))["valid_frac"] == 0


def test_clip_quality_and_cloudlike():
    img = paint(textured(10, 10), 0, 0, 5, 10, (230, 232, 235))
    img = paint(img, 5, 0, 10, 2, (0, 0, 0))
    q = pw.clip_quality(img)
    assert q["cloudlike"] == pytest.approx(0.5) and q["nodata"] == pytest.approx(0.1)
    assert pw.is_cloudlike(200, 200, 210) and not pw.is_cloudlike(230, 120, 100) and not pw.is_cloudlike(150, 150, 150)


def test_heatmap_is_transparent_where_unchanged():
    a = textured(20, 20)
    b = paint(a, 5, 5, 15, 15, (220, 40, 40))
    c = pw.compare(a, b)
    w, h, ch, d = pw.png_decode(pw.heatmap_png(c, 20, 20, 30))
    assert ch == 4 and d[3] == 0 and d[(10 * 20 + 10) * 4 + 3] > 100


# ------------------------------------------------------------------ mocked STAC
def item(id_, dt_, cloud, grid="MGRS-11SQV", **assets):
    return {"id": id_, "properties": {"datetime": dt_ + "T18:00:00Z", "eo:cloud_cover": cloud, "grid:code": grid, "proj:epsg": 32611},
            "assets": {"visual": {"href": f"https://example.invalid/{id_}/TCI.tif"}}}


def test_stac_search_follows_next_pages():
    pages = [{"features": [item("a", "2026-01-03", 5)], "links": [{"rel": "next", "body": {"page": 2}}]},
             {"features": [item("b", "2026-02-03", 9)], "links": []}]
    seen = []

    def post(url, body):
        seen.append(body)
        return pages[len(seen) - 1]

    got = pw.stac_search(10.0, 20.0, 800, ["2026-01", "2026-02"], post=post)
    assert [g["id"] for g in got] == ["a", "b"] and len(seen) == 2
    assert seen[0]["collections"] == ["sentinel-2-l2a"] and seen[0]["datetime"] == "2026-01-01T00:00:00Z/2026-02-28T23:59:59Z"
    assert len(seen[0]["bbox"]) == 4


def test_pick_tile_and_clearest_scene_per_month():
    items = [item("a1", "2026-01-03", 40), item("a2", "2026-01-13", 3), item("a3", "2026-02-03", 8),
             item("b1", "2026-01-05", 1, grid="MGRS-11SQU")]
    months = ["2026-01", "2026-02"]
    assert pw.pick_tile(items, months) == "MGRS-11SQV"   # covers both months
    c = pw.candidates_by_month(items, months, "MGRS-11SQV")
    assert [i["id"] for i in c["2026-01"]] == ["a2", "a1"] and [i["id"] for i in c["2026-02"]] == ["a3"]
    assert pw.pick_tile([], months) is None


def test_window_fits_uses_metadata():
    lat, lon = 36.0157, -114.7375
    e, n = pw.utm_forward(lat, lon, 32611)
    it = item("x", "2026-01-01", 1)
    it["assets"]["visual"]["proj:shape"] = [10980, 10980]
    it["assets"]["visual"]["proj:transform"] = [10, 0, e - 5000, 0, -10, n + 5000]
    assert pw.window_fits(it, lat, lon, 80)
    it["assets"]["visual"]["proj:transform"] = [10, 0, e - 100, 0, -10, n + 5000]   # point 10 px from the left edge, box needs 40
    assert not pw.window_fits(it, lat, lon, 80)
    assert pw.window_fits({"properties": {}, "assets": {}}, lat, lon, 80)   # no metadata: reader decides


# ------------------------------------------------------------------ range reader on a synthetic COG
@pytest.mark.parametrize("predictor", [1, 2])
def test_cog_window_matches_source_pixels_across_tiles(predictor):
    W, H = 70, 50
    src = bytes((x * 7 + y * 13 + c * 29) % 256 for y in range(H) for x in range(W) for c in range(3))
    blob = build_tiff(W, H, 32, 16, src, 500000.0, 4000000.0, predictor=predictor)
    log = []
    cog = pw.Cog("u", rng=fake_range(blob, log))
    assert (cog.width, cog.height, cog.tw, cog.th, cog.x0, cog.y0) == (70, 50, 32, 16, 500000.0, 4000000.0)
    assert cog.colrow(500100.0, 3999800.0) == (10.0, 20.0)
    col0, row0, w, h = 25, 10, 30, 30   # spans 2x3 tiles
    got = cog.window(col0, row0, w, h)
    want = b"".join(src[((row0 + y) * W + col0) * 3:((row0 + y) * W + col0 + w) * 3] for y in range(h))
    assert got == want
    with pytest.raises(pw.WatchError):
        cog.window(60, 0, 20, 10)


def test_cog_rejects_unsupported():
    with pytest.raises(pw.WatchError):
        pw.Cog("u", rng=fake_range(b"RIFFxxxxxxxxxxxxxxxx" * 4))
    with pytest.raises(pw.WatchError):
        pw.Cog("u", rng=fake_range(b"II+\x00" + b"\x00" * 100))   # BigTIFF


def test_fetch_clip_centres_on_the_point():
    lat, lon = 36.0157, -114.7375
    e, n = pw.utm_forward(lat, lon, 32611)
    W = H = 64
    x0, y0 = e - 320.0, n + 320.0   # point at pixel (32, 32)
    src = bytes((x * 4 + c) % 256 if (x, y) != (32, 32) else 255 for y in range(H) for x in range(W) for c in range(3))
    blob = build_tiff(W, H, 32, 32, src, x0, y0)
    it = item("s", "2026-01-01", 1)
    img = pw.fetch_clip(it, lat, lon, 16, cog_factory=lambda url: pw.Cog(url, rng=fake_range(blob)))
    assert (img.w, img.h) == (16, 16)
    cx = 32 - 8                                      # clip column 0 = source column 24
    assert img.data[0] == (cx * 4) % 256
    assert img.data[(8 * 16 + 8) * 3] == 255          # the marked pixel sits at the clip centre


# ------------------------------------------------------------------ whole scan with fakes
def fake_world(frames):
    """frames: {month: Image}. Search returns one scene per month; clip returns that month's image."""
    items = [item(f"S_{m}", f"{m}-10", 2) for m in frames]
    return (lambda lat, lon, size, months: items), (lambda it, lat, lon, px: frames[it["properties"]["datetime"][:7]])


def test_scan_flags_the_month_with_unusual_change():
    base = textured(40, 40)
    frames = {"2026-01": base, "2026-02": base, "2026-03": base, "2026-04": paint(base, 5, 5, 35, 35, (210, 60, 40)), "2026-05": paint(base, 5, 5, 35, 35, (210, 60, 40))}
    search, clip = fake_world(frames)
    out = pw.scan(10, 20, list(frames), 400, search=search, clip=clip)
    s = out["result"]["summary"]
    assert s["flagged_months"] == ["2026-04"]
    assert s["first_vs_last"]["changed_frac"] > 0.3 and s["months_with_clip"] == 5
    assert "2026-01..2026-04" in out["heat"] and "2026-01..2026-05" in out["heat"]


def test_scan_skips_cloudy_months_and_reports_missing():
    base = textured(40, 40)
    frames = {"2026-01": base, "2026-02": solid(40, 40, (240, 240, 240)), "2026-04": base}
    search, clip = fake_world(frames)
    out = pw.scan(10, 20, ["2026-01", "2026-02", "2026-03", "2026-04"], 400, search=search, clip=clip)
    st = {r["month"]: r["status"] for r in out["result"]["months"]}
    assert st == {"2026-01": "ok", "2026-02": "cloudy", "2026-03": "no-scene", "2026-04": "ok"}
    assert out["result"]["summary"]["flagged_months"] == []
    assert out["result"]["months"][3]["vs_prev"]["against"] == "2026-01"


def test_scan_falls_back_to_next_candidate_when_best_is_cloudy():
    base = textured(40, 40)
    best, second = item("best", "2026-01-05", 1), item("second", "2026-01-15", 6)
    clips = {"best": solid(40, 40, (240, 240, 240)), "second": base}
    out = pw.scan(10, 20, ["2026-01"], 400, search=lambda *a: [best, second], clip=lambda it, *a: clips[it["id"]])
    r = out["result"]["months"][0]
    assert r["status"] == "ok" and r["scene"] == "second" and r["tried"] == 2


def test_scan_records_clip_errors_with_cause():
    def clip(it, *a):
        raise pw.WatchError("boom")
    out = pw.scan(10, 20, ["2026-01"], 400, search=lambda *a: [item("x", "2026-01-05", 1)], clip=clip)
    r = out["result"]["months"][0]
    assert r["status"] == "error" and "boom" in r["note"] and out["result"]["summary"]["errors"] == ["2026-01"]


def test_scan_fails_when_no_scene_contains_box():
    far = item("far", "2026-01-05", 1)
    far["assets"]["visual"].update({"proj:transform": [10, 0, 0, 0, -10, 0], "proj:shape": [10980, 10980]})
    with pytest.raises(pw.WatchError):
        pw.scan(36.0, -114.7, ["2026-01"], 800, search=lambda *a: [far], clip=lambda *a: None)


# ------------------------------------------------------------------ report + CLI
def run_cli(monkeypatch, tmp_path, frames, extra=()):
    search, clip = fake_world(frames)
    real = pw.scan
    monkeypatch.setattr(pw, "scan", lambda *a, **k: real(*a, search=search, clip=clip, **k))
    out = tmp_path / "r.html"
    rc = pw.main(["10.5,20.5", "--end", "2026-03", "--months", "3", "--size", "400", "--out", str(out), "--json", str(tmp_path / "r.json"), "--quiet", *extra])
    return rc, out


def test_cli_writes_report_and_json(monkeypatch, tmp_path, capsys):
    base = textured(40, 40)
    rc, out = run_cli(monkeypatch, tmp_path, {"2026-01": base, "2026-02": base, "2026-03": paint(base, 0, 0, 40, 40, (200, 50, 50))})
    assert rc == 0
    page = out.read_text()
    assert "<script src" not in page and "http://" not in page.split("</style>")[1].split("<h2")[0] and "cannot see" in page
    js = json.loads((tmp_path / "r.json").read_text())
    assert js["summary"]["months_with_clip"] == 3 and js["query"]["lat"] == 10.5 and js["tool"] == "property-watch"
    assert "wrote" in capsys.readouterr().out


def test_cli_fails_loudly_when_search_fails(monkeypatch, tmp_path, capsys):
    def boom(*a, **k):
        raise pw.WatchError("STAC down")
    monkeypatch.setattr(pw, "scan", boom)
    rc = pw.main(["10.5,20.5", "--out", str(tmp_path / "r.html"), "--quiet"])
    assert rc == 1 and "FAILED: STAC down" in capsys.readouterr().err and not (tmp_path / "r.html").exists()


def test_cli_fails_when_fewer_than_two_clips(monkeypatch, tmp_path, capsys):
    rc, _ = run_cli(monkeypatch, tmp_path, {"2026-03": textured(40, 40)})
    assert rc == 1 and "FAILED" in capsys.readouterr().err


def test_cli_nonzero_when_a_month_errors(monkeypatch, tmp_path, capsys):
    base = textured(40, 40)
    items = [item(f"S_{m}", f"{m}-10", 2) for m in ("2026-01", "2026-02", "2026-03")]

    def clip(it, *a):
        if it["id"] == "S_2026-02":
            raise pw.WatchError("range read failed HTTP 503")
        return base
    real = pw.scan
    monkeypatch.setattr(pw, "scan", lambda *a, **k: real(*a, search=lambda *x: items, clip=clip, **k))
    rc = pw.main(["10.5,20.5", "--end", "2026-03", "--months", "3", "--size", "400", "--out", str(tmp_path / "r.html"), "--quiet"])
    assert rc == 2 and "FAILED: 1 month(s)" in capsys.readouterr().err


def test_cli_geocode_failure_is_failed_line(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(pw, "http_get", lambda *a, **k: type("R", (), {"status_code": 200, "json": lambda self: []})())
    rc = pw.main(["Nowhereville 99999", "--out", str(tmp_path / "r.html"), "--quiet"])
    assert rc == 1 and "FAILED: Nominatim found nothing" in capsys.readouterr().err


def test_geocode_parses_nominatim(monkeypatch):
    monkeypatch.setattr(pw, "http_get", lambda *a, **k: type("R", (), {"status_code": 200, "json": lambda self: [{"lat": "1.5", "lon": "-2.5", "display_name": "X"}]})())
    assert pw.geocode("x") == (1.5, -2.5, "X")


def test_report_escapes_the_place_name():
    base = textured(40, 40)
    search, clip = fake_world({"2026-01": base, "2026-02": base})
    out = pw.scan(10, 20, ["2026-01", "2026-02"], 400, search=search, clip=clip)
    page = pw.render_report(out, "<script>alert(1)</script>")
    assert "<script>alert(1)</script>" not in page and "&lt;script&gt;alert(1)" in page
