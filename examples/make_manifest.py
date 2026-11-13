#!/usr/bin/env python3
"""Rebuild examples/MANIFEST.json (sha256 + size of each saved-run file, the exact command, the numbers) from examples/*/result.json."""
import hashlib, json, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import property_watch as pw  # noqa: E402

RUNS = {
    "richland-parish": "python3 property_watch.py --lat 32.5289 --lon -91.6361 --size 2000 --end 2026-10 --out examples/richland-parish/report.html --json examples/richland-parish/result.json",
}
out = {"tool": "property-watch", "version": pw.__version__, "python": sys.version.split()[0], "runs": []}
for name, cmd in RUNS.items():
    res = json.loads((HERE / name / "result.json").read_text())
    files = {f.name: {"bytes": f.stat().st_size, "sha256": hashlib.sha256(f.read_bytes()).hexdigest()} for f in sorted((HERE / name).iterdir()) if f.is_file()}
    out["runs"].append({
        "name": name, "command": cmd, "window": [res["months"][0]["month"], res["months"][-1]["month"]], "tile": res["tile"],
        "box_m": res["size_m"], "pixels": res["size_px"], "scenes": {m["month"]: m["scene"] for m in res["months"]},
        "months_without_clip": {m["month"]: m["status"] for m in res["months"] if m["status"] != "ok"},
        "flagged_months": res["summary"]["flagged_months"], "typical_month_to_month_changed_frac": res["summary"]["typical_changed_frac"],
        "first_vs_last_changed_frac": (res["summary"]["first_vs_last"] or {}).get("changed_frac"), "files": files,
        "imagery": "Contains modified Copernicus Sentinel data, via Element84 Earth Search (sentinel-2-l2a). Re-running later can pick different scenes if new ones are added.",
    })
(HERE / "MANIFEST.json").write_text(json.dumps(out, indent=2) + "\n")
print("wrote examples/MANIFEST.json")
