# property-watch

Give it an address (or lat/lon) and a 12-month window. It pulls the clearest free Sentinel-2 image of that spot for each month, compares the months pixel by pixel, flags the months where a lot changed, and writes **one self-contained HTML report** with a before/after slider and a change heatmap.

One Python file, standard library plus `requests`. Imagery is the free [Element84 Earth Search](https://element84.com/earth-search/) STAC API (Sentinel-2 L2A on AWS); geocoding is free [Nominatim](https://nominatim.org/). No account, no API key.

## What it can and cannot see (read this first)

Sentinel-2 true-colour pixels are **10 m x 10 m**. So this is good for **large** changes: a new warehouse or data-centre roof, a cleared lot, graded or paved ground, a big new road, a burned or logged patch. It **cannot** see anything a few pixels wide: **a pool, a shed, a fence, a car, a solar panel, a tree removal on one lot.** It does not detect pools and nothing in the report should be read as if it did.

Colour also changes for reasons that are not construction: crops and trees through the seasons, water level, snow, shadows, a hazy day. A flag means "look at this month", not "a building changed". Cloudy clips are rejected and the tool tries the next-clearest scene of that month; if none is clear the month is shown as `cloudy` and skipped.

## Quickstart

```bash
git clone https://github.com/MegaStretch23/property-watch
cd property-watch
python3 -m pip install requests          # or: uv run property_watch.py ...  (PEP 723 header installs it)

# an address (geocoded with Nominatim) ...
python3 property_watch.py "Hoover Dam" --end 2026-10 --out report.html

# ... or coordinates, a 2 km box, plus the numbers as JSON
python3 property_watch.py --lat 32.5289 --lon -91.6361 --size 2000 --end 2026-10 --out report.html --json result.json
```

Open `report.html` in a browser: pick a pair of months, drag the slider, tick "change heatmap". Only your box is fetched from each 110 km scene, with HTTP range requests (measured: about 2 MB per month for a 2 km box, because whole 1024-pixel tiles are downloaded); a 12-month run took under a minute here.

Options: `--months 12` (window length), `--size 800` (box edge in metres), `--threshold 30` (colour difference that counts as changed, 0-255), `--flag 0.05` (share of pixels changed that flags a month), `--outlier 1.5`, `--max-cloud 0.2`. Exit codes: `0` ok, `1` failed (geocoder, STAC or another request failed, or fewer than two usable months), `2` some month could not be read (listed in the output). Failures print a `FAILED:` line on stderr.

## How the comparison works

1. STAC search for Sentinel-2 L2A scenes over the box for each month; keep the one MGRS tile that covers the most months so every clip shares one pixel grid.
2. Per month take the lowest `eo:cloud_cover` scene, read only the box out of its true-colour `TCI.tif` (a cloud-optimised GeoTIFF; the small reader in this file handles deflate tiles), and re-check clouds on the clip itself.
3. For two clips: mask no-data and cloud-like pixels, subtract the box-wide median colour shift (so a hazier day does not light everything up), take the largest RGB difference per pixel, smooth it over 3x3 pixels, call a pixel changed at or above the threshold. Score = changed pixels / valid pixels.
4. A month is **flagged** when its score against the previous usable clip is above `--flag` **and** at least `--outlier` times that place's own median month-to-month score, so farmland that always churns does not flag every month. The report also gives first clip vs last clip, which are about 11 months apart and so roughly the same season.

## Sample output (saved run)

`examples/richland-parish/` is a real run on a public location, the Meta Richland Parish data-centre site in Louisiana (coordinates from Nominatim; a 2 km box, window 2025-11 to 2026-10, Python 3.12.3). Files, scene ids, command and sha256 sums are in [`examples/MANIFEST.json`](examples/MANIFEST.json); open [`examples/richland-parish/report.html`](examples/richland-parish/report.html) after cloning (it is about 2.4 MB).

```
11 of 12 months with a clip (2026-06 was cloudy and skipped)
typical month-to-month change      41.2 % of valid pixels   (surrounding farmland and ground work churn every month)
flagged month                      2026-04 (65.9 %, above the 1.5x typical bar)
first clip (2025-11) vs last (2026-10)   59.8 % of valid pixels changed
```

In the report the November clip and the October clip show visibly different roofs and re-graded ground inside the box. The month-to-month number is high because the site is surrounded by fields and active earthworks; that is why the flag is relative to the place's own typical change and why the first-vs-last figure is the one to read first. Re-running later can choose different scenes if new ones are added, so numbers can differ slightly.

## Tests

`python3 -m pytest` runs 30 offline tests: synthetic arrays for the difference logic (brightness shift, block change, noise, cloud and no-data masks), a mocked paginated STAC response, a synthetic tiled GeoTIFF for the range reader, UTM projection against reference values, and the CLI failure paths. No network needed.

## Related

I also run hosted data tools on Apify and sell small source packs on Gumroad. They are unrelated to satellite imagery; property-watch is free and has no paid version.

- [SEC Form 4 insider cluster scanner](https://apify.com/johnveds15/sec-form4-insider-clusters?utm_source=github&utm_medium=readme&utm_campaign=property-watch) on Apify: pay per event, $0.0004 per Form 4 filing parsed and $0.02 per insider cluster detected (prices read from the Apify API on 2026-11-13). Its free open-source core is [sec-form4-clusters](https://github.com/MegaStretch23/sec-form4-clusters).
- [pdfguard](https://vedovino.gumroad.com/l/gvvcsh?utm_source=github&utm_medium=readme&utm_campaign=property-watch): finds hidden prompt-injection text in PDFs before your LLM reads them (Python, local), Gumroad source pack, $19.
- [blame-ai](https://vedovino.gumroad.com/l/smhyy?utm_source=github&utm_medium=readme&utm_campaign=property-watch): AI-assist ledger for your git history (Python, local, no dependencies), Gumroad source pack, $19.

## Data credits and license

Code: MIT. Imagery: contains modified Copernicus Sentinel data, served by Element84 Earth Search from the open Sentinel-2 L2A archive on AWS. Geocoding: data (c) OpenStreetMap contributors (ODbL) via Nominatim; follow its [usage policy](https://operations.osmfoundation.org/policies/nominatim/) (one request per run here, identifying User-Agent).
