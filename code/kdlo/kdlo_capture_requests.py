#!/usr/bin/env python3
"""Resolve the CDSE burst metadata for the recorded KDLO measurements.

`onboarder.insar_measurements.burst_id` stores the **CDSE OData burst Id** (not a
tower-internal UUID): ``d716e5d7-...`` resolves to
``S1A-SLC-20230504T124402-024434-IW1-VH-381391``. The Bursts endpoint is public
(no token required), so the exact burst of every measurement — including its
polarisation, which varies by date — can be pinned without guessing.

Input  : ``--ids`` file with ``date<TAB>burst_id<TAB>orbit_direction`` lines,
         dumped from the live service with ``psql``:

             docker exec tower-postgres psql -U tower -d tower -A -F'\t' -t -c \
               "SELECT m.acquisition_ts::date, m.burst_id, m.orbit_direction
                FROM onboarder.insar_measurements m
                WHERE m.asset_id='93b18457-2e4d-4d20-a444-be1c4733423b'
                ORDER BY m.acquisition_ts" > ids.tsv

Output : ``data/kdlo/strips/capture_requests.json`` — the manifest the capture
         binary consumes (one entry per measurement, in acquisition order).

Run:  python3 code/kdlo/kdlo_capture_requests.py --ids /tmp/kdlo_ids.tsv
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time
import urllib.parse
import urllib.request

BURSTS_URL = "https://catalogue.dataspace.copernicus.eu/odata/v1/Bursts"
SITE_LAT = 44.96580177546652
SITE_LON = -97.58975377350735
STRIP_WIDTH = 11
STRIP_LINES = 400
DEFAULT_OUT = pathlib.Path(__file__).resolve().parents[2] / "data/kdlo/strips"


def fetch_burst(burst_id, tries=5, timeout=45):
    """One ``Bursts`` record by Id (the endpoint is public)."""
    query = BURSTS_URL + "?$filter=" + urllib.parse.quote(f"Id eq '{burst_id}'")
    last = None
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(query + "&$top=1", timeout=timeout) as fh:
                payload = json.load(fh)
            values = payload.get("value") or []
            return values[0] if values else None
        except Exception as exc:                       # transient 429/5xx/network
            last = exc
            time.sleep(3.0 * (attempt + 1))
    raise SystemExit(f"burst {burst_id}: {last}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ids", required=True,
                    help="TSV: date<TAB>burst_id<TAB>orbit_direction")
    ap.add_argument("--out", default=str(DEFAULT_OUT / "capture_requests.json"))
    args = ap.parse_args()

    entries = []
    lines = [l for l in pathlib.Path(args.ids).read_text().splitlines() if l.strip()]
    for i, line in enumerate(lines):
        # `psql -F'\t'` writes the two characters `\` `t` unless the shell
        # unescapes them, so accept both spellings of the separator.
        fields = line.replace("\\t", "\t").split("\t")
        date, burst_id, orbit_dir = (p.strip() for p in fields[:3])
        burst = fetch_burst(burst_id)
        if burst is None:
            raise SystemExit(f"burst {burst_id} ({date}) not found in the catalogue")
        entries.append({
            "date": date,
            "burst_id": burst_id,
            "orbit_direction": orbit_dir or burst.get("OrbitDirection"),
            "relative_orbit": burst.get("RelativeOrbitNumber"),
            "subswath": burst.get("SwathIdentifier"),
            "polarisation": burst.get("PolarisationChannels"),
            "burst_name": burst.get("Name"),
            "parent_product_id": burst.get("ParentProductId"),
            "parent_product_name": burst.get("ParentProductName"),
            "content_start": (burst.get("ContentDate") or {}).get("Start"),
            "lines_per_burst": burst.get("LinesPerBurst"),
            "lines": burst.get("Lines"),
            "samples_per_burst": burst.get("SamplesPerBurst"),
            "byte_offset": burst.get("ByteOffset"),
        })
        print(f"[{i + 1}/{len(lines)}] {date} {burst.get('Name')}")

    manifest = {
        "site": {"latitude": SITE_LAT, "longitude": SITE_LON,
                 "name": "KDLO-TV Tower Garden City SD (OSM 357089613)"},
        "strip": {"width": STRIP_WIDTH, "lines": STRIP_LINES,
                  "row_major": True,
                  "format": "u32 width, u32 height, 2*width*height little-endian f64 (re, im)",
                  "decoder": "decode_slc_rect_at_location(annotation, tiff, lat, lon, 11, 400)"},
        "source": {
            "endpoint": BURSTS_URL,
            "auth": "public (no token)",
            "measurement_table": "onboarder.insar_measurements",
            "asset_id": "93b18457-2e4d-4d20-a444-be1c4733423b",
        },
        "n_entries": len(entries),
        "entries": entries,
    }
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(manifest, indent=1) + "\n")
    pols = sorted({e["polarisation"] for e in entries})
    print(f"wrote {len(entries)} entries -> {out} (polarisations: {pols})")


if __name__ == "__main__":
    main()
