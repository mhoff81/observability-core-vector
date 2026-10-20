#!/usr/bin/env python3
"""HALTERN — the burst id behind every committed window (``data/haltern/reference/haltern_bursts.json``).

Why this file exists
--------------------
The committed window index (``data/haltern/haltern_windows_index.txt``) is the export of
``onboarder.insar_window_samples`` and has **eight** pipe-delimited fields:

    measurement_id|asset_id|asset_name|segment_index|acquisition_ts|burst_id|w|h

The **6th field is the CDSE burst id** — the very column the database stores under
that name. This package labelled it ``request_id`` (``haltern_ocv_masks.py``,
``read_index``), which collides with a *different* column of the measurement
extract: the pipeline's own request UUID (``haltern_measurements_full.txt``, field 5 —
**one single** distinct value in this record, and no value of it occurs in the
index). The label is corrected in this package and the fact behind it is pinned
here.

What the ids identify (the point of pinning them)
-------------------------------------------------
Every id resolves against the **public** CDSE OData catalogue to the SLC burst the
window was cut from, which gives the record a pass geometry it did not have
before: sub-swath, polarisation, relative orbit, azimuth time, the parent product
and the S3 object the pixels live in. That is what turns the ASC/DESC cross-cut of
§A2 from "two pass geometries" into four: **both** directions are two sub-swaths
of the IW mode on two relative orbits (ascending IW1/rel. orbit 15 and IW3/rel.
orbit 88, descending IW1/rel. orbit 37 and IW3/rel. orbit 139) at two times of
day — and it is the reason that block may only be read as a
*description*. ``fig_haltern_ocv_orbits.py`` reads this file and cross-checks
every id's ``OrbitDirection`` against the measurement extract's own
``orbit_direction`` column.

Two steps, like ``haltern_ocv_masks.py`` separates ``--extract``
------------------------------------------------------------
    --write   query the CDSE catalogue and (re)write the committed file [network]
    --check   re-verify the committed file offline against the committed index
              and the measurement extract (the default)

Resolving needs the live catalogue — a mutable external input — so the committed
file is pinned **structurally** (the identity facts below), never re-derived
offline. It carries no wall-clock timestamp: every field in it is the catalogue's
own datum, so two resolutions of an unchanged index are byte-identical.

    python3 code/haltern/haltern_bursts_resolve.py            # check (offline)
    python3 code/haltern/haltern_bursts_resolve.py --write    # resolve + write (network)
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import haltern_ocv_masks as masks  # noqa: E402

DATA = os.path.join(HERE, os.pardir, os.pardir, "data", "haltern")
OUT = os.path.join(DATA, "reference", "haltern_bursts.json")

ODATA = "https://catalogue.dataspace.copernicus.eu/odata/v1/Bursts"

# The catalogue fields this package uses, and nothing else. ``GeoFootprint`` is
# kept: it is the burst's ground geometry, i.e. the evidence that the window's
# burst covers the site, and ``S3Path`` is the object a re-decode would read.
KEEP = ["Id", "Name", "SwathIdentifier", "PolarisationChannels", "OrbitDirection",
        "RelativeOrbitNumber", "OperationalMode", "PlatformSerialIdentifier",
        "AzimuthTime", "BeginningDateTime", "EndingDateTime", "LinesPerBurst",
        "SamplesPerBurst", "AbsoluteBurstId", "ByteOffset", "S3Path",
        "ParentProductId", "ParentProductName", "GeoFootprint"]


# Structural identity of the committed file: these facts come from the live
# catalogue and cannot be re-derived from the record, so they are asserted
# instead — the same shape as the anchor's ``OSM_EXPECTED``.
EXPECTED = {
    "n_bursts": 120,
    "operational_mode": "IW",
    "swaths": ["IW1", "IW3"],
    "relative_orbits": [15, 37, 88, 139],
    "orbit_directions": ["ASCENDING", "DESCENDING"],
    "polarisations": ["VH", "VV"],
    # The headline of the pass geometry: **both** pass directions of this record
    # are *two* sub-swaths of the IW mode on *two* relative orbits — ascending is
    # IW1 on relative orbit 15 (~17:16 UT) and IW3 on relative orbit 88
    # (~17:24 UT), descending is IW1 on relative orbit 37 (~05:49 UT) and IW3 on
    # relative orbit 139 (~05:41 UT) — so the export's daily strongest-window
    # pick may alternate sub-swaths within one direction.
    "swath_per_orbit": {"ASCENDING": ["IW1", "IW3"], "DESCENDING": ["IW1", "IW3"]},
    "relative_orbit_per_orbit": {"ASCENDING": [15, 88], "DESCENDING": [37, 139]},
}


def index_bursts():
    """``burst_id -> [measurement_id, ...]`` of the committed index."""
    out = {}
    for mid, meta in masks.read_index().items():
        out.setdefault(meta["burst_id"], []).append(mid)
    return out


def resolve(burst_ids):
    """Query the catalogue for every id and keep the fields in :data:`KEEP`."""
    import requests

    recs = {}
    for i, bid in enumerate(sorted(burst_ids), 1):
        params = {"$filter": f"Id eq '{bid}'", "$top": "1"}
        last, value = None, []
        for attempt in range(3):
            try:
                r = requests.get(ODATA, params=params, timeout=60)
                r.raise_for_status()
                value = r.json().get("value") or []
                break
            except Exception as exc:                      # noqa: BLE001
                last = exc
                time.sleep(1.5 * (attempt + 1))
        if not value:
            raise SystemExit(f"{bid}: not resolvable against {ODATA} ({last})")
        rec = value[0]
        missing = [k for k in KEEP if k not in rec]
        if missing:
            raise SystemExit(f"{bid}: catalogue record lacks {missing}")
        recs[bid] = {k: rec[k] for k in KEEP}
        print(f"  [{i}/{len(burst_ids)}] {bid} -> {rec['Name']} "
              f"({rec['SwathIdentifier']} {rec['PolarisationChannels']} "
              f"{rec['OrbitDirection']}, rel. orbit {rec['RelativeOrbitNumber']}, "
              f"{rec['AzimuthTime']})")
        time.sleep(0.3)
    return recs


def write():
    """Query the catalogue, rewrite the committed file and check it."""
    ids = index_bursts()
    print(f"resolving {len(ids)} burst ids against {ODATA} …")
    recs = resolve(ids)
    payload = {
        "source": ODATA,
        "query": "$filter=Id eq '<burst_id>' (one query per distinct index id)",
        "index": "data/haltern/haltern_windows_index.txt",
        "index_field": "6th field (burst_id) of data/haltern/haltern_windows_index.txt",
        "index_sha256": masks._sha256(masks.INDEX_PATH),
        "collected_by": "code/haltern/haltern_bursts_resolve.py --write",
        "note": ("resolved from the live CDSE catalogue — a mutable external "
                 "input, so this file is pinned structurally by --check and by "
                 "the pin block of figs. A-E, never re-derived offline"),
        "n_bursts": len(recs),
        "bursts": recs,
    }
    os.makedirs(os.path.dirname(os.path.abspath(OUT)), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(payload, fh, indent=1, sort_keys=True)
        fh.write("\n")
    print(f"written: {os.path.relpath(OUT)} ({len(recs)} bursts)")
    return payload


def check(payload=None):
    """Re-verify the committed file offline. Returns ``(problems, facts)``."""
    if payload is None:
        with open(OUT) as fh:
            payload = json.load(fh)
    problems, recs = [], payload["bursts"]
    ids = index_bursts()
    meas = masks.read_measurements()

    # 1) this file and the committed index describe the same windows ---------
    if payload.get("index_sha256") != masks._sha256(masks.INDEX_PATH):
        problems.append("index_sha256: the committed index changed since the "
                        "resolution was written — re-run --write")
    if set(recs) != set(ids):
        problems.append("burst id set differs from the index: "
                        f"{sorted(set(recs) ^ set(ids))}")
    if len(recs) != EXPECTED["n_bursts"]:
        problems.append(f"n_bursts {len(recs)} != {EXPECTED['n_bursts']}")

    # 2) the identity facts of the pass geometry ----------------------------
    got = {key: sorted({r[field] for r in recs.values()}) for key, field in (
        ("swaths", "SwathIdentifier"), ("relative_orbits", "RelativeOrbitNumber"),
        ("orbit_directions", "OrbitDirection"),
        ("polarisations", "PolarisationChannels"))}
    for key, value in got.items():
        if value != EXPECTED[key]:
            problems.append(f"{key} {value} != {EXPECTED[key]}")
    if sorted({r["OperationalMode"] for r in recs.values()}) != \
            [EXPECTED["operational_mode"]]:
        problems.append("operational modes != "
                        f"[{EXPECTED['operational_mode']}]")
    # > the pass geometry of this record is a sub-swath, not only a direction
    swath_per_orbit = {o: sorted({r["SwathIdentifier"] for r in recs.values()
                                  if r["OrbitDirection"] == o})
                       for o in got["orbit_directions"]}
    relorb_per_orbit = {o: sorted({r["RelativeOrbitNumber"] for r in recs.values()
                                   if r["OrbitDirection"] == o})
                        for o in got["orbit_directions"]}
    if swath_per_orbit != EXPECTED["swath_per_orbit"]:
        problems.append(f"swath per orbit {swath_per_orbit} != "
                        f"{EXPECTED['swath_per_orbit']}")
    if relorb_per_orbit != EXPECTED["relative_orbit_per_orbit"]:
        problems.append(f"relative orbit per orbit {relorb_per_orbit} != "
                        f"{EXPECTED['relative_orbit_per_orbit']}")

    # 3) every burst's orbit direction against the measurement extract -------
    n_agrees = n_with_meas = 0
    for bid, mids in ids.items():
        rec = recs.get(bid)
        if rec is None:
            continue
        seen = {meas[m].get("orbit_direction") for m in mids if m in meas}
        seen.discard(None)
        if not seen:
            continue
        n_with_meas += 1
        if len(seen) > 1:
            problems.append(f"{bid}: its windows disagree on orbit_direction "
                            f"{sorted(seen)}")
        elif rec["OrbitDirection"] not in seen:
            problems.append(f"{bid}: catalogue {rec['OrbitDirection']} != "
                            f"measurement {sorted(seen)}")
        else:
            n_agrees += 1
    if n_agrees != n_with_meas:
        problems.append(f"orbit direction agrees on {n_agrees} of "
                        f"{n_with_meas} bursts")

    facts = {
        "n_bursts": len(recs),
        "swath_per_orbit": swath_per_orbit,
        "relative_orbit_per_orbit": relorb_per_orbit,
        "polarisations": got["polarisations"],
        "orbit_direction_agrees_with_measurements": {
            "n_bursts": n_with_meas, "n_agrees": n_agrees},
    }
    return problems, facts


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--write", action="store_true",
                    help="query the CDSE catalogue and rewrite the committed file")
    args = ap.parse_args(argv)
    problems, facts = check(write() if args.write else None)
    if problems:
        print("deviations:", file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)
        raise SystemExit(1)
    print(f"ok: {os.path.relpath(OUT)} — {facts['n_bursts']} bursts, swath per "
          f"orbit {facts['swath_per_orbit']}, relative orbit per orbit "
          f"{facts['relative_orbit_per_orbit']}, polarisations "
          f"{facts['polarisations']}; the catalogue's orbit direction agrees "
          f"with the measurement extract on "
          f"{facts['orbit_direction_agrees_with_measurements']['n_agrees']} of "
          f"{facts['n_bursts']} bursts")


if __name__ == "__main__":
    main()
