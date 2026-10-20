#!/usr/bin/env python3
"""WERDER — the bursts under every committed acquisition date (``data/werder_havel/reference/werder_havel_bursts.json``).

The converse of Haltern
-----------------------
Haltern's export committed a burst id per window, so ``haltern_bursts_resolve.py``
asks the catalogue about 120 *ids*. This record committed the opposite: 991 rows,
28 columns, **no burst id and no coordinate** (its ``request_id`` is the track tag —
``werder_havel_ocv_masks.py``). What it does commit is a bridge outline — 704.0426 m
x 35.9292 m, OSM way 413381418 (``werder_havel_aoi.py``) — and 491 acquisition dates,
and its own ``track`` column assigns every one of those dates to exactly one of four
pass geometries:

    asc_146_iw1    126 dates    ASCENDING   rel. orbit 146   IW1
    asc_44_iw3     118          ASCENDING   rel. orbit  44   IW3
    desc_168_iw1   124          DESCENDING  rel. orbit 168   IW1
    desc_95_iw3    123          DESCENDING  rel. orbit  95   IW3

(the ``ref_`` and ``refout_`` tags re-measure the same four geometries — eight tags,
four geometries, 126+118+124+123 = 491 dates). So this file asks **the outline
itself**, with one OData query per acquisition date:

    /odata/v1/Bursts?$filter=OData.CSC.Intersects(area=geography'SRID=4326;<ring>')
        and ContentDate/Start ge <date>T00:00:00.000Z
        and ContentDate/Start lt <date>T23:59:59.999Z

— 491 date queries, which is why the answer can be, and is, committed whole.

What the answer is (and why it is checkable)
--------------------------------------------
Each date answers with the bursts of *its own* pass geometry, in the IW SLC dual-pol
pair (VV and VH: two catalogue records per frame — distinct ``Id``, one shared
``BurstId``), and where the 704 m outline straddles a burst frame boundary it answers
with *two frames*: on the ascending IW3 geometry the outline sits on the boundary of
the 092669/092670 pair (2.75 s apart in azimuth) on **every one** of its 118 dates,
and the other three geometries answer with one frame on every one of theirs. The
committed layer is therefore 118 x 2 + 126 + 124 + 123 = **609 frames in 1218
catalogue records** over the 491 dates, and ``n_frames`` per date is a pinned fact per
geometry (2, 1, 1, 1).

None of this could be read off the record — which is exactly why the record's own
columns are the test of it, and ``--check`` runs that test date by date: every burst
of a date must carry the ``orbit_direction``, ``relative_orbit``, ``subswath`` and
platform the extract states for it, its ``AzimuthTime`` must sit within 5 s of the
extract's ``acquisition_ts`` (the layer's worst offset is 2.0 s), and its footprint
must contain the committed site point — which all 1218 of them do, each of them also
containing all 10 vertices of the committed outline (the geometric meaning of the
``/Intersects()`` that asked for it).

What this file does *not* claim
-------------------------------
That a burst here is the one the pipeline measured. The export committed no burst id,
so what is pinned is the **covering** burst layer of the record — its pass geometry,
its polarisations, its frames and its object paths — not a reconstruction of the
pipeline's own choice inside that layer. The pixel audit of Phase 3
(``werder_havel_geometry.json``) is where that distinction is priced.

Payload of a burst (the number the download step prices)
--------------------------------------------------------
The catalogue commits no file size for a burst, but a burst SLC is CS16 (complex
int16, 4 bytes per sample), so

    payload_bytes = 4 x LinesPerBurst x SamplesPerBurst

is what one burst of this layer weighs: **128.9 .. 158.0 MB** over the committed
layer — ~129..140 MB on the two IW1 geometries (1502 x 22 961 on the ascending one,
1507 x 22 774 on the descending one) and ~155..158 MB on the two IW3 geometries
(1517 x 25 632 / 1524 x 25 521), with the second platform (S1C, 42 records — 21 of
the dates — per track) trimming a few lines and widening the IW3 samples slightly. It
is a **definition**, like ``WINDOW_SIZE`` in ``werder_havel_ocv_masks.py``, not a
measurement: ``--write`` prints the range over the whole record and ``--check``
re-derives it from the committed file.

Two steps, like the rest of the family
--------------------------------------
    --write     query the CDSE catalogue and (re)write the committed file [network]
    --check     re-verify the committed file offline (the default)
    --probe N   query the N first dates and print, touch nothing [network]

Resolving needs the live catalogue — a mutable external input — so the committed file
is pinned **structurally** (the identity facts below), never re-derived offline, and it
carries no wall-clock timestamp: every field in it is either the catalogue's own datum
or the record's own ``track`` tag for that date.

    python3 code/werder_havel/werder_havel_bursts_resolve.py           # check (offline)
    python3 code/werder_havel/werder_havel_bursts_resolve.py --write   # resolve + write
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import werder_havel_aoi as aoi  # noqa: E402
import werder_havel_ocv_masks as masks  # noqa: E402

DATA = os.path.join(HERE, os.pardir, os.pardir, "data", "werder_havel")
OUT = os.path.join(DATA, "reference", "werder_havel_bursts.json")

ODATA = "https://catalogue.dataspace.copernicus.eu/odata/v1/Bursts"
COLLECTED_BY = "code/werder_havel/werder_havel_bursts_resolve.py --write"
QUERY = ("$filter=OData.CSC.Intersects(area=geography'SRID=4326;<the committed AOI "
         "ring>) and ContentDate/Start ge <date>T00:00:00.000Z and ContentDate/Start "
         "lt <date>T23:59:59.999Z, $top=100 (one query per acquisition date)")

NOTE = ("derived from the live CDSE catalogue as the Intersects() of the committed "
        "OSM outline (the record carries no burst id): a mutable external input, so "
        "this file is pinned structurally by --check and by the pin block of figs. "
        "A-E, never re-derived offline")

# The catalogue fields this package commits, and nothing else. ``GeoFootprint`` is
# the burst's ground geometry — the evidence that it covers the committed outline —
# ``ContentDate`` is the field the per-date query filtered on, and ``S3Path`` is the
# object a re-decode would read.
KEEP = ["Id", "Name", "BurstId", "AbsoluteBurstId", "SwathIdentifier",
        "PolarisationChannels", "OrbitDirection", "RelativeOrbitNumber",
        "OperationalMode", "PlatformSerialIdentifier", "AzimuthTime", "ContentDate",
        "BeginningDateTime", "EndingDateTime", "LinesPerBurst", "SamplesPerBurst",
        "DatatakeID", "ByteOffset", "S3Path", "ParentProductId", "ParentProductName",
        "ParentProductType", "GeoFootprint"]

# CS16: a burst SLC sample is two int16, i.e. 4 bytes (a definition, not a read).
BYTES_PER_SAMPLE = 4

# Structural identity of the committed file. The geometry facts come from the live
# catalogue and cannot be re-derived from the record, so they are asserted instead —
# the same shape as the anchor's ``OSM_EXPECTED`` and Haltern's ``EXPECTED``.
EXPECTED = {
    "n_dates": 491,
    "n_queries": 491,
    "operational_mode": "IW",
    "swaths": ["IW1", "IW3"],
    "relative_orbits": [44, 95, 146, 168],
    "orbit_directions": ["ASCENDING", "DESCENDING"],
    "polarisations": ["VH", "VV"],
    "platforms": ["A", "C"],
    # the four pass geometries of the record, as the extract's own ``track`` tags
    # name them — the key is ``<orbit_direction>/<relative_orbit>/<subswath>``
    "track_of_geometry": {
        "ASCENDING/146/IW1": "asc_146_iw1",
        "ASCENDING/44/IW3": "asc_44_iw3",
        "DESCENDING/168/IW1": "desc_168_iw1",
        "DESCENDING/95/IW3": "desc_95_iw3",
    },
    "n_dates_per_track": {"asc_146_iw1": 126, "asc_44_iw3": 118,
                          "desc_168_iw1": 124, "desc_95_iw3": 123},
    # the headline of the resolution: the outline straddles a burst frame boundary on
    # ascending IW3 only, so that geometry carries two frames on *every* one of its
    # dates and the other three carry one — 118 x 2 + 126 + 124 + 123 = 609 frames
    # over 491 dates (values are the distinct counts seen over the track's dates)
    "frames_per_track_per_date": {"asc_146_iw1": [1], "asc_44_iw3": [2],
                                  "desc_168_iw1": [1], "desc_95_iw3": [1]},
    "n_frames": 609,
    "n_burst_records": 1218,      # the IW SLC dual-pol pair (VV, VH) per frame
    # the test the record's own columns make of the catalogue's answer
    "azimuth_time_max_offset_s": 5.0,
    "site_point_in_every_footprint": True,
    # the ``$top`` the probe fixed: no date of this record has more than four
    # catalogue records (2 frames x 2 polarisations), and ``$count`` is what says so
    "top": 100,
}


# ---------------------------------------------------------------------------
# The record's own tags, and the geometry they name
# ---------------------------------------------------------------------------
def geometry_key(orbit_direction, relative_orbit, subswath):
    """``'ASCENDING/146/IW1'`` — the key both the record and the catalogue are read on."""
    return f"{orbit_direction}/{int(relative_orbit)}/{subswath}"


def record_dates():
    """``{date: row}`` of the record, one row per acquisition date.

    The four target tracks partition the 491 dates (the ``ref_``/``refout_`` rows
    repeat them), so a date maps to exactly one row once the reference rows are
    dropped — and the ambiguity is *raised*, not resolved silently, if that ever
    stops being true.
    """
    rows = {}
    for row in masks.read_measurements().values():
        if row["role"] != "target":
            continue
        date = row["date"]
        if date in rows and rows[date]["track"] != row["track"]:
            raise SystemExit(f"{date}: two target tracks in the record: "
                             f"{rows[date]['track']} and {row['track']}")
        rows[date] = row
    return rows


def track_geometries(rows):
    """``{track: geometry_key}`` of the extract's eight tags (four geometries)."""
    out = {}
    for row in masks.read_measurements().values():
        key = geometry_key(row["orbit_direction"], row["relative_orbit"],
                           row["subswath"])
        tag = row["track"]
        if out.setdefault(tag, key) != key:
            raise SystemExit(f"{tag}: the record gives it two geometries")
    return out


def aoi_ring():
    """The committed outline as ``[[lon, lat], ...]`` (the query's polygon)."""
    with open(aoi.OUT) as fh:
        return json.load(fh)["aoi"]["polygon_lonlat"]


def aoi_polygon(ring):
    """``POLYGON((lon lat, ...))`` of a GeoJSON-order ring."""
    return "POLYGON((" + ", ".join(f"{lon} {lat}" for lon, lat in ring) + "))"


def date_filter(date):
    """The ``$filter`` of one acquisition date of this record."""
    ring = aoi_polygon(aoi_ring())
    return (f"OData.CSC.Intersects(area=geography'SRID=4326;{ring}') "
            f"and ContentDate/Start ge {date}T00:00:00.000Z "
            f"and ContentDate/Start lt {date}T23:59:59.999Z")


# ---------------------------------------------------------------------------
# The two offline readings of a burst record
# ---------------------------------------------------------------------------
def _parse_iso(stamp):
    """``datetime`` of a catalogue ISO stamp (second resolution is enough here)."""
    return datetime.strptime(stamp[:19], "%Y-%m-%dT%H:%M:%S").replace(
        tzinfo=timezone.utc)


def payload_bytes(rec):
    """``4 x LinesPerBurst x SamplesPerBurst`` — the CS16 weight of one burst."""
    return BYTES_PER_SAMPLE * int(rec["LinesPerBurst"]) * int(rec["SamplesPerBurst"])


def footprint_ring(rec):
    """The burst footprint's outer ring as ``[(lon, lat), ...]``."""
    return [(p[0], p[1]) for p in rec["GeoFootprint"]["coordinates"][0]]


def point_in_ring(lat, lon, ring):
    """Ray casting in lon/lat — enough for a point against a burst quadrilateral."""
    inside = False
    n = len(ring)
    for i in range(n):
        x1, y1 = ring[i]
        x2, y2 = ring[(i + 1) % n]
        if (y1 > lat) != (y2 > lat):
            x_cross = x1 + (lat - y1) * (x2 - x1) / (y2 - y1)
            if lon < x_cross:
                inside = not inside
    return inside


# ---------------------------------------------------------------------------
# Resolve: one query per acquisition date [network]
# ---------------------------------------------------------------------------
def resolve(dates, top=None, quiet=False):
    """Query the catalogue for every date and keep the fields in :data:`KEEP`.

    One query per date, retried three times like the family's other resolvers, and
    a **hard failure** if a date answers with nothing: this record's dates are its
    own acquisition days, so an empty answer is a finding about the catalogue (or
    about the query) and never a row to skip silently.
    """
    import requests

    top = EXPECTED["top"] if top is None else top
    out = {}
    for i, date in enumerate(sorted(dates), 1):
        params = {"$filter": date_filter(date), "$top": str(top), "$count": "true"}
        last, body = None, {}
        for attempt in range(3):
            try:
                r = requests.get(ODATA, params=params, timeout=120)
                r.raise_for_status()
                body = r.json()
                break
            except Exception as exc:                      # noqa: BLE001
                last = exc
                time.sleep(1.5 * (attempt + 1))
        value = body.get("value") or []
        if not value:
            raise SystemExit(f"{date}: no burst on {ODATA} ({last}) — filter: "
                             f"{params['$filter']}")
        if body.get("@odata.count") not in (None, len(value)):
            raise SystemExit(f"{date}: $count {body['@odata.count']} != {len(value)} "
                             f"returned at $top={top}")
        recs = {}
        for rec in value:
            missing = [k for k in KEEP if k not in rec]
            if missing:
                raise SystemExit(f"{date}: catalogue record lacks {missing}")
            recs[rec["Id"]] = {k: rec[k] for k in KEEP}
        out[date] = recs
        if not quiet:
            frames = sorted({r["BurstId"] for r in recs.values()})
            sizes = sorted({round(payload_bytes(r) / 1e6, 1) for r in recs.values()})
            print(f"  [{i}/{len(dates)}] {date}: {len(recs)} records, "
                  f"{len(frames)} frame(s), {sizes} MB — "
                  f"{recs[sorted(recs)[0]]['Name']}")
    return out


def write():
    """Query the catalogue, rewrite the committed file and check it."""
    rows = record_dates()
    tracks = track_geometries(rows)
    print(f"resolving {len(rows)} acquisition dates against {ODATA} …")
    per_date = resolve(rows)
    dates = {}
    for date, recs in per_date.items():
        frames = {r["BurstId"] for r in recs.values()}
        dates[date] = {"track": rows[date]["track"], "n_frames": len(frames),
                       "bursts": recs}
    payload = {
        "source": ODATA,
        "query": QUERY,
        "aoi": "data/werder_havel/reference/werder_havel_aoi.json",
        "aoi_sha256": masks._sha256(aoi.OUT),
        "aoi_polygon_role": "the /Intersects() polygon of every query above",
        "record": "data/werder_havel/werder_havel_measurements_full.txt",
        "record_sha256": masks._sha256(masks.MEAS_PATH),
        "collected_by": COLLECTED_BY,
        "note": NOTE,
        "bytes_per_sample": BYTES_PER_SAMPLE,
        "n_dates": len(dates),
        "n_queries": len(dates),
        "n_frames": sum(d["n_frames"] for d in dates.values()),
        "n_burst_records": sum(len(d["bursts"]) for d in dates.values()),
        "track_geometries": tracks,
        "dates": dates,
    }
    os.makedirs(os.path.dirname(os.path.abspath(OUT)), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(payload, fh, indent=1, sort_keys=True)
        fh.write("\n")
    print(f"written: {os.path.relpath(OUT)} ({len(dates)} dates, "
          f"{payload['n_frames']} frames, {payload['n_burst_records']} records)")
    return payload


# ---------------------------------------------------------------------------
# Check (offline)
# ---------------------------------------------------------------------------
def check(payload=None):
    """Re-verify the committed file offline. Returns ``(problems, facts)``."""
    problems = []
    rows = record_dates()
    tracks = track_geometries(rows)
    if payload is None:
        if not os.path.exists(OUT):
            return [f"{os.path.relpath(OUT)}: missing — run --write"], {}
        with open(OUT) as fh:
            payload = json.load(fh)
    dates = payload.get("dates") or {}

    # 1) this file and the two committed inputs it was derived from ------------
    if payload.get("record_sha256") != masks._sha256(masks.MEAS_PATH):
        problems.append("record_sha256: the committed extract changed since the "
                        "resolution was written — re-run --write")
    if payload.get("aoi_sha256") != masks._sha256(aoi.OUT):
        problems.append("aoi_sha256: the committed outline changed since the "
                        "resolution was written — re-run --write")
    if payload.get("collected_by") != COLLECTED_BY:
        problems.append(f"collected_by {payload.get('collected_by')!r} != "
                        f"{COLLECTED_BY!r}")
    # the file's own statement of what was asked and how it may be read is part of
    # the pin: editing the filter text or the mutable-input caveat must not pass
    if payload.get("query") != QUERY:
        problems.append("query: the committed filter text is not this module's — "
                        "the file and the code that writes it disagree")
    if payload.get("note") != NOTE:
        problems.append("note: the committed provenance caveat is not this "
                        "module's")
    if payload.get("aoi") != "data/werder_havel/reference/werder_havel_aoi.json" \
            or payload.get("aoi_polygon_role") != \
            "the /Intersects() polygon of every query above":
        problems.append("aoi/aoi_polygon_role: the committed pointer to the query "
                        "polygon is not this module's")
    if payload.get("source") != ODATA:
        problems.append(f"source {payload.get('source')!r} != {ODATA!r}")
    if payload.get("bytes_per_sample") != BYTES_PER_SAMPLE:
        problems.append(f"bytes_per_sample {payload.get('bytes_per_sample')!r} != "
                        f"{BYTES_PER_SAMPLE}")
    if set(dates) != set(rows):
        problems.append("the committed dates differ from the record's: "
                        f"{sorted(set(dates) ^ set(rows))[:8]}")
    if payload.get("n_dates") != len(dates) or \
            payload.get("n_dates") != EXPECTED["n_dates"]:
        problems.append(f"n_dates {payload.get('n_dates')} != {len(dates)} / "
                        f"{EXPECTED['n_dates']}")
    if payload.get("n_queries") != EXPECTED["n_queries"]:
        problems.append(f"n_queries {payload.get('n_queries')} != "
                        f"{EXPECTED['n_queries']}")
    if payload.get("track_geometries") != tracks:
        problems.append("track_geometries: the committed tag->geometry map is not "
                        "what the extract now says")

    # 2) the record's own columns against the catalogue's answer, date by date --
    with open(aoi.OUT) as fh:
        site = json.load(fh)["site_point"]
    ring = aoi_ring()
    n_records = n_frames = 0
    per_track_frames = {}
    payload_mb = []
    az_max = 0.0
    n_outside = 0
    min_ring_in = len(ring)

    for date in sorted(dates):
        entry = dates[date]
        recs = entry["bursts"]
        row = rows.get(date)
        if row is None:
            continue
        track = entry.get("track")
        if track != row["track"]:
            problems.append(f"{date}: track {track!r} != the extract's "
                            f"{row['track']!r}")
            continue
        want_geometry = tracks[track]
        want_pol = set()
        for rid, rec in recs.items():
            if rid != rec["Id"]:
                problems.append(f"{date}: key {rid} != Id {rec['Id']}")
            key = geometry_key(rec["OrbitDirection"], rec["RelativeOrbitNumber"],
                               rec["SwathIdentifier"])
            if key != want_geometry:
                problems.append(f"{date}: catalogue {key} != extract {want_geometry}")
            if rec["OperationalMode"] != EXPECTED["operational_mode"]:
                problems.append(f"{date}: operational mode "
                                f"{rec['OperationalMode']}")
            if "S1" + rec["PlatformSerialIdentifier"] != row["platform"]:
                problems.append(f"{date}: platform "
                                f"S1{rec['PlatformSerialIdentifier']} != extract "
                                f"{row['platform']}")
            if rec["ContentDate"]["Start"][:10] != date:
                problems.append(f"{date}: ContentDate/Start "
                                f"{rec['ContentDate']['Start']} is not that day")
            want_pol.add(rec["PolarisationChannels"])
            offset = abs((_parse_iso(rec["AzimuthTime"])
                          - _parse_iso(row["acquisition_ts"])).total_seconds())
            az_max = max(az_max, offset)
            if offset > EXPECTED["azimuth_time_max_offset_s"]:
                problems.append(f"{date}: AzimuthTime {rec['AzimuthTime']} is "
                                f"{offset:.1f} s from acquisition_ts "
                                f"{row['acquisition_ts']}")
            fp = footprint_ring(rec)
            if fp[0] != fp[-1]:
                problems.append(f"{date}: the burst footprint is not a closed ring")
            if not point_in_ring(site["lat"], site["lon"], fp):
                n_outside += 1
                problems.append(f"{date}: the site point is outside the footprint "
                                f"of {rec['Name']}")
            min_ring_in = min(min_ring_in, sum(1 for lon, lat in ring
                                               if point_in_ring(lat, lon, fp)))
            payload_mb.append(payload_bytes(rec) / 1e6)
        if want_pol != set(EXPECTED["polarisations"]):
            problems.append(f"{date}: polarisations {sorted(want_pol)} != "
                            f"{EXPECTED['polarisations']}")
        frames = len({r["BurstId"] for r in recs.values()})
        if entry.get("n_frames") != frames:
            problems.append(f"{date}: n_frames {entry.get('n_frames')} != "
                            f"{frames} distinct BurstId in bursts")
        n_frames += frames
        n_records += len(recs)
        per_track_frames.setdefault(track, set()).add(frames)



    # 3) the identity facts of the layer, per track and per geometry ------------
    per_track_frames = {t: sorted(v) for t, v in per_track_frames.items()}
    if per_track_frames != EXPECTED["frames_per_track_per_date"]:
        problems.append(f"frames per date per track {per_track_frames} != "
                        f"{EXPECTED['frames_per_track_per_date']}")
    n_dates_per_track = {t: sum(1 for d in dates.values() if d.get("track") == t)
                         for t in sorted(EXPECTED["n_dates_per_track"])}
    if n_dates_per_track != EXPECTED["n_dates_per_track"]:
        problems.append(f"dates per track {n_dates_per_track} != "
                        f"{EXPECTED['n_dates_per_track']}")
    # the extract's own tag column, folded onto the four geometries: the committed
    # map must be exactly it (the reference tags of a geometry carry its target tag)
    got_map = {}
    for tag, key in sorted(tracks.items()):
        target = tag.replace("refout_", "").replace("ref_", "")
        if got_map.setdefault(key, target) != target:
            problems.append(f"{key}: the extract maps it to two target tags")
    if got_map != EXPECTED["track_of_geometry"]:
        problems.append(f"track of geometry {got_map} != "
                        f"{EXPECTED['track_of_geometry']}")
    allrecs = [r for d in dates.values() for r in d["bursts"].values()]
    got = {key: sorted({r[field] for r in allrecs}) for key, field in (
        ("swaths", "SwathIdentifier"), ("relative_orbits", "RelativeOrbitNumber"),
        ("orbit_directions", "OrbitDirection"),
        ("polarisations", "PolarisationChannels"),
        ("platforms", "PlatformSerialIdentifier"))}
    for key, value in got.items():
        if value != EXPECTED[key]:
            problems.append(f"{key} {value} != {EXPECTED[key]}")
    if payload.get("n_frames") != n_frames or n_frames != EXPECTED["n_frames"]:
        problems.append(f"n_frames {payload.get('n_frames')} / {n_frames} != "
                        f"{EXPECTED['n_frames']}")
    if payload.get("n_burst_records") != n_records or \
            n_records != EXPECTED["n_burst_records"]:
        problems.append(f"n_burst_records {payload.get('n_burst_records')} / "
                        f"{n_records} != {EXPECTED['n_burst_records']}")

    facts = {
        "n_dates": len(dates),
        "n_queries": payload.get("n_queries"),
        "n_frames": n_frames,
        "n_burst_records": n_records,
        "n_dates_per_track": n_dates_per_track,
        "frames_per_track_per_date": per_track_frames,
        "payload_mb_range": ([round(min(payload_mb), 1), round(max(payload_mb), 1)]
                             if payload_mb else None),
        "azimuth_time_max_offset_s": round(az_max, 3),
        "site_point_in_every_footprint": n_outside == 0,
        "min_ring_vertices_in_footprint": min_ring_in,
        "ring_vertices": len(ring),
        "swaths": got["swaths"],
        "platforms": got["platforms"],
        "relative_orbits": got["relative_orbits"],
    }
    return problems, facts


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--write", action="store_true",
                    help="query the CDSE catalogue and rewrite the committed file")
    ap.add_argument("--check", action="store_true",
                    help="re-verify the committed file offline (the default)")
    ap.add_argument("--probe", type=int, metavar="N",
                    help="query the N first dates and print, touch nothing")
    ap.add_argument("--json", action="store_true", help="print the facts and exit")
    args = ap.parse_args(argv)

    if args.probe:
        rows = record_dates()
        print(f"probing {args.probe} of {len(rows)} dates against {ODATA} …")
        resolve(sorted(rows)[:args.probe])
        return 0

    problems, facts = check(write() if args.write else None)
    if args.json:
        print(json.dumps(facts, indent=1, sort_keys=True))
    if problems:
        print(f"CHECK FAILED — {len(problems)} problem(s)")
        for p in problems[:40]:
            print(f"  - {p}")
        return 1
    if args.json:
        return 0
    print(f"OK — {os.path.relpath(OUT)}: {facts['n_dates']} dates in "
          f"{facts['n_queries']} catalogue queries -> {facts['n_frames']} burst "
          f"frames in {facts['n_burst_records']} records; frames per date "
          f"{facts['frames_per_track_per_date']}; dates per track "
          f"{facts['n_dates_per_track']}")
    print(f"     payload {facts['payload_mb_range'][0]} .. "
          f"{facts['payload_mb_range'][1]} MB per burst (CS16); |AzimuthTime - "
          f"acquisition_ts| <= {facts['azimuth_time_max_offset_s']} s; the site "
          f"point lies inside every footprint "
          f"{facts['site_point_in_every_footprint']} (at least "
          f"{facts['min_ring_vertices_in_footprint']} of the "
          f"{facts['ring_vertices']} outline vertices in each); platforms "
          f"{facts['platforms']}, relative orbits {facts['relative_orbits']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

