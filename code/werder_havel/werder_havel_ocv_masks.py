#!/usr/bin/env python3
"""WERDER OCV — the record readers and the echo-mask rule of the Havelbrücke Werder record.

Havelbrücke Werder (A 10, Berliner Ring, western ring) crosses the Potsdamer Havel /
**Großer Zernsee** between the Phöben and Leest junctions at 52.41753 N 12.91467 E
(Werder (Havel), Potsdam-Mittelmark, Brandenburg). It is a **704.043 m x 35.929 m**
closed OSM outline: an eleven-span, continuous steel girder bridge, deck on top
(``bridge:structure=beam``, ``material=steel``), whose documented span division is
2 x 55 m end spans + 9 x 66 m inner spans = 704.0 m — 0.04 m (0.006 %) shorter than
the map outline this package anchors to.

This site is the family's **positive control**. Its superstructures were renewed in
two steps — *Richtungsfahrbahn Nauen 1997*, *Richtungsfahrbahn Werder 2004*, both
**before** the Sentinel-1 archive — and the one documented state change *inside* the
archive is the renewal of the expansion joints, of more than 1800 m of guardrails, of
the mastic asphalt and of the corrosion protection of the steel superstructure,
begun 2024-07-01 and finished 2025-09-05 (Autobahn GmbH). Nothing here failed:
``event_kind`` is a documented change, never a precursor, and the value of the
package is the *positive* direction of the channel.

What the record carries
-----------------------
``data/werder_havel/werder_havel_measurements_full.txt`` is the export of
``onboarder.measurements`` for this asset: **991 rows**, 28 pipe-delimited columns,
one row per ``(track, date)`` — the strata are

  * ``track`` — eight Sentinel-1 tags: ``asc_146_iw1``, ``asc_44_iw3``,
    ``desc_168_iw1``, ``desc_95_iw3`` and the ``ref_`` / ``refout_`` pair of each of
    the two IW1 tracks (the reference and the outside-reference chains);
  * ``role`` — ``target`` (491), ``reference`` (250), ``reference_outside`` (250);
  * ``state_window`` — ``W2_bestand`` (241 dates, 2022-07-02 .. 2024-06-29),
    ``W3_sanierung`` (189, 2024-07-03 .. 2025-09-04), ``W4_post`` (61,
    2025-09-05 .. 2025-12-04), over **491 distinct acquisition dates**.

The record carries **no coordinates and no catalogue id**: none of the 28 columns is
a latitude, a longitude or a burst id, and ``request_id`` is the *track tag* (eight
distinct values, i.e. one per track — not the pipeline's request UUID). The anchor of
this package is therefore the bridge outline itself
(``code/werder_havel/werder_havel_aoi.py``, OSM way 413381418) and the burst layer is
*derived* from that area (``code/werder_havel/werder_havel_bursts_resolve.py``),
which is the opposite of Haltern, whose export committed a burst id per window.

The echo mask
-------------
The rule below is the project-wide tower/bridge echo definition — identical to
``haltern_ocv_masks.echo_mask`` and to the backend's
``sarsolve::coherence::coherence_mast_echo``:

    intensity = |z|^2
    peak      = max(intensity)
    median    = np.median(intensity)
    mask      = intensity >= 0.30 * peak  and  intensity >= 5.0 * median
    A         = count(mask)                    (no echo if A < 2)
    gamma2    = |sum(z[mask])|^2 / (A * sum(|z[mask]|^2)),  clamped to [0, 1]

The window of this record is **80 x 80 = 6400 complex samples** per acquisition
(``WINDOW_SIZE``) — 131x the 49-sample payload of Haltern and YWF, and ~250 MB of
complex JSON across the 991 rows, which is why that payload is the one **external**
input of this package: it is streamed on demand (``read_windows``, ``read_window_ids``)
and everything this package commits is derived from it. The export's own
``coherence_masked_pixels`` column is **640 in every one of the 991 rows**, i.e.
exactly 10 % of 6400: that column is a **definition** (the brightest decile of the
window), *not* a count of the mask the rule above produces. The two are kept apart
here and wherever the number is used; ``werder_havel_aoi_geometry.py`` is where the
constant is priced against the structure, and ``--json`` below prints the census of
the record's own columns.

Usage
-----
  python3 code/werder_havel/werder_havel_ocv_masks.py --check   # record invariants
  python3 code/werder_havel/werder_havel_ocv_masks.py --json    # the record census
Command lines are always run from the **repository root**.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, os.pardir, os.pardir, "data", "werder_havel")

MEAS_PATH = os.path.join(DATA, "werder_havel_measurements_full.txt")
SEGMENTS_PATH = os.path.join(DATA, "werder_havel_segments.txt")
WINDOWS_PATH = os.path.join(DATA, "werder_havel_windows_full.txt")
REFERENCE = os.path.join(DATA, "reference")
AOI_PATH = os.path.join(REFERENCE, "werder_havel_aoi.json")
BURSTS_PATH = os.path.join(REFERENCE, "werder_havel_bursts.json")
GEOMETRY_PATH = os.path.join(REFERENCE, "werder_havel_geometry.json")

WINDOW_SIZE = 80                     # the backend's window size of this record
WINDOW_PX = WINDOW_SIZE * WINDOW_SIZE

PEAK_FRAC = 0.30
MEDIAN_MULT = 5.0
MIN_N_MASKED = 2
RULE_TOL = 1e-9

# The asset, exactly as all 991 rows of the export name it.
ASSET = "werder"
ASSET_ID = "40f420c6-c38e-5eab-8f2b-ac00d66707d0"
ASSET_NAME = "Havelbrücke Werder (A 10)"
ASSET_TYPE = "bridge"

# The documented event of this record (Autobahn GmbH, carried by the Step-0
# hierarchy document of the external chain). The record's own strata are cut on
# exactly these two dates, which is the corroboration werder_havel_aoi.py commits.
EVENT_START = "2024-07-01"
EVENT_END = "2025-09-05"
EVENT_SOURCE = ("Autobahn GmbH via havel_werder_analysis/havel_werder_asset_hierarchy.md "
                "(Step-0 A0 document) and havel_werder_asset.json")

# Structural identity of the committed record: the facts every later step depends
# on, asserted rather than trusted (the family's ``RECORD_EXPECTED`` pattern).
RECORD_EXPECTED = {
    "rows": 991,
    "columns": 28,
    "n_dates": 491,
    "n_tracks": 8,
    "asset_id": ASSET_ID,
    "asset_name": ASSET_NAME,
    "asset_type": ASSET_TYPE,
    "first_date": "2022-07-02",
    "last_date": "2025-12-04",
    "platforms": ["S1A", "S1C"],
    "roles": {"target": 491, "reference": 250, "reference_outside": 250},
    "tracks": {"asc_146_iw1": 126, "asc_44_iw3": 118, "desc_168_iw1": 124,
               "desc_95_iw3": 123, "ref_asc_146_iw1": 126, "ref_desc_168_iw1": 124,
               "refout_asc_146_iw1": 126, "refout_desc_168_iw1": 124},
    "state_windows": {"W2_bestand": 241, "W3_sanierung": 189, "W4_post": 61},
    "window_size": WINDOW_SIZE,
    "masked_pixels": 640,       # = 0.10 * WINDOW_PX, a definition (module docstring)
    "orbit_direction_per_track": {"asc_146_iw1": "ASCENDING",
                                  "asc_44_iw3": "ASCENDING",
                                  "desc_168_iw1": "DESCENDING",
                                  "desc_95_iw3": "DESCENDING",
                                  "ref_asc_146_iw1": "ASCENDING",
                                  "ref_desc_168_iw1": "DESCENDING",
                                  "refout_asc_146_iw1": "ASCENDING",
                                  "refout_desc_168_iw1": "DESCENDING"},
    "relative_orbit_per_track": {"asc_146_iw1": 146, "asc_44_iw3": 44,
                                 "desc_168_iw1": 168, "desc_95_iw3": 95,
                                 "ref_asc_146_iw1": 146, "ref_desc_168_iw1": 168,
                                 "refout_asc_146_iw1": 146, "refout_desc_168_iw1": 168},
    "subswath_per_track": {"asc_146_iw1": "IW1", "asc_44_iw3": "IW3",
                           "desc_168_iw1": "IW1", "desc_95_iw3": "IW3",
                           "ref_asc_146_iw1": "IW1", "ref_desc_168_iw1": "IW1",
                           "refout_asc_146_iw1": "IW1", "refout_desc_168_iw1": "IW1"},
    # every reference chain re-measures the *same* dates on the same relative
    # orbit, so the whole record is 491 acquisitions, each seen eight ways
    "reference_tracks_share_target_dates": True,
}

WINDOW_NOTE = ("one 80x80 complex window per acquisition at the site point; the export "
               "committed no latitude/longitude, no segment_index and no burst id (the "
               "anchor is the OSM bridge outline, see werder_havel_aoi.py), and a "
               "window carries no pixel-spacing spec")


# ---------------------------------------------------------------------------
# The echo mask (port of carola_ocv_masks / haltern_ocv_masks)
# ---------------------------------------------------------------------------
def largest_component(mask):
    """``(n_components, largest_size)`` of a 2-D boolean mask, 8-connectivity.

    Dependency-free BFS, verbatim from ``carola_ocv_masks.largest_component``.
    """
    h, w = mask.shape
    seen = np.zeros_like(mask, dtype=bool)
    sizes = []
    for r0 in range(h):
        for c0 in range(w):
            if not mask[r0, c0] or seen[r0, c0]:
                continue
            stack, size = [(r0, c0)], 0
            seen[r0, c0] = True
            while stack:
                r, c = stack.pop()
                size += 1
                for dr in (-1, 0, 1):
                    for dc in (-1, 0, 1):
                        rr, cc = r + dr, c + dc
                        if (0 <= rr < h and 0 <= cc < w and mask[rr, cc]
                                and not seen[rr, cc]):
                            seen[rr, cc] = True
                            stack.append((rr, cc))
            sizes.append(size)
    return len(sizes), (max(sizes) if sizes else 0)


def echo_mask(z, peak_frac=PEAK_FRAC, median_mult=MEDIAN_MULT,
              min_n_masked=MIN_N_MASKED):
    """The echo mask of one complex window, and the numbers it is built from.

    The rule is the project-wide one (module docstring). ``has_echo`` is the
    ``A >= min_n_masked`` gate; the caller decides what a missing echo means
    (this package never imputes one).
    """
    z = np.asarray(z, dtype=complex)
    inten = z.real ** 2 + z.imag ** 2
    peak = float(inten.max()) if inten.size else 0.0
    median = float(np.median(inten)) if inten.size else 0.0
    if peak <= 0.0:
        mask = np.zeros(inten.shape, dtype=bool)
    else:
        mask = (inten >= peak_frac * peak) & (inten >= median_mult * median)
    n_masked = int(mask.sum())
    pr, pc = np.unravel_index(int(np.argmax(inten)), inten.shape)
    n_comp, largest = largest_component(mask)
    sel = z[mask]
    return {
        "peak": peak,
        "median": median,
        "n_masked": n_masked,
        "has_echo": n_masked >= min_n_masked,
        "max_masked": (float(inten[mask].max()) if n_masked else 0.0),
        "min_masked": (float(inten[mask].min()) if n_masked else 0.0),
        "ss_re": float(sel.real.sum()),
        "ss_im": float(sel.imag.sum()),
        "ss_abs2": float((sel.real ** 2 + sel.imag ** 2).sum()),
        "peak_row": int(pr),
        "peak_col": int(pc),
        "n_components": n_comp,
        "largest_component": largest,
        "mask": mask,
    }


def gamma2_of(ss_re, ss_im, ss_abs2, n_masked):
    """``|sum(z)|^2 / (A * sum(|z|^2))`` of a mask, clamped to [0, 1]."""
    if n_masked < 1 or ss_abs2 <= 0.0:
        return 0.0
    den = n_masked * ss_abs2
    return float(min(max((ss_re * ss_re + ss_im * ss_im) / den, 0.0), 1.0))


# ---------------------------------------------------------------------------
# The committed record
# ---------------------------------------------------------------------------
def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_measurements(path=MEAS_PATH):
    """``{id: row}`` of the measurement extract; ``id`` is ``<track>_<date>``.

    Every field stays a **string**: the export writes empty fields, and turning
    those into ``0`` is the one thing the family loaders must never do (an empty
    ``coherence_gamma2`` is a missing value, not a zero).
    """
    with open(path) as fh:
        keys = fh.readline().rstrip("\n").split("|")
        out = {}
        for line in fh:
            parts = line.rstrip("\n").split("|")
            if len(parts) != len(keys):
                continue
            row = dict(zip(keys, parts))
            out[row["id"]] = row
    return out


def read_window_ids(path=WINDOWS_PATH):
    """``{measurement_id: (w, h)}`` of the window file, without its payload.

    The window file of this record is ~250 MB of complex JSON, so the payload is
    read a line at a time and only on demand (``read_windows``).
    """
    out = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("|", 3)
            if len(parts) == 4:
                out[parts[0]] = (int(parts[1]), int(parts[2]))
    return out


def read_windows(path=WINDOWS_PATH, ids=None):
    """Yield ``(measurement_id, complex ndarray (h, w))`` of the window file.

    ``ids`` restricts the yield set. Malformed lines are skipped: the reader is
    also the acceptance test (``read_window_ids`` counts what is there).
    """
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("|", 3)
            if len(parts) != 4:
                continue
            if ids is not None and parts[0] not in ids:
                continue
            w, h = int(parts[1]), int(parts[2])
            try:
                arr = json.loads(parts[3])
            except json.JSONDecodeError:
                continue
            z = np.asarray([complex(a[0], a[1]) for a in arr], dtype=complex)
            yield parts[0], z.reshape(h, w)


def read_segments(path=SEGMENTS_PATH):
    """The committed segment list, as ``(keys, rows)``.

    Kept deliberately dumb (like ``haltern_ocv_masks``' reader): the file is an
    export, and what it *does not* carry is a fact of the record, not an error.
    """
    with open(path) as fh:
        keys = fh.readline().rstrip("\n").split("|")
        rows = [dict(zip(keys, line.rstrip("\n").split("|")))
                for line in fh if line.strip()]
    return keys, rows


# ---------------------------------------------------------------------------
# The census of the record's own columns, and its acceptance test
# ---------------------------------------------------------------------------
def census(meas=None):
    """The record's own facts, as a JSON tree (the input of every later step)."""
    meas = meas if meas is not None else read_measurements()
    rows = list(meas.values())
    keys = list(rows[0]) if rows else []
    rel = os.path.join("data", "werder_havel")
    dates = sorted({r["date"] for r in rows})
    windows = {}
    for w in sorted({r["state_window"] for r in rows}):
        ds = sorted({r["date"] for r in rows if r["state_window"] == w})
        windows[w] = {"n_dates": len(ds), "first_date": ds[0], "last_date": ds[-1]}
    per_track = {t: sorted({r["date"] for r in rows if r["track"] == t})
                 for t in sorted({r["track"] for r in rows})}
    n_masked = sorted({int(r["coherence_masked_pixels"]) for r in rows
                       if r["coherence_masked_pixels"]})
    return {
        "asset": ASSET,
        "asset_id": ASSET_ID,
        "asset_name": ASSET_NAME,
        "asset_type": ASSET_TYPE,
        "measurements": f"{rel}/werder_havel_measurements_full.txt",
        "measurements_sha256": _sha256(MEAS_PATH),
        "rows": len(rows),
        "columns": len(keys),
        "column_names": keys,
        "n_dates": len(dates),
        "first_date": dates[0],
        "last_date": dates[-1],
        "n_tracks": len(per_track),
        "tracks": {t: len(d) for t, d in per_track.items()},
        "roles": dict(sorted(Counter(r["role"] for r in rows).items())),
        "state_windows": {w: v["n_dates"] for w, v in windows.items()},
        "state_window_ranges": windows,
        "platforms": sorted({r["platform"] for r in rows}),
        "orbit_direction_per_track": {t: sorted({r["orbit_direction"] for r in rows
                                                 if r["track"] == t})[0]
                                      for t in sorted(per_track)},
        "relative_orbit_per_track": {t: int(sorted({int(r["relative_orbit"])
                                                    for r in rows
                                                    if r["track"] == t})[0])
                                     for t in sorted(per_track)},
        "subswath_per_track": {t: sorted({r["subswath"] for r in rows
                                          if r["track"] == t})[0]
                               for t in sorted(per_track)},
        "reference_tracks_share_target_dates": all(
            set(per_track.get(ref, [])) == set(per_track.get(tgt, []))
            for ref, tgt in (("ref_asc_146_iw1", "asc_146_iw1"),
                             ("refout_asc_146_iw1", "asc_146_iw1"),
                             ("ref_desc_168_iw1", "desc_168_iw1"),
                             ("refout_desc_168_iw1", "desc_168_iw1"))),
        "masked_pixels": (n_masked[0] if len(n_masked) == 1 else None),
        "masked_pixels_values": n_masked,
        "masked_pixels_frac_of_window": (round(n_masked[0] / float(WINDOW_PX), 6)
                                         if len(n_masked) == 1 else None),
        "window_size": WINDOW_SIZE,
        "window_px": WINDOW_PX,
        "empty_fields": {k: sum(1 for r in rows if r.get(k, "") == "") for k in keys},
        "segments_sha256": _sha256(SEGMENTS_PATH),
    }


def check(meas=None):
    """Re-verify the committed record offline. Returns ``(problems, facts)``."""
    facts = census(meas)
    problems = []
    for key, want in RECORD_EXPECTED.items():
        got = facts.get(key)
        if got != want:
            problems.append(f"{key}: {got!r} != {want!r}")

    # the one column that is a *definition*, not a measurement
    if facts["masked_pixels_values"] != [RECORD_EXPECTED["masked_pixels"]]:
        problems.append("coherence_masked_pixels is not one committed value: "
                        f"{facts['masked_pixels_values'][:8]}")
    elif 10 * RECORD_EXPECTED["masked_pixels"] != WINDOW_PX:
        problems.append(f"{RECORD_EXPECTED['masked_pixels']} is not 10 % of {WINDOW_PX}")

    # the strata are contiguous and they are cut on the documented event
    order = ["W2_bestand", "W3_sanierung", "W4_post"]
    r = facts["state_window_ranges"]
    if sorted(r) != order:
        problems.append(f"state windows {sorted(r)} != {order}")
    else:
        if not (r[order[0]]["last_date"] < r[order[1]]["first_date"]
                <= r[order[2]]["first_date"]):
            problems.append(f"state windows are not ordered in time: {r}")
        if r[order[2]]["first_date"] != EVENT_END:
            problems.append(f"W4_post starts {r[order[2]]['first_date']}, the "
                            f"documented reopening is {EVENT_END}")
        if not (r[order[0]]["last_date"] < EVENT_START < r[order[1]]["first_date"]):
            problems.append(f"the documented start {EVENT_START} does not fall "
                            f"between W2_bestand and W3_sanierung: "
                            f"{r[order[0]]['last_date']} .. {r[order[1]]['first_date']}")
    return problems, facts


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true",
                    help="re-verify the committed record (the default)")
    ap.add_argument("--json", action="store_true",
                    help="print the census of the record's own columns and exit")
    args = ap.parse_args(argv)
    if args.json:
        print(json.dumps(census(), indent=1, sort_keys=True))
        return 0
    problems, facts = check()
    if problems:
        print("deviations:", file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)
        raise SystemExit(1)
    print(f"ok: {facts['rows']} rows x {facts['columns']} columns, "
          f"{facts['n_tracks']} tracks, {facts['n_dates']} acquisitions "
          f"{facts['first_date']} .. {facts['last_date']}; roles {facts['roles']}; "
          f"state windows {facts['state_windows']}; platforms "
          f"{facts['platforms']}; coherence_masked_pixels "
          f"{facts['masked_pixels']} = "
          f"{facts['masked_pixels_frac_of_window']:.2f} x {WINDOW_PX} "
          f"(a definition, not a mask)")
    return 0


if __name__ == "__main__":
    sys.exit(main())


