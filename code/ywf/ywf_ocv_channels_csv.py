#!/usr/bin/env python3
"""YWF OCV-Paper — build ``data/ywf/ywf_ocv_channels.csv``.

Why this generator exists
-------------------------
The YWF record is a committed *window* export: ``data/ywf/ywf_windows_full.txt``
holds the 7x7 complex payload of every committed chip (57 windows, 38 kB,
committed in full — this site needs no separate extraction step). The mask layer
is the committed cache ``data/ywf/ywf_windows_mask_cache.txt``, written by
``ywf_ocv_masks.py --extract`` and re-verified here.

One complication is specific to this site: ``ywf_ocv_masks.echo_mask`` returns
``None`` for a window without an echo, so the cache holds **only the 17 windows
that carry an echo** (16 pre, 1 post) out of the record's 57, and its manifest
records the 16 dropped ones in ``n_without_echo``. A channel table that copied
the cache row for row would therefore hide the fact that 16 of the site's 33
unique chips show *nothing*, which is the whole point of this site's headline
(echo coverage). The CSV is consequently written one row per **unique chip** of
the window file — the same de-duplicated reference set the mask layer uses — and
``echo_mask_present`` marks the 17 that carry an echo; the build then verifies
that the cache's 17 rows are exactly those 17, statistic by statistic.

What the generator does
-----------------------
  1. reads the committed window payloads (``core.chip_table``) and the
     committed measurement extract ``ywf_measurements_full.txt``;
  2. recomputes the six mask dimensions of every unique chip with
     ``ywf_ocv_masks.vector_from_row`` — ``A``, ``D``, ``F``, ``S`` and the
     masked coherence ``gamma2`` from the re-derived mask, and ``P`` from the
     mask frequency map of the 17 echo chips;
  3. derives ``segment`` (mast section), ``section_index``, ``state``,
     ``pre_collapse``, ``day``, ``month``, ``year``, ``season``, ``echo_mode``
     and the segment-duplication provenance
     (``payload_group_size``, ``payload_shared_with_segments``);
  4. writes ``data/ywf/ywf_ocv_channels.csv`` (one row per unique chip) plus
     ``ywf_ocv_channels_meta.json`` (file sha256s, the mask rule, the guard
     counts, the definitions and the summary of every channel per state).

The build **refuses to write** if the mask layer cannot be verified: the
re-derived mask of every chip must satisfy the mask rule (checked in
``vector_from_row`` — the threshold on every masked pixel, the peak among the
masked pixels, ``A >= 2``), and every row of the committed cache must match a
unique chip of the window file in ``A``, the peak, the median, the mask pixels
and the component counts.

The coverage of that mask is summarised twice: per state (the site's headline)
and per **orbit direction** (ascending / descending). The orbit block is the
record's only second cross-cut — every unique chip is the same mast section,
``segment_resolved=False`` — and it is descriptive: both orbits straddle the
collapse date and the single post-event echo is an ascending pass. The build
cross-checks both summaries against the chip table, so a row whose
``orbit_direction`` disagrees with its measurement or whose
``echo_mask_present`` disagrees with the per-orbit count is a hard failure.

A note on ``coherence_masked_pixels`` and on the duplicated segments
--------------------------------------------------------------------
The measurement extract carries the pipeline's own mask column
``coherence_masked_pixels``. It is **not** the mask of these chips: the pipeline
computed it on the pre-fix *asset-point* chip of the whole-scene window (75,
133 ... px) while this package recomputes the echo mask of the committed 7x7
window (2..7 px, 49 px in total). The CSV keeps both (``coherence_masked_pixels``
as pipeline provenance, ``A``/``gamma2_mask`` as the recomputed mask) and reports
how many of them coincide.

24 of the record's 57 windows are byte-identical duplicates that differ only in
the request's ``segment_index`` (2 vs. 3) on one date — the record is **not**
segment-resolved. The CSV is de-duplicated on ``(asset, day, md5(payload))``
(first occurrence wins), so no duplicated chip can inflate a state, and the two
provenance columns state the group size (``2``) and the sharing segments
(``[2, 3]``) of every affected row.

Derived columns (all recomputable from the row, no lookups):
    segment        mast-section label of the request (from the window index)
    section_index  0..4 of that request (from the segments table)
    state          ``pre-collapse (healthy)`` | ``post-collapse`` (collapse date)
    pre_collapse   state == pre
    day/month/year acquisition day and its parts
    season         winter_Oct-Feb | summer_Mar-Jul | shoulder_Aug-Sep
    echo_mode      compact (<=2 px) | intermediate (3..5) | distributed (>=6)
    echo_mask_present  the chip carries an echo mask at all (A >= 2)
    payload_group_size              byte-identical windows sharing this payload
    payload_shared_with_segments    the segment indices of those windows
    A/D/F/S/P      echo-mask geometry of the chip (see ``ywf_ocv_masks``)
    gamma2_mask    |sum z|^2 / (A * sum |z|^2) over the mask

Usage:
  python3 code/ywf/ywf_ocv_channels_csv.py                 # build + verify
  python3 code/ywf/ywf_ocv_channels_csv.py --verify-only   # verify only
Command lines are always run from the **repository root**.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ywf_ocv_core as core    # noqa: E402
import ywf_ocv_masks as masks  # noqa: E402
import ywf_ocv_stats as st     # noqa: E402

DEFAULT_OUT = core.CSV_PATH
DEFAULT_META = core.CSV_META_PATH
TOL = 1e-9


def fmt(v):
    if v is None:
        return ""
    if isinstance(v, bool):
        return "1" if v else "0"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        return f"{v:.17g}"
    return str(v)


def row_of(rec, src):
    """One CSV row of a unique chip of the window file.

    ``rec`` comes from ``core.chip_table``: it carries the payload, the
    re-derived echo mask (``None`` when the chip has no echo), the mast-section
    label/index of the request and the segment-duplication provenance. ``src``
    is the measurement extract row, the source of the pipeline's own columns.
    """
    v = core.vector_of(rec)
    day = rec["date"]
    month = int(day[5:7]) if len(day) >= 7 else None
    row = {"id": rec["mid"]}
    for c in core.SOURCE_TEXT:
        row[c] = src.get(c, "")
    for c in core.SOURCE_NUMERIC:
        row[c] = core._num(src.get(c))
    for c in core.MASK_COLUMNS:          # every mask column starts empty
        row[c] = None
    row.update({
        "segment": rec["segment"],
        "section_index": rec["section_index"],
        "state": rec["state"],
        "pre_collapse": rec["state"] == core.PRE,
        "day": day, "month": month,
        "year": int(day[:4]) if len(day) >= 4 else None,
        "season": st.season_of({"month": month}),
        "echo_mode": core.echo_mode_of(v["A"]) if v else None,
        "echo_mask_present": v is not None,
        "payload_group_size": rec["payload_group_size"],
        "payload_shared_with_segments": json.dumps(
            rec["payload_shared_with_segments"]),
        "gamma2_export": core._num(src.get("coherence_gamma2")),
    })
    if v is not None:
        row.update({
            "A": v["A"], "D": v["D"], "F": v["F"], "S": v["S"],
            "P": rec.get("P"), "gamma2_mask": v["gamma2"],
            "mask_bbox_area": v["bbox_area"], "mask_row_span": v["row_span"],
            "mask_col_span": v["col_span"],
            "mask_centroid_row": v["centroid_row"],
            "mask_centroid_col": v["centroid_col"],
            "mask_peak_row": v["peak_row"], "mask_peak_col": v["peak_col"],
            "mask_n_components": v["n_components"],
            "mask_largest_component": v["largest_component"],
            "mask_fragmentation": v["fragmentation"],
        })
    return row


def cache_report(unique, cache):
    """Cross-check the committed mask cache against the re-derived masks.

    The cache of this site holds **only** the chips that carry an echo (see the
    module docstring), so the check has two parts: the cache's row count equals
    the number of chips with an echo, and every sufficient statistic of every
    cache row is reproduced by the payload. Any difference is a hard failure —
    the CSV is built from the payloads, and it may only be written when the
    committed mask layer says the same thing.
    """
    by_id = {r["mid"]: r for r in unique}
    echo = [r for r in unique if r["mask"] is not None]
    problems = []
    if len(cache) != len(echo):
        problems.append(f"cache rows {len(cache)} != chips with an echo "
                        f"{len(echo)}")
    for c in cache:
        rec = by_id.get(c["mid"])
        if rec is None:
            problems.append(f"{c['mid']}: cache row is not a unique chip")
            continue
        m = rec["mask"]
        if m is None:
            problems.append(f"{c['mid']}: cache row has no echo mask")
            continue
        checks = [("n_masked", c["n_masked"], m["A"]),
                  ("peak", c["peak"], m["peak"]),
                  ("median", c["median"], m["median"]),
                  ("max_masked", c["max_masked"], m["max_masked"]),
                  ("min_masked", c["min_masked"], m["min_masked"]),
                  ("n_components", c["n_components"], m["n_components"]),
                  ("largest_component", c["largest_component"],
                   m["largest_component"]),
                  ("ss_re", c["ss_re"], m["ss_re"]),
                  ("ss_im", c["ss_im"], m["ss_im"]),
                  ("ss_abs2", c["ss_abs2"], m["ss_abs2"])]
        for col, got, want in checks:
            if abs(float(got) - float(want)) > TOL * max(1.0, abs(float(want))):
                problems.append(f"{c['mid']}.{col}: cache {got!r} != "
                                f"recomputed {want!r}")
        if list(c["px"]) != list(m["px"]):
            problems.append(f"{c['mid']}.px: cache mask pixels differ from "
                            f"the recomputed mask")
        if c.get("P") is None or abs(c["P"] - rec["P"]) > TOL:
            problems.append(f"{c['mid']}.P: cache {c.get('P')!r} != "
                            f"{rec['P']!r}")
    return {"n_cache_rows": len(cache), "n_chips_with_echo": len(echo),
            "problems": problems}


def build(args):
    """Write the CSV + metadata of every unique chip of the window file."""
    manifest = core.load_manifest(args.manifest)
    meas = core.load_measurements(args.meas)
    seg_table = core.load_segments(args.segments)
    windows, unique = core.chip_table(args.windows)
    if not unique:
        raise SystemExit("ERROR: the committed window file is empty — "
                         "data/ywf/ywf_windows_full.txt holds the payloads")

    cache = masks.load_cache(args.cache)
    rep = cache_report(unique, cache)
    if rep["problems"]:
        for p in rep["problems"][:10]:
            print("PROBLEM:", p)
        raise SystemExit("ERROR: the committed mask cache does not match the "
                         "recomputed chips — run ywf_ocv_masks.py --extract "
                         "first")

    rows, n_export_cmp, n_export_mismatch = [], 0, 0
    for rec in unique:
        src = meas.get(rec["mid"], {})
        row = row_of(rec, src)
        exp = row["coherence_masked_pixels"]
        if exp is not None and row["A"] is not None:
            n_export_cmp += 1
            if int(exp) != int(row["A"]):
                n_export_mismatch += 1
        rows.append(row)

    rows.sort(key=lambda r: (r["day"], r["asset_name"] or "",
                             r["segment"] or "", r["id"]))
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=core.CSV_COLUMNS)
        w.writeheader()
        for r in rows:
            w.writerow({c: fmt(r.get(c)) for c in core.CSV_COLUMNS})

    meta = build_meta(rows, manifest, args, windows, unique, rep)
    meta["n_coherence_masked_pixels_compared"] = n_export_cmp
    meta["n_coherence_masked_pixels_differing"] = n_export_mismatch
    with open(args.meta, "w") as fh:
        json.dump(meta, fh, indent=1, sort_keys=True)
    print(f"written: {os.path.relpath(args.out, core.DATA)} "
          f"({len(rows)} rows x {len(core.CSV_COLUMNS)} columns)")
    print(f"  {len(windows)} committed windows -> {len(unique)} unique chips "
          f"({len(windows) - len(unique)} duplicated segments de-duplicated)")
    print(f"  echo mask: {rep['n_chips_with_echo']} of {len(unique)} unique "
          f"chips, all re-verified against the committed cache")
    bo = core.coverage_by_orbit(unique)["per_orbit"]
    print("  orbit split: " + ", ".join(
        f"{core.ORBIT_SHORT[o]} {bo[o]['n_with_echo']}/{bo[o]['n_chips']} with "
        f"an echo" for o in core.ORBITS))
    print(f"  coherence_masked_pixels is a different-generation column: "
          f"{n_export_cmp - n_export_mismatch}/{n_export_cmp} agree by chance "
          f"(it was computed on the asset-point chip, not on this window)")
    return rows, meta


def build_meta(rows, manifest, args, windows, unique, cache_rep):
    """The provenance metadata of the channel table."""
    by_state = st.by_state(rows)
    seg_counts = {}
    for r in rows:
        seg_counts[r["segment"]] = seg_counts.get(r["segment"], 0) + 1
    return {
        "generator": "code/ywf/ywf_ocv_channels_csv.py",
        "site": core.SITE,
        "event": core.COLLAPSE_DATE,
        "out_csv": os.path.basename(args.out),
        "row_unit": ("one row per unique chip (asset, day, md5(payload)) of the "
                     "committed window file; the 17 chips that carry an echo are "
                     "the rows of the committed mask cache"),
        "inputs": {
            "windows": {"file": os.path.basename(args.windows),
                        "sha256": core.sha256(args.windows)},
            "mask_cache": {"file": os.path.basename(args.cache),
                           "sha256": core.sha256(args.cache)},
            "mask_cache_manifest": {"file": os.path.basename(args.manifest),
                                    "sha256": core.sha256(args.manifest)},
            "measurements": {"file": os.path.basename(args.meas),
                             "sha256": core.sha256(args.meas)},
            "segments": {"file": os.path.basename(args.segments),
                         "sha256": core.sha256(args.segments)},
        },
        "mask_rule": manifest.get("mask_rule"),
        "mask_window": f"{manifest.get('window_size')}x"
                       f"{manifest.get('window_size')}",
        "mask_columns": list(core.MASK_COLUMNS),
        "echo_mode_cuts": core.ECHO_MODE_CUTS,
        "state_definition": f"{core.PRE} = day < {core.COLLAPSE_DATE} <= "
                            f"{core.POST}",
        "states": list(core.STATES),
        "ocv": list(core.OCV),
        "features": list(core.FEATURES),
        "sections": list(core.SECTIONS),
        "n_rows": len(rows),
        "n_columns": len(core.CSV_COLUMNS),
        "n_windows": len(windows),
        "n_unique_chips": len(unique),
        "n_duplicate_windows": len(windows) - len(unique),
        "n_pre_collapse": len(by_state[core.PRE]),
        "n_post_collapse": len(by_state[core.POST]),
        "n_echo_masks": sum(1 for r in rows if r.get("A") is not None),
        "n_rows_without_echo_mask": sum(1 for r in rows
                                        if r.get("A") is None),
        "n_P_set": sum(1 for r in rows if r.get("P") is not None),
        "sections_counts": seg_counts,
        "n_coherence_masked_pixels_matched": sum(
            1 for r in rows if r.get("coherence_masked_pixels") is not None
            and r.get("A") is not None
            and int(r["coherence_masked_pixels"]) == int(r["A"])),
        "cache_cross_check": {k: v for k, v in cache_rep.items()
                              if k != "problems"},
        "mask_cache_manifest": manifest,
        "orbits": list(core.ORBITS),
        "n_by_orbit": {core.ORBIT_SHORT[o]: sum(
            1 for r in rows if core.orbit_of(r) == o) for o in core.ORBITS},
        "coverage": core.coverage_block(unique),
        "coverage_by_orbit": core.coverage_by_orbit(unique),
        "summary_by_state": core.summary_by_state(rows),
        "date_range": [min(r["day"] for r in rows), max(r["day"] for r in rows)],
    }


def verify(args):
    """Re-derive every CSV cell from the committed window payloads."""
    rows = core.load_csv(args.out)
    meta = core.load_csv_meta(args.meta)
    windows, unique = core.chip_table(args.windows)
    by_id = {r["mid"]: r for r in unique}
    cache = masks.load_cache(args.cache)
    seg_table = core.load_segments(args.segments)
    rep = cache_report(unique, cache)
    problems, warnings = list(rep["problems"]), []
    n_cells = 0

    with open(args.out) as fh:
        header = fh.readline().rstrip("\n").split(",")
    if header != core.CSV_COLUMNS:
        problems.append("CSV header != CSV_COLUMNS")
    if len(header) != len(set(header)):
        problems.append("CSV header has duplicate column names")

    if len(rows) != meta["n_rows"]:
        problems.append(f"meta.n_rows {meta['n_rows']} != CSV rows {len(rows)}")
    if len(rows) != len(unique):
        problems.append(f"CSV rows {len(rows)} != unique chips {len(unique)}")

    n_masks = n_unmasked = n_p_set = 0
    n_export_cmp = n_export_match = 0
    for r in rows:
        rec = by_id.get(r["id"])
        if rec is None:
            problems.append(f"{r['id']}: not a unique chip of the window file")
            continue
        v = core.vector_of(rec)
        if v is None:
            n_unmasked += 1
            for col in ("A", "D", "F", "S", "P", "gamma2_mask",
                        "mask_bbox_area", "mask_fragmentation"):
                if r.get(col) is not None:
                    problems.append(f"{r['id']}.{col} set although the chip "
                                    f"has no echo mask")
        else:
            n_masks += 1
            checks = [("A", v["A"], 0.0), ("D", v["D"], TOL),
                      ("F", v["F"], 0.0), ("S", v["S"], TOL),
                      ("P", rec.get("P"), TOL),
                      ("gamma2_mask", v["gamma2"], TOL),
                      ("mask_bbox_area", v["bbox_area"], 0.0),
                      ("mask_row_span", v["row_span"], 0.0),
                      ("mask_col_span", v["col_span"], 0.0),
                      ("mask_centroid_row", v["centroid_row"], TOL),
                      ("mask_centroid_col", v["centroid_col"], TOL),
                      ("mask_n_components", v["n_components"], 0.0),
                      ("mask_largest_component", v["largest_component"], 0.0)]
            for col, want, tol in checks:
                got = r.get(col)
                n_cells += 1
                if got is None or abs(float(got) - float(want)) > tol:
                    problems.append(f"{r['id']}.{col}: {got!r} != {want!r}")
            n_p_set += 1
            if not (0.0 <= r["P"] <= 1.0):
                problems.append(f"{r['id']}.P out of [0, 1]: {r['P']!r}")
        # derived fields must be recomputable from the record alone
        if r["day"] != rec["date"]:
            problems.append(f"{r['id']}.day {r['day']!r} != {rec['date']!r}")
        if r["state"] != core.state_of(rec["date"]):
            problems.append(f"{r['id']}.state != state_of(day)")
        if bool(r.get("pre_collapse")) != (r["state"] == core.PRE):
            problems.append(f"{r['id']}.pre_collapse != state == pre")
        if core.orbit_of(r) != core.orbit_of(rec):
            problems.append(f"{r['id']}.orbit_direction {core.orbit_of(r)!r} != "
                            f"measurement {core.orbit_of(rec)!r}")
        if core.orbit_of(r) not in core.ORBITS:
            warnings.append(f"{r['id']}: orbit_direction {core.orbit_of(r)!r} "
                            f"is not one of {core.ORBITS}")
        if bool(r.get("echo_mask_present")) != (v is not None):
            problems.append(f"{r['id']}.echo_mask_present != (A is not None)")
        if r["segment"] != rec["segment"]:
            problems.append(f"{r['id']}.segment {r['segment']!r} != "
                            f"window-index label {rec['segment']!r}")
        if r.get("section_index") is None or \
                int(r["section_index"]) != rec["section_index"]:
            problems.append(f"{r['id']}.section_index {r.get('section_index')!r}"
                            f" != {rec['section_index']!r}")
        seg_rec = seg_table.get((r["asset_id"], rec["section_index"]))
        if seg_rec is None:
            warnings.append(f"{r['id']}: mast section {rec['section_index']} is "
                            f"not in the segments file")
        elif seg_rec["label"] != r["segment"]:
            problems.append(f"{r['id']}.segment {r['segment']!r} != segment "
                            f"table {seg_rec['label']!r}")
        if (r.get("echo_mode") or None) != (core.echo_mode_of(v["A"]) if v else None):
            problems.append(f"{r['id']}.echo_mode != echo_mode_of(A)")
        if r["season"] != st.season_of({"month": r["month"]}):
            problems.append(f"{r['id']}.season != season_of(month)")
        if int(r["payload_group_size"]) != rec["payload_group_size"]:
            problems.append(f"{r['id']}.payload_group_size != window group")
        if r["payload_shared_with_segments"] != json.dumps(
                rec["payload_shared_with_segments"]):
            problems.append(f"{r['id']}.payload_shared_with_segments != window "
                            f"group")
        if r.get("gamma2_mask") is not None and not (0.0 <= r["gamma2_mask"] <= 1.0):
            problems.append(f"{r['id']}.gamma2_mask out of [0, 1]")
        exp = r.get("coherence_masked_pixels")
        if exp is not None and r.get("A") is not None:
            n_export_cmp += 1
            if int(exp) == int(r["A"]):
                n_export_match += 1

    by_state = st.by_state(rows)
    for s, recs in by_state.items():
        if not recs:
            warnings.append(f"state {s} is empty")

    # the orbit cross-cut: the CSV's own echo flag per orbit must reproduce the
    # chip table's mask, and the summary must be the committed one
    by_orbit = {o: [r for r in rows if core.orbit_of(r) == o]
                for o in core.ORBITS}
    cov_orbit = core.coverage_by_orbit(unique)
    n_orbit_rows = sum(len(v) for v in by_orbit.values())
    if n_orbit_rows != len(rows):
        problems.append(f"CSV rows with a known orbit {n_orbit_rows} != CSV "
                        f"rows {len(rows)}")
    if n_orbit_rows != len(unique):
        problems.append(f"CSV rows with a known orbit {n_orbit_rows} != unique "
                        f"chips {len(unique)}")
    for o in core.ORBITS:
        want = cov_orbit["per_orbit"][o]
        n_echo = sum(1 for r in by_orbit[o] if r.get("A") is not None)
        if not by_orbit[o]:
            warnings.append(f"orbit {o} is empty")
        if len(by_orbit[o]) != want["n_chips"]:
            problems.append(f"{o}: CSV rows {len(by_orbit[o])} != chip table "
                            f"{want['n_chips']}")
        if n_echo != want["n_with_echo"]:
            problems.append(f"{o}: CSV rows with an echo {n_echo} != chip table "
                            f"{want['n_with_echo']}")
    if meta.get("coverage_by_orbit") != cov_orbit:
        problems.append("coverage_by_orbit recomputed != meta.coverage_by_orbit")

    if n_cells < 5 * len(rows):
        warnings.append(f"only {n_cells} cells compared")

    return {
        "ok": not problems,
        "n_rows": len(rows), "n_columns": len(header),
        "n_cells_checked": n_cells,
        "n_windows": len(windows), "n_unique_chips": len(unique),
        "n_echo_masks": n_masks,
        "n_rows_without_echo_mask": n_unmasked,
        "n_P_set": n_p_set,
        "n_pre_collapse": len(by_state[core.PRE]),
        "n_post_collapse": len(by_state[core.POST]),
        "n_by_orbit": {core.ORBIT_SHORT[o]: len(by_orbit[o])
                       for o in core.ORBITS},
        "coverage_by_orbit": cov_orbit,
        "sections_counts": meta.get("sections_counts"),
        "mask_rule": meta.get("mask_rule"),
        "n_cache_rows": rep["n_cache_rows"],
        "n_chips_with_echo": rep["n_chips_with_echo"],
        "n_coherence_masked_pixels_compared": n_export_cmp,
        "n_coherence_masked_pixels_matched": n_export_match,
        "n_coherence_masked_pixels_mismatched": n_export_cmp - n_export_match,
        "problems": problems[:60], "n_problems": len(problems),
        "warnings": warnings,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=DEFAULT_OUT, help="target CSV")
    ap.add_argument("--meta", default=DEFAULT_META,
                    help="provenance metadata (JSON)")
    ap.add_argument("--windows", default=core.WINDOWS_PATH,
                    help="committed window payloads (the CSV's source)")
    ap.add_argument("--cache", default=core.CACHE_PATH,
                    help="committed window cache (the mask layer)")
    ap.add_argument("--manifest", default=core.MANIFEST_PATH)
    ap.add_argument("--meas", default=core.MEAS_PATH)
    ap.add_argument("--segments", default=core.SEGMENTS_PATH)
    ap.add_argument("--verify-only", action="store_true",
                    help="write nothing, only check the existing CSV")
    args = ap.parse_args()

    if not args.verify_only:
        build(args)

    rep = verify(args)
    print(json.dumps({k: v for k, v in rep.items() if k != "problems"},
                     indent=1, sort_keys=True))
    if rep["problems"]:
        for p in rep["problems"][:60]:
            print("PROBLEM:", p)
        raise SystemExit(f"VERIFICATION FAILED ({rep['n_problems']} problems)")
    print("VERIFICATION OK")


if __name__ == "__main__":
    main()
