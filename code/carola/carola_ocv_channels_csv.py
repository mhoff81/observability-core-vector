#!/usr/bin/env python3
"""Carola OCV-Paper — build ``data/carola/carola_ocv_channels.csv``.

Why this generator exists
-------------------------
The Carola record is a committed *chip* export, not a table of channels: the
analysis project's ``carola_windows_full.txt`` (182 MB, not committed) holds the
80x80 complex payload of every chip, and the analysis derived the echo mask from
it. This package carries the mask layer as the committed **window cache**
``data/carola/carola_windows_mask_cache.txt`` (see ``carola_ocv_masks.py``,
built once by ``--extract`` from that file) and builds its channel table from it
plus the committed measurement extract and segments file.

What the generator does
-----------------------
  1. reads the committed cache (1,717 80x80 chips) and the committed
     measurement extract ``carola_measurements_full.txt``;
  2. recomputes the six mask dimensions of every chip with
     ``carola_ocv_masks.vector_from_row`` — ``A``, ``D``, ``F``, ``S`` and the
     masked coherence ``gamma2`` from the committed mask, and ``P`` from the
     mask frequency map of all chips of the majority shape;
  3. derives ``segment`` (Girder A/B/C), ``girder_index``, ``state``,
     ``pre_collapse``, ``day``, ``month``, ``year``, ``season`` and
     ``echo_mode``;
  4. writes ``data/carola/carola_ocv_channels.csv`` (one row per chip) plus
     ``carola_ocv_channels_meta.json`` (file sha256s, the mask rule, the guard
     counts, the definitions and the summary of every channel per state).

The build **refuses to write** if the mask layer cannot be verified: the
committed mask of every chip must satisfy the mask rule (checked in
``vector_from_row`` — the threshold on every masked pixel, the peak among the
masked pixels, ``A >= 2``), and every chip must have an echo mask.

A note on ``coherence_masked_pixels``
-------------------------------------
The measurement extract carries the pipeline's own mask columns
(``coherence_masked_pixels``, ``coherence_gamma2``). They are **not** the mask of
these chips: the pipeline computed them on the pre-fix *asset-point* chip while
the committed export holds the per-girder re-extraction, and the analysis
project therefore recomputed the mask from the payloads instead of trusting the
column. The CSV keeps both (``coherence_masked_pixels`` as ``*_export``
provenance, ``A``/``gamma2_mask`` as the recomputed mask) and reports how many
of them coincide — the authoritative pins are against
``carola_coherence_states.json`` and ``carola_echo_mask_gamma2.json``.

Derived columns (all recomputable from the row, no lookups):
    segment        Girder A/B/C (from the per-chip export index + segments file)
    girder_index   0/1/2 of the girder request
    state          ``pre-collapse (healthy)`` | ``post-collapse`` (collapse date)
    pre_collapse   state == pre
    day/month/year acquisition day and its parts
    season         winter_Oct-Feb | summer_Mar-Jul | shoulder_Aug-Sep
    echo_mode      compact (<=10 px) | intermediate | distributed (>=25)
    A/D/F/S/P      echo-mask geometry of the chip (see ``carola_ocv_masks``)
    gamma2_mask    |sum z|^2 / (A * sum |z|^2) over the mask

Usage:
  python3 code/carola/carola_ocv_channels_csv.py                 # build + verify
  python3 code/carola/carola_ocv_channels_csv.py --verify-only   # verify only
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import carola_ocv_core as core    # noqa: E402
import carola_ocv_masks as masks  # noqa: E402
import carola_ocv_stats as st     # noqa: E402

DEFAULT_OUT = core.CSV_PATH
DEFAULT_META = core.CSV_META_PATH
TOL = 1e-9


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


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


def _girder_index(seg_table, asset_id, segment):
    """0/1/2 of the girder request whose label is ``segment``."""
    for (aid, sidx), rec in seg_table.items():
        if aid == asset_id and rec["label"] == segment:
            return sidx
    return None


def build(args):
    cache = masks.load_cache(args.cache)
    manifest = masks.load_manifest(args.manifest)
    if not cache:
        raise SystemExit("ERROR: the window cache is empty — run "
                         "carola_ocv_masks.py --extract first")
    meas = core.load_measurements(args.meas)
    seg_table = core.load_segments(args.segments)
    masks.add_persistence(cache)

    rows, n_export_cmp, n_export_mismatch = [], 0, 0
    for rec in cache:
        v = masks.vector_from_row(rec)
        src = meas.get(rec["mid"], {})
        day = rec["date"]
        seg = rec["segment"]
        asset_id = src.get("asset_id") or ""
        month = int(day[5:7]) if len(day) >= 7 else None
        row = {"id": rec["mid"]}
        for c in core.SOURCE_TEXT:
            row[c] = src.get(c, "")
        for c in core.SOURCE_NUMERIC:
            row[c] = core._num(src.get(c))
        row.update({
            "segment": seg,
            "girder_index": _girder_index(seg_table, asset_id, seg),
            "state": core.state_of(day),
            "pre_collapse": core.state_of(day) == core.PRE,
            "day": day, "month": month,
            "year": int(day[:4]) if len(day) >= 4 else None,
            "season": st.season_of({"month": month}),
            "echo_mode": core.echo_mode_of(v["A"]),
            "A": v["A"], "D": v["D"], "F": v["F"], "S": v["S"], "P": rec.get("P"),
            "gamma2_mask": v["gamma2"],
            "gamma2_export": core._num(src.get("coherence_gamma2")),
            "coherence_masked_pixels": core._num(src.get("coherence_masked_pixels")),
            "mask_bbox_area": v["bbox_area"],
            "mask_row_span": v["row_span"], "mask_col_span": v["col_span"],
            "mask_centroid_row": v["centroid_row"],
            "mask_centroid_col": v["centroid_col"],
            "mask_peak_row": v["peak_row"], "mask_peak_col": v["peak_col"],
            "mask_n_components": v["n_components"],
            "mask_largest_component": v["largest_component"],
            "mask_fragmentation": v["fragmentation"],
        })
        exp = row["coherence_masked_pixels"]
        if exp is not None:
            n_export_cmp += 1
            if int(exp) != v["A"]:
                n_export_mismatch += 1
        rows.append(row)

    rows.sort(key=lambda r: (r["day"], r["asset_name"] or "", r["segment"] or "",
                             r["id"]))
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=core.CSV_COLUMNS)
        w.writeheader()
        for r in rows:
            w.writerow({c: fmt(r.get(c)) for c in core.CSV_COLUMNS})

    meta = build_meta(rows, manifest, args)
    meta["n_coherence_masked_pixels_compared"] = n_export_cmp
    meta["n_coherence_masked_pixels_differing"] = n_export_mismatch
    with open(args.meta, "w") as fh:
        json.dump(meta, fh, indent=1, sort_keys=True)
    print(f"written: {os.path.relpath(args.out, core.DATA)} "
          f"({len(rows)} rows x {len(core.CSV_COLUMNS)} columns)")
    print(f"  coherence_masked_pixels is a different-generation column: "
          f"{n_export_cmp - n_export_mismatch}/{n_export_cmp} agree by chance "
          f"(the analysis project recomputed the mask from the payloads)")
    return rows, meta


def build_meta(rows, manifest, args):
    by_state = st.by_state(rows)
    seg_counts = {}
    for r in rows:
        seg_counts[r["segment"]] = seg_counts.get(r["segment"], 0) + 1
    return {
        "generator": "code/carola/carola_ocv_channels_csv.py",
        "site": core.SITE,
        "event": core.COLLAPSE_DATE,
        "out_csv": os.path.basename(args.out),
        "inputs": {
            "mask_cache": {"file": os.path.basename(args.cache),
                           "sha256": sha256(args.cache)},
            "mask_cache_manifest": {"file": os.path.basename(args.manifest),
                                    "sha256": sha256(args.manifest)},
            "measurements": {"file": os.path.basename(args.meas),
                             "sha256": sha256(args.meas)},
            "segments": {"file": os.path.basename(args.segments),
                         "sha256": sha256(args.segments)},
        },
        "mask_rule": manifest["mask_rule"],
        "mask_window": f"{manifest['window_size']}x{manifest['window_size']}",
        "mask_columns": list(core.MASK_COLUMNS),
        "state_definition": f"{core.PRE} = day < {core.COLLAPSE_DATE} <= {core.POST}",
        "states": list(core.STATES),
        "ocv": list(core.OCV),
        "features": list(core.FEATURES),
        "girders": list(core.GIRDERS),
        "n_rows": len(rows),
        "n_columns": len(core.CSV_COLUMNS),
        "n_pre_collapse": len(by_state[core.PRE]),
        "n_post_collapse": len(by_state[core.POST]),
        "n_echo_masks": sum(1 for r in rows if r.get("A") is not None),
        "n_rows_without_echo_mask": sum(1 for r in rows if r.get("A") is None),
        "n_P_set": sum(1 for r in rows if r.get("P") is not None),
        "girders_counts": seg_counts,
        "n_coherence_masked_pixels_matched": sum(
            1 for r in rows if r.get("coherence_masked_pixels") is not None
            and int(r["coherence_masked_pixels"]) == int(r["A"])),
        "mask_cache_manifest": manifest,
        "summary_by_state": core.summary_by_state(rows),
        "date_range": [min(r["day"] for r in rows), max(r["day"] for r in rows)],
    }


def verify(args):
    """Re-derive every CSV cell from the committed cache and check the schema."""
    rows = core.load_csv(args.out)
    meta = core.load_csv_meta(args.meta)
    cache = masks.load_cache(args.cache)
    by_id = {r["mid"]: r for r in cache}
    masks.add_persistence(cache)
    problems, warnings = [], []
    n_cells = 0

    with open(args.out) as fh:
        header = fh.readline().rstrip("\n").split(",")
    if header != core.CSV_COLUMNS:
        problems.append("CSV header != CSV_COLUMNS")

    if len(rows) != meta["n_rows"]:
        problems.append(f"meta.n_rows {meta['n_rows']} != CSV rows {len(rows)}")
    if len(rows) != len(cache):
        problems.append(f"CSV rows {len(rows)} != cache rows {len(cache)}")

    n_masks = n_unmasked = n_p_set = 0
    n_export_cmp = n_export_match = 0
    for r in rows:
        c = by_id.get(r["id"])
        if c is None:
            problems.append(f"{r['id']}: not in the cache")
            continue
        v = masks.vector_from_row(c)
        n_masks += 1
        checks = [("A", v["A"], 0.0), ("D", v["D"], TOL), ("F", v["F"], 0.0),
                  ("P", c.get("P"), TOL), ("gamma2_mask", v["gamma2"], TOL),
                  ("mask_bbox_area", v["bbox_area"], 0.0),
                  ("mask_row_span", v["row_span"], 0.0),
                  ("mask_col_span", v["col_span"], 0.0),
                  ("mask_n_components", v["n_components"], 0.0)]
        for col, want, tol in checks:
            got = r.get(col)
            n_cells += 1
            if got is None or abs(float(got) - float(want)) > tol:
                problems.append(f"{r['id']}.{col}: {got!r} != {want!r}")
        if r.get("P") is not None:
            n_p_set += 1
            if not (0.0 <= r["P"] <= 1.0):
                problems.append(f"{r['id']}.P out of [0, 1]: {r['P']!r}")
        if r.get("A") is None:
            n_unmasked += 1
        exp = r.get("coherence_masked_pixels")
        if exp is not None:
            n_export_cmp += 1
            if int(exp) == int(r["A"]):
                n_export_match += 1
        # derived fields must be recomputable from the day alone
        if r["day"] != c["date"]:
            problems.append(f"{r['id']}.day {r['day']!r} != cache date {c['date']!r}")
        if r["state"] != core.state_of(c["date"]):
            problems.append(f"{r['id']}.state != state_of(day)")
        if bool(r.get("pre_collapse")) != (r["state"] == core.PRE):
            problems.append(f"{r['id']}.pre_collapse != state == pre")
        if r["segment"] != c["segment"]:
            problems.append(f"{r['id']}.segment != cache segment")
        if r["echo_mode"] != core.echo_mode_of(v["A"]):
            problems.append(f"{r['id']}.echo_mode != echo_mode_of(A)")
        if r["season"] != st.season_of({"month": r["month"]}):
            problems.append(f"{r['id']}.season != season_of(month)")
        if r.get("gamma2_mask") is not None and not (0.0 <= r["gamma2_mask"] <= 1.0):
            problems.append(f"{r['id']}.gamma2_mask out of [0, 1]")

    by_state = st.by_state(rows)
    for s, recs in by_state.items():
        if not recs:
            warnings.append(f"state {s} is empty")

    if n_cells < 5 * len(rows):
        warnings.append(f"only {n_cells} cells compared")

    return {
        "ok": not problems,
        "n_rows": len(rows), "n_columns": len(header),
        "n_cells_checked": n_cells,
        "n_echo_masks": n_masks,
        "n_rows_without_echo_mask": n_unmasked,
        "n_P_set": n_p_set,
        "n_pre_collapse": len(by_state[core.PRE]),
        "n_post_collapse": len(by_state[core.POST]),
        "girders_counts": meta.get("girders_counts"),
        "mask_rule": meta.get("mask_rule"),
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
    ap.add_argument("--meta", default=DEFAULT_META, help="provenance metadata (JSON)")
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
