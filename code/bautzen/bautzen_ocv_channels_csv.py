#!/usr/bin/env python3
"""Bautzen OCV-Paper — enumerate the four rect series into one channel table.

Why this generator exists
-------------------------
The mask dimensions ``A`` (area), ``D`` (density), ``F`` (fragmentation),
``S`` (centroid<->deck-band-peak shift) and ``P`` (persistence) and the deck-band
``gamma^2`` are functions of the *deck-edge mask* of one date. The committed
reference ``bautzen_deck_edge_mask.json`` stores the mask statistics per date but
not the OCV vector; ``bautzen_deck_edge_state.json`` stores ``gamma^2`` and the
state metrics but again not a flat per-date table. This generator recomputes all
of it from the four committed rect files, with exactly the thresholds and the
order of operations of ``bautzen_ocv_masks.py``, and writes the result as
*columns* into ``data/bautzen/bautzen_ocv_channels.csv``. Afterwards that CSV
(plus its meta JSON) is the only input of the figure pipeline.

Unlike the LUMO generator this one needs no external cache: the four rect files
are committed next to the CSV (see ``data/README.md``), so a fresh clone rebuilds
the table offline. The committed ``mask_pixels`` column carries the mask itself
(the masked pixels, row-major), which is what lets the figure pipeline rebuild
per-state persistence without touching the rect files.

Derivation (identical to the origin, see ``bautzen_ocv_masks.py``):
  gamma2 = |sum z|^2 / (N sum |z|^2) on the deck-edge mask, floor-corrected
  A = n_masked
  D = n_masked / bbox_area
  F = number of 8-neighbourhood components
  S = ||mask centroid - deck band peak|| in px
  P = mean pixel frequency across the dates of the same series

Usage:
  python3 code/bautzen/bautzen_ocv_channels_csv.py                # build + verify
  python3 code/bautzen/bautzen_ocv_channels_csv.py --verify-only  # verify only
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import bautzen_ocv_core as core    # noqa: E402
import bautzen_ocv_masks as masks  # noqa: E402

DEFAULT_CSV = os.path.join(core.DATA, "bautzen_ocv_channels.csv")
DEFAULT_META = os.path.join(core.DATA, "bautzen_ocv_channels_meta.json")
REF_MASK = os.path.join(core.REF, core.REFERENCE_FILES["deck_edge_mask"])
REF_STATE = os.path.join(core.REF, core.REFERENCE_FILES["deck_edge_state"])

CSV_COLUMNS = [
    "site", "series", "rect_file", "acquisition_date", "orbit_direction",
    "month", "state",
    "band_contrast", "peak_row", "peak_row_offset", "n_masked",
    "n_masked_house", "n_components", "largest_component", "has_mask",
    "centroid_row", "centroid_col", "band_peak_r", "band_peak_c",
    "A", "D", "F", "S", "P", "P_frac_ge_50pct", "P_frac_ge_30pct",
    "gamma2_band_raw", "gamma2_band_floor_corrected", "n_band",
    "gamma2_house_raw", "gamma2_house_floor_corrected", "n_house",
    "mask_rows", "mask_cols", "bbox_area", "row_span", "col_span",
    "mask_pixels",
]
MASK_JSON_KEYS = ["anchor_row", "peak_row", "peak_row_offset", "n_masked",
                  "band_contrast", "n_masked_house", "n_components",
                  "largest_component", "has_mask", "centroid_row",
                  "centroid_col", "date"]
STATE_JSON_KEYS = ["date", "state", "band_contrast", "peak_row",
                   "peak_row_offset", "n_masked", "has_mask", "centroid_row",
                   "gamma2_band_raw", "gamma2_band_floor_corrected", "n_band",
                   "gamma2_house_raw", "gamma2_house_floor_corrected", "n_house"]


def _fmt(v):
    if v is None:
        return ""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, float):
        return repr(v)
    return str(v)


def _pixels_str(mask):
    rr, cc = np.nonzero(mask)
    return ";".join(f"{int(r)}:{int(c)}" for r, c in zip(rr, cc))


def build(args):
    """One mask per date of every series; writes the CSV and returns the meta."""
    rows, series_block, sources = [], {}, {}
    problems_guard = []   # cross-check channel path vs. analysis path (see below)
    for lab in core.SERIES_ORDER:
        fn = core.SERIES[lab]
        path = os.path.join(core.DATA, fn)
        if not os.path.exists(path):
            raise SystemExit(f"ERROR: rect file missing: {path}")
        sources[fn] = core.sha256(path)
        geom = core.geom_of(lab)
        recs = masks.load_complex(path)

        # Channel side: one enriched record per date (mask + vector + gamma^2).
        infos = []
        for date, z in recs:
            inten = np.abs(z) ** 2
            info, full = masks.date_mask(inten)
            info["date"] = date
            info["series"] = lab
            info["state"] = core.state_of(date, geom)
            # gamma^2 on the deck band mask (structurally attributed, deck-edge
            # statement) ...
            raw_b, corr_b, n_b = masks.gamma2(z[full])
            info["gamma2_band_raw"] = raw_b
            info["gamma2_band_floor_corrected"] = corr_b
            info["n_band"] = n_b
            # ... and on the signal-based house/peak mask (no deck-edge statement)
            raw_h, corr_h, n_h = masks.gamma2(z[masks.house_mask(inten)])
            info["gamma2_house_raw"] = raw_h
            info["gamma2_house_floor_corrected"] = corr_h
            info["n_house"] = n_h
            info["mask_rows"], info["mask_cols"] = (int(s) for s in full.shape)
            info["mask"] = full
            info["shape"] = full.shape
            infos.append(info)
        infos.sort(key=lambda d: d["date"])
        # A date without a deck-edge pixel has no anchor, no peak and no
        # centroid. The committed mask reference stores those as an explicit
        # ``null`` instead of omitting the key, so the same convention is used
        # here: a key is never absent, only empty. Without this the verifier
        # could not tell "not computed" from "computed as nothing".
        for d in infos:
            for k in MASK_JSON_KEYS:
                d.setdefault(k, None)
        persistence = masks.add_persistence(infos)

        # Mask side: the verbatim series summary + union maps (these are the
        # numbers of the committed mask reference). NOTE the input: analyse_series
        # is the port of the inventory script and takes *intensities*, never the
        # complex chips (a complex median sorts lexicographically and would give a
        # different mask).
        recs_i = [(date, np.abs(z) ** 2) for date, z in recs]
        per_date, summary, union50, union30, counts = masks.analyse_series(recs_i)
        # Cross-check: the mask of the analysis path must be the mask of the
        # channel path, date by date (two independent traversals, one result).
        side = {d["date"]: d for d in per_date}
        for d in infos:
            p = side[d["date"]]
            for k in p:
                if k not in d and p[k] is None:
                    continue          # absent here == explicit null there
                if not _close(d.get(k), p[k]):
                    problems_guard.append(
                        f"{lab} {d['date']}: channel path {k} = {d.get(k)!r} "
                        f"!= analysis path {p[k]!r}")
        series_block[lab] = {
            "file": fn, "geom": geom, "n_dates": len(infos),
            "n_by_state": {s: sum(1 for d in infos if d["state"] == s)
                           for s in core.STATES},
            "dates_by_state": {s: [d["date"] for d in infos if d["state"] == s]
                               for s in core.STATES},
            "per_date": infos, "summary": summary, "union50": union50,
            "union30": union30, "counts": counts, "persistence": persistence,
        }
        for d in infos:
            row = {
                "site": core.SITE, "series": lab, "rect_file": fn,
                "acquisition_date": d["date"],
                "orbit_direction": geom,
                "month": core.month_of(d["date"]), "state": d["state"],
                "mask_pixels": _pixels_str(d["mask"]),
            }
            for col in CSV_COLUMNS:
                if col in row:
                    continue
                row[col] = d.get(col)
            rows.append({c: _fmt(row.get(c)) for c in CSV_COLUMNS})
        print(f"  {lab}: {len(infos)} dates, "
              f"{series_block[lab]['n_by_state']}, "
              f"deck edge detectable = {summary['deck_edge_detectable']}")

    series_block["_cross_check"] = problems_guard
    if problems_guard:
        for p in problems_guard:
            print(f"  CROSS-CHECK PROBLEM: {p}")

    meta = {
        "generator": "code/bautzen/bautzen_ocv_channels_csv.py",
        "out_csv": os.path.basename(args.out),
        "site": core.SITE,
        "reference_source": {
            "mask": "data/bautzen/reference/" + core.REFERENCE_FILES["deck_edge_mask"],
            "state": "data/bautzen/reference/" + core.REFERENCE_FILES["deck_edge_state"],
        },
        "rect_files": sorted(sources),
        "rect_file_sha256": sources,
        "mask_definition": {
            "layer": "deck_edge",
            "center_row": masks.CENTER_ROW, "center_col": masks.CENTER_COL,
            "anchor_half": masks.ANCHOR_HALF, "col_half": masks.COL_HALF,
            "band_half": masks.BAND_HALF, "k_local": masks.K_LOCAL,
            "min_n_masked": masks.MIN_N_MASKED,
            "contrast_detect": masks.CONTRAST_DETECT,
            "house_peak_frac": masks.GAMMA2_PEAK_FRAC,
            "house_median_mult": masks.GAMMA2_MEDIAN_MULT,
            "house_min_masked": masks.GAMMA2_MIN_MASKED,
            "low_n_gate": masks.LOW_N_GATE,
        },
        "dimensions": {
            "gamma2": "gamma2_band_floor_corrected: |sum z|^2/(N sum|z|^2) on the "
                      "deck-edge mask, (N*raw-1)/(N-1), no clipping",
            "A": "n_masked (deck-edge pixels of the anchored band)",
            "D": "n_masked / bbox_area",
            "F": "number of 8-neighbourhood components of the mask",
            "S": "distance mask centroid <-> deck band peak (px) -- the global "
                 "chip peak is clutter at this site, see bautzen_ocv_masks.py",
            "P": "mean pixel frequency inside the own mask, frequency map of the "
                 "dates of the same series",
        },
        "states": {s: {"window": w} for s, w in core.STATE_WINDOWS.items()},
        "series": {lab: {k: v for k, v in blk.items()
                         if k not in ("counts", "per_date")
                         and not k.startswith("_")}
                   for lab, blk in series_block.items()
                   if not lab.startswith("_")},
        "channel_columns": CSV_COLUMNS,
    }
    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    with open(args.meta, "w") as fh:
        json.dump(meta, fh, indent=1, sort_keys=True, default=float)
    return meta, series_block, rows


def _close(a, b, tol=1e-9):
    """Structural comparison of a recomputed value with a committed one.

    Scalars must agree within ``tol``; strings exactly; ``None`` only with
    ``None``. Lists/tuples and dicts (the reference stores a few, e.g. pixel
    lists or per-block sub-dicts) are compared element by element, so the
    verifier can walk a reference row key by key without having to know which
    key holds what -- and without ever calling float() on a container.
    """
    if isinstance(a, (list, tuple)) or isinstance(b, (list, tuple)):
        return (isinstance(a, (list, tuple)) and isinstance(b, (list, tuple))
                and len(a) == len(b)
                and all(_close(x, y, tol) for x, y in zip(a, b)))
    if isinstance(a, dict) or isinstance(b, dict):
        return (isinstance(a, dict) and isinstance(b, dict)
                and set(a) == set(b)
                and all(_close(a[k], b[k], tol) for k in a))
    if a is None or b is None:
        return a is b or a == b
    if isinstance(a, bool) or isinstance(b, bool):
        return bool(a) == bool(b)
    if isinstance(a, str) or isinstance(b, str):
        return a == b
    try:
        return abs(float(a) - float(b)) <= tol
    except (TypeError, ValueError):
        return a == b


def verify(args, series_block, rows):
    """Hard self-check of the written table against the two committed references.

    ok == False => abort (exit code != 0). Checks, in this order:
      1. the CSV itself (header, one row per date, A/D/P consistency),
      2. the mask reference (config, per-date mask fields, series summary and
         both union maps, pixel-exact),
      3. the state reference (state assignment, per-date gamma^2, the per-state
         metrics incl. the union maps and the Jaccard overlaps).
    """
    problems = list(series_block.get("_cross_check") or [])
    warnings = []
    mask_ref = json.load(open(REF_MASK))
    state_ref = json.load(open(REF_STATE))

    # ---- 1. the CSV ------------------------------------------------------
    with open(args.out, newline="") as fh:
        reader = csv.DictReader(fh)
        fields = list(reader.fieldnames or [])
        got = list(reader)
    if fields != CSV_COLUMNS:
        problems.append(f"CSV header differs from CSV_COLUMNS: {fields}")
    if len(got) != len(rows):
        problems.append(f"{len(got)} CSV rows instead of {len(rows)}")
    for r in got:
        if r["A"] != r["n_masked"]:
            problems.append(f"A != n_masked at {r['series']} {r['acquisition_date']}")
            break
    for r in got:
        if r["D"] == "" or r["bbox_area"] == "":
            continue
        if float(r["bbox_area"]) > 0 and float(r["D"]) != float(r["A"]) / float(r["bbox_area"]):
            problems.append(f"D != A/bbox_area at {r['series']} {r['acquisition_date']}")
            break
    for r in got:
        if r["P"] and not (0.0 <= float(r["P"]) <= 1.0):
            problems.append(f"P outside [0,1] at {r['series']} {r['acquisition_date']}")
            break

    # ---- 2. the mask reference ------------------------------------------
    if mask_ref["config"] != masks.CONFIG:
        problems.append("mask reference config != masks.CONFIG")
    for lab in core.SERIES_ORDER:
        ref = mask_ref["series"][lab]
        blk = series_block[lab]
        if ref["file"] != blk["file"]:
            problems.append(f"{lab}: file {blk['file']} != {ref['file']}")
        if len(ref["per_date"]) != blk["n_dates"]:
            problems.append(f"{lab}: {blk['n_dates']} dates != {len(ref['per_date'])}")
        mine = {d["date"]: d for d in blk["per_date"]}
        for rd in ref["per_date"]:
            md = mine.get(rd["date"])
            if md is None:
                problems.append(f"{lab}: date {rd['date']} missing")
                continue
            # The reference dictates the key set: every key it stores must be
            # present here with the same value. (Dates without a deck-edge pixel
            # carry no centroid/peak, so the reference rows are not all equally
            # wide -- iterating over MASK_JSON_KEYS would raise a KeyError.)
            for k in rd:
                if k not in md:
                    problems.append(f"{lab} {rd['date']}: mask.{k} missing "
                                    f"(committed {rd[k]!r})")
                elif not _close(md[k], rd[k]):
                    problems.append(f"{lab} {rd['date']}: mask.{k} = {md[k]!r} "
                                    f"!= committed {rd[k]!r}")
        for k in ref["summary"]:
            if not _close(blk["summary"].get(k), ref["summary"][k]):
                problems.append(f"{lab}: summary.{k} = {blk['summary'].get(k)!r} "
                                f"!= committed {ref['summary'][k]!r}")
        for key in ("union50", "union30"):
            mine_px = [[r, c] for r, c, _f in blk[key]]
            ref_px = [[r, c] for r, c, _f in ref[key]]
            if mine_px != ref_px:
                problems.append(f"{lab}: {key} pixels differ from the committed list")
            elif ([f for _r, _c, f in blk[key]]
                  != [f for _r, _c, f in ref[key]]):
                problems.append(f"{lab}: {key} frequencies differ (pixel order equal)")
    return _verify_state(state_ref, series_block, problems, warnings, got)


def _verify_state(state_ref, series_block, problems, warnings, got):
    """Part 3 of ``verify``: the state reference (see its docstring)."""
    for lab in core.SERIES_ORDER:
        ref = state_ref["series"][lab]
        blk = series_block[lab]
        for k in ("file", "geom", "n_dates"):
            if not _close(blk[k], ref[k]):
                problems.append(f"{lab}: {k} = {blk[k]!r} != committed {ref[k]!r}")
        if blk["n_by_state"] != ref["n_by_state"]:
            problems.append(f"{lab}: n_by_state {blk['n_by_state']} "
                            f"!= {ref['n_by_state']}")
        if blk["dates_by_state"] != ref["dates_by_state"]:
            problems.append(f"{lab}: dates_by_state differ from the committed file")
        mine = {d["date"]: d for d in blk["per_date"]}
        for rd in ref["per_date"]:
            md = mine.get(rd["date"])
            if md is None:
                problems.append(f"{lab}: state date {rd['date']} missing")
                continue
            for k in rd:
                if k not in md:
                    problems.append(f"{lab} {rd['date']}: {k} missing "
                                    f"(committed {rd[k]!r})")
                elif not _close(md[k], rd[k]):
                    problems.append(f"{lab} {rd['date']}: {k} = {md[k]!r} "
                                    f"!= committed {rd[k]!r}")
        # Per-state metrics, rebuilt from the per-date masks of this series.
        for stt in list(core.STATES) + ["ALL"]:
            sub = blk["per_date"] if stt == "ALL" else [
                d for d in blk["per_date"] if d["state"] == stt]
            m = masks.state_metrics(sub, masks.counts_of(sub))
            refm = core.translate(ref["metrics"][stt])
            for k in refm:
                if k in ("u50", "u30"):
                    continue
                if not _close(m.get(k), refm[k]):
                    problems.append(f"{lab} {stt}: metrics.{k} = {m.get(k)!r} "
                                    f"!= committed {refm[k]!r}")
            for k in ("u50", "u30"):
                if ([[r, c] for r, c, _f in m[k]]
                        != [[r, c] for r, c, _f in refm[k]]):
                    problems.append(f"{lab} {stt}: {k} pixels differ "
                                    f"from the committed list")
        # Jaccard overlaps between the >= 50 % unions of the states.
        sets = {}
        for stt in core.STATES:
            sub = [d for d in blk["per_date"] if d["state"] == stt]
            u50, _u30 = masks.unions(masks.counts_of(sub), len(sub))
            sets[stt] = {(r, c) for r, c, _f in u50}
        for key, want in ref["jaccard"].items():
            a, b = key.split("|")
            # An empty union (both states without a single mask pixel) is
            # undefined, not 0 -- the committed reference stores null there.
            union = sets[a] | sets[b]
            got_j = (len(sets[a] & sets[b]) / len(union)) if union else None
            if not _close(got_j, want):
                problems.append(f"{lab}: jaccard {key} = {got_j!r} "
                                f"!= committed {want!r}")
    if not got:
        warnings.append("the written table has no rows")
    return {
        "ok": not problems,
        "n_rows": len(got),
        "n_series": len(core.SERIES_ORDER),
        "n_dates": {lab: series_block[lab]["n_dates"] for lab in core.SERIES_ORDER},
        "problems": problems,
        "warnings": warnings,
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=DEFAULT_CSV, help="target CSV (enumerated)")
    ap.add_argument("--meta", default=DEFAULT_META, help="provenance metadata (JSON)")
    ap.add_argument("--verify-only", action="store_true",
                    help="only check (the files are rewritten with identical content)")
    args = ap.parse_args(argv)

    meta, series_block, rows = build(args)
    print(f"written: {os.path.relpath(args.out)} ({len(rows)} rows)")
    print(f"written: {os.path.relpath(args.meta)}")
    rep = verify(args, series_block, rows)
    print(json.dumps(rep, indent=1, sort_keys=True))
    if not rep["ok"]:
        raise SystemExit("VERIFICATION FAILED")
    print("VERIFICATION OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())



