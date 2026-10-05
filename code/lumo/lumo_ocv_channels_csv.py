#!/usr/bin/env python3
"""LUMO OCV-Paper — copy ``lumo_channels.csv`` and extend it with A/D/F/S/P.

Why this generator exists
-------------------------
The mask dimensions ``A`` (area), ``D`` (density), ``F`` (fragmentation),
``S`` (centroid<->peak shift) and ``P`` (persistence) are functions of the
*echo mask* of one overpass. No committed JSON of this stack stores them per
overpass: ``lumo_echo_mask_vector.json`` / ``lumo_emv_multivariate_cv.json``
store only aggregates (``per_dimension``), ``lumo_tower_coherence_states.json``
only ``n_masked`` and ``peak_intensity``. The raw source is the machine-local
burst cache ``monthly_bursts/<month>/<key>/{meta.json,strip.bin}`` (not
committed, ~100 MB, a pipeline artefact).

This generator reads the cache **once**, computes the mask geometry with the
same thresholds and definitions as ``analyze_lumo_echo_mask_vector.py`` and
writes the result as *columns* into ``data/lumo/lumo_ocv_channels.csv``. Afterwards
that CSV is the only input of the figure pipeline — the cache is never needed
again.

Derivation (identical to the origin, see ``lumo_ocv_masks.py``):
  A = n_masked
  D = n_masked / bbox_area
  F = number of 8-neighbourhood components
  S = ||centroid - peak|| in px
  P = mean pixel frequency across the 177 dates of the majority chip shape
      (the 1 date with a deviating shape gets P = empty)

Only the 178 ``coherence`` rows get values; the remaining rows of the copied
CSV stay empty in the new columns.

Usage:
  python3 code/lumo/lumo_ocv_channels_csv.py                # build + verify
  python3 code/lumo/lumo_ocv_channels_csv.py --verify-only  # verify only
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, os.pardir, os.pardir, "data", "lumo")
REF = os.path.join(DATA, "reference")

sys.path.insert(0, HERE)
import lumo_ocv_masks as masks  # noqa: E402

DEFAULT_CSV = os.path.join(DATA, "lumo_channels.csv")
DEFAULT_OUT = os.path.join(DATA, "lumo_ocv_channels.csv")
DEFAULT_META = os.path.join(DATA, "lumo_ocv_channels_meta.json")
DEFAULT_BURST_DIR = os.path.abspath(os.path.join(
    HERE, os.pardir, os.pardir, os.pardir, os.pardir,
    "lumo_dam6_analysis", "monthly_bursts"))

COH_JSON = os.path.join(REF, "lumo_tower_coherence_states.json")
EMV_JSON = os.path.join(REF, "lumo_echo_mask_vector.json")

# Appended columns (order = documentation):
#   A/D/F/S/P            -> the EMV vector used by the figure pipeline
#   mask_*/bbox_*/peak_* -> auditable ingredients (sizes in px) that allow
#                           recomputing D, F and S directly from the CSV
NEW_COLS = ["A", "D", "F", "S", "P",
            "mask_rows", "mask_cols", "bbox_area",
            "centroid_r", "centroid_c", "peak_r", "peak_c"]
DIMS = ["A", "D", "F", "S", "P"]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_committed():
    with open(COH_JSON) as fh:
        return json.load(fh)


def committed_keys(coh):
    return {(b["date"], b["orbit"]) for b in coh["bursts"]}


def walk_cache(burst_dir, keys):
    """One vector per cache directory of a committed overpass (cache order)."""
    rows, n_skipped = [], 0
    for month in sorted(os.listdir(burst_dir)):
        md = os.path.join(burst_dir, month)
        if not os.path.isdir(md):
            continue
        for key in sorted(os.listdir(md)):
            bd = os.path.join(md, key)
            meta_path = os.path.join(bd, "meta.json")
            strip_path = os.path.join(bd, "strip.bin")
            if not (os.path.isfile(meta_path) and os.path.isfile(strip_path)):
                continue
            with open(meta_path) as fh:
                meta = json.load(fh)
            if (meta.get("acquisition_date"), meta.get("orbit")) not in keys:
                continue
            try:
                z, _w, _h = masks.read_complex_bin(strip_path)
                vec = masks.vector_from_strip(z)
            except Exception as exc:  # pragma: no cover - cache integrity
                print(f"skip {month}/{key}: {exc}", file=sys.stderr)
                n_skipped += 1
                continue
            if vec is None:
                n_skipped += 1
                continue
            vec.update({"date": meta.get("acquisition_date"), "orbit": meta.get("orbit"),
                        "polarisation": meta.get("polarisation"),
                        "burst_id": meta.get("burst_id")})
            rows.append(vec)
    return rows, n_skipped


def committed_rows(raw, coh):
    """Restrict to the 178 committed overpasses, VV preferred, in COH-JSON order.

    Verbatim rebuild of ``analyze_lumo_echo_mask_vector.committed_rows`` +
    ``analyze_lumo_tower_coherence.dedup_overpasses``. The order deliberately is
    the order of the committed JSON, because bootstrap CIs (seed 7) and
    Jonckheere permutations depend on the order.
    """
    by_key = {}
    for r in raw:
        by_key.setdefault((r["date"], r["orbit"]), []).append(r)

    def is_vv(b):
        return str(b.get("polarisation", "")).lower() == "vv"

    rows, n_mismatch, missing = [], 0, []
    for ref_row in coh["bursts"]:
        cands = by_key.get((ref_row["date"], ref_row["orbit"]))
        if not cands:
            missing.append((ref_row["date"], ref_row["orbit"]))
            continue
        vv = [c for c in cands if is_vv(c)]
        c = dict(vv[0] if vv else cands[0])
        if c["A"] != ref_row["n_masked"]:
            n_mismatch += 1
        c["damage_label"] = ref_row["damage_label"]
        c["committed_gamma2"] = ref_row["gamma2"]
        c["committed_peak_intensity"] = ref_row["peak_intensity"]
        c["committed_burst_id"] = ref_row["burst_id"]
        rows.append(c)
    return rows, {"n_checked": len(rows), "n_mismatch": n_mismatch, "missing": missing}


def fmt_num(v):
    if v is None:
        return ""
    if isinstance(v, float):
        return repr(v)
    return str(v)


def build(args):
    coh = load_committed()
    keys = committed_keys(coh)
    raw, n_skipped = walk_cache(args.burst_dir, keys)
    rows, check = committed_rows(raw, coh)
    persistence_meta = masks.add_persistence(rows)
    if check["missing"]:
        raise SystemExit(f"ERROR: cache incomplete for {len(check['missing'])} "
                         f"committed overpasses: {check['missing'][:3]}")

    by_bid = {r["committed_burst_id"]: r for r in rows}
    by_keymap = {(r["date"], r["orbit"]): r for r in rows}

    with open(args.csv, newline="") as fh:
        reader = csv.DictReader(fh)
        fieldnames = list(reader.fieldnames)
        src_rows = list(reader)
    for col in NEW_COLS:
        if col in fieldnames:
            raise SystemExit(f"ERROR: column {col} already exists in {args.csv}")

    n_filled, n_by_bid, n_by_key, n_unmatched = 0, 0, 0, 0
    for r in src_rows:
        if not (r.get("coherence_gamma2") or "").strip():
            for col in NEW_COLS:
                r[col] = ""
            continue
        m = by_bid.get(r.get("burst_id"))
        if m is not None:
            n_by_bid += 1
        else:
            m = by_keymap.get((r.get("acquisition_date"), r.get("orbit_direction")))
            if m is not None:
                n_by_key += 1
        if m is None:
            n_unmatched += 1
            for col in NEW_COLS:
                r[col] = ""
            continue
        r["A"] = fmt_num(m["A"])
        r["D"] = fmt_num(m["D"])
        r["F"] = fmt_num(int(m["F"]))
        r["S"] = fmt_num(m["S"])
        r["P"] = fmt_num(m.get("P"))
        r["mask_rows"] = fmt_num(int(m["shape"][0]))
        r["mask_cols"] = fmt_num(int(m["shape"][1]))
        r["bbox_area"] = fmt_num(int(m["bbox_area"]))
        r["centroid_r"] = fmt_num(m["centroid_row"])
        r["centroid_c"] = fmt_num(m["centroid_col"])
        r["peak_r"] = fmt_num(int(m["peak_row"]))
        r["peak_c"] = fmt_num(int(m["peak_col"]))
        n_filled += 1

    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames + NEW_COLS)
        writer.writeheader()
        writer.writerows(src_rows)

    n_by_state = {lab: sum(1 for r in rows if r["damage_label"] == lab)
                  for lab in ("healthy", "DAM 3", "DAM 4", "DAM 6")}
    meta = {
        "generator": "code/lumo/lumo_ocv_channels_csv.py",
        "source_csv": os.path.relpath(args.csv, DATA),
        "source_csv_sha256": sha256(args.csv),
        "out_csv": os.path.relpath(args.out, DATA),
        "burst_dir": args.burst_dir,
        "mask_definition": {"peak_frac": masks.PEAK_FRAC,
                            "median_mult": masks.MEDIAN_MULT,
                            "min_n_masked": masks.MIN_N_MASKED},
        "dimensions": {
            "A": "n_masked (identical to coherence_masked_pixels)",
            "D": "n_masked / bbox_area",
            "F": "number of 8-neighbourhood components of the mask",
            "S": "distance mask centroid <-> peak pixel (px)",
            "P": ("mean pixel frequency inside the own mask; "
                  "frequency map from the dates of the majority chip shape"),
        },
        "persistence": persistence_meta,
        "n_cache_dirs_matched": len(raw),
        "n_cache_skipped": n_skipped,
        "n_overpasses": len(rows),
        "n_by_state": n_by_state,
        "n_rows_filled_by_burst_id": n_by_bid,
        "n_rows_filled_by_date_orbit": n_by_key,
        "n_rows_unmatched": n_unmatched,
        "cross_check_n_masked_vs_committed": {"n_checked": check["n_checked"],
                                              "n_mismatch": check["n_mismatch"]},
        "committed_reference": "data/lumo/reference/lumo_tower_coherence_states.json",
        "created_rows": n_filled,
    }
    with open(args.meta, "w") as fh:
        json.dump(meta, fh, indent=1, sort_keys=True)
    return meta, rows


def _f(x):
    x = (x or "").strip()
    return float(x) if x else None


def verify(args, cache_dir=None):
    """Hard self-check of the extended CSV (ok=False => abort)."""
    problems, warnings = [], []
    coh = load_committed()
    with open(args.out, newline="") as fh:
        reader = csv.DictReader(fh)
        fields = list(reader.fieldnames)
        all_rows = list(reader)
    rows = [r for r in all_rows if (r.get("coherence_gamma2") or "").strip()]

    for col in NEW_COLS:
        if col not in fields:
            problems.append(f"missing column {col}")
    if len(rows) != len(coh["bursts"]):
        problems.append(f"{len(rows)} coherence rows instead of {len(coh['bursts'])}")
    if [r["burst_id"] for r in rows] != [b["burst_id"] for b in coh["bursts"]]:
        problems.append("order of the coherence rows != committed COH-JSON")
    n_filled_elsewhere = sum(
        1 for r in all_rows
        if not (r.get("coherence_gamma2") or "").strip()
        and any((r.get(c) or "").strip() for c in NEW_COLS))
    if n_filled_elsewhere:
        problems.append(f"{n_filled_elsewhere} non-coherence rows carry EMV values")

    n_p, n_odd, n_bad_dim, n_geom = 0, 0, 0, 0
    max_peak_dev, max_gamma_dev, n_peak_dev = 0.0, 0.0, 0
    for r, b in zip(rows, coh["bursts"]):
        A = _f(r.get("coherence_masked_pixels"))
        A_new = _f(r.get("A"))
        D, F, S, P = _f(r.get("D")), _f(r.get("F")), _f(r.get("S")), _f(r.get("P"))
        bbox = _f(r.get("bbox_area"))
        mr, mc = _f(r.get("mask_rows")), _f(r.get("mask_cols"))
        if A_new != A:
            n_bad_dim += 1
        if A != float(b["n_masked"]):
            problems.append(f"A != committed n_masked at {b['burst_id']}")
        g = _f(r.get("coherence_gamma2"))
        max_gamma_dev = max(max_gamma_dev, abs((g or 0.0) - b["gamma2"]))
        pk = _f(r.get("peak_intensity"))
        if pk is not None and pk != b["peak_intensity"]:
            n_peak_dev += 1
            max_peak_dev = max(max_peak_dev, abs(pk - b["peak_intensity"]))
        if bbox is None or A is None or D is None or abs(D - A / bbox) > 1e-12:
            n_geom += 1
        if F is None or F < 1.0 or S is None or S < 0.0:
            n_bad_dim += 1
        if (mr, mc) != (400.0, 11.0):
            n_odd += 1
            if P is not None:
                problems.append(f"P set for a non-majority shape at {b['burst_id']}")
        elif P is None or not (0.0 <= P <= 1.0):
            problems.append(f"P missing/outside [0,1] at {b['burst_id']}")
        else:
            n_p += 1
    if n_bad_dim:
        problems.append(f"{n_bad_dim} rows with implausible A/D/F/S values")
    if n_geom:
        problems.append(f"{n_geom} rows without D == A/bbox_area")
    if n_p != 177 or n_odd != 1:
        problems.append(f"P distribution {n_p} set / {n_odd} excluded "
                        f"(expected 177 / 1)")
    if max_gamma_dev > 1e-12:
        problems.append(f"gamma2 deviates from the committed JSON (max {max_gamma_dev:g})")
    if n_peak_dev:
        warnings.append(f"{n_peak_dev} peak_intensity rounding differences "
                        f"(max {max_peak_dev:g}) vs. committed JSON")

    n_cache_checked = 0
    if cache_dir and os.path.isdir(cache_dir):
        keys = committed_keys(coh)
        raw, _ = walk_cache(cache_dir, keys)
        rec, _ = committed_rows(raw, coh)
        masks.add_persistence(rec)
        by_bid = {x["committed_burst_id"]: x for x in rec}
        for r, b in zip(rows, coh["bursts"]):
            m = by_bid.get(b["burst_id"])
            if m is None:
                problems.append(f"cache recomputation missing for {b['burst_id']}")
                continue
            n_cache_checked += 1
            for col in DIMS:
                if _f(r.get(col)) != m.get(col):
                    problems.append(f"cache recomputation {col} != CSV at {b['burst_id']}")
                    break
    return {
        "ok": not problems,
        "n_rows": len(all_rows), "n_coherence_rows": len(rows),
        "n_P_set": n_p, "n_P_excluded_other_shape": n_odd,
        "max_gamma2_deviation": max_gamma_dev,
        "n_peak_intensity_roundtrip_diff": n_peak_dev,
        "n_cache_rows_reexamined": n_cache_checked,
        "problems": problems, "warnings": warnings,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", default=DEFAULT_CSV, help="source CSV (copy from espoo)")
    ap.add_argument("--out", default=DEFAULT_OUT, help="target CSV (extended)")
    ap.add_argument("--meta", default=DEFAULT_META, help="provenance metadata (JSON)")
    ap.add_argument("--burst-dir", default=DEFAULT_BURST_DIR,
                    help="burst cache monthly_bursts/ (only for build/recompute)")
    ap.add_argument("--verify-only", action="store_true",
                    help="write nothing, only check the existing target CSV")
    ap.add_argument("--no-cache-verify", action="store_true",
                    help="skip the cache recomputation of the masks")
    args = ap.parse_args()

    if not args.verify_only:
        if not os.path.isdir(args.burst_dir):
            raise SystemExit(
                f"ERROR: burst cache not found: {args.burst_dir}\n"
                f"        (needed only once; --verify-only checks without cache)")
        meta, _rows = build(args)
        print(f"written: {os.path.relpath(args.out)} "
              f"({meta['created_rows']} rows extended, "
              f"{meta['n_rows_unmatched']} without cache match)")
        print(f"written: {os.path.relpath(args.meta)}")

    rep = verify(args, cache_dir=(None if args.no_cache_verify else args.burst_dir))
    print(json.dumps(rep, indent=1, sort_keys=True))
    if not rep["ok"]:
        raise SystemExit("VERIFICATION FAILED")
    print("VERIFICATION OK")


if __name__ == "__main__":
    main()


