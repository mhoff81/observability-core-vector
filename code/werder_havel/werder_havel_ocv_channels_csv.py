#!/usr/bin/env python3
"""WERDER OCV — the committed mask cache and the support channel table.

Why this generator exists
-------------------------
``data/werder_havel/werder_havel_measurements_full.txt`` is the export of
``onboarder.measurements``: 28 columns, one row per ``(track, date)``, and **no
mask dimension**. The mask layer of this site is the project-wide echo mask of
the 80x80 complex window of every acquisition, and that payload is the one
*external* input of the package (~250 MB of complex JSON, ``WINDOW_NOTE`` in
``werder_havel_ocv_masks.py``). Exactly as Carola does with
``carola_ocv_masks.py --extract`` + ``carola_ocv_channels_csv.py``, this file
splits the work in two:

  * ``--extract`` streams the payload once and writes the committed **window
    mask cache** (``support/werder_havel_windows_mask_cache.txt``): the echo
    mask of every chip as pixel coordinates plus the sufficient statistics the
    six dimensions need. The manifest beside it records the source file's
    sha256. Nothing downstream reads the 250 MB file again.
  * the default mode builds ``support/werder_havel_ocv_channels.csv`` from that
    cache and the committed measurement extract, and ``--check`` re-verifies
    the committed table offline.

The six dimensions are the registered ones, in the family's spelling
(``carola_ocv_masks.vector_from_row`` / ``ocv_reconfigure``):

    A = |mask|                                   (echo gate: A >= 2)
    D = A / bbox_area                            (mask occupancy)
    F = 8-connected components of the mask       (fragmentation, count)
    S = ||centroid - brightest pixel||  [px]     (mask shift; identity scale)
    P = mean pixel frequency inside the own mask (echo persistence)
    gamma2 = |sum z|^2 / (A * sum |z|^2)         (masked coherence consistency)

``P`` needs a population of masks of the same shape, so it is computed **within
a patch**: the export's ``target`` rows are the bridge's own chip, while
``reference`` and ``reference_outside`` are two *different* patches (a
reference object elsewhere). Their masks are never mixed into the bridge's
frequency map -- a frequency map over three objects would be a property of the
export, not of this asset. The trajectory plate reads the ``target`` rows only.

What this table is not
----------------------
The export's own ``coherence_masked_pixels`` (640 in every row, i.e. 10 % of
6400) and ``coherence_gamma2`` are carried as provenance. They are **not** the
mask of these chips: 640 is a definition (the brightest decile of the window),
not a count of the echo mask, and the two columns are never conflated
(``A``/``gamma2_mask`` are recomputed from the payload). The meta file reports
how often the two happen to coincide.

Usage
-----
  python3 code/werder_havel/werder_havel_ocv_channels_csv.py            # build + check
  python3 code/werder_havel/werder_havel_ocv_channels_csv.py --check    # check only
  python3 code/werder_havel/werder_havel_ocv_channels_csv.py --extract \\
      --windows /home/projects/havel_werder_analysis/corrected2/werder_windows_full.txt
  python3 code/werder_havel/werder_havel_ocv_channels_csv.py --json     # census

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
sys.path.insert(0, HERE)
import werder_havel_ocv_masks as masks  # noqa: E402

DATA = masks.DATA
SUPPORT = os.path.join(DATA, "support")
CACHE_PATH = os.path.join(SUPPORT, "werder_havel_windows_mask_cache.txt")
MANIFEST_PATH = CACHE_PATH + ".manifest.json"
OUT_PATH = os.path.join(SUPPORT, "werder_havel_ocv_channels.csv")
META_PATH = os.path.join(SUPPORT, "werder_havel_ocv_channels_meta.json")

SCHEMA = "werder-ocv-channels-v1"
WINDOW_NOTE = masks.WINDOW_NOTE

#: The cache columns. ``px`` is ``row,col;row,col;...`` -- exactly the encoding
#: ``carola_ocv_masks._px_encode`` uses, so the two caches stay comparable.
CACHE_COLUMNS = ["mid", "date", "role", "track", "w", "h", "peak", "median",
                 "max_masked", "min_masked", "ss_re", "ss_im", "ss_abs2",
                 "n_masked", "peak_row", "peak_col", "n_components",
                 "largest_component", "px"]

#: The table's columns, in file order. Dimensions first, then provenance, then
#: the record's own co-variates (panel c1 of the trajectory plate).
CSV_COLUMNS = [
    "id", "date", "month", "year", "role", "track", "state_window",
    "orbit_direction", "relative_orbit", "subswath", "platform",
    "acquisition_ts",
    "A", "D", "F", "S", "P", "gamma2_mask", "echo", "patch",
    "temperature_c", "precipitation_mm", "shift_x", "shift_y",
    "coherence_gamma2", "coherence_masked_pixels", "phase_scatterer_count",
    "phase_coherence", "displacement_los_m",
]

#: Columns that must never be empty on any row (the six dimensions and the
#: identifiers). An empty co-variate is a missing value and stays empty.
REQUIRED = ["id", "date", "month", "year", "role", "track", "state_window",
            "orbit_direction", "relative_orbit", "subswath", "platform",
            "A", "D", "F", "S", "P", "gamma2_mask", "echo", "patch"]


# ---------------------------------------------------------------------------
# Formatting
# ---------------------------------------------------------------------------
def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fmt(value):
    """Write a value into the pipe-delimited table (empty stays empty)."""
    if value is None:
        return ""
    if isinstance(value, bool):
        return "1" if value else "0"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return "%.17g" % value
    return str(value)


def _px_decode(text):
    if not text:
        return []
    out = []
    for part in text.split(";"):
        r, c = part.split(",")
        out.append((int(r), int(c)))
    return out


def _px_encode(px):
    return ";".join("%d,%d" % (r, c) for r, c in px)


# ---------------------------------------------------------------------------
# The mask cache
# ---------------------------------------------------------------------------
def extract(windows_path=None, out_path=CACHE_PATH,
            manifest_path=MANIFEST_PATH, meas_path=masks.MEAS_PATH):
    """Stream the window payload once, write the committed mask cache.

    Every mask is checked against the rule before it is written (the same
    re-verification ``vector_from_row`` performs on read), and the manifest
    records the source path, its sha256, the shape census and the guard counts.
    """
    windows_path = windows_path or masks.WINDOWS_PATH
    if not os.path.exists(windows_path):
        raise SystemExit("ERROR: no window payload at %s -- pass --windows"
                         % windows_path)
    meas = masks.read_measurements(meas_path)
    rows, shapes, n_skipped, n_unknown = [], Counter(), 0, 0
    for mid, z in masks.read_windows(windows_path):
        rec = meas.get(mid)
        if rec is None:
            n_unknown += 1
            continue
        stats = masks.echo_mask(z)
        if not stats["has_echo"]:
            n_skipped += 1
        h, w = z.shape
        shapes["%dx%d" % (w, h)] += 1
        mask = stats["mask"]
        px = list(zip(np.nonzero(mask)[0].tolist(),
                      np.nonzero(mask)[1].tolist()))
        rows.append({
            "mid": mid, "date": rec["date"], "role": rec["role"],
            "track": rec["track"], "w": w, "h": h,
            "peak": stats["peak"], "median": stats["median"],
            "max_masked": stats["max_masked"], "min_masked": stats["min_masked"],
            "ss_re": stats["ss_re"], "ss_im": stats["ss_im"],
            "ss_abs2": stats["ss_abs2"], "n_masked": stats["n_masked"],
            "peak_row": stats["peak_row"], "peak_col": stats["peak_col"],
            "n_components": stats["n_components"],
            "largest_component": stats["largest_component"],
            "px": px,
        })
    rows.sort(key=lambda r: (r["date"], r["track"], r["mid"]))
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write("# werder windows mask cache v1 — the echo mask of every "
                 "committed chip\n")
        fh.write("# rule: intensity >= %.2f * max(intensity) AND intensity >= "
                 "%.1f * np.median(intensity); min_n_masked=%d\n"
                 % (masks.PEAK_FRAC, masks.MEDIAN_MULT, masks.MIN_N_MASKED))
        fh.write("# extractor: code/werder_havel/werder_havel_ocv_channels_csv.py"
                 " --extract\n")
        fh.write("|".join(CACHE_COLUMNS) + "\n")
        for r in rows:
            fh.write("|".join([
                r["mid"], r["date"], r["role"], r["track"], str(r["w"]),
                str(r["h"]), repr(r["peak"]), repr(r["median"]),
                repr(r["max_masked"]), repr(r["min_masked"]), repr(r["ss_re"]),
                repr(r["ss_im"]), repr(r["ss_abs2"]), str(r["n_masked"]),
                str(r["peak_row"]), str(r["peak_col"]),
                str(r["n_components"]), str(r["largest_component"]),
                _px_encode(r["px"]),
            ]) + "\n")
    manifest = {
        "schema": SCHEMA,
        "generator": "code/werder_havel/werder_havel_ocv_channels_csv.py "
                     "--extract",
        "out_cache": os.path.basename(out_path),
        "out_cache_sha256": _sha256(out_path),
        "source_windows": os.path.abspath(windows_path),
        "source_windows_sha256": _sha256(windows_path),
        "source_measurements": os.path.relpath(meas_path,
                                               os.path.dirname(DATA)),
        "source_measurements_sha256": _sha256(meas_path),
        "rows": len(rows),
        "skipped_no_echo": n_skipped,
        "skipped_unknown_id": n_unknown,
        "shapes": dict(sorted(shapes.items())),
        "rule": {
            "peak_frac": masks.PEAK_FRAC, "median_mult": masks.MEDIAN_MULT,
            "min_n_masked": masks.MIN_N_MASKED,
            "text": "intensity >= 0.30 * max(intensity) AND intensity >= "
                    "5.0 * np.median(intensity); A = count(mask)",
        },
        "window_note": WINDOW_NOTE,
        "roles": dict(sorted(Counter(r["role"] for r in rows).items())),
        "min_A": min(r["n_masked"] for r in rows),
        "max_A": max(r["n_masked"] for r in rows),
    }
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=1, sort_keys=True)
        fh.write("\n")
    print("wrote %s (%d rows) and %s"
          % (os.path.relpath(out_path), len(rows),
             os.path.relpath(manifest_path)))
    return manifest


def load_cache(path=CACHE_PATH):
    """The committed cache as rows with ``mask_rows``/``mask_cols`` arrays."""
    rows = []
    header = None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("|")
            if header is None:
                header = parts
                continue
            d = dict(zip(header, parts))
            d["px"] = _px_decode(d["px"])
            d["w"], d["h"] = int(d["w"]), int(d["h"])
            d["n_masked"] = int(d["n_masked"])
            d["peak_row"] = int(d["peak_row"])
            d["peak_col"] = int(d["peak_col"])
            d["n_components"] = int(d["n_components"])
            d["largest_component"] = int(d["largest_component"])
            for key in ("peak", "median", "max_masked", "min_masked",
                        "ss_re", "ss_im", "ss_abs2"):
                d[key] = float(d[key])
            d["mask_rows"] = np.array([r for r, _ in d["px"]], dtype=int)
            d["mask_cols"] = np.array([c for _, c in d["px"]], dtype=int)
            rows.append(d)
    return rows


def vector_from_row(row):
    """The six mask dimensions of one cache row, rule re-verified.

    Verbatim the arithmetic of ``carola_ocv_masks.vector_from_row``: the mask
    rule is re-checked on the committed pixels (threshold on the faintest masked
    pixel, the brightest masked pixel *is* the peak, ``A >= 2``), then
    ``A``/``D``/``F``/``S``/``gamma2`` are derived from the committed set.
    """
    px = row["px"]
    A = len(px)
    if A != row["n_masked"]:
        raise ValueError("%s: n_masked %d != |px| %d"
                         % (row["mid"], row["n_masked"], A))
    if A < masks.MIN_N_MASKED:
        raise ValueError("%s: A=%d < %d"
                         % (row["mid"], A, masks.MIN_N_MASKED))
    peak, median = row["peak"], row["median"]
    lo_thr = max(masks.PEAK_FRAC * peak, masks.MEDIAN_MULT * median)
    if row["min_masked"] * (1.0 + masks.RULE_TOL) < lo_thr:
        raise ValueError("%s: committed mask violates the rule "
                         "(min|z|^2=%r < %r)"
                         % (row["mid"], row["min_masked"], lo_thr))
    if abs(row["max_masked"] - peak) > masks.RULE_TOL * max(peak, 1e-30):
        raise ValueError("%s: peak %r is not a masked pixel (max_masked=%r)"
                         % (row["mid"], peak, row["max_masked"]))
    rows = np.array([r for r, _ in px], dtype=float)
    cols = np.array([c for _, c in px], dtype=float)
    row_span = int(rows.max() - rows.min() + 1)
    col_span = int(cols.max() - cols.min() + 1)
    bbox_area = int(row_span * col_span)
    cen_r, cen_c = float(rows.mean()), float(cols.mean())
    return {
        "A": A,
        "D": (A / bbox_area) if bbox_area > 0 else None,
        "F": row["n_components"],
        "S": float(np.hypot(cen_r - row["peak_row"], cen_c - row["peak_col"])),
        "gamma2": masks.gamma2_of(row["ss_re"], row["ss_im"],
                                  row["ss_abs2"], A),
        "bbox_area": bbox_area,
        "row_span": row_span,
        "col_span": col_span,
        "centroid_row": cen_r,
        "centroid_col": cen_c,
        "peak_row": row["peak_row"],
        "peak_col": row["peak_col"],
        "n_components": row["n_components"],
        "largest_component": row["largest_component"],
        "fragmentation": (1.0 - row["largest_component"] / float(max(A, 1))),
    }


def add_persistence(rows):
    """``P`` per patch: mean pixel frequency inside the own mask.

    A pixel-exact frequency map is only comparable across chips of the same
    shape. Every Werder chip is 80x80, so the shape filter is a no-op here and
    is asserted rather than assumed. The map is computed **inside a role**
    (``target`` / ``reference`` / ``reference_outside``): the export's reference
    chains measure different patches, and a frequency map over three objects
    would describe the export, not the bridge.
    """
    shapes = Counter((r["w"], r["h"]) for r in rows)
    if len(shapes) != 1:
        raise SystemExit("ERROR: the cache holds %d chip shapes: %s"
                         % (len(shapes), dict(shapes)))
    h, w = next(iter(shapes))
    out = {}
    for role in sorted({r["role"] for r in rows}):
        group = [r for r in rows if r["role"] == role]
        freq = np.zeros((h, w), dtype=float)
        for r in group:
            freq[r["mask_rows"], r["mask_cols"]] += 1.0
        freq /= float(len(group))
        for r in group:
            values = freq[r["mask_rows"], r["mask_cols"]]
            r["P"] = float(values.mean())
            r["P_frac_ge_50pct"] = float((values >= 0.5).mean())
            r["P_frac_ge_30pct"] = float((values >= 0.3).mean())
        out[role] = {"n_chips": len(group),
                     "median_P": float(np.median([r["P"] for r in group])),
                     "min_P": float(min(r["P"] for r in group)),
                     "max_P": float(max(r["P"] for r in group))}
    return {"shape": [w, h], "patches": out,
            "note": "P is computed within a role (a patch), never across the "
                    "export's three objects"}


# ---------------------------------------------------------------------------
# The channel table
# ---------------------------------------------------------------------------
def _number(text):
    text = (text or "").strip()
    return float(text) if text else None


def build(cache_path=CACHE_PATH, out_path=OUT_PATH, meta_path=META_PATH,
          meas_path=masks.MEAS_PATH):
    """Recompute the six dimensions from the cache and write the table + meta."""
    cache = load_cache(cache_path)
    if not cache:
        raise SystemExit("ERROR: the mask cache is empty -- run --extract first")
    persistence = add_persistence(cache)
    meas = masks.read_measurements(meas_path)
    rows, n_export_gamma2, export_pixels = [], 0, Counter()
    for rec in cache:
        source = meas[rec["mid"]]
        vec = vector_from_row(rec)
        export_gamma2 = _number(source.get("coherence_gamma2"))
        export_px = _number(source.get("coherence_masked_pixels"))
        if export_px is not None:
            export_pixels[int(export_px)] += 1
        if export_gamma2 is not None:
            n_export_gamma2 += 1
        rows.append({
            "id": rec["mid"],
            "date": rec["date"],
            "month": rec["date"][:7],
            "year": int(rec["date"][:4]),
            "role": rec["role"],
            "track": rec["track"],
            "state_window": source.get("state_window", "").strip(),
            "orbit_direction": source.get("orbit_direction", "").strip().upper(),
            "relative_orbit": source.get("relative_orbit", "").strip(),
            "subswath": source.get("subswath", "").strip(),
            "platform": source.get("platform", "").strip(),
            "acquisition_ts": source.get("acquisition_ts", "").strip(),
            "A": vec["A"], "D": vec["D"], "F": vec["F"], "S": vec["S"],
            "P": rec["P"], "gamma2_mask": vec["gamma2"],
            "echo": vec["A"] >= masks.MIN_N_MASKED,
            "patch": rec["role"],
            "temperature_c": _number(source.get("temperature_c")),
            "precipitation_mm": _number(source.get("precipitation_mm")),
            "shift_x": _number(source.get("shift_x")),
            "shift_y": _number(source.get("shift_y")),
            "coherence_gamma2": export_gamma2,
            "coherence_masked_pixels": export_px,
            "phase_scatterer_count": _number(source.get("phase_scatterer_count")),
            "phase_coherence": _number(source.get("phase_coherence")),
            "displacement_los_m": _number(source.get("displacement_los_m")),
        })
    rows.sort(key=lambda r: (r["date"], r["track"]))
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write("|".join(CSV_COLUMNS) + "\n")
        for row in rows:
            fh.write("|".join(fmt(row.get(column))
                              for column in CSV_COLUMNS) + "\n")
    meta = build_meta(rows, persistence, cache_path, out_path, meas_path,
                      n_export_gamma2, export_pixels)
    with open(meta_path, "w") as fh:
        json.dump(meta, fh, indent=1, sort_keys=True)
        fh.write("\n")
    print("wrote %s (%d rows x %d columns) and %s"
          % (os.path.relpath(out_path), len(rows), len(CSV_COLUMNS),
             os.path.relpath(meta_path)))
    return meta, rows


def build_meta(rows, persistence, cache_path, out_path, meas_path,
               n_export_gamma2, export_pixels):
    """The meta file: provenance, the rule, the census and the guard counts."""
    target = [r for r in rows if r["role"] == "target"]
    return {
        "schema": SCHEMA,
        "generator": "code/werder_havel/werder_havel_ocv_channels_csv.py",
        "asset": masks.ASSET, "asset_id": masks.ASSET_ID,
        "asset_name": masks.ASSET_NAME, "asset_type": masks.ASSET_TYPE,
        "table": os.path.relpath(out_path, os.path.dirname(DATA)),
        "table_sha256": _sha256(out_path),
        "mask_cache": os.path.relpath(cache_path, os.path.dirname(DATA)),
        "mask_cache_sha256": _sha256(cache_path),
        "measurements": os.path.relpath(meas_path, os.path.dirname(DATA)),
        "measurements_sha256": _sha256(meas_path),
        "rows": len(rows), "columns": len(CSV_COLUMNS),
        "column_names": CSV_COLUMNS,
        "roles": dict(sorted(Counter(r["role"] for r in rows).items())),
        "target": {
            "rows": len(target),
            "n_dates": len({r["date"] for r in target}),
            "n_tracks": len({r["track"] for r in target}),
            "tracks": dict(sorted(Counter(r["track"] for r in target).items())),
            "orbit_directions": dict(sorted(
                Counter(r["orbit_direction"] for r in target).items())),
            "first_date": target[0]["date"], "last_date": target[-1]["date"],
            "echo": sum(1 for r in target if r["echo"]),
        },
        "mask_rule": {
            "peak_frac": masks.PEAK_FRAC, "median_mult": masks.MEDIAN_MULT,
            "min_n_masked": masks.MIN_N_MASKED,
            "text": "intensity >= 0.30 * max(intensity) AND intensity >= "
                    "5.0 * np.median(intensity); A = count(mask); "
                    "gamma2 = |sum z|^2 / (A * sum |z|^2), clamped to [0, 1]",
            "window": "80 x 80 complex samples (%d px)" % masks.WINDOW_PX,
            "window_note": WINDOW_NOTE,
        },
        "dimensions": {
            "A": "count of the echo mask (px)",
            "D": "A / bbox_area of the mask",
            "F": "8-connected components of the mask",
            "S": "distance from the mask centroid to its brightest pixel (px)",
            "P": "mean pixel frequency inside the own mask, per patch",
            "gamma2_mask": "masked coherence consistency, clamped to [0, 1]",
        },
        "persistence": persistence,
        "export_columns": {
            "coherence_gamma2": "the export's own coherence column — provenance, "
                                "never the mask rule's gamma2",
            "coherence_masked_pixels": "the export's own definition (the "
                                       "brightest decile of the window), not a "
                                       "count of the echo mask",
            "masked_pixels_values": dict(sorted(export_pixels.items())),
            "rows_with_export_gamma2": n_export_gamma2,
            "n_A_equal_export_masked_pixels": sum(
                1 for r in rows if r["coherence_masked_pixels"] == r["A"]),
        },
        "event": {"start": masks.EVENT_START, "end": masks.EVENT_END,
                  "source": masks.EVENT_SOURCE},
        "state_windows": {w: sum(1 for r in rows if r["state_window"] == w)
                          for w in sorted({r["state_window"] for r in rows})},
        "notes": [
            "one row per (track, date) of the committed measurement extract; "
            "the mask dimensions are recomputed from the committed mask cache",
            "P is computed within a patch (role): the reference chains measure "
            "different objects and are never folded into the bridge's frequency "
            "map",
            "echo = 1 wherever A >= 2; a chip without an echo keeps its "
            "geometry dimensions but carries no behaviour channel",
        ],
    }


# ---------------------------------------------------------------------------
# The acceptance test
# ---------------------------------------------------------------------------
def _read_table(path=OUT_PATH):
    with open(path) as fh:
        keys = fh.readline().rstrip("\n").split("|")
        rows = [dict(zip(keys, line.rstrip("\n").split("|")))
                for line in fh if line.strip()]
    return keys, rows


def check(cache_path=CACHE_PATH, out_path=OUT_PATH, meta_path=META_PATH,
          meas_path=masks.MEAS_PATH):
    """Re-verify the committed table and cache offline.

    Returns ``(problems, facts)``: nothing here trusts the meta file -- every
    dimension is recomputed from the committed pixels, the rule is re-checked
    per chip, and the sha256 chain (payload -> cache -> table -> meta) is
    followed one link at a time.
    """
    problems = []
    for path in (cache_path, out_path, meta_path, MANIFEST_PATH):
        if not os.path.exists(path):
            raise SystemExit("ERROR: missing %s" % path)
    with open(meta_path) as fh:
        meta = json.load(fh)
    with open(MANIFEST_PATH) as fh:
        manifest = json.load(fh)
    if _sha256(out_path) != meta["table_sha256"]:
        problems.append("the table's sha256 does not match the meta file")
    if _sha256(cache_path) != meta["mask_cache_sha256"]:
        problems.append("the cache's sha256 does not match the meta file")
    if manifest["out_cache_sha256"] != meta["mask_cache_sha256"]:
        problems.append("the manifest and the meta disagree on the cache")

    cache = load_cache(cache_path)
    persistence = add_persistence(cache)
    keys, rows = _read_table(out_path)
    if keys != CSV_COLUMNS:
        problems.append("columns %r != %r" % (keys, CSV_COLUMNS))
    if len(rows) != len(cache):
        problems.append("%d table rows != %d cache rows"
                        % (len(rows), len(cache)))
    if len(rows) != meta["rows"]:
        problems.append("%d table rows != %d in the meta"
                        % (len(rows), meta["rows"]))

    by_id = {r["mid"]: r for r in cache}
    for row in rows:
        rec = by_id.get(row["id"])
        if rec is None:
            problems.append("%s: not in the cache" % row["id"])
            continue
        vec = vector_from_row(rec)
        if row["date"] != rec["date"]:
            problems.append("%s: date %r != %r"
                            % (row["id"], row["date"], rec["date"]))
        if row["role"] != rec["role"] or row["track"] != rec["track"]:
            problems.append("%s: role/track differ from the cache" % row["id"])
        want = {"A": vec["A"], "F": vec["F"], "D": vec["D"], "S": vec["S"],
                "gamma2_mask": vec["gamma2"], "P": rec["P"]}
        for dim, target in want.items():
            got = float(row[dim])
            if abs(got - target) > 1e-12 * max(1.0, abs(target)):
                problems.append("%s: %s %r != %r"
                                % (row["id"], dim, got, target))
        if float(row["D"]) <= 0.0 or float(row["D"]) > 1.0:
            problems.append("%s: D=%r outside (0, 1]" % (row["id"], row["D"]))
        if row["echo"] != ("1" if vec["A"] >= masks.MIN_N_MASKED else "0"):
            problems.append("%s: echo flag disagrees with A" % row["id"])
        if row["state_window"] not in ("W2_bestand", "W3_sanierung",
                                       "W4_post"):
            problems.append("%s: state_window %r"
                            % (row["id"], row["state_window"]))
        for column in REQUIRED:
            if row.get(column, "") == "":
                problems.append("%s: empty required column %s"
                                % (row["id"], column))
                break

    target = [r for r in rows if r["role"] == "target"]
    if len(target) != masks.RECORD_EXPECTED["roles"]["target"]:
        problems.append("target rows %d != %d"
                        % (len(target),
                           masks.RECORD_EXPECTED["roles"]["target"]))
    if len({r["date"] for r in target}) != len(target):
        problems.append("the target chain is not one row per date")
    if any(float(r["gamma2_mask"]) <= 0.0 for r in target):
        problems.append("a target chip carries gamma2_mask <= 0 "
                        "(log undefined)")
    if meta["persistence"] != persistence:
        problems.append("the meta's persistence block does not recompute")
    facts = {"rows": len(rows), "target": len(target),
             "cache_rows": len(cache), "first": target[0]["date"],
             "last": target[-1]["date"],
             "orbit": dict(sorted(Counter(r["orbit_direction"]
                                          for r in target).items())),
             "meta": meta}
    return problems, facts


def census(**kwargs):
    """The table's census: what the committed file holds, checked."""
    problems, facts = check(**kwargs)
    meta = facts["meta"]
    return {
        "table": meta["table"], "table_sha256": meta["table_sha256"],
        "rows": facts["rows"], "columns": meta["columns"],
        "roles": meta["roles"], "target": meta["target"],
        "state_windows": meta["state_windows"],
        "persistence": meta["persistence"],
        "export": meta["export_columns"],
        "problems": problems,
    }


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--extract", action="store_true",
                    help="rebuild the committed mask cache from the window "
                         "payload (the one external input)")
    ap.add_argument("--windows", default=None,
                    help="the window payload to stream with --extract "
                         "(default: %s)" % masks.WINDOWS_PATH)
    ap.add_argument("--check", action="store_true",
                    help="re-verify the committed table (the default mode)")
    ap.add_argument("--json", action="store_true",
                    help="print the census of the committed table and exit")
    ap.add_argument("--cache", default=CACHE_PATH)
    ap.add_argument("--out", default=OUT_PATH)
    ap.add_argument("--meta", default=META_PATH)
    ap.add_argument("--meas", default=masks.MEAS_PATH)
    args = ap.parse_args(argv)
    if args.extract:
        extract(args.windows, args.cache, args.cache + ".manifest.json",
                args.meas)
        build(args.cache, args.out, args.meta, args.meas)
        return 0
    if args.json:
        print(json.dumps(census(cache_path=args.cache, out_path=args.out,
                                meta_path=args.meta, meas_path=args.meas),
                         indent=1, sort_keys=True))
        return 0
    if not args.check:
        build(args.cache, args.out, args.meta, args.meas)
    problems, facts = check(cache_path=args.cache, out_path=args.out,
                            meta_path=args.meta, meas_path=args.meas)
    if problems:
        print("deviations:", file=sys.stderr)
        for problem in problems:
            print("  - %s" % problem, file=sys.stderr)
        raise SystemExit(1)
    meta = facts["meta"]
    print("ok: %d rows (%d target, %s) %s .. %s; %d columns; echo %d/%d; "
          "patches %s"
          % (facts["rows"], facts["target"], facts["orbit"], facts["first"],
             facts["last"], meta["columns"], meta["target"]["echo"],
             facts["target"],
             {k: v["n_chips"] for k, v in meta["persistence"]["patches"].items()}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
