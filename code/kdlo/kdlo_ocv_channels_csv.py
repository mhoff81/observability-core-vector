#!/usr/bin/env python3
"""KDLO OCV-Paper — build ``data/kdlo/kdlo_ocv_channels.csv`` from committed files.

Why this generator exists
-------------------------
The KDLO record is two committed JSON files (``fig_kdlo_1_raw.json``,
``fig_kdlo_2_raw.json``) holding the 60 API measurement rows of the two figure
windows, plus four committed reference JSONs. Nine of the channels the KDLO
project reported are **not** in those rows, because the API never returned them:
they exist only in the pipeline's own Postgres table
``onboarder.insar_measurements``:

    phase_snr_db, phase_observable, phase_unobservable_reason, snow_depth_m,
    snowfall_cm, precipitation_mm, wind_gust_ms, intensity, incidence_angle_deg

This generator therefore does two things:

1. ``--fetch-db-extra`` (run **once**, needs the local ``tower-postgres``):
   a read-only ``SELECT`` of those nine columns for the 60 measurement ids of
   the two raw files, written to ``data/kdlo/kdlo_db_extra.json`` with the exact
   SQL and its sha256. Nothing else is taken from the database — no geometry, no
   masks, no time series, no writes.
2. the default build: merge the two raw files with the vendored extraction,
   derive ``part/state/month/year/season/echo_mode``, recompute the echo mask
   and write ``data/kdlo/kdlo_ocv_channels.csv`` + ``kdlo_ocv_channels_meta.json``.
   From then on the CSV is the only input of the figure pipeline; the database
   is never needed again, and the build is reproducible offline.
3. the mask layer. The API never returned the *mask* of an overpass, only its
   size (``coherence_masked_pixels``), so the KDLO record had no mask geometry
   at all. The masks are now recomputed from the committed full-dwell strip
   cache ``data/kdlo/strips/`` — the exact ``11 x 400`` complex strips the
   pipeline decoded, written by
   ``tower/backend/src/cli/analysis/kdlo_strip_capture.rs`` — with the rule of
   the pipeline itself (``kdlo_ocv_masks.py``). ``A`` must therefore reproduce
   ``coherence_masked_pixels`` and ``gamma2_mask`` must reproduce
   ``coherence_gamma2``; the generator refuses to write a CSV where that fails.
   The strips *are* committed, so this is a second half of the mask layer, not a
   new non-committed dependency.

Because step 1 is the *only* non-committed input, the build refuses to run
without ``kdlo_db_extra.json``, and the verification compares the nine DB
columns against the committed channel JSONs (their medians/quartiles are
committed facts): if the extraction came from a different database state than
the one the KDLO project used, that comparison fails.

Derived columns (all recomputable from the raw row, no lookups):
    part       part1/part2, from which raw file the row came
    state      epoch: pre | rebuild | single (see ``kdlo_ocv_core``)
    month      1..12 of the acquisition day
    year       YYYY of the acquisition day
    season     winter_Oct-Feb | summer_Mar-Jul | shoulder_Aug-Sep
    echo_mode  compact (<=10 masked px) | intermediate | distributed (>=25)
    A/D/F/S/P  echo-mask geometry of the strip (see ``kdlo_ocv_masks``)

Usage:
  python3 code/kdlo/kdlo_ocv_channels_csv.py --fetch-db-extra   # once, needs DB
  python3 code/kdlo/kdlo_ocv_channels_csv.py                    # build + verify
  python3 code/kdlo/kdlo_ocv_channels_csv.py --verify-only      # verify only
  python3 code/kdlo/kdlo_ocv_channels_csv.py --verify-only --no-strip-verify
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, os.pardir, os.pardir, "data", "kdlo")

sys.path.insert(0, HERE)
import kdlo_ocv_core as core  # noqa: E402
import kdlo_ocv_masks as masks  # noqa: E402

DEFAULT_OUT = core.CSV_PATH
DEFAULT_META = core.CSV_META_PATH
DEFAULT_DB_EXTRA = core.DB_EXTRA_PATH
DEFAULT_STRIPS = masks.STRIPS_DIR
DEFAULT_MANIFEST = masks.MANIFEST_PATH

DB_CONTAINER = "tower-postgres"
DB_USER = "tower"
DB_NAME = "tower"
DB_TABLE = "onboarder.insar_measurements"
TOL = 1e-12


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch_sql(ids):
    """The exact read-only SELECT used by ``--fetch-db-extra`` (one statement)."""
    cols = ["m.id", "m.request_id::text",
            "to_char(m.acquisition_ts, 'YYYY-MM-DD')"]
    cols += [f"m.{c}" for c in core.DB_EXTRA_COLUMNS]
    return ("select " + ", ".join(cols) + f" from {DB_TABLE} m where m.id in ("
            + ", ".join(f"'{i}'" for i in sorted(ids))
            + ") order by m.acquisition_ts, m.id;")


def run_sql(sql):
    out = subprocess.run(
        ["docker", "exec", DB_CONTAINER, "psql", "-U", DB_USER, "-d", DB_NAME,
         "-t", "-A", "-F", "|", "-c", sql],
        capture_output=True, text=True, timeout=120)
    if out.returncode != 0:
        raise SystemExit(f"ERROR: psql failed ({DB_CONTAINER}): {out.stderr.strip()}")
    return [ln.split("|") for ln in out.stdout.splitlines() if ln.strip()]


def _num_or_text(text):
    """One psql text cell -> JSON value (``''`` is NULL, ``t``/``f`` a bool)."""
    if text == "":
        return None
    if text in ("t", "f"):
        return text == "t"
    try:
        return float(text)
    except ValueError:
        return text


def fetch_db_extra(args):
    """The one-time, read-only extraction -> ``data/kdlo/kdlo_db_extra.json``."""
    by_part, prov = core.load_raw_parts()
    lookup = {}
    for p in core.PART_ORDER:
        for r in by_part[p]:
            lookup[str(r["id"])] = (p, str(r["acquisition_ts"])[:10])
    ids = sorted(lookup)
    sql = fetch_sql(ids)
    rows = run_sql(sql)
    keys = ["id", "request_id", "day"] + list(core.DB_EXTRA_COLUMNS)

    out, missing, moved = {}, [], []
    for vals in rows:
        if len(vals) != len(keys):
            raise SystemExit(f"ERROR: unexpected column count in DB answer: {vals!r}")
        d = dict(zip(keys, vals))
        d = {k: _num_or_text(v) for k, v in d.items()}
        if d["id"] not in lookup:
            raise SystemExit(f"ERROR: DB returned an unknown measurement id {d['id']}")
        part, day = lookup[d["id"]]
        if d["day"] != day:
            moved.append((d["id"], d["day"], day))
        d["part"] = f"part{part}"
        out[d["id"]] = d
    for i in ids:
        if i not in out:
            missing.append(i)
    if moved:
        raise SystemExit(f"ERROR: {len(moved)} measurements moved day in the DB "
                         f"(first: {moved[0]})")
    if missing:
        raise SystemExit(f"ERROR: {len(missing)} of {len(ids)} ids not in the DB "
                         f"(first: {missing[0]})")

    got = run_sql(f"select count(*) from {DB_TABLE} where asset_id = "
                  f"'{prov[1]['asset_id']}'")
    total = int(got[0][0]) if got else -1
    doc = {
        "generator": "code/kdlo/kdlo_ocv_channels_csv.py --fetch-db-extra",
        "created_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "database": f"{DB_USER}@{DB_CONTAINER}/{DB_NAME}",
        "table": DB_TABLE,
        "read_only": True,
        "access": "docker exec psql, one SELECT, no writes",
        "sql": sql,
        "sql_sha256": hashlib.sha256(sql.encode()).hexdigest(),
        "asset_id": prov[1]["asset_id"],
        "columns": list(core.DB_EXTRA_COLUMNS),
        "n_ids_requested": len(ids),
        "n_rows": len(out),
        "n_asset_rows_in_db": total,
        "n_asset_rows_outside_figure_windows": total - len(out),
        "rows": {k: out[k] for k in sorted(out)},
    }
    with open(args.db_extra, "w") as fh:
        json.dump(doc, fh, indent=1, sort_keys=True)
    print(f"written: {os.path.relpath(args.db_extra, DATA)} "
          f"({len(out)}/{len(ids)} ids, {total} rows for the asset in the DB)")
    return doc



def attach_masks(rows, strips_dir=DEFAULT_STRIPS, manifest=DEFAULT_MANIFEST):
    """Add the echo-mask columns of every row (in place) from the strip cache.

    Returns ``(meta, problems)`` — the strip-side metadata of
    ``kdlo_ocv_masks.load_vectors`` plus the acceptance findings
    (``A == coherence_masked_pixels``, ``gamma2_mask == coherence_gamma2``).
    Rows whose day has no strip keep the mask columns empty and are reported.
    """
    vec, meta = masks.load_vectors(strips_dir, manifest)
    by_day = {v["date"]: v for v in vec}
    problems, missing = [], []
    for r in rows:
        day = core.day_of(r)
        v = by_day.get(day)
        if v is None:
            missing.append(day)
            for col in core.MASK_COLUMNS:
                r[col] = None
            continue
        shape = v["shape"]
        r["A"] = v["A"]
        r["D"] = v["D"]
        r["F"] = v["F"]
        r["S"] = v["S"]
        r["P"] = v.get("P")
        r["mask_rows"], r["mask_cols"] = int(shape[0]), int(shape[1])
        r["bbox_area"] = v["bbox_area"]
        r["row_span"], r["col_span"] = v["row_span"], v["col_span"]
        r["centroid_r"], r["centroid_c"] = v["centroid_row"], v["centroid_col"]
        r["peak_r"], r["peak_c"] = v["peak_row"], v["peak_col"]
        r["peak_intensity"] = v["peak"]
        r["gamma2_mask"] = v["gamma2"]

        a_rec = core._num(r.get("coherence_masked_pixels"))
        g_rec = core._num(r.get("coherence_gamma2"))
        if a_rec is None or int(a_rec) != v["A"]:
            problems.append(f"{day}: A={v['A']} != coherence_masked_pixels={a_rec}")
        if g_rec is None or abs(v["gamma2"] - g_rec) > TOL:
            problems.append(f"{day}: gamma2_mask={v['gamma2']!r} != "
                            f"coherence_gamma2={g_rec!r}")
    if missing:
        problems.append(f"{len(missing)} CSV day(s) without a strip: {missing[:5]}")
    meta = dict(meta)
    meta["n_attached"] = len(rows) - len(missing)
    meta["days_without_strip"] = missing
    return meta, problems


def strip_cache_meta(args, strip_meta):
    """Provenance block of the mask layer for the CSV metadata."""
    return {
        "dir": os.path.relpath(args.strips_dir, DATA),
        "manifest": os.path.relpath(args.manifest, DATA),
        "manifest_sha256": sha256(args.manifest),
        "capture": ("tower/backend/src/cli/analysis/kdlo_strip_capture.rs — "
                    "CopernicusClient::download_burst_strip_ranged"),
        "decoder": "decode_slc_rect_at_location(annotation, tiff, lat, lon, 11, 400)",
        "n_strips": strip_meta["n_strips"],
        "n_with_mask": strip_meta["n_with_mask"],
        "n_without_mask": strip_meta["n_without_mask"],
        "n_attached_to_csv_rows": strip_meta["n_attached"],
        "days_without_strip": strip_meta["days_without_strip"],
        "majority_shape": strip_meta["majority_shape"],
        "n_majority_shape": strip_meta["n_majority_shape"],
        "n_excluded_other_shape": strip_meta["n_excluded_other_shape"],
        "mask_rule": {
            "intensity": "|z|^2 of the 11 x 400 full-dwell strip samples",
            "peak": "max(intensity)",
            "median": "sorted(intensity)[len // 2] — the *upper* median the Rust "
                      "pipeline uses (np.median is a diagnostic only)",
            "mask": "intensity >= 0.30 * peak and intensity >= 5.0 * median",
            "min_n_masked": masks.MIN_N_MASKED,
            "gamma2": "|sum(z[mask])|^2 / (A * sum(|z[mask]|^2)), clamped to [0,1]",
        },
        "acceptance": {
            "A == coherence_masked_pixels": "exact (integer)",
            "gamma2_mask == coherence_gamma2": f"abs(diff) <= {TOL:g}",
            "shape": f"{strip_meta['majority_shape']}",
        },
    }


def build(args):
    """Merge raw + db-extra, write the CSV and its provenance metadata."""
    if not os.path.exists(args.db_extra):
        raise SystemExit(
            f"ERROR: {os.path.relpath(args.db_extra, DATA)} not found — the nine "
            f"DB-only columns cannot be invented; run --fetch-db-extra once "
            f"against the local tower-postgres.")
    by_part, prov = core.load_raw_parts()
    extra, extra_prov = core.db_extra_rows(args.db_extra)
    rows = core.merge_rows(by_part, extra)
    if len(rows) != sum(prov[p]["n_rows"] for p in core.PART_ORDER):
        raise SystemExit("ERROR: merged row count != sum of the raw parts")

    strip_meta, strip_problems = attach_masks(rows, args.strips_dir, args.manifest)
    if strip_problems:
        raise SystemExit("ERROR: echo-mask layer incomplete: "
                         + "; ".join(strip_problems[:5]))

    core.write_csv(rows, args.out)

    meta = {
        "generator": "code/kdlo/kdlo_ocv_channels_csv.py",
        "out_csv": os.path.relpath(args.out, DATA),
        "created_rows": len(rows),
        "columns": list(core.CSV_COLUMNS),
        "derived_columns": list(core.DERIVED_COLUMNS),
        "source_columns": list(core.SOURCE_COLUMNS),
        "db_columns": list(core.DB_EXTRA_COLUMNS),
        "mask_columns": list(core.MASK_COLUMNS),
        "raw_parts": {
            str(p): {"file": prov[p]["file"], "sha256": prov[p]["sha256"],
                     "asset_id": prov[p]["asset_id"], "window": prov[p]["window"],
                     "n_rows": prov[p]["n_rows"]}
            for p in core.PART_ORDER},
        "db_extra": {"file": extra_prov["file"], "sha256": extra_prov["sha256"],
                     "source": extra_prov["source"], "columns": extra_prov["columns"],
                     "n_rows": extra_prov["n_rows"],
                     "fetched_at": extra_prov["fetched_at"]},
        "reference_files": {
            name: {"file": fname, "sha256": sha256(os.path.join(core.REF, fname))}
            for name, fname in sorted(core.REFERENCE_FILES.items())},
        "definitions": {
            "part": "part1/part2 — which committed raw file the row came from",
            "state": ("epoch from the acquisition day: pre = 2022-06-02..2022-12-11, "
                      "rebuild = 2023-05-04..2024-09-19, single = 2022-12-23 "
                      "(between the windows, in no state statistic)"),
            "month": "month number of the acquisition day",
            "year": "year of the acquisition day",
            "season": ("winter_Oct-Feb / summer_Mar-Jul / shoulder_Aug-Sep "
                       "(pipeline season_of split, verbatim)"),
            "echo_mode": ("compact A <= 10 px, distributed A >= 25 px, "
                          "intermediate in between (pipeline discriminator)"),
            "A": ("echo-mask area (n_masked) — recomputed from the committed strip "
                  "cache and equal to coherence_masked_pixels"),
            "D": "echo-mask density D = A / bbox_area",
            "F": "echo-mask fragmentation F = number of 8-connected components",
            "S": "echo-mask shift S = ||mask centroid - peak pixel||  [px]",
            "P": ("echo-mask persistence P = mean pixel frequency of the own mask "
                  "across the 64 dates of the majority strip shape"),
            "peak_intensity": "max |z|^2 of the full-dwell strip",
            "gamma2_mask": ("masked coherence recomputed from the strip; must equal "
                            "coherence_gamma2"),
        },
        "strip_cache": strip_cache_meta(args, strip_meta),
        "state_windows": {k: list(v) for k, v in core.STATE_WINDOWS.items()},
        "singleton_day": core.SINGLETON_DAY,
        "n_by_state": core.n_by_state(rows),
        "n_by_part": {key: sum(1 for r in rows if r["part"] == key)
                      for key in ("part1", "part2")},
        "n_by_echo_mode": {m: sum(1 for r in rows if r["echo_mode"] == m)
                           for m in ("compact", "intermediate", "distributed")},
        "n_states_rows": len(core.state_rows(rows)),
    }
    with open(args.meta, "w") as fh:
        json.dump(meta, fh, indent=1, sort_keys=True)
    print(f"written: {os.path.relpath(args.out, DATA)} ({len(rows)} rows, "
          f"{len(core.CSV_COLUMNS)} columns, "
          f"{strip_meta['n_attached']} echo masks from the strip cache)")
    print(f"written: {os.path.relpath(args.meta, DATA)}")
    return meta, rows



def _diff(a, b, tol, path="", out=None):
    """Recursive dict/list diff; floats compare with ``tol`` (exact if tol=0)."""
    if out is None:
        out = []
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            if k not in a:
                out.append((path + "/" + k, "<missing>", b[k]))
            elif k not in b:
                out.append((path + "/" + k, a[k], "<missing>"))
            else:
                _diff(a[k], b[k], tol, path + "/" + k, out)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            out.append((path, f"len={len(a)}", f"len={len(b)}"))
        else:
            for i, (x, y) in enumerate(zip(a, b)):
                _diff(x, y, tol, f"{path}[{i}]", out)
    elif (isinstance(a, (int, float)) and isinstance(b, (int, float))
          and not isinstance(a, bool) and not isinstance(b, bool)):
        if abs(a - b) > tol:
            out.append((path, a, b))
    elif a != b:
        out.append((path, a, b))
    return out


def _count_numbers(x):
    """How many numeric leaves a committed JSON block contains (audit count)."""
    if isinstance(x, dict):
        return sum(_count_numbers(v) for v in x.values())
    if isinstance(x, list):
        return sum(_count_numbers(v) for v in x)
    if isinstance(x, (int, float)) and not isinstance(x, bool):
        return 1
    return 0


def verify(args):
    """Hard self-check of the extended CSV (``ok=False`` => abort)."""
    problems, warnings = [], []
    n_cells = n_agg = 0

    if not os.path.exists(args.out):
        return {"ok": False, "problems": [f"{args.out} does not exist"],
                "warnings": [], "n_rows": 0}

    with open(args.out, newline="") as fh:
        reader = csv.DictReader(fh)
        fields = list(reader.fieldnames or [])
        csv_rows = [dict(r) for r in reader]

    # (a) structure -------------------------------------------------------
    if fields != list(core.CSV_COLUMNS):
        problems.append("CSV header != CSV_COLUMNS (order matters)")
    if len(csv_rows) != 60:
        problems.append(f"{len(csv_rows)} rows in the CSV instead of 60")
    ids = [r.get("id") or "" for r in csv_rows]
    if len(set(ids)) != len(ids):
        problems.append("duplicate measurement ids in the CSV")

    # (b) re-derivation from the committed raw files + vendored DB columns --
    by_part, prov = core.load_raw_parts()
    extra, extra_prov = core.db_extra_rows(args.db_extra)
    merged = core.merge_rows(by_part, extra)
    if len(merged) != len(csv_rows):
        problems.append(f"merged {len(merged)} rows vs {len(csv_rows)} CSV rows")
    strip_meta, strip_problems = attach_masks(merged, args.strips_dir, args.manifest)
    problems.extend(strip_problems)
    for i, (m, r) in enumerate(zip(merged, csv_rows)):
        for col in core.CSV_COLUMNS:
            n_cells += 1
            want = core._csv_cell(m.get(col))
            got = r.get(col, "")
            if got != want:
                problems.append(f"row {i + 1} ({m['part']} {core.day_of(m)}) column "
                                f"{col}: CSV {got!r} != derived {want!r}")
    if [core._csv_cell(m.get("id")) for m in merged] != ids:
        problems.append("row order of the CSV != the merged committed order")

    # (c) derived columns recomputed from the CSV's own day ----------------
    for r in csv_rows:
        day = (r.get("acquisition_ts") or "")[:10]
        if r.get("state") != core.state_of(day):
            problems.append(f"state of {day} is {r.get('state')!r}")
        if r.get("month") != str(int(day[5:7])) or r.get("year") != day[:4]:
            problems.append(f"month/year of {day} are {r.get('month')}/{r.get('year')}")
        season = core.st.season_of({"month": int(day[5:7])})
        if r.get("season") != season:
            problems.append(f"season of {day} is {r.get('season')!r} != {season!r}")
        px = core._num(r.get("coherence_masked_pixels"))
        if r.get("echo_mode") != core.echo_mode_of(px):
            problems.append(f"echo_mode of {day} is {r.get('echo_mode')!r} for A={px}")
        want_part = "part1" if day <= "2022-12-23" else "part2"
        if r.get("part") != want_part:
            problems.append(f"part of {day} is {r.get('part')!r}")

    # (d) the two committed channel JSONs must be reproduced from the CSV ----
    ref = core.load_reference()
    rows_all = core.load_csv(args.out)
    for p in core.PART_ORDER:
        reference = (core.part_stats(core.rows_of_part(rows_all, 1), 1, prov=prov)
                     if p == 2 else None)
        got = core.part_stats(core.rows_of_part(rows_all, p), p, prov=prov,
                              reference=reference)
        want = ref[f"kdlo_{p}_channels"]
        n_agg += _count_numbers(want)
        d = _diff(got, want, TOL)
        if d:
            for path, g, w in d[:10]:
                problems.append(f"part{p} channels JSON {path}: derived {g!r} != "
                                f"committed {w!r}")

    # (e) epoch accounting -------------------------------------------------
    counts = core.n_by_state(rows_all)
    if counts != core.N_BY_STATE_EXPECTED:
        problems.append(f"n_by_state {counts} != {core.N_BY_STATE_EXPECTED}")
    if len(core.state_rows(rows_all)) != sum(core.N_BY_STATE_EXPECTED.values()):
        problems.append("state_rows count != 59")
    singleton = [r for r in rows_all if r["state"] == "single"]
    if [core.day_of(r) for r in singleton] != [core.SINGLETON_DAY]:
        problems.append("the singleton acquisition is not exactly 2022-12-23")
    for label in core.STATE_LIST:
        lo, hi = core.STATE_WINDOWS[label]
        days = [core.day_of(r) for r in rows_all if r["state"] == label]
        if not all(lo <= d <= hi for d in days):
            problems.append(f"a {label} row lies outside {lo}..{hi}")

    # (f) provenance of the one non-committed input ------------------------
    meta = core.load_csv_meta(args.meta)
    if extra_prov["sha256"] != meta["db_extra"]["sha256"]:
        problems.append("db_extra sha256 in the meta != the vendored file")
    if extra_prov["sha256"] != sha256(args.db_extra):
        problems.append("db_extra sha256 of the vendored file changed")
    if extra_prov["columns"] != list(core.DB_EXTRA_COLUMNS):
        problems.append("vendored DB columns != DB_EXTRA_COLUMNS")
    if extra_prov["n_rows"] != len(rows_all):
        problems.append("vendored DB rows != CSV rows")
    n_missing = sum(1 for r in rows_all if r["id"] not in extra)
    if n_missing:
        problems.append(f"{n_missing} CSV rows have no vendored DB row")
    for p in core.PART_ORDER:
        if meta["raw_parts"][str(p)]["sha256"] != prov[p]["sha256"]:
            problems.append(f"raw part {p} sha256 in the meta != the committed file")
        if sha256(core.RAW_PATHS[p]) != prov[p]["sha256"]:
            problems.append(f"raw part {p} sha256 changed on disk")
    for name, fname in core.REFERENCE_FILES.items():
        if meta["reference_files"][name]["sha256"] != sha256(
                os.path.join(core.REF, fname)):
            problems.append(f"reference file {fname} sha256 changed")
    if meta["columns"] != fields:
        problems.append("meta columns != CSV header")
    if meta["n_by_state"] != counts:
        problems.append("meta n_by_state != recomputed")
    if meta["n_states_rows"] != len(core.state_rows(rows_all)):
        problems.append("meta n_states_rows != recomputed")

    # (g) the DB extraction must agree with the committed catalogue facts ---
    avail = ref["s1_availability"]
    doc = core.load_db_extra(args.db_extra)
    if doc["n_asset_rows_in_db"] - doc["n_rows"] != \
            avail["pipeline"]["measured_outside_figure_windows"]:
        problems.append("DB asset rows minus figure rows != the availability file's "
                        "measured_outside_figure_windows")
    if doc["asset_id"] != avail["pipeline"]["asset_id"]:
        problems.append("vendored DB asset_id != the availability file's asset_id")
    if doc["sql_sha256"] != hashlib.sha256(fetch_sql(
            [r["id"] for r in rows_all]).encode()).hexdigest():
        problems.append("vendored DB SQL is not the SQL this generator would run")

    # (h) the echo-mask layer: the CSV's geometry must be the geometry of the
    #     strip the pipeline decoded -------------------------------------
    n_masks, n_unmasked = 0, 0
    n_p_set = 0
    for r in merged:
        day = core.day_of(r)
        if r.get("A") is None:
            n_unmasked += 1
            problems.append(f"{day}: no echo mask (no strip in the cache)")
            continue
        n_masks += 1
        a_rec = core._num(r.get("coherence_masked_pixels"))
        g_rec = core._num(r.get("coherence_gamma2"))
        if a_rec is None or int(a_rec) != int(r["A"]):
            problems.append(f"{day}: A={r['A']} != coherence_masked_pixels={a_rec}")
        g_mask = core._num(r.get("gamma2_mask"))
        if g_rec is None or g_mask is None or abs(g_mask - g_rec) > TOL:
            problems.append(f"{day}: gamma2_mask={g_mask!r} != coherence_gamma2={g_rec!r}")
        bbox = core._num(r.get("bbox_area"))
        if not bbox or abs(core._num(r.get("D")) - core._num(r["A"]) / bbox) > TOL:
            problems.append(f"{day}: D != A / bbox_area")
        if (r.get("mask_rows"), r.get("mask_cols")) != (400, 11):
            problems.append(f"{day}: mask shape {r.get('mask_rows')}x"
                            f"{r.get('mask_cols')} != 400x11")
        p = core._num(r.get("P"))
        if p is None or not (0.0 <= p <= 1.0):
            problems.append(f"{day}: P={p!r} not in [0, 1]")
        else:
            n_p_set += 1
        if int(r.get("row_span") or 0) > 400 or int(r.get("col_span") or 0) > 11:
            problems.append(f"{day}: mask span {r.get('row_span')}x"
                            f"{r.get('col_span')} outside the strip")
        if core._num(r.get("F")) is not None and not (1 <= core._num(r["F"]) <= r["A"]):
            problems.append(f"{day}: F={r.get('F')} outside [1, A]")
    if n_masks != len(merged):
        problems.append(f"{n_masks} of {len(merged)} rows carry an echo mask")

    # (i) integrity of the committed strip cache --------------------------
    n_strips_hashed = 0
    if not args.no_strip_verify:
        for e in masks.load_manifest(args.manifest)["entries"]:
            path = os.path.join(args.strips_dir, e["strip_file"])
            if not os.path.exists(path):
                problems.append(f"strip cache: {e['strip_file']} is missing")
                continue
            n_strips_hashed += 1
            if sha256(path) != e["strip_sha256"]:
                problems.append(f"strip cache: {e['strip_file']} sha256 changed")

    # (j) provenance of the mask layer in the metadata --------------------
    if meta.get("mask_columns") != list(core.MASK_COLUMNS):
        problems.append("meta mask_columns != MASK_COLUMNS")
    sc = meta.get("strip_cache") or {}
    if sc.get("manifest_sha256") != sha256(args.manifest):
        problems.append("meta strip_cache.manifest_sha256 != the manifest on disk")
    if sc.get("n_strips") != len(masks.load_manifest(args.manifest)["entries"]):
        problems.append("meta strip_cache.n_strips != the manifest's entry count")
    if sc.get("n_attached_to_csv_rows") != n_masks:
        problems.append("meta strip_cache.n_attached_to_csv_rows != masks in the CSV")
    if sc.get("mask_rule", {}).get("min_n_masked") != masks.MIN_N_MASKED:
        problems.append("meta strip_cache.mask_rule.min_n_masked != the module's")

    return {
        "ok": not problems,
        "n_rows": len(csv_rows), "n_columns": len(fields),
        "n_cells_checked": n_cells,
        "n_aggregate_values_checked": n_agg,
        "n_echo_masks": n_masks,
        "n_rows_without_echo_mask": n_unmasked,
        "n_strips_hashed": n_strips_hashed,
        "n_strips_with_mask": strip_meta["n_with_mask"],
        "mask_shape": strip_meta["majority_shape"],
        "n_P_set": n_p_set,
        "problems": problems, "warnings": warnings,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=DEFAULT_OUT, help="target CSV")
    ap.add_argument("--meta", default=DEFAULT_META, help="provenance metadata (JSON)")
    ap.add_argument("--db-extra", default=DEFAULT_DB_EXTRA,
                    help="vendored one-time DB extraction (JSON)")
    ap.add_argument("--fetch-db-extra", action="store_true",
                    help="run the one-time read-only SELECT against tower-postgres")
    ap.add_argument("--verify-only", action="store_true",
                    help="write nothing, only check the existing CSV")
    ap.add_argument("--strips-dir", default=DEFAULT_STRIPS,
                    help="committed full-dwell strip cache (the echo-mask layer)")
    ap.add_argument("--manifest", default=DEFAULT_MANIFEST,
                    help="strip cache manifest (per-strip sha256 pins)")
    ap.add_argument("--no-strip-verify", action="store_true",
                    help="skip re-hashing the committed strips")
    args = ap.parse_args()

    if args.fetch_db_extra:
        fetch_db_extra(args)
        if args.verify_only:
            return

    if not args.verify_only:
        build(args)

    rep = verify(args)
    print(json.dumps(rep, indent=1, sort_keys=True))
    if not rep["ok"]:
        raise SystemExit("VERIFICATION FAILED")
    print("VERIFICATION OK")


if __name__ == "__main__":
    main()
