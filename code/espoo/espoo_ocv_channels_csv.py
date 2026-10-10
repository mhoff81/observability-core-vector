#!/usr/bin/env python3
"""Espoo OCV-Paper — build ``data/espoo/espoo_ocv_channels.csv``.

Why this generator exists
-------------------------
The LUMO package of this repository extends the committed overpass table with the
five mask dimensions ``A, D, F, S, P``, and its generator has to open a 100 MB
machine-local burst cache to do it. Espoo needs none of that, and the reason is a
property of the site's record, not a shortcut:

* the committed ``espoo_channels.csv`` **is** the site's own export — the same 150
  acquisitions that ``espoo_mast_observability.json`` commits under
  ``acquisitions``, in the same order, with exactly the same 27 fields;
* the export delivers exactly **one** mask-geometry quantity,
  ``coherence_masked_pixels``. The mask raster is not in the export, so there is
  no ``D`` (density), ``F`` (fragmentation), ``S`` (centroid <-> peak shift) or
  ``P`` (persistence) to compute: a generator that produced them would be
  inventing an echo-mask layer this record does not have.

The generator therefore does the honest amount of work: it **verifies** the copy
against the committed JSON field by field, and adds the bookkeeping the figure
pipeline needs but the export leaves out.

Added columns (the whole of the extension):

  A                = ``coherence_masked_pixels`` — the OCV name of the one
                     measured mask quantity. Copied, not recomputed, and proven
                     equal row by row against the committed column *and* against
                     the committed reference statistics. ``D``/``F``/``S``/``P``
                     are deliberately absent: no mask raster, no echo-mask layer.
  acquisition_date = ``YYYY-MM-DD`` of ``acquisition_ts`` — the column exists but
                     is empty on all 150 rows, so it is filled in place
  month            = ``YYYY-MM`` of ``acquisition_ts`` — same, empty, filled
  doy              = day of year, from ``acquisition_ts`` (appended)
  phase_available  = 1 where the phase ladder of the record delivered a number
                     (``phase_coherence`` present), 0 otherwise (appended)

Nothing else is added and nothing is dropped: the table keeps its 150 rows in its
committed order, because the bootstrap CIs (fixed seeds) and the permutation
tests depend on that order.

Usage:
  python3 code/espoo/espoo_ocv_channels_csv.py                # build + verify
  python3 code/espoo/espoo_ocv_channels_csv.py --verify-only  # verify only
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, os.pardir, os.pardir, "data", "espoo")
REF = os.path.join(DATA, "reference")

sys.path.insert(0, HERE)
import espoo_ocv_core as core  # noqa: E402

DEFAULT_SRC = core.RAW_CSV_PATH
DEFAULT_OUT = core.CSV_PATH
DEFAULT_META = core.CSV_META_PATH

OBS_JSON = os.path.join(REF, "espoo_mast_observability.json")
CHAN_JSON = os.path.join(REF, "fig_espoo_channels.json")

# The 27 fields ``espoo_mast_observability.json`` commits per acquisition — in
# the order the file lists them. The verifier compares the CSV against all of
# them; a column the export adds later shows up as a pin failure.
COMMITTED_FIELDS = [
    "id", "acquisition_ts", "pass_label", "orbit_direction",
    "displacement_los_m", "brightness_ratio", "measured_frequency_hz",
    "baseline_frequency_hz", "frequency_drop_pct", "wind_speed_ms",
    "temperature_c", "traffic_load_label", "status",
    "phase_detected_frequency_hz", "phase_coherence", "phase_status",
    "coherence_gamma2", "coherence_masked_pixels", "ndvi", "ndvi_scene_date",
    "ndvi_confound_flagged", "sub_aperture_modulation", "coherent_sum_amplitude",
    "registration_status", "phase_rms_rad", "condition_label", "campaign_label",
]

# The export carries the acquisition timestamps but leaves the two columns the
# pipeline wants **empty** on all 150 rows. The generator fills those in place
# rather than appending a second spelling of the same thing.
FILL_COLS = ["month", "acquisition_date"]

# Appended columns (order = documentation). ``A`` is the OCV name of the one
# measured mask quantity; ``doy`` and ``phase_available`` are bookkeeping the
# export does not have a slot for at all.
NEW_COLS = ["A", "doy", "phase_available"]

# ``fig_espoo_channels.json`` states the pooled and the per-orbit statistics of
# five record columns. The generator recomputes them from the source CSV: this is
# what proves the *copy* is the committed table and not merely a similar one.
CHANNEL_STATS_COLS = {
    "coherence_gamma2": "gamma2",
    "coherence_masked_pixels": "A",
    "phase_coherence": "phase_coherence",
    "phase_rms_rad": "phase_rms_rad",
    "phase_snr_db": "phase_snr_db",
}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _num(x):
    x = ("" if x is None else str(x)).strip()
    try:
        return float(x) if x else None
    except ValueError:
        return None


def _same(a, b, tol=1e-12):
    """Compare a CSV string against a committed JSON scalar."""
    if b is None:
        return a is None or str(a).strip() == ""
    if isinstance(b, bool):
        return str(a).strip().lower() == str(b).lower()
    if isinstance(b, (int, float)):
        av = _num(a)
        if av is not None:
            return abs(av - float(b)) <= tol * max(1.0, abs(float(b)))
    return str(a).strip() == str(b)


def read_source(path):
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        fieldnames = list(reader.fieldnames)
        rows = list(reader)
    return fieldnames, rows


def verify_copy(rows, obs, fieldnames):
    """Row by row, field by field: the CSV **is** the committed acquisitions."""
    committed = obs["acquisitions"]
    if len(rows) != len(committed):
        raise SystemExit(f"ERROR: {len(rows)} CSV rows vs "
                         f"{len(committed)} committed acquisitions")
    n_checked = 0
    for i, (got, want) in enumerate(zip(rows, committed)):
        for f in COMMITTED_FIELDS:
            if not _same(got.get(f), want.get(f)):
                raise SystemExit(
                    f"ERROR: row {i} ({got.get('id')}) field {f}: "
                    f"csv {got.get(f)!r} != committed {want.get(f)!r}")
            n_checked += 1
    # The two columns this generator fills must arrive empty — otherwise it
    # would silently overwrite a committed value.
    for col in FILL_COLS:
        if col not in fieldnames:
            raise SystemExit(f"ERROR: source table has no {col!r} column")
        n_filled_in_src = sum(1 for r in rows if (r.get(col) or "").strip())
        if n_filled_in_src:
            raise SystemExit(f"ERROR: source column {col!r} is not empty "
                             f"({n_filled_in_src} values)")
    for col in NEW_COLS:
        if col in fieldnames:
            raise SystemExit(f"ERROR: column {col!r} already exists in the "
                             f"source table")
    n_coherence = sum(1 for r in rows
                      if (r.get("coherence_gamma2") or "").strip())
    return {"n_rows": len(rows), "n_fields_checked": n_checked,
            "n_coherence_rows": n_coherence,
            "n_other_rows": len(rows) - n_coherence,
            "filled_columns": list(FILL_COLS)}


def derive_row(row):
    """The filled + appended columns of one source row."""
    ts = dt.datetime.fromisoformat(row["acquisition_ts"])
    return {
        "A": row["coherence_masked_pixels"],
        "acquisition_date": ts.date().isoformat(),
        "month": ts.strftime("%Y-%m"),
        "doy": str(ts.timetuple().tm_yday),
        "phase_available": "1" if (row.get("phase_coherence") or "").strip() else "0",
    }


def _median(vals):
    vals = sorted(v for v in vals if v is not None)
    if not vals:
        return None
    n = len(vals)
    return vals[n // 2] if n % 2 else 0.5 * (vals[n // 2 - 1] + vals[n // 2])


def channel_stats(rows):
    """The five blocks of ``fig_espoo_channels.json``, recomputed from the CSV."""
    def block(col, sel=None):
        vals = [_num(r.get(col)) for r in rows
                if sel is None or r["orbit_direction"] == sel]
        vals = [v for v in vals if v is not None]
        return {"n": len(vals), "median": _median(vals)}

    out = {"n_acquisitions": len(rows)}
    for col, key in CHANNEL_STATS_COLS.items():
        out[key] = block(col)
        for orb, suffix in (("ASCENDING", "_ascending"),
                            ("DESCENDING", "_descending")):
            out[key + suffix] = block(col, orb)
    return out


def verify_channels(rows, chan):
    """Recomputed pooled/per-orbit statistics vs. the committed channel JSON."""
    got = channel_stats(rows)
    n_checked = 0

    def check_block(alias, committed):
        nonlocal n_checked
        want = chan[committed]
        g = got[alias]
        for kk, wv in want.items():
            gv = g[kk]
            if isinstance(wv, float) and isinstance(gv, float):
                ok = abs(gv - wv) <= 1e-12 * max(1.0, abs(wv))
            else:
                ok = gv == wv
            if not ok:
                raise SystemExit(f"ERROR: {committed}.{kk}: recomputed {gv!r} "
                                 f"!= committed {wv!r}")
            n_checked += 1

    if chan["n_acquisitions"] != got["n_acquisitions"]:
        raise SystemExit("ERROR: n_acquisitions disagrees")
    n_checked += 1
    for col, alias in CHANNEL_STATS_COLS.items():
        for suffix in ("", "_ascending", "_descending"):
            check_block(alias + suffix, col + suffix)
    # The committed LUMO reference scale the verdict is read against — the same
    # constants the classification layer of ``espoo_ocv_stats`` uses, so this
    # proves the two files agree on the scale.
    lumo = {"ascending_gamma2": core.st.LUMO_ASC_G2_MEDIAN,
            "descending_gamma2": core.st.LUMO_DESC_G2_MEDIAN}
    for k, want in chan["lumo_reference"].items():
        if k not in lumo or abs(lumo[k] - want) > 1e-12:
            raise SystemExit(f"ERROR: lumo_reference.{k}: package constant "
                             f"{lumo.get(k)!r} != committed {want!r}")
        n_checked += 1
    # The monthly composition of the record (both orbits, all 25 months).
    monthly = {m: {"ASCENDING": 0, "DESCENDING": 0} for m in chan["monthly"]}
    for r in rows:
        monthly[r["month"]][r["orbit_direction"]] += 1
    for m, want in chan["monthly"].items():
        for orb in ("ASCENDING", "DESCENDING"):
            if monthly[m][orb] != want[orb]:
                raise SystemExit(f"ERROR: monthly.{m}.{orb}: recomputed "
                                 f"{monthly[m][orb]} != committed {want[orb]}")
            n_checked += 1
    dates = sorted(r.get("acquisition_date") or r["date"] for r in rows)
    for key, want in (("date_first", dates[0]), ("date_last", dates[-1])):
        if chan[key] != want:
            raise SystemExit(f"ERROR: {key}: recomputed {want!r} != "
                             f"committed {chan[key]!r}")
        n_checked += 1
    # ``A`` is a copy of the committed column: prove the two agree value by
    # value, so the statistics above cannot pass on a half-copied table.
    a_vals = [_num(r.get("A")) for r in rows]
    src_vals = [_num(r.get("coherence_masked_pixels")) for r in rows]
    if a_vals != src_vals:
        raise SystemExit("ERROR: A and coherence_masked_pixels disagree")
    n_checked += len(rows)
    return n_checked, monthly


def build(args):
    fieldnames, rows = read_source(args.csv)
    with open(OBS_JSON) as fh:
        obs = json.load(fh)
    with open(CHAN_JSON) as fh:
        chan = json.load(fh)

    copy_check = verify_copy(rows, obs, fieldnames)

    # Fill the empty columns and append the new ones, then verify the channel
    # statistics — those are computed from the derived table, i.e. from exactly
    # what the figure pipeline will read.
    out_names = fieldnames + NEW_COLS
    for r in rows:
        d = derive_row(r)
        for col in FILL_COLS + NEW_COLS:
            r[col] = d[col]

    n_chan_checks, monthly = verify_channels(rows, chan)

    if args.verify_only:
        print(f"verify only: {copy_check['n_rows']} rows, "
              f"{copy_check['n_fields_checked']} copied fields, "
              f"{n_chan_checks} channel checks — all ok")
        return 0

    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=out_names)
        writer.writeheader()
        writer.writerows(rows)

    # Read the artifact back and check it is what the pipeline will load.
    loaded = core.load_csv(args.out)
    meta = {
        "generator": "code/espoo/espoo_ocv_channels_csv.py",
        "source_csv": os.path.relpath(args.csv, os.path.join(HERE, os.pardir, os.pardir)),
        "output_csv": os.path.relpath(args.out, os.path.join(HERE, os.pardir, os.pardir)),
        "source_sha256": sha256(args.csv),
        "output_sha256": sha256(args.out),
        "committed_acquisitions_json": "data/espoo/reference/espoo_mast_observability.json",
        "committed_channels_json": "data/espoo/reference/fig_espoo_channels.json",
        "filled_columns": FILL_COLS,
        "new_columns": NEW_COLS,
        "mask_definition": {
            "A": "coherence_masked_pixels — copied from the committed record, not "
                 "recomputed: the export carries no mask raster, so A is the "
                 "pipeline's own mask size and there is no D/F/S/P of this site.",
            "acquisition_date": "YYYY-MM-DD of acquisition_ts (column was empty)",
            "month": "YYYY-MM of acquisition_ts (column was empty)",
            "doy": "day of year of acquisition_ts",
            "phase_available": "1 where phase_coherence is present, else 0",
        },
        "derived_checks": {
            "rows_read": copy_check["n_rows"],
            "copied_fields_verified_against_committed_json":
                copy_check["n_fields_checked"],
            "coherence_rows": copy_check["n_coherence_rows"],
            "rows_left_empty": copy_check["n_other_rows"],
            "channel_stats_verified": n_chan_checks,
            "columns_filled_in_place": len(FILL_COLS),
            "columns_appended": len(NEW_COLS),
            "n_months": len(monthly),
            "loaded_by_core": len(loaded),
            "A_equals_committed_masked_pixels": True,
        },
        "monthly": monthly,
        "n_by_orbit": core.n_by_state(loaded),
        "date_range": [min(r["date"] for r in loaded),
                       max(r["date"] for r in loaded)],
    }
    with open(args.meta, "w", encoding="utf-8") as fh:
        json.dump(meta, fh, indent=1, sort_keys=True)
        fh.write("\n")

    print(f"wrote {os.path.relpath(args.out)}: {len(rows)} rows, "
          f"+{len(NEW_COLS)} columns ({', '.join(NEW_COLS)})")
    print(f"      {copy_check['n_fields_checked']} copied fields verified, "
          f"{n_chan_checks} channel checks verified")
    print(f"wrote {os.path.relpath(args.meta)}")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Build data/espoo/espoo_ocv_channels.csv from the committed "
                    "site record table + the committed reference JSONs.")
    ap.add_argument("--csv", default=DEFAULT_SRC, help="source record table")
    ap.add_argument("--out", default=DEFAULT_OUT, help="extended table to write")
    ap.add_argument("--meta", default=DEFAULT_META, help="metadata JSON to write")
    ap.add_argument("--verify-only", action="store_true",
                    help="run every check, write nothing")
    args = ap.parse_args(argv)
    return build(args)


if __name__ == "__main__":
    raise SystemExit(main())
