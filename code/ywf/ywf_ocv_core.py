#!/usr/bin/env python3
"""YWF OCV-Paper — data layer, constants and the loaders of the committed tables.

Everything the other modules need, in one place: paths, the analysis channels,
the mask vector, the site constants, the event and the loaders.

The site
--------
Yeongdeok Wind Farm (Changpo Wind Power Complex, Samgye-ri, 36.4236 N
129.4203 E) is Unit 21, a steel monopole wind-turbine tower that **collapsed on
2026-02-02**, which divides the record into the two states of the *same*
structure

    pre-collapse (healthy)   day <  2026-02-02  <=  post-collapse

The record straddles the event, so the season moves with the state (the pre side
is Aug-Jan, the post side Feb-Apr) and every contrast below carries a co-variate
control beside it. The record is small: 57 committed 7x7 windows on 32 dates,
33 of them unique chips — and only **6** unique chips after the event, of which
**1** carries an echo at all. The collapse contrast of this site is therefore
reported as an **echo-coverage** statement with an explicit power caveat, never
as a damage verdict.

The export
----------
Per acquisition the pipeline stores one 7x7 complex window per **mast-section
request** (``segment_index`` 0..4 -> ``turbine-mast-section-0`` .. ``-4``, see
``data/ywf/ywf_segments.txt``), but the committed record only holds segments 2
and 3, and 24 of its 57 windows are byte-identical duplicates that differ only
in that request index: both requests decoded the **same** chip at the asset
point, so the record is **not segment-resolved** (see ``ywf_ocv_masks``). This
package consequently carries **one** mask layer, the echo mask of the whole
7x7 window — a per-mast-section layer would be dishonest.

The mask layer
--------------
  * ``echo mask`` — the whole 7x7 window, the layer of the 6D vector
    ``x = [gamma2, P, D, A, F, S]``, recomputed from the committed payload
    ``data/ywf/ywf_windows_full.txt`` (38 kB, committed in full) and re-verified
    from the committed cache ``data/ywf/ywf_windows_mask_cache.txt``;
  * the mask rule is the project-wide echo definition (``0.30 x peak`` **and**
    ``5 x np.median``, ``min_n_masked = 2``), identical to Carola's.

``A`` is *not* the pipeline's own ``coherence_masked_pixels``: that column was
computed on the pre-fix asset-point chip of the whole-scene window (75, 133, ...
px), while this package recomputes the echo mask of the committed 7x7 window
(2..7 px). The CSV keeps both and reports how many rows coincide; see
``code/ywf/README.md``.

Reference file (``data/ywf/reference/``)
----------------------------------------
  ywf_ocv_findings.json   the committed analysis-layer reference: the
                          echo-coverage table and its Fisher test, the mask
                          dimensions per state, the pipeline-column and
                          co-variate controls, the timeline and the record's
                          provenance facts. This site has **no upstream analysis
                          project** to pin against (the record *is* the site's
                          own export), so the reference is regenerated
                          deterministically by ``fig_ywf_ocv_paper.py
                          --write-reference`` and pinned on every later run.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import ywf_ocv_masks as masks  # noqa: E402
import ywf_ocv_stats as st     # noqa: E402

import numpy as np  # noqa: E402

DATA = os.path.join(HERE, os.pardir, os.pardir, "data", "ywf")
FIGDIR = os.path.join(HERE, os.pardir, os.pardir, "figures", "ywf")
REFERENCE = os.path.join(DATA, "reference")

CACHE_PATH = masks.CACHE_PATH
MANIFEST_PATH = masks.MANIFEST_PATH
MEAS_PATH = masks.MEAS_PATH
INDEX_PATH = masks.INDEX_PATH
SEGMENTS_PATH = masks.SEGMENTS_PATH
WINDOWS_PATH = masks.DEFAULT_WINDOWS

CSV_PATH = os.path.join(DATA, "ywf_ocv_channels.csv")
CSV_META_PATH = os.path.join(DATA, "ywf_ocv_channels_meta.json")
REF_FINDINGS = os.path.join(REFERENCE, "ywf_ocv_findings.json")

SITE = "Yeongdeok Wind Farm (Samgye-ri)"
SITE_SHORT = "YWF"
ASSET_ID = "af683c7e-416a-4395-86cc-094ee83d9497"
ASSET_NAME = "Yeongdeok Wind Turbine (Samgye-ri)"
COLLAPSE_DATE = "2026-02-02"
EVENT = COLLAPSE_DATE

STATES = list(st.STATE_ORDER)
PRE, POST = STATES
STATE_SHORT = dict(st.STATE_SHORT)

OCV = ["gamma2", "P", "D"]
FEATURES = ["gamma2", "P", "D", "A", "F", "S"]
DIMS = list(FEATURES)

# The mast-section requests of the export. Section 0/1/4 never decoded a chip
# in this record; the constants stay complete so the reader sees the whole
# segment table of ``ywf_segments.txt``.
SECTIONS = [f"turbine-mast-section-{i}" for i in range(5)]
SEGMENT_INDICES = {lab: i for i, lab in enumerate(SECTIONS)}
CHANNEL_AXIS = {
    "gamma2": "gamma^2  (masked coherence)",
    "P": "P  (mask persistence)",
    "D": "D  (mask density A / bbox)",
    "A": "A  (echo pixels n_masked)",
    "F": "F  (8-connected components)",
    "S": "S  (||centroid - peak||, px)",
    "wind_speed_ms": "wind (m/s)",
    "temperature_c": "temperature (C)",
    "precipitation_mm": "precip (mm)",
    "wind_gust_ms": "gust (m/s)",
    "brightness_ratio": "brightness ratio",
    "intensity": "intensity",
    "coherence": "coherence",
    "coherence_gamma2": "gamma2 (export)",
    "cumulative_mm": "cumulative d_LOS (mm)",
    "structural_frequency_hz": "structural freq (Hz)",
}
CHANNEL_TITLE = {
    "gamma2": "gamma2", "P": "P", "D": "D", "A": "A", "F": "F", "S": "S",
    "wind_speed_ms": "wind (m/s)", "temperature_c": "temperature (C)",
    "precipitation_mm": "precip (mm)", "wind_gust_ms": "gust (m/s)",
    "brightness_ratio": "brightness ratio", "intensity": "intensity",
    "coherence": "coherence", "coherence_gamma2": "gamma2 (export)",
    "cumulative_mm": "cumulative d_LOS (mm)",
    "structural_frequency_hz": "structural freq (Hz)",
}

WEATHER = ["wind_speed_ms", "temperature_c", "precipitation_mm", "wind_gust_ms",
           "relative_humidity_pct"]
# The pipeline's own columns beside the weather co-variates — the controls the
# figures use to show *what else* moves with the state.
CONTROLS = ["wind_speed_ms", "temperature_c", "brightness_ratio", "intensity",
            "coherence", "coherence_gamma2", "cumulative_mm",
            "structural_frequency_hz"]
PIPELINE_CHANNELS = ["intensity", "brightness_ratio", "coherence",
                     "coherence_gamma2", "coherence_masked_pixels",
                     "cumulative_mm", "structural_frequency_hz",
                     "fused_frequency_hz"]

MODELS = {
    "gamma2": ["gamma2"],
    "PD": ["P", "D"],
    "gamma2PD": ["gamma2", "P", "D"],
    "x_6d": ["gamma2", "P", "D", "A", "F", "S"],
    "weather": ["wind_speed_ms", "temperature_c", "precipitation_mm",
                "wind_gust_ms"],
}
HEADLINE_PAIR = "x_6d"
LAMBDA_GRID = [0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9]
HEADLINE_LAMBDA = 0.5
RNG_SEED = 7
N_BOOT = 10000
N_PERM = 1000

MASK_COLUMNS = ["A", "D", "F", "S", "P", "gamma2_mask", "gamma2_export",
                "mask_bbox_area", "mask_row_span", "mask_col_span",
                "mask_centroid_row", "mask_centroid_col",
                "mask_peak_row", "mask_peak_col", "mask_n_components",
                "mask_largest_component", "mask_fragmentation"]

SOURCE_NUMERIC = ["intensity", "displacement_los_m", "brightness_ratio",
                  "measured_frequency_hz", "baseline_frequency_hz",
                  "frequency_drop_pct", "wind_speed_ms", "temperature_c",
                  "precipitation_mm", "wind_gust_ms", "relative_humidity_pct",
                  "incidence_angle_deg", "amplitude_expected_m",
                  "amplitude_measured_m", "thermal_frequency_adjusted_hz",
                  "coherence", "cumulative_mm", "structural_frequency_hz",
                  "dynamic_traffic_frequency_hz", "phase_detected_frequency_hz",
                  "phase_coherence", "fused_frequency_hz", "fused_alpha",
                  "fusion_confidence", "phase_snr_db", "phase_scatterer_count",
                  "phase_modal_consistency", "phase_quality_score",
                  "phase_rms_rad", "coherence_masked_pixels",
                  "coherence_gamma2", "sub_aperture_modulation",
                  "coherent_sum_amplitude"]
SOURCE_TEXT = ["asset_id", "asset_name", "asset_type", "request_id",
               "segment_index", "acquisition_ts", "pass_label",
               "orbit_direction", "condition_label", "campaign_label",
               "traffic_load_label", "status", "phase_observable",
               "phase_unobservable_reason", "fusion_source", "fusion_route",
               "registration_status", "ndvi_scene_date"]
SOURCE_COLUMNS = ["id"] + SOURCE_TEXT + SOURCE_NUMERIC
# ``section_index`` is the girder-index analogue of the bridge packages; the
# last three derived columns are YWF-specific provenance of the duplicated
# segment payloads (the property that makes this record not segment-resolved).
DERIVED_COLUMNS = ["segment", "section_index", "state", "pre_collapse", "day",
                   "month", "year", "season", "echo_mode", "echo_mask_present",
                   "payload_group_size", "payload_shared_with_segments"]
CSV_COLUMNS = SOURCE_COLUMNS + DERIVED_COLUMNS + MASK_COLUMNS


def _num(v):
    """CSV/JSON cell -> float or None."""
    if v is None or v == "" or v == "null" or v == "None":
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _bool(v):
    return str(v) in ("1", "True", "true")


def day_of(row):
    """Acquisition day (YYYY-MM-DD) of a CSV row."""
    d = row.get("day")
    return d if d else (row.get("acquisition_ts") or "")[:10]


def state_of(day):
    """The collapse alphabet: ``pre-collapse (healthy)`` before the event."""
    return PRE if str(day) < COLLAPSE_DATE else POST


def echo_mode_of(a):
    """Echo-mode bucket of a masked-pixel count, on the site's own scale.

    The project-wide cuts (compact <= 10 px, distributed >= 25 px) belong to the
    80x80 bridge chips; a YWF window holds 49 px in total, so those cuts would
    put every echo in one bucket. This site keeps the three names on its own
    scale — 2 px, 3..5 px, 6+ px of 49 — and the metadata states the cut used.
    """
    a = _num(a)
    if a is None:
        return None
    if a <= 2:
        return "compact"
    if a <= 5:
        return "intermediate"
    return "distributed"


ECHO_MODE_CUTS = {"compact": "A <= 2", "intermediate": "3 <= A <= 5",
                  "distributed": "A >= 6", "window_px": 49}


def payload_digest(payload):
    """The cache's payload digest: ``md5`` of the committed payload text."""
    return hashlib.md5(payload.encode()).hexdigest()[:12]


def read_windows(path=WINDOWS_PATH):
    """``[(mid, w, h, payload_text)]`` of the committed window file."""
    out = []
    if not os.path.isfile(path):
        return out
    for line in open(path).read().splitlines():
        if not line.strip():
            continue
        mid, w, h, payload = line.split("|", 3)
        out.append((mid, int(w), int(h), payload))
    return out


def chip_table(path=WINDOWS_PATH):
    """``(windows, unique)`` — the committed 7x7 windows as chip records.

    Every record carries the acquisition metadata of the index, the payload
    digest, the decoded ``z`` chip (``(h, w)`` complex) and the re-derived echo
    mask (``masks.echo_mask``, ``None`` when the chip has no echo). The
    ``payload_*`` fields are the segment-duplication provenance: windows of one
    ``(asset, day)`` whose payloads are byte identical share a group, and the
    record is *not* segment-resolved (see the module docstring).

    ``windows`` keeps every committed line in file order; ``unique`` keeps the
    first occurrence of each ``(asset, day, payload)`` — the mask layer's
    reference set, i.e. the rows of ``ywf_windows_mask_cache.txt``.
    """
    index = masks.read_index(INDEX_PATH)
    meas = masks.read_measurements(MEAS_PATH)
    groups = {}
    for mid, _, _, payload in read_windows(path):
        meta = index.get(mid)
        if meta is None:
            continue
        groups.setdefault((meta["asset"], meta["date"]), []).append(
            (segment_index_of(meta["segment"]), payload_digest(payload)))

    windows, seen = [], {}
    for mid, w, h, payload in read_windows(path):
        meta = index.get(mid)
        if meta is None:
            continue
        digest = payload_digest(payload)
        group = {s for s, d in groups[(meta["asset"], meta["date"])]
                 if d == digest}
        samples = [complex(p[0], p[1]) for p in json.loads(payload)]
        z = np.array(samples, dtype=complex).reshape(h, w)
        m = meas.get(mid) or {}
        res = masks.echo_mask(z)
        sect = segment_index_of(meta["segment"])
        rec = {
            "mid": mid, "asset": meta["asset"], "asset_id": meta["asset_id"],
            "segment": section_label_of(sect) or meta["segment"],
            "section_index": sect,
            "date": meta["date"], "ts": meta["ts"],
            "request_id": meta["request_id"], "w": w, "h": h,
            "digest": digest, "z": z,
            "px": res["px"] if res else [],
            "orbit": m.get("orbit_direction"), "pass_label": m.get("pass_label"),
            "state": state_of(meta["date"]),
            "payload_group_size": len(group),
            "payload_shared_with_segments": sorted(group),
            "mask": res,
        }
        key = (meta["asset"], meta["date"], digest)
        rec["first_of_group"] = key not in seen
        seen.setdefault(key, mid)
        windows.append(rec)
    unique = [r for r in windows if r["first_of_group"]]
    add_persistence(windows, unique)
    return windows, unique


def add_persistence(windows, unique):
    """``P`` of every window, from the de-duplicated reference set.

    ``P`` is the mean pixel frequency inside the own mask over the unique chips
    (``ywf_ocv_masks.add_persistence``), so the two byte-identical windows of a
    duplicated segment pair necessarily get the same ``P``: the duplication
    cannot inflate the persistence of a pixel.
    """
    ref = [r for r in unique if r["mask"] is not None]
    masks.add_persistence(ref)
    by_group = {(r["date"], r["digest"]): r["P"] for r in ref}
    for r in windows:
        r["P"] = (None if r["mask"] is None
                  else by_group.get((r["date"], r["digest"])))
        r["P_reference_chips"] = len(ref)
    return windows


def section_index_of(label):
    """0..4 of a mast-section label (``turbine-mast-section-N``)."""
    return SEGMENT_INDICES.get(label)


def section_label_of(index):
    """The label of a mast-section request index (0..4) of the export."""
    idx = segment_index_of(index)
    return SECTIONS[idx] if idx is not None and 0 <= idx < len(SECTIONS) else None


def segment_index_of(value):
    """The numeric mast-section request index of a ``segment`` cell.

    The committed window index stores the request's raw ``segment_index`` (the
    cache rows and the CSV rows both carry it as text), so the label of a row is
    ``SECTIONS`` looked up with this number.
    """
    try:
        return int(float(str(value).strip()))
    except (TypeError, ValueError):
        return None


def cache_shaped(rec):
    """A cache row of a chip record — the input of ``masks.vector_from_row``."""
    row = dict(rec["mask"] or {})
    row.update({"mid": rec["mid"], "asset": rec["asset"],
                "segment": rec["segment"], "date": rec["date"]})
    return row


def vector_of(rec):
    """The derived 6D vector of a chip record (``None`` when it has no echo)."""
    if rec["mask"] is None:
        return None
    v = masks.vector_from_row(cache_shaped(rec))
    v["P"] = rec.get("P")
    return v


# ---------------------------------------------------------------------------
# Loaders of the committed tables
# ---------------------------------------------------------------------------
def load_cache(path=CACHE_PATH):
    """The committed window cache (the mask layer's committed input)."""
    return masks.load_cache(path)


def load_manifest(path=MANIFEST_PATH):
    return masks.load_manifest(path)


def load_measurements(path=MEAS_PATH):
    """``measurement_id -> column dict`` of the committed measurement extract."""
    return masks.read_measurements(path)


def load_index(path=INDEX_PATH):
    """``measurement_id -> index row`` of the committed window index."""
    return masks.read_index(path)


def load_segments(path=SEGMENTS_PATH):
    """``(asset_id, segment_index) -> row`` of the committed segment table."""
    out = {}
    if not os.path.isfile(path):
        return out
    lines = open(path).read().splitlines()
    if not lines:
        return out
    keys = lines[0].split("|")
    for line in lines[1:]:
        p = line.split("|")
        if len(p) != len(keys):
            continue
        d = dict(zip(keys, p))
        d["segment_index"] = int(d["segment_index"])
        d["fundamental_hz"] = _num(d.get("fundamental_hz"))
        out[(d["asset_id"], d["segment_index"])] = d
    return out


def load_csv(path=CSV_PATH):
    """The generated channel table, numeric cells as floats."""
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        for c in SOURCE_NUMERIC + MASK_COLUMNS + ["payload_group_size"]:
            r[c] = _num(r.get(c))
        r["gamma2"] = r["gamma2_mask"]
        r["section_index"] = _num(r.get("section_index"))
        r["pre_collapse"] = _bool(r.get("pre_collapse"))
        r["echo_mask_present"] = _bool(r.get("echo_mask_present"))
    return rows


def load_csv_meta(path=CSV_META_PATH):
    with open(path) as fh:
        return json.load(fh)


def load_reference(path):
    with open(path) as fh:
        return json.load(fh)


def load_findings(path=REF_FINDINGS):
    """The committed analysis-layer reference (``data/ywf/reference/``)."""
    return load_reference(path)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Small aggregate helpers the figure script shares with the generator
# ---------------------------------------------------------------------------
def row_brief(r):
    """The columns of one chip that enter every figure and every pin."""
    return {
        "id": r.get("id"), "day": day_of(r), "state": r.get("state"),
        "segment": r.get("segment"), "orbit": r.get("orbit_direction"),
        "gamma2": _num(r.get("gamma2_mask")), "P": _num(r.get("P")),
        "D": _num(r.get("D")), "A": _num(r.get("A")), "F": _num(r.get("F")),
        "S": _num(r.get("S")),
        "wind_speed_ms": _num(r.get("wind_speed_ms")),
        "temperature_c": _num(r.get("temperature_c")),
        "precipitation_mm": _num(r.get("precipitation_mm")),
        "wind_gust_ms": _num(r.get("wind_gust_ms")),
        "coherence": _num(r.get("coherence")),
    }


def summ(rows, key):
    """``st._stats`` over one channel — the CSV's own summary shape."""
    return st._stats([r.get(key) for r in rows])


def feature_column(k):
    """The CSV column that carries feature ``k`` (``gamma2`` is stored as
    ``gamma2_mask``; the pipeline's own ``gamma2_export`` lives beside it)."""
    return "gamma2_mask" if k == "gamma2" else k


def summary_by_state(rows):
    """``{state_short: {feature: st._stats}}`` — the metadata's summary block."""
    return {STATE_SHORT[s]: {k: summ([r for r in rows if r.get("state") == s],
                                     feature_column(k))
                             for k in FEATURES}
            for s in STATES}


def by_section(rows):
    """Rows per mast-section label, in ``SECTIONS`` order (unknowns appended)."""
    keys = list(SECTIONS)
    for r in rows:
        s = r.get("segment")
        if s and s not in keys:
            keys.append(s)
    return {g: [r for r in rows if r.get("segment") == g] for g in keys}


# The section plays the role of Carola's girder; the alias keeps the figure
# scripts of the six packages readable as one diff.
by_girder = by_section


def fisher_p(a, b, c, d):
    """Two-sided Fisher exact p of the 2x2 table ``[[a, b], [c, d]]``."""
    if st.sps is None:
        return None
    return float(st.sps.fisher_exact([[a, b], [c, d]])[1])


def coverage_block(chips, event=COLLAPSE_DATE):
    """Echo coverage of the *unique* chips per state, plus its Fisher test.

    ``chips`` is the ``unique`` list of :func:`chip_table`. A chip "carries an
    echo" when the mask rule holds on it (``A >= 2``); the block is the site's
    collapse signature, and the Fisher test states how far 6 post-event chips
    can carry it.
    """
    table, block = {}, {}
    for s in STATES:
        sub = [c for c in chips if c["state"] == s]
        n_echo = sum(1 for c in sub if c["mask"] is not None)
        table[s] = [n_echo, len(sub) - n_echo]
        block[STATE_SHORT[s]] = {
            "n_chips": len(sub), "n_with_echo": n_echo,
            "n_without_echo": len(sub) - n_echo,
            "echo_rate": (n_echo / len(sub)) if sub else None,
        }
    order = [STATE_SHORT[s] for s in STATES]
    return {
        "event": event,
        "window_px": 49,
        "n_windows": len(chips),
        "per_state": block,
        "fisher": {
            "table_rows": order, "table_cols": ["echo", "no_echo"],
            "table": [table[s] for s in STATES],
            "p_two_sided": fisher_p(table[PRE][0], table[PRE][1],
                                    table[POST][0], table[POST][1]),
        },
    }


# ---------------------------------------------------------------------------
# CLI — a read-only brief of the committed layer
# ---------------------------------------------------------------------------
def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--chips", action="store_true",
                    help="one line per unique chip with its 6D vector")
    args = ap.parse_args(argv)

    windows, unique = chip_table()
    cov = coverage_block(unique)
    manifest = load_manifest()
    print(f"{SITE} — event {COLLAPSE_DATE}")
    print(f"  {len(windows)} committed windows on "
          f"{len({c['date'] for c in windows})} dates, "
          f"{len(unique)} unique chips")
    for k in ("pre", "post"):
        b = cov["per_state"][k]
        print(f"  {k:>4}: {b['n_chips']:3d} unique chips, "
              f"{b['n_with_echo']:3d} with an echo "
              f"({(b['echo_rate'] or 0.0) * 100:.1f} %)")
    p = cov["fisher"]["p_two_sided"]
    print(f"  Fisher exact (echo vs. state): p = "
          f"{'n/a (scipy missing)' if p is None else f'{p:.4f}'}")
    print(f"  not segment-resolved: "
          f"{manifest.get('n_deduplicated_same_payload')} duplicated windows, "
          f"segment_resolved={manifest.get('segment_resolved')}")
    if args.chips:
        print()
        print(f"{'date':10s} {'seg':>3s} {'A':>2s} {'D':>5s} {'F':>2s} "
              f"{'S':>5s} {'P':>5s} {'gamma2':>6s}  echo_mode")
        for c in sorted(unique, key=lambda r: (r["date"], r["segment"])):
            v = vector_of(c)
            if v is None:
                print(f"{c['date']:10s} {c['segment']:>3s}  - no echo")
                continue
            print(f"{c['date']:10s} {c['segment']:>3s} {v['A']:2d} "
                  f"{v['D']:5.2f} {v['F']:2d} {v['S']:5.2f} {v['P']:5.3f} "
                  f"{v['gamma2']:6.3f}  {echo_mode_of(v['A'])}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
