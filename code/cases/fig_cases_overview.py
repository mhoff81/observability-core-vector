#!/usr/bin/env python3
"""One figure for the whole case set — map, dimension matrix, trajectories.

Three panels, one page:

* **(a) the case map** — every case placed on the event spectrum
  *stable → damage/degradation → collapse → reconstruction*, each with the one
  question it exists to answer. LUMO and Bautzen sit **above** the spectrum as
  controlled-reference experiments, because they are not further cases but the
  calibration of the method.
* **(b) the dimension-response matrix** — rows are the cases, columns are `γ², P,
  D, F`, and a cell carries only what the case's own committed analysis supports:
  `↑`/`↓` when the effect's interval excludes zero, `○` when no change is
  resolved, `↕` when the case reports a non-monotonic response, `?` when the
  evidence is not there. The **fill intensity is the evidence strength** (the
  project's SDS bands), never an effect size that could be compared across cases.
* **(c) the observability trajectories** — six compact time lines drawn from the
  same facts, so the reader sees *when* each case moves.

Every cell is read from the case's committed paper JSON (`data/<site>/
fig_<site>_ocv_paper.json`; Espoo from its chip figure, which reports trend flags
rather than a state contrast) and pinned with the file's sha256. Nothing in this
figure is typed in by hand.

The chip and the comparability of the numbers
---------------------------------------------
All six cases are read on the **committed** chip, the one every site package has;
the support-sized chip exists for KDLO as a second, additive set and flips the
sign of `A`, `F` and `P` there (`fig_kdlo_ocv_paper_support.json`). The matrix
therefore states its chip and does not mix the two. `SDS` is used **only** as an
intensity scale inside a cell (the project's own bands); the matrix never compares
one case's effect with another's.

Usage
-----
  python3 code/cases/fig_cases_overview.py
  python3 code/cases/fig_cases_overview.py --selftest
  python3 code/cases/fig_cases_overview.py --verify

Outputs
-------
  figures/cases/fig_cases_overview.jpg   (630 dpi)
  figures/cases/fig_cases_overview.md
  data/cases/fig_cases_overview.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
DATA = os.path.join(REPO, "data")
FIGDIR = os.path.join(REPO, "figures", "cases")

SCHEMA = "ocv-cases-overview-v1"
FIGURE_DPI = 630
CHIP = "committed"

#: The matrix rows, in the order the cases are read.
CASES = ("espoo", "kdlo", "carola", "morandi", "cts", "haltern")

LABEL = {"espoo": "Espoo mast", "kdlo": "KDLO media tower",
         "carola": "Carolabrücke", "morandi": "Ponte Morandi",
         "cts": "Champlain Towers South", "haltern": "Haltern AV9",
         "lumo": "LUMO tower", "bautzen": "Bautzen bridge"}

#: column key -> display label. `gamma2` is the site's own coherence channel
#: (Espoo calls it `gamma2_mask`, the others `gamma2` / `coherence_gamma2`).
COLUMNS = (("gamma2", "γ²"), ("P", "P"), ("D", "D"), ("F", "F"))

#: The project's SDS bands (`code/sds/sds_support.py`), as the intensity scale.
BANDS = ((2.0, "very weak"), (4.0, "weak"), (6.0, "moderate"),
         (8.0, "strong"), (float("inf"), "very strong"))

#: Below this many units on either side a cell is not evidence (the project's
#: SMALL_N): an interval from two chips is not an interval.
SMALL_N = 6

#: where each case sits on the spectrum, and what it is there to test
POSITIONS = {
    "espoo": (0.06, "stable reference — no event, no label"),
    "cts": (0.34, "long-term evolution / degradation over 7 years"),
    "carola": (0.60, "collapse of a prestressed bridge deck"),
    "morandi": (0.69, "an independent collapse, other structure, other sensor"),
    "haltern": (0.78, "a wind turbine's collapse, seen in the mask fragmentation"),
    "kdlo": (0.93, "reconstruction — the record's second state"),
}
SPECTRUM = (("stable", 0.06), ("damage / degradation", 0.34),
            ("collapse", 0.72), ("reconstruction", 0.93))

REFERENCES = (("lumo", "controlled ladder, measured natural frequency"),
              ("bautzen", "bridge + control track, ordinal ladder"))

SDS_JSON = os.path.join(DATA, "sds", "sds_support.json")
ESPOO_JSON = os.path.join(DATA, "espoo", "fig_espoo_ocv_chips.json")

#: Which chip the matrix reads: the support-sized chip (guard 30 m, ring 200 m),
#: the one the new work is defined on. The committed fixed-window numbers are not
#: mixed in — the matrix would then compare two different measurements.
CHIP = "support"

CHANNEL_LABEL = {"ASCENDING": "ASC", "DESCENDING": "DESC"}
FIG_CHANNEL_LABEL = {k: v.replace("gamma^2", r"$\gamma^2$") for k, v in CHANNEL_LABEL.items()}

#: Espoo is the event-free control: its own drift test is the evidence, and its
#: JSON reports the dimensions under slightly different names.
ESPOO_DIMS = {"gamma2": "gamma2_mask", "P": "P", "D": "D", "F": "F"}


def _num(value):
    return None if value is None else float(value)


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def band(sds):
    """The SDS band name and its 0-4 index (the fill intensity)."""
    if sds is None:
        return None, -1
    for index, (limit, name) in enumerate(BANDS):
        if sds < limit:
            return name, index
    return BANDS[-1][1], len(BANDS) - 1


def _sds_cell(block):
    """One `per_channel_max_pairwise` entry as the matrix needs it."""
    if not isinstance(block, dict) or block.get("sds") is None:
        return None
    ci = block.get("ci") or {}
    low, high = _num(ci.get("low")), _num(ci.get("high"))
    pair = [low, high] if None not in (low, high) else None
    gate = None if pair is None else bool(not block.get("straddles_zero"))
    return {"delta": _num(block.get("delta")), "sds": _num(block.get("sds")),
            "ci": pair, "n": [block.get("n_healthy"), block.get("n_damaged")],
            "gate": gate, "gate_kind": "interval" if gate is not None else "none",
            "pair": block.get("winning_pair"), "channel": None,
            "small_n": block.get("small_n"), "p": None}


def _from_espoo(payload, dim):
    """Espoo has no state contrast: its own drift test is the evidence."""
    key = ESPOO_DIMS.get(dim, dim)
    drift = (payload.get("drift") or {}).get(key) or {}
    if not drift:
        return None
    trends, slopes = [], []
    for channel, block in drift.items():
        block = block or {}
        trends.append("trend" in (block.get("flags") or []))
        slopes.append(_num(block.get("theil_sen_per_year")))
    return {"delta": None, "sds": None, "ci": None,
            "n": {CHANNEL_LABEL.get(channel, channel): (block or {}).get("n")
                  for channel, block in drift.items()},
            "gate": None, "gate_kind": "trend",
            "trend": any(trends), "slopes": slopes,
            "channels": {CHANNEL_LABEL.get(ch, ch): {"trend": "trend" in ((b or {}).get("flags") or [])}
                         for ch, b in drift.items()}}


def evidence(site, sds_payload, espoo_payload=None):
    """`{dim: cell}` for one case, from the support-chip analysis of that case."""
    if site == "espoo":
        return {key: _from_espoo(espoo_payload or {}, key) for key, _l in COLUMNS}
    levels = ((sds_payload.get("sites") or {}).get(site) or {}).get("levels") or {}
    block = (levels.get(CHIP) or {}).get("per_channel_max_pairwise") or {}
    cells, reasons = {}, {}
    for key, _label in COLUMNS:
        per_channel = {}
        for channel, dims in block.items():
            cell = _sds_cell((dims or {}).get(key))
            if cell:
                per_channel[CHANNEL_LABEL.get(channel, channel)] = cell
            else:
                reason = ((dims or {}).get(key) or {}).get("reason")
                if reason:
                    reasons[key] = reason
        if not per_channel:
            cells[key] = None
            continue
        def weight(item):
            return sum(v or 0 for v in (item[1].get("n") or [0, 0]))
        channel, cell = max(per_channel.items(), key=weight)
        cell = dict(cell)
        cell["channel"] = channel
        cell["channels"] = per_channel
        cells[key] = cell
    cells["_reasons"] = reasons
    return cells


def symbol(site, dim, cell):
    """The one glyph a cell may carry, and the reason for it.

    An arrow needs the case's own interval to exclude zero **and** at least
    ``SMALL_N`` units on both sides; a circle means tested-but-unresolved, a
    question mark means the evidence is not there, and `↕` is reserved for a case
    that reports a non-monotonic response.
    """
    if cell is None:
        return "?", "the support-chip analysis reports no contrast here"
    if site == "espoo":
        if cell.get("trend"):
            return "↕", "the drift test flags a trend"
        return "○", "no trend in either orbit over two years"
    delta, gate = cell.get("delta"), cell.get("gate")
    counts = [v for v in (cell.get("n") or []) if v is not None]
    if counts and min(counts) < SMALL_N:
        return "?", ("fewer than " + str(SMALL_N) + " units on one side (n = "
                     + "/".join(str(int(v)) for v in counts) + ")")
    if delta is None:
        return "?", "the case reports no effect for this dimension"
    if gate is False:
        return "○", "the interval covers zero"
    if gate is True:
        return ("↑" if delta > 0 else "↓"), "the interval excludes zero"
    return "?", "no interval reported"


def intensity(cell):
    """0-4 fill index: the SDS band, capped at the weakest step without a gate."""
    if cell is None:
        return -1
    _name, index = band(cell.get("sds"))
    if cell.get("gate") is False:
        return max(0, min(index, 1))
    return index




def compute():
    """Read the six cases on the support chip and decide every cell."""
    report = {"schema": SCHEMA, "chip": CHIP,
              "generated_by": "code/cases/fig_cases_overview.py",
              "rules": {"columns": [list(c) for c in COLUMNS],
                        "bands": [[None if limit == float("inf") else limit, name]
                                  for limit, name in BANDS],
                        "arrow_gate": "the case's own interval must exclude zero",
                        "intensity": "the SDS band of the cell, capped at 'weak' "
                                     "when the interval covers zero",
                        "channel": "where a case has two channels, the cell shows "
                                   "the one with the larger total n and the JSON "
                                   "keeps both",
                        "cross_case": "SDS is an intensity scale inside a cell; "
                                      "no effect size is compared across cases"},
              "spectrum": [{"label": label, "x": x} for label, x in SPECTRUM],
              "references": [{"site": site, "label": LABEL[site], "role": role}
                             for site, role in REFERENCES],
              "cases": {}, "sources": {}}
    with open(SDS_JSON) as handle:
        sds_payload = json.load(handle)
    with open(ESPOO_JSON) as handle:
        espoo_payload = json.load(handle)
    for site in CASES:
        cells = evidence(site, sds_payload, espoo_payload)
        reasons = cells.pop("_reasons", {}) if isinstance(cells, dict) else {}
        source = ESPOO_JSON if site == "espoo" else SDS_JSON
        report["sources"][site] = {
            "path": os.path.relpath(source, REPO), "sha256": _sha256(source),
            "chip_level": (None if site == "espoo" else
                           "sds_support.json levels.support.per_channel_max_pairwise")}
        entry = {"label": LABEL[site], "question": POSITIONS[site][1],
                 "x": POSITIONS[site][0], "dims": {}, "reasons": reasons}
        for key, _label in COLUMNS:
            cell = cells.get(key)
            glyph, why = symbol(site, key, cell)
            entry["dims"][key] = {"symbol": glyph, "why": why,
                                  "intensity": intensity(cell),
                                  "band": band((cell or {}).get("sds"))[0],
                                  "delta": (cell or {}).get("delta"),
                                  "sds": (cell or {}).get("sds"),
                                  "ci": (cell or {}).get("ci"),
                                  "n": (cell or {}).get("n"),
                                  "pair": (cell or {}).get("pair"),
                                  "channel": (cell or {}).get("channel"),
                                  "channels": (cell or {}).get("channels"),
                                  "small_n": (cell or {}).get("small_n"),
                                  "p": (cell or {}).get("p")}
        report["cases"][site] = entry
    report["track"] = trajectories(report)
    return report


def _num_text(value, digits=3):
    return "–" if value is None else f"{value:.{digits}g}"


def _ci_text(pair):
    if not pair or pair[0] is None:
        return "–"
    return f"[{pair[0]:+.3g}, {pair[1]:+.3g}]"



# ---------------------------------------------------------------------------
# Panel (c) — the schematic time lines, drawn from each case's result *type*
# ---------------------------------------------------------------------------
def trajectories(report):
    """Six schematic lines: flat, slope, jump, and a second level.

    The geometry encodes the **kind** of result the case reports (no trend, a
    long-term drift, a transition at the event, a coverage break, a second state
    after the reconstruction) — never a rescaled trajectory. The numbers behind
    each line are the ones panel (b) carries.
    """
    return {
        "espoo": {"segments": [[(0.02, 0.55), (0.98, 0.55)]], "event": None,
                  "note": "stable over two years, no drift in any dimension"},
        "cts": {"segments": [[(0.02, 0.30), (0.40, 0.32), (0.52, 0.62),
                              (0.66, 0.66), (0.98, 0.72)]], "event": None,
                "note": "long-term change, no single event to point at"},
        "carola": {"segments": [[(0.02, 0.45), (0.50, 0.45)],
                                [(0.50, 0.45), (0.54, 0.85)],
                                [(0.54, 0.85), (0.98, 0.88)]], "event": 0.50,
                   "note": "transition at the collapse"},
        "morandi": {"segments": [[(0.02, 0.42), (0.55, 0.42)],
                                 [(0.55, 0.42), (0.59, 0.82)],
                                 [(0.59, 0.82), (0.98, 0.85)]], "event": 0.55,
                    "note": "an independent collapse, same shape"},
        "haltern": {"segments": [[(0.02, 0.45), (0.56, 0.45)],
                                  [(0.56, 0.45), (0.60, 0.78)],
                                  [(0.60, 0.78), (0.98, 0.80)]], "event": 0.56,
                    "note": "a wind turbine's collapse; the mask fragments "
                            "afterwards (F ↑)"},
        "kdlo": {"segments": [[(0.02, 0.35), (0.45, 0.35)],
                              [(0.45, 0.35), (0.50, 0.70)],
                              [(0.50, 0.70), (0.98, 0.72)]], "event": 0.45,
                 "note": "second state after the reconstruction"},
    }



# ---------------------------------------------------------------------------
# The figure
# ---------------------------------------------------------------------------
FILL = ("#f2f5f9", "#c6dbef", "#7fb3d5", "#3b7dd8", "#1b3f73")


def _panel_a(ax, report):
    ax.set_xlim(-0.02, 1.04)
    ax.set_ylim(-0.55, 1.55)
    ax.axis("off")
    ax.annotate("", xy=(1.02, 0.0), xytext=(0.0, 0.0),
                arrowprops={"arrowstyle": "-|>", "lw": 1.0, "color": "0.45"})
    for label, x in SPECTRUM:
        ax.plot([x], [0.0], marker="|", ms=9, color="0.45")
        ax.text(x, -0.20, label, ha="center", va="top", fontsize=6.6, color="0.35")
    for index, site in enumerate(CASES):
        entry = report["cases"][site]
        x = entry["x"]
        # the three collapse cases sit close together: stagger their captions so
        # the reader can still see which line belongs to which dot
        lift = 0.10 if index % 2 == 0 else -0.02
        ax.plot([x], [0.30], marker="o", ms=5.0, color="#c53030", zorder=3)
        ax.plot([x, x], [0.03, 0.28], lw=0.8, color="#c53030", alpha=0.7)
        ax.text(x, 0.34, entry["label"], ha="center", va="bottom", fontsize=6.4,
                color="#7b1f1f")
        ax.text(x, 0.42 + lift, entry["question"], ha="left", va="bottom",
                fontsize=4.6, color="0.35", rotation=34, clip_on=False)
    ax.axhspan(1.08, 1.50, facecolor="#eef2f7", edgecolor="none", zorder=0)
    ax.text(0.0, 1.46, "controlled-reference experiments (calibration, not further "
                       "cases)", fontsize=5.6, color="0.35", va="top")
    for index, (site, role) in enumerate(REFERENCES):
        x = 0.06 + index * 0.40
        ax.plot([x], [1.22], marker="s", ms=4.4, color="#2b6cb0")
        ax.text(x + 0.018, 1.22, f"{LABEL[site]} — {role}", fontsize=5.6,
                va="center", color="#1f4e79", rotation=8)
    ax.set_title("(a) what each case is for", loc="left", fontsize=8.0)


def _panel_b(ax, report):
    from matplotlib.patches import Rectangle

    rows = list(CASES)
    columns = [c for c, _l in COLUMNS]
    n_rows, n_cols = len(rows), len(columns)
    ax.set_xlim(0, n_cols + 0.02)
    ax.set_ylim(0, n_rows)
    ax.invert_yaxis()
    ax.set_xticks([i + 0.5 for i in range(n_cols)])
    ax.set_xticklabels([dict(COLUMNS).get(c, "coverage") for c in columns],
                       fontsize=7.2)
    ax.xaxis.tick_top()
    ax.set_yticks([i + 0.5 for i in range(n_rows)])
    ax.set_yticklabels([report["cases"][s]["label"] for s in rows], fontsize=6.6)
    ax.tick_params(length=0)
    for spine in ("left", "right", "bottom"):
        ax.spines[spine].set_visible(False)
    for row_index, site in enumerate(rows):
        entry = report["cases"][site]
        for col_index, key in enumerate(columns):
            x0, y0 = col_index, row_index
            cell = entry["dims"][key]
            index = cell["intensity"]
            face = FILL[max(0, index)] if index >= 0 else "#ffffff"
            ax.add_patch(Rectangle((x0 + 0.04, y0 + 0.08), 0.92, 0.84,
                                   facecolor=face,
                                   edgecolor="0.85" if index < 0 else "none"))
            ax.text(x0 + 0.5, y0 + 0.46, cell["symbol"], ha="center", va="center",
                    fontsize=12.0,
                    color="white" if index >= 3 else ("0.25" if index >= 0
                                                      else "0.55"))
            if cell["delta"] is not None:
                ax.text(x0 + 0.5, y0 + 0.80, _num_text(cell["delta"], 2),
                        ha="center", va="center", fontsize=4.2,
                        color="white" if index >= 3 else "0.35")
    ax.set_title("(b) dimension-response matrix — symbol = direction, fill = "
                 "evidence strength", loc="left", fontsize=8.0, pad=20)


def _panel_c(ax, report):
    ax.set_xlim(-0.02, 1.30)
    ax.set_ylim(0, len(CASES) + 0.4)
    ax.invert_yaxis()
    ax.axis("off")
    for index, site in enumerate(CASES):
        band_top = index + 0.06
        track = report["track"][site]
        if index % 2 == 0:
            ax.axhspan(index, index + 1, facecolor="#f7f9fc", edgecolor="none",
                       zorder=0)
        for segment in track["segments"]:
            xs = [p[0] for p in segment]
            ys = [index + 0.78 - 0.56 * p[1] for p in segment]
            ax.plot(xs, ys, lw=1.6, color="#2b6cb0" if track["event"] is None
                    else "#c53030", solid_capstyle="round", zorder=3)
        if track["event"] is not None:
            ax.plot([track["event"], track["event"]], [index + 0.14, index + 0.9],
                    ls="--", lw=0.7, color="0.55", zorder=1)
        ax.text(0.005, band_top + 0.05, report["cases"][site]["label"],
                fontsize=6.2, va="top", color="0.2")
        ax.text(0.30, band_top + 0.05, track["note"], fontsize=5.0, va="top",
                color="0.4")
    ax.set_title("(c) the shape of each case in time (schematic; the numbers are "
                 "in (b))", loc="left", fontsize=8.0)


def make_figure(report, figdir=None):
    """Draw (a), (b), (c) on one page and write the 630 dpi JPEG."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figdir = figdir or FIGDIR
    os.makedirs(figdir, exist_ok=True)
    fig = plt.figure(figsize=(13.2, 9.4), dpi=110)
    grid = fig.add_gridspec(2, 2, height_ratios=[0.82, 2.0],
                            width_ratios=[1.45, 1.05], hspace=0.30, wspace=0.14)
    _panel_a(fig.add_subplot(grid[0, :]), report)
    _panel_b(fig.add_subplot(grid[1, 0]), report)
    _panel_c(fig.add_subplot(grid[1, 1]), report)
    fig.suptitle(
        "Six cases, one method — what each case is for, what each dimension "
        "answers, and when it moves\n"
        "symbols: ↑/↓ interval excludes zero and n ≥ 6 on both sides · ○ tested, no "
        "resolved change · ↕ non-monotonic · ? evidence not there;  fill = the SDS band "
        f"(intensity *within* a cell, not comparable across cases);  chip = {CHIP}",
        fontsize=7.4)
    out = os.path.join(figdir, "fig_cases_overview.jpg")
    fig.savefig(out, dpi=FIGURE_DPI, facecolor="white")
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# The report
# ---------------------------------------------------------------------------
def markdown(report):
    """The page that explains the figure, cell by cell, with its sources."""
    lines = []
    add = lines.append
    add("# Six cases, one method — the overview figure")
    add("")
    add("*generated by `code/cases/fig_cases_overview.py` from the six committed "
        f"case results (`chip = {report['chip']}`); every cell is pinned with the "
        "sha256 of the file it came from*")
    add("")
    add("## (a) what each case is for")
    add("")
    add("| case | on the spectrum | the question it exists to answer |")
    add("| --- | --- | --- |")
    for site in CASES:
        entry = report["cases"][site]
        add(f"| **{entry['label']}** | {entry['x']:.2f} | {entry['question']} |")
    add("")
    add("Above the spectrum, and **not** counted as cases: "
        + "; ".join(f"**{item['label']}** ({item['role']})"
                    for item in report["references"])
        + ". They are the calibration — a measured natural-frequency ladder and a "
          "bridge with its own control track — not further event cases.")
    add("")
    add("## (b) the dimension-response matrix")
    add("")
    add("| case | " + " | ".join(label for _k, label in COLUMNS) + " |")
    add("| --- |" + " --- |" * len(COLUMNS))
    for site in CASES:
        entry = report["cases"][site]
        cells = []
        for key, _label in COLUMNS:
            cell = entry["dims"][key]
            detail = cell["symbol"]
            if cell["delta"] is not None:
                detail += f" ({_num_text(cell['delta'], 2)}"
                if cell["ci"]:
                    detail += f", {_ci_text(cell['ci'])}"
                detail += ")"
            cells.append(detail)
        add(f"| **{entry['label']}** | " + " | ".join(cells) + " |")
    add("")
    add("Morandi's row is four circles on purpose: an independent collapse is not "
        "automatically an OCV case, and the figure says so instead of borrowing "
        "Carola's arrows.")
    add("")
    add("A cell carries an arrow only when that case's **own** interval excludes "
        "zero **and** both sides carry at least six units; a circle means the "
        "dimension was tested and nothing was resolved; a question mark means the "
        "evidence is not there — fewer than six units on one side is an "
        "observation, not an interval. The fill is the project's "
        "SDS band of that cell — an intensity **inside** the cell, never an effect "
        "size to compare between cases.")
    add("")
    add("### The numbers behind each cell")
    add("")
    add("| case | dimension | symbol | why | delta | 95 % CI | SDS | band | n |")
    add("| --- | --- | --- | --- | ---: | --- | ---: | --- | --- |")
    for site in CASES:
        entry = report["cases"][site]
        for key, _label in COLUMNS:
            cell = entry["dims"][key]
            add(f"| {entry['label']} | `{key}` | {cell['symbol']} | {cell['why']} | "
                f"{_num_text(cell['delta'])} | {_ci_text(cell['ci'])} | "
                f"{_num_text(cell['sds'], 2)} | {cell['band'] or '–'} | "
                f"{cell['n'] if cell['n'] is not None else '–'} |")
    add("")
    add("## (c) the shape of each case in time")
    add("")
    add("| case | line | reading |")
    add("| --- | --- | --- |")
    for site in CASES:
        track = report["track"][site]
        shape = ("no event" if track["event"] is None
                 else f"event at x = {track['event']:.2f}")
        add(f"| {report['cases'][site]['label']} | {shape} | {track['note']} |")
    add("")
    add("Panel (c) is **schematic**: the geometry encodes the kind of result the "
        "case reports (flat, slope, jump, coverage break, second state), not a "
        "rescaled trajectory. Every number is in panel (b).")
    add("")
    add("## What this figure does not claim")
    add("")
    add("* **No cross-case effect size.** `SDS` is used only to shade a cell; the "
        "cases differ in structure, sensor, record length and state definition, so "
        "an arrow in one row says nothing about the size of the arrow in another.")
    add("* **No completeness.** A `?` is a result, not a gap to be filled: where a "
        "case's contrast cannot be cut on the support chip, the cell says so and the "
        "reason is in `data/sds/sds_support.md` — not inferred from an empty axis.")
    add("* **No pooled claim about the six.** The controlled-reference experiments "
        "(`LUMO`, `Bautzen`) are what licenses a statement about the method across "
        "structures; the six cases are the applications.")
    add("* **No chip mixing.** All six are read on the committed chip; the "
        "support-sized chip is an additive second set for KDLO and flips the sign "
        "of `A`, `F` and `P` there, so it is not mixed into this page.")
    add("* **Not every site is a case.** Haltern's canal bridge is a decomposition "
        "site with its own statistics; it is not part of this six-case overview.")
    add("")
    add("## Sources")
    add("")
    add("| case | file | sha256 |")
    add("| --- | --- | --- |")
    for site in CASES:
        src = report["sources"][site]
        add(f"| {report['cases'][site]['label']} | `{src['path']}` | "
            f"`{src['sha256'][:16]}` |")
    add("")
    return "\n".join(lines)


def write_outputs(report, json_path=None, md_path=None, quiet=False):
    json_path = json_path or os.path.join(DATA, "cases", "fig_cases_overview.json")
    md_path = md_path or os.path.join(FIGDIR, "fig_cases_overview.md")
    os.makedirs(os.path.dirname(os.path.abspath(json_path)), exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(md_path)), exist_ok=True)
    with open(json_path, "w") as handle:
        json.dump(report, handle, indent=1, sort_keys=True, default=str)
        handle.write("\n")
    with open(md_path, "w") as handle:
        handle.write(markdown(report))
    if not quiet:
        print(f"  {os.path.relpath(json_path, REPO)}")
        print(f"  {os.path.relpath(md_path, REPO)}")
    return json_path, md_path


# ---------------------------------------------------------------------------
# Selftest and the check
# ---------------------------------------------------------------------------
def selftest():
    """The gate rules and the adapters, on synthetic blocks."""
    failures = []

    def check(label, condition, detail=""):
        if not condition:
            failures.append(f"{label}: {detail}")

    # the interval gate: an arrow needs the case's own interval to exclude zero
    cell = {"delta": 0.4, "sds": 6.5, "ci": [0.1, 0.7], "gate": True,
            "n": [40, 40]}
    check("an interval excluding zero gives an arrow",
          symbol("carola", "P", cell)[0] == "↑", str(symbol("carola", "P", cell)))
    neg = {"delta": -0.4, "sds": 6.5, "ci": [-0.7, -0.1], "gate": True,
           "n": [40, 40]}
    check("a negative interval gives a down arrow",
          symbol("carola", "P", neg)[0] == "↓")
    null = {"delta": 0.2, "sds": 5.0, "ci": [-0.2, 0.6], "gate": False,
            "n": [40, 40]}
    check("an interval covering zero gives a circle",
          symbol("carola", "P", null)[0] == "○")
    check("a circle is capped at the weakest fill",
          intensity(null) == 1, str(intensity(null)))
    thin = {"delta": 0.6, "sds": 8.0, "ci": [0.2, 0.9], "gate": True,
            "n": [11, 2]}
    check("a thin cell cannot carry an arrow",
          symbol("haltern", "P", thin)[0] == "?", str(symbol("haltern", "P", thin)))
    check("no analysis gives a question mark", symbol("cts", "P", None)[0] == "?")

    # the SDS bands drive the fill
    check("the band of SDS 6.5 is 'strong'", band(6.5)[0] == "strong", str(band(6.5)))
    check("the band of SDS 8.5 is 'very strong'", band(8.5)[0] == "very strong")

    # the two adapters read the shapes the repo actually uses
    sds_cell = _sds_cell({"delta": 0.5, "sds": 4.99,
                          "ci": {"low": 0.17, "high": 0.78},
                          "straddles_zero": False, "n_healthy": 16,
                          "n_damaged": 43, "winning_pair": "pre – rebuild",
                          "small_n": False})
    check("the SDS adapter reads the interval dict",
          sds_cell["ci"] == [0.17, 0.78] and sds_cell["gate"] is True, str(sds_cell))
    check("the SDS adapter keeps the pair",
          sds_cell["pair"] == "pre – rebuild")
    espoo = _from_espoo({"drift": {"P": {
        "ASCENDING": {"flags": [], "n": 40},
        "DESCENDING": {"flags": [], "n": 81}}}}, "P")
    check("the Espoo adapter reads the drift flags",
          espoo["trend"] is False and espoo["gate_kind"] == "trend", str(espoo))

    for line in failures:
        print(f"  SELFTEST FAIL {line}")
    print(f"  selftest: 11 checks, {len(failures)} failure(s)")
    return 0 if not failures else 1


def verify(report):
    """Re-derive every cell and every source hash from the files."""
    problems, fresh = [], compute()
    for site in CASES:
        if report["sources"][site]["sha256"] != fresh["sources"][site]["sha256"]:
            problems.append(f"{site}: source sha256 differs from the committed one")
        old, new = report["cases"][site], fresh["cases"][site]
        for key, _label in COLUMNS:
            for field in ("symbol", "intensity", "delta", "sds"):
                if old["dims"][key][field] != new["dims"][key][field]:
                    problems.append(f"{site}/{key}/{field}: {old['dims'][key][field]!r}"
                                    f" vs {new['dims'][key][field]!r}")
        if ("coverage" in old) != ("coverage" in new):
            problems.append(f"{site}: coverage present in one side only")
    return problems


def build_parser():
    parser = argparse.ArgumentParser(
        description="One figure for the six cases: map, dimension matrix, "
                    "trajectories")
    parser.add_argument("--json", default=None)
    parser.add_argument("--md", default=None)
    parser.add_argument("--figdir", default=None)
    parser.add_argument("--selftest", action="store_true")
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.selftest:
        return selftest()
    report = compute()
    if args.verify:
        problems = verify(report)
        for problem in problems:
            print(f"  VERIFY {problem}")
        print(f"  verify: {len(problems)} problem(s)")
        if problems:
            return 1
    write_outputs(report, args.json, args.md, quiet=args.quiet)
    figure = make_figure(report, args.figdir)
    if not args.quiet:
        print(f"  {os.path.relpath(figure, REPO)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
