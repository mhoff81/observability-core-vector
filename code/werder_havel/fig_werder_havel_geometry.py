#!/usr/bin/env python3
"""WERDER HAVEL — the geometry figure: the grid, the deck in pixels, and the fill ratio.

Figure: ``figures/werder_havel/fig_werder_havel_geometry.png`` (+ ``.md``)
Result: ``data/werder_havel/fig_werder_havel_geometry.json``

The record's window is 80 x 80 px per acquisition and its exported decile is 640 px —
both of which are *pixel* statements, and the export carries no pixel-spacing spec. The
derived grid of ``werder_havel_geometry.py`` (metres per pixel off the committed burst
footprints) is what they are priced against; this figure draws that division, so the
reader sees the deck, the window and the grid in the same picture instead of three
numbers in three files.

Panels (2 x 2, 183 mm wide):

  a  the deck, the committed queried ring and one 80 x 80 window per pass geometry, in
     metres: the deck is 704.043 x 35.929 m at bearing 59.59 deg and the windows are
     249..302 m across range against the 548..678 m the deck needs along range — i.e.
     *not one of the four geometries contains the structure*;
  b  metres per pixel per pass geometry (range and azimuth, with the orbit tube's own
     min..max spread): 3.77 x 15.45 m/px on IW1 against 3.12 x 15.30 m/px on IW3, a
     21 % difference in the range spacing;
  c  the deck in pixels: 434/436 px on IW1 and 530/535 px on IW3 against the exported
     640-px decile and the 6400-px window, with the count the external 0.86 claim
     implies (550.4 px) dashed beside them;
  d  the fill ratio: the derived deck share of the decile on each geometry against the
     external claim 0.86 (reproduced to 2.8-3.7 % on IW3, missed by 21 % on IW1).

What this figure is not. It reads no SLC, no mask and no coherence number: every
quantity is the committed geometry of ``data/werder_havel/reference/werder_havel_
geometry.json``, re-derived here by importing ``werder_havel_geometry`` and asserting
that the committed file still matches (``geo.check()``, a stop, not a warning). The 0.86
is the external chain's claim — an input outside this repository, mutable, quoted with
its path in the geometry file and never re-measured from here; the figure prices it, it
does not endorse it. The four reference/outside windows of the record are *not* drawn and
their grid *not* claimed: they measure other windows outside the queried ring (see the
geometry file's ``reference_tracks_note``).

Usage:
  python3 code/werder_havel/fig_werder_havel_geometry.py
  python3 code/werder_havel/fig_werder_havel_geometry.py --figdir /tmp/x --json /tmp/y.json --md /tmp/z.md
  python3 code/werder_havel/fig_werder_havel_geometry.py --no-figures
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Polygon  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import werder_havel_geometry as geo  # noqa: E402

ROOT = os.path.abspath(os.path.join(HERE, os.pardir, os.pardir))
FIGDIR = os.path.join(ROOT, "figures", "werder_havel")
DATA = os.path.join(ROOT, "data", "werder_havel")
NAME = "fig_werder_havel_geometry"
OUT_JSON = os.path.join(DATA, NAME + ".json")
OUT_MD = os.path.join(FIGDIR, NAME + ".md")

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "mathtext.fontset": "stix",
    "font.size": 8.5,
    "axes.titlesize": 9.0,
    "axes.labelsize": 8.5,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7.0,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linewidth": 0.4,
    "figure.dpi": 110,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
})
MM = 1.0 / 25.4
COL_WIDTH = 183 * MM

# The four *covering* pass geometries of the record, and a colour per sub-swath: the
# IW1/IW3 split is the whole content of panels (b)-(d).
TRACKS = list(geo.BURSTS_EXPECTED["tracks"])
SWATH_COLOR = {"IW1": "#1f4e79", "IW3": "#b5651d"}
CLAIM_COLOR = "#7b1e1e"
DECILE_COLOR = "#3f7d3f"


# ---------------------------------------------------------------------------
# The committed geometry, in the picture's own frame
# ---------------------------------------------------------------------------
def ring_local(polygon_lonlat, ref):
    """The committed ring in the anchor's flat metres, east/north of the site point.

    Reuses the module's own constants (``M_PER_DEG_LAT`` / ``M_PER_DEG_LON`` at
    ``LAT0``) so the picture and the metres per pixel are the *same* projection, and
    the anchor's ``site_point`` - the arithmetic mean of the ring's nine distinct
    nodes - as the origin, so the ring is drawn about the point the windows are
    centred on.
    """
    kx = geo.M_PER_DEG_LON * math.cos(math.radians(geo.LAT0))
    ky = geo.M_PER_DEG_LAT
    out = []
    for lon, lat in polygon_lonlat:
        p = ((lon - ref["lon"]) * kx, (lat - ref["lat"]) * ky)
        if not out or (abs(p[0] - out[0][0]) > 1e-9 or abs(p[1] - out[0][1]) > 1e-9):
            out.append(p)
    return out


def rect_corners(cx, cy, along_a, along_b, bearing_deg):
    """Corners of a rectangle centred at ``(cx, cy)``.

    ``along_a`` runs along the bearing - a compass bearing of the *range* axis, so
    ``(sin, cos)`` in east/north metres - and ``along_b`` across it.
    """
    rad = math.radians(bearing_deg)
    ux, uy = math.sin(rad), math.cos(rad)
    vx, vy = uy, -ux
    a, b = along_a / 2.0, along_b / 2.0
    return [(cx + a * ux + b * vx, cy + a * uy + b * vy),
            (cx - a * ux + b * vx, cy - a * uy + b * vy),
            (cx - a * ux - b * vx, cy - a * uy - b * vy),
            (cx + a * ux - b * vx, cy + a * uy - b * vy)]


def _ring_extent(ring, bearing_deg):
    """``(span, width, centroid)`` of the ring along and across a bearing, in metres."""
    rad = math.radians(bearing_deg)
    ux, uy = math.sin(rad), math.cos(rad)
    a = [(x * ux + y * uy) for x, y in ring]
    b = [(y * ux - x * uy) for x, y in ring]
    return {"span_m": max(a) - min(a), "width_m": max(b) - min(b),
            "centroid_m": [sum(p[0] for p in ring) / len(ring),
                           sum(p[1] for p in ring) / len(ring)]}


def compute():
    """Everything the figure and the report need, from the committed geometry file."""
    payload = geo.derive()
    problems, _ = geo.check()
    if problems:
        raise SystemExit("the committed geometry file no longer re-derives, so this "
                         "figure refuses to draw it:\n  - " + "\n  - ".join(problems))
    a = geo.load_aoi()
    sp = a["site_point"]
    ring = ring_local(a["aoi"]["polygon_lonlat"], sp)
    exp = geo.masks.RECORD_EXPECTED
    tracks = [t for t in TRACKS if t in payload["chip_grid"]]
    return {
        "payload": payload,
        "site": payload["site"],
        "asset_id": payload["asset_id"],
        "asset_name": payload["asset_name"],
        "sources": payload["sources"],
        "sp": sp,
        "ring": ring,
        "ring_bbox_m": [min(p[0] for p in ring), min(p[1] for p in ring),
                        max(p[0] for p in ring), max(p[1] for p in ring)],
        # the ring's own extent *along* the fitted span and across it, for the check
        # that the committed ring is the fitted deck and not another polygon
        "ring_extent_m": _ring_extent(ring, payload["deck"]["span_bearing_deg"]),
        "site_point_vs_area_centroid_m": a["site_point_vs_area_centroid_m"],
        "site_point_kind": a["site_point_kind"],
        "tracks": tracks,
        "sub_swath": exp["subswath_per_track"],
        "orbit": exp["orbit_direction_per_track"],
        "rel_orbit": exp["relative_orbit_per_track"],
        "event": {"start": geo.masks.EVENT_START, "end": geo.masks.EVENT_END,
                  "source": geo.masks.EVENT_SOURCE},
        "record": {"rows": exp["rows"], "columns": exp["columns"],
                   "n_dates": exp["n_dates"], "first_date": exp["first_date"],
                   "last_date": exp["last_date"],
                   "state_windows": exp["state_windows"]},
        "window": payload["window"],
        "deck": payload["deck"],
        "grid": payload["chip_grid"],
        "grid_note": payload["chip_grid_note"],
        "vs": payload["deck_vs_pixel_grid"],
        "vs_note": payload["deck_vs_pixel_grid_note"],
        "win": payload["deck_vs_window"],
        "win_note": payload["deck_vs_window_note"],
        "claim": payload["fill_ratio_claim"],
        "corr": payload["area_corrections"],
        "ref_note": payload["reference_tracks_note"],
        "notes": payload["notes"],
    }


# ---------------------------------------------------------------------------
# The panels
# ---------------------------------------------------------------------------
def save_fig(fig, name, out_dir=None):
    """Write one 600 dpi PNG (the binaries stay out of git)."""
    out_dir = out_dir or FIGDIR
    os.makedirs(out_dir, exist_ok=True)
    png = os.path.join(out_dir, name + ".png")
    fig.savefig(png)
    plt.close(fig)
    print(f"  {os.path.relpath(png)}")
    return png


def _deck_patch(ax, d, **kw):
    ax.add_patch(Polygon(rect_corners(0.0, 0.0, d["span_m"], d["width_m"],
                                      d["span_bearing_deg"]),
                         closed=True, fill=False, **kw))


def panel_a(ax, res):
    """The deck, the queried ring and one 80 x 80 window per pass geometry, in metres."""
    d, grid, win = res["deck"], res["grid"], res["win"]
    ax.add_patch(Polygon(res["ring"], closed=True, fc="#d6d6d6", ec="#222222", lw=0.9,
                         zorder=4, label="the committed queried ring"))
    _deck_patch(ax, d, ec="#222222", lw=0.7, ls=(0, (3, 2)), zorder=5,
                label=f"the fitted rectangle {d['span_m']:.3f} x {d['width_m']:.3f} m")
    for t in res["tracks"]:
        s = grid[t]
        ax.add_patch(Polygon(rect_corners(0.0, 0.0, win[t]["window_range_m"],
                                          win[t]["window_azimuth_m"],
                                          s["range_axis_bearing_deg"]),
                             closed=True, fill=False, lw=1.0, zorder=6,
                             ec=SWATH_COLOR[s["sub_swath"][0]],
                             label=f"80 x 80 window, {t} ({s['sub_swath'][0]})"))
    ax.plot([0.0], [0.0], marker="+", ms=5, mew=0.9, color="#000000", zorder=7,
            label="site point (window centre)")
    half = 1.06 * max(v["window_azimuth_m"] for v in win.values()) / 2.0
    ax.set_xlim(-half, half)
    ax.set_ylim(-half, half)
    ax.set_aspect("equal")
    ax.set_xlabel("east of the site point [m]")
    ax.set_ylabel("north of the site point [m]")
    nr = sorted(v["deck_span_over_window_range"] for v in win.values())
    ax.text(0.02, 0.02,
            f"window {min(v['window_range_m'] for v in win.values()):.0f}-"
            f"{max(v['window_range_m'] for v in win.values()):.0f} m across range\n"
            f"deck needs {nr[0]:.2f}-{nr[-1]:.2f} of that range span\n"
            "the deck is 704 m long by 36 m wide, the window 249-302 m\n"
            "across range: the window samples the deck's middle",
            transform=ax.transAxes, fontsize=6.6, va="bottom", ha="left",
            bbox=dict(fc="white", ec="#999999", lw=0.4, alpha=0.9, pad=2.0))
    ax.set_title("(a) the deck, the ring and the record's window, in metres", loc="left")
    ax.legend(loc="upper left", fontsize=6.0, framealpha=0.9)


def panel_b(ax, res):
    """Metres of ground per pixel, per pass geometry, with the orbit tube's spread."""
    grid = res["grid"]
    x = np.arange(len(res["tracks"]))
    bw = 0.34
    for i, t in enumerate(res["tracks"]):
        s = grid[t]
        col = SWATH_COLOR[s["sub_swath"][0]]
        for off, key, lo_key, alpha in (
                (-bw / 2.0, "range_m_per_px", "range_m_per_px_spread", 1.0),
                (bw / 2.0, "azimuth_m_per_px", "azimuth_m_per_px_spread", 0.45)):
            v = s[key]
            lo, hi = s[lo_key]
            ax.bar(i + off, v, width=bw, color=col, alpha=alpha, edgecolor="k", lw=0.3)
            ax.errorbar([i + off], [v], yerr=[[v - lo], [hi - v]], fmt="none",
                        ecolor="k", elinewidth=0.6, capsize=1.8)
        ax.text(i, max(s["azimuth_m_per_px_spread"][1], s["range_m_per_px_spread"][1]) + 0.5,
                f"{s['px_area_m2']:.1f} m2/px", ha="center", fontsize=6.8)
    ax.set_xticks(x)
    ax.set_xticklabels([t.replace("_", " ") + f"\n{grid[t]['sub_swath'][0]}"
                        f" x{grid[t]['n_frames']}" for t in res["tracks"]], fontsize=6.6)
    ax.set_ylim(0, 19.5)
    ax.set_ylabel("metres of ground per pixel")
    ax.set_title("(b) the grid of the record, off the committed footprints", loc="left")
    ax.legend(handles=[
        Line2D([], [], color=SWATH_COLOR["IW1"], lw=4, label="IW1"),
        Line2D([], [], color=SWATH_COLOR["IW3"], lw=4, label="IW3"),
        Line2D([], [], color="#555555", lw=4, label="range spacing (solid)"),
        Line2D([], [], color="#555555", lw=4, alpha=0.45, label="azimuth spacing (pale)"),
    ], loc="upper right", fontsize=6.2, framealpha=0.9)




def panel_c(ax, res):
    """The deck in pixels against the exported decile and the claim's own count."""
    vs, claim = res["vs"], res["claim"]
    masked = res["window"]["masked_px"]
    wpx = res["window"]["px"]
    implied_px = res["deck"]["area_m2"] / (claim["claim"] * masked)
    x = np.arange(len(res["tracks"]))
    cols = [SWATH_COLOR[res["grid"][t]["sub_swath"][0]] for t in res["tracks"]]
    px = [vs[t]["deck_area_px"] for t in res["tracks"]]
    ax.bar(x, px, width=0.5, color=cols, edgecolor="k", lw=0.4)
    for i, t in enumerate(res["tracks"]):
        ax.text(i, px[i] + 16,
                f"{px[i]:.1f} px\n= {vs[t]['deck_frac_of_masked_decile']:.3f} of the decile",
                ha="center", fontsize=6.6)
    ax.axhline(masked, color=DECILE_COLOR, lw=1.0, ls="--")
    ax.text(len(x) - 0.45, masked + 8, f"the exported decile, {masked} px",
            color=DECILE_COLOR, ha="right", fontsize=6.8)
    ax.axhline(implied_px, color=CLAIM_COLOR, lw=0.9, ls=":")
    ax.text(-0.55, implied_px + 8, f"the count the {claim['claim']} claim implies, "
            f"{implied_px:.1f} px", color=CLAIM_COLOR, ha="left", fontsize=6.6)
    ax.set_xticks(x)
    ax.set_xticklabels([t.replace("_", " ") for t in res["tracks"]], fontsize=6.6)
    ax.set_xlim(-0.6, len(x) - 0.4)
    ax.set_ylim(0, 760)
    ax.set_ylabel("pixels of ground the deck covers")
    ax.set_title("(c) the deck in pixels (deck_area_px = area / px area)", loc="left")
    ax.text(0.02, 0.03,
            f"the window is {wpx} px; the deck is "
            f"{min(v['deck_frac_of_window'] for v in vs.values()) * 100:.1f}-"
            f"{max(v['deck_frac_of_window'] for v in vs.values()) * 100:.1f} % of it",
            transform=ax.transAxes, fontsize=6.6, va="bottom",
            bbox=dict(fc="white", ec="#999999", lw=0.4, alpha=0.9, pad=2.0))


def panel_d(ax, res):
    """The derived fill ratio against the external 0.86 claim."""
    vs, claim, per = res["vs"], res["claim"], res["claim"]["per_track"]
    x = np.arange(len(res["tracks"]))
    cols = [SWATH_COLOR[res["grid"][t]["sub_swath"][0]] for t in res["tracks"]]
    der = [vs[t]["deck_frac_of_masked_decile"] for t in res["tracks"]]
    ax.bar(x, der, width=0.5, color=cols, edgecolor="k", lw=0.4)
    ax.axhline(claim["claim"], color=CLAIM_COLOR, lw=1.1, ls="--")
    ax.text(len(x) - 0.45, claim["claim"] + 0.008,
            f"the external claim {claim['claim']:.2f}", color=CLAIM_COLOR,
            ha="right", fontsize=6.8)
    for i, t in enumerate(res["tracks"]):
        ax.text(i, der[i] + 0.008,
                f"{der[i]:.3f}\n{per[t]['derived_minus_claim_pct']:+.2f} %",
                ha="center", fontsize=6.6)
    ax.set_xticks(x)
    ax.set_xticklabels([t.replace("_", " ") + f"\n{res['grid'][t]['sub_swath'][0]}"
                        for t in res["tracks"]], fontsize=6.6)
    ax.set_xlim(-0.6, len(x) - 0.4)
    ax.set_ylim(0, 0.98)
    ax.set_ylabel("deck area / the 640-px decile")
    ax.set_title("(d) the fill ratio: derived against the claim", loc="left")
    g3 = sorted(per[t]["derived_minus_claim_pct"]
                for t in res["tracks"] if res["grid"][t]["sub_swath"][0] == "IW3")
    g1 = sorted(per[t]["derived_minus_claim_pct"]
                for t in res["tracks"] if res["grid"][t]["sub_swath"][0] == "IW1")
    ax.text(0.02, 0.03,
            f"IW3 {g3[0]:+.1f}..{g3[-1]:+.1f} % of the claim, "
            f"IW1 {g1[0]:+.1f}..{g1[-1]:+.1f} %:\n"
            "the claim is an IW3 reading, priced, never adopted",
            transform=ax.transAxes, fontsize=6.6, va="bottom",
            bbox=dict(fc="white", ec="#999999", lw=0.4, alpha=0.9, pad=2.0))


def fig_geometry(res, out_dir=None):
    """The 2 x 2 plate: the deck, the grid, the deck in pixels and the claim."""
    fig, axes = plt.subplots(2, 2, figsize=(COL_WIDTH, COL_WIDTH * 1.02))
    panel_a(axes[0][0], res)
    panel_b(axes[0][1], res)
    panel_c(axes[1][0], res)
    panel_d(axes[1][1], res)
    d = res["deck"]
    fig.suptitle(f"{res['asset_name']} - the geometry of the record: a "
                 f"{d['span_m']:.1f} m x {d['width_m']:.1f} m deck on the committed "
                 f"burst grid", fontsize=9.0)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.98))
    return save_fig(fig, NAME, out_dir)


# ---------------------------------------------------------------------------
# The teeth: what the figure asserts, and never draws
# ---------------------------------------------------------------------------
def self_check(res):
    """Every claim the figure makes, as a boolean - the drawing is downstream of this.

    These are not plot assertions: each one is a statement the plate's caption, axes
    or annotations make, checked against the committed geometry (which ``compute`` has
    already re-derived) and, where the arithmetic is the content, against the two
    factors it was built from. A failure is a stop.
    """
    grid, vs, win = res["grid"], res["vs"], res["win"]
    claim, w, d, corr = res["claim"], res["window"], res["deck"], res["corr"]
    tr = res["tracks"]
    checks = []

    def add(name, ok, detail=""):
        checks.append({"check": name, "ok": bool(ok), "detail": detail})

    problems, _ = geo.check()
    add("the committed geometry re-derives field by field", not problems,
        f"geo.check(): {len(problems)} problem(s)" if problems else "geo.check(): no problem")

    add("the window is 80 x 80 = 6400 px and its decile 640 px = 10 % of it",
        w["side_px"] == 80 and w["px"] == w["side_px"] ** 2 == 6400
        and w["masked_px"] == 640 and abs(w["masked_frac_of_window"] - 0.1) < 1e-9,
        f"{w['side_px']} x {w['side_px']} = {w['px']} px, decile {w['masked_px']} px")

    add("no pass geometry of the record contains the deck in one window",
        all(win[t]["deck_contained"] is False for t in tr),
        "deck spans " + ", ".join(f"{win[t]['deck_span_over_window_range']:.2f}"
                                  for t in tr) + " x the window's range extent")

    sides = [win[t]["square_side_px_to_contain"] for t in tr]
    add("containing the deck needs 147..220 px per side", 147 <= min(sides) and max(sides) <= 220,
        f"{min(sides)}..{max(sides)} px per side, the window is {w['side_px']}")

    dev = max(abs(grid[t]["px_area_m2"] - grid[t]["range_m_per_px"] * grid[t]["azimuth_m_per_px"])
              for t in tr)
    add("px_area_m2 is range x azimuth on every geometry", dev <= 0.02,
        f"max |diff| {dev:.4f} m2 over the four geometries")

    dev = max(abs(vs[t]["deck_area_px"] - d["area_m2"] / grid[t]["px_area_m2"]) for t in tr)
    add("deck_area_px = deck area / px area on every geometry", dev <= 0.05,
        f"max |diff| {dev:.4f} px")

    r3 = [claim["per_track"][t]["derived_minus_claim_pct"] for t in tr
          if grid[t]["sub_swath"][0] == "IW3"]
    r1 = [claim["per_track"][t]["derived_minus_claim_pct"] for t in tr
          if grid[t]["sub_swath"][0] == "IW1"]
    add("the external claim is reproduced within 4 % on IW3",
        max(abs(v) for v in r3) <= 4.0, ", ".join(f"{v:+.2f} %" for v in r3))
    add("the external claim is missed by more than 15 % on IW1",
        min(abs(v) for v in r1) >= 15.0, ", ".join(f"{v:+.2f} %" for v in r1))

    add("IW1 resolves the deck coarser than IW3 on both orbits",
        all(grid[a]["px_area_m2"] > grid[b]["px_area_m2"]
            for a, b in (("asc_146_iw1", "asc_44_iw3"), ("desc_168_iw1", "desc_95_iw3"))),
        f"{min(grid[t]['px_area_m2'] for t in tr):.1f}.."
        f"{max(grid[t]['px_area_m2'] for t in tr):.1f} m2/px")

    add("the claim's implied grid is finer than every derived grid",
        claim["implied_px_area_m2"] < min(grid[t]["px_area_m2"] for t in tr),
        f"claim implies {claim['implied_px_area_m2']:.2f} m2/px, derived "
        f"{min(grid[t]['px_area_m2'] for t in tr):.2f}.."
        f"{max(grid[t]['px_area_m2'] for t in tr):.2f} m2/px")

    wsp = [100.0 * abs(1.0 - vs[t]["deck_frac_of_masked_decile_wgs84_area"]
                       / vs[t]["deck_frac_of_masked_decile"]) for t in tr]
    add("the WGS84 reading of the deck's area moves the ratio 0.87 %",
        all(0.80 <= v <= 0.95 for v in wsp),
        f"{min(wsp):.2f}..{max(wsp):.2f} %")

    sl = list(corr["ground_range_not_slant"]["deck_frac_if_slant_range"].values())
    add("slant range would push every ratio above 1", all(v > 1.0 for v in sl),
        f"{min(sl):.3f}..{max(sl):.3f} on the four geometries")

    dd = [corr["complex_sample_not_two_numbers"]["deck_frac_if_decile_were_1280"][t]
          for t in tr]
    add("counting the payload's numbers as pixels would halve every ratio",
        all(abs(2.0 * dd[i] - vs[tr[i]]["deck_frac_of_masked_decile"]) <= 0.01
            for i in range(len(tr))),
        f"{min(dd):.3f}..{max(dd):.3f} against the committed 0.68..0.84")

    cen = res["ring_extent_m"]["centroid_m"]
    add("the ring's own nodes are centred on the site point the window is centred on",
        max(abs(cen[0]), abs(cen[1])) <= 1.0, f"centroid {cen[0]:.2f} / {cen[1]:.2f} m")

    ext = res["ring_extent_m"]
    add("the queried ring is this deck, not another polygon",
        abs(ext["span_m"] - d["span_m"]) <= 0.01 * d["span_m"]
        and abs(ext["width_m"] - d["width_m"]) <= 0.10 * d["width_m"],
        f"ring {ext['span_m']:.2f} x {ext['width_m']:.2f} m along that bearing against the "
        f"fitted {d['span_m']:.3f} x {d['width_m']:.3f} m")
    return checks



# ---------------------------------------------------------------------------
# The report, the JSON and the command line
# ---------------------------------------------------------------------------
def _tbl(rows, head):
    out = ["| " + " | ".join(head) + " |", "|" + "|".join(["---"] * len(head)) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def report_tables(res):
    """The four tables of the ``.md``: the grid, the window, the pixels, the corrections."""
    grid, vs, win, claim = res["grid"], res["vs"], res["win"], res["claim"]
    corr, tr = res["corr"], res["tracks"]

    grid_rows = [[t, grid[t]["sub_swath"][0], res["orbit"][t][:3].lower(),
                  "{:.3f} ({:.3f}-{:.3f})".format(grid[t]["range_m_per_px"],
                                                 *grid[t]["range_m_per_px_spread"]),
                  "{:.3f} ({:.3f}-{:.3f})".format(grid[t]["azimuth_m_per_px"],
                                                  *grid[t]["azimuth_m_per_px_spread"]),
                  f"{grid[t]['px_area_m2']:.2f}", grid[t]["n_frames"],
                  len(grid[t]["burst_ids"])] for t in tr]

    win_rows = [[t, "{:.0f} x {:.0f}".format(win[t]["window_range_m"],
                                             win[t]["window_azimuth_m"]),
                 "{:.2f} x {:.2f}".format(win[t]["deck_span_over_window_range"],
                                          win[t]["deck_span_over_window_azimuth"]),
                 "{} / {}".format(win[t]["side_px_to_contain_range"],
                                  win[t]["side_px_to_contain_azimuth"]),
                 win[t]["square_side_px_to_contain"],
                 "no" if not win[t]["deck_contained"] else "YES"] for t in tr]

    px_rows = [[t, f"{vs[t]['deck_area_px']:.1f}",
                f"{vs[t]['deck_frac_of_masked_decile']:.3f}",
                f"{vs[t]['deck_frac_of_masked_decile_wgs84_area']:.3f}",
                f"{100.0 * vs[t]['deck_frac_of_window']:.2f} %",
                f"{claim['per_track'][t]['derived_minus_claim_pct']:+.2f} %"] for t in tr]

    cs = corr["complex_sample_not_two_numbers"]["deck_frac_if_decile_were_1280"]
    sl = corr["ground_range_not_slant"]["deck_frac_if_slant_range"]
    wg = corr["flat_constants_not_wgs84"]["deck_frac_if_wgs84_area"]
    corr_rows = [
        ["12 800 payload numbers read as 12 800 px (the window twice)",
         " / ".join(f"{cs[t]:.3f}" for t in tr)],
        ["slant range x azimuth instead of ground",
         " / ".join(f"{sl[t]:.3f}" for t in tr)],
        ["the deck's WGS84 area instead of the flat constants",
         " / ".join(f"{wg[t]:.3f}" for t in tr)],
        ["**the committed reading**",
         " / ".join(f"**{vs[t]['deck_frac_of_masked_decile']:.3f}**" for t in tr)],
    ]
    return grid_rows, win_rows, px_rows, corr_rows


def report(res):
    """The ``.md`` beside the plate: the same numbers in text, with their caveats."""
    grid, vs, win, claim = res["grid"], res["vs"], res["win"], res["claim"]
    w, d, rec, ev = res["window"], res["deck"], res["record"], res["event"]
    tr = res["tracks"]
    grid_rows, win_rows, px_rows, corr_rows = res["tables"]

    rng_lo = min(grid[t]["range_m_per_px"] for t in tr)
    rng_hi = max(grid[t]["range_m_per_px"] for t in tr)
    azi_lo = min(grid[t]["azimuth_m_per_px"] for t in tr)
    azi_hi = max(grid[t]["azimuth_m_per_px"] for t in tr)
    pa_lo = min(grid[t]["px_area_m2"] for t in tr)
    pa_hi = max(grid[t]["px_area_m2"] for t in tr)
    dpx = [vs[t]["deck_area_px"] for t in tr]
    w_lo = min(win[t]["window_range_m"] for t in tr)
    w_hi = max(win[t]["window_range_m"] for t in tr)
    a_lo = min(win[t]["window_azimuth_m"] for t in tr)
    a_hi = max(win[t]["window_azimuth_m"] for t in tr)
    ov_lo = min(win[t]["deck_span_over_window_range"] for t in tr)
    ov_hi = max(win[t]["deck_span_over_window_range"] for t in tr)
    sd = [win[t]["square_side_px_to_contain"] for t in tr]
    c_lo = min(vs[t]["deck_frac_of_masked_decile"] for t in tr)
    c_hi = max(vs[t]["deck_frac_of_masked_decile"] for t in tr)
    fr_lo = 100.0 * min(vs[t]["deck_frac_of_window"] for t in tr)
    fr_hi = 100.0 * max(vs[t]["deck_frac_of_window"] for t in tr)
    r3 = " and ".join(f"{claim['per_track'][t]['derived_minus_claim_pct']:+.2f} % ({t})"
                      for t in tr if grid[t]["sub_swath"][0] == "IW3")
    r1 = " and ".join(f"{claim['per_track'][t]['derived_minus_claim_pct']:+.2f} % ({t})"
                      for t in tr if grid[t]["sub_swath"][0] == "IW1")

    L = [f"# {NAME}\n"]
    L.append(f"**Site.** {res['site']} - the committed OSM outline of `{res['asset_name']}` "
             f"(`{res['asset_id']}`), {d['span_m']:.3f} m x {d['width_m']:.3f} m at bearing "
             f"{d['span_bearing_deg']:.2f} deg, {d['area_source']}.\n")
    L.append("**Record.** {rows} rows x {cols} columns, {nd} acquisition dates "
             "{f}..{l}; state windows {sw}; the documented change {s}..{e} ({src}).\n".format(
                 rows=rec["rows"], cols=rec["columns"], nd=rec["n_dates"], f=rec["first_date"],
                 l=rec["last_date"],
                 sw=", ".join(f"{k} {v}" for k, v in rec["state_windows"].items()),
                 s=ev["start"], e=ev["end"], src=ev["source"]))
    L.append(f"**Figure.** `figures/werder_havel/{NAME}.png`, 2 x 2 panels at 183 mm, "
             f"600 dpi. **Result.** `data/werder_havel/{NAME}.json`.\n")
    L.append(f"**In one line.** The record's window is 80 x 80 = {w['px']} px and its "
             f"exported decile {w['masked_px']} px; the four covering geometries put a pixel "
             f"at {pa_lo:.2f}-{pa_hi:.2f} m2 of ground ({rng_lo:.2f}-{rng_hi:.2f} m in range "
             f"by {azi_lo:.2f}-{azi_hi:.2f} m in azimuth), so the deck is {min(dpx):.0f}-"
             f"{max(dpx):.0f} px of the decile ({c_lo:.2f}-{c_hi:.2f} of it) - and the window, "
             f"{w_lo:.0f}-{w_hi:.0f} m wide in range by {a_lo:.0f}-{a_hi:.0f} m in azimuth, "
             f"holds {ov_lo:.2f}-{ov_hi:.2f} x the range span the deck needs: not one of the "
             f"four covers the structure, and the deck is only {fr_lo:.1f}-{fr_hi:.1f} % of "
             f"the window.\n")
    L.append("## What this figure is\n")
    L.append("A division, drawn. `werder_havel_geometry.py` derives metres per pixel from the "
             "committed burst footprints - the family's off-the-footprint method: the longest "
             "footprint edge is the range axis, collinear pieces are summed, and "
             "`SamplesPerBurst` / `LinesPerBurst` divide the extents, with no product-spacing "
             "lookup. This figure reads that result, re-derives it and stops if the committed "
             "file no longer matches (`geo.check()`), then draws the four quantities a reader "
             "has to hold at once: the deck in metres, the grid in metres per pixel, the deck "
             "as a count of pixels, and the fill ratio the external chain claims for it.\n")
    L.append("## What this figure is not\n")
    L.append(f"It is not a measurement: no SLC, no mask and no coherence value is read here, "
             f"and the only external input is the {claim['claim']:.2f} claim, quoted with its "
             "path, priced on every geometry and never adopted as this package's number. It is "
             f"not the echo mask either - the {w['masked_px']} px is the export's own "
             "`coherence_masked_pixels` (the brightest decile of the window by definition), "
             "not the peak/median rule of `werder_havel_ocv_masks.py`. And it claims nothing "
             "for the record's four reference/outside chains: they measure other windows "
             "outside the queried ring, so their own metres per pixel is not derivable from "
             "the committed footprints (see `reference_tracks_note`).\n")
    L.append("## The four panels\n")
    L.append(f"**a - the deck, the ring and the window, in metres.** The grey ring is the "
             f"committed outline (`werder_havel_aoi.json`, {len(res['ring'])} distinct nodes, "
             f"drawn about the site point - the mean of those nodes, the point the record's "
             f"one window per acquisition is taken about, "
             f"`werder_havel_ocv_masks.WINDOW_NOTE`; that mean sits "
             f"{res['site_point_vs_area_centroid_m']:.1f} m off the ring's own area "
             f"centroid). The dashed rectangle is the fitted {d['span_m']:.3f} x "
             f"{d['width_m']:.3f} m deck every number in this file uses; the four coloured "
             f"rectangles are the same 80 x 80 window on each covering geometry's own axes "
             f"(its range axis at "
             f"{min(grid[t]['range_axis_bearing_deg'] for t in tr):.2f}-"
             f"{max(grid[t]['range_axis_bearing_deg'] for t in tr):.2f} deg of bearing, the "
             f"longest footprint edge). The window is {w_lo:.0f}-{w_hi:.0f} m across range "
             f"against the {d['span_m']:.1f} m deck: it samples the deck's middle, and "
             f"containing the bridge would need {min(sd)}-{max(sd)} px per side, not "
             f"{w['side_px']}. The ring's own extent along that bearing is "
             f"{res['ring_extent_m']['span_m']:.2f} x "
             f"{res['ring_extent_m']['width_m']:.2f} m against the fitted "
             f"{d['span_m']:.3f} x {d['width_m']:.3f} m - the map outline is not a rectangle, "
             f"so the fitted width is the audit's and the ring keeps its own.\n")
    L.append(f"**b - the grid of the record.** Per covering geometry, the range spacing (solid) "
             f"and the azimuth spacing (pale), the orbit tube's own min..max spread as "
             f"whiskers: {rng_lo:.2f}-{rng_hi:.2f} m by {azi_lo:.2f}-{azi_hi:.2f} m, i.e. "
             f"{pa_lo:.2f}-{pa_hi:.2f} m2 per pixel. The IW3 geometries resolve range "
             f"{100.0 * (1.0 - min(grid[t]['range_m_per_px'] for t in tr) / max(grid[t]['range_m_per_px'] for t in tr)):.0f} % "
             f"finer than the IW1 ones; azimuth is within "
             f"{100.0 * (azi_hi - azi_lo) / azi_hi:.1f} % across all four.\n")
    L.append(f"**c - the deck in pixels.** `deck_area_px = area_m2 / (range_m_per_px x "
             f"azimuth_m_per_px)`: {min(dpx):.1f}-{max(dpx):.1f} px against the exported "
             f"{w['masked_px']}-px decile (green) and against the {claim['claim']:.2f} claim's "
             f"own {claim['implied_px_area_m2']:.3f} m2/px, i.e. "
             f"{d['area_m2'] / (claim['claim'] * w['masked_px']):.1f} px (dotted). The window "
             f"is {w['px']} px; the deck is {fr_lo:.1f}-{fr_hi:.1f} % of it.\n")
    L.append(f"**d - the fill ratio.** The derived deck share of the decile on each geometry "
             f"against the external {claim['claim']:.2f}: below it by {r3} on the IW3 "
             f"geometries and by {r1} on the IW1 ones. The claim is an IW3 reading, "
             f"reproduced to within "
             f"{max(abs(claim['per_track'][t]['derived_minus_claim_pct']) for t in tr if grid[t]['sub_swath'][0] == 'IW3'):.2f} % "
             f"there and missed by more than "
             f"{min(abs(claim['per_track'][t]['derived_minus_claim_pct']) for t in tr if grid[t]['sub_swath'][0] == 'IW1'):.0f} % "
             f"on the two IW1 ones: the same deck and the same decile, read on a 21 % "
             f"coarser range spacing.\n")
    L.append("## The numbers the panels draw\n")
    L.append("**The grid (panel b), off the committed footprints.**\n")
    L.append(_tbl(grid_rows, ["track", "swath", "orbit", "range m/px (spread)",
                              "azimuth m/px (spread)", "m2/px", "frames", "bursts"]))
    L.append("\n**The window against the deck (panel a).** `deck_span_over_window_*` is how "
             "many window spans the deck needs; `side_px_to_contain_*` the square window that "
             "*would* hold it.\n")
    L.append(_tbl(win_rows, ["track", "window m (range x azimuth)",
                             "deck spans (range x azimuth)", "side px to contain (r / a)",
                             "square px to contain", "contains?"]))
    L.append(f"\n**The deck in pixels (panels c and d).** The decile is {w['masked_px']} px = "
             f"{w['masked_frac_of_window']:.0%} of the {w['px']}-px window by definition; the "
             f"WGS84 column is the flat-vs-WGS84 sensitivity carried beside the committed "
             f"reading, never instead of it.\n")
    L.append(_tbl(px_rows, ["track", "deck px", "deck px / decile",
                            "same, on the WGS84 area", "deck / window", "vs claim"]))
    L.append("\n**The three corrections the ratio can be mis-read through** (committed on "
             "every geometry; the last row is panel d's own arithmetic).\n")
    L.append(_tbl(corr_rows, ["reading"] + list(tr)))
    L.append(f"\nThe claim's implied grid is {claim['implied_px_area_m2']:.3f} m2/px, "
             f"{min(grid[t]['px_area_m2'] for t in tr) - claim['implied_px_area_m2']:.3f} m2 "
             f"finer than the finest geometry of the record. {claim['implied_grid_note']}\n")
    L.append(f"The claim is on the record as a quote, with its own verdict from the geometry "
             f"audit: {claim['verdict']}\n")
    L.append("## Caveats\n")
    L.append(f"- **The claim is a mutable external input.** `{claim['source']}` - the quote is "
             f"committed here, the file is not, and nothing in this package re-measures it.\n")
    L.append("- **One grid per geometry, not per acquisition.** The spreading in panel b is the "
             "orbit tube's own; the committed grid is its median, so the ratio moves by that "
             "spread and not by the median alone.\n")
    L.append("- **No pixel layer is claimed.** The record's windows carry no spacing spec. The "
             "derived grid is the covering layer's, at the site; the four reference/outside "
             "chains are excluded for that reason (`reference_tracks_note`).\n")
    L.append("- **This package's echo mask is a different object.** "
             "`werder_havel_ocv_masks.echo_mask` (peak/median, with 640 px per mask as its own "
             "definition) is not the export's decile column. The figure prices the export's "
             "column because that is what the claim is about.\n")
    L.append("## Reproduce\n")
    L.append("```sh\npython3 code/werder_havel/werder_havel_geometry.py --check\n"
             "python3 code/werder_havel/fig_werder_havel_geometry.py\n```\n")
    L.append("The figure writes `figures/werder_havel/fig_werder_havel_geometry.png` and "
             "`data/werder_havel/fig_werder_havel_geometry.json`. It reads the committed "
             "`data/werder_havel/reference/werder_havel_geometry.json` - re-derived, not "
             "copied - and `data/werder_havel/reference/werder_havel_aoi.json` for the ring of "
             "panel a. No network, no SLC, no window payload, offline.\n")

    return "\n".join(L)


# ---------------------------------------------------------------------------
# The result JSON and the command line
# ---------------------------------------------------------------------------
def build_json(res, checks):
    """The committed geometry, plus what this figure drew and asserted over it."""
    out = dict(res["payload"])
    passed = sum(1 for c in checks if c["ok"])
    out["figure"] = {
        "name": NAME,
        "png": f"figures/werder_havel/{NAME}.png",
        "md": f"figures/werder_havel/{NAME}.md",
        "width_mm": round(COL_WIDTH / MM, 1),
        "dpi": 600,
        "generated": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "panels": [
            "a the deck, the committed ring and one 80x80 window per pass geometry, in metres",
            "b metres per pixel per pass geometry, with the orbit tube's own spread",
            "c the deck in pixels (deck_area_px = area / px area) against the 640-px decile "
            "and the count the claim implies",
            "d the fill ratio of panel c against the external claim",
        ],
        "reads": [
            "data/werder_havel/reference/werder_havel_geometry.json - re-derived by importing "
            "werder_havel_geometry, not copied, and geo.check() must pass",
            "data/werder_havel/reference/werder_havel_aoi.json - the ring of panel a only",
        ],
        "ring_extent_m": res["ring_extent_m"],
        "ring_note": (f"the ring drawn in panel a: {len(res['ring'])} distinct nodes, "
                      f"{res['site_point_kind']} as its origin/centre - the point the "
                      f"record's one window per acquisition is taken about"),
        "self_check": {"n_checks": len(checks), "n_passed": passed, "checks": checks},
    }
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--figdir", default=FIGDIR, help="directory for the PNG")
    ap.add_argument("--json", default=OUT_JSON, help="path of the result JSON")
    ap.add_argument("--md", default=OUT_MD, help="path of the report")
    ap.add_argument("--no-figures", action="store_true",
                    help="write the JSON and the .md only (no matplotlib output)")
    args = ap.parse_args(argv)

    t0 = time.time()
    res = compute()
    res["tables"] = report_tables(res)
    checks = self_check(res)
    bad = [c for c in checks if not c["ok"]]
    print(f"{NAME}: {len(checks) - len(bad)}/{len(checks)} self-checks pass")
    for c in checks:
        print(f"  [{'ok' if c['ok'] else 'FAIL'}] {c['check']}: {c['detail']}")
    if bad:
        raise SystemExit(f"{len(bad)} self-check(s) failed - nothing written")

    if not args.no_figures:
        fig_geometry(res, args.figdir)

    for path in (args.json, args.md):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(args.json, "w", encoding="utf-8") as fh:
        json.dump(build_json(res, checks), fh, indent=2, ensure_ascii=False)
        fh.write("\n")
    with open(args.md, "w", encoding="utf-8") as fh:
        fh.write(report(res))
    print(f"  {os.path.relpath(args.json, ROOT)}")
    print(f"  {os.path.relpath(args.md, ROOT)}")
    print(f"  {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()


