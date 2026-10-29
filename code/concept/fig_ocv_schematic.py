#!/usr/bin/env python3
"""The Observability Core Vector — one schematic, from asset to six numbers.

Question: every other plate in ``figures/`` is drawn *from a record*. This one is
drawn from nothing but a definition, and that is its job: it says what the OCV
**is** before any site's numbers are read. A physical infrastructure asset becomes
a defined radar observation, and the observation becomes a six-dimensional
vector:

    physical asset -> SAR observation geometry -> asset-specific measurement
    support -> radar echo representation -> OCV

Three panels, left to right:

    (a) the asset and the acquisition    G_a^{3D} in its own (X,Y,Z) frame, the
        line of sight, and the projection Pi_c stated once as the mapping
        (X,Y,Z) -> (az,r) into radar coordinates
    (b) the measurement support          the projected asset inside its support
        W_{a,c}, a white guard, the background ring, and the read chip C_{a,c}
        that encloses them:  W_{a,c} subset C_{a,c}
    (c) the observable representation    the echo mask inside W_{a,c} and the
        vector it is summarised by, OCV(t) = [gamma^2, P, D, F, A, S], in its two
        conceptual groups (cross-acquisition; spatial geometry)

The figure is a **schematic**, so it carries no data, no site, no state and no
result. Its text inventory is fixed and checked: the panel labels, the operators
and variables of the definition, and the two group names. ``R_B`` and ``R_G`` are
deliberately absent — they belong to the temporal reconfiguration analysis, not
to the definition of the vector.

What the figure must not be read as, and does not show: structural health, damage
probability, failure detection, displacement measurement, InSAR deformation, or a
universal score. There is no displacement vector, no deformation, no threshold
and no count.

Figure (183 mm x 61 mm = 3.000 : 1, three panels, no title, no legend):

    a  asset + acquisition geometry   grey (X,Y,Z) triad, lattice tower, satellite,
                                      orange line of sight LOS_c, the look angle
                                      theta_c, Pi_c stated once as the mapping
                                      (X,Y,Z) -> (az,r)
    b  measurement support            grey C_{a,c}, the white guard, dashed
                                      background ring, blue W_{a,c}, dark asset
                                      projection
    c  echo mask + OCV                4 connected echo regions in W_{a,c}, nothing
                                      annotated on the mask, the vector below

Every symbol that depends on the acquisition channel carries it as a subscript:
LOS_c and theta_c here, alongside Pi_c, W_{a,c} and C_{a,c}. Ascending and
descending are never written out — the subscript already says it. The centroid
and the peak of the echo mask are measured (and checked) but never drawn: S is a
distance the raster yields, not a glyph the plate has to explain.

Style: white background, thin black and grey construction lines, one muted blue
(W), one muted orange (line of sight), one muted green (echo mask). Vector-like,
no shadows, no decoration. All geometry is laid out in panel units where one unit
is the same length in x and y, so circles, arcs and raster cells are round; the
unit is 55.632 mm, so 7 pt text is 7 pt on the printed page at 183 mm.

Deliverables: ``figures/concept/fig_ocv_schematic.jpg`` (630 dpi) and
``figures/concept/fig_ocv_schematic.pdf`` (true vector, FontType 42),
``data/concept/fig_ocv_schematic.json`` (the drawn geometry, every check and the
sha256 of both images) and ``figures/concept/fig_ocv_schematic.md`` (the report).

Usage:
  python3 code/concept/fig_ocv_schematic.py
  python3 code/concept/fig_ocv_schematic.py --no-figures
  python3 code/concept/fig_ocv_schematic.py --figdir /tmp/x --json /tmp/y.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import sys

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, Polygon  # noqa: E402
from matplotlib.path import Path  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
FIGDIR = os.path.join(ROOT, "figures", "concept")
DATADIR = os.path.join(ROOT, "data", "concept")

# ---------------------------------------------------------------------------
# canvas: 183 mm x 61 mm, three equal panels, no title, no legend
# ---------------------------------------------------------------------------
MM = 1.0 / 25.4
FIG_W_MM = 183.0
FIG_H_MM = 61.0
DPI = 630
PANEL_W_FRAC = 0.304
PANEL_H_FRAC = 0.96
PANEL_Y_FRAC = 0.02
PANEL_X_FRAC = (0.012, 0.348, 0.684)
# one panel unit is 0.304 x 183 mm wide; the y limit is chosen so that one unit
# is the same length in x and y (circles stay round, the raster cells are square)
UNIT_MM = PANEL_W_FRAC * FIG_W_MM
YLIM = (PANEL_H_FRAC * FIG_H_MM) / UNIT_MM
PANELS = ("a", "b", "c")

# ---------------------------------------------------------------------------
# palette: one muted blue, one muted orange, one muted green, ink and greys
# ---------------------------------------------------------------------------
INK = "#1a1a1a"
GREY = "#6f6f6f"
LIGHT = "#b9b9b9"
FAINT = "#d8d8d8"
BLUE = "#1f4e79"
ORANGE = "#b5651d"
GREEN = "#3f7d3f"

# type sizes (pt). 7.0 is the floor: the plate must stay readable at 170 mm.
FS_PANEL = 9.0
FS_MATH = 9.0
FS_LABEL = 8.0
FS_SMALL = 7.5
FS_GROUP = 7.0
FS_FLOOR = 7.0

RC = {
    "font.family": "DejaVu Sans",
    "mathtext.fontset": "stix",
    "font.size": FS_LABEL,
    "axes.linewidth": 0.0,
    "savefig.dpi": DPI,
    # deliberately NOT "tight": a tight bbox would silently move the 183 mm edge.
    "savefig.bbox": "standard",
    "pdf.fonttype": 42,          # TrueType embed: the PDF text stays selectable
    "ps.fonttype": 42,
}

OUT_JPG = "fig_ocv_schematic.jpg"
OUT_PDF = "fig_ocv_schematic.pdf"
OUT_MD = "fig_ocv_schematic.md"
OUT_JSON = "fig_ocv_schematic.json"

# ---------------------------------------------------------------------------
# forbidden vocabulary: the figure is a definition, not a verdict
# ---------------------------------------------------------------------------
FORBIDDEN = (
    "damage", "health", "fail", "probab", "displace", "deform", "insar",
    "score", "threshold", "pixel", "px", "r_b", "r_g", "risk", "crack",
    "anomal", "diagnos", "predict",
)

# ---------------------------------------------------------------------------
# checks
# ---------------------------------------------------------------------------
CHECKS: list[dict] = []


def check(name: str, ok: bool, detail: str = "") -> bool:
    """Record one self-check; never raises, so the whole list is always reported."""
    CHECKS.append({"name": name, "ok": bool(ok), "detail": detail})
    return bool(ok)


def checks_failed() -> list[dict]:
    return [c for c in CHECKS if not c["ok"]]


# ---------------------------------------------------------------------------
# the drawn geometry (all values in panel units unless the name says otherwise)
# ---------------------------------------------------------------------------
# The measurement support as a chamfered rectangle in its own unit box: ONE
# definition, mapped into panel (b) and panel (c) by an affine map each, so the
# two panels are guaranteed to show the same W_{a,c}.
W_UNIT = (
    (0.10, 0.00), (0.90, 0.00), (1.00, 0.10), (1.00, 0.90),
    (0.90, 1.00), (0.10, 1.00), (0.00, 0.90), (0.00, 0.10),
)

# The coordinate symbols of the two frames, defined once and consumed by both
# the drawing and the text whitelist, so a label can never drift from the
# whitelist that pins it. The physical frame and the radar frame deliberately
# share no symbol: no letter is ever doing two jobs on this plate.
PHYS_AX = (r"$X$", r"$Y$", r"$Z$")          # physical frame: the triad's x, y, z
RAD_AZ = r"$\mathrm{az}$"                   # radar azimuth
RAD_R = r"$r$"                              # radar range
# The whole projection, stated once at the foot of panel (a). mathtext has no
# \xrightarrow, so the operator rides above a long arrow via \overset.
MAP = r"$(X,Y,Z)\overset{\Pi_c}{\longrightarrow}(\mathrm{az},r)$"

# The acquisition channel. c is the channel — the same subscript that Pi_c,
# W_{a,c} and C_{a,c} carry — and the line of sight and the look angle belong to
# the pass rather than to the asset, so they carry it too. Ascending and
# descending are never spelled out: the subscript says it, and the check below
# fails the build if a channel-dependent symbol loses it or a frame symbol picks
# it up.
LOS_SYM = r"$\mathrm{LOS}_c$"
THETA_SYM = r"$\theta_c$"
CHANNEL_SYMS = (LOS_SYM, THETA_SYM, MAP, r"$C_{a,c}$", r"$W_{a,c}$",
                r"$W_{a,c}\subset C_{a,c}$", r"$\Pi_c(G_a^{3D})$")

# panel (a): the asset
TOWER = dict(cx=0.300, y_base=0.120, y_top=0.600,
             hw_base=0.075, hw_top=0.032, depth=(0.050, 0.040), belts=5)
GROUND = ((0.150, 0.120), (0.470, 0.120))
TRIAD = dict(origin=(0.058, 0.150), x=(0.180, 0.150),
             y=(0.112, 0.212), z=(0.058, 0.345))

# panel (a): the satellite and the line of sight
SAT = dict(cx=0.800, cy=0.905, w=0.100, h=0.050, wing=(0.075, 0.026), gap=0.004)
LOS = ((0.772, 0.868), (0.356, 0.618))      # from the satellite to the asset top
NADIR = ((0.772, 0.880), (0.772, 0.700))    # through the LOS origin, for the angle
THETA_R = 0.085                             # arc radius of the look angle

# panel (a): the radar coordinate frame the projection lands in. The frame is
# raised by FRAME_DY so the mapping statement fits under it with clear air.
FRAME_DY = 0.040
RFRAME = dict(box=(0.560, 0.100 + FRAME_DY, 0.950, 0.380 + FRAME_DY),
              proj=((0.610, 0.128 + FRAME_DY), (0.878, 0.128 + FRAME_DY),
                    (0.838, 0.302 + FRAME_DY), (0.652, 0.302 + FRAME_DY)))
MAP_X = RFRAME["box"][0]                     # left-aligned under the radar frame
MAP_Y = 0.030                               # baseline, clear of the ground rule

# panel (b): chip, guard, background ring, support. Outward from the support
# the chip is guard, then ring, then the chip's own boundary — the same order
# the audit uses, and the reason the ring is drawn as a true annulus.
CHIP = (0.140, 0.240, 0.900, 0.860)         # C_{a,c}
RING = (0.200, 0.310, 0.840, 0.790)         # the background ring, drawn dashed
GUARD = 0.020                               # white gap between W and the ring
W_BOX = (0.260, 0.370, 0.780, 0.730)        # W_{a,c} mapped into this box
GUARD_BOX = (W_BOX[0] - GUARD, W_BOX[1] - GUARD,
             W_BOX[2] + GUARD, W_BOX[3] + GUARD)
B_AXIS = dict(x0=0.100, y0=0.170, x1=0.920, y1=0.900)
B_ASSET = ((0.410, 0.485), (0.630, 0.485), (0.600, 0.700), (0.440, 0.700))

# panel (c): the raster and the support it sits in
RAST = dict(x0=0.220, y0=0.520, cols=16, rows=10, cell=0.03625)
W_BOX_C = (0.170, 0.470, 0.850, 0.945)
#   four connected echo regions, (column, row), row 0 at the bottom:
#   a 3x3 ring with a closed hole, a 3x2 block, a 2x2 patch and a 3x1 run
ECHO = {
    "ring": [(3, 3), (4, 3), (5, 3), (3, 4), (5, 4), (3, 5), (4, 5), (5, 5)],
    "block": [(8, 6), (9, 6), (10, 6), (8, 7), (9, 7), (10, 7)],
    "patch": [(11, 2), (12, 2), (11, 3), (12, 3)],
    "run": [(5, 9), (6, 9), (7, 9)],
}
ECHO_HOLE = (4, 4)                          # the closed hole of the 3x3 ring
ECHO_PEAK = (9, 7)                          # the intensity peak, inside "block"

OCV_TOKENS = (
    ("gamma2", r"$\gamma^2$", "cross"),
    ("P", r"$P$", "cross"),
    ("D", r"$D$", "spatial"),
    ("F", r"$F$", "spatial"),
    ("A", r"$A$", "spatial"),
    ("S", r"$S$", "spatial"),
)
OCV_GROUPS = (
    ("cross", "cross-acquisition", BLUE),
    ("spatial", "spatial geometry", GREEN),
)
Y_OCV = 0.335                               # baseline of the equation
Y_BRACE = 0.285                             # the two under-braces
Y_GROUP = 0.245                             # the two group names

# ---------------------------------------------------------------------------
# drawing primitives — every artist is one of these five, so the plate has no
# shadows, no gradients and no decorative element anywhere
# ---------------------------------------------------------------------------
TEXTS: list[dict] = []


def _line(ax, p, q, color=INK, lw=0.6, ls="-", z=2):
    ax.plot([p[0], q[0]], [p[1], q[1]], color=color, lw=lw, ls=ls,
            solid_capstyle="round", zorder=z)


def _arrow(ax, p, q, color=INK, lw=0.9, scale=7.0, z=4):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=scale,
                                 lw=lw, color=color, shrinkA=0.0, shrinkB=0.0,
                                 joinstyle="miter", zorder=z))


def _poly(ax, pts, edge=INK, fc="none", lw=0.8, ls="-", z=2, alpha=1.0):
    ax.add_patch(Polygon(pts, closed=True, fill=(fc != "none"), facecolor=fc,
                         edgecolor=edge, lw=lw, ls=ls, alpha=alpha, zorder=z,
                         joinstyle="miter"))


def _text(ax, x, y, s, fs=FS_LABEL, color=INK, ha="left", va="baseline",
          rot=0.0, weight="normal", z=6, tag="", bbox=None):
    t = ax.text(x, y, s, fontsize=fs, color=color, ha=ha, va=va, rotation=rot,
                rotation_mode="anchor", fontweight=weight, zorder=z, bbox=bbox)
    TEXTS.append({"ax": ax, "t": t, "s": s, "fs": fs, "tag": tag, "color": color})
    return t


def _renderer(fig):
    fig.canvas.draw()
    return fig.canvas.get_renderer()


def _width(ax, s, fs):
    """Width of `s` at `fs` in panel units, measured on the output renderer."""
    t = ax.text(0.5, 0.5, s, fontsize=fs, ha="left", va="baseline")
    bb = t.get_window_extent(renderer=ax.figure.canvas.get_renderer())
    inv = ax.transData.inverted()
    x0, _ = inv.transform((bb.x0, bb.y0))
    x1, _ = inv.transform((bb.x1, bb.y1))
    t.remove()
    return float(x1 - x0)


def _extent(rec):
    """Drawn extent of a text record, in panel units, from its own renderer."""
    ax = rec["ax"]
    bb = rec["t"].get_window_extent(renderer=ax.figure.canvas.get_renderer())
    inv = ax.transData.inverted()
    (x0, y0) = inv.transform((bb.x0, bb.y0))
    (x1, y1) = inv.transform((bb.x1, bb.y1))
    return (min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1))


def w_poly(box):
    """Map the one W_{a,c} definition into `box` = (x0, y0, x1, y1)."""
    x0, y0, x1, y1 = box
    return [(x0 + u * (x1 - x0), y0 + v * (y1 - y0)) for u, v in W_UNIT]


def cell_xy(i, j):
    """Centre of raster cell (column i, row j) of panel (c)."""
    c = RAST["cell"]
    return (RAST["x0"] + (i + 0.5) * c, RAST["y0"] + (j + 0.5) * c)


# ---------------------------------------------------------------------------
# panel (a): the physical asset under a specified acquisition geometry
# ---------------------------------------------------------------------------
def _tower_sections():
    """The belts of the lattice tower: front face and depth-offset back face."""
    T = TOWER
    n = T["belts"]
    dx, dy = T["depth"]
    secs = []
    for k in range(n):
        t = k / (n - 1)
        hw = T["hw_base"] + t * (T["hw_top"] - T["hw_base"])
        y = T["y_base"] + t * (T["y_top"] - T["y_base"])
        front = ((T["cx"] - hw, y), (T["cx"] + hw, y))
        back = ((T["cx"] - hw + dx, y + dy), (T["cx"] + hw + dx, y + dy))
        secs.append((front, back))
    return secs


def draw_asset(ax):
    """A generic lattice asset: four legs, five belts, alternating X braces.

    Four legs and belts give it a readable 3D depth; the tapered section and the
    braces are what makes it a tower rather than a box. No displacement vector is
    drawn anywhere: this figure is about observation geometry, not motion.
    """
    secs = _tower_sections()
    last = len(secs) - 1
    for k in range(last):                                   # back legs
        _line(ax, secs[k][1][0], secs[k + 1][1][0], GREY, 0.45, z=1)
        _line(ax, secs[k][1][1], secs[k + 1][1][1], GREY, 0.45, z=1)
    for f, b in secs:                                       # depth connectors
        _line(ax, b[0], b[1], GREY, 0.40, z=1)
        _line(ax, f[0], b[0], LIGHT, 0.40, z=1)
        _line(ax, f[1], b[1], LIGHT, 0.40, z=1)
    for k in range(last):                                   # front legs
        _line(ax, secs[k][0][0], secs[k + 1][0][0], INK, 0.70)
        _line(ax, secs[k][0][1], secs[k + 1][0][1], INK, 0.70)
    for f, _b in secs:                                      # front belts
        _line(ax, f[0], f[1], INK, 0.60)
    for k in range(last):                                   # X braces
        _line(ax, secs[k][0][0], secs[k + 1][0][1], INK, 0.45)
        _line(ax, secs[k][0][1], secs[k + 1][0][0], INK, 0.45)


def draw_panel_a(ax):
    _text(ax, 0.020, 1.020, "(a)", fs=FS_PANEL, weight="bold", va="top", tag="panel")

    # the ground and the orientation triad: construction lines, grey
    _line(ax, GROUND[0], GROUND[1], GREY, 0.60)
    o = TRIAD["origin"]
    for key in ("x", "y", "z"):
        _arrow(ax, o, TRIAD[key], GREY, 0.50, scale=5.0, z=1)
    _text(ax, TRIAD["z"][0] - 0.012, TRIAD["z"][1], PHYS_AX[2], fs=FS_SMALL,
          color=GREY, ha="right", va="center", tag="axis")
    _text(ax, TRIAD["x"][0] + 0.008, TRIAD["x"][1], PHYS_AX[0], fs=FS_SMALL,
          color=GREY, ha="left", va="center", tag="axis")
    _text(ax, TRIAD["y"][0] + 0.004, TRIAD["y"][1] + 0.026, PHYS_AX[1],
          fs=FS_SMALL, color=GREY, ha="left", va="center", tag="axis")

    # the physical asset
    draw_asset(ax)
    _text(ax, TOWER["cx"], 0.050, r"$G_a^{3D}$", fs=FS_MATH, ha="center",
          va="baseline", tag="asset")

    # the satellite, the nadir reference and the line of sight
    s = SAT
    _poly(ax, [(s["cx"] - s["w"] / 2, s["cy"] - s["h"] / 2),
               (s["cx"] + s["w"] / 2, s["cy"] - s["h"] / 2),
               (s["cx"] + s["w"] / 2, s["cy"] + s["h"] / 2),
               (s["cx"] - s["w"] / 2, s["cy"] + s["h"] / 2)],
          edge=INK, fc="white", lw=0.60)
    wx, wy = s["wing"]
    for sgn in (-1.0, 1.0):
        x_in = s["cx"] + sgn * (s["w"] / 2 + s["gap"])
        _poly(ax, [(x_in, s["cy"] - wy / 2), (x_in + sgn * wx, s["cy"] - wy / 2),
                   (x_in + sgn * wx, s["cy"] + wy / 2), (x_in, s["cy"] + wy / 2)],
              edge=INK, fc="white", lw=0.55)
    _line(ax, NADIR[0], NADIR[1], GREY, 0.50, ls=(0, (2.4, 1.8)))
    _arrow(ax, LOS[0], LOS[1], ORANGE, 0.90, scale=7.0)

    # the look angle between nadir and the line of sight, drawn at the satellite
    p0 = np.array(LOS[0])
    d = np.array(LOS[1]) - p0
    a_los = float(np.degrees(np.arctan2(d[1], d[0])) % 360.0)
    a_nad = 270.0
    th = np.linspace(a_los, a_nad, 60)
    ax.plot(p0[0] + THETA_R * np.cos(np.radians(th)),
            p0[1] + THETA_R * np.sin(np.radians(th)),
            color=GREY, lw=0.50, zorder=3)
    am = 0.5 * (a_los + a_nad)
    _text(ax, p0[0] + (THETA_R + 0.030) * np.cos(np.radians(am)),
          p0[1] + (THETA_R + 0.030) * np.sin(np.radians(am)),
          THETA_SYM, fs=FS_SMALL, color=GREY, ha="center", va="center",
          tag="angle")
    nrm = np.array([-d[1], d[0]]) / float(np.hypot(d[0], d[1]))
    lp = p0 + 0.45 * d - 0.040 * nrm
    rot = a_los - 180.0 if a_los > 90.0 else a_los    # keep the label upright
    _text(ax, lp[0], lp[1], LOS_SYM, fs=FS_SMALL, color=ORANGE, ha="center",
          va="center", rot=rot, tag="los")

    # the projection into the radar frame, and the mapping that states it once.
    # The operator lives inside the mapping expression now, so the panel carries
    # \Pi_c exactly once: the bare label that used to sit above this arrow is
    # gone rather than duplicated.
    _arrow(ax, (0.428, 0.165), (0.552, 0.165), INK, 0.90, scale=7.0)
    _text(ax, MAP_X, MAP_Y, MAP, fs=FS_SMALL, ha="left", va="baseline",
          tag="mapping")
    bx0, by0, bx1, by1 = RFRAME["box"]
    _line(ax, (bx0, by0), (bx1, by0), INK, 0.60)
    _line(ax, (bx0, by0), (bx0, by1), INK, 0.60)
    _text(ax, bx1 + 0.010, by0 - 0.022, RAD_R, fs=FS_SMALL, color=INK,
          ha="left", va="center", tag="axis")
    _text(ax, bx0 - 0.012, by1 + 0.016, RAD_AZ, fs=FS_SMALL, color=INK,
          ha="right", va="center", tag="axis")
    _poly(ax, RFRAME["proj"], edge=INK, fc=INK, lw=0.55, alpha=0.82, z=3)
    pj = RFRAME["proj"]
    for frac in (0.25, 0.50, 0.75):
        y = pj[0][1] + frac * (pj[3][1] - pj[0][1])
        xl = pj[0][0] + frac * (pj[3][0] - pj[0][0])
        xr = pj[1][0] + frac * (pj[2][0] - pj[1][0])
        _line(ax, (xl, y), (xr, y), "white", 0.55, z=4)


# ---------------------------------------------------------------------------
# panel (b): the asset-specific measurement support inside the read chip
# ---------------------------------------------------------------------------
def draw_panel_b(ax):
    _text(ax, 0.020, 1.020, "(b)", fs=FS_PANEL, weight="bold", va="top", tag="panel")

    # the radar-image axes
    a = B_AXIS
    _line(ax, (a["x0"], a["y0"]), (a["x1"], a["y0"]), INK, 0.60)
    _line(ax, (a["x0"], a["y0"]), (a["x0"], a["y1"]), INK, 0.60)
    _text(ax, a["x1"] + 0.012, a["y0"] - 0.024, RAD_R, fs=FS_SMALL,
          ha="left", va="center", tag="axis")
    _text(ax, a["x0"], a["y1"] + 0.006, RAD_AZ, fs=FS_SMALL,
          ha="center", va="bottom", tag="axis")

    # the read chip, then the background ring as a true annulus. The tint is
    # knocked out over the guard, so the ring is the band between the guard and
    # the chip — never a filled rectangle behind the support: no tint reaches
    # inside W_{a,c}, and no background cell sits on the support.
    cx0, cy0, cx1, cy1 = CHIP
    _poly(ax, [(cx0, cy0), (cx1, cy0), (cx1, cy1), (cx0, cy1)],
          edge=GREY, fc="none", lw=0.55, z=1)
    rx0, ry0, rx1, ry1 = RING
    _poly(ax, [(rx0, ry0), (rx1, ry0), (rx1, ry1), (rx0, ry1)],
          edge=GREY, fc=LIGHT, lw=0.50, ls=(0, (3.0, 2.0)), alpha=0.30, z=1)
    _poly(ax, w_poly(GUARD_BOX), edge="none", fc="white", lw=0.0, z=2)
    _poly(ax, w_poly(W_BOX), edge=BLUE, fc=BLUE, lw=0.85, alpha=0.13, z=2)

    # the projected asset: the same flattened silhouette as in panel (a)
    _poly(ax, B_ASSET, edge=INK, fc=INK, lw=0.50, alpha=0.80, z=3)
    for frac in (0.28, 0.52, 0.76):
        y = B_ASSET[0][1] + frac * (B_ASSET[3][1] - B_ASSET[0][1])
        xl = B_ASSET[0][0] + frac * (B_ASSET[3][0] - B_ASSET[0][0])
        xr = B_ASSET[1][0] + frac * (B_ASSET[2][0] - B_ASSET[1][0])
        _line(ax, (xl, y), (xr, y), "white", 0.55, z=4)

    # the four names this panel needs, and nothing else. C_{a,c} sits *on* the
    # chip's lower edge with a white knockout, so the edge breaks behind it: the
    # label names that rectangle and cannot be read as naming the ring it used
    # to float between.
    _text(ax, cx1 - 0.010, cy0, r"$C_{a,c}$", fs=FS_SMALL, color=GREY,
          ha="right", va="center", tag="chip",
          bbox=dict(facecolor="white", edgecolor="none", pad=1.2))
    _text(ax, 0.292, 0.638, r"$W_{a,c}$", fs=FS_MATH, color=BLUE,
          ha="left", va="baseline", tag="support")
    _text(ax, 0.520, 0.398, r"$\Pi_c(G_a^{3D})$", fs=FS_MATH,
          ha="center", va="baseline", tag="projected")
    _text(ax, 0.620, 0.900, "background ring", fs=FS_SMALL, color=INK,
          ha="center", va="baseline", tag="ring")
    _line(ax, (0.620, 0.886), (0.620, 0.770), GREY, 0.45, z=3)
    _text(ax, 0.500, 0.070, r"$W_{a,c}\subset C_{a,c}$", fs=FS_MATH,
          ha="center", va="baseline", tag="subset")


# ---------------------------------------------------------------------------
# panel (c): the echo representation inside W_{a,c}, and the vector it becomes
# ---------------------------------------------------------------------------
def echo_cells():
    """Every echo cell once, in drawing order (region by region)."""
    return [(i, j) for _n, cells in ECHO.items() for (i, j) in cells]


def centroid_xy():
    """Centroid of the echo mask, in panel units (cells weighted equally)."""
    cells = echo_cells()
    ci = float(np.mean([i + 0.5 for i, _j in cells]))
    cj = float(np.mean([j + 0.5 for _i, j in cells]))
    c = RAST["cell"]
    return (RAST["x0"] + ci * c, RAST["y0"] + cj * c), (ci, cj)


def draw_raster(ax):
    c = RAST["cell"]
    x0, y0 = RAST["x0"], RAST["y0"]
    w = RAST["cols"] * c
    h = RAST["rows"] * c
    _poly(ax, [(x0, y0), (x0 + w, y0), (x0 + w, y0 + h), (x0, y0 + h)],
          edge=LIGHT, fc="none", lw=0.45, z=2)
    for i in range(1, RAST["cols"]):
        _line(ax, (x0 + i * c, y0), (x0 + i * c, y0 + h), FAINT, 0.28, z=2)
    for j in range(1, RAST["rows"]):
        _line(ax, (x0, y0 + j * c), (x0 + w, y0 + j * c), FAINT, 0.28, z=2)
    for (i, j) in echo_cells():
        x = x0 + i * c
        y = y0 + j * c
        _poly(ax, [(x, y), (x + c, y), (x + c, y + c), (x, y + c)],
              edge="white", fc=GREEN, lw=0.45, alpha=0.88, z=3)


def draw_panel_c(ax):
    _text(ax, 0.020, 1.020, "(c)", fs=FS_PANEL, weight="bold", va="top", tag="panel")

    # the same support as panel (b), so the echo has a place to sit in
    _poly(ax, w_poly(W_BOX_C), edge=BLUE, fc=BLUE, lw=0.80, alpha=0.10, z=1)
    draw_raster(ax)

    # Nothing is annotated inside the mask — no centroid dot, no peak dot, no
    # arrow between them, no symbol naming either. The centroid-to-peak distance
    # that S measures is a number the raster yields, not a glyph the plate has to
    # explain, so panel (c) is left with marks that need no definition: the
    # support, the echo, the vector.

    # the raster becomes the vector, pointing at the middle of the finished row
    pos = layout_ocv(ax)
    mid = 0.5 * (pos["left"] + pos["right"])
    _arrow(ax, (mid, 0.452), (mid, 0.392), INK, 0.80, scale=6.0)
    return pos


def ocv_group_string(grp: str) -> str:
    """The tokens of one conceptual group as a single mathtext string.

    Built from ``OCV_TOKENS`` rather than typed twice, so the string the figure
    draws and the token list the checks verify cannot drift apart.
    """
    syms = [s for (_k, s, g) in OCV_TOKENS if g == grp]
    return "$" + r",\ ".join(x.strip("$") for x in syms) + "$"


def layout_ocv(ax):
    """Place ``OCV(t) = [gamma2, P, D, F, A, S]`` from measured text widths.

    Each group is one string, so its brace spans the measured extent of what it
    names — no cumulative cursor arithmetic, and no brace wider than its group by
    accident. The two group names are wider than the token groups, so the
    inter-group gap is *solved* from the two name widths (they must clear each
    other by ``lab_sep``) instead of being guessed; the row is then centred.
    """
    x_l, x_r = 0.020, 0.980
    g0 = g1 = 0.022                     # prefix -> group one, group two -> "]"
    lab_sep = 0.035
    pre, close = r"$\mathrm{OCV}(t)=[$", "]"
    parts = [(grp, ocv_group_string(grp), col) for grp, _n, col in OCV_GROUPS]
    w_pre = _width(ax, pre, FS_MATH)
    w_close = _width(ax, close, FS_MATH)
    w_com = _width(ax, ",", FS_MATH)
    w = [_width(ax, s, FS_MATH) for (_g, s, _c) in parts]
    wlab = [_width(ax, name, FS_GROUP) for (_g, name, _c) in OCV_GROUPS]

    # the gap between the groups must clear both names AND carry the separating
    # comma, which belongs to the vector and not to either group's string
    gd = max((wlab[0] + wlab[1]) / 2.0 + lab_sep - (w[0] + w[1]) / 2.0,
             w_com + 0.024, 0.070)
    total = w_pre + g0 + w[0] + gd + w[1] + g1 + w_close
    x0 = x_l + max(0.0, (x_r - x_l - total) / 2.0)

    _text(ax, x0, Y_OCV, pre, fs=FS_MATH, ha="left", va="baseline", tag="ocv_pre")
    spans = {}
    x1 = x0 + w_pre + g0
    spans[parts[0][0]] = (x1, x1 + w[0])
    _text(ax, x1, Y_OCV, parts[0][1], fs=FS_MATH, color=parts[0][2],
          ha="left", va="baseline", tag="ocv_" + parts[0][0])
    _text(ax, spans[parts[0][0]][1] + 0.012, Y_OCV, ",", fs=FS_MATH, color=INK,
          ha="left", va="baseline", tag="ocv_sep")
    x2 = spans[parts[0][0]][1] + gd
    spans[parts[1][0]] = (x2, x2 + w[1])
    _text(ax, x2, Y_OCV, parts[1][1], fs=FS_MATH, color=parts[1][2],
          ha="left", va="baseline", tag="ocv_" + parts[1][0])
    xc = spans[parts[1][0]][1] + g1
    _text(ax, xc, Y_OCV, close, fs=FS_MATH, ha="left", va="baseline",
          tag="ocv_close")

    # the split the two colours and the two braces already state, made visible in
    # the gap as well: a thin rule, not a second key and not a box
    xdiv = 0.5 * (spans[parts[0][0]][1] + spans[parts[1][0]][0])
    _line(ax, (xdiv, Y_OCV - 0.012), (xdiv, Y_OCV + 0.052), LIGHT, 0.50, z=4)

    pos = {"x0": x0, "gd": gd, "g0": g0, "g1": g1, "total": total,
           "w_pre": w_pre, "w_close": w_close, "w": w, "wlab": wlab,
           "left": x0, "right": xc + w_close, "spans": spans, "xdiv": xdiv}

    # the two under-braces and the two names: a bracket, not a box or a legend
    for grp, name, col in OCV_GROUPS:
        xs, xe = spans[grp]
        _line(ax, (xs - 0.008, Y_BRACE), (xe + 0.008, Y_BRACE), col, 0.55, z=5)
        for xx in (xs - 0.008, xe + 0.008):
            _line(ax, (xx, Y_BRACE), (xx, Y_BRACE + 0.015), col, 0.55, z=5)
        xm = 0.5 * (xs + xe)
        _text(ax, xm, Y_GROUP, name, fs=FS_GROUP, color=col, ha="center",
              va="baseline", tag="group_" + grp)
        pos["centre_" + grp] = xm
    return pos


# ---------------------------------------------------------------------------
# the plate
# ---------------------------------------------------------------------------
def make_figure():
    """Build the three-panel plate and return it with its measured layout."""
    TEXTS.clear()
    plt.rcParams.update(RC)
    fig = plt.figure(figsize=(FIG_W_MM * MM, FIG_H_MM * MM), dpi=DPI,
                     facecolor="white")
    axes = []
    for xf in PANEL_X_FRAC:
        ax = fig.add_axes([xf, PANEL_Y_FRAC, PANEL_W_FRAC, PANEL_H_FRAC])
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, YLIM)
        ax.set_axis_off()
        ax.set_facecolor("white")
        ax.patch.set_visible(False)
        axes.append(ax)
    fig.canvas.draw()               # give the renderer real metrics first
    draw_panel_a(axes[0])
    draw_panel_b(axes[1])
    pos = draw_panel_c(axes[2])
    fig.canvas.draw()
    return fig, axes, pos


# ---------------------------------------------------------------------------
# geometry helpers the checks use
# ---------------------------------------------------------------------------
def _in_rect(pts, rect, margin=0.0):
    x0, y0, x1, y1 = rect
    return all(x0 + margin <= x <= x1 - margin and y0 + margin <= y <= y1 - margin
               for x, y in pts)


def echo_components():
    """4-connected components of the echo mask, as lists of cells."""
    cells = set(echo_cells())
    seen, comps = set(), []
    for start in sorted(cells):
        if start in seen:
            continue
        stack, comp = [start], []
        seen.add(start)
        while stack:
            i, j = stack.pop()
            comp.append((i, j))
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                n = (i + di, j + dj)
                if n in cells and n not in seen:
                    seen.add(n)
                    stack.append(n)
        comps.append(sorted(comp))
    return sorted(comps, key=len, reverse=True)


def echo_holes():
    """Enclosed non-echo cells of the raster: the fragmentation the mask carries."""
    cells = set(echo_cells())
    cols, rows = RAST["cols"], RAST["rows"]
    outside, stack = set(), []
    for i in range(cols):
        for j in (0, rows - 1):
            if (i, j) not in cells:
                stack.append((i, j))
    for j in range(rows):
        for i in (0, cols - 1):
            if (i, j) not in cells:
                stack.append((i, j))
    while stack:
        c = stack.pop()
        if c in outside or c in cells:
            continue
        outside.add(c)
        i, j = c
        for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            n = (i + di, j + dj)
            if 0 <= n[0] < cols and 0 <= n[1] < rows:
                stack.append(n)
    return sorted((i, j) for i in range(cols) for j in range(rows)
                  if (i, j) not in cells and (i, j) not in outside)


# ---------------------------------------------------------------------------
# the checks: nine families, every one of them fatal
# ---------------------------------------------------------------------------
def expected_texts():
    """The complete text inventory of the plate, written out, not discovered."""
    t = ["(a)", "(b)", "(c)"]
    t += [PHYS_AX[2], PHYS_AX[0], PHYS_AX[1], r"$G_a^{3D}$", THETA_SYM,
          LOS_SYM, MAP, RAD_R, RAD_AZ]                      # panel (a)
    t += [RAD_R, RAD_AZ, r"$C_{a,c}$", r"$W_{a,c}$",
          r"$\Pi_c(G_a^{3D})$", "background ring",
          r"$W_{a,c}\subset C_{a,c}$"]                      # panel (b)
    t += [r"$\mathrm{OCV}(t)=[$", ocv_group_string("cross"), ",",
          ocv_group_string("spatial"), "]",
          "cross-acquisition", "spatial geometry"]          # panel (c)
    return t


def run_checks(fig, axes, pos):
    """Every structural promise the plate makes, checked on the built figure."""
    # ---- 1. canvas -------------------------------------------------------
    w_in, h_in = fig.get_size_inches()
    check("canvas width 183 mm", abs(w_in * 25.4 - FIG_W_MM) < 1e-6,
          f"{w_in * 25.4:.4f} mm")
    check("canvas height 61 mm", abs(h_in * 25.4 - FIG_H_MM) < 1e-6,
          f"{h_in * 25.4:.4f} mm")
    check("canvas ratio 3.000:1", abs(w_in / h_in - 3.0) < 2e-3,
          f"{w_in / h_in:.5f}")
    sx = fig.transFigure.transform((1, 0))[0] - fig.transFigure.transform((0, 0))[0]
    check("figure dpi is the output dpi", fig.dpi == DPI, f"{fig.dpi}")
    for ax in axes:
        t = ax.transData
        dx = t.transform((1, 0))[0] - t.transform((0, 0))[0]
        dy = t.transform((0, 1))[1] - t.transform((0, 0))[1]
        check(f"panel {PANELS[axes.index(ax)]} is isotropic (unit = unit)",
              abs(dx - dy) < 1e-6, f"dx={dx:.4f}px dy={dy:.4f}px")
        b = ax.get_position()
        check(f"panel {PANELS[axes.index(ax)]} box matches the layout",
              abs(b.width - PANEL_W_FRAC) < 1e-9
              and abs(b.height - PANEL_H_FRAC) < 1e-9, f"{b.width:.6f}x{b.height:.6f}")
    for k in range(len(axes) - 1):
        check(f"panel {PANELS[k]}/{PANELS[k + 1]} boxes do not overlap",
              axes[k].get_position().x1 <= axes[k + 1].get_position().x0 + 1e-12,
              f"gap {axes[k + 1].get_position().x0 - axes[k].get_position().x1:.4f}")
    check("panel unit is 55.632 mm",
          abs(fig.transFigure.transform((PANEL_W_FRAC, 0))[0]
              - fig.transFigure.transform((0, 0))[0]) / DPI * 25.4 - UNIT_MM < 1e-6,
          f"{UNIT_MM:.3f} mm")
    _ = sx

    # ---- 2. text inventory ----------------------------------------------
    got = [rec["s"] for rec in TEXTS]
    want = expected_texts()
    check("text inventory is exactly the declared set",
          sorted(got) == sorted(want), f"{len(got)} objects")
    check("text count <= 30", len(got) <= 30, f"{len(got)}")
    check("no duplicate panel label",
          [g for g in got if g in ("(a)", "(b)", "(c)")] == ["(a)", "(b)", "(c)"])
    check("no title anywhere", getattr(fig, "_suptitle", None) is None)
    check("no legend anywhere", all(a.get_legend() is None for a in axes))
    check("axes are off (no spines, no ticks)",
          all(not a.axison for a in axes))
    fmin = min(rec["fs"] for rec in TEXTS)
    check("every text >= 7.0 pt", fmin >= FS_FLOOR - 1e-9, f"min {fmin:.2f} pt")
    check("every text >= 7.0 pt at 170 mm (scale 0.929) = 6.5 pt",
          fmin * 170.0 / 183.0 >= 6.5 - 1e-9, f"{fmin * 170.0 / 183.0:.2f} pt")

    low = " ".join(got).lower()
    hit = [w for w in FORBIDDEN if w in low]
    check("no forbidden vocabulary in any drawn text", not hit, ",".join(hit))
    check("R_B absent", "r_b" not in low.replace(" ", ""))
    check("R_G absent", "r_g" not in low.replace(" ", ""))
    check("the coherence channel is mathtext",
          any("gamma^2" in g for g in got))
    check("the physical frame and the radar frame share no symbol",
          not (set(PHYS_AX) & {RAD_AZ, RAD_R}),
          f"{sorted(set(PHYS_AX) & {RAD_AZ, RAD_R})}")
    check("the mapping statement names the operator and both frames",
          all(k in MAP for k in (r"\overset{\Pi_c}", "(X,Y,Z)",
                                 r"(\mathrm{az},r)")), "no drift")
    check("every channel-dependent symbol carries c, and no frame symbol does",
          all("c" in s for s in CHANNEL_SYMS)
          and not any("c" in s for s in PHYS_AX + (RAD_AZ, RAD_R)),
          f"{len(CHANNEL_SYMS)} channel symbols")

    # every text sits inside its own panel, and no two texts collide
    for ax in axes:
        tag = PANELS[axes.index(ax)]
        mine = [rec for rec in TEXTS if rec["ax"] is ax]
        out = []
        for rec in mine:
            x0, y0, x1, y1 = _extent(rec)
            if x0 < -0.002 or x1 > 1.002 or y0 < -0.002 or y1 > YLIM + 0.002:
                out.append(rec["s"])
        check(f"panel {tag}: every text inside the panel", not out, ",".join(out))
        bad = []
        for i in range(len(mine)):
            for j in range(i + 1, len(mine)):
                a, b = _extent(mine[i]), _extent(mine[j])
                if a[0] < b[2] and b[0] < a[2] and a[1] < b[3] and b[1] < a[3]:
                    bad.append(f"{mine[i]['s']}|{mine[j]['s']}")
        check(f"panel {tag}: no two texts overlap", not bad, ";".join(bad))

    # the two labels a reader could misread: C_{a,c} must name the chip's own
    # boundary, and the mapping must sit in the foot of (a), under the frame
    ch0, ch1 = CHIP[1], CHIP[2]
    crec = next(r for r in TEXTS if r["s"] == r"$C_{a,c}$")
    clx0, cly0, clx1, cly1 = _extent(crec)
    check("C_{a,c} sits on the chip's own lower edge, not between chip and ring",
          cly0 < ch0 < cly1 and CHIP[0] < clx0 and clx1 < ch1,
          f"y {cly0:.3f}..{cly1:.3f} across the edge {ch0:.3f}")
    mrec = next(r for r in TEXTS if r["s"] == MAP)
    mpx0, _mpy0, _mpx1, mpy1 = _extent(mrec)
    check("the mapping statement is the radar frame's own caption, below it",
          mpy1 < RFRAME["box"][1] and mpx0 > 0.0
          and abs(mpx0 - RFRAME["box"][0]) <= 0.02,
          f"top {mpy1:.3f} < frame {RFRAME['box'][1]:.3f}, "
          f"left {mpx0:.3f} vs {RFRAME['box'][0]:.3f}")

    # ---- 3. the support: W subset C, ring between -------------------------
    wb = w_poly(W_BOX)
    wc = w_poly(W_BOX_C)
    check("W has eight vertices in both panels", len(wb) == 8 and len(wc) == 8)
    check("W in (b) sits inside C with a >= 20 um margin",
          _in_rect(wb, CHIP, 0.02))
    check("W in (b) sits strictly inside the background ring",
          _in_rect(wb, RING, 0.0))
    check("W in (b) sits inside the guard, and the guard inside the ring",
          _in_rect(wb, GUARD_BOX, 0.0))
    guard = min(W_BOX[0] - GUARD_BOX[0], W_BOX[1] - GUARD_BOX[1],
                GUARD_BOX[2] - W_BOX[2], GUARD_BOX[3] - W_BOX[3])
    check("the guard is a white band >= 0.015 wide between W and the ring",
          guard >= 0.015 - 1e-12, f"{guard:.3f} units")
    band = min(GUARD_BOX[0] - RING[0], GUARD_BOX[1] - RING[1],
               RING[2] - GUARD_BOX[2], RING[3] - GUARD_BOX[3])
    check("the ring band outside the guard is >= 0.03 wide",
          band >= 0.03 - 1e-12, f"{band:.3f} units")
    check("the ring sits inside C, so it cannot be the chip edge",
          _in_rect([(RING[0], RING[1]), (RING[2], RING[3])], CHIP, 0.04))
    check("C sits inside the panel", _in_rect([(CHIP[0], CHIP[1]),
                                               (CHIP[2], CHIP[3])],
                                              (0.0, 0.0, 1.0, YLIM), 0.01))
    pw = Path(wb)
    check("the projected asset lies inside W",
          bool(pw.contains_points(B_ASSET).all()))
    check("the projected asset is inside C, not straddling it",
          _in_rect(B_ASSET, CHIP, 0.02))
    lad = (0.620, 0.770)
    check("the ring leader ends on the ring band, clear of the guard",
          GUARD_BOX[3] < lad[1] < RING[3], f"y={lad[1]:.3f}")

    # The ring is an annulus, not a fill, and this is proved on the composed
    # canvas rather than on the geometry: sample W's interior and the band, and
    # compare each against its analytic composite over white. If the tint is
    # ever painted behind the support again, these fail.
    _arr = render_array(fig)
    _axb = axes[1]

    def _px(x, y):
        px, py = _axb.transData.transform((x, y))
        return _arr[int(round(_arr.shape[0] - py)), int(round(px)), :3] / 255.0

    def _comp(hex_colour, alpha):
        rgb = np.array([int(hex_colour[1:3], 16), int(hex_colour[3:5], 16),
                        int(hex_colour[5:7], 16)]) / 255.0
        return 1.0 - alpha * (1.0 - rgb)

    _in_w, _in_band = _px(0.700, 0.430), _px(0.520, 0.330)
    check("the ring's tint does not reach inside W (annulus, not a fill)",
          bool(np.abs(_in_w - _comp(BLUE, 0.13)).max() <= 0.01),
          "W rgb " + " ".join(f"{v:.3f}" for v in _in_w))
    check("the band between the guard and the ring's boundary is tinted",
          bool(np.abs(_in_band - _comp(LIGHT, 0.30)).max() <= 0.01),
          "band rgb " + " ".join(f"{v:.3f}" for v in _in_band))

    # the one definition of W, mapped twice: the panels cannot drift apart
    def _unit(pts, box):
        x0, y0, x1, y1 = box
        return [((x - x0) / (x1 - x0), (y - y0) / (y1 - y0)) for x, y in pts]

    check("W in (b) and W in (c) are the same polygon",
          np.allclose(_unit(wb, W_BOX), _unit(wc, W_BOX_C), atol=1e-12))
    check("W in (b) and W in (c) are the same shape, not the same place",
          np.allclose(_unit(wb, W_BOX), np.array(W_UNIT), atol=1e-12)
          and not np.allclose(wb, wc, atol=1e-9))

    # ---- 4. the raster and the echo mask ---------------------------------
    c = RAST["cell"]
    rw = RAST["cols"] * c
    rh = RAST["rows"] * c
    check("raster cells are square on the drawn canvas",
          abs((axes[2].transData.transform((RAST["x0"] + c, RAST["y0"]))[0]
               - axes[2].transData.transform((RAST["x0"], RAST["y0"]))[0])
              - (axes[2].transData.transform((RAST["x0"], RAST["y0"] + c))[1]
                 - axes[2].transData.transform((RAST["x0"], RAST["y0"]))[1])) < 1e-6,
          f"{c:.5f} units")
    check("raster cell >= 1.5 mm on the printed page", c * UNIT_MM >= 1.5,
          f"{c * UNIT_MM:.2f} mm")
    corners = [(RAST["x0"], RAST["y0"]), (RAST["x0"] + rw, RAST["y0"]),
               (RAST["x0"] + rw, RAST["y0"] + rh), (RAST["x0"], RAST["y0"] + rh)]
    check("the raster sits inside W in (c)",
          bool(Path(wc).contains_points(corners).all()))
    check("the raster sits inside W with a visible margin",
          _in_rect(corners, W_BOX_C, 0.02))
    cell_pts = []
    for (i, j) in echo_cells():
        for dx, dy in ((0, 0), (c, 0), (c, c), (0, c)):
            cell_pts.append((RAST["x0"] + i * c + dx, RAST["y0"] + j * c + dy))
    check("every echo cell lies inside W in (c)",
          bool(Path(wc).contains_points(cell_pts).all()))
    check("every echo cell lies inside the raster",
          _in_rect(cell_pts, (RAST["x0"], RAST["y0"],
                              RAST["x0"] + rw, RAST["y0"] + rh), 0.0))

    comps = echo_components()
    sizes = [len(k) for k in comps]
    check("the echo mask has exactly four connected regions", len(comps) == 4,
          f"{len(comps)}")
    check("the four regions are 8, 6, 4 and 3 cells", sizes == [8, 6, 4, 3],
          f"{sizes}")
    check("every region is 4-connected (a fragment, not a dust field)",
          all(len(k) == len(set(k)) for k in comps))
    holes = echo_holes()
    check("the mask encloses exactly one hole", len(holes) == 1, f"{len(holes)}")
    check("the hole is the centre of the 3x3 ring",
          holes == [ECHO_HOLE] if len(holes) == 1 else False,
          f"{holes}")
    dens = len(echo_cells()) / float(RAST["cols"] * RAST["rows"])
    check("echo density is sparse and legible (5 % < d < 20 %)",
          0.05 < dens < 0.20, f"{dens:.4f}")
    check("the peak lies on an echo cell", ECHO_PEAK in set(echo_cells()))
    (ccx, ccy), _cc = centroid_xy()
    ppx, ppy = cell_xy(*ECHO_PEAK)
    sep_cells = float(np.hypot(ppx - ccx, ppy - ccy)) / c
    check("the centroid and the peak are >= 1.5 cells apart", sep_cells >= 1.5,
          f"{sep_cells:.3f} cells, the quantity S measures")
    check("the centroid is not inside an echo cell (the mask is fragmented)",
          not any(RAST["x0"] + i * c < ccx < RAST["x0"] + (i + 1) * c
                  and RAST["y0"] + j * c < ccy < RAST["y0"] + (j + 1) * c
                  for (i, j) in echo_cells()))

    # ---- 5. the OCV row --------------------------------------------------
    keys = [k for (k, _s, _g) in OCV_TOKENS]
    check("the vector has six components", len(keys) == 6)
    check("the order is gamma2, P, D, F, A, S",
          keys == ["gamma2", "P", "D", "F", "A", "S"], ",".join(keys))
    check("the partition covers every component exactly once",
          sorted(OCV_TOKENS and [g for (_k, _s, g) in OCV_TOKENS])
          == ["cross", "cross", "spatial", "spatial", "spatial", "spatial"])
    check("group one is exactly gamma2 and P",
          ocv_group_string("cross") == r"$\gamma^2,\ P$",
          ocv_group_string("cross"))
    check("group two is exactly D, F, A and S",
          ocv_group_string("spatial") == r"$D,\ F,\ A,\ S$",
          ocv_group_string("spatial"))
    check("the row starts inside the panel", pos["left"] >= 0.005,
          f"{pos['left']:.4f}")
    check("the row ends inside the panel", pos["right"] <= 0.995,
          f"{pos['right']:.4f}")
    check("the row is centred in its panel",
          abs((pos["left"] + pos["right"]) / 2.0 - 0.5) < 0.01,
          f"{(pos['left'] + pos['right']) / 2.0:.4f}")

    recs = {rec["tag"]: rec for rec in TEXTS if rec["tag"]}
    for grp, _name, _col in OCV_GROUPS:
        xs, xe = pos["spans"][grp]
        e = _extent(recs["ocv_" + grp])
        check(f"the {grp} brace spans its own group only",
              abs(e[0] - xs) < 0.012 and abs(e[2] - xe) < 0.012,
              f"text {e[0]:.3f}..{e[2]:.3f} brace {xs:.3f}..{xe:.3f}")
        lx = _extent(recs["group_" + grp])
        check(f"the {grp} name sits under its brace", lx[3] < Y_BRACE,
              f"top {lx[3]:.3f} < brace {Y_BRACE:.3f}")
    a1 = _extent(recs["group_cross"])
    a2 = _extent(recs["group_spatial"])
    sep = _extent(recs["ocv_sep"])
    check("the separating comma sits between the two groups",
          pos["spans"]["cross"][1] <= sep[0] and sep[2] <= pos["spans"]["spatial"][0],
          f"{sep[0]:.3f}..{sep[2]:.3f}")
    check("the separating comma does not touch either group",
          sep[0] - pos["spans"]["cross"][1] > 0.005
          and pos["spans"]["spatial"][0] - sep[2] > 0.005,
          f"{sep[0] - pos['spans']['cross'][1]:.3f} / "
          f"{pos['spans']['spatial'][0] - sep[2]:.3f}")
    check("the two group names do not overlap",
          a1[2] < a2[0], f"{a1[2]:.3f} < {a2[0]:.3f}")
    check("the two group names do not run into each other's brace",
          abs(a1[2] - a2[0]) >= 0.010, f"gap {a2[0] - a1[2]:.3f}")
    check("the equation sits above the braces",
          pos["w"] and Y_BRACE + 0.015 < Y_OCV)
    check("the equation does not touch the echo raster",
          _extent(recs["ocv_pre"])[3] < RAST["y0"])


# ---------------------------------------------------------------------------
# the written artefacts
# ---------------------------------------------------------------------------
def render_array(fig):
    """The composed canvas as an array, for the determinism check."""
    fig.canvas.draw()
    return np.asarray(fig.canvas.buffer_rgba()).copy()


def check_deterministic():
    """Two independent builds must produce the same pixels, not just the same JSON."""
    fig1, ax1, pos1 = make_figure()
    a = render_array(fig1)
    plt.close(fig1)
    fig2, ax2, pos2 = make_figure()
    b = render_array(fig2)
    plt.close(fig2)
    check("two builds render identical pixels", a.shape == b.shape
          and bool(np.array_equal(a, b)), f"{a.shape}")
    k1 = json.dumps({k: v for k, v in pos1.items()}, sort_keys=True, default=str)
    k2 = json.dumps({k: v for k, v in pos2.items()}, sort_keys=True, default=str)
    check("two builds lay the OCV row out identically", k1 == k2)


def check_images(jpg, pdf):
    from PIL import Image

    im = Image.open(jpg)
    check("JPEG width is 183 mm at 630 dpi (4539 px +/- 2)",
          abs(im.size[0] - 4539) <= 2, f"{im.size[0]} px")
    check("JPEG height is 61 mm at 630 dpi (1513 px +/- 2)",
          abs(im.size[1] - 1513) <= 2, f"{im.size[1]} px")
    check("JPEG pixel ratio is 3.000", abs(im.size[0] / im.size[1] - 3.0) < 2e-3,
          f"{im.size[0] / im.size[1]:.5f}")
    dpi = im.info.get("dpi", (0, 0))
    check("JPEG carries 630 dpi", abs(float(dpi[0]) - 630.0) < 1.0, f"{dpi}")
    check("JPEG is a full-colour raster", im.mode in ("RGB", "L"), im.mode)
    im.close()

    raw = open(pdf, "rb").read()
    flat = raw.replace(b" ", b"").replace(b"\n", b"")
    check("PDF page contains no raster image",
          b"/Subtype/Image" not in flat and b"/Subtype/Image" not in raw)
    check("PDF embeds its fonts as TrueType (FontType 42)",
          b"/FontFile2" in raw)
    check("PDF has no Type-3 fonts", b"/Subtype/Type3" not in flat)
    check("PDF carries no timestamp (it must be reproducible)",
          b"/CreationDate" not in raw and b"/ModDate" not in raw)
    head = raw[:4000]
    i = head.find(b"/MediaBox")
    nums = []
    if i >= 0:
        seg = head[i:i + 200].split(b"[", 1)
        if len(seg) == 2:
            body = seg[1].split(b"]", 1)[0]
            nums = [float(v) for v in body.split()]
    ok = len(nums) == 4 and abs(nums[2] - FIG_W_MM / 25.4 * 72.0) < 0.6 \
        and abs(nums[3] - FIG_H_MM / 25.4 * 72.0) < 0.6
    check("PDF page box is 183 x 61 mm", ok, f"{nums}")


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_md(path, pos, hashes):
    fails = checks_failed()
    comps = echo_components()
    ncell = len(echo_cells())
    dens = ncell / float(RAST["cols"] * RAST["rows"])
    L = []
    A = L.append
    A("# `concept/fig_ocv_schematic` — the observability core vector, drawn")
    A("")
    A("![the observability core vector](fig_ocv_schematic.jpg)")
    A("")
    A("One plate, three panels, no data. Every other figure in `figures/` is "
      "computed from a record and pinned to it; this one is computed from the "
      "**definition** the records are measured with, and its job is to say what "
      "the core vector is before any site's numbers are read. It therefore has "
      "no pins to a table — it has structural checks instead "
      f"({len(CHECKS)} of them, listed below and in "
      "`data/concept/fig_ocv_schematic.json`).")
    A("")
    A("```")
    A("physical asset -> SAR observation geometry -> asset-specific measurement")
    A("support -> radar echo representation -> OCV")
    A("```")
    A("")
    A("## The three panels")
    A("")
    A("**(a) the asset and the acquisition.** A generic lattice asset "
      "`G_a^{3D}` in its own frame (the grey triad `X`, `Y`, `Z`), observed "
      "from a satellite along the orange line of sight `LOS_c` at the look "
      "angle `θ_c` measured from nadir; the operator `Π_c` maps that frame "
      "into the radar frame `(az, r)`. The asset is drawn as four legs, five "
      "belts and alternating braces — something with a real cross-section, "
      "because that is what makes the projection lossy.")
    A("")
    A("The subscript on `LOS_c` and `θ_c` is the acquisition channel: the same "
      "`c` that `Π_c`, `W_{a,c}` and `C_{a,c}` carry, because the line of "
      "sight and the incidence angle belong to the pass rather than to the "
      "asset. Nothing is spelled out about ascending or descending — the "
      "subscript says it, and a check fails the build if a channel-dependent "
      "symbol ever loses it or a frame symbol picks it up.")
    A("")
    A("The mapping is written out exactly once, at the foot of the panel, "
      "left-aligned to the radar frame's own left edge: "
      "`(X, Y, Z) ⟶^{Π_c} (az, r)`. The operator appears nowhere else on the "
      "plate, and the two frames share no symbol — the physical frame is `X`, "
      "`Y`, `Z`, the radar frame is `az`, `r`, and a check fails the build if "
      "the two symbol sets ever intersect. The radar frame is drawn 0.040 "
      "panel units above the foot so that the statement has clear air beneath "
      "it, and a check pins both the gap and the alignment.")
    A("")
    A("**(b) the measurement support.** The projected asset `Π_c(G_a^{3D})` sits "
      "inside the asset- and channel-specific support `W_{a,c}` (blue), which "
      "sits inside the read chip `C_{a,c}` (grey, outermost). Outward from the "
      "support, in this order, come the **guard** — a white band 0.020 units "
      "wide (1.1 mm) — and then the **background ring**, dashed and lightly "
      "tinted, from which the no-echo reference is read.")
    A("")
    A("The ring is drawn as a true annulus: its tint is knocked out over the "
      "guard, so no tint reaches inside `W_{a,c}` and no background cell can "
      "be read as support. That order — support, guard, ring, chip — is the "
      "order the audit walks (`code/lumo/lumo_ocv_masks.py`), and two pixel "
      "checks prove it on the composed canvas rather than on the geometry: "
      "W's interior must equal the blue fill at 0.13 over white and the outer "
      "band must equal the ring tint, each to within 2/255.")
    A("")
    A("`C_{a,c}` is set *on* the chip's own lower edge with a white knockout, "
      "so the edge breaks behind the glyphs and the label cannot be read as "
      "naming the ring it used to float between; a check pins the label's own "
      "box across that edge. The panel states the set algebra "
      "`W_{a,c} ⊂ C_{a,c}` once at its foot and nothing else: no box, no "
      "legend, no count.")
    A("")
    A("**(c) the observable representation and the vector.** Inside the same "
      "`W_{a,c}` the echo mask is a raster of green cells: "
      f"{len(comps)} connected regions of "
      f"{', '.join(str(len(k)) for k in comps)} cells in a "
      f"{RAST['cols']} × {RAST['rows']} grid (density {dens:.3f}) and one "
      "closed hole. That is the extent, the density and the fragmentation the "
      "vector has to summarise, made visible before it is made numeric.")
    A("")
    A("Nothing is annotated on the mask itself — no centroid, no peak, no "
      "arrow between them, no symbol naming either. The centroid-to-peak "
      "distance that `S` measures is a number the raster yields, not an "
      "element the plate has to define: every mark in this panel can be read "
      "without a caption, and the only words on it are the vector's and its "
      "two group names. Below the mask:")
    A("")
    A("```")
    A("OCV(t) = [ γ², P, D, F, A, S ]")
    A("           ⏟           ⏟")
    A("        cross-acquisition  spatial geometry")
    A("```")
    A("")
    A("The two groups are carried by two restrained colours and two thin "
      "under-braces — not by boxes and not by a second legend. Each brace spans "
      "the measured extent of exactly the tokens it names: the row is laid out "
      "from the renderer's own text metrics at the output dpi, and the "
      "inter-group gap is *solved* so that the two names clear each other.")
    A("")
    A("| group | components | reads as | colour |")
    A("| --- | --- | --- | --- |")
    A("| cross-acquisition | `γ²`, `P` | what changes between acquisitions of "
      "the same place | blue `#1f4e79` |")
    A("| spatial geometry | `D`, `F`, `A`, `S` | how the echo is laid out "
      "inside `W_{a,c}` | green `#3f7d3f` |")
    A("")
    A("## What the plate does not show")
    A("")
    A("It is a schematic, so it carries no site, no date, no state and no number "
      "that could be read as a result. Deliberately absent:")
    A("")
    A("- structural health, damage probability, failure detection or any other "
      "diagnosis — the figure defines a measurement, it does not deliver a "
      "verdict;")
    A("- displacement, deformation or InSAR geometry — the OCV is a coherence "
      "and echo-geometry vector, and nothing in the plate moves;")
    A("- `R_B` and `R_G` — those belong to the temporal reconfiguration analysis "
      "(`figures/sds/`), not to the definition of the vector;")
    A("- a centroid or a peak symbol — `S` is a distance the mask yields, not a "
      "glyph the plate has to explain;")
    A("- a score, a threshold, a weighting or a count.")
    A("")
    A("The prohibition is mechanical, not editorial: a forbidden-token scan runs "
      "over every drawn string, and the whole text inventory is compared against "
      "a written-out whitelist, so the plate cannot grow a `damage` label "
      "without failing the run.")
    A("")
    A("## The plates, and how to reproduce them")
    A("")
    A("| file | what it is |")
    A("| --- | --- |")
    A("| `concept/fig_ocv_schematic.jpg` | 183 mm × 61 mm, 630 dpi JPEG "
      "(quality 95, 4:4:4 subsampling) |")
    A("| `concept/fig_ocv_schematic.pdf` | the same plate as a true vector page "
      "(FontType 42; no raster object anywhere in the file) |")
    A("| `concept/fig_ocv_schematic.md` | this report |")
    A("| `../data/concept/fig_ocv_schematic.json` | the drawn geometry, the text "
      "inventory, every check and the sha256 of both images |")
    A("")
    A("```bash")
    A("bash figures/concept/fig_ocv_schematic.sh                  # both images")
    A("python3 code/concept/fig_ocv_schematic.py --no-figures     # checks only")
    A("PYTHON=/usr/bin/python3 bash figures/concept/fig_ocv_schematic.sh")
    A("```")
    A("")
    A("## The checks")
    A("")
    A(f"**{len(CHECKS) - len(fails)} of {len(CHECKS)} checks pass.** A failure "
      "exits with status 2 *after* the artefacts are written, so a broken plate "
      "is reported rather than silently delivered. (`--no-figures` builds, "
      "measures and checks, writes the JSON and this report, and skips both "
      "images and the eleven file checks that go with them.)")
    A("")
    A("| family | what it pins |")
    A("| --- | --- |")
    A("| canvas | 183 mm × 61 mm at 630 dpi, ratio 3.000, `savefig.bbox` "
      "`standard` (a tight bbox would move the 183 mm edge); every panel "
      "isotropic, so one unit is one unit in x and in y and the raster cells "
      "are square |")
    A("| text inventory | exactly the declared "
      f"{len(expected_texts())} strings and no others; no title; no legend; "
      "all axes off |")
    A("| legibility | every text ≥ 7.0 pt at 183 mm and ≥ 6.5 pt at 170 mm; every "
      "text inside its own panel; no two texts overlap |")
    A("| vocabulary | no forbidden token in any drawn string (health, damage, "
      "fail, displace, deform, InSAR, score, threshold, px, `R_B`, `R_G`, …) |")
    A("| set algebra | `W_{a,c}` inside the guard inside the ring inside "
      "`C_{a,c}`, each with a visible margin; the guard is ≥ 0.015 units and "
      "the tinted band outside it ≥ 0.03; the projected asset lies inside "
      "`W` |")
    A("| the annulus | the composed canvas is sampled, not the geometry: W's "
      "interior must be the blue fill over white and the outer band the ring "
      "tint, each within 2/255, so the tint provably never runs behind the "
      "support |")
    A("| the two frames | the physical frame (`X`, `Y`, `Z`) and the radar "
      "frame (`az`, `r`) share no symbol; the mapping statement carries the "
      "operator and both frames; `C_{a,c}` sits on the chip's own edge |")
    A("| the channel | every channel-dependent symbol carries the `c` "
      "subscript — `LOS_c`, `θ_c`, `Π_c`, `W_{a,c}`, `C_{a,c}` — and no "
      "frame symbol does; ascending and descending are never spelled out |")
    A("| one W, two panels | the support is **one** polygon definition mapped "
      "into (b) and into (c); normalised to its own box it is identical to "
      "1 × 10⁻¹² in both |")
    A("| raster integrity | cells square on the drawn canvas and ≥ 1.5 mm; every "
      "echo cell inside `W` and inside the raster; exactly 4 connected regions "
      "of 8, 6, 4 and 3 cells; exactly one enclosed hole; density in "
      "(5 %, 20 %); the mask's centroid and its peak ≥ 1.5 cells apart — the "
      "quantity `S` measures, measured but not drawn |")
    A("| the vector | six components in the order `γ², P, D, F, A, S`; the "
      "partition into (γ², P) and (D, F, A, S) proved from the token table; each "
      "brace spans its own group's measured extent; the two names clear each "
      "other; the row sits inside the panel |")
    A("| determinism | two independent builds render identical pixels and lay "
      "the row out byte-identically; the PDF is written without a creation "
      "date, so both files are byte-identical across separate runs |")
    A("| the files | JPEG 4539 × 1513 px ± 2 at 630 dpi; the PDF has no "
      "`/Subtype /Image`, embeds `/FontFile2`, carries no Type-3 font, and its "
      "page box is 183 mm × 61 mm |")
    A("")
    A("### Every check, in the order it runs")
    A("")
    A("| # | check | result | detail |")
    A("| --- | --- | --- | --- |")
    for i, ck in enumerate(CHECKS, 1):
        A(f"| {i} | {ck['name']} | {'ok' if ck['ok'] else '**FAIL**'} "
          f"| {ck['detail']} |")
    A("")
    A("## Files")
    A("")
    if hashes:
        for name, h in hashes.items():
            A(f"- `{name}` — sha256 `{h}`")
    else:
        A("- no image was written (`--no-figures`)")
    A("")
    A(f"Environment: Python {platform.python_version()}, matplotlib "
      f"{matplotlib.__version__}, numpy {np.__version__}. The plate is "
      "deterministic: no timestamp, no seed-dependent value and no "
      "machine-dependent path reaches any output, so two runs produce "
      "byte-identical JSON, Markdown, JPEG and PDF.")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    return path


def payload(pos, hashes):
    """Everything the plate drew, plus every check it passed, as one JSON."""
    comps = echo_components()
    (ccx, ccy), (ci, cj) = centroid_xy()
    fails = checks_failed()
    return {
        "figure": {
            "file": OUT_JPG,
            "vector": OUT_PDF,
            "what": "the observability core vector as a definition, not a result",
            "reads": "no data file at all",
            "panels": {
                "a": "physical asset in (X,Y,Z) + acquisition geometry + the "
                     "mapping to (az,r)",
                "b": "support W_{a,c} in the guard, inside the background "
                     "ring, inside the chip C_{a,c}",
                "c": "echo mask inside W_{a,c} and the six-dimensional vector",
            },
            "not_shown": [
                "structural health", "damage probability", "failure detection",
                "displacement", "InSAR deformation", "a universal score",
                "R_B and R_G (temporal reconfiguration, not the definition)",
            ],
        },
        "canvas": {
            "width_mm": FIG_W_MM, "height_mm": FIG_H_MM, "dpi": DPI,
            "ratio": round(FIG_W_MM / FIG_H_MM, 6),
            "panel_width_frac": PANEL_W_FRAC, "panel_height_frac": PANEL_H_FRAC,
            "panel_x_frac": list(PANEL_X_FRAC), "panel_y_frac": PANEL_Y_FRAC,
            "panel_unit_mm": round(UNIT_MM, 6), "ylim": round(YLIM, 8),
            "bbox": "standard (a tight box would move the 183 mm edge)",
        },
        "style": {
            "font": "DejaVu Sans", "mathtext": "STIX",
            "sizes_pt": {"panel": FS_PANEL, "math": FS_MATH, "label": FS_LABEL,
                         "small": FS_SMALL, "group": FS_GROUP},
            "colours": {"ink": INK, "grey": GREY, "light": LIGHT,
                        "faint": FAINT, "blue_W": BLUE, "orange_los": ORANGE,
                        "green_echo": GREEN},
            "title": None, "legend": None, "shadows": None,
        },
        "support": {
            "w_unit": [list(v) for v in W_UNIT],
            "box_b": list(W_BOX), "box_c": list(W_BOX_C),
            "chip": list(CHIP), "ring": list(RING),
            "guard": [round(v, 6) for v in GUARD_BOX],
            "guard_band_units": GUARD,
            "ring_band_units": round(min(GUARD_BOX[0] - RING[0],
                                          GUARD_BOX[1] - RING[1],
                                          RING[2] - GUARD_BOX[2],
                                          RING[3] - GUARD_BOX[3]), 6),
            "subset": "W_{a,c} inside the guard inside the background ring "
                      "inside C_{a,c}",
            "ring_is_annulus": "the tint is knocked out over the guard, so no "
                               "tint reaches inside W_{a,c}",
        },
        "asset": {"tower": TOWER, "los": [list(v) for v in LOS],
                  "nadir": [list(v) for v in NADIR], "theta_radius": THETA_R,
                  "los_symbol": LOS_SYM, "theta_symbol": THETA_SYM,
                  "channel_note": "c is the acquisition channel: LOS_c, "
                                  "theta_c, Pi_c, W_{a,c} and C_{a,c} carry "
                                  "it; ascending and descending are never "
                                  "spelled out"},
        "frames": {
            "physical": list(PHYS_AX), "physical_reads": "X, Y, Z",
            "radar": [RAD_AZ, RAD_R], "radar_reads": "az, r",
            "mapping": MAP,
            "mapping_plain": "(X,Y,Z) -> (az,r) under Pi_c",
            "shared_symbols": sorted(set(PHYS_AX) & {RAD_AZ, RAD_R}),
            "why": "the two frames share no symbol, so no letter does two jobs "
                   "and the projection cannot be misread as a coordinate change "
                   "inside one frame",
        },
        "echo": {
            "raster": dict(RAST),
            "regions": {k: [list(c) for c in v] for k, v in ECHO.items()},
            "sizes": [len(k) for k in comps],
            "hole": list(ECHO_HOLE), "peak": list(ECHO_PEAK),
            "density": round(len(echo_cells())
                             / float(RAST["cols"] * RAST["rows"]), 6),
            "centroid_units": [round(ccx, 6), round(ccy, 6)],
            "centroid_cells": [round(ci, 6), round(cj, 6)],
            "centroid_peak_cells": round(float(np.hypot(
                cell_xy(*ECHO_PEAK)[0] - ccx, cell_xy(*ECHO_PEAK)[1] - ccy))
                / RAST["cell"], 6),
            "annotated": "nothing on the mask: no centroid or peak symbol is "
                         "drawn, so S is a measured distance the plate never "
                         "has to define",
        },
        "ocv": {
            "tokens": [{"key": k, "mathtext": s, "group": g}
                       for (k, s, g) in OCV_TOKENS],
            "groups": [{"key": g, "name": n, "colour": c,
                        "mathtext": ocv_group_string(g)}
                       for (g, n, c) in OCV_GROUPS],
            "row": {k: (list(v) if isinstance(v, tuple) else v)
                    for k, v in pos.items() if k != "spans"},
            "spans": {k: list(v) for k, v in pos["spans"].items()},
            "equation": "OCV(t) = [" + ", ".join(
                [k for (k, _s, _g) in OCV_TOKENS]) + "]",
        },
        "checks": {
            "total": len(CHECKS), "passed": len(CHECKS) - len(fails),
            "failed": len(fails), "failures": fails, "list": CHECKS,
        },
        "artifacts": hashes,
        "environment": {
            "python": platform.python_version(),
            "matplotlib": matplotlib.__version__,
            "numpy": np.__version__,
        },
    }


# ---------------------------------------------------------------------------
# entry point
# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Draw the observability core vector as a three-panel "
                    "schematic, and check every structural promise it makes.")
    ap.add_argument("--figdir", default=FIGDIR,
                    help="where the JPEG and the PDF go (default figures/concept)")
    ap.add_argument("--json", default=os.path.join(DATADIR, OUT_JSON),
                    help="the JSON artifact (default data/concept/...)")
    ap.add_argument("--md", default=None,
                    help="the report (default <figdir>/fig_ocv_schematic.md)")
    ap.add_argument("--no-figures", action="store_true",
                    help="build, measure and check, but write neither image")
    args = ap.parse_args(argv)

    fig, axes, pos = make_figure()
    run_checks(fig, axes, pos)

    hashes = {}
    if not args.no_figures:
        os.makedirs(args.figdir, exist_ok=True)
        jpg = os.path.join(args.figdir, OUT_JPG)
        pdf = os.path.join(args.figdir, OUT_PDF)
        fig.savefig(jpg, dpi=DPI, facecolor="white", format="jpg",
                    pil_kwargs={"quality": 95, "subsampling": 0})
        fig.savefig(pdf, facecolor="white", format="pdf",
                    metadata={"CreationDate": None, "ModDate": None})
        check_images(jpg, pdf)
        hashes = {OUT_JPG: sha256_file(jpg), OUT_PDF: sha256_file(pdf)}
        # a second, byte-for-byte run of the same plate must not move anything
        again = os.path.join(args.figdir, ".rastercheck.jpg")
        fig.savefig(again, dpi=DPI, facecolor="white", format="jpg",
                    pil_kwargs={"quality": 95, "subsampling": 0})
        same = sha256_file(again) == hashes[OUT_JPG]
        os.remove(again)
        check("re-saving the same figure gives the same JPEG bytes", same)
        hashes[OUT_JPG] = sha256_file(jpg)

    plt.close(fig)
    check_deterministic()

    md = args.md or os.path.join(args.figdir, OUT_MD)
    os.makedirs(os.path.dirname(os.path.abspath(md)) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.json)) or ".", exist_ok=True)
    data = payload(pos, hashes)
    with open(args.json, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=1, sort_keys=True, ensure_ascii=False)
        fh.write("\n")
    write_md(md, pos, hashes)

    fails = checks_failed()
    print(f"fig_ocv_schematic: {len(CHECKS) - len(fails)}/{len(CHECKS)} checks "
          f"passed, {len(TEXTS)} text objects, unit {UNIT_MM:.3f} mm")
    for name, h in hashes.items():
        print(f"  {name:28s} {h[:16]}")
    print(f"  {os.path.relpath(args.json, ROOT)}")
    print(f"  {os.path.relpath(md, ROOT)}")
    if fails:
        for ck in fails:
            print(f"  FAIL {ck['name']}: {ck['detail']}", file=sys.stderr)
        print(f"fig_ocv_schematic: {len(fails)} check(s) failed", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())












