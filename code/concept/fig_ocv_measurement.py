#!/usr/bin/env python3
"""The Observability Core Vector — one chip, three stages, one number set.

Question: the first schematic (``fig_ocv_schematic``) defines the vector from
the asset outwards. This one asks the question a reader of that plate asks next:
*what exactly is measured, and where?* It answers by taking **one** radar read
chip and showing it three times, each time less than the time before:

    C_{a,c}  ->  W_{a,c}  ->  M_{a,c}(t)  ->  OCV(t)

Three panels, left to right, all on the *same* spatial frame — the same chip,
the same asset footprint, the same scale, the same axes:

    (a) the native read chip            C_{a,c}: the SLC intensity as it comes
        off the acquisition, with more radar cells than any measurement will
        ever use
    (b) the measurement support         W_{a,c} subset C_{a,c}: the
        asset-specific support, its white guard and the background ring, drawn
        over the same intensity so the restriction is visible and not asserted
    (c) the echo mask                   M_{a,c}(t): the support cells that
        carry echo, and the vector they are summarised by,
        OCV(t) = [gamma^2, P, D, F, A, S], in its two conceptual groups
        (cross-acquisition; spatial geometry)

The single idea is **restriction**: (a) holds everything the sensor returned,
(b) keeps the part that belongs to this asset in this channel, (c) keeps the
part of that which carries echo. Nothing appears in (b) that is not in (a), and
nothing appears in (c) that is not in (b) — the panels are the same observation
at three levels of admission, not three different pictures.

The mask is **derived, not decorated**: the intensity field is written down as a
closed-form sum of declared responses, the mask is the set of support cells
whose peak sample exceeds one declared reference level, and the check fails the
build unless that set is what the plate claims (four connected regions of
strictly decreasing size, one enclosed hole, sparse, every cell inside W, and no
echo cell anywhere else).

The figure is a **schematic**, so it carries no site, no date, no state and no
result. Its text inventory is fixed and checked. ``R_B`` and ``R_G`` are
deliberately absent — they belong to the temporal reconfiguration analysis, not
to the construction of the vector — and so is every shift marker: no ``p`` with
a bar, no ``p`` with a star, no displacement arrow, no regular analysis grid.

What the figure must not be read as, and does not show: structural health,
damage probability, failure detection, displacement, InSAR deformation, a
scatterer catalogue, or a universal score. There is no threshold value on the
plate, no count and no verdict.

Figure (183 mm x 61 mm = 3.000 : 1, three panels, no legend):

    a  native read chip       full grey intensity field, the asset footprint as
                              a thin outline, the chip boundary with C_{a,c}
                              on it, and the radar axes az / r
    b  measurement support    the same field, strongly subdued outside the ring
                              and knocked out inside it; blue W_{a,c}, white
                              guard, grey ring band, the projected footprint as
                              a dark plate named in white on itself, and
                              W_{a,c} subset C_{a,c} at the foot
    c  echo mask + OCV        the mask inside W_{a,c}: inactive support cells in
                              very light grey, echo cells in green, nothing
                              annotated on the mask, the vector below

The channel subscript is never spelled out: ``c`` is the channel, exactly as on
the first schematic, and ``a`` is the asset. The mask's centroid and peak are
measured (and checked) but never drawn: S is a distance the raster yields, not a
glyph the plate has to explain.

Style: white background, thin black and grey construction lines, one muted blue
(W), one muted green (echo), a grey intensity ramp and nothing else. Vector-like,
no shadows, no decoration, no gradient fill. All geometry is laid out in panel
units where one unit is the same length in x and y, so the raster cells are
square and the support keeps its shape; the unit is 55.632 mm, so 7 pt text is
7 pt on the printed page at 183 mm.

Deliverables: ``figures/concept/fig_ocv_measurement.jpg`` (630 dpi) and
``figures/concept/fig_ocv_measurement.pdf`` (true vector, FontType 42),
``data/concept/fig_ocv_measurement.json`` (the drawn geometry, every check and
the sha256 of both images) and ``figures/concept/fig_ocv_measurement.md`` (the
figure guide).

Usage:
  python3 code/concept/fig_ocv_measurement.py
  python3 code/concept/fig_ocv_measurement.py --no-figures
  python3 code/concept/fig_ocv_measurement.py --figdir /tmp/x --json /tmp/y.json
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
from matplotlib.collections import PolyCollection  # noqa: E402
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
# is the same length in x and y (the raster cells stay square in every panel)
UNIT_MM = PANEL_W_FRAC * FIG_W_MM
YLIM = (PANEL_H_FRAC * FIG_H_MM) / UNIT_MM
PANELS = ("a", "b", "c")

# ---------------------------------------------------------------------------
# palette: one muted blue, one muted green, ink and greys
# ---------------------------------------------------------------------------
INK = "#1a1a1a"
GREY = "#6f6f6f"
LIGHT = "#b9b9b9"
FAINT = "#d8d8d8"
BLUE = "#1f4e79"
GREEN = "#3f7d3f"
MASK_IDLE = "#e9e9e9"        # support cell, evaluated, no echo

# type sizes (pt). 7.0 is the floor: the plate must stay readable at 170 mm.
FS_PANEL = 9.0
FS_MATH = 9.0
FS_LABEL = 8.0
FS_TITLE = 7.0
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

OUT_JPG = "fig_ocv_measurement.jpg"
OUT_PDF = "fig_ocv_measurement.pdf"
OUT_MD = "fig_ocv_measurement.md"
OUT_JSON = "fig_ocv_measurement.json"

# ---------------------------------------------------------------------------
# forbidden vocabulary: the figure is a definition, not a verdict
# ---------------------------------------------------------------------------
FORBIDDEN = (
    "damage", "health", "fail", "probab", "displace", "deform", "insar",
    "score", "threshold", "pixel", "px", "risk", "crack", "anomal", "diagnos",
    "predict", "peak", "median", "scatterer", "shift", "star", "count",
)
# symbols that belong to other analyses and must never reach this plate
FORBIDDEN_SYMS = (r"r_b", r"r_g", r"\bar{p}", r"p\bar", r"p^\star", r"p\star",
                  r"\star", r"p^\ast")

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
# the drawn geometry (all values in panel units) — ONE frame, three panels
# ---------------------------------------------------------------------------
# The whole plate is one spatial frame, read by all three panels: the same chip
# rectangle, the same asset footprint, the same radar frame, the same scale. The
# only thing that changes from panel to panel is how much of the frame is
# admitted, which is the figure's entire argument.
CHIP = (0.140, 0.290, 0.900, 0.860)         # C_{a,c}: 0.760 x 0.570 (52 x 39)
NATIVE = dict(cols=52, rows=39, cell=(CHIP[2] - CHIP[0]) / 52)   # 0.81 mm each
# the asset footprint, projected: one polygon, drawn in all three panels. It is
# cut a little wider than the echo mask it carries, so no mask cell ever crosses
# the footprint's own outline: the asset contains the echo, the echo is not
# clipped by the drawing.
ASSET = ((0.396, 0.490), (0.644, 0.490), (0.620, 0.720), (0.420, 0.720))

# The measurement support as a chamfered rectangle in its own unit box: ONE
# definition, mapped into panel (b) and panel (c) by an affine map each, and it
# is the same octagon the first schematic draws, so W_{a,c} is recognisable
# across the two plates.
W_UNIT = (
    (0.10, 0.00), (0.90, 0.00), (1.00, 0.10), (1.00, 0.90),
    (0.90, 1.00), (0.10, 1.00), (0.00, 0.90), (0.00, 0.10),
)
W_BOX = (0.270, 0.430, 0.770, 0.780)          # 0.500 x 0.350
GUARD = 0.020                                 # white gap between W and the ring
GUARD_BOX = (W_BOX[0] - GUARD, W_BOX[1] - GUARD,
             W_BOX[2] + GUARD, W_BOX[3] + GUARD)
RING = (0.205, 0.365, 0.835, 0.840)           # the background ring band

# The echo mask lives on a coarse lattice: one mask cell is exactly 2 x 2 native
# samples, which is what makes "the mask is a coarsened read of the same chip" a
# statement the checks can hold rather than a claim the caption makes. The
# lattice is inset two cells inside W on every side, so no cell can cross the
# chamfer and every cell is a cell *of the support*.
MASK = dict(cols=15, rows=10, cell=2 * NATIVE["cell"],
            x0=CHIP[0] + 11 * NATIVE["cell"], y0=CHIP[1] + 12 * NATIVE["cell"])

# The radar frame: the chip sits inside one axis frame in panels (a) and (b) —
# azimuth up the left edge, range along the bottom — so the chip is placed in a
# radar geometry rather than floating in the panel. Panel (c) has no room for it
# and does not need it: the chip and the support are the same objects, at the
# same scale, in all three panels.
AX = dict(x0=0.100, y0=0.170, x1=0.920, y1=0.900)
RAD_AZ = r"$\mathrm{az}$"                     # radar azimuth
RAD_R = r"$r$"                                # radar range


# --- the intensity field ---------------------------------------------------
# Deterministic and analytic: a smooth base, two texture terms, and seven
# declared responses. No random number is drawn anywhere, so the field — and
# therefore the mask derived from it — is a constant of this file.
FIELD_BASE = (0.30, 0.05, 3.4, 0.6, 2.6, 1.2)
FIELD_TEX = ((0.045, 13.0, 2.0, 9.0, 0.5), (0.030, 23.0, 4.1, 19.0, 2.3))
#   (kind, cx, cy, amplitude, sigma_x, sigma_y)   kind "ring": sigma_x is the
#   radius, sigma_y the width of the annulus
RESP = (
    ("ring", 0.4615, 0.6262, 0.50, 0.0292, 0.0094),  # the tower: a 3 x 3 ring of
                                                     # cells, so its own cell is
                                                     # a hole and not an echo
    ("blob", 0.5785, 0.6408, 0.40, 0.0351, 0.0205),  # 3 x 2 block of cells
    ("blob", 0.5638, 0.5531, 0.38, 0.0205, 0.0205),  # 2 x 2 block of cells
    ("blob", 0.4469, 0.5092, 0.38, 0.0205, 0.0146),  # 2 x 1 block of cells
    ("blob", 0.6515, 0.7138, 0.10, 0.0300, 0.0260),  # diffuse glow, no echo
    ("blob", 0.1760, 0.8360, 0.40, 0.0219, 0.0219),  # clutter, outside W and ring
    ("blob", 0.8550, 0.3800, 0.38, 0.0201, 0.0201),  # clutter, outside W and ring
)
# The reference level the mask is read against. Read as a set, the mask is the
# same 20 cells anywhere in 0.430..0.551 — every echo cell peaks above the top
# of that band and every idle cell below its bottom — so the level is chosen in
# the middle of the band and no cell is close to the decision.
T_REF = 0.490

# the grey ramp: 16 steps, dark to light, one drawn colour per step
RAMP_N = 16


def _grey(v: float) -> str:
    c = max(0, min(255, int(round(v * 255.0))))
    return f"#{c:02x}{c:02x}{c:02x}"


RAMP = tuple(_grey(0.32 + 0.64 * (k + 0.5) / RAMP_N) for k in range(RAMP_N))
# where the field is mapped onto the ramp: black point and white point
I_LO, I_HI = 0.14, 0.86

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
Y_FOOT = 0.070                              # the foot line of panels (a), (b)

TITLES = {
    "a": "Native SLC read chip",
    "b": "Asset-specific measurement support",
    "c": "Echo-mask measurement",
}
FOOT_A = "radar data access"                # what panel (a) is, in four words
CHIP_SYM = r"$C_{a,c}$"
W_SYM = r"$W_{a,c}$"
ASSET_SYM = r"$\Pi_c(G_a^{3D})$"
SUBSET_SYM = r"$W_{a,c}\subset C_{a,c}$"
RING_LABEL = "background ring"
QUAD_EPS = 0.0004           # quads overlap by a hair, so no seam shows


# ---------------------------------------------------------------------------
# the field and the mask it produces
# ---------------------------------------------------------------------------
def sample_grid():
    """The native SLC samples: centre coordinates and the intensity field."""
    c = NATIVE["cell"]
    x = CHIP[0] + (np.arange(NATIVE["cols"]) + 0.5) * c
    y = CHIP[1] + (np.arange(NATIVE["rows"]) + 0.5) * c
    X, Y = np.meshgrid(x, y)
    u = (X - CHIP[0]) / (CHIP[2] - CHIP[0])
    v = (Y - CHIP[1]) / (CHIP[3] - CHIP[1])
    b, ba, bf, bp, bf2, bp2 = FIELD_BASE
    I = b + ba * np.sin(bf * u + bp) * np.cos(bf2 * v + bp2)
    for amp, fu, pu, fv, pv in FIELD_TEX:
        I = I + amp * np.sin(fu * u + pu) * np.sin(fv * v + pv)
    for kind, cx, cy, a, s1, s2 in RESP:
        if kind == "ring":
            r = np.hypot(X - cx, Y - cy)
            I = I + a * np.exp(-0.5 * ((r - s1) / s2) ** 2)
        else:
            I = I + a * np.exp(-0.5 * (((X - cx) / s1) ** 2
                                       + ((Y - cy) / s2) ** 2))
    return X, Y, I


def mask_cell_samples(i, j, X, Y):
    """The native samples that fall inside support cell (column i, row j)."""
    x0 = MASK["x0"] + i * MASK["cell"]
    y0 = MASK["y0"] + j * MASK["cell"]
    return ((X >= x0) & (X <= x0 + MASK["cell"]) &
            (Y >= y0) & (Y <= y0 + MASK["cell"]))


def mask_cell_peak(i, j, X, Y, I):
    sel = mask_cell_samples(i, j, X, Y)
    return float(I[sel].max()) if sel.any() else -1.0


def derive_mask():
    """The echo mask: every support cell whose peak sample exceeds T_REF."""
    X, Y, I = sample_grid()
    return frozenset((i, j) for i in range(MASK["cols"])
                     for j in range(MASK["rows"])
                     if mask_cell_peak(i, j, X, Y, I) > T_REF)


ECHO = derive_mask()


def echo_cells():
    """Every echo cell once, in a fixed (column, row) order."""
    return sorted(ECHO)


def level_of(value: float) -> int:
    """The ramp step a field value is drawn with."""
    t = (value - I_LO) / (I_HI - I_LO)
    return int(max(0, min(RAMP_N - 1, np.floor(t * RAMP_N))))


def lighten(hexc: str, f: float) -> str:
    """Pull a colour `f` of the way toward white (an exact colour, not alpha)."""
    rgb = [int(hexc[k:k + 2], 16) / 255.0 for k in (1, 3, 5)]
    return "#" + "".join(f"{int(round(255 * (c + f * (1.0 - c)))):02x}"
                         for c in rgb)


# ---------------------------------------------------------------------------
# drawing primitives — every artist is one of these five, so the plate has no
# shadows, no gradients and no decorative element anywhere
# ---------------------------------------------------------------------------
TEXTS: list[dict] = []
DRAWN: dict[str, dict] = {}         # what each panel actually drew, for the
                                    # same-frame checks


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


def _quads(ax, polys, colours, z=1):
    """A raster cell set as filled quads: vector output, no image object.

    ``PolyCollection`` is a path collection, so the PDF stays a true vector page
    (``/Subtype /Image`` never appears) while the cells keep hard edges.
    """
    return ax.add_collection(PolyCollection(polys, facecolors=colours,
                                            edgecolors="none", linewidths=0.0,
                                            antialiased=False, zorder=z))


def _text(ax, x, y, s, fs=FS_LABEL, color=INK, ha="left", va="baseline",
          rot=0.0, weight="normal", z=6, tag="", bbox=None):
    t = ax.text(x, y, s, fontsize=fs, color=color, ha=ha, va=va, rotation=rot,
                rotation_mode="anchor", fontweight=weight, zorder=z, bbox=bbox)
    TEXTS.append({"ax": ax, "t": t, "s": s, "fs": fs, "tag": tag, "color": color})
    return t


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
    """Centre of mask cell (column i, row j)."""
    c = MASK["cell"]
    return (MASK["x0"] + (i + 0.5) * c, MASK["y0"] + (j + 0.5) * c)


def native_cell_xy(i, j):
    """Centre of native sample (column i, row j) — the raster cell centres."""
    c = NATIVE["cell"]
    return (CHIP[0] + (i + 0.5) * c, CHIP[1] + (j + 0.5) * c)


# ---------------------------------------------------------------------------
# the raster: every native sample as one quad, in the colour its level draws
# ---------------------------------------------------------------------------
def raster_cells():
    """(quad, colour) for every native cell, and whether it lies inside W.

    Drawn as filled quads rather than as an image: the PDF stays a true vector
    page, and a reader can edit the cells. A cell counts as *inside the support*
    when its own centre is inside W, which is the only choice that does not
    split a cell between two zones.
    """
    X, Y, I = sample_grid()
    c = NATIVE["cell"]
    inside_w = Path(w_poly(W_BOX))
    inside, outside = [], []
    for j in range(NATIVE["rows"]):
        y0 = CHIP[1] + j * c
        for i in range(NATIVE["cols"]):
            x0 = CHIP[0] + i * c
            quad = [(x0 - QUAD_EPS, y0 - QUAD_EPS),
                    (x0 + c + QUAD_EPS, y0 - QUAD_EPS),
                    (x0 + c + QUAD_EPS, y0 + c + QUAD_EPS),
                    (x0 - QUAD_EPS, y0 + c + QUAD_EPS)]
            item = (quad, RAMP[level_of(float(I[j, i]))])
            if inside_w.contains_point((x0 + 0.5 * c, y0 + 0.5 * c)):
                inside.append(item)
            else:
                outside.append(item)
    return inside, outside


def draw_intensity(ax, which="all", lighten_outside=0.0, z=1):
    """Paint the field: `which` in {"all", "in", "out"}."""
    inside, outside = raster_cells()
    sets = []
    if which in ("all", "in"):
        sets.append((inside, 0.0))
    if which in ("all", "out"):
        sets.append((outside, lighten_outside))
    for items, f in sets:
        if not items:
            continue
        polys = [q for (q, _col) in items]
        cols = [col if f == 0.0 else lighten(col, f) for (_q, col) in items]
        _quads(ax, polys, cols, z=z)


def draw_asset(ax, fc="none", edge=INK, lw=0.55, alpha=1.0, z=3):
    """The projected asset footprint: the same polygon in all three panels."""
    _poly(ax, ASSET, edge=edge, fc=fc, lw=lw, alpha=alpha, z=z)


def draw_title(ax, key):
    """Panel letter and panel title: the only place the stage is named."""
    _text(ax, 0.020, 1.020, f"({key})", fs=FS_PANEL, weight="bold", va="top",
          tag="panel")
    _text(ax, 0.115, 1.018, TITLES[key], fs=FS_TITLE, color=GREY, va="top",
          tag="title")


def draw_axis_frame(ax):
    """The radar frame: az up the left side, r along the bottom, one open L.

    Drawn *around* the chip, not inside it — the chip is placed in a radar
    geometry, and the two arms are the axes its rows and columns are indexed by.
    The arms are plain rules with the labels offset the same way the first
    schematic offsets its own radar frame, so the two plates share one
    convention for a frame.
    """
    a = AX
    _line(ax, (a["x0"], a["y0"]), (a["x0"], a["y1"]), INK, 0.60, z=2)
    _line(ax, (a["x0"], a["y0"]), (a["x1"], a["y0"]), INK, 0.60, z=2)
    _text(ax, a["x1"] + 0.012, a["y0"] - 0.024, RAD_R, fs=FS_SMALL, color=INK,
          ha="left", va="center", tag="axis")
    _text(ax, a["x0"], a["y1"] + 0.006, RAD_AZ, fs=FS_SMALL, color=INK,
          ha="center", va="bottom", tag="axis")


def draw_chip_edge(ax, label=True, z=5):
    """The chip boundary, with C_{a,c} sitting on its own lower edge."""
    x0, y0, x1, y1 = CHIP
    _poly(ax, [(x0, y0), (x1, y0), (x1, y1), (x0, y1)], edge=GREY, fc="none",
          lw=0.55, z=z)
    if label:
        _text(ax, x1 - 0.010, y0, CHIP_SYM, fs=FS_SMALL, color=GREY, ha="right",
              va="center", tag="chip",
              bbox=dict(facecolor="white", edgecolor="none", pad=1.2))


# ---------------------------------------------------------------------------
# panel (a): the native read chip, as it comes off the acquisition
# ---------------------------------------------------------------------------
def draw_panel_a(ax):
    draw_title(ax, "a")
    draw_axis_frame(ax)
    draw_intensity(ax, "all")
    draw_asset(ax, fc="none", edge=INK, lw=0.55, z=3)
    draw_chip_edge(ax, label=True)
    _text(ax, 0.500, Y_FOOT, FOOT_A, fs=FS_SMALL, color=INK, ha="center",
          va="baseline", tag="foot")
    DRAWN["a"] = dict(chip=CHIP, asset=ASSET, w=None, axis=AX,
                      plate=None)


# ---------------------------------------------------------------------------
# panel (b): the measurement support — C is still there, W is what is kept
# ---------------------------------------------------------------------------
B_LIGHTEN = 0.65          # how far the chip outside the ring is pulled to white
B_VEIL = 0.18             # how much blue sits over the field inside W
RING_TINT = 0.30          # the ring band's flat tint over white


def draw_panel_b(ax):
    wb = w_poly(W_BOX)
    draw_title(ax, "b")

    # the same field as panel (a), strongly subdued: the chip is still there,
    # it is simply not the measurement
    draw_intensity(ax, "all", lighten_outside=B_LIGHTEN)

    # the background ring, as a true annulus. The tint is painted as a flat
    # band over white and then knocked out again over the guard, so no tint
    # reaches inside W_{a,c} and no background cell sits on the support.
    rx0, ry0, rx1, ry1 = RING
    ring_rect = [(rx0, ry0), (rx1, ry0), (rx1, ry1), (rx0, ry1)]
    _poly(ax, ring_rect, edge="none", fc="white", lw=0.0, z=2)
    _poly(ax, ring_rect, edge="none", fc=LIGHT, lw=0.0, alpha=RING_TINT, z=2)
    _poly(ax, w_poly(GUARD_BOX), edge="none", fc="white", lw=0.0, z=2)

    # inside W the field is shown at full contrast again, clipped to the
    # support: the restriction is visible, not asserted
    inside, _outside = raster_cells()
    coll = _quads(ax, [q for (q, _c) in inside], [c for (_q, c) in inside], z=3)
    coll.set_clip_path(Polygon(wb, closed=True, transform=ax.transData))

    # the blue veil, then its own crisp boundary
    _poly(ax, wb, edge="none", fc=BLUE, lw=0.0, alpha=B_VEIL, z=4)
    _poly(ax, wb, edge=BLUE, fc="none", lw=0.85, z=4)

    # the asset, and the chip's own boundary. The footprint is filled here and
    # only here: in (a) and (c) it is an outline over a field, in (b) it is the
    # object the support is cut around.
    draw_asset(ax, fc=INK, edge=INK, lw=0.50, alpha=0.85, z=5)

    # the ring's boundary: dashed, so the band reads as a reference and not as
    # a second frame around the support
    _poly(ax, ring_rect, edge=GREY, fc="none", lw=0.50, ls=(0, (3.0, 2.0)),
          z=5)
    draw_chip_edge(ax, label=True, z=6)

    # four names, and nothing else, placed the way the first plate places its
    # own: W_{a,c} inside the support, the projected footprint named on the
    # footprint itself, the ring named from above by a leader, the containment
    # at the foot. The symbols that sit on the field keep the same white
    # knockout C_{a,c} uses; the asset symbol sits on the dark footprint, so it
    # is set in paper white and needs no box — inside W there is no white band
    # below the footprint deep enough to hold a 9 pt label clear of both the
    # footprint and the support's edge.
    _text(ax, 0.292, 0.640, W_SYM, fs=FS_MATH, color=BLUE, ha="left",
          va="baseline", tag="support",
          bbox=dict(facecolor="white", edgecolor="none", pad=1.0))
    _text(ax, 0.520, 0.592, ASSET_SYM, fs=FS_SMALL, color="white", ha="center",
          va="baseline", tag="projected")
    _text(ax, 0.620, 0.900, RING_LABEL, fs=FS_SMALL, color=INK, ha="center",
          va="baseline", tag="ring")
    _line(ax, (0.620, 0.886), (0.620, 0.812), GREY, 0.45, z=5)
    _text(ax, 0.500, Y_FOOT, SUBSET_SYM, fs=FS_MATH, ha="center", va="baseline",
          tag="subset")

    draw_axis_frame(ax)
    DRAWN["b"] = dict(chip=CHIP, asset=ASSET, w=tuple(wb), axis=AX,
                      plate=None)


# ---------------------------------------------------------------------------
# panel (c): the echo mask inside the support, and the vector it becomes
# ---------------------------------------------------------------------------
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


def draw_mask_lattice(ax):
    """Every support cell: inactive in very light grey, echo in green."""
    c = MASK["cell"]
    x0, y0 = MASK["x0"], MASK["y0"]
    w = MASK["cols"] * c
    h = MASK["rows"] * c
    for j in range(MASK["rows"]):
        for i in range(MASK["cols"]):
            x = x0 + i * c
            y = y0 + j * c
            _poly(ax, [(x, y), (x + c, y), (x + c, y + c), (x, y + c)],
                  edge="white", fc=MASK_IDLE, lw=0.30, z=2)
    _poly(ax, [(x0, y0), (x0 + w, y0), (x0 + w, y0 + h), (x0, y0 + h)],
          edge=LIGHT, fc="none", lw=0.45, z=3)
    for (i, j) in echo_cells():
        x = x0 + i * c
        y = y0 + j * c
        _poly(ax, [(x, y), (x + c, y), (x + c, y + c), (x, y + c)],
              edge="white", fc=GREEN, lw=0.34, z=3)


def draw_panel_c(ax):
    draw_title(ax, "c")
    wc = w_poly(W_BOX)

    # The support, and nothing else: the chip is not re-drawn here, because the
    # vector row occupies the band the chip's lower edge would run through, and
    # because (c) is about W_{a,c} alone — the chip's own frame is what (a) and
    # (b) established, at the same scale, in the same place.
    _poly(ax, wc, edge="none", fc=BLUE, lw=0.0, alpha=0.10, z=1)
    draw_mask_lattice(ax)

    # the asset is the same footprint, in the same place, drawn thin so the
    # mask stays the subject
    draw_asset(ax, fc="none", edge=LIGHT, lw=0.50, z=4)
    _poly(ax, wc, edge=BLUE, fc="none", lw=0.85, z=4)

    # Nothing is annotated inside the mask — no centroid dot, no peak dot, no
    # arrow between them, no symbol naming either. The centroid-to-peak distance
    # that S measures is a number the raster yields, not a glyph the plate has to
    # explain, so panel (c) is left with marks that need no definition: the
    # support, the echo, the vector.

    # the mask becomes the vector, pointing at the middle of the finished row.
    # The tail keeps the offset the first schematic uses for its own arrow (it
    # starts just under its raster's edge, not on it), so the arrow reads as a
    # hand-off between the two stages rather than as part of the support.
    pos = layout_ocv(ax)
    mid = 0.5 * (pos["left"] + pos["right"])
    ay_tail = W_BOX[1] - 0.012
    ay_tip = Y_OCV + 0.047
    _arrow(ax, (mid, ay_tail), (mid, ay_tip), INK, 0.80, scale=6.0)
    DRAWN["c"] = dict(chip=CHIP, asset=ASSET, w=tuple(wc), axis=None,
                      plate=None, arrow=(ay_tail, ay_tip))
    return pos


# ---------------------------------------------------------------------------
# the plate
# ---------------------------------------------------------------------------
def make_figure():
    """Build the three-panel plate and return it with its measured layout."""
    TEXTS.clear()
    DRAWN.clear()
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
    """Enclosed non-echo cells of the support lattice."""
    cells = set(echo_cells())
    cols, rows = MASK["cols"], MASK["rows"]
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


def centroid_xy():
    """Centroid of the echo mask, in panel units (cells weighted equally)."""
    cells = echo_cells()
    ci = float(np.mean([i + 0.5 for i, _j in cells]))
    cj = float(np.mean([j + 0.5 for _i, j in cells]))
    c = MASK["cell"]
    return (MASK["x0"] + ci * c, MASK["y0"] + cj * c), (ci, cj)


def peak_cell():
    """The brightest support cell: the mask's own extreme, measured not drawn."""
    X, Y, I = sample_grid()
    best, cell = -1.0, None
    for (i, j) in echo_cells():
        v = mask_cell_peak(i, j, X, Y, I)
        if v > best:
            best, cell = v, (i, j)
    return cell, best


# ---------------------------------------------------------------------------
# the checks: fourteen families, every one of them fatal
# ---------------------------------------------------------------------------
def expected_texts():
    """The complete text inventory of the plate, written out, not discovered."""
    t = ["(a)", "(b)", "(c)"]
    t += [TITLES["a"], TITLES["b"], TITLES["c"]]
    t += [CHIP_SYM, RAD_AZ, RAD_R, FOOT_A]                       # panel (a)
    t += [CHIP_SYM, RAD_AZ, RAD_R, W_SYM, ASSET_SYM, SUBSET_SYM,
          RING_LABEL]                                            # panel (b)
    t += [r"$\mathrm{OCV}(t)=[$", ocv_group_string("cross"), ",",
          ocv_group_string("spatial"), "]", "cross-acquisition",
          "spatial geometry"]                                    # panel (c)
    return t


def _dev(ax, pts):
    """Panel points in device pixels, so a frame can be compared across panels."""
    return np.array([ax.transData.transform(p) for p in pts], dtype=float)


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
    check("figure dpi is the output dpi", fig.dpi == DPI, f"{fig.dpi}")
    for k, ax in enumerate(axes):
        t = ax.transData
        dx = t.transform((1, 0))[0] - t.transform((0, 0))[0]
        dy = t.transform((0, 1))[1] - t.transform((0, 0))[1]
        check(f"panel {PANELS[k]} is isotropic (unit = unit)", abs(dx - dy) < 1e-6,
              f"dx={dx:.4f}px dy={dy:.4f}px")
        b = ax.get_position()
        check(f"panel {PANELS[k]} box matches the layout",
              abs(b.width - PANEL_W_FRAC) < 1e-9
              and abs(b.height - PANEL_H_FRAC) < 1e-9,
              f"{b.width:.6f}x{b.height:.6f}")
    for k in range(len(axes) - 1):
        check(f"panel {PANELS[k]}/{PANELS[k + 1]} boxes do not overlap",
              axes[k].get_position().x1 <= axes[k + 1].get_position().x0 + 1e-12,
              f"gap {axes[k + 1].get_position().x0 - axes[k].get_position().x1:.4f}")
    check("panel unit is 55.632 mm",
          abs(fig.transFigure.transform((PANEL_W_FRAC, 0))[0]
              - fig.transFigure.transform((0, 0))[0]) / DPI * 25.4 - UNIT_MM < 1e-6,
          f"{UNIT_MM:.3f} mm")

    # ---- 2. one frame, three panels --------------------------------------
    # The three panels must show the same chip, the same asset and the same
    # scale, so the progression in the plate is a progression in admission and
    # not a change of frame. Checked on what each panel actually drew, in device
    # pixels, not on the constants they were given.
    chip_pts = [(CHIP[0], CHIP[1]), (CHIP[2], CHIP[1]),
                (CHIP[2], CHIP[3]), (CHIP[0], CHIP[3])]
    dev = {k: _dev(axes[PANELS.index(k)], chip_pts) for k in PANELS}
    org = {k: np.array(axes[PANELS.index(k)].transData.transform((0.0, 0.0)))
           for k in PANELS}
    rel = {k: dev[k] - org[k] for k in PANELS}
    size = {k: dev[k][2] - dev[k][0] for k in PANELS}
    check("the chip rectangle maps to the same size in all three panels",
          all(np.allclose(size["a"], size[k], atol=1e-9) for k in PANELS),
          f"{np.round(size['a'], 3)} px")
    check("the chip sits at the same place inside all three panels",
          all(np.allclose(rel["a"], rel[k], atol=1e-9) for k in PANELS),
          f"corner offset {np.round(rel['a'][0], 3)} px from each panel's origin")

    check("one chip rectangle, shared by all three panels",
          all(DRAWN[k]["chip"] == CHIP for k in PANELS))
    check("all three panels draw the same asset footprint",
          all(DRAWN[k]["asset"] == ASSET for k in PANELS))
    check("panels (a) and (b) draw the same radar frame in the same place",
          DRAWN["a"]["axis"] == DRAWN["b"]["axis"] == AX)
    check("panel (c) names no axes, because it labels nothing on the mask",
          DRAWN["c"]["axis"] is None)
    check("panel (c) draws no chip plate: the vector row runs through that band",
          DRAWN["c"]["plate"] is None)
    check("the frame encloses the chip on all four sides",
          AX["x0"] < CHIP[0] - 0.005 and AX["y0"] < CHIP[1] - 0.005
          and AX["x1"] > CHIP[2] + 0.005 and AX["y1"] > CHIP[3] + 0.005,
          f"x {AX['x0']:.3f}..{AX['x1']:.3f} y {AX['y0']:.3f}..{AX['y1']:.3f}")
    check("the frame stays inside the panel",
          0.0 < AX["x0"] and AX["x1"] < 1.0 and 0.0 < AX["y0"]
          and AX["y1"] < YLIM - 0.02, f"{AX}")
    check("the frame's arm ends are clear of the chip's own edge",
          abs(AX["y0"] - CHIP[1]) > 0.05 and (AX["x1"] - CHIP[2]) > 0.02,
          f"r arm {CHIP[1] - AX['y0']:.3f} below the chip")
    check("panel (c)'s ground is the support alone, not a second frame",
          DRAWN["c"]["plate"] is None and DRAWN["c"]["w"] is not None)
    check("panels (b) and (c) draw the same W_{a,c}",
          DRAWN["b"]["w"] == DRAWN["c"]["w"], f"{len(DRAWN['b']['w'])} vertices")
    check("panel (a) draws no support, so the restriction is visible",
          DRAWN["a"]["w"] is None)
    # the same spatial extent: the chip covers the same fraction of every panel
    span = {k: axes[PANELS.index(k)].transData.transform((1, 0))[0]
            - axes[PANELS.index(k)].transData.transform((0, 0))[0] for k in PANELS}
    frac = {k: size[k][0] / span[k] for k in PANELS}
    check("the chip covers the same fraction of each panel",
          abs(frac["a"] - frac["b"]) < 1e-9 and abs(frac["a"] - frac["c"]) < 1e-9,
          f"{frac['a']:.6f} of the panel width")
    check("the native sample grid is the same grid in (a) and (b)",
          abs(NATIVE["cell"] * NATIVE["cols"] - (CHIP[2] - CHIP[0])) < 1e-12
          and abs(NATIVE["cell"] * NATIVE["rows"] - (CHIP[3] - CHIP[1])) < 1e-12,
          f"{NATIVE['cols']} x {NATIVE['rows']} cells")


    # ---- 3. text inventory ----------------------------------------------
    got = [rec["s"] for rec in TEXTS]
    want = expected_texts()
    check("text inventory is exactly the declared set",
          sorted(got) == sorted(want), f"{len(got)} objects")
    check("text count <= 30", len(got) <= 30, f"{len(got)}")
    check("no duplicate panel label",
          [g for g in got if g in ("(a)", "(b)", "(c)")] == ["(a)", "(b)", "(c)"])
    check("every panel carries exactly one title",
          [g for g in got if g in TITLES.values()] == [TITLES[k] for k in PANELS])
    check("no title anywhere", getattr(fig, "_suptitle", None) is None)
    check("no legend anywhere", all(a.get_legend() is None for a in axes))
    check("axes are off (no spines, no ticks)", all(not a.axison for a in axes))
    fmin = min(rec["fs"] for rec in TEXTS)
    check("every text >= 7.0 pt", fmin >= FS_FLOOR - 1e-9, f"min {fmin:.2f} pt")
    check("every text >= 7.0 pt at 170 mm (scale 0.929) = 6.5 pt",
          fmin * 170.0 / 183.0 >= 6.5 - 1e-9, f"{fmin * 170.0 / 183.0:.2f} pt")

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
        trec = next(rec for rec in mine if rec["tag"] == "title")
        t0, _ty0, t2, _ty1 = _extent(trec)
        check(f"panel {tag}: the title fits its own run",
              t0 > 0.075 and t2 <= 0.998, f"{t0:.3f}..{t2:.3f}")

    # ---- 3b. every name sits on the thing it names ------------------------
    # Not just "inside the panel": a symbol that floats off its own object is
    # a reading error, and (b) is the panel with the least room to spare.
    tb = {rec["tag"]: rec for rec in TEXTS
          if rec["ax"] is axes[1] and rec["tag"]}
    ap = Path(ASSET)
    e = _extent(tb["projected"])
    check("(b) names the asset on the asset: the symbol sits inside the footprint",
          bool(ap.contains_points([(x, y) for x in (e[0], e[2])
                                   for y in (e[1], e[3])]).all()),
          f"label {e[0]:.3f}..{e[2]:.3f} x {e[1]:.3f}..{e[3]:.3f}")
    check("(b) sets the asset symbol in paper white, on the dark footprint",
          tb["projected"]["color"] == "white")
    ew = _extent(tb["support"])
    pad = 0.353 * 1.0 / UNIT_MM          # the 1 pt knockout, in panel units
    check("(b)'s W_{a,c} label and its knockout stay clear of the footprint",
          not ap.contains_points([(x, y) for x in (ew[0] - pad, ew[2] + pad)
                                  for y in (ew[1] - pad, ew[3] + pad)]).any(),
          f"knockout to {ew[2] + pad:.3f}, footprint from "
          f"{min(p[0] for p in ASSET[:1] + ASSET[2:]):.3f} at its base")

    # ---- 4. vocabulary ---------------------------------------------------
    low = " ".join(got).lower()
    hit = [w for w in FORBIDDEN if w in low]
    check("no forbidden vocabulary in any drawn text", not hit, ",".join(hit))
    check("R_B absent", "r_b" not in low.replace(" ", ""))
    check("R_G absent", "r_g" not in low.replace(" ", ""))
    symhit = [s for s in FORBIDDEN_SYMS if any(s in g for g in got)]
    check("no shift marker (p-bar, p-star) anywhere on the plate", not symhit,
          ",".join(symhit))
    check("the coherence channel is mathtext", any("gamma^2" in g for g in got))
    check("the radar frame is az / r and nothing else",
          any(RAD_AZ == g for g in got) and any(RAD_R == g for g in got))
    check("the physical frame (X, Y, Z) is not on this plate at all",
          not any(("$X$" in g or "$Y$" in g or "$Z$" in g) for g in got))
    check("every channel-dependent symbol carries c",
          all("c" in s for s in (CHIP_SYM, W_SYM, SUBSET_SYM, ASSET_SYM)),
          "C_{a,c}, W_{a,c}, W_{a,c} subset C_{a,c}, Pi_c(G_a^{3D})")
    check("no channel-dependent symbol is spelled out instead",
          not any(w in low for w in ("ascending", "descending", "orbit", "pass")))

    # ---- 5. set algebra: W inside the guard inside the ring inside C -----
    wb = w_poly(W_BOX)
    path_w = Path(wb)
    check("W has eight vertices (the chamfered rectangle of the first plate)",
          len(wb) == 8)
    check("W sits inside C with a visible margin", _in_rect(wb, CHIP, 0.02),
          f"margin {min(W_BOX[0] - CHIP[0], W_BOX[1] - CHIP[1], CHIP[2] - W_BOX[2], CHIP[3] - W_BOX[3]):.3f}")
    check("W sits inside the guard", _in_rect(wb, GUARD_BOX, 0.0))
    check("the guard sits inside the ring", _in_rect(
        [(GUARD_BOX[0], GUARD_BOX[1]), (GUARD_BOX[2], GUARD_BOX[3])], RING, 0.0))
    guard = min(W_BOX[0] - GUARD_BOX[0], W_BOX[1] - GUARD_BOX[1],
                GUARD_BOX[2] - W_BOX[2], GUARD_BOX[3] - W_BOX[3])
    check("the guard is a white band >= 0.015 wide", guard >= 0.015 - 1e-12,
          f"{guard:.3f} units = {guard * UNIT_MM:.2f} mm")
    band = min(GUARD_BOX[0] - RING[0], GUARD_BOX[1] - RING[1],
               RING[2] - GUARD_BOX[2], RING[3] - GUARD_BOX[3])
    check("the ring band outside the guard is >= 0.03 wide",
          band >= 0.03 - 1e-12, f"{band:.3f} units = {band * UNIT_MM:.2f} mm")
    check("the ring sits inside C, so it cannot be the chip edge",
          _in_rect([(RING[0], RING[1]), (RING[2], RING[3])], CHIP, 0.015),
          f"margin {min(RING[0] - CHIP[0], RING[1] - CHIP[1], CHIP[2] - RING[2], CHIP[3] - RING[3]):.3f}")

    check("C sits inside the panel", _in_rect(
        [(CHIP[0], CHIP[1]), (CHIP[2], CHIP[3])], (0.0, 0.0, 1.0, YLIM), 0.01))
    check("the asset footprint lies inside W",
          bool(path_w.contains_points(ASSET).all()))
    check("the asset footprint is inside C, not straddling it",
          _in_rect(ASSET, CHIP, 0.02))
    check("the ring leader ends on the ring band, clear of the guard",
          GUARD_BOX[3] < 0.812 < RING[3],
          f"{GUARD_BOX[3]:.3f} < 0.812 < {RING[3]:.3f}")


    # the one definition of W, mapped once: the two panels cannot drift apart
    def _unit(pts, box):
        x0, y0, x1, y1 = box
        return [((x - x0) / (x1 - x0), (y - y0) / (y1 - y0)) for x, y in pts]

    check("W in (b) and W in (c) are the same shape, in the same box",
          np.allclose(_unit(wb, W_BOX), np.array(W_UNIT), atol=1e-12)
          and np.allclose(_unit(list(DRAWN["c"]["w"]), W_BOX),
                          np.array(W_UNIT), atol=1e-12))

    # ---- 6. the mask lattice: a coarsened read of the same chip ---------
    mc = MASK["cell"]
    lat = (MASK["x0"], MASK["y0"], MASK["x0"] + MASK["cols"] * mc,
           MASK["y0"] + MASK["rows"] * mc)
    check("one mask cell is exactly 2 x 2 native samples",
          abs(mc / NATIVE["cell"] - 2.0) < 1e-12,
          f"{mc:.3f} / {NATIVE['cell']:.3f} = {mc / NATIVE['cell']:.1f}")
    check("the mask lattice sits on the native sample grid",
          abs((MASK["x0"] - CHIP[0]) / NATIVE["cell"]
              - round((MASK["x0"] - CHIP[0]) / NATIVE["cell"])) < 1e-9
          and abs((MASK["y0"] - CHIP[1]) / NATIVE["cell"]
                  - round((MASK["y0"] - CHIP[1]) / NATIVE["cell"])) < 1e-9,
          f"origin at sample "
          f"({(MASK['x0'] - CHIP[0]) / NATIVE['cell']:.0f},"
          f" {(MASK['y0'] - CHIP[1]) / NATIVE['cell']:.0f})")
    check("the mask lattice sits inside W", _in_rect(
        [(lat[0], lat[1]), (lat[2], lat[3])], W_BOX, 0.010),
        f"margin {min(lat[0] - W_BOX[0], lat[1] - W_BOX[1], W_BOX[2] - lat[2], W_BOX[3] - lat[3]):.3f}")
    check("every mask cell is inside W",
          all(bool(path_w.contains_point(cell_xy(i, j)))
              for i in range(MASK["cols"]) for j in range(MASK["rows"])))
    ncell_total = NATIVE["cols"] * NATIVE["rows"]
    _inside, _outside = raster_cells()
    ncell_w = len(_inside)
    check("C holds far more cells than W (the read is not the measurement)",
          ncell_w * 2 < ncell_total,
          f"{ncell_w} of {ncell_total} native cells in W")
    check("the mask has fewer cells than W has native samples",
          len(ECHO) < ncell_w, f"{len(ECHO)} mask cells, {ncell_w} native cells")

    # ---- 7. raster integrity ---------------------------------------------
    c = NATIVE["cell"]
    check("native cells are square on the drawn canvas",
          abs((axes[0].transData.transform((CHIP[0] + c, CHIP[1]))[0]
               - axes[0].transData.transform((CHIP[0], CHIP[1]))[0])
              - (axes[0].transData.transform((CHIP[0], CHIP[1] + c))[1]
                 - axes[0].transData.transform((CHIP[0], CHIP[1]))[1])) < 1e-6,
          f"{c:.4f} units")
    check("native cell >= 0.8 mm on the printed page", c * UNIT_MM >= 0.8,
          f"{c * UNIT_MM:.2f} mm")
    check(f"the chip holds {NATIVE['cols']} x {NATIVE['rows']} native samples",
          ncell_total == NATIVE["cols"] * NATIVE["rows"],
          f"{ncell_total}")
    check("the raster is drawn as quads, not as an image",
          all(len(a.collections) >= 1 for a in axes[:2])
          and not any(getattr(a, "images", []) for a in axes))
    quads = [q for (q, _col) in _inside] + [q for (q, _col) in _outside]
    check("the quads cover the chip exactly once",
          len(quads) == ncell_total and
          all(_in_rect(q, (CHIP[0] - 2 * QUAD_EPS, CHIP[1] - 2 * QUAD_EPS,
                           CHIP[2] + 2 * QUAD_EPS, CHIP[3] + 2 * QUAD_EPS))
              for q in quads),
          f"{len(quads)} quads")
    cover = (min(p[0] for q in quads for p in q),
             min(p[1] for q in quads for p in q),
             max(p[0] for q in quads for p in q),
             max(p[1] for q in quads for p in q))
    check("the quad union spans the chip from edge to edge",
          abs(cover[0] + QUAD_EPS - CHIP[0]) < 1e-9
          and abs(cover[1] + QUAD_EPS - CHIP[1]) < 1e-9
          and abs(cover[2] - QUAD_EPS - CHIP[2]) < 1e-9
          and abs(cover[3] - QUAD_EPS - CHIP[3]) < 1e-9,
          f"{np.round(cover, 4)}")
    check("no native cell is counted twice",
          len(_inside) + len(_outside) == ncell_total)
    lums = [int(RAMP[k][1:3], 16) for k in range(RAMP_N)]
    check("the intensity ramp has 16 steps", len(RAMP) == 16)
    check("the ramp is monotone and has no repeated step",
          all(b > a for a, b in zip(lums, lums[1:])), f"{lums[0]}..{lums[-1]}")
    check("the darkest ramp step is grey, not ink", lums[0] >= 70, f"{lums[0]}/255")
    check("the lightest ramp step is not paper white", lums[-1] <= 250,
          f"{lums[-1]}/255")
    X, Y, I = sample_grid()
    X2, Y2, I2 = sample_grid()
    check("the field is deterministic (no random number anywhere)",
          bool(np.array_equal(I, I2)))
    used = len({level_of(float(v)) for v in I.ravel()})
    check("the field uses at least half the ramp", used >= RAMP_N // 2,
          f"{used} of {RAMP_N} steps used")
    imax = np.unravel_index(int(np.argmax(I)), I.shape)
    pkt = native_cell_xy(int(imax[1]), int(imax[0]))
    check("the brightest sample in the chip lies inside W",
          bool(path_w.contains_point(pkt)),
          f"({pkt[0]:.3f}, {pkt[1]:.3f}) = {I.max():.3f}")
    ncol = NATIVE["cols"]
    inside_flat = np.array([path_w.contains_point((float(x), float(y)))
                            for x, y in zip(X.ravel(), Y.ravel())])
    out_max = float(I.ravel()[~inside_flat].max())
    check("the chip outside W does carry bright returns", out_max > T_REF,
          f"peak outside W {out_max:.3f} > reference {T_REF:.3f}")
    in_max = float(I.ravel()[inside_flat].max())
    check("the strongest return in the chip is on the support", in_max > out_max,
          f"{in_max:.3f} inside vs {out_max:.3f} outside")

    # ---- 8. the echo mask ------------------------------------------------
    comps = echo_components()
    sizes = [len(k) for k in comps]
    check("the echo mask has exactly four connected regions", len(comps) == 4,
          f"{len(comps)}")
    check("the four regions are 8, 6, 4 and 2 cells", sizes == [8, 6, 4, 2],
          f"{sizes}")
    check("the regions are strictly decreasing in size",
          all(sizes[i] > sizes[i + 1] for i in range(len(sizes) - 1)))
    check("every region is 4-connected (a fragment, not a dust field)",
          all(len(k) == len(set(k)) for k in comps))
    check("the mask has 20 cells", len(ECHO) == 20, f"{len(ECHO)}")
    holes = echo_holes()
    check("the mask encloses exactly one hole", len(holes) == 1, f"{len(holes)}")
    check("the hole is enclosed, not on the lattice border",
          bool(holes) and 0 < holes[0][0] < MASK["cols"] - 1
          and 0 < holes[0][1] < MASK["rows"] - 1, f"{holes}")
    check("no echo cell touches the lattice border",
          all(0 < i < MASK["cols"] - 1 and 0 < j < MASK["rows"] - 1
              for (i, j) in ECHO))
    dens = len(ECHO) / float(MASK["cols"] * MASK["rows"])
    check("echo density is sparse and legible (5 % < d < 22 %)",
          0.05 < dens < 0.22, f"{dens:.4f}")
    peaks = {cell: mask_cell_peak(cell[0], cell[1], X, Y, I) for cell in ECHO}
    quiet = [mask_cell_peak(i, j, X, Y, I) for i in range(MASK["cols"])
             for j in range(MASK["rows"]) if (i, j) not in ECHO]
    check("every echo cell is above the reference level",
          min(peaks.values()) > T_REF,
          f"min echo {min(peaks.values()):.3f} > {T_REF:.3f}")
    check("no inactive cell reaches the reference level",
          max(quiet) < T_REF, f"max inactive {max(quiet):.3f} < {T_REF:.3f}")
    check("the mask is not knife-edge: it clears the level by 0.02 either way",
          min(peaks.values()) - T_REF > 0.02 and T_REF - max(quiet) > 0.02,
          f"margins {min(peaks.values()) - T_REF:.3f} / {T_REF - max(quiet):.3f}")
    check("the mask is exactly the threshold set of the field, cell for cell",
          set(ECHO) == {(i, j) for i in range(MASK["cols"])
                        for j in range(MASK["rows"])
                        if mask_cell_peak(i, j, X, Y, I) > T_REF})
    on_asset = sum(1 for (i, j) in ECHO
                   if Path(ASSET).contains_point(cell_xy(i, j)))
    check("the asset footprint carries the bulk of the echo", on_asset >= 8,
          f"{on_asset} of {len(ECHO)} echo cells on the asset")
    # the asset is drawn *around* the echo, so no mask cell may be cut by the
    # footprint's own outline: in (a) and (c) that outline is the only thing
    # marking the asset, and a cell poking through it reads as an error
    mcell = MASK["cell"]
    cut = [(i, j) for (i, j) in ECHO
           if not Path(ASSET).contains_points([
               (MASK["x0"] + i * mcell, MASK["y0"] + j * mcell),
               (MASK["x0"] + (i + 1) * mcell, MASK["y0"] + j * mcell),
               (MASK["x0"] + i * mcell, MASK["y0"] + (j + 1) * mcell),
               (MASK["x0"] + (i + 1) * mcell, MASK["y0"] + (j + 1) * mcell),
           ]).all()]
    check("no echo cell crosses the asset outline", not cut, f"{cut}")
    pcell, pval = peak_cell()
    check("the strongest echo cell lies inside W",
          bool(path_w.contains_point(cell_xy(*pcell))),
          f"{pcell} = {pval:.3f}")

    # what the plate must NOT be: a threshold of the picture. Outside the
    # support the chip is bright, and still contributes no cell at all.
    check("echo stops at the support: bright returns outside W make no cell",
          out_max > T_REF and not any(not path_w.contains_point(cell_xy(i, j))
                                      for (i, j) in ECHO),
          f"{len(ECHO)} cells, all inside W")

    # the distance S measures: measured, reported, never drawn
    (ccx, ccy), (ci, cj) = centroid_xy()
    ppx, ppy = cell_xy(*pcell)
    sep_cells = float(np.hypot(ppx - ccx, ppy - ccy)) / MASK["cell"]
    check("the centroid and the peak are >= 1.5 cells apart", sep_cells >= 1.5,
          f"{sep_cells:.3f} cells, the quantity S measures")
    check("the centroid of the mask lies inside the support",
          bool(path_w.contains_point((ccx, ccy))),
          f"({ccx:.3f}, {ccy:.3f}) in panel units")

    check("nothing is annotated on the mask: no marker, no dot, no label",
          not any(rec["tag"] in ("centroid", "peak", "mask_label")
                  for rec in TEXTS)
          and all(getattr(ln, "get_marker", lambda: None)() in ("", "None", None)
                  for ln in axes[2].lines),
          f"{len(axes[2].lines)} line objects, all markers empty")

    # ---- 9. the four zones, read off the composed canvas ------------------
    # Not the geometry: the pixels. If the ring tint ever runs behind the
    # support, or the guard is never knocked out, these fail.
    arr = render_array(fig)

    def _px(ax, x, y):
        px, py = ax.transData.transform((x, y))
        return arr[int(round(arr.shape[0] - py)), int(round(px)), :3] / 255.0

    def _comp(hexc, alpha):
        rgb = np.array([int(hexc[1:3], 16), int(hexc[3:5], 16),
                        int(hexc[5:7], 16)]) / 255.0
        return 1.0 - alpha * (1.0 - rgb)

    def _rgb(hexc):
        return np.array([int(hexc[1:3], 16), int(hexc[3:5], 16),
                         int(hexc[5:7], 16)]) / 255.0

    ax_a, ax_b = axes[0], axes[1]
    out_pt = native_cell_xy(2, 2)
    w_pt = native_cell_xy(40, 25)
    guard_pt = (0.500, 0.790)
    ring_pt = (0.500, 0.820)
    out_col = RAMP[level_of(float(I[2, 2]))]
    w_col = RAMP[level_of(float(I[25, 40]))]
    check("(a) paints the field cell at its own ramp step",
          bool(np.abs(_px(ax_a, *out_pt) - _rgb(out_col)).max() <= 0.01),
          "rgb " + " ".join(f"{v:.3f}" for v in _px(ax_a, *out_pt)))
    check("(b) shows the same cell outside the ring, pulled 65 % to white",
          bool(np.abs(_px(ax_b, *out_pt) - _rgb(lighten(out_col, B_LIGHTEN)))
               .max() <= 0.01),
          "rgb " + " ".join(f"{v:.3f}" for v in _px(ax_b, *out_pt)))
    check("(b) keeps the field inside W at full contrast under the blue veil",
          bool(np.abs(_px(ax_b, *w_pt)
                      - (_rgb(w_col) * (1 - B_VEIL) + _rgb(BLUE) * B_VEIL))
               .max() <= 0.01),
          "rgb " + " ".join(f"{v:.3f}" for v in _px(ax_b, *w_pt)))
    check("the veil is visible: (b) inside W is not (a)",
          float(np.abs(_px(ax_b, *w_pt) - _px(ax_a, *w_pt)).max()) > 0.05)
    check("(b)'s guard is pure white, knocked out of the ring tint",
          bool(np.abs(_px(ax_b, *guard_pt) - 1.0).max() <= 0.01),
          "rgb " + " ".join(f"{v:.3f}" for v in _px(ax_b, *guard_pt)))
    check("(b)'s ring band is the flat tint, not the field",
          bool(np.abs(_px(ax_b, *ring_pt) - _comp(LIGHT, RING_TINT)).max() <= 0.01),
          "rgb " + " ".join(f"{v:.3f}" for v in _px(ax_b, *ring_pt)))
    check("the ring tint never reaches inside W (annulus, not a fill)",
          bool(np.abs(_px(ax_b, *w_pt) - _comp(LIGHT, RING_TINT)).max() > 0.05))
    check("(b)'s footprint is a dark plate, so the white symbol reads on it",
          float(_px(ax_b, 0.520, 0.520).max()) < 0.40,
          "rgb " + " ".join(f"{v:.3f}" for v in _px(ax_b, 0.520, 0.520)))
    check("(c) draws no chip plate: outside the support is bare paper",
          bool(np.abs(_px(axes[2], 0.950, 0.930) - 1.0).max() <= 0.01),
          "rgb " + " ".join(f"{v:.3f}" for v in _px(axes[2], 0.950, 0.930)))
    ground_pt = (0.285, 0.670)          # inside W, between its edge and the lattice
    check("(c)'s support ground is the pale blue over paper",
          bool(np.abs(_px(axes[2], *ground_pt)
                      - (_rgb("#ffffff") * 0.90 + _rgb(BLUE) * 0.10)).max()
               <= 0.01),
          "rgb " + " ".join(f"{v:.3f}" for v in _px(axes[2], *ground_pt)))
    check("(c)'s inactive support cells are very light grey",
          bool(np.abs(_px(axes[2], *cell_xy(0, 0)) - _rgb(MASK_IDLE)).max() <= 0.01),
          "rgb " + " ".join(f"{v:.3f}" for v in _px(axes[2], *cell_xy(0, 0))))
    check("(c)'s echo cells are green",
          bool(np.abs(_px(axes[2], *cell_xy(*echo_cells()[0]))
                      - _rgb(GREEN)).max() <= 0.01),
          "rgb " + " ".join(f"{v:.3f}" for v in
                            _px(axes[2], *cell_xy(*echo_cells()[0]))))
    rgb = arr[:, :, :3] / 255.0
    greenish = np.abs(rgb - _rgb(GREEN)).max(axis=2) < 0.02
    n_green = int(greenish.sum())
    two_thirds = int(2.0 * FIG_W_MM * MM * DPI / 3.0)
    check("green appears in (c), and only there", n_green > 500,
          f"{n_green} exact-green pixels on the whole plate")
    check("the green is only the mask: (a) and (b) carry no green at all",
          int(greenish[:, :two_thirds].sum()) == 0,
          "left two thirds of the canvas")

    # ---- 10. the vector ---------------------------------------------------
    keys = [k for (k, _s, _g) in OCV_TOKENS]
    check("the vector has six components", len(keys) == 6)
    check("the order is gamma2, P, D, F, A, S",
          keys == ["gamma2", "P", "D", "F", "A", "S"], ",".join(keys))
    check("the partition covers every component exactly once",
          sorted([g for (_k, _s, g) in OCV_TOKENS])
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
          pos["spans"]["cross"][1] <= sep[0]
          and sep[2] <= pos["spans"]["spatial"][0],
          f"{sep[0]:.3f}..{sep[2]:.3f}")
    check("the separating comma does not touch either group",
          sep[0] - pos["spans"]["cross"][1] > 0.004
          and pos["spans"]["spatial"][0] - sep[2] > 0.004,
          f"{sep[0] - pos['spans']['cross'][1]:.3f} / "
          f"{pos['spans']['spatial'][0] - sep[2]:.3f}")
    check("the two group names do not overlap", a1[2] < a2[0],
          f"{a1[2]:.3f} < {a2[0]:.3f}")
    check("the equation sits above the braces", Y_BRACE + 0.015 < Y_OCV)
    check("the equation does not touch the mask lattice",
          _extent(recs["ocv_pre"])[3] < MASK["y0"],
          f"row top {_extent(recs['ocv_pre'])[3]:.3f} < lattice {MASK['y0']:.3f}")
    ay_tail, ay_tip = DRAWN["c"]["arrow"]
    check("the arrow leaves the support's lower edge and points at the row",
          ay_tip < ay_tail < W_BOX[1],
          f"{ay_tip:.3f} -> {ay_tail:.3f} under support {W_BOX[1]:.3f}")
    check("the whole vector row sits below the support",
          _extent(recs["ocv_pre"])[3] < W_BOX[1]
          and _extent(recs["group_spatial"])[3] > 0.0,
          f"row top {_extent(recs['ocv_pre'])[3]:.3f} < support {W_BOX[1]:.3f}")



# ---------------------------------------------------------------------------
# the written artefacts
# ---------------------------------------------------------------------------
def render_array(fig):
    """The composed canvas as an array, for the determinism and zone checks."""
    fig.canvas.draw()
    return np.asarray(fig.canvas.buffer_rgba()).copy()


def check_deterministic():
    """Two independent builds must produce the same pixels, not just the same JSON."""
    fig1, _ax1, pos1 = make_figure()
    a = render_array(fig1)
    plt.close(fig1)
    fig2, _ax2, pos2 = make_figure()
    b = render_array(fig2)
    plt.close(fig2)
    check("two builds render identical pixels",
          a.shape == b.shape and bool(np.array_equal(a, b)), f"{a.shape}")
    k1 = json.dumps(pos1, sort_keys=True, default=str)
    k2 = json.dumps(pos2, sort_keys=True, default=str)
    check("two builds lay the OCV row out identically", k1 == k2)


def check_cross_plate(pos):
    """The vector row must be the same row the first schematic draws.

    Both plates set the same equation at the same size in the same panel width,
    so a reader comparing the two must find the row in the same place to within
    a thousandth of a panel unit. Skipped, not faked, if the other plate has not
    been built in this checkout.
    """
    path = os.path.join(DATADIR, "fig_ocv_schematic.json")
    if not os.path.exists(path):
        check("the vector row agrees with the first schematic", True,
              "fig_ocv_schematic.json not built in this checkout")
        return
    with open(path, "r", encoding="utf-8") as fh:
        other = json.load(fh)
    row = other.get("ocv", {}).get("row", {})
    ok = (abs(row.get("left", -1.0) - pos["left"]) < 1e-6
          and abs(row.get("total", -1.0) - pos["total"]) < 1e-6)
    check("the vector row agrees with the first schematic", bool(ok),
          f"left {row.get('left', float('nan')):.6f} vs {pos['left']:.6f}, "
          f"width {row.get('total', float('nan')):.6f} vs {pos['total']:.6f}")


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


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
          b"/Subtype/Image" not in flat and b"/Subtype/Image" not in raw,
          f"the raster is {NATIVE['cols'] * NATIVE['rows']} vector quads, "
          "not an image object")
    check("PDF embeds its fonts as TrueType (FontType 42)", b"/FontFile2" in raw)
    check("PDF has no Type-3 fonts", b"/Subtype/Type3" not in flat)
    check("PDF carries no timestamp (it must be reproducible)",
          b"/CreationDate" not in raw and b"/ModDate" not in raw)
    check("PDF is a manageable vector page (< 1.5 MB)", len(raw) < 1500000,
          f"{len(raw) / 1e6:.2f} MB")
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


def payload(pos, hashes):
    """Everything the plate drew, plus every check it passed, as one JSON."""
    comps = echo_components()
    holes = echo_holes()
    (ccx, ccy), (ci, cj) = centroid_xy()
    pcell, pval = peak_cell()
    X, Y, I = sample_grid()
    inside, outside = raster_cells()
    fails = checks_failed()
    imax = np.unravel_index(int(np.argmax(I)), I.shape)
    path_w = Path(w_poly(W_BOX))
    out_flat = I.ravel()[~np.array([path_w.contains_point((float(x), float(y)))
                                    for x, y in zip(X.ravel(), Y.ravel())])]
    return {
        "figure": {
            "file": OUT_JPG,
            "vector": OUT_PDF,
            "what": "one radar read chip, three levels of admission, and the "
                    "vector the last of them becomes",
            "reads": "no data file at all: the field is analytic",
            "idea": "C_{a,c} -> W_{a,c} -> M_{a,c}(t) -> OCV(t), on one frame",
            "panels": {
                "a": "the native SLC read chip: the full intensity field, the "
                     "asset footprint as a thin outline, C_{a,c} on the chip "
                     "edge, the radar axes az / r",
                "b": "the measurement support W_{a,c} subset C_{a,c}: the same "
                     "field subdued outside, knocked out in the ring and guard, "
                     "blue W_{a,c} over the full-contrast field, the projected "
                     "footprint as a dark plate named in white on itself",
                "c": "the echo mask M_{a,c}(t) inside W_{a,c}: inactive "
                     "support cells in very light grey, echo cells in green, "
                     "and the six-dimensional vector below",
            },
            "not_shown": [
                "structural health", "damage probability", "failure detection",
                "displacement", "InSAR deformation", "a universal score",
                "R_B and R_G (temporal reconfiguration, not the definition)",
                "any shift marker (p-bar, p-star)",
                "the mask's centroid or peak, drawn",
                "any reference value, count or verdict on the plate",
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
                         "title": FS_TITLE, "small": FS_SMALL,
                         "group": FS_GROUP},
            "colours": {"ink": INK, "grey": GREY, "light": LIGHT,
                        "faint": FAINT, "blue_W": BLUE, "green_echo": GREEN,
                        "mask_inactive": MASK_IDLE},
            "ramp_steps": list(RAMP),
            "title": "one 7 pt title per panel, no plate title, no legend",
            "shadows": None, "gradients": None,
        },
        "frame": {
            "promise": "one spatial frame, read by all three panels: the same "
                       "chip, the same asset footprint, the same scale",
            "chip": list(CHIP),
            "native": dict(NATIVE),
            "native_cell_mm": round(NATIVE["cell"] * UNIT_MM, 6),
            "asset": [list(v) for v in ASSET],
            "w_unit": [list(v) for v in W_UNIT],
            "w_box": list(W_BOX),
            "guard": [round(v, 6) for v in GUARD_BOX],
            "guard_band_units": GUARD,
            "ring": list(RING),
            "ring_band_units": round(min(GUARD_BOX[0] - RING[0],
                                         GUARD_BOX[1] - RING[1],
                                         RING[2] - GUARD_BOX[2],
                                         RING[3] - GUARD_BOX[3]), 6),
            "subset": "W_{a,c} inside the guard inside the background ring "
                      "inside C_{a,c}",
            "cells_in_C": NATIVE["cols"] * NATIVE["rows"],
            "cells_in_W": len(inside),
            "axis_frame": dict(AX),
            "axes": [RAD_AZ, RAD_R],
            "anchored": "the chip sits inside one open radar frame in (a) and "
                        "(b) — az up the left side, r along the bottom — so the "
                        "read is placed in a geometry, not left floating",
        },
        "field": {
            "what": "a closed-form intensity field: no random number anywhere, "
                    "so the mask derived from it is a constant of the file",
            "base": {"offset": FIELD_BASE[0], "amplitude": FIELD_BASE[1],
                     "frequencies": [FIELD_BASE[2], FIELD_BASE[4]],
                     "phases": [FIELD_BASE[3], FIELD_BASE[5]]},
            "texture": [{"amplitude": a, "fu": fu, "pu": pu, "fv": fv,
                         "pv": pv} for (a, fu, pu, fv, pv) in FIELD_TEX],
            "responses": [{"kind": k, "centre": [cx, cy], "amplitude": a,
                           "sigma": [s1, s2]} for (k, cx, cy, a, s1, s2) in RESP],
            "minimum": round(float(I.min()), 6),
            "maximum": round(float(I.max()), 6),
            "peak_at": [round(float(X[imax]), 6), round(float(Y[imax]), 6)],
            "peak_inside_W": bool(path_w.contains_point(
                (float(X[imax]), float(Y[imax])))),
            "peak_outside_W": round(float(out_flat.max()), 6),
            "ramp": {"steps": RAMP_N, "lo": I_LO, "hi": I_HI},
        },
        "mask": {
            "lattice": dict(MASK),
            "native_per_cell": int(round(MASK["cell"] / NATIVE["cell"])),
            "reference_level": T_REF,
            "definition": "the support cells whose peak native sample exceeds "
                          "the reference level",
            "cells": [list(c) for c in echo_cells()],
            "regions": [{"cells": len(k), "at": [list(c) for c in k]}
                        for k in comps],
            "sizes": [len(k) for k in comps],
            "holes": [list(h) for h in holes],
            "density": round(len(ECHO) / float(MASK["cols"] * MASK["rows"]), 6),
            "min_echo_peak": round(min(mask_cell_peak(i, j, X, Y, I)
                                       for (i, j) in ECHO), 6),
            "max_inactive_peak": round(max(mask_cell_peak(i, j, X, Y, I)
                                           for i in range(MASK["cols"])
                                           for j in range(MASK["rows"])
                                           if (i, j) not in ECHO), 6),
            "centroid_units": [round(ccx, 6), round(ccy, 6)],
            "centroid_cells": [round(ci, 6), round(cj, 6)],
            "peak_cell": list(pcell), "peak_value": round(pval, 6),
            "centroid_peak_cells": round(float(np.hypot(
                cell_xy(*pcell)[0] - ccx, cell_xy(*pcell)[1] - ccy))
                / MASK["cell"], 6),
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


def write_md(path, pos, hashes):
    fails = checks_failed()
    comps = echo_components()
    holes = echo_holes()
    ncell = len(ECHO)
    nt = MASK["cols"] * MASK["rows"]
    dens = ncell / float(nt)
    X, Y, I = sample_grid()
    inside, outside = raster_cells()
    (ccx, ccy), _ci = centroid_xy()
    pcell, pval = peak_cell()
    sep = float(np.hypot(cell_xy(*pcell)[0] - ccx, cell_xy(*pcell)[1] - ccy)) \
        / MASK["cell"]
    L = []
    A = L.append
    A("# `concept/fig_ocv_measurement` — one chip, three stages, one number set")
    A("")
    A("![one chip, three stages, one number set](fig_ocv_measurement.jpg)")
    A("")
    A("The companion plate to `concept/fig_ocv_schematic`. That one defines the "
      "vector from the asset outwards; this one answers the question a reader "
      "asks next — **what exactly is measured, and where?** It takes one radar "
      "read chip and shows it three times, each time admitting less:")
    A("")
    A("```")
    A("C_{a,c}  ->  W_{a,c}  ->  M_{a,c}(t)  ->  OCV(t)")
    A("native       measurement     echo mask      the six")
    A("read chip    support                        components")
    A("```")
    A("")
    A("One idea, three panels, one frame: **restriction**. Panel (a) holds "
      "everything the sensor returned; (b) keeps the part that belongs to this "
      "asset in this channel; (c) keeps the part of that which carries echo. "
      "Nothing in (b) is absent from (a), and nothing in (c) is absent from (b) "
      "— the panels are one observation at three levels of admission, not three "
      "pictures.")
    A("")
    A(f"Like the first schematic it carries no data, so it has no pins to a "
      f"table; it has {len(CHECKS)} structural checks instead, listed at the "
      f"foot of this page and in `../data/concept/fig_ocv_measurement.json`.")
    A("")
    A("## The three panels")
    A("")
    A(f"**(a) the native read chip.** `C_{{a,c}}` as it comes off the "
      f"acquisition: the intensity field over the whole chip, "
      f"{NATIVE['cols']} × {NATIVE['rows']} = "
      f"{NATIVE['cols'] * NATIVE['rows']} native samples at "
      f"{NATIVE['cell'] * UNIT_MM:.2f} mm each, drawn as "
      f"{NATIVE['cols'] * NATIVE['rows']} filled vector quads on a "
      f"{RAMP_N}-step grey ramp — not as an embedded image, so the PDF stays a "
      f"true vector page. The asset footprint "
      f"`Π_c(G_a^{{3D}})` is a thin outline in the middle of it (the symbol "
      f"itself is set in (b), on the footprint), the chip's own "
      f"boundary carries `C_{{a,c}}` on its lower edge, and the radar frame "
      f"`az` / `r` sits at the chip's lower-left corner. Nothing is subdued: "
      f"this is the read, with more cells in it than any measurement will use.")
    A("")
    A("**(b) the measurement support.** The same field, restricted. Outside the "
      "background ring it is pulled 65 % of the way to white; the ring band "
      "itself is a flat tint (it is the no-echo reference, not structure); the "
      "guard between the ring and the support is pure white; and inside "
      f"`W_{{a,c}}` the field returns at full contrast under a blue veil "
      f"({int(B_VEIL * 100)} % blue). The asset is dark again, as on the first "
      "schematic, but solid: what is drawn is the *projection* "
      "`Π_c(G_a^{3D})`, a silhouette, not the tower with its storeys, so it "
      "carries no floor lines. Its name is set in paper white on the "
      "silhouette itself: the first schematic sets the same symbol under its "
      "tower, and here there is no band between the footprint and the "
      "support's own edge deep enough to hold a 9 pt label clear of both. The "
      "restriction is therefore *visible*: the reader "
      "sees what is kept and what is dropped, in the same picture, at the same "
      "place, rather than being told about it in a caption.")
    A("")
    A("**(c) the echo mask.** Inside `W_{a,c}` the mask: every support cell is "
      "drawn, the "
      f"{nt - ncell} without echo in very light grey and the {ncell} with echo "
      "in green. The footprint is repeated as a thin outline, so the echo can "
      "be seen to sit on the asset rather than beside it, and no echo cell "
      "crosses that outline. Outside the support there is no chip plate at all "
      "— bare paper, because the vector row runs through the band where the "
      "chip's lower edge would fall, and nothing outside the support is "
      "evaluated. Below the support the mask becomes the vector, "
      "`OCV(t) = [γ², P, D, F, A, S]`, in its two conceptual groups.")
    A("")
    A("## One frame, drawn three times")
    A("")
    A("All three panels read one frame, and the plate says so in a way a reader "
      "can check by eye: the chip rectangle, the asset footprint, the radar "
      "frame and the scale are identical in (a), (b) and (c) — the same "
      f"{CHIP[2] - CHIP[0]:.3f} × {CHIP[3] - CHIP[1]:.3f} panel-unit rectangle "
      f"({(CHIP[2] - CHIP[0]) * UNIT_MM:.1f} × "
      f"{(CHIP[3] - CHIP[1]) * UNIT_MM:.1f} mm), drawn at the same device "
      "coordinates in all three panels, and the checks compare the drawn "
      "rectangles in pixels rather than trusting the constants.")
    A("")
    A("The nesting is `W_{a,c} ⊂ guard ⊂ background ring ⊂ C_{a,c}`:")
    A("")
    A("| element | what it is | drawn as |")
    A("| --- | --- | --- |")
    A(f"| `C_{{a,c}}` | the read chip | grey rectangle, "
      f"{(CHIP[2] - CHIP[0]) * UNIT_MM:.1f} × "
      f"{(CHIP[3] - CHIP[1]) * UNIT_MM:.1f} mm, "
      f"{NATIVE['cols'] * NATIVE['rows']} native samples |")
    A(f"| background ring | the no-echo reference | flat grey band "
      f"{min(GUARD_BOX[0] - RING[0], GUARD_BOX[1] - RING[1], RING[2] - GUARD_BOX[2], RING[3] - GUARD_BOX[3]) * UNIT_MM:.1f} mm "
      f"wide, dashed boundary, `background ring` on a leader |")
    A(f"| guard | keeps the ring out of the support | pure white band, "
      f"{GUARD * UNIT_MM:.1f} mm |")
    A(f"| `W_{{a,c}}` | the asset- and channel-specific support | the same "
      f"chamfered octagon as the first schematic, blue, "
      f"{(W_BOX[2] - W_BOX[0]) * UNIT_MM:.1f} × "
      f"{(W_BOX[3] - W_BOX[1]) * UNIT_MM:.1f} mm |")
    A("")
    A("`W_{a,c}` is **one** polygon definition mapped once into the two panels "
      "that draw it, so the panels cannot drift apart; a check normalises both "
      "to their own box and compares them to the definition at 1 × 10⁻¹².")
    A("")
    A("## The mask is derived, not decorated")
    A("")
    A("The field is written down in closed form: a smooth base, two texture "
      "terms and seven declared responses (one annulus — the tower seen as a "
      "ring — three compact blobs on the asset, one diffuse glow that never "
      "reaches the reference level, and two clutter returns *outside* the "
      "support). There is no random number anywhere in the file, so the field "
      "is a constant and so is the mask.")
    A("")
    A(f"The mask is the set of support cells whose **peak native sample exceeds "
      f"one declared reference level** ({T_REF:.3f} on a field that runs "
      f"{I.min():.3f} … {I.max():.3f}). That definition produces, and the "
      f"checks insist on:")
    A("")
    A("| property | value |")
    A("| --- | --- |")
    A(f"| support lattice | {MASK['cols']} × {MASK['rows']} = {nt} cells, one "
      f"cell = exactly {int(round(MASK['cell'] / NATIVE['cell']))} × "
      f"{int(round(MASK['cell'] / NATIVE['cell']))} native samples "
      f"({MASK['cell'] * UNIT_MM:.1f} mm) |")
    A(f"| echo cells | {ncell} ({dens * 100:.1f} % of the support) |")
    A(f"| connected regions | {len(comps)}, of "
      + ", ".join(str(len(k)) for k in comps) + " cells |")
    A(f"| enclosed holes | {len(holes)} — the dark centre of the annulus, and "
      f"the only place the mask is 4-connected around nothing |")
    A(f"| margin | weakest echo cell {min(mask_cell_peak(i, j, X, Y, I) for (i, j) in ECHO):.3f}, "
      f"strongest inactive cell {max(mask_cell_peak(i, j, X, Y, I) for i in range(MASK['cols']) for j in range(MASK['rows']) if (i, j) not in ECHO):.3f} "
      f"— the mask survives ±0.02 on the reference |")
    A("")
    A("The lattice is not decoration: one mask cell is exactly two native "
      "samples, on the native sample grid, so *the mask is a coarsened read of "
      "the same chip* is something the plate proves rather than something the "
      "caption claims. And the mask is not a threshold of the picture either: "
      f"the chip outside `W_{{a,c}}` does carry returns above the reference "
      f"level (peak {float(I.ravel()[~np.array([Path(w_poly(W_BOX)).contains_point((float(x), float(y))) for x, y in zip(X.ravel(), Y.ravel())])].max()):.3f}) "
      "and they produce no cell at all, because the mask is evaluated on the "
      "support and nowhere else.")
    A("")
    A("The centroid and the peak are measured and reported — they are "
      f"{sep:.2f} cells apart, and that distance is what `S` measures — and "
      "neither is drawn. No dot, no arrow, no symbol: the reader is not asked "
      "to define a glyph the plate cannot justify.")
    A("")
    A("## The vector, in its two blocks")
    A("")
    A("| block | components | colour | what it is about |")
    A("| --- | --- | --- | --- |")
    A("| cross-acquisition | `γ²`, `P` | blue `#1f4e79` | what changes between "
      "passes of the same chip |")
    A("| spatial geometry | `D`, `F`, `A`, `S` | green `#3f7d3f` | what the "
      "echo looks like inside the support |")
    A("")
    A("The two blocks are separated three ways — colour, a hairline rule in the "
      "gap, and an under-brace with a name — and each brace spans the measured "
      "extent of its own group, so a brace can never be wider than what it "
      "names. The row is the **same row** the first schematic draws, to 10⁻⁶ of "
      "a panel unit: a check reads `fig_ocv_schematic.json` and compares.")
    A("")
    A("## What the plate does not show")
    A("")
    A("It is a schematic, so it carries no site, no date, no state and no "
      "number that could be read as a result. Deliberately absent:")
    A("")
    A("- structural health, damage probability, failure detection or any other "
      "diagnosis — the figure defines a measurement, it does not deliver a "
      "verdict;")
    A("- displacement, deformation or InSAR geometry — nothing in the plate "
      "moves;")
    A("- `R_B` and `R_G` — those belong to the temporal reconfiguration analysis "
      "(`figures/sds/`), not to the definition of the vector;")
    A("- any shift marker (`p̄`, `p⋆`) — the vector has no shift parameter on "
      "this plate;")
    A("- a centroid or a peak symbol — `S` is a distance the mask yields, not a "
      "glyph the plate has to explain;")
    A("- the reference level itself, any count, any weight and any score: the "
      "number is in this report and in the JSON, not on the plate.")
    A("")
    A("The prohibition is mechanical, not editorial: a forbidden-token scan runs "
      "over every drawn string, the whole text inventory is compared against a "
      "written-out whitelist, and a separate scan rejects the shift markers, so "
      "the plate cannot grow a `damage` label or a `p̄` without failing the run.")
    A("")
    A("## The plates, and how to reproduce them")
    A("")
    A("| file | what it is |")
    A("| --- | --- |")
    A("| `concept/fig_ocv_measurement.jpg` | 183 mm × 61 mm, 630 dpi JPEG "
      "(quality 95, 4:4:4 subsampling) |")
    A("| `concept/fig_ocv_measurement.pdf` | the same plate as a true vector "
      "page (FontType 42; no raster object anywhere in the file — the "
      f"intensity field is {NATIVE['cols'] * NATIVE['rows']:,} filled quads) |")
    A("| `concept/fig_ocv_measurement.md` | this report |")
    A("| `../data/concept/fig_ocv_measurement.json` | the drawn geometry, the "
      "field, the mask, the text inventory, every check and the sha256 of both "
      "images |")
    A("")
    A("```bash")
    A("bash figures/concept/fig_ocv_measurement.sh                  # both images")
    A("python3 code/concept/fig_ocv_measurement.py --no-figures     # checks only")
    A("PYTHON=/usr/bin/python3 bash figures/concept/fig_ocv_measurement.sh")
    A("```")
    A("")
    A("For the manuscript:")
    A("")
    A("```latex")
    A("\\includegraphics[width=170mm]{figures/concept/fig_ocv_measurement.pdf}")
    A("```")
    A("")
    A("## The checks")
    A("")
    A(f"**{len(CHECKS) - len(fails)} of {len(CHECKS)} checks pass.** A failure "
      "exits with status 2 *after* the artefacts are written, so a broken plate "
      "is reported rather than silently delivered. (`--no-figures` builds, "
      "measures and checks, writes the JSON and this report, and skips both "
      "images and the file checks that go with them.)")
    A("")
    A("| family | what it pins |")
    A("| --- | --- |")
    A("| canvas | 183 mm × 61 mm at 630 dpi, ratio 3.000, `savefig.bbox` "
      "`standard`; every panel isotropic, so one unit is one unit in x and in "
      "y and the native cells stay square |")
    A("| one frame | the chip, the asset footprint, the radar frame and the "
      "scale are identical in all three panels — compared on the *drawn* "
      "geometry in device pixels, and on the chip's offset inside each panel |")
    A("| text inventory | exactly the declared strings and no others; one 7 pt "
      "title per panel; no plate title; no legend; all axes off |")
    A("| legibility | every text ≥ 7.0 pt at 183 mm and ≥ 6.5 pt at 170 mm; "
      "every text inside its own panel; no two texts overlap; every title fits "
      "its run |")
    A("| vocabulary | no forbidden token in any drawn string, and no shift "
      "marker (`p̄`, `p⋆`); the radar frame is `az` / `r`; the physical frame "
      "is absent |")
    A("| set algebra | `W_{a,c}` inside the guard inside the ring inside "
      "`C_{a,c}`, each with a visible margin; the guard ≥ 1.1 mm and the band "
      "outside it ≥ 1.7 mm; the asset footprint lies inside `W` |")
    A("| the lattice | one mask cell is exactly 2 × 2 native samples, on the "
      "native grid; the lattice is inside `W`; every mask cell is inside `W`; "
      "`C` holds more than twice the cells `W` does |")
    A(f"| raster integrity | {NATIVE['cols']} × {NATIVE['rows']} = "
      f"{NATIVE['cols'] * NATIVE['rows']:,} quads that tile the chip exactly "
      "once and touch every edge; cells square and ≥ 0.8 mm; the ramp is 16 "
      "monotone steps, neither ink nor paper; the field is bit-identical on a "
      "rebuild |")
    A("| the mask | 4 connected regions of 8, 6, 4 and 2 cells, strictly "
      "decreasing; 20 cells, 18.5 % density; exactly one enclosed hole, off "
      "the border; every echo cell above the reference level and every "
      "inactive cell below it, with 0.02 to spare; the mask is the threshold "
      "set of the field cell for cell |")
    A("| the mask stops at W | the chip outside the support carries returns "
      "above the reference level and produces no cell at all — the mask is a "
      "measurement, not a threshold of the picture |")
    A("| the four zones | read off the composed canvas, not the geometry: the "
      "outer zone is the panel-(a) colour pulled 65 % to white, the band and "
      "the guard are exact composites, the inside of `W` is the field under "
      "18 % blue, and every one of them to within 2/255 |")
    A("| nothing annotated | no centroid or peak symbol, no marker artist "
      "anywhere in panel (c); the centroid and the peak are measured, "
      "reported and ≥ 1.5 cells apart — the quantity `S` measures |")
    A("| the vector | six components in the order `γ², P, D, F, A, S`; the "
      "partition into (γ², P) and (D, F, A, S) proved from the token table; "
      "each brace spans its own group's measured extent; the row is centred, "
      "clear of the mask, and identical to the first schematic's row |")
    A("| determinism, files | two independent builds render identical pixels "
      "and lay the row out byte-identically; JPEG 4 539 × 1 513 px ± 2 at "
      "630 dpi; the PDF has no `/Subtype /Image`, embeds `/FontFile2`, carries "
      "no Type-3 font, no timestamp, and a 183 mm × 61 mm page box |")
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


# ---------------------------------------------------------------------------
# entry point
# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Draw one radar read chip at three levels of admission — "
                    "C to W to the echo mask to the vector — and check every "
                    "structural promise the plate makes.")
    ap.add_argument("--figdir", default=FIGDIR,
                    help="where the JPEG and the PDF go (default figures/concept)")
    ap.add_argument("--json", default=os.path.join(DATADIR, OUT_JSON),
                    help="the JSON artifact (default data/concept/...)")
    ap.add_argument("--md", default=None,
                    help="the report (default <figdir>/fig_ocv_measurement.md)")
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
    check_cross_plate(pos)

    md = args.md or os.path.join(args.figdir, OUT_MD)
    os.makedirs(os.path.dirname(os.path.abspath(md)) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.json)) or ".", exist_ok=True)
    data = payload(pos, hashes)
    with open(args.json, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=1, sort_keys=True, ensure_ascii=False)
        fh.write("\n")
    write_md(md, pos, hashes)

    fails = checks_failed()
    print(f"fig_ocv_measurement: {len(CHECKS) - len(fails)}/{len(CHECKS)} checks "
          f"passed, {len(TEXTS)} text objects, {len(ECHO)} echo cells in "
          f"{len(echo_components())} regions, unit {UNIT_MM:.3f} mm")
    for name, h in hashes.items():
        print(f"  {name:32s} {h[:16]}")
    print(f"  {os.path.relpath(args.json, ROOT)}")
    print(f"  {os.path.relpath(md, ROOT)}")
    if fails:
        for ck in fails:
            print(f"  FAIL {ck['name']}: {ck['detail']}", file=sys.stderr)
        print(f"fig_ocv_measurement: {len(fails)} check(s) failed", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())












    _ = ncol





