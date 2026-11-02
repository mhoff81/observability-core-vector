"""The one typography and the one print width of this repository's figures.

Every figure under ``figures/`` is printed at the same width, **170 mm** — the
text width it is set in — and every figure therefore has to carry the same
point sizes *at that width*. Before this module existed each generator carried
its own ``COL_WIDTH = 183 * MM`` and its own ``rcParams`` block, so the same
class of panel came out between 3.2 pt and 13.4 pt depending on how far
``bbox_inches="tight"`` let a long title or a wide legend push the canvas past
183 mm: the type is authored in points, so a canvas that grows *after* the
layout shrinks every glyph with it. The rule this module encodes is the one the
print check needs:

  * the **canvas** is ``FINAL_WIDTH_MM * MM`` inches wide and is not exceeded —
    a figure whose tight bounding box is wider than the canvas is a defect in
    the layout, not something to absorb by scaling;
  * the **type** is a fixed, small ladder: panel title 8.8 pt, axis label
    7.6 pt, tick and legend 6.8 pt, a ``suptitle`` at 9.0 pt, annotations at
    6.2 pt, all hung off ``font.size`` 8.0;
  * the **pen** is a small ladder too (:data:`LW`): a drawn series 1.4 pt, a
    plane path or a rule 1.0, the grid 0.6, a marker edge 0.9.  It is *not* in
    :func:`rcparams` -- a generator that wants it passes it as an override, so
    the plates already written keep their pens;
  * a **frame shared by several plates** is a constant, not a copy
    (:data:`TRAJECTORY`): the six timeline plates of this repository are the
    same 170 x 156 mm template and are read from one table, so they cannot
    drift apart the way six hard-coded ``figsize`` calls do;
  * a drawn figure is **measured after it is written** (``measure``) and a
    deviation from 170 ± 0.5 mm aborts the run, exactly like a failed pin.

Nothing here imports matplotlib at module level: the constants and the image
reader are stdlib-only, so ``audit_figure_style.py`` can check the sources and
the artefacts on an interpreter without matplotlib.

Usage (in a generator):

    import figure_style as style
    style.apply()                          # the rcParams block, once
    fig = plt.figure(figsize=(style.COL_WIDTH, 0.9 * style.COL_WIDTH))
    ...
    leg = fig.legend(...)                  # a figure artist, under the grid
    fig.tight_layout(rect=style.grid_rect(fig, leg))
    style.assert_fits(fig)                 # canvas is 170 mm, tight box fits
    style.save(fig, path)                  # write at that canvas, then measure

``save`` writes the raster and then, beside it and with the same stem, a ``.pdf``
and an ``.eps`` -- and measures *all three*.  The PDF is the vector copy for
LaTeX and keeps the transparency; PostScript has no alpha, so the EPS resolves
it per ``eps=`` (``"blend"`` composites each translucent artist over white,
``"raster"`` takes the PDF through ``pdftops`` -- see :func:`blend_alpha` and
:func:`eps_from_pdf`).

A generator reaches this module the way the site packages reach each other: the
site directory is already on the path (``sys.path.insert(0, HERE)``), so
``code/`` is one line further up --

    HERE = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(1, os.path.dirname(HERE))      # code/figure_style.py

-- or, from a script with no preamble yet, the ``REPO`` form of
``sds/sds_support.py``::

    HERE = os.path.dirname(os.path.abspath(__file__))
    REPO = os.path.dirname(os.path.dirname(HERE))
    for path in (os.path.join(REPO, "code"),):
        if path not in sys.path:
            sys.path.insert(0, path)

Usage (as a check):

    python3 code/figure_style.py            # print the ladder, measure files

The CLI measures whatever it is given: a raster (PNG/JPEG, width from its
declared density) or a vector container (PDF ``/MediaBox``, EPS
``%%HiResBoundingBox``), which carries its width in points and has no density at
all.
"""
from __future__ import annotations

import os
import re
import shutil
import struct
import subprocess
import sys
import tempfile

MM = 1.0 / 25.4

#: The width every figure is finally set at, in millimetres.
FINAL_WIDTH_MM = 170.0

#: That width in inches, to be handed to ``figsize``.
COL_WIDTH = FINAL_WIDTH_MM * MM

#: A saved figure may deviate this much from :data:`FINAL_WIDTH_MM` (mm).
WIDTH_TOL_MM = 0.5

#: The point-size ladder.  One value per role, used by *every* figure.
FS = {
    "base": 8.0,        # font.size -- the anchor everything else hangs off
    "title": 8.8,       # axes.titlesize, the panel title
    "suptitle": 9.0,    # figure.titlesize, the one-line title above the grid
    "label": 7.6,       # axes.labelsize, the axis label
    "tick": 6.8,        # xtick/ytick.labelsize
    "legend": 6.8,      # legend.fontsize
    "annotation": 6.2,  # text/annotate inside a panel -- a floor, not a value
}

#: The band each role has to fall in, in points.  ``None`` is an open end.
BANDS = {
    "title": (8.5, 9.0),
    "suptitle": (8.5, 9.0),
    "label": (7.0, 8.0),
    "tick": (6.5, 7.0),
    "legend": (6.5, 7.0),
    "annotation": (6.0, None),
}

#: rcParams key -> role.  The keys :func:`rcparams` writes and the audit reads.
RC_ROLES = {
    "font.size": "base",
    "axes.titlesize": "title",
    "figure.titlesize": "suptitle",
    "axes.labelsize": "label",
    "xtick.labelsize": "tick",
    "ytick.labelsize": "tick",
    "legend.fontsize": "legend",
}

#: Roles the type check fails on.  ``base`` is a fixed anchor and annotations
#: are only held to a floor, so neither is a band.
CHECKED_ROLES = ("title", "suptitle", "label", "tick", "legend")

#: The pen ladder, in points.  Deliberately thin: five widths, each a *kind* of
#: ink rather than a per-plate taste.  A generator that wants it hands it to
#: :func:`apply` as an override, so no plate already written moves::
#:
#:     plt = style.apply(dpi=FIG_DPI, **{"axes.linewidth": style.LW["spine"],
#:                                       "grid.linewidth": style.LW["grid"]})
#:
#: ``series`` is a monthly series line, ``plane`` a trajectory path inside a
#: plane panel, ``rule`` the event and zero lines, ``grid`` the mesh, ``spine``
#: the axes frame, ``marker_edge`` the outline of a drawn month, ``bar_edge``
#: the outline of a hollow bar, and ``band`` a filled span or envelope, which
#: has no outline at all.
LW = {
    "series": 1.4,
    "plane": 1.0,
    "rule": 1.0,
    "grid": 0.6,
    "spine": 1.0,
    "marker_edge": 0.9,
    "bar_edge": 0.9,
    "band": 0.0,
}

#: The marker size of a drawn month, in points, per class of panel.
MS = {"month": 3.4, "plane": 3.6, "legend": 4.2}

#: The frame of the six OCV **timeline** plates -- one table, two packages.
#:
#: ``cases/ocv_trajectory.py`` (the five case sites) and
#: ``werder_havel/fig_werder_havel_ocv_trajectory.py`` (the interval site) draw
#: the same 3 x 4 template, so they must also deliver the same physical plate:
#: **170 x 156 mm**, raster plus ``.pdf`` and ``.eps``, nothing tight-cropped.
#: Every number a reader could compare between two of those plates lives here
#: and nowhere else; a per-generator copy of ``figsize`` is exactly how the six
#: came out between 190.2 and 193.3 mm wide before.
#:
#: The margins are fractions of the canvas.  ``bottom`` holds the rotated date
#: ticks and the plate's one legend; ``top`` holds the panel tag of the first
#: row; ``wspace``/``hspace`` are matplotlib's, in units of the axes size, and
#: are wide enough for the y label of the next column and the date ticks of the
#: row below.
TRAJECTORY = {
    "width_mm": FINAL_WIDTH_MM,
    "height_mm": 156.0,
    "grid": {
        "left": 0.080,
        "right": 0.938,
        "top": 0.966,
        "bottom": 0.170,
        "wspace": 0.55,
        "hspace": 0.55,
    },
    #: The band the legend is allowed to occupy, in mm, measured up from the
    #: canvas' bottom edge.  The legend is re-flowed until it fits this width
    #: *and* this band, so it can neither widen the plate nor grow into the
    #: date ticks of the bottom row.
    "legend": {"rows": 2, "band_mm": 15.0, "gap_mm": 2.0, "margin_mm": 4.0},
}

#: Every figure is written at this resolution (600 dpi = the print pipeline).
FIGURE_DPI = 600

#: The vector containers written beside each raster, under the same stem.
VECTOR_FORMATS = ("pdf", "eps")

#: How an EPS resolves the alpha PostScript cannot express.
#:
#: ``"blend"``  -- composite each translucent artist over white (vector text,
#:                 colours as the PDF shows them on white, overlaps do not
#:                 accumulate);
#: ``"raster"`` -- take the PDF through poppler's ``pdftops``: the plate stays
#:                 pixel-identical, but the page becomes a 300 dpi bitmap;
#: ``"opaque"`` -- leave matplotlib's own flattening (tints at full strength).
EPS_MODES = ("blend", "raster", "opaque")

#: Points per millimetre, the other way round: vector boxes are given in points.
PT_TO_MM = 25.4 / 72.0


def rcparams(dpi=FIGURE_DPI, **overrides):
    """The complete ``rcParams`` block of this repository's figures.

    Only typography, output resolution and font embedding are set: nothing here
    touches line widths, markers or colours, so a generator that switches to
    this block keeps its drawn geometry and changes only its type.
    """
    rc = {
        "font.family": "DejaVu Sans",
        "font.size": FS["base"],
        "axes.titlesize": FS["title"],
        "axes.labelsize": FS["label"],
        "xtick.labelsize": FS["tick"],
        "ytick.labelsize": FS["tick"],
        "legend.fontsize": FS["legend"],
        "figure.titlesize": FS["suptitle"],
        "savefig.dpi": dpi,
        "figure.dpi": 100,
        # The canvas is saved whole: ``None`` is how matplotlib >= 3.11 spells
        # "the figure, not its bounding box" (older versions call it
        # "standard").  ``"tight"`` here would crop the canvas to the drawn
        # content and shrink every glyph with it -- the defect this module
        # exists to stop -- so the block names the default explicitly.
        "savefig.bbox": None,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
    rc.update(overrides)
    return rc


def apply(dpi=FIGURE_DPI, **overrides):
    """Push :func:`rcparams` into matplotlib and return ``pyplot``."""
    import matplotlib.pyplot as plt
    plt.rcParams.update(rcparams(dpi=dpi, **overrides))
    return plt



# ---------------------------------------------------------------------------
# Measuring what was written: PNG pHYs and JPEG JFIF/SOF, stdlib only
# ---------------------------------------------------------------------------
def _png_size(data):
    """(width, height, dpi) of a PNG, dpi from the ``pHYs`` chunk (or None)."""
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError("not a PNG")
    width, height = struct.unpack(">II", data[16:24])
    dpi = None
    i = 8
    while i + 8 <= len(data):
        length = struct.unpack(">I", data[i:i + 4])[0]
        kind = data[i + 4:i + 8]
        if kind == b"pHYs":
            x_ppm, _y_ppm, unit = struct.unpack(">IIB", data[i + 8:i + 17])
            if unit == 1 and x_ppm:
                dpi = x_ppm * 0.0254
        if kind in (b"IDAT", b"IEND"):
            break
        i += 12 + length
    return width, height, dpi


def _jpeg_size(data):
    """(width, height, dpi) of a JPEG, dpi from JFIF ``APP0`` (or None)."""
    i, dpi, size = 2, None, None
    while i + 3 < len(data):
        if data[i] != 0xFF:
            i += 1
            continue
        marker = data[i + 1]
        if marker in (0xD8, 0xD9) or 0xD0 <= marker <= 0xD7 or marker == 0xFF:
            i += 2 if marker != 0xFF else 1
            continue
        length = struct.unpack(">H", data[i + 2:i + 4])[0]
        if marker == 0xE0 and data[i + 4:i + 9] == b"JFIF\x00":
            unit = data[i + 11]
            x_density = struct.unpack(">H", data[i + 12:i + 14])[0]
            if unit == 1:
                dpi = float(x_density)          # dots per inch
            elif unit == 2:
                dpi = x_density * 2.54          # dots per centimetre
        elif 0xC0 <= marker <= 0xCF and marker not in (0xC4, 0xC8, 0xCC):
            height, width = struct.unpack(">HH", data[i + 5:i + 9])
            size = (width, height)
        if marker == 0xDA:                       # start of scan: metadata done
            break
        i += 2 + length
    if size is None:
        raise ValueError("no SOF marker")
    return size[0], size[1], dpi


def _pdf_size(data):
    """(width_pt, height_pt) of the first page, from ``/MediaBox`` or ``/CropBox``.

    A PDF states its page size in points, not in density and pixels: the box is
    what the page *is*, so the print width follows from it with no dpi involved.
    """
    for key in (b"/MediaBox", b"/CropBox"):
        found = re.findall(key + rb"\s*\[\s*(-?[\d.]+)\s+(-?[\d.]+)\s+"
                                  rb"(-?[\d.]+)\s+(-?[\d.]+)\s*\]", data)
        if found:
            x0, y0, x1, y1 = (float(v) for v in found[0])
            return x1 - x0, y1 - y0
    raise ValueError("no /MediaBox in the PDF")


def _eps_size(data):
    """(width_pt, height_pt) from ``%%HiResBoundingBox``, else ``%%BoundingBox``.

    The HiRes box is the exact one (a matplotlib EPS writes
    ``0 0 481.889764 481.889764`` for 170 mm); the plain box is the same area
    rounded to whole points, so it is only a fallback.  ``%%BoundingBox: (atend)``
    carries no numbers and is skipped by requiring four of them.
    """
    for key in (b"%%HiResBoundingBox:", b"%%BoundingBox:"):
        found = re.findall(key + rb"[ \t]*(-?[\d.]+)[ \t]+(-?[\d.]+)[ \t]+"
                                  rb"(-?[\d.]+)[ \t]+(-?[\d.]+)", data)
        if found:
            x0, y0, x1, y1 = (float(v) for v in found[0])
            if x1 > x0 and y1 > y0:
                return x1 - x0, y1 - y0
    raise ValueError("no bounding box in the EPS")


def measure(path):
    """What a written figure really is: width/height in mm, dpi, pixels.

    Returns a dict with ``path``, ``format``, ``vector``, ``px``, ``dpi``, ``mm``
    and ``height_mm``.  For a raster ``dpi`` is the declared density and ``px``
    the pixel size; ``dpi`` is ``None`` when the file carries no density and
    ``mm`` is then ``None`` -- a raster without a declared density has no print
    width, which is a defect in itself.  A vector container has no density by
    construction: its box *is* its size, so ``dpi`` and ``px`` are ``None`` and
    ``mm`` is the box in millimetres.
    """
    with open(path, "rb") as fh:
        data = fh.read()
    vector = False
    if data[:8] == b"\x89PNG\r\n\x1a\n":
        width, height, dpi = _png_size(data)
        fmt = "png"
    elif data[:2] == b"\xff\xd8":
        width, height, dpi = _jpeg_size(data)
        fmt = "jpeg"
    elif data[:5] == b"%PDF-":
        width, height = _pdf_size(data)
        dpi, fmt, vector = None, "pdf", True
    elif data[:4] == b"%!PS":
        width, height = _eps_size(data)
        dpi, fmt, vector = None, "eps", True
    else:
        raise ValueError("%s: neither PNG, JPEG, PDF nor EPS" % path)
    if vector:
        return {"path": path, "format": fmt, "vector": True, "px": None,
                "dpi": None, "mm": width * PT_TO_MM, "height_mm": height * PT_TO_MM}
    return {
        "path": path,
        "format": fmt,
        "vector": False,
        "px": (width, height),
        "dpi": dpi,
        "mm": (width / dpi * 25.4) if dpi else None,
        "height_mm": (height / dpi * 25.4) if dpi else None,
    }


def effective_pt(pt, width_mm, at_mm=FINAL_WIDTH_MM):
    """What *pt* becomes when a figure *width_mm* wide is set at *at_mm*."""
    if not width_mm:
        return None
    return pt * (at_mm / width_mm)


def width_mm(path):
    """The saved width of *path* in millimetres (``None`` without a density)."""
    return measure(path)["mm"]


def assert_width(path, expected_mm=FINAL_WIDTH_MM, tol=WIDTH_TOL_MM):
    """Abort unless *path* is *expected_mm* ± *tol* wide.  Returns the measure.

    A miss means the layout let something (a long panel title, a legend, a
    colour bar) grow the canvas past the print width -- the type is then set
    smaller than declared, which is exactly what this module exists to stop.
    """
    info = measure(path)
    got = info["mm"]
    if got is None:
        raise RuntimeError("%s: no declared density, so the print width is "
                           "unknown -- save with an explicit dpi" % path)
    if abs(got - expected_mm) > tol:
        raise RuntimeError(
            "%s: %.1f mm wide, expected %.1f ± %.1f mm -- the layout overflows "
            "the print width, so every point size shrinks by x%.3f"
            % (path, got, expected_mm, tol, expected_mm / got))
    return info


def tight_size_mm(fig):
    """(canvas_mm, tight_mm) of a drawn figure, before it is written.

    The fit guard: if ``tight_mm`` exceeds ``canvas_mm`` the save would widen
    the canvas and shrink the type, so the layout has to be fixed first.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    box = fig.get_tightbbox(renderer)
    return fig.get_size_inches()[0] * 25.4, box.width * 25.4


def assert_fits(fig, tol=WIDTH_TOL_MM):
    """Abort unless the drawn figure is the print width and fits its canvas.

    Two failures, both layout defects rather than write-time surprises: a canvas
    that is not ``FINAL_WIDTH_MM`` (the generator kept its own ``COL_WIDTH``),
    and a tight bounding box wider than the canvas (the text -- a long title, a
    wide legend, a colour bar -- does not fit, so printing would have to shrink
    every glyph).  Returns ``(canvas_mm, tight_mm)``.
    """
    canvas, tight = tight_size_mm(fig)
    if abs(canvas - FINAL_WIDTH_MM) > tol:
        raise RuntimeError(
            "the canvas is %.1f mm wide, not the print width %.1f ± %.1f mm -- "
            "use style.COL_WIDTH as the figure width"
            % (canvas, FINAL_WIDTH_MM, tol))
    if tight > canvas + tol:
        raise RuntimeError(
            "the drawn figure overflows its canvas: tight %.1f mm > canvas "
            "%.1f mm -- shorten the text or re-cut the grid, do not widen the "
            "canvas" % (tight, canvas))
    return canvas, tight


# ---------------------------------------------------------------------------
# The EPS sidecar: PostScript has no alpha, so it has to be resolved
# ---------------------------------------------------------------------------
def _declared_alpha(artist):
    """The alpha *artist* carries, or ``None`` when it carries none.

    A per-mark alpha array (a ``Collection``) is not one number and is left to
    the colours: matplotlib folds it into them, so the colour walk sees it
    there.  ``float()`` is the test: it accepts neither ``None`` nor an array.
    """
    getter = getattr(artist, "get_alpha", None)
    try:
        return float(getter()) if callable(getter) else None
    except (TypeError, ValueError):
        return None


def _over_background(rgba, alpha, background):
    """``alpha*C + (1-alpha)*background`` per channel, opaque."""
    return tuple(alpha * c + (1.0 - alpha) * b
                 for c, b in zip(rgba[:3], background)) + (1.0,)


def _marks_of(artist, alpha, to_rgba):
    """[(attribute, index, colour, effective_alpha)] of *artist*.

    The effective alpha is the artist's own times the colour's; *index* is the
    position in a per-mark colour array (a ``Collection``) and ``None`` for a
    single colour.
    """
    marks = []
    for attribute in ("facecolor", "edgecolor", "color"):
        getter = getattr(artist, "get_" + attribute, None)
        if not callable(getter):
            continue
        try:
            value = getter()
        except Exception:                       # an artist without that property
            continue
        rows = getattr(value, "ndim", 1) == 2
        for index, color in enumerate(list(value) if rows else [value]):
            if isinstance(color, str) and color in ("none", "auto"):
                continue
            try:
                rgba = to_rgba(color)
            except (ValueError, TypeError):
                continue
            effective = rgba[3] * alpha
            if 0.0 < effective < 1.0:
                marks.append((attribute, index if rows else None, rgba,
                              effective))
    return marks


def _rgba_marks(artist, to_rgba):
    """``(alpha, marks)`` of *artist*, read with its alpha suspended.

    A ``Patch`` and a ``Collection`` fold their own alpha into the colour they
    return (``set_alpha`` re-derives the face colour from
    ``_original_facecolor``), so reading them as they stand would count the
    alpha twice -- 0.35 *and* 0.35 where the file carries 0.35 once.  The alpha
    is suspended for the read and put back afterwards, so the artist is left
    exactly as it was found.

    An artist with no alpha of its own is still inspected: ``grid.alpha`` and a
    colour written as ``(r, g, b, a)`` put the transparency in the *colour*, and
    the PS backend flattens those just the same.
    """
    declared = _declared_alpha(artist)
    alpha = 1.0 if declared is None else declared
    if alpha < 1.0:
        artist.set_alpha(1.0)
    try:
        marks = _marks_of(artist, alpha, to_rgba)
    finally:
        if alpha < 1.0:
            artist.set_alpha(declared)
    return alpha, marks


def alpha_artists(fig):
    """Every artist of *fig* whose marks are drawn below full opacity.

    Two artists are reachable only by name: ``Text.get_children`` and
    ``Annotation.get_children`` return ``[]``, so a label's box and an
    annotation's arrow are in no artist tree and ``fig.findobj`` never sees
    them, although the PS backend draws them like any other patch.  Missing the
    arrow is not cosmetic -- its alpha is in the arrow, so the EPS would carry
    the raw colour where the PDF shows a tint.

    This is what :func:`blend_alpha` rewrites and what :func:`save` reports, so
    a plate with 68 translucent artists says so instead of silently tinting.
    Read-only: every artist keeps the alpha it was found with.
    """
    from matplotlib.colors import to_rgba
    from matplotlib.text import Text

    found = []
    for art in fig.findobj():
        alpha, marks = _rgba_marks(art, to_rgba)
        if marks:
            found.append((art, alpha, marks))
    for text in fig.findobj(Text):
        for patch in (text.get_bbox_patch(),
                      getattr(text, "arrow_patch", None)):
            if patch is None:
                continue
            alpha, marks = _rgba_marks(patch, to_rgba)
            if marks:
                found.append((patch, alpha, marks))
    return found


def blend_alpha(fig, background=(1.0, 1.0, 1.0)):
    """Resolve every translucent artist to the colour it has over *background*.

    matplotlib's PS backend does **not** blend a translucent artist into an EPS,
    it writes the *unmixed* colour opaque: a marker at alpha 0.55 comes out at
    full strength while the PDF next to it shows the tint.  This walk does the
    mixing instead -- ``alpha*C + (1-alpha)*background`` per channel, then alpha
    1 -- so the EPS carries the colour the PDF shows on the white page.  The
    arithmetic is checkable in the file: matplotlib writes colours as plain
    ``setrgbcolor`` text, so a composite can be asserted by value.

    What it cannot do is accumulate: two translucent marks drawn on top of each
    other become two opaque marks at their single-layer composite, so a plate
    whose alpha *is* the information (an overplot density cloud) has to be
    written with ``eps="raster"`` instead of ``eps="blend"``.  For the same
    reason it does not look *through* a mark: a translucent white label box over
    a line becomes the white box, where the PDF still shows the line at 22 %.

    Returns ``[(type_name, [(attribute, before, after)])]`` -- one entry per
    artist, so the caller can report what was flattened.
    """
    blended = []
    for art, _alpha, marks in alpha_artists(fig):
        art.set_alpha(1.0)              # the setters below must not re-fold it
        changes = []
        for attribute, index, rgba, effective in marks:
            after = _over_background(rgba, effective, background)
            setter = getattr(art, "set_" + attribute)
            if index is None:
                setter(after)
            else:
                colors = list(getattr(art, "get_" + attribute)())
                if len(colors) == 1:            # one colour for every mark
                    setter([after])
                else:
                    colors[index] = after
                    setter(colors)
            changes.append((attribute, rgba, after))
        blended.append((type(art).__name__, changes))
    return blended


def _normalise_eps(path):
    """Drop the two EPS header lines that change from run to run.

    ``%%CreationDate`` is a timestamp and ``%%Title`` is the *output file name*
    (poppler writes the source PDF's path there), so two runs of the same plate
    never produce the same bytes.  Both are optional in DSC, and the PDF side
    sets ``metadata={"CreationDate": None, "ModDate": None}`` for the same
    reason, so the two containers are treated alike.  Returns the number of
    lines dropped.
    """
    with open(path, "rb") as fh:
        lines = fh.read().decode("latin-1").splitlines(True)
    kept = [line for line in lines
            if not line.startswith(("%%CreationDate:", "%%Title:"))]
    with open(path, "wb") as fh:
        fh.write("".join(kept).encode("latin-1"))
    return len(lines) - len(kept)


def eps_from_pdf(pdf_path, eps_path, timeout=180):
    """Turn our own PDF into an EPS with poppler's ``pdftops``.

    The faithful route, and the only one that keeps the plate pixel-identical to
    the PDF: PostScript cannot express alpha, so ``pdftops`` resolves it by
    rasterising what it cannot express -- the *whole page*, at 300 dpi, inside
    an EPS wrapper of exactly the right size.  The type is then a bitmap too,
    which is why this is not the default.

    Returns ``True`` when the file was written, ``False`` when ``pdftops`` is
    not installed -- the caller then falls back to the vector route rather than
    failing a run over a missing external tool.
    """
    exe = shutil.which("pdftops")
    if exe is None:
        return False
    subprocess.run([exe, "-eps", "-level3", pdf_path, eps_path],
                   check=True, timeout=timeout,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return os.path.exists(eps_path)


def save_eps(fig, path, eps="blend", source=None, expect_mm=FINAL_WIDTH_MM,
             tol=WIDTH_TOL_MM, timeout=180, verify=True):
    """Write *fig* as EPS at *path*, normalise it, and hold it to the width.

    The EPS sidecar without the rest of :func:`save`.  A generator that writes
    its own raster and PDF -- because its canvas is *not* the print width, as
    with the 183 mm concept plates -- borrows exactly this container instead of
    copying the alpha handling.  *expect_mm* is the width that plate is authored
    at (:data:`FINAL_WIDTH_MM` for a print plate, its own canvas for a page
    plate), so the measurement stays a real check rather than being skipped.

    *source* is the PDF the ``eps="raster"`` conversion reads; it is written to a
    scratch directory when the caller has none.  With ``verify=False`` nothing is
    measured.  Returns the measure of the written file, or ``None``.
    """
    if eps not in EPS_MODES:
        raise ValueError("eps=%r: expected one of %s"
                         % (eps, ", ".join(EPS_MODES)))
    path = os.path.abspath(path)
    note = eps
    if eps == "raster":
        scratch = None
        if source is None:              # the PDF is what this mode reads
            scratch = tempfile.mkdtemp(prefix="figure_style_")
            source = os.path.join(scratch, "plate.pdf")
            fig.savefig(source, format="pdf", facecolor="white",
                        metadata={"CreationDate": None, "ModDate": None})
        try:
            if not eps_from_pdf(source, path, timeout=timeout):
                note = "blend (no pdftops)"
        finally:
            if scratch:
                shutil.rmtree(scratch, ignore_errors=True)
    if not note.startswith("raster"):
        if note.startswith("blend"):
            blended = blend_alpha(fig)
            note += ", %d artist(s) blended" % len(blended)
            left = alpha_artists(fig)
            if left:
                raise RuntimeError(
                    "eps='blend' left %d artist(s) translucent (%s): the "
                    "PS backend would write their raw colour opaque -- "
                    "write this plate with eps='raster'"
                    % (len(left), ", ".join(sorted({
                        type(art).__name__ for art, _a, _m in left}))))
        fig.savefig(path, format="eps", facecolor="white")
    _normalise_eps(path)
    if not verify:
        return None
    info = assert_width(path, expected_mm=expect_mm, tol=tol)
    print("  eps %.2f mm  %s" % (info["mm"], note))
    return info


def save(fig, path, dpi=FIGURE_DPI, verify=True, vector=VECTOR_FORMATS,
         eps="blend"):
    """Write *path* at the declared canvas, then its vector sidecars.

    The canvas is **not** cropped to the drawn content.  ``bbox_inches="tight"``
    makes the file as wide as the layout happens to be, which is exactly the
    defect this module exists to stop: the reader then gets the ladder at a size
    nobody authored.  ``assert_fits`` runs first, so a layout that needs more
    than the canvas is refused instead of cropped; the width assertion then
    checks every file that was written.

    Beside the raster, under the same stem, the containers named in *vector*:

    * ``.pdf`` -- the plate as vector, transparency intact and byte-reproducible
      (no ``/CreationDate``, no ``/ModDate``); this is the copy for LaTeX;
    * ``.eps`` -- the same plate in PostScript, which has no alpha; how that is
      resolved is *eps*, one of :data:`EPS_MODES` (see :func:`blend_alpha` and
      :func:`eps_from_pdf`).  The sidecar itself is written by :func:`save_eps`,
      which a generator can also call on its own when its canvas is not the
      print width.

    Every container is measured, and a deviation from 170 +/- 0.5 mm aborts the
    run: a sidecar that is cropped or padded is as much a defect as a cropped
    raster.  Returns the raster path.
    """
    if eps not in EPS_MODES:
        raise ValueError("eps=%r: expected one of %s"
                         % (eps, ", ".join(EPS_MODES)))
    assert_fits(fig)
    path = os.path.abspath(path)
    stem = os.path.splitext(path)[0]
    fig.savefig(path, dpi=dpi, facecolor="white")
    pdf = None
    if "pdf" in vector:
        pdf = stem + ".pdf"
        fig.savefig(pdf, format="pdf", facecolor="white",
                    metadata={"CreationDate": None, "ModDate": None})
        if verify:
            info = assert_width(pdf)
            print("  %s %.2f mm  %5.0f kB  vector"
                  % (info["format"], info["mm"],
                     os.path.getsize(pdf) / 1024.0))
    if "eps" in vector:
        save_eps(fig, stem + ".eps", eps=eps, source=pdf, verify=verify)
    return path


def grid_rect(fig, bottom=None, top=None, pad_mm=2.0):
    """The ``tight_layout(rect=...)`` for a fixed canvas.

    ``tight_layout`` already reserves the height of a figure suptitle, so the
    top of the rect has to stay ``1.0``.  A rect that ends below the suptitle --
    the shape every generator carried while ``bbox_inches="tight"`` trimmed the
    canvas -- counts the caption twice and leaves a white band the height of the
    caption between the caption and the first row of panels; the cropping used
    to hide that band, a fixed canvas shows it.

    The rest is the part ``tight_layout`` does *not* know about: a
    ``fig.legend`` and a ``fig.text`` are figure artists and are in no axes'
    decorations, so they have to be declared, or they are drawn across the axis
    labels::

        leg = fig.legend(..., loc="lower center")
        fig.tight_layout(rect=style.grid_rect(fig, bottom=leg))
        fig.tight_layout(rect=style.grid_rect(fig, top=note))   # a header strip

    *bottom* and *top* are such an artist, an explicit height in millimetres, or
    ``None``.  Returns ``(0.0, bottom_fraction, 1.0, top_fraction)``.
    """
    height_mm = fig.get_size_inches()[1] * 25.4

    def edge(artist, which):
        if artist is None:
            return None
        if not hasattr(artist, "get_window_extent"):
            return float(artist)
        fig.canvas.draw()
        box = artist.get_window_extent(fig.canvas.get_renderer())
        return getattr(box, which) * 25.4 / fig.dpi

    low = edge(bottom, "y1")
    low = pad_mm if low is None else low + pad_mm
    high = edge(top, "y0")
    high = height_mm - pad_mm if high is None else high - pad_mm
    return (0.0, min(1.0, max(0.0, low / height_mm)), 1.0,
            min(1.0, max(0.0, high / height_mm)))


def _main(argv):
    """Print the ladder, then measure the files named on the command line."""
    for role in ("base", "title", "suptitle", "label", "tick", "legend",
                 "annotation"):
        band = BANDS.get(role)
        if band and band[1]:
            spec = "%.1f-%.1f pt" % band
        elif band:
            spec = ">= %.1f pt" % band[0]
        else:
            spec = "anchor"
        print("%-11s %5.1f pt   %s" % (role, FS[role], spec))
    print("%-11s %5.1f mm   (%.3f in, tol %.1f mm)"
          % ("width", FINAL_WIDTH_MM, COL_WIDTH, WIDTH_TOL_MM))
    for arg in argv[1:]:
        info = measure(arg)
        got = info["mm"]
        if info["vector"]:
            how = "%s box" % info["format"]
        else:
            how = "dpi=%s" % (("%.0f" % info["dpi"]) if info["dpi"] else "?")
        print("%-58s %7.1f mm  %s"
              % (os.path.relpath(arg), got if got else float("nan"), how))
    return 0


if __name__ == "__main__":
    sys.exit(_main(sys.argv))

