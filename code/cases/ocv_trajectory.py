#!/usr/bin/env python3
"""One standard OCV **trajectory** plate per case site, on the support-sized chip.

The question, and why it needs its own template
-----------------------------------------------
Every other figure of this repository answers a *contrast* question ("does the
state differ?"). This one answers the **timeline** question the six cases are
actually read for: *where in the record does the vector move, and does it move
where the event is?* It is deliberately the **same plate for every site** -- same
twelve panels, same standardisation, same axes logic -- because the comparison a
reader makes is between plates, and a per-site design would make that comparison a
comparison of designs.

The reference timeline is defined **case by case** (see
``code/cases/ocv_trajectory_spec.md``), and the post-event tail is **drawn**, not
cut off: a trajectory that stops at the event cannot show whether the vector
returned, stayed or kept moving.

The five rules the template enforces
------------------------------------
* **R1 -- one cut.** The plate reads the **support-sized** chip
  (``data/<site>/support/<site>_ocv_channels.csv``), the same cut the shared
  orbit plate and the SDS reconfiguration read. The committed chip is never mixed
  in: two chips are two measurements.
* **R2 -- the case's own reference, and the window is chosen per asset.** The
  reference timeline is the case's own reference state (``pre`` /
  ``pre-collapse (healthy)`` / *all* for the event-free control) **restricted to
  the window the registry declares for that site**
  (``SITES[...]["reference_window"]``, inclusive calendar months; ``None`` means
  "the state's own span"). The post tail is every chip of the post state. A
  window is a **choice**, so the plate draws what the choice leaves out: a
  healthy chip outside the window is a third class -- drawn, in the reference
  frame, and standardising nothing.
* **R3 -- per-block standardisation.** The two response blocks of
  ``ocv_reconfigure.BLOCKS`` are standardised on **their own** reference sets:
  the behaviour block (``gamma2, P``) against the reference chips that **carry an
  echo** (the two channels only exist where the mask rule holds), the geometry
  block (``A, D, F, S``) against **every** reference chip. Pooled over the orbits,
  except Espoo, whose only contrast is the orbit geometry and which is therefore
  standardised **per orbit**.
* **R4 -- a robust z.** Transform first (log for ``gamma2, P, A, D, F``, identity
  for ``S``, exactly ``ocv_reconfigure.TRANSFORM``), then
  ``z = (x - median_ref) / (1.4826 * MAD_ref)``, falling back to the IQR and then
  to the standard deviation if the MAD is zero, and clipping at ``+-6``. A block
  magnitude is the **RMS** of its standardised dimensions, as
  ``ocv_reconfigure.magnitude`` defines it.
* **R5 -- no interpolation.** A month with no chip is a hole: the line breaks
  there and never spans it, and neither does the event. The monthly median is
  drawn where a month holds a chip, with its ``n``.
* **R6 -- the outside class is tested too.** Healthy chips of the reference state
  that the chosen window leaves out are drawn (grey, hollow) and read against the
  reference with the same interval, so the plate never shows a chip it does not
  also measure. Where the window is the state's own span, the block says so and
  the report prints the reason instead of the table.

The twelve panels
-----------------
``a1`` z(gamma2) -- ``a2`` z(P) -- ``a3`` behaviour plane -- ``a4`` R_B
``b1`` z(A) -- ``b2`` z(F) -- ``b3`` geometry plane (``D``, ``S`` annotated)
``b4`` R_G -- ``c1`` weather -- ``c2`` observation cadence -- ``c3`` echo
coverage -- ``c4`` the trajectory matrix (all six dimensions, one row each).

What it does not do
-------------------
No classifier, no damage verdict, no cross-site pooling, no significance claim on
the curve. Both intervals of the plate are off the curve: the **contrast** (post
minus reference) and R6's **outside class** (the excluded healthy chips minus the
reference), each from a seeded circular moving-block bootstrap over the chip
series (``BOOT`` draws, ``BLOCK`` chips, seed ``SEED``) -- the same bootstrap
family as the orbit plate, so the interval means the same thing.

Outputs
-------
``figures/<site>/fig_<site>_ocv_trajectory.jpg``   630 dpi
``figures/<site>/fig_<site>_ocv_trajectory.md``    the report
``data/<site>/fig_<site>_ocv_trajectory.json``     every plotted number

``--verify`` re-derives every number from the table and fails if the JSON and the
table disagree; ``--selftest`` checks the rules above without touching a file.

Usage:
  python3 code/cases/ocv_trajectory.py --site carola
  python3 code/cases/ocv_trajectory.py --all
  python3 code/cases/ocv_trajectory.py --selftest
  python3 code/cases/ocv_trajectory.py --site carola --verify
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import os
import random
import sys
import tempfile

import matplotlib

matplotlib.use("Agg")

import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, os.pardir, os.pardir))
DATA = os.path.join(REPO, "data")

SCHEMA = "ocv-trajectory-v2"
FIGURE_DPI = 630
BOOT = 2000
BLOCK = 3
SEED = 7
Z_CLIP = 6.0
MM = 1.0 / 25.4
COL_WIDTH = 183 * MM

#: The frozen pins of the held-out assessment. They are constants, not choices
#: made after looking at 2024: the seed, the block length and the quantile
#: levels the plate already used, and a persistence rule fixed before the
#: held-out numbers were read. ``--selftest`` asserts them, so a later edit
#: cannot silently move the criterion.
QUANTILES = (5, 50, 95)
DECISION_QUANTILE = 95
PERSIST_MIN = 2
SUMMER_2024 = ("2024-06-01", "2024-08-31")

#: The names the report and the plate use for the three frozen classes. The
#: *analysis* names stay ``reference`` / ``outside`` / ``post`` in the JSON.
CLASS_LABEL = {
    "reference": "frozen pre-event reference",
    "outside": "held-out pre-event",
    "post": "post-event",
}
CLASS_NOTE = {
    "reference": "the frozen pre-event reference: the chips the standardisation "
                 "rests on (R2/R3)",
    "outside": "held-out pre-event: reference-state chips after the reference "
               "window and strictly before the event; drawn in the reference "
               "frame, standardising nothing",
    "post": "post-event: every chip of the registry's post state(s)",
}

#: The response blocks, verbatim from ``ocv_reconfigure.BLOCKS`` (asserted by
#: ``--selftest`` whenever that module can be imported, so the two can never drift
#: apart silently).
BLOCKS = {
    "behaviour": ("gamma2", "P"),
    "geometry": ("A", "D", "F", "S"),
}
GEOMETRY_PLANE = ("A", "F")
DIMS = ("gamma2", "P", "A", "D", "F", "S")
DIM_LABEL = {"gamma2": r"$\gamma^2$", "P": "P", "A": "A", "D": "D",
             "F": "F", "S": "S"}
DIM_LONG = {"gamma2": "coherence consistency", "P": "echo persistence",
            "A": "mask extent", "D": "mask occupancy",
            "F": "mask fragmentation", "S": "mask shift"}
#: The same names in plain markdown (the figure uses mathtext, a table does not).
DIM_MD = {"gamma2": "γ²", "P": "P", "A": "A", "D": "D", "F": "F", "S": "S"}

#: Which block a dimension belongs to; a dimension that is in neither block would
#: be a registry error.
BLOCK_OF = {dim: name for name, dims in BLOCKS.items() for dim in dims}


#: The six sites. ``reference`` is the state the standardisation is taken from
#: (``None`` means "no state column": every chip is the reference, the event-free
#: control), ``reference_window`` the **chosen** slice of that state as an
#: inclusive range of calendar months (``None`` = the state's own span, first
#: chip to last), ``post`` are the states drawn as the post-event tail and
#: ``excluded`` the states the registry removes before any statistic, with the
#: reason it gives.
SITES = {
    "carola": {
        "label": "Carolabrücke, Dresden",
        "event": "2024-09-11",
        "reference": "pre-collapse (healthy)",
        # The reference of record: the first sixteen healthy months. The rest of
        # the healthy state (2024-01 .. 2024-09, 63 chips) is the record's own
        # pre-collapse drift -- drawn, but standardising nothing (R2).
        "reference_window": ("2022-09", "2023-12"),
        "states": ("pre-collapse (healthy)", "post-collapse"),
        "post": ("post-collapse",),
        "excluded": {},
        "weather": ("temperature_c", "wind_speed_ms"),
        "note": "28 months of reference before the collapse, 24 months of tail "
                "after it: the longest pair of the six, and the only one whose "
                "post side is longer than many whole records here.",
    },
    "cts": {
        "label": "Champlain Towers South, Surfside",
        "event": "2021-06-24",
        "reference": "pre-collapse (healthy)",
        "reference_window": None,
        "states": ("pre-collapse (healthy)", "post-collapse"),
        "post": ("post-collapse",),
        "excluded": {},
        "weather": ("temperature_c", "precipitation_mm"),
        "note": "one orbit only (ASCENDING, relative orbit 48): the orbit columns "
                "are empty on the other direction, so R3's pooling is the whole "
                "record here. The reference is the site's six-year pre-trend.",
    },
    "morandi": {
        "label": "Ponte Morandi, Genova",
        "event": "2018-08-14",
        "reference": "pre-collapse (healthy)",
        "reference_window": None,
        "states": ("pre-collapse (healthy)", "post-collapse"),
        "post": ("post-collapse",),
        "excluded": {},
        "weather": ("temperature_c", "wind_speed_ms"),
        "note": "nine post chips on thirteen days: the tail is drawn and its "
                "contrast is reported with the interval it can carry, which is "
                "wide by construction.",
    },
    "kdlo": {
        "label": "KDLO media tower, Garden City SD",
        "event": "2022-12-11",
        "reference": "pre",
        "reference_window": None,
        "states": ("pre", "single", "rebuild"),
        "post": ("rebuild",),
        "excluded": {"single": "one acquisition at the collapse itself "
                               "(2022-12-23): one unit is not a state"},
        "weather": ("temperature_c", "wind_speed_ms"),
        "note": "the tail is the **rebuild**, not a collapse aftermath: the "
                "contrast answers a different question, and the plate says so in "
                "its own subtitle rather than in a footnote.",
    },
    "haltern": {
        "label": "Haltern AV9 wind turbine",
        "event": "2021-09-29",
        "reference": "pre-collapse (healthy)",
        "reference_window": None,
        "states": ("pre-collapse (healthy)", "post-collapse"),
        "post": ("post-collapse",),
        "excluded": {},
        "weather": ("temperature_c", "wind_speed_ms"),
        "note": "13 of the record's 120 dates are refused at the raster edge, all "
                "ascending: the cut's own manifest carries them, and the plate "
                "does not pretend the missing 11 % are months without data.",
    },
    "espoo": {
        "label": "Espoo Kurttila mast (event-free control)",
        "event": None,
        "reference": None,
        "reference_window": None,
        "states": ("all",),
        "post": (),
        "excluded": {},
        "weather": ("temperature_c", "wind_speed_ms"),
        "note": "no event to anchor on and no state column: the whole record is "
                "the reference *and* the trajectory, standardised **per orbit** "
                "(the only contrast the record carries), with no contrast row.",
    },
}

#: The transform per dimension: ``ocv_reconfigure.TRANSFORM`` says exactly this
#: for the sites it registers (Carola, Morandi, Haltern, KDLO) -- log for every
#: strictly positive dimension, identity for the non-negative shift ``S``.
TRANSFORM = {dim: ("identity" if dim == "S" else "log") for dim in DIMS}


def table_path(site):
    return os.path.join(DATA, site, "support", site + "_ocv_channels.csv")


def out_paths(site):
    return {"json": os.path.join(DATA, site, "fig_%s_ocv_trajectory.json" % site),
            "md": os.path.join(REPO, "figures", site,
                               "fig_%s_ocv_trajectory.md" % site),
            "jpg": os.path.join(REPO, "figures", site,
                                "fig_%s_ocv_trajectory.jpg" % site)}




# ---------------------------------------------------------------------------
# Small statistics, spelled out (no numpy, no scipy: every number in the JSON is
# reproducible on a bare interpreter, and the moving-block bootstrap below is the
# only randomness in the file)
# ---------------------------------------------------------------------------


def median(values):
    """Median of a non-empty sequence (mean of the two middles when even)."""
    ordered = sorted(values)
    n = len(ordered)
    if n == 0:
        raise ValueError("median of an empty sequence")
    mid = n // 2
    if n % 2:
        return ordered[mid]
    return 0.5 * (ordered[mid - 1] + ordered[mid])


def percentile(values, q):
    """Linear-interpolation percentile (the type-7 rule numpy defaults to)."""
    ordered = sorted(values)
    n = len(ordered)
    if n == 0:
        raise ValueError("percentile of an empty sequence")
    if n == 1:
        return ordered[0]
    pos = (n - 1) * (q / 100.0)
    low = int(math.floor(pos))
    high = min(low + 1, n - 1)
    frac = pos - low
    return ordered[low] * (1.0 - frac) + ordered[high] * frac


def stdev(values):
    """Sample standard deviation (``n - 1``), the last fallback of the scale."""
    n = len(values)
    if n < 2:
        return 0.0
    mean = sum(values) / n
    return math.sqrt(sum((v - mean) ** 2 for v in values) / (n - 1))


def scale_of(values):
    """``(median, scale, source)`` -- the robust centre and scale of a series.

    MAD first (``1.4826 * MAD``, consistent for a normal sample), then the IQR
    divided by 1.349, then the sample standard deviation, then 1.0 with the source
    named ``degenerate`` -- a reference that does not move at all cannot
    standardise anything, and the plate says so instead of dividing by zero.
    """
    centre = median(values)
    mad = median([abs(v - centre) for v in values])
    scale = 1.4826 * mad
    source = "mad"
    if scale <= 0.0:
        scale = (percentile(values, 75) - percentile(values, 25)) / 1.349
        source = "iqr"
    if scale <= 0.0:
        scale = stdev(values)
        source = "sd"
    if scale <= 0.0:
        scale = 1.0
        source = "degenerate"
    return centre, scale, source


def apply_transform(value, how):
    """``log`` or ``identity``; a non-positive value under ``log`` is an error."""
    if value is None:
        return None
    if how == "log":
        if value <= 0.0:
            raise ValueError("log transform of a non-positive value: %r" % (value,))
        return math.log(value)
    return value


def month_key(day):
    return (day.year, day.month)


def month_date(key):
    """The 15th of the month -- a tick mark, never a data coordinate."""
    return dt.date(key[0], key[1], 15)


def next_month(key):
    return (key[0] + 1, 1) if key[1] == 12 else (key[0], key[1] + 1)


def month_text(key):
    return "%04d-%02d" % key


def parse_day(text):
    """The calendar day of an ISO timestamp in any of the shapes the tables use."""
    text = (text or "").strip()
    if not text:
        return None
    if text.endswith("Z"):
        text = text[:-1] + "+00:00"
    try:
        stamp = dt.datetime.fromisoformat(text)
    except ValueError:
        try:
            return dt.date.fromisoformat(text[:10])
        except ValueError:
            return None
    return stamp.date()


def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 16), b""):
            digest.update(chunk)
    return digest.hexdigest()


def round_floats(value, digits=10):
    """Round every float of a nested structure, so two runs compare exactly."""
    if isinstance(value, float):
        if value != value or value in (float("inf"), float("-inf")):
            return None
        return round(value, digits)
    if isinstance(value, dict):
        return {k: round_floats(v, digits) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [round_floats(v, digits) for v in value]
    return value


# ---------------------------------------------------------------------------
# Read, partition, standardise
# ---------------------------------------------------------------------------

#: The family's threshold for "this chip carries an echo at all" (``A >= 2``); a
#: chip below it has no behaviour channels, which is exactly why R3 gives the
#: behaviour block its own reference set.
ECHO_MIN = 2.0
ORBITS = ("ASCENDING", "DESCENDING")
ORBIT_SHORT = {"ASCENDING": "ASC", "DESCENDING": "DESC"}
CHANNEL_COLUMN = {"gamma2": "gamma2_mask"}


def rms(values):
    """Root mean square of the available standardised responses.

    ``ocv_reconfigure.magnitude`` verbatim: the block magnitude divides by the
    number of dimensions it actually had, so a block of three is not
    automatically larger than a block of two.
    """
    usable = [v for v in values if v is not None]
    if not usable:
        return None
    return math.sqrt(sum(v * v for v in usable) / len(usable))


def reference_label(site):
    """The state (or the word ``all``) the standardisation is taken from."""
    spec = SITES[site]
    return spec["reference"] if spec["reference"] else "all"


#: Sentinel for "nothing was overridden before this one".
_MISSING = object()

#: Reference windows substituted for the duration of one build, keyed by site.
#: Empty in every ordinary run: the registry is then the only source.
WINDOW_OVERRIDE = {}

#: The rule that selected an overridden window, keyed by site -- a sentence that
#: is carried into the report, never a number that computes anything.
WINDOW_RULE = {}


def reference_window(site):
    """The frozen reference window in force for ``site`` (R2).

    Normally the registry's own choice, which is what the committed plates say.
    ``window_override`` substitutes another for the length of one build; the
    prospective runner freezes the window its stated rule selected that way, so
    the registry -- and therefore the descriptive plate -- keeps saying what it
    always said, and this accessor is the only place the two can differ.
    """
    if site in WINDOW_OVERRIDE:
        return WINDOW_OVERRIDE[site]
    return SITES[site]["reference_window"]


def window_rule(site):
    """The sentence that selected the window in force, or ``None``."""
    return WINDOW_RULE.get(site)


class window_override(object):
    """Build one site against ``window`` (and ``rule``), then restore the registry.

    A context manager, not a set-and-forget: a caller that wants another window
    has it for exactly one ``with`` block, and a failing build cannot leave the
    registry altered for the next case in the same process.
    """

    def __init__(self, site, window, rule=None):
        self.site = site
        self.window = window
        self.rule = rule
        self.previous = _MISSING
        self.previous_rule = _MISSING

    def __enter__(self):
        self.previous = WINDOW_OVERRIDE.get(self.site, _MISSING)
        self.previous_rule = WINDOW_RULE.get(self.site, _MISSING)
        WINDOW_OVERRIDE[self.site] = self.window
        if self.rule:
            WINDOW_RULE[self.site] = self.rule
        return self.window

    def __exit__(self, *exception):
        if self.previous is _MISSING:
            WINDOW_OVERRIDE.pop(self.site, None)
        else:
            WINDOW_OVERRIDE[self.site] = self.previous
        if self.previous_rule is _MISSING:
            WINDOW_RULE.pop(self.site, None)
        else:
            WINDOW_RULE[self.site] = self.previous_rule
        return False


def within(day, window):
    """True when ``day`` falls in an inclusive calendar-**month** window (R2).

    ``window`` is ``(first, last)`` as ``"YYYY-MM"``, or ``None`` for "no
    restriction". The month is the unit because a month is the plate's own bin
    (R5): a window picked by the day could cut a month in half, and then a month
    would sit on both sides of the choice the plate is supposed to show.
    """
    if window is None:
        return True
    key = month_text(month_key(day))
    return window[0] <= key <= window[1]


def reference_state(site, chip):
    """True for a chip of the reference *state* -- the chosen window not applied."""
    reference = SITES[site]["reference"]
    return True if reference is None else chip["state"] == reference


def in_reference(site, chip):
    """True for a reference chip: the reference state inside the chosen window.

    The window is per asset (R2). Every statistic of the plate -- the centre and
    the scale, the z of every chip, the blocks, the weather reference and the
    reference side of the contrast -- is taken on this set and on no other.
    """
    return reference_state(site, chip) and within(
        chip["day"], reference_window(site))


def event_day(site):
    """The documented event as a date, or ``None`` on the event-free control."""
    text = SITES[site]["event"]
    return dt.date.fromisoformat(text) if text else None


def outside_chips(site, chips):
    """The OUTSIDE class: held-out pre-event chips.

    A chip is held out when it belongs to the reference **state** but not to the
    frozen reference **window**, and it is acquired **strictly before** the
    documented event. The strict day rule keeps the class on one side of the
    event: a reference-state chip on or after the event date belongs to neither
    class, is counted in the audit and is never silently absorbed into OUTSIDE.

    The window, where the asset has one, also bounds the class from the left: a
    reference-state chip acquired *before* the frozen cycle was never held out
    from anything, so it is the audit's ``prefix`` bucket and never OUTSIDE.
    """
    event = event_day(site)
    frozen = reference_window(site)
    held = []
    for chip in chips:
        if not reference_state(site, chip) or in_reference(site, chip):
            continue
        if event is not None and chip["day"] >= event:
            continue
        if frozen is not None and month_text(month_key(chip["day"])) < frozen[0]:
            continue
        held.append(chip)
    return held


def chip_class(site, chip):
    """``reference`` / ``outside`` / ``prefix`` / ``post`` / ``unassigned``.

    The three classes are fixed before any statistic is computed: REFERENCE is
    the reference state inside the frozen window, OUTSIDE what that window leaves
    out on the pre-event side, POST the registry's post state(s). The event day
    is handled explicitly: a post-labelled chip on the event day is POST (the
    registry's own label), a reference-state chip on or after the event is
    *unassigned* and reported, never absorbed.
    """
    if in_reference(site, chip):
        return "reference"
    if reference_state(site, chip) and not in_reference(site, chip):
        event = event_day(site)
        frozen = reference_window(site)
        if event is None or chip["day"] < event:
            if frozen is not None and month_text(
                    month_key(chip["day"])) < frozen[0]:
                return "prefix"
            return "outside"
        return "post" if chip["state"] in SITES[site]["post"] else "unassigned"
    if chip["state"] in SITES[site]["post"]:
        return "post"
    return "unassigned"


def reference_split(site, chips):
    """``(reference, outside, state)`` -- the choice, what it leaves out, the state.

    ``outside`` is the third class of R2: chips of the reference *state* that the
    chosen window leaves out. They are drawn in the reference frame and they
    standardise nothing -- that is the whole difference between them and the
    reference, and the plate says so in its legend.
    """
    window = reference_window(site)
    state = [chip for chip in chips if reference_state(site, chip)]
    ref = [chip for chip in state if within(chip["day"], window)]
    out = outside_chips(site, chips)
    return ref, out, state


def group_label(site, chip):
    """The reference group a chip is standardised with (R3).

    Pooled over the orbits everywhere except Espoo, whose only contrast *is* the
    pass geometry -- there the two orbits are their own reference sets.
    """
    return chip["orbit"] if site == "espoo" else "pooled"


def load_chips(site):
    """Every chip of the site's support cut, with the audit of what was dropped."""
    path = table_path(site)
    if not os.path.exists(path):
        raise SystemExit("%s: no support table %s"
                         % (site, os.path.relpath(path, REPO)))
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle))
    excluded = {}
    refused = []
    chips = []
    for row in rows:
        state = (row.get("state") or "").strip() or None
        if state is not None and state in SITES[site]["excluded"]:
            excluded[state] = excluded.get(state, 0) + 1
            continue
        day = parse_day(row.get("acquisition_ts") or row.get("date"))
        if day is None:
            refused.append({"id": row.get("id"), "reason": "no acquisition date"})
            continue
        orbit = (row.get("orbit_direction") or "").strip().upper()
        if orbit not in ORBITS:
            refused.append({"id": row.get("id"),
                            "reason": "orbit_direction %r" % (orbit,)})
            continue
        values = {}
        for dim in DIMS:
            raw = (row.get(CHANNEL_COLUMN.get(dim, dim)) or "").strip()
            values[dim] = float(raw) if raw else None
        weather = {}
        for column in SITES[site]["weather"]:
            raw = (row.get(column) or "").strip()
            weather[column] = float(raw) if raw else None
        chips.append({"day": day, "orbit": orbit, "state": state,
                      "values": values, "weather": weather,
                      "echo": (values.get("A") or 0.0) >= ECHO_MIN})
    chips.sort(key=lambda chip: (chip["day"], chip["orbit"]))
    audit = {"rows": len(rows), "kept": len(chips), "excluded": excluded,
             "refused": refused}
    return chips, audit


def reference_groups(site, chips):
    """The reference sets R3 standardises against, as ``{label: [chips]}``."""
    refs = [chip for chip in chips if in_reference(site, chip)]
    if site == "espoo":
        # Espoo's only contrast *is* the pass geometry (ASCENDING n = 31,
        # DESCENDING n = 44), so its reference is per orbit and never pooled.
        return {orbit: [c for c in refs if c["orbit"] == orbit] for orbit in ORBITS}
    return {"pooled": refs}


def standardise(site, chips):
    """Centre and scale per dimension x reference group, and the z of every chip."""
    stats = {}
    for label, group in reference_groups(site, chips).items():
        for dim in DIMS:
            block = BLOCK_OF[dim]
            if block == "behaviour":
                usable = [c for c in group
                          if c["echo"] and c["values"][dim] is not None]
                rule = "echo-bearing reference chips"
            else:
                usable = [c for c in group if c["values"][dim] is not None]
                rule = "all reference chips"
            if not usable:
                raise SystemExit("%s: no reference chip carries %s (%s)"
                                 % (site, dim, label))
            values = [apply_transform(c["values"][dim], TRANSFORM[dim])
                      for c in usable]
            centre, scale, source = scale_of(values)
            stats[(dim, label)] = {
                "transform": TRANSFORM[dim], "rule": rule, "n": len(values),
                "centre": centre, "scale": scale, "source": source,
                "z_min": (min(values) - centre) / scale,
                "z_max": (max(values) - centre) / scale,
            }
    weather = weather_reference(site, chips)
    clipped = 0
    for chip in chips:
        label = group_label(site, chip)
        z = {}
        for dim in DIMS:
            cell = stats[(dim, label)]
            raw = chip["values"][dim]
            if raw is None:
                z[dim] = None
                continue
            value = (apply_transform(raw, cell["transform"]) - cell["centre"]) \
                / cell["scale"]
            if value > Z_CLIP or value < -Z_CLIP:
                clipped += 1
                value = max(-Z_CLIP, min(Z_CLIP, value))
            z[dim] = value
        chip["z"] = z
        chip["R_B"] = rms([z[d] for d in BLOCKS["behaviour"]])
        chip["R_G"] = rms([z[d] for d in BLOCKS["geometry"]])
        weather_z = {}
        for column, per_group in weather.items():
            cell = per_group.get(label)
            raw = chip["weather"].get(column)
            if cell is None or raw is None:
                weather_z[column] = None
                continue
            weather_z[column] = max(-Z_CLIP, min(Z_CLIP,
                                                (raw - cell["centre"])
                                                / cell["scale"]))
        chip["weather_z"] = weather_z
    return stats, weather, clipped


def weather_reference(site, chips):
    """The same robust z for the record's own co-variates (identity scale).

    Panel c1 draws them, so they are standardised exactly like a dimension --
    against the reference chips, on the identity scale they come in.
    """
    out = {}
    for column in SITES[site]["weather"]:
        per_group = {}
        for label, group in reference_groups(site, chips).items():
            values = [c["weather"][column] for c in group
                      if c["weather"].get(column) is not None]
            if not values:
                continue
            centre, scale, source = scale_of(values)
            per_group[label] = {"transform": "identity", "n": len(values),
                                "distinct": len(set(values)),
                                "centre": centre, "scale": scale,
                                "source": source}
        if per_group:
            out[column] = per_group
    return out



def split_event_month(site, chips):
    """The event month that carries chips on both sides, else ``None``.

    A calendar month is the plate's bin (R5), but a bin may not combine
    observations from both sides of the event. The event month is split into a
    ``pre`` and a ``post`` bin exactly where the record holds chips on both
    sides; a month that holds only one side stays a plain month.
    """
    event = event_day(site)
    if event is None:
        return None
    month = month_key(event)
    days = [chip["day"] for chip in chips if month_key(chip["day"]) == month]
    if any(day < event for day in days) and any(day >= event for day in days):
        return month
    return None


def bin_key(entry):
    """``(year, month, side)`` -- side is ``whole`` / ``pre`` / ``post``."""
    return (int(entry["month"][:4]), int(entry["month"][5:]),
            entry.get("part") or "whole")


def bin_date(entry):
    """The x anchor of a bin: 5th (pre), 25th (post), 15th (whole)."""
    year, month, side = bin_key(entry)
    day = 5 if side == "pre" else 25 if side == "post" else 15
    return dt.date(year, month, day)


def follows(left, right):
    """True when ``right`` is the next occupied bin after ``left`` in time."""
    lkey, rkey = bin_key(left), bin_key(right)
    if lkey[:2] == rkey[:2]:
        return lkey[2] == "pre" and rkey[2] == "post"
    year, month = (lkey[0] + 1, 1) if lkey[1] == 12 else (lkey[0], lkey[1] + 1)
    return rkey[0] == year and rkey[1] == month


def bin_value(entry, key):
    """The value of one quantity in one monthly bin (z for a dimension)."""
    return entry["z"][key] if key in DIMS else entry[key]


def monthly(site, chips):
    """One row per calendar-month bin that holds a chip, in date order.

    The event month is split where the record holds chips on both sides of the
    event (``split_event_month``), so no bin combines pre- and post-event
    acquisitions. Each row carries the census of its bin, the monthly medians the
    plate draws, and the class composition the report tabulates.
    """
    split = split_event_month(site, chips)
    event = event_day(site)
    groups = {}
    for chip in chips:
        key = month_key(chip["day"])
        part = None
        if split is not None and key == split:
            part = "pre" if chip["day"] < event else "post"
        groups.setdefault((key, part), []).append(chip)
    order = {"whole": 0, "pre": 0, "post": 1}
    series = []
    for key, part in sorted(groups, key=lambda item: (
            item[0][0], item[0][1], order[item[1] or "whole"])):
        group = groups[(key, part)]
        label = month_text(key) + ("" if part is None else " (%s)" % part)
        entry = {"month": month_text(key), "label": label, "part": part,
                 "anchor": bin_date({"month": month_text(key),
                                     "part": part}).isoformat(),
                 "n": len(group),
                 "n_reference": sum(1 for c in group
                                    if chip_class(site, c) == "reference"),
                 "n_outside": sum(1 for c in group
                                  if chip_class(site, c) == "outside"),
                 "n_post": sum(1 for c in group
                               if chip_class(site, c) == "post"),
                 "n_unassigned": sum(1 for c in group
                                     if chip_class(site, c) == "unassigned"),
                 "n_prefix": sum(1 for c in group
                                 if chip_class(site, c) == "prefix"),
                 "n_ascending": sum(1 for c in group
                                    if c["orbit"] == "ASCENDING"),
                 "n_descending": sum(1 for c in group
                                     if c["orbit"] == "DESCENDING"),
                 "n_echo": sum(1 for c in group if c["echo"]),
                 "z": {}, "weather_z": {}}
        for dim in DIMS:
            values = [c["z"][dim] for c in group if c["z"][dim] is not None]
            entry["z"][dim] = median(values) if values else None
        for name in ("R_B", "R_G"):
            values = [c[name] for c in group if c[name] is not None]
            entry[name] = median(values) if values else None
        for column in group[0]["weather_z"]:
            values = [c["weather_z"][column] for c in group
                      if c["weather_z"].get(column) is not None]
            entry["weather_z"][column] = median(values) if values else None
        # Only a rule-selected window can leave a bin before it, and a key that
        # is always zero would change every committed JSON for nothing.
        if not entry["n_prefix"]:
            del entry["n_prefix"]
        series.append(entry)
    return series


def quantile_block(values):
    """The empirical q05 / q50 / q95 summary of one reference distribution."""
    return {"n": len(values), "q05": percentile(values, 5),
            "q50": percentile(values, 50), "q95": percentile(values, 95),
            "min": min(values), "max": max(values)}


def reference_envelope(site, chips, months):
    """The reference distributions, at CHIP and at MONTHLY aggregation.

    A monthly held-out state may only be compared with a monthly reference state
    (the like-with-like rule), so the decision distribution is the distribution
    of the **monthly medians of the reference months** -- never of single
    reference chips and never the reference maximum. ``q05``/``q50``/``q95`` are
    empirical type-7 quantiles of those monthly states.
    """
    def series(key):
        pairs = []
        for month in months:
            if not month["n"] or month["n_reference"] != month["n"]:
                continue
            value = bin_value(month, key)
            if value is not None:
                pairs.append({"month": month["label"], "value": value})
        return pairs

    monthly_dist = {}
    for key in list(DIMS) + ["R_B", "R_G"]:
        pairs = series(key)
        if not pairs:
            continue
        block = quantile_block([pair["value"] for pair in pairs])
        block["values"] = pairs
        block["role"] = ("decision distribution: a magnitude has no sign, so only "
                         "an exceedance above q95 counts"
                         if key in ("R_B", "R_G") else
                         "directional reference distribution: outside "
                         "[q05, q95] is a departure from the reference's own "
                         "monthly range")
        monthly_dist[key] = block

    chip_dist = {}
    reference_chips = [c for c in chips if in_reference(site, c)]
    for key in list(DIMS) + ["R_B", "R_G"]:
        if key in DIMS:
            values = [c["z"][key] for c in reference_chips
                      if c["z"][key] is not None]
        else:
            values = [c[key] for c in reference_chips if c[key] is not None]
        if values:
            chip_dist[key] = quantile_block(values)
    labels = [pair["month"] for pair in
              monthly_dist.get("gamma2", {}).get("values", [])]
    return {"aggregation": "monthly median of each calendar-month bin",
            "quantiles": list(QUANTILES),
            "decision_quantile": DECISION_QUANTILE,
            "decision_level": "monthly",
            "persistence_min_months": PERSIST_MIN,
            "n_reference_months": len(labels),
            "months": labels,
            "monthly": monthly_dist,
            "chip": chip_dist,
            "chip_note": "descriptive only: a single chip is never compared "
                         "with a monthly state in this analysis"}


def envelope_position(entries, key, distribution, one_sided=False):
    """Where a class's monthly states sit in the reference's monthly envelope."""
    values = [(entry["label"], bin_value(entry, key)) for entry in entries]
    values = [(label, value) for label, value in values if value is not None]
    above = [label for label, value in values if value > distribution["q95"]]
    below = [label for label, value in values if value < distribution["q05"]]
    level = median([value for _label, value in values]) if values else None
    if level is None:
        position = "no monthly state"
    elif level > distribution["q95"]:
        position = "above q95"
    elif level < distribution["q05"]:
        position = "below q05"
    else:
        position = "inside q05..q95"
    return {"level": level, "q05": distribution["q05"],
            "q50": distribution["q50"], "q95": distribution["q95"],
            "position": position, "months_above_q95": above,
            "months_below_q05": below, "n_months": len(values),
            "n_above_q95": len(above), "n_below_q05": len(below),
            "max_value": max([value for _label, value in values]) if values else None,
            "min_value": min([value for _label, value in values]) if values else None,
            "outside_envelope": above if one_sided else sorted(above + below),
            "criterion": ("above q95 (one-sided: a magnitude has no sign)"
                          if one_sided else
                          "outside [q05, q95] (two-sided signed departure)")}


def persistence_runs(entries, key, distribution):
    """Consecutive held-out months above the monthly reference q95."""
    threshold = distribution["q95"]
    flags = [(entry, bin_value(entry, key)) for entry in entries]
    flags = [(entry, value) for entry, value in flags if value is not None]
    runs, current, previous = [], [], None
    for entry, value in flags:
        above = value > threshold
        if above and current and previous is not None and follows(previous, entry):
            current.append(entry)
        elif above:
            current = [entry]
        else:
            if current:
                runs.append(current)
            current = []
        previous = entry
    if current:
        runs.append(current)
    blocks = [{"labels": [e["label"] for e in run], "start": run[0]["label"],
               "end": run[-1]["label"], "length": len(run),
               "values": [bin_value(e, key) for e in run]} for run in runs]
    longest = max(blocks, key=lambda block: block["length"]) if blocks else None
    return {"quantity": key, "threshold": threshold,
            "threshold_quantile": DECISION_QUANTILE, "n_months": len(flags),
            "n_above_q95": sum(1 for _e, value in flags if value > threshold),
            "months_above_q95": [e["label"] for e, value in flags
                                 if value > threshold],
            "runs": blocks, "longest_run": longest,
            "sustained": bool(longest and longest["length"] >= PERSIST_MIN),
            "persistence_min_months": PERSIST_MIN}


def directional_runs(entries):
    """Runs of the joint behaviour pattern ``z(gamma2) > 0 and z(P) < 0``.

    Descriptive only: run labels and lengths, no threshold, no test, and the run
    is never turned into a significance statistic.
    """
    pattern = "z(gamma2) > 0 and z(P) < 0"
    flags = []
    for entry in entries:
        gamma2, persistence = entry["z"]["gamma2"], entry["z"]["P"]
        flags.append((entry, bool(gamma2 is not None and persistence is not None
                                  and gamma2 > 0.0 and persistence < 0.0)))
    runs, current, previous = [], [], None
    for entry, matches in flags:
        if matches and current and previous is not None and follows(previous, entry):
            current.append(entry)
        elif matches:
            current = [entry]
        else:
            if current:
                runs.append(current)
            current = []
        previous = entry
    if current:
        runs.append(current)
    blocks = [{"labels": [e["label"] for e in run], "start": run[0]["label"],
               "end": run[-1]["label"], "length": len(run)} for run in runs]
    longest = max(blocks, key=lambda block: block["length"]) if blocks else None
    return {"pattern": pattern, "runs": blocks, "longest_run": longest,
            "n_months": sum(1 for _e, matches in flags if matches),
            "months": [e["label"] for e, matches in flags if matches],
            "note": "descriptive only: no statistic is derived from these runs"}


def resampled_medians(values, draws, block, rng):
    """``draws`` circular moving-block bootstrap medians of an ordered series.

    Blocks of ``block`` consecutive chips are drawn with wrap-around, so the
    short-range correlation of a 6- to 12-day cadence survives the resampling --
    the same construction the orbit plate uses, which is why the interval means
    the same thing in both.
    """
    n = len(values)
    if n == 0:
        return []
    out = []
    for _ in range(draws):
        sample = []
        while len(sample) < n:
            start = rng.randrange(n)
            for step in range(block):
                sample.append(values[(start + step) % n])
        out.append(median(sample[:n]))
    return out


def bootstrap_cell(reference_values, other_values, rng):
    """``median(other) - median(reference)`` with R5's interval, in reference scales.

    ``other`` is the class being read -- the tail, or R2's outside class -- and
    ``reference`` the frame it is read in. Both sides are resampled from ``rng``,
    reference first and then other, so a caller that seeds once draws in exactly
    the order it did when this was a closure of ``contrast()`` and the contrast's
    committed numbers cannot move when the same interval is reused below.
    """
    observed = median(other_values) - median(reference_values)
    ref_draws = resampled_medians(reference_values, BOOT, BLOCK, rng)
    other_draws = resampled_medians(other_values, BOOT, BLOCK, rng)
    deltas = [o - r for o, r in zip(other_draws, ref_draws)]
    low, high = percentile(deltas, 2.5), percentile(deltas, 97.5)
    return {"delta": observed, "ci": [low, high],
            "excludes_zero": bool(low > 0.0 or high < 0.0),
            "direction": ("up" if observed > 0 else
                          "down" if observed < 0 else "flat")}


def contrast(site, chips):
    """Post minus reference per dimension and per block, in reference scales."""
    spec = SITES[site]
    if not spec["post"]:
        return {"applies": False, "draws": BOOT, "block": BLOCK, "seed": SEED,
                "reason": "the record labels no event: there is no post side to "
                          "contrast the reference with, and the plate draws the "
                          "trajectory alone"}
    ref = [c for c in chips if in_reference(site, c)]
    post = [c for c in chips if c["state"] in spec["post"]]
    rng = random.Random(SEED)
    result = {"applies": True, "draws": BOOT, "block": BLOCK, "seed": SEED,
              "states": list(spec["post"]), "n_reference": len(ref),
              "n_post": len(post), "dims": {}}

    # R5's interval, one resample per side (``bootstrap_cell``) -- the same
    # construction the outside class below is read with.
    for dim in DIMS:
        ref_values = [c["z"][dim] for c in ref if c["z"][dim] is not None]
        post_values = [c["z"][dim] for c in post if c["z"][dim] is not None]
        if not ref_values or not post_values:
            continue
        entry = bootstrap_cell(ref_values, post_values, rng)
        entry["n_reference"] = len(ref_values)
        entry["n_post"] = len(post_values)
        result["dims"][dim] = entry
    for name, dim in (("R_B", "gamma2"), ("R_G", "A")):
        ref_values = [c[name] for c in ref if c[name] is not None]
        post_values = [c[name] for c in post if c[name] is not None]
        if not ref_values or not post_values:
            continue
        entry = bootstrap_cell(ref_values, post_values, rng)
        entry["n_reference"] = len(ref_values)
        entry["n_post"] = len(post_values)
        entry["block"] = BLOCK_OF[dim]
        result["dims"][name] = entry
    return result


def outside_contrast(site, chips, months=None, envelope=None, stats=None):
    """The OUTSIDE class, read against the frozen reference.

    OUTSIDE is the third class of R2: reference-state chips the chosen window
    leaves out, on the pre-event side only (``outside_chips``). It standardises
    nothing -- the centre and the scale are the reference's -- so the plate reads
    it *in the reference's own frame*. This block reports, for every channel and
    for both block magnitudes, the signed median departure ``r`` (in transform
    units), the standardised departure ``z`` (``r`` over the reference's own
    robust scale), the seeded moving-block interval, the direction, and where the
    class sits in the **monthly** reference envelope.

    ``applies`` is false where a site has nothing outside -- R2's window is the
    state's own span there -- and the report then prints ``reason`` in place of
    the tables, so every site's report keeps the same sections.
    """
    spec = SITES[site]
    ref = [c for c in chips if in_reference(site, c)]
    outside = outside_chips(site, chips)
    months = monthly(site, chips) if months is None else months
    envelope = (reference_envelope(site, chips, months)
                if envelope is None else envelope)
    if stats is None:
        stats, _weather, _clipped = standardise(site, chips)
    held = [m for m in months if m["n"] and m["n_outside"] == m["n"]]
    if not outside:
        return {"applies": False, "draws": BOOT, "block": BLOCK, "seed": SEED,
                "n_reference": len(ref), "n_outside": 0, "n_intervals": 0,
                "intervals": [], "dims": {}, "persistence": {},
                "directional": {"runs": [], "longest_run": None},
                "reason": ("the reference window is the state's own span (R2), so "
                           "no healthy chip sits outside it and there is nothing "
                           "for the plate to test"
                           if spec["reference_window"] is None else
                           "the chosen reference window leaves no healthy chip "
                           "outside it, so there is nothing for the plate to test")}
    rng = random.Random(SEED)
    label = group_label(site, outside[0])
    result = {"applies": True, "draws": BOOT, "block": BLOCK, "seed": SEED,
              "class": "outside", "label": CLASS_LABEL["outside"],
              "state": spec["reference"], "window": window(outside),
              "months": len(held), "n_intervals": len(held),
              "intervals": [m["label"] for m in held],
              "n_reference": len(ref), "n_outside": len(outside), "dims": {},
              "criterion": ("a held-out month is an unusual departure only when "
                            "its monthly value exceeds the monthly reference q95 "
                            "of the same quantity (magnitudes, one-sided) or "
                            "leaves [q05, q95] (signed z); a run of at least %d "
                            "consecutive such months is sustained" % PERSIST_MIN)}

    # The bootstrap is the plate's own: one resample per side, reference first,
    # in the same order as the contrast, so the committed intervals cannot move.
    for dim in DIMS:
        ref_values = [c["z"][dim] for c in ref if c["z"][dim] is not None]
        out_values = [c["z"][dim] for c in outside if c["z"][dim] is not None]
        if not ref_values or not out_values:
            continue
        entry = bootstrap_cell(ref_values, out_values, rng)
        cell = stats[(dim, label)]
        out_transformed = [apply_transform(c["values"][dim], cell["transform"])
                           for c in outside if c["values"][dim] is not None]
        departure = median(out_transformed) - cell["centre"]
        entry.update({
            "n_reference": len(ref_values), "n_outside": len(out_values),
            "n_chips": len(out_values), "n_intervals": len(held),
            "intervals": [m["label"] for m in held],
            "r": departure,
            "r_unit": "transformed units (%s)" % cell["transform"],
            "z": departure / cell["scale"],
            "ci_r": [entry["ci"][0] * cell["scale"],
                     entry["ci"][1] * cell["scale"]],
            "reference": {"transform": cell["transform"],
                          "centre": cell["centre"], "scale": cell["scale"],
                          "source": cell["source"], "n": cell["n"],
                          "group": label},
            "reference_comparison": envelope_position(held, dim,
                                                      envelope["monthly"][dim]),
        })
        result["dims"][dim] = entry

    for name, dim in (("R_B", "gamma2"), ("R_G", "A")):
        ref_values = [c[name] for c in ref if c[name] is not None]
        out_values = [c[name] for c in outside if c[name] is not None]
        if not ref_values or not out_values:
            continue
        entry = bootstrap_cell(ref_values, out_values, rng)
        comparison = envelope_position(held, name, envelope["monthly"][name],
                                       one_sided=True)
        entry.update({
            "n_reference": len(ref_values), "n_outside": len(out_values),
            "n_chips": len(out_values), "n_intervals": len(held),
            "intervals": [m["label"] for m in held],
            "block": BLOCK_OF[dim], "value": median(out_values),
            "monthly_median": comparison["level"],
            "reference_comparison": comparison,
        })
        result["dims"][name] = entry

    result["persistence"] = {
        name: persistence_runs(held, name, envelope["monthly"][name])
        for name in ("R_B", "R_G")}
    result["directional"] = directional_runs(held)
    result["reference_envelope_months"] = envelope["n_reference_months"]
    return result


#: The audit's fourth bucket is not one of the three frozen classes, so it is
#: named here and not in ``CLASS_LABEL``: the registry's own windows leave no
#: chip before them, and a label no committed report carries would change all
#: six of them for nothing.
PREFIX_LABEL = "before the frozen reference"
PREFIX_NOTE = ("reference-state chips acquired before the frozen window: they "
               "standardise nothing and are held out from nothing")


COMPARE_ORDER = ("gamma2", "P", "R_B", "A", "D", "F", "S", "R_G")
COMPARE_LABEL = {"gamma2": "z(gamma2)", "P": "z(P)", "R_B": "R_B",
                 "A": "z(A)", "D": "z(D)", "F": "z(F)", "S": "z(S)",
                 "R_G": "R_G"}


def class_members(site, chips):
    """``{class: [chips]}`` for the three frozen classes (+ any unassigned)."""
    groups = {"reference": [], "outside": [], "prefix": [], "post": [],
              "unassigned": []}
    for chip in chips:
        groups[chip_class(site, chip)].append(chip)
    return groups


def class_comparison(site, chips, months, envelope):
    """The eight quantities, once per frozen class.

    Each cell is the median over the chips of that class, in the reference's own
    scales: the reference column standardises to zero by construction, the
    held-out column is the departure the assessment reads, and the post column is
    the tail the contrast already carries.
    """
    groups = class_members(site, chips)
    result = {"order": list(COMPARE_ORDER), "labels": dict(COMPARE_LABEL),
              "units": "z (reference scales); R is a per-chip RMS in the same "
                       "units", "classes": {}}
    for name in ("reference", "outside", "post"):
        group = groups[name]
        values = {}
        for key in COMPARE_ORDER:
            if key in DIMS:
                numbers = [c["z"][key] for c in group
                           if c["z"][key] is not None]
            else:
                numbers = [c[key] for c in group if c[key] is not None]
            values[key] = median(numbers) if numbers else None
        if name == "reference":
            bins = [m for m in months if m["n_reference"] == m["n"] and m["n"]]
        elif name == "outside":
            bins = [m for m in months if m["n_outside"] == m["n"] and m["n"]]
        else:
            bins = [m for m in months if m.get("n_post", 0) == m["n"] and m["n"]]
        result["classes"][name] = {"label": CLASS_LABEL[name],
                                    "note": CLASS_NOTE[name],
                                    "n_chips": len(group),
                                    "n_intervals": len(bins),
                                    "intervals": [m["label"] for m in bins],
                                    "values": values}
    return result


def summer_block(site, chips, months, envelope, outside):
    """The Summer-2024 subset of OUTSIDE -- post-hoc and descriptive only.

    The window is fixed in ``SUMMER_2024`` and is a strict subset of OUTSIDE. It
    is never used to define a threshold, never tested, and never an input to the
    verdict: it exists because a reader will ask what the last summer before the
    event looked like, and the answer must be labelled as a post-hoc description
    rather than passed off as a pre-registered contrast.
    """
    if not outside.get("applies"):
        return {"applies": False,
                "reason": "no held-out pre-event class on this site"}
    start, end = (dt.date.fromisoformat(text) for text in SUMMER_2024)
    first, last = SUMMER_2024[0][:7], SUMMER_2024[1][:7]
    held = [c for c in chips
            if chip_class(site, c) == "outside" and start <= c["day"] <= end]
    bins = [m for m in months if m["n"] and m["n_outside"] == m["n"]
            and first <= m["month"] <= last]

    def medians(group):
        out = {}
        for key in COMPARE_ORDER:
            if key in DIMS:
                numbers = [c["z"][key] for c in group
                           if c["z"][key] is not None]
            else:
                numbers = [c[key] for c in group if c[key] is not None]
            out[key] = median(numbers) if numbers else None
        return out

    months_pattern = []
    for month in bins:
        gamma2, persistence = month["z"]["gamma2"], month["z"]["P"]
        months_pattern.append({
            "month": month["label"], "gamma2": gamma2, "P": persistence,
            "R_B": month["R_B"], "R_G": month["R_G"],
            "pattern": bool(gamma2 is not None and persistence is not None
                            and gamma2 > 0.0 and persistence < 0.0)})
    above = {name: [label for label in
                    outside["persistence"][name]["months_above_q95"]
                    if first <= label[:7] <= last]
             for name in ("R_B", "R_G")}
    return {"applies": True, "descriptive_only": True, "post_hoc": True,
            "inferential": False, "window": list(SUMMER_2024),
            "n_chips": len(held), "n_intervals": len(bins),
            "intervals": [m["label"] for m in bins],
            "medians": medians(held), "months": months_pattern,
            "months_above_q95": above,
            "note": "post-hoc descriptive subset of OUTSIDE; no interval, no "
                    "p-value, no new statistic and no verdict input"}


def environment_context(site, months, outside):
    """Descriptive coincidence of held-out departures with the co-variates.

    The environment is reported *beside* the OCV numbers, never inside them: no
    weather correction, no attribution, no causal claim. An exceedance that
    coincides with a weather departure is still an OCV statement.
    """
    if not outside.get("applies"):
        return {"applies": False}
    held = [m for m in months if m["n"] and m["n_outside"] == m["n"]]
    rows = []
    for month in held:
        rows.append({
            "month": month["label"], "R_B": month["R_B"],
            "R_G": month["R_G"], "gamma2": month["z"]["gamma2"],
            "P": month["z"]["P"], "weather_z": dict(month["weather_z"]),
            "above_q95": {name: bool(month[name] is not None
                                      and month[name] >
                                      outside["persistence"][name]["threshold"])
                          for name in ("R_B", "R_G")}})
    coincident = [row["month"] for row in rows
                  if any(row["above_q95"].values())]
    return {"applies": True, "variables": list(SITES[site]["weather"]),
            "months": rows, "coincident_months": coincident,
            "note": "The co-variates are reported beside the OCV quantities, "
                    "not inside them: no weather correction is applied and no "
                    "causal attribution is made"}


ALLOWED_VERDICTS = (
    "no unusual departure resolved",
    "isolated unusual departure",
    "sustained unusual behaviour departure",
    "sustained unusual geometry departure",
    "sustained multidimensional observability departure",
)


def verdict_word(behaviour, geometry):
    """The verdict, derived from the two block-magnitude persistence rules."""
    if behaviour["sustained"] and geometry["sustained"]:
        return ALLOWED_VERDICTS[4]
    if behaviour["sustained"]:
        return ALLOWED_VERDICTS[2]
    if geometry["sustained"]:
        return ALLOWED_VERDICTS[3]
    if behaviour["n_above_q95"] or geometry["n_above_q95"]:
        return ALLOWED_VERDICTS[1]
    return ALLOWED_VERDICTS[0]

def assessment(site, chips, months, envelope, outside, comparison, summer,
               environment):
    """The six questions, each answered from the frozen classes and nothing else.

    The verdict is derived, never chosen: it is a function of the two block
    magnitudes' persistence rules, and nothing else -- not the summer subset, not
    the signed channels, not the environment. The wording comes from
    ``ALLOWED_VERDICTS``.
    """
    if not outside.get("applies"):
        return {"applies": False,
                "reason": "there is no held-out pre-event class on this site: "
                          "the reference window is the reference state's own "
                          "span, so there is nothing to assess against it"}
    behaviour = outside["persistence"]["R_B"]
    geometry = outside["persistence"]["R_G"]
    verdict = verdict_word(behaviour, geometry)
    basis = {name: {
        "threshold": outside["persistence"][name]["threshold"],
        "threshold_quantile": DECISION_QUANTILE,
        "n_above_q95": outside["persistence"][name]["n_above_q95"],
        "months_above_q95": outside["persistence"][name]["months_above_q95"],
        "longest_run": outside["persistence"][name]["longest_run"],
        "sustained": outside["persistence"][name]["sustained"]}
        for name in ("R_B", "R_G")}

    def span(run):
        if not run:
            return "no run"
        return "%s to %s" % (run["start"], run["end"])

    def labels(items):
        return ", ".join(items) if items else "none"

    gamma2 = outside["dims"].get("gamma2", {}).get("reference_comparison", {})
    persistence_dim = outside["dims"].get("P", {}).get("reference_comparison", {})
    pattern_run = outside["directional"].get("longest_run")
    directional = bool(gamma2.get("outside_envelope")
                       or persistence_dim.get("outside_envelope")
                       or (pattern_run and pattern_run["length"] >= PERSIST_MIN))
    pieces = []
    if persistence_dim.get("months_below_q05"):
        pieces.append("z(P) lies below its monthly reference envelope in %d of %d "
                      "held-out months (%s)"
                      % (persistence_dim["n_below_q05"],
                         persistence_dim["n_months"],
                         labels(persistence_dim["months_below_q05"])))
    if persistence_dim.get("months_above_q95"):
        pieces.append("z(P) lies above it in %s"
                      % labels(persistence_dim["months_above_q95"]))
    if gamma2.get("months_above_q95"):
        pieces.append("z(gamma2) lies above it in %s"
                      % labels(gamma2["months_above_q95"]))
    if gamma2.get("months_below_q05"):
        pieces.append("z(gamma2) lies below it in %s"
                      % labels(gamma2["months_below_q05"]))
    if pattern_run and pattern_run["length"] >= PERSIST_MIN:
        pieces.append("the joint pattern z(gamma2) > 0 and z(P) < 0 holds for %d "
                      "consecutive months (%s)"
                      % (pattern_run["length"], span(pattern_run)))
    if directional:
        answer = ("Yes - directional reconfiguration is visible: " + "; ".join(pieces)
                  + ". This is a statement about direction, not about the size of "
                    "the displacement, and it is not an anomaly verdict by "
                    "itself.")
    else:
        answer = ("No - the signed monthly medians stay inside the monthly "
                  "reference envelope and the joint pattern does not persist.")
    q1 = {"question": "Is there directional behaviour reconfiguration in the "
                       "held-out period?",
          "answer": answer,
          "evidence": {"gamma2_months_above_q95": gamma2.get("months_above_q95", []),
                       "gamma2_months_below_q05": gamma2.get("months_below_q05", []),
                       "P_months_above_q95":
                           persistence_dim.get("months_above_q95", []),
                       "P_months_below_q05":
                           persistence_dim.get("months_below_q05", []),
                       "months_outside_envelope": sorted(
                           set(gamma2.get("outside_envelope", []))
                           | set(persistence_dim.get("outside_envelope", []))),
                       "pattern_runs": outside["directional"].get("runs", []),
                       "longest_pattern_run": pattern_run}}

    def magnitude_question(question, key):
        cell = outside["dims"].get(key, {}).get("reference_comparison", {})
        runs = outside["persistence"].get(key, {})
        longest = runs.get("longest_run")
        if cell.get("months_above_q95"):
            answer = ("Yes - %d of %d held-out months %s above the monthly "
                      "reference q95 of `%s` (%s); the longest consecutive run is "
                      "%d month(s) (%s), so the departure is %s."
                      % (runs.get("n_above_q95", 0), runs.get("n_months", 0),
                         "lies" if runs.get("n_above_q95", 0) == 1 else "lie",
                         key, labels(cell["months_above_q95"]),
                         longest["length"] if longest else 0, span(longest),
                         "sustained under the pre-specified rule"
                         if runs.get("sustained") else
                         "not sustained under the pre-specified rule"))
        else:
            answer = ("No - no held-out month exceeds the monthly reference q95 "
                      "of `%s`: the highest monthly value is %s against a "
                      "threshold of %s."
                      % (key, number(cell.get("max_value")),
                         number(cell.get("q95"))))
        return {"question": question, "answer": answer,
                "evidence": {"threshold": cell.get("q95"),
                             "months_above_q95": cell.get("months_above_q95", []),
                             "n_above_q95": cell.get("n_above_q95", 0),
                             "longest_run": runs.get("longest_run"),
                             "sustained": runs.get("sustained", False)}}

    q2 = magnitude_question("Does the held-out behaviour magnitude (R_B) exceed "
                            "the monthly reference envelope?", "R_B")
    q3 = magnitude_question("Does the held-out geometry magnitude (R_G) exceed "
                            "the monthly reference envelope?", "R_G")
    q4 = {"question": "Is any held-out departure temporally persistent?",
          "answer": ("`R_B`: %d of %d months above the q95 threshold, longest run "
                     "%d month(s) (%s); `R_G`: %d of %d months above the q95 "
                     "threshold, longest run %d month(s) (%s). Persistence is "
                     "decided by the fixed rule (at least %d consecutive months), "
                     "not by the shape of the curve."
                     % (behaviour["n_above_q95"], behaviour["n_months"],
                        (behaviour["longest_run"] or {"length": 0})["length"],
                        span(behaviour["longest_run"]),
                        geometry["n_above_q95"], geometry["n_months"],
                        (geometry["longest_run"] or {"length": 0})["length"],
                        span(geometry["longest_run"]), PERSIST_MIN)),
          "evidence": {"R_B": basis["R_B"], "R_G": basis["R_G"]}}

    groups = class_members(site, chips)
    start, end = (dt.date.fromisoformat(text) for text in SUMMER_2024)
    summer_group = [c for c in groups["outside"] if start <= c["day"] <= end]
    rest = [c for c in groups["outside"] if not (start <= c["day"] <= end)]

    def medians(group):
        out = {}
        for key in COMPARE_ORDER:
            if key in DIMS:
                numbers = [c["z"][key] for c in group
                           if c["z"][key] is not None]
            else:
                numbers = [c[key] for c in group if c[key] is not None]
            out[key] = median(numbers) if numbers else None
        return out

    summer_medians, rest_medians = medians(summer_group), medians(rest)
    q5 = {"question": "Does the last summer before the event differ from the "
                       "rest of the held-out period?",
          "answer": ("Post-hoc, descriptive only: the Summer-2024 subset (%d chips, "
                     "%d monthly intervals, %s to %s) has median z(gamma2) = %s, "
                     "z(P) = %s, R_B = %s, R_G = %s; the remaining held-out chips "
                     "(%d) have %s, %s, %s, %s. The difference is not tested, no "
                     "interval is attached to it, and it is not an input to the "
                     "verdict."
                     % (len(summer_group), len(summer["intervals"]),
                        SUMMER_2024[0], SUMMER_2024[1],
                        signed(summer_medians["gamma2"]),
                        signed(summer_medians["P"]),
                        number(summer_medians["R_B"]),
                        number(summer_medians["R_G"]), len(rest),
                        signed(rest_medians["gamma2"]), signed(rest_medians["P"]),
                        number(rest_medians["R_B"]), number(rest_medians["R_G"]))),
          "evidence": {"window": list(SUMMER_2024), "descriptive_only": True,
                       "post_hoc": True, "n_chips": len(summer_group),
                       "n_intervals": len(summer["intervals"]),
                       "n_remainder": len(rest), "summer": summer_medians,
                       "remainder": rest_medians, "difference_tested": False}}

    coincident = environment.get("coincident_months", [])
    if coincident:
        rows = [row for row in environment.get("months", [])
                if row["month"] in coincident]
        detail = "; ".join(
            "%s: %s" % (row["month"], ", ".join(
                "%s z = %s" % (column, signed(value))
                for column, value in row["weather_z"].items()
                if value is not None))
            for row in rows)
        answer = ("The month(s) above the magnitude envelope - %s - coincide with "
                  "environmental variation (%s). This is a temporal coincidence "
                  "only: no weather correction is applied and no causal "
                  "attribution is made." % (", ".join(coincident), detail))
    else:
        answer = ("No held-out month exceeds the magnitude envelope, so there is "
                  "no held-out departure to compare with the environment.")
    q6 = {"question": "Do the pronounced held-out departures coincide with "
                       "environmental variation?",
          "answer": answer,
          "evidence": {"variables": environment.get("variables", []),
                       "coincident_months": coincident,
                       "months": environment.get("months", []),
                       "note": environment.get("note")}}

    return {"applies": True, "criterion": outside["criterion"],
            "questions": {
                "directional_behaviour_reconfiguration": q1,
                "behaviour_magnitude": q2,
                "geometry_magnitude": q3,
                "temporal_persistence": q4,
                "summer_2024": q5,
                "environmental_coincidence": q6},
            "verdict": verdict, "allowed_verdicts": list(ALLOWED_VERDICTS),
            "verdict_inputs": ["R_B", "R_G"], "basis": basis,
            "class_medians": comparison["classes"],
            "summer_is_descriptive_only": True,
            "not_claiming": ["damage", "collapse prediction",
                             "failure precursor", "causation",
                             "a new significance statistic",
                             "a directional score"]}


PANELS = (
    ("a1", "behaviour", "z(gamma2) -- monthly median, with the monthly reference envelope"),
    ("a2", "behaviour", "z(P) -- monthly median, with the monthly reference envelope"),
    ("a3", "behaviour", "the behaviour plane (z(gamma2), z(P))"),
    ("a4", "behaviour", "R_B -- the behaviour block magnitude, with its monthly "
                        "reference envelope"),
    ("b1", "geometry", "z(A) -- monthly median, with the monthly reference envelope"),
    ("b2", "geometry", "z(F) -- monthly median, with the monthly reference envelope"),
    ("b3", "geometry", "the geometry plane (z(A), z(F))"),
    ("b4", "geometry", "R_G -- the geometry block magnitude, with its monthly "
                       "reference envelope"),
    ("c1", "context", "environmental context (weather, descriptive only)"),
    ("c2", "context", "the observation cadence per orbit"),
    ("c3", "context", "observation availability: solid = an echo, hollow = a "
                      "mask without one"),
    ("c4", "context", "the trajectory matrix: six dimensions x every month"),
)


def window(chips):
    """``[first, last]`` day of a chip list, as ISO strings."""
    if not chips:
        return None
    return [chips[0]["day"].isoformat(), chips[-1]["day"].isoformat()]


def build(site, chips=None):
    """The whole report of one site -- every number the plate draws.

    ``chips`` substitutes the support cut: the prospective runner re-derives a
    report with the post-event chips deleted (a hindsight guard) and compares
    the frozen reference and the held-out class with the published ones. Nothing
    that is published is ever built this way.
    """
    spec = SITES[site]
    if chips is None:
        chips, audit = load_chips(site)
    else:
        chips = list(chips)
        audit = {"rows": len(chips), "kept": len(chips), "excluded": {},
                 "refused": [], "substituted": True}
    if not chips:
        raise SystemExit("%s: the support table holds no chip" % site)
    stats, weather, clipped = standardise(site, chips)
    months = monthly(site, chips)
    envelope = reference_envelope(site, chips, months)
    con = contrast(site, chips)
    out_con = outside_contrast(site, chips, months, envelope, stats)
    comparison = class_comparison(site, chips, months, envelope)
    summer = summer_block(site, chips, months, envelope, out_con)
    environment = environment_context(site, months, out_con)
    verdict = assessment(site, chips, months, envelope, out_con, comparison,
                         summer, environment)
    ref_chips, outside_class, state_chips = reference_split(site, chips)
    members = class_members(site, chips)
    post_chips = members["post"]
    unassigned = members["unassigned"]
    labels = sorted({label for _dim, label in stats})

    def census(group):
        return {"window": window(group), "n_chips": len(group),
                "n_ascending": sum(1 for c in group
                                   if c["orbit"] == "ASCENDING"),
                "n_descending": sum(1 for c in group
                                    if c["orbit"] == "DESCENDING"),
                "n_echo": sum(1 for c in group if c["echo"]),
                "months": len({month_key(c["day"]) for c in group})}

    report = {
        "schema": SCHEMA,
        "site": site,
        "label": spec["label"],
        "note": spec["note"],
        "event": spec["event"],
        "table": os.path.relpath(table_path(site), REPO),
        "sha256": sha256_of(table_path(site)),
        "cut": {"rows": audit["rows"], "kept": audit["kept"],
                "excluded": audit["excluded"], "refused": audit["refused"],
                "window": window(chips)},
        "reference": dict(
            census(ref_chips),
            state=reference_label(site),
            chosen=list(reference_window(site)) if reference_window(site)
            else None,
            policy=("per orbit" if site == "espoo" else "pooled over the orbits"),
            groups=labels,
            # The whole reference state, and what the choice leaves out of it --
            # both counts are censuses, never statistics (R2).
            span=census(state_chips),
            outside=census(outside_class) if outside_class else {
                "window": None, "n_chips": 0, "n_ascending": 0,
                "n_descending": 0, "n_echo": 0, "months": 0},
            choice=(("the window the stated rule selected for this asset (%s)"
                     % window_rule(site))
                    if window_rule(site) else
                    "the registry's window for this asset"
                    if reference_window(site) else
                    "no sub-window: the state's own span"),
            label=CLASS_LABEL["reference"],
            note=CLASS_NOTE["reference"],
        ),
        "post": dict(census(post_chips), states=list(spec["post"]),
                     label=CLASS_LABEL["post"], note=CLASS_NOTE["post"])
        if post_chips else None,
        "classes": {
            "labels": dict(CLASS_LABEL), "notes": dict(CLASS_NOTE),
            "event": spec["event"],
            "event_rule": ("OUTSIDE is strictly before %s; POST is the "
                           "registry's post state(s); a chip on the event day "
                           "is assigned by its registry state, and a "
                           "reference-state chip on or after the event is "
                           "reported as unassigned"
                           % (spec["event"] or "the (absent) event")),
            "reference": dict(census(ref_chips),
                              label=CLASS_LABEL["reference"]),
            "outside": dict(census(outside_class),
                            label=CLASS_LABEL["outside"])
            if outside_class else {
                "label": CLASS_LABEL["outside"], "window": None,
                "n_chips": 0, "n_ascending": 0, "n_descending": 0,
                "n_echo": 0, "months": 0},
            "post": dict(census(post_chips), label=CLASS_LABEL["post"])
            if post_chips else {
                "label": CLASS_LABEL["post"], "window": None,
                "n_chips": 0, "n_ascending": 0, "n_descending": 0,
                "n_echo": 0, "months": 0},
            "unassigned": [{"day": c["day"].isoformat(),
                            "orbit": c["orbit"], "state": c["state"]}
                           for c in unassigned],
            "n_unassigned": len(unassigned),
            "reconciliation": {"kept": len(chips),
                               "reference": len(ref_chips),
                               "outside": len(outside_class),
                               "post": len(post_chips),
                               "unassigned": len(unassigned)},
        },
        "reference_envelope": envelope,
        "class_comparison": comparison,
        "summer_2024": summer,
        "environment_context": environment,
        "assessment": verdict,
        "standardisation": {
            label: {dim: dict(stats[(dim, label)]) for dim in DIMS}
            for label in labels
        },
        "weather_reference": {
            column: {label: dict(cell) for label, cell in per_group.items()}
            for column, per_group in weather.items()
        },
        "clipped": clipped,
        "bootstrap": {"draws": BOOT, "block": BLOCK, "seed": SEED,
                      "kind": "circular moving-block, one resample per group"},
        "months": months,
        "contrast": con,
        "outside": out_con,
        "figure": {"dpi": FIGURE_DPI, "panels": [
            {"id": pid, "block": block, "shows": shows}
            for pid, block, shows in PANELS]},
    }
    # The prospective additions are emitted only when the build actually has
    # them: the registry's own windows leave no chip before them, so every
    # committed report, and therefore every committed plate, stays identical.
    if members["prefix"]:
        report["classes"]["prefix"] = dict(census(members["prefix"]),
                                           label=PREFIX_LABEL)
        report["classes"]["reconciliation"]["prefix"] = len(members["prefix"])
    if window_rule(site):
        report["reference"]["choice_rule"] = window_rule(site)
    return report


COLOUR = {"behaviour": "#2b6cb0", "geometry": "#2f855a", "context": "#6b46c1",
          "asc": "#2b6cb0", "desc": "#dd6b20", "reference": "#4a5568",
          "post": "#c53030", "mixed": "#b7791f", "band": "#e2e8f0",
          "outside": "#a0aec0", "envelope": "#cbd5e0", "prefix": "#b7c0c9"}

#: The marker of each month class. ``outside`` (R2's third class) is a hollow
#: grey circle -- the same shape as a reference month, a different colour and an
#: empty face, so a reader never mistakes a chip the choice left out for one the
#: standardisation rests on.
MARKER = {"ref": "o", "outside": "o", "prefix": "o", "post": "s",
          "mixed": "D"}


def class_colour(cls, colour):
    """The colour a month class is drawn in: its block's, or grey if outside.

    A month before the frozen window is the holder grey of the held-out class in
    a lighter tone: same shape, same empty face -- it belongs to the audit, not
    to the standardisation.
    """
    if cls == "prefix":
        return COLOUR["prefix"]
    return COLOUR["outside"] if cls == "outside" else colour

WEATHER_LABEL = {"temperature_c": "T [C]", "wind_speed_ms": "wind [m/s]",
                 "precipitation_mm": "precipitation [mm/h]",
                 "relative_humidity_pct": "relative humidity [%]"}


def month_class(entry):
    """``ref`` / ``prefix`` / ``outside`` / ``post`` / ``mixed`` -- its place.

    ``outside`` is R2's third class: a month of the reference *state* that the
    chosen reference window leaves out. It is drawn in the reference frame (grey,
    hollow) but it standardises nothing, so it is neither ``ref`` nor ``post``.
    ``prefix`` is the audit's fourth bucket: a reference-state month entirely
    before the frozen window. Precedence is ``ref`` > ``prefix`` > ``post`` >
    ``outside`` > ``mixed``; a month can only be ``mixed`` if the window rule is
    ever relaxed to day granularity.
    """
    if entry["n_reference"] == entry["n"]:
        return "ref"
    if entry.get("n_prefix", 0) == entry["n"]:
        return "prefix"
    if entry.get("n_outside", 0) == entry["n"]:
        return "outside"
    if entry.get("n_post", 0) == entry["n"]:
        return "post"
    return "mixed"


def packed(fig, parts, fontsize, width_in, sep=" · "):
    """``parts`` joined by ``sep``, broken into the fewest lines that fit ``width_in``.

    A header written as one line is wider than the plate and, because the plate is
    saved with ``bbox_inches="tight"``, the *file* is then widened to hold the header
    -- the plate stops being the size the figure says it is (the Carola plate came out
    16.7 in wide for a 7.2 in column). The lines are measured on the figure's own
    renderer at ``fontsize``, so the break falls where the text actually ends instead
    of where a character count guesses it does.

    The break is not greedy. The parts are distributed over ``k`` lines, ``k`` as small
    as fits, so that the longest line is as short as it can be; between equally long
    choices the last line is given as much as it can take, because a header whose last
    line is a lonely date reads as an accident rather than as a design.
    """
    try:
        renderer = fig.canvas.get_renderer()
    except AttributeError:                    # a backend without an Agg renderer
        renderer = None
    probe = fig.text(0.0, 0.0, "", fontsize=fontsize)
    cache = {}

    def width(first, last):
        """The width in inches of ``parts[first:last]`` joined on one line."""
        if (first, last) not in cache:
            text = sep.join(parts[first:last])
            if renderer is None:
                cache[(first, last)] = 0.5 * fontsize / 72.0 * len(text)
            else:
                probe.set_text(text)
                cache[(first, last)] = (probe.get_window_extent(renderer).width
                                        / fig.dpi)
        return cache[(first, last)]

    def splits(lines):
        """For every prefix of ``parts``, the best way to spread it over ``lines``."""
        table = [{} for _ in range(len(parts) + 1)]
        table[0][0] = ((0.0, 0.0), 0)
        for end in range(1, len(parts) + 1):
            for count in range(1, lines + 1):
                for cut in range(count - 1, end):
                    if count - 1 not in table[cut]:
                        continue
                    here = width(cut, end)
                    worst = max(table[cut][count - 1][0][0], here)
                    trial = ((worst, -here), cut)
                    if count not in table[end] or trial < table[end][count]:
                        table[end][count] = trial
        return table

    lines = len(parts)
    for wanted in range(1, len(parts) + 1):
        table = splits(wanted)
        if wanted in table[len(parts)] and table[len(parts)][wanted][0][0] <= width_in:
            lines = wanted
            break
    out, end = [], len(parts)
    while lines:
        cut = table[end][lines][1]
        out.append(sep.join(parts[cut:end]))
        end, lines = cut, lines - 1
    probe.remove()
    return chr(10).join(reversed(out))


def fit_legend(fig, handles, width_in, centre, gap_in=0.07):
    """``fig.legend`` hung under the panels, in as few rows as fit ``width_in``.

    The same trap as the header: one row of long labels is wider than the plate, and
    ``bbox_inches="tight"`` then widens the *file* to hold the key instead of the key
    being re-flowed to the plate. Columns are dropped one at a time -- widest row
    first -- until the legend, measured by the figure's own renderer, fits, and the rows
    are centred on ``centre`` (the panel block's own centre, not the canvas'). The key
    hangs from the bottom of the panels' own ink, tick labels included, so a second row
    grows downwards into the tight box instead of upwards into the month labels.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    lowest = min(ax.get_tightbbox(renderer).y0 for ax in fig.axes)
    anchor = (centre, (lowest - gap_in * fig.dpi) / fig.bbox.height)
    legend = None
    for ncol in range(len(handles), 0, -1):
        if legend is not None:
            legend.remove()
        legend = fig.legend(handles=handles, loc="upper center", ncol=ncol,
                            frameon=False, bbox_to_anchor=anchor)
        fig.canvas.draw()
        if (legend.get_window_extent(fig.canvas.get_renderer()).width / fig.dpi
                <= width_in):
            break
    return legend


def draw(report, path, caption=None):
    """The 3 x 4 plate. Every ordinate drawn here is a number of ``report``.

    ``caption`` replaces the first header line, which is all the prospective run
    changes: its report is drawn from a rule-selected window, and the header
    says so rather than calling it a choice. It is not called ``title`` because
    ``draw`` has a nested ``title`` of its own and the parameter would be
    shadowed by it.

    One design for all six sites: the frozen reference months filled, the
    held-out pre-event months hollow grey circles, the post-event months hollow
    coloured squares, the event a red dashed line, the frozen reference window a
    grey band, and the monthly reference envelope (q05-q95 of the reference
    months' monthly medians) a light band behind each magnitude panel.
    """
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "mathtext.fontset": "stix",
        "font.size": 7.4,
        "axes.titlesize": 7.4,
        "axes.labelsize": 6.9,
        "xtick.labelsize": 6.2,
        "ytick.labelsize": 6.2,
        "legend.fontsize": 6.2,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.22,
        "grid.linewidth": 0.4,
        "figure.dpi": 110,
        "savefig.dpi": FIGURE_DPI,
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,
    })
    months = report["months"]
    keys = [bin_key(month) for month in months]
    xs = [bin_date(month) for month in months]
    classes = [month_class(month) for month in months]
    envelopes = report["reference_envelope"]["monthly"]
    event = dt.date.fromisoformat(report["event"]) if report["event"] else None
    ref_window = report["reference"]["window"]
    fig = plt.figure(figsize=(COL_WIDTH, 0.92 * COL_WIDTH))
    grid = fig.add_gridspec(3, 4, hspace=0.80, wspace=0.62,
                            left=0.068, right=0.975, top=0.875, bottom=0.145)
    where = {}
    for row in range(3):
        for col in range(4):
            where["abc"[row] + str(col + 1)] = fig.add_subplot(grid[row, col])

    def linked(index):
        """True when bin i and bin i + 1 are adjacent in time and same side."""
        if classes[index] == "mixed" or classes[index + 1] == "mixed":
            return False
        left, right = keys[index], keys[index + 1]
        if left[:2] == right[:2]:
            return False                      # the event splits this month
        year, month = (left[0] + 1, 1) if left[1] == 12 else (left[0], left[1] + 1)
        return right[0] == year and right[1] == month

    def decorate(ax, zero=True, envelope_key=None):
        cell = envelopes.get(envelope_key) if envelope_key else None
        if cell:
            ax.axhspan(cell["q05"], cell["q95"], color=COLOUR["envelope"],
                       alpha=0.55, zorder=0, lw=0)
        span = report["reference"]["span"]["window"]
        if span and report["reference"]["outside"]["n_chips"]:
            ax.axvspan(dt.date.fromisoformat(span[0]),
                       dt.date.fromisoformat(span[1]),
                       color=COLOUR["outside"], alpha=0.30, zorder=0, lw=0)
        if ref_window:
            ax.axvspan(dt.date.fromisoformat(ref_window[0]),
                       dt.date.fromisoformat(ref_window[1]),
                       color=COLOUR["band"], alpha=0.75, zorder=0, lw=0)
        if event:
            ax.axvline(event, color=COLOUR["post"], lw=0.8, ls=(0, (4, 2)),
                       zorder=1)
        if zero:
            ax.axhline(0.0, color="#4a5568", lw=0.6, ls=":", zorder=1)
        ax.xaxis.set_major_locator(mdates.AutoDateLocator(minticks=3,
                                                          maxticks=6))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
        for label in ax.get_xticklabels():
            label.set_rotation(30)
            label.set_horizontalalignment("right")

    def title(ax, ident, text):
        ax.set_title("%s  ·  %s" % (ident, text), loc="left", pad=3.0,
                     fontsize=7.2)

    def series(ax, values, ident, text, colour, ylabel, envelope_key=None):
        """One monthly series: solid while reference, dashed once post."""
        decorate(ax, envelope_key=envelope_key)
        for i in range(len(xs) - 1):
            if not linked(i):
                continue
            style = ("-" if classes[i] in ("ref", "outside", "prefix")
                     else (0, (3, 2)))
            ax.plot(xs[i:i + 2], values[i:i + 2],
                    color=class_colour(classes[i], colour), lw=1.0,
                    ls=style, zorder=2)
        for i, (x, value) in enumerate(zip(xs, values)):
            if value is None:
                continue
            ax.plot([x], [value], marker=MARKER[classes[i]], ms=2.6,
                    mfc=colour if classes[i] == "ref" else "white",
                    mec=class_colour(classes[i], colour), mew=0.7, zorder=3)
        ax.set_ylabel(ylabel)
        title(ax, ident, text)
        return ax

    def plane(ax, ident, px, py, cx, cy, text):
        """The two dimensions of a block, month by month -- no date axis here."""
        ax.axhline(0.0, color="#4a5568", lw=0.5, ls=":")
        ax.axvline(0.0, color="#4a5568", lw=0.5, ls=":")
        for i in range(len(px) - 1):
            if linked(i):
                ax.plot(px[i:i + 2], py[i:i + 2], color="#a0aec0", lw=0.6,
                        zorder=1)
        for i in range(len(px)):
            ax.plot([px[i]], [py[i]], marker=MARKER[classes[i]], ms=2.8,
                    mfc=COLOUR["reference"] if classes[i] == "ref" else "white",
                    mec=(COLOUR["post"] if classes[i] == "post" else
                         class_colour(classes[i], COLOUR["reference"])),
                    mew=0.7, zorder=3)
        ax.set_xlabel(cx)
        ax.set_ylabel(cy)
        title(ax, ident, text)
        return ax

    def line_only(ax, values, colour):
        for i in range(len(xs) - 1):
            if not linked(i):
                continue
            ax.plot(xs[i:i + 2], values[i:i + 2],
                    color=class_colour(classes[i], colour), lw=1.0,
                    ls="-" if classes[i] in ("ref", "outside") else (0, (3, 2)),
                    zorder=2)
        for i, value in enumerate(values):
            if value is None:
                continue
            ax.plot([xs[i]], [value],
                    marker=MARKER[classes[i]],
                    ms=2.4, mfc=colour if classes[i] == "ref" else "white",
                    mec=class_colour(classes[i], colour), mew=0.6, zorder=3)

    dim_series = {dim: [m["z"][dim] for m in months] for dim in DIMS}
    block_series = {name: [m[name] for m in months] for name in ("R_B", "R_G")}
    ref = report["reference"]
    post = report["post"]

    # --- row a: the behaviour block ---------------------------------------
    series(where["a1"], dim_series["gamma2"], "a1",
           "z(%s)" % DIM_LABEL["gamma2"], COLOUR["behaviour"],
           "z(%s)" % DIM_LABEL["gamma2"], "gamma2")
    series(where["a2"], dim_series["P"], "a2", "z(P)", COLOUR["behaviour"],
           "z(P)", "P")
    plane(where["a3"], "a3", dim_series["gamma2"], dim_series["P"],
          "z(%s)" % DIM_LABEL["gamma2"], "z(P)", "behaviour plane")
    series(where["a4"], block_series["R_B"], "a4", "R_B", COLOUR["behaviour"],
           "R_B", "R_B")

    # --- row b: the geometry block ----------------------------------------
    series(where["b1"], dim_series["A"], "b1", "z(A)", COLOUR["geometry"],
           "z(A)", "A")
    series(where["b2"], dim_series["F"], "b2", "z(F)", COLOUR["geometry"],
           "z(F)", "F")
    plane(where["b3"], "b3", dim_series["A"], dim_series["F"],
          "z(A)", "z(F)", "geometry plane")
    series(where["b4"], block_series["R_G"], "b4", "R_G", COLOUR["geometry"],
           "R_G", "R_G")

    # --- row c: the record's own context ----------------------------------
    ax = where["c1"]
    columns = list(months[0]["weather_z"]) if months else []
    decorate(ax, zero=True)
    for index, column in enumerate(columns[:2]):
        line_only(ax, [m["weather_z"].get(column) for m in months],
                  [COLOUR["context"], COLOUR["asc"]][index])
    ax.set_ylabel("z (own reference)")
    for index, column in enumerate(columns[:2]):
        ax.text(0.03, 0.96 - 0.11 * index, WEATHER_LABEL.get(column, column),
                transform=ax.transAxes, fontsize=5.4, va="top",
                color=[COLOUR["context"], COLOUR["asc"]][index],
                bbox=dict(facecolor="white", alpha=0.7, lw=0, pad=1.0))
    title(ax, "c1", "environmental context")

    ax = where["c2"]
    decorate(ax, zero=False)
    ax.bar(xs, [m["n_ascending"] for m in months], width=20,
           color=COLOUR["asc"], lw=0)
    ax.bar(xs, [m["n_descending"] for m in months], width=20,
           bottom=[m["n_ascending"] for m in months], color=COLOUR["desc"], lw=0)
    ax.set_ylabel("chips per month")
    tallest = max([m["n"] for m in months] or [1])
    ax.set_ylim(0, tallest * 1.32)
    ax.legend(handles=[Patch(facecolor=COLOUR["asc"], label="ASC"),
                       Patch(facecolor=COLOUR["desc"], label="DESC")],
              loc="upper left", frameon=False, handlelength=1.0,
              borderaxespad=0.4, labelspacing=0.35, fontsize=5.2)
    title(ax, "c2", "cadence")

    ax = where["c3"]
    decorate(ax, zero=False)
    echo = [m["n_echo"] for m in months]
    plain = [m["n"] - m["n_echo"] for m in months]
    ax.bar(xs, echo, width=20, color=COLOUR["behaviour"], lw=0)
    ax.bar(xs, plain, width=20, bottom=echo, facecolor="white",
           edgecolor=COLOUR["behaviour"], lw=0.6)
    ax.set_ylabel("chips per month")
    ax.set_ylim(0, max([m["n"] for m in months] or [1]) * 1.12)
    title(ax, "c3", "observation availability")

    ax = where["c4"]
    grid_values = [[(m["z"][dim] if m["z"][dim] is not None else float("nan"))
                    for m in months] for dim in DIMS]
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#f7fafc")
    image = ax.imshow(grid_values, aspect="auto", cmap=cmap, vmin=-3.0,
                      vmax=3.0, interpolation="nearest")
    ax.set_yticks(range(len(DIMS)))
    ax.set_yticklabels([DIM_LABEL[dim] for dim in DIMS])
    step = max(1, int(math.ceil(len(months) / 5.0)))
    ticks = list(range(0, len(months), step))
    ax.set_xticks(ticks)
    ax.set_xticklabels([months[i]["label"] for i in ticks], rotation=30,
                       ha="right")
    ax.grid(False)
    if event:
        column = sum(1 for key in keys if key[:2] < (event.year, event.month))
        ax.axvline(column - 0.5, color=COLOUR["post"], lw=0.8, ls=(0, (4, 2)))
    bar = fig.colorbar(image, ax=ax, fraction=0.05, pad=0.02)
    bar.set_label("z", fontsize=6.2)
    bar.ax.tick_params(labelsize=5.8)
    # Two lines, not one: c4 is the narrowest column and its one-line title is wider
    # than the panel, which -- with ``bbox_inches="tight"`` -- widened the whole file
    # to fit the title. The break is where the column ends.
    ax.set_title("c4  ·  trajectory matrix" + chr(10) + "z per dimension × month",
                 loc="left", pad=3.0, fontsize=7.2, linespacing=1.15)

    # --- the plate's own header and legend --------------------------------
    prefix_class = report["classes"].get("prefix")
    head = caption or ("%s — OCV trajectory on the support-sized chip"
                     % report["label"])
    parts = ["n = %d chips (%s → %s)"
             % (report["cut"]["kept"], report["cut"]["window"][0],
                report["cut"]["window"][1]),
             "frozen pre-event reference%s: %d chips, %s → %s"
             % (" (rule-selected %s → %s)"
                % (ref["chosen"][0], ref["chosen"][1])
                if ref["chosen"] and ref.get("choice_rule") else
                " (chosen %s → %s)" % (ref["chosen"][0], ref["chosen"][1])
                if ref["chosen"] else "",
                ref["n_chips"], ref["window"][0], ref["window"][1])]
    if ref["outside"]["n_chips"]:
        parts.append("held-out pre-event: %d chips, %s → %s"
                     % (ref["outside"]["n_chips"],
                        ref["outside"]["window"][0],
                        ref["outside"]["window"][1]))
    if prefix_class and prefix_class["n_chips"]:
        parts.append(("before the frozen reference: %d chips, %s → %s")
                     % (prefix_class["n_chips"],
                        prefix_class["window"][0], prefix_class["window"][1]))
    if post:
        parts.append("post-event %s: %d chips, %s → %s"
                     % ("/".join(post["states"]), post["n_chips"],
                        post["window"][0], post["window"][1]))
    if event:
        parts.append("event %s" % report["event"])
    # The header is packed to the panel block's own width and starts at the block's own
    # left edge, measured rather than guessed: flush with the leftmost ink of the block
    # (the y labels of column one) and ending at the right spine of the last column
    # (a4/b4 span it in full; c4's is narrowed by its colourbar). Written as one line
    # the header is wider than the plate, and ``bbox_inches="tight"`` would then grow
    # the *file* to fit the header instead of the header to fit the plate.
    renderer = fig.canvas.get_renderer()
    panel_left = min(ax.get_tightbbox(renderer).x0
                     for ax in where.values()) / fig.bbox.width
    panel_right = max(ax.get_position().x1 for ident, ax in where.items()
                      if ident.endswith("4"))
    fig.suptitle(packed(fig, [head] + parts, 8.0,
                        (panel_right - panel_left) * COL_WIDTH),
                 fontsize=8.0, x=panel_left, ha="left", y=0.995, va="top",
                 linespacing=1.30)
    handles = [
        Line2D([], [], color=COLOUR["reference"], marker="o", ms=3.2, lw=0.9,
               label="frozen pre-event reference month"),
        Patch(facecolor=COLOUR["band"], alpha=0.75,
              label="frozen reference window"),
        Patch(facecolor=COLOUR["envelope"], alpha=0.55,
              label="reference envelope (q05–q95)"),
    ]
    if ref["outside"]["n_chips"]:
        handles.append(
            Line2D([], [], color=COLOUR["outside"], marker="o", ms=3.2, lw=0.9,
                   mfc="white", label="held-out pre-event month"))
    if prefix_class and prefix_class["n_chips"]:
        handles.append(
            Line2D([], [], color=COLOUR["prefix"], marker="o", ms=3.2, lw=0.9,
                   mfc="white", label="month before the frozen reference"))
    if post:
        handles.append(
            Line2D([], [], color=COLOUR["post"], marker="s", ms=3.2, lw=0.9,
                   ls=(0, (3, 2)), mfc="white", label="post-event month"))
        if any(month_class(m) == "mixed" for m in months):
            handles.append(
                Line2D([], [], color=COLOUR["mixed"], marker="D", ms=3.2, lw=0,
                       mfc="white", label="month on both sides"))
    if event:
        handles.append(Line2D([], [], color=COLOUR["post"], lw=0.8,
                              ls=(0, (4, 2)), label="the event"))
    # The key is the plate's only legend, so it is re-flowed until it fits the column
    # rather than allowed to widen the file: widest row first, then more rows.
    fit_legend(fig, handles, (panel_right - panel_left) * COL_WIDTH,
               (panel_left + panel_right) / 2.0)
    fig.savefig(path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    return path


def number(value, digits=3):
    """``value`` with ``digits`` significant digits, or an en dash."""
    if value is None:
        return "–"
    if value == int(value) and abs(value) < 1e12:
        return "%d" % int(value)
    return "%.*g" % (digits, value)


def signed(value, digits=2):
    return "–" if value is None else "%+.*f" % (digits, value)


def reference_lead(ref):
    """The sentence that says where the frozen window comes from.

    One place, so the report cannot describe a rule-selected window without
    saying which rule selected it, and cannot describe a chosen one as a rule.
    """
    rule = ref.get("choice_rule")
    if rule:
        return ("The reference window is selected by a **stated rule** (%s), "
                "which reads acquisition dates and per-month counts only and "
                "is frozen before any held-out number is computed;" % rule)
    return "The reference window is a **stated choice per asset** (R2):"


def report_md(report):
    """The site's report -- the same numbers as the plate, never a new one."""
    site = report["site"]
    ref = report["reference"]
    post = report["post"]
    out = []
    add = out.append
    add("# %s — OCV trajectory on the support-sized chip" % report["label"])
    add("")
    add("**Site** %s · **event** %s · **reference state** `%s` (%s) · **table** `%s` "
        "(sha256 `%s`)"
        % (report["label"], report["event"] or "none", ref["state"],
           ref["policy"], report["table"], report["sha256"][:16]))
    add("")
    add("```")
    add("standardisation  z = (t(x) - median_ref) / (1.4826 * MAD_ref),  |z| <= %g" % Z_CLIP)
    add("                 t = log for gamma2, P, A, D, F; identity for S")
    add("blocks           behaviour = %s   geometry = %s   (ocv_reconfigure.BLOCKS)"
        % (list(BLOCKS["behaviour"]), list(BLOCKS["geometry"])))
    add("reference per block  behaviour: echo-bearing reference chips; "
        "geometry: all reference chips (R3)")
    add("bootstrap        %d circular moving-block draws, block %d, seed %d (R5)"
        % (BOOT, BLOCK, SEED))
    add("```")
    add("")
    add("## The window")
    add("")
    add("| part | state | window | chips | ASC | DESC | echo-bearing | months |")
    add("| --- | --- | --- | --- | --- | --- | --- | --- |")
    add("| frozen pre-event reference (chosen) | `%s` | %s → %s | **%d** | %d | %d | %d | %d |"
        % (ref["state"], ref["window"][0], ref["window"][1], ref["n_chips"],
           ref["n_ascending"], ref["n_descending"], ref["n_echo"],
           ref["months"]))
    outside = ref["outside"]
    if outside["n_chips"]:
        add("| held-out pre-event | `%s` | %s → %s | %d | %d | %d | %d | %d |"
            % (ref["state"], outside["window"][0], outside["window"][1],
               outside["n_chips"], outside["n_ascending"],
               outside["n_descending"], outside["n_echo"], outside["months"]))
    prefix_class = report["classes"].get("prefix")
    if prefix_class and prefix_class["n_chips"]:
        add(("| before the frozen reference | `%s` | %s → %s | %d | %d "
             "| %d | %d | %d |")
            % (ref["state"], prefix_class["window"][0],
               prefix_class["window"][1], prefix_class["n_chips"],
               prefix_class["n_ascending"], prefix_class["n_descending"],
               prefix_class["n_echo"], prefix_class["months"]))
    if ref["chosen"]:
        add("| frozen pre-event state | `%s` | %s → %s | %d | %d | %d | %d | %d |"
            % (ref["state"], ref["span"]["window"][0],
               ref["span"]["window"][1], ref["span"]["n_chips"],
               ref["span"]["n_ascending"], ref["span"]["n_descending"],
               ref["span"]["n_echo"], ref["span"]["months"]))
    if post:
        add("| tail | `%s` | %s → %s | **%d** | %d | %d | %d | %d |"
            % ("/".join(post["states"]), post["window"][0],
               post["window"][1], post["n_chips"], post["n_ascending"],
               post["n_descending"], post["n_echo"], post["months"]))
    add("")
    if ref["chosen"]:
        add("%s `%s` → `%s`, inclusive calendar months. The `%s` state "
            "holds **%d** chips in all (%s → %s); the choice keeps **%d** "
            "and leaves **%d** chips outside it on the pre-event side "
            "(%s → %s) — drawn in grey, standardising nothing."
            % (reference_lead(ref), ref["chosen"][0], ref["chosen"][1],
               ref["state"], ref["span"]["n_chips"],
               ref["span"]["window"][0], ref["span"]["window"][1],
               ref["n_chips"], outside["n_chips"], outside["window"][0],
               outside["window"][1]))
    else:
        add("The reference window is the state's own span (R2, `None`): no "
            "sub-window is chosen, so all **%d** chips of `%s` (%s → %s) "
            "standardise, and the plate draws nothing outside it."
            % (ref["span"]["n_chips"], ref["state"],
               ref["span"]["window"][0], ref["span"]["window"][1]))
    add("")
    add("The cut holds **%d** of the table's **%d** rows (%d refused, %s excluded); "
        "the plate draws %d monthly bins from %s to %s."
        % (report["cut"]["kept"], report["cut"]["rows"],
           len(report["cut"]["refused"]),
           ", ".join("`%s` × %d" % (k, v) for k, v in
                     sorted(report["cut"]["excluded"].items())) or "none",
           len(report["months"]), report["months"][0]["month"],
           report["months"][-1]["month"]))
    add("")
    if report["cut"]["excluded"]:
        for state, reason in sorted(SITES[site]["excluded"].items()):
            add("* `%s` — %s." % (state, reason))
        add("")
    add("### The three frozen classes")
    add("")
    add("| class | label in this report | chips | ASC | DESC | echo-bearing | monthly intervals |")
    add("| --- | --- | --- | --- | --- | --- | --- |")
    classes = report["classes"]
    for name in ("reference", "outside", "post"):
        cell = classes[name]
        if not cell["n_chips"]:
            add("| `%s` | %s | 0 | – | – | – | 0 |" % (name, cell["label"]))
            continue
        add("| `%s` | %s | **%d** | %d | %d | %d | %d |"
            % (name, cell["label"], cell["n_chips"], cell["n_ascending"],
               cell["n_descending"], cell["n_echo"], cell["months"]))
    add("")
    add(classes["event_rule"] + ".")
    add("")
    reconciliation = classes["reconciliation"]
    counted = ["**%d** reference" % reconciliation["reference"],
               "**%d** held-out pre-event" % reconciliation["outside"]]
    if reconciliation.get("prefix"):
        counted.append("**%d** before the frozen reference"
                       % reconciliation["prefix"])
    counted += ["**%d** post-event" % reconciliation["post"],
                "**%d** unassigned" % reconciliation["unassigned"]]
    add("The classes reconcile with the cut: **%d** kept = %s. The event month "
        "is split into a pre and a post bin wherever the record holds chips on "
        "both sides, so no monthly aggregation combines the two."
        % (reconciliation["kept"], " + ".join(counted)))
    add("")
    add("## The reference (R3 / R4)")
    add("")
    for label in ref["groups"]:
        add("**%s**%s" % (label, " — per orbit, never pooled"
                          if label != "pooled" else ""))
        add("")
        add("| dimension | block | standardised against | n | t | median | scale | from | z range |")
        add("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
        for dim in DIMS:
            cell = report["standardisation"][label][dim]
            add("| `%s` | %s | %s | %d | %s | %s | %s | %s | %s … %s |"
                % (dim, BLOCK_OF[dim], cell["rule"], cell["n"],
                   cell["transform"], number(cell["centre"], 4),
                   number(cell["scale"], 4), cell["source"],
                   signed(cell["z_min"]), signed(cell["z_max"])))
        add("")
    if report["weather_reference"]:
        add("### Environmental context (panel c1)")
        add("")
        add("Standardised on the same reference, on their own (identity) scale.")
        add("")
        add("| column | group | n | distinct | median | scale | from |")
        add("| --- | --- | --- | --- | --- | --- | --- |")
        thin = {}
        for column in sorted(report["weather_reference"]):
            per_group = report["weather_reference"][column]
            for label in sorted(per_group):
                cell = per_group[label]
                add("| `%s` | %s | %d | %d | %s | %s | %s |"
                    % (column, label, cell["n"], cell["distinct"],
                       number(cell["centre"], 4), number(cell["scale"], 4),
                       cell["source"]))
                if cell["distinct"] <= 12 and cell["n"] > cell["distinct"]:
                    previous = thin.get(column)
                    pair = (cell["distinct"], cell["n"])
                    thin[column] = pair if previous is None \
                        else (min(previous[0], pair[0]), max(previous[1], pair[1]))
        add("")
        if thin:
            add("A co-variate that carries only a handful of distinct values is a "
                "climatology, not a per-acquisition measurement — %s — so panel c1 "
                "steps between those values and cannot resolve anything finer. It "
                "is drawn for exactly that reason: to keep the reader from reading "
                "the other panels' wiggles as weather."
                % "; ".join("`%s` (%d values over up to %d chips)" % (column, pair[0],
                                                                     pair[1])
                            for column, pair in sorted(thin.items())))
            add("")
    if report["clipped"]:
        add("**%d** of the standardised values were clipped at ±%g; they are "
            "counted here and drawn at the clip, never dropped."
            % (report["clipped"], Z_CLIP))
        add("")
    add("## The reference distribution and the monthly envelope")
    add("")
    envelope = report["reference_envelope"]
    add("The decision distribution is built at **monthly** aggregation: each of "
        "the **%d** reference months contributes one state (the monthly median of "
        "its chips), and the envelope is the empirical q05/q50/q95 of those "
        "states. A held-out month is compared only with this monthly distribution "
        "-- never with the chip-level reference and never with the reference "
        "maximum, which is reported below for transparency and is **not** a "
        "decision boundary. With %d monthly states the quantiles are empirical "
        "order statistics, not confidence bounds."
        % (envelope["n_reference_months"], envelope["n_reference_months"]))
    add("")
    add("| quantity | monthly ref states | q05 | q50 | q95 | min | max | used for |")
    add("| --- | --- | --- | --- | --- | --- | --- | --- |")
    for key in list(DIMS) + ["R_B", "R_G"]:
        cell = envelope["monthly"].get(key)
        if not cell:
            continue
        add("| `%s` | %d | %s | %s | %s | %s | %s | %s |"
            % (key, cell["n"], signed(cell["q05"]), signed(cell["q50"]),
               signed(cell["q95"]), signed(cell["min"]), signed(cell["max"]),
               "above q95 = unusual magnitude" if key in ("R_B", "R_G")
               else "outside [q05, q95] = signed departure"))
    add("")
    add("The chip-level distributions are stored beside them in the JSON "
        "(`reference_envelope.chip`) and are **descriptive only**: nothing in the "
        "held-out assessment compares a monthly state with a single chip.")
    add("")
    add("## The monthly trajectory")
    add("")
    add("| month | n | ASC | DESC | echo | " +
        " | ".join("z(%s)" % DIM_MD[d] for d in DIMS) +
        " | R_B | R_G |")
    add("| --- | --- | --- | --- | --- | " + " | ".join("---" for _ in DIMS) +
        " | --- | --- |")
    for month in report["months"]:
        add("| %s | %d | %d | %d | %d | %s | %s | %s |"
            % (month["label"], month["n"], month["n_ascending"],
               month["n_descending"], month["n_echo"],
               " | ".join(signed(month["z"][d]) for d in DIMS),
               signed(month["R_B"]), signed(month["R_G"])))
    add("")
    add("A block magnitude is computed **per chip** (the RMS of that chip's "
        "standardised dimensions, R4) and the panel and the table show the "
        "**median of those chip values** — so `R_B` is the typical joint "
        "displacement of a single overpass. It is therefore *not* the RMS of the "
        "two columns beside it, and it does not vanish when those two medians are "
        "near zero: a month whose chips scatter widely in opposite directions has "
        "a zero median in each dimension and a large magnitude per chip.")
    add("")
    add("## Frozen reference / held-out pre-event / post-event")
    add("")
    add("Each cell is the median over the chips of that class, in the reference's "
        "own scales. The reference column is zero by construction (the "
        "standardisation rests on its median); it is printed so the other two "
        "columns are read against the frame they are measured in.")
    add("")
    comparison = report["class_comparison"]
    add("| quantity | frozen pre-event reference | held-out pre-event | post-event |")
    add("| --- | --- | --- | --- |")
    for key in comparison["order"]:
        cells = []
        for name in ("reference", "outside", "post"):
            value = comparison["classes"][name]["values"].get(key)
            cells.append("–" if value is None else signed(value))
        add("| `%s` | %s | %s | %s |" % (key, cells[0], cells[1], cells[2]))
    add("")
    add("n = **%d** / **%d** / **%d** chips over **%d** / **%d** / **%d** monthly "
        "intervals."
        % (comparison["classes"]["reference"]["n_chips"],
           comparison["classes"]["outside"]["n_chips"],
           comparison["classes"]["post"]["n_chips"],
           comparison["classes"]["reference"]["n_intervals"],
           comparison["classes"]["outside"]["n_intervals"],
           comparison["classes"]["post"]["n_intervals"]))
    add("")
    add("## The contrast (post minus reference, in reference scales)")
    add("")
    con = report["contrast"]
    if con.get("applies"):
        add("Circular moving-block bootstrap, **%d** draws, block **%d**, seed "
            "**%d**, one resample per group; the reference side holds **%d** "
            "chips and the tail **%d**."
            % (con["draws"], con["block"], con["seed"], con["n_reference"],
               con["n_post"]))
        add("")
        add("| quantity | n_ref | n_post | delta z | 95 % CI | verdict |")
        add("| --- | --- | --- | --- | --- | --- |")
        for key in list(DIMS) + ["R_B", "R_G"]:
            cell = con["dims"].get(key)
            if not cell:
                continue
            add("| `%s` | %d | %d | %s | %s … %s | %s |"
                % (key, cell["n_reference"], cell["n_post"],
                   signed(cell["delta"]), signed(cell["ci"][0]),
                   signed(cell["ci"][1]),
                   "**excludes zero**" if cell["excludes_zero"]
                   else "spans zero"))
        add("")
        add("The interval is on the **difference of medians**, not on the curve: "
            "a month that moves inside it is not a finding, and a month that "
            "leaves it is not yet one either — %s."
            % ("the tail is the rebuild, so the contrast answers a different "
               "question than the collapse plates of the other sites"
               if site == "kdlo" else
               "this record carries one reference and one tail, not a series "
               "of states"))
    else:
        add(con["reason"] + ".")
    add("")
    add("## The held-out pre-event class (OUTSIDE)")
    add("")
    oc = report["outside"]
    if oc.get("applies"):
        add("OUTSIDE is the third frozen class: **%d** chips of `%s` that the "
            "reference window leaves out on the pre-event side (%s → %s), over "
            "**%d** monthly intervals (%s). They standardise nothing -- the centre "
            "and the scale are the reference's -- so every number below reads them "
            "in the frozen reference's own frame. The interval is the plate's own "
            "bootstrap (**%d** draws, block **%d**, seed **%d**, one resample per "
            "side)."
            % (oc["n_outside"], oc["state"], oc["window"][0], oc["window"][1],
               oc["n_intervals"], ", ".join(oc["intervals"]), oc["draws"],
               oc["block"], oc["seed"]))
        add("")
        add("**Signed channels** -- `r` is the signed median departure in "
            "transform units, `z = r / scale_ref` the same departure in reference "
            "scales, and the interval is on the difference of medians.")
        add("")
        add("| quantity | n chips | n int | r | z | 95 % CI (z) | q05 | q95 | monthly level | position | months outside envelope |")
        add("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
        for dim in DIMS:
            cell = oc["dims"].get(dim)
            if not cell:
                continue
            positioning = cell["reference_comparison"]
            add("| `%s` | %d | %d | %s | %s | %s … %s | %s | %s | %s | %s | %s |"
                % (dim, cell["n_chips"], cell["n_intervals"],
                   signed(cell["r"], 3), signed(cell["z"]),
                   signed(cell["ci"][0]), signed(cell["ci"][1]),
                   signed(positioning["q05"]), signed(positioning["q95"]),
                   signed(positioning["level"]), positioning["position"],
                   ", ".join(positioning["outside_envelope"]) or "none"))
        add("")
        add("**Block magnitudes** -- a magnitude has no sign, so the envelope test "
            "is one-sided: only a monthly value above q95 counts as an exceedance. "
            "A value below q05 is reported as a *lower* magnitude, never as an "
            "anomaly. `position` is where the class's median monthly state sits; "
            "the last-but-one column lists the individual months.")
        add("")
        add("| quantity | n chips | n int | median R | delta | 95 % CI | q95 | monthly level | position | months above q95 | longest run | sustained |")
        add("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
        for name in ("R_B", "R_G"):
            cell = oc["dims"].get(name)
            runs = oc["persistence"].get(name, {})
            if not cell:
                continue
            positioning = cell["reference_comparison"]
            longest = runs.get("longest_run")
            add("| `%s` | %d | %d | %s | %s | %s … %s | %s | %s | %s | %s | %s | %s |"
                % (name, cell["n_chips"], cell["n_intervals"],
                   number(cell["value"]), signed(cell["delta"]),
                   signed(cell["ci"][0]), signed(cell["ci"][1]),
                   signed(positioning["q95"]), signed(positioning["level"]),
                   positioning["position"],
                   ", ".join(positioning["months_above_q95"]) or "none",
                   ("%d (%s → %s)" % (longest["length"], longest["start"],
                                       longest["end"])) if longest else "none",
                   "yes" if runs.get("sustained") else "no"))
        add("")
        add("**Criterion (fixed before the numbers were read).** " + oc["criterion"]
            + ".")
        add("")
        directional = oc["directional"]
        pattern_run = directional.get("longest_run")
        add("**Directional pattern (descriptive only).** The joint pattern `%s` "
            "holds in **%d** held-out months; its longest run is **%d** month(s) "
            "(%s). No statistic is derived from these runs."
            % (directional["pattern"], directional.get("n_months", 0),
               (pattern_run or {"length": 0})["length"],
               ("%s → %s" % (pattern_run["start"], pattern_run["end"]))
               if pattern_run else "no run"))
        add("")
    else:
        add(oc["reason"] + ".")
        add("")
    add("## Held-out pre-event assessment")
    add("")
    assessment = report["assessment"]
    if assessment.get("applies"):
        add("The three statements are kept apart throughout:")
        add("")
        add("* **Directional reconfiguration** -- a signed monthly median leaves "
            "the monthly reference envelope; a statement about direction, not size.")
        add("* **Reconfiguration magnitude** -- R_B or R_G above the monthly "
            "reference q95; a statement about the joint size of the displacement, "
            "with no direction.")
        add("* **Unusual reference departure** -- only the magnitude statement, and "
            "only when the pre-specified criterion is met; a run of at least %d "
            "consecutive months makes it sustained." % PERSIST_MIN)
        add("")
        for question in assessment["questions"].values():
            add("**Q. %s**" % question["question"])
            add("")
            add(question["answer"])
            add("")
        add("**Verdict — %s.** It is derived from the R_B and R_G persistence "
            "rules alone (`verdict_inputs`: %s); the Summer-2024 subset and the "
            "environment are not inputs. It does not claim damage, a collapse "
            "prediction, a failure precursor or causation."
            % (assessment["verdict"],
               ", ".join("`%s`" % key for key in assessment["verdict_inputs"])))
        add("")
    else:
        add(assessment["reason"] + ".")
        add("")
    summer = report["summer_2024"]
    add("### Summer 2024 — post-hoc descriptive subset")
    add("")
    if summer.get("applies"):
        add("**Post-hoc and descriptive only -- no interval, no p-value, no new "
            "statistic, no verdict input.** The window `%s` → `%s` was fixed after "
            "OUTSIDE was defined, and it is a strict subset of it: **%d** chips "
            "over **%d** monthly intervals (%s). Medians: z(γ²) %s, z(P) %s, "
            "R_B %s, R_G %s; months above the monthly q95: R_B %s, R_G %s."
            % (summer["window"][0], summer["window"][1], summer["n_chips"],
               summer["n_intervals"], ", ".join(summer["intervals"]),
               signed(summer["medians"]["gamma2"]), signed(summer["medians"]["P"]),
               number(summer["medians"]["R_B"]),
               number(summer["medians"]["R_G"]),
               ", ".join(summer["months_above_q95"]["R_B"]) or "none",
               ", ".join(summer["months_above_q95"]["R_G"]) or "none"))
        add("")
    else:
        add(summer["reason"] + ".")
        add("")
    environment = report["environment_context"]
    add("### Environmental context of the held-out departures")
    add("")
    if environment.get("applies"):
        add(environment["note"] + ".")
        add("")
        variables = environment["variables"]
        add("| month | R_B | R_G | z(γ²) | z(P) | "
            + " | ".join(WEATHER_LABEL.get(column, column)
                         for column in variables) + " | above q95 |")
        add("| --- | --- | --- | --- | --- | "
            + " | ".join("---" for _ in variables) + " | --- |")
        for row in environment["months"]:
            add("| %s | %s | %s | %s | %s | %s | %s |"
                % (row["month"], number(row["R_B"]), number(row["R_G"]),
                   signed(row["gamma2"]), signed(row["P"]),
                   " | ".join(signed(row["weather_z"].get(column))
                             for column in variables),
                   ", ".join(name for name, flag in row["above_q95"].items()
                             if flag) or "none"))
        add("")
        if environment["coincident_months"]:
            add("The month(s) above the magnitude envelope -- %s -- coincide with "
                "the co-variation recorded above. This is a temporal coincidence, "
                "not an explanation: OCV is not weather-corrected in this "
                "repository." % ", ".join(environment["coincident_months"]))
            add("")
    else:
        add("No held-out class, so there is no held-out departure to compare "
            "with the environment.")
        add("")
    add("## The panels")
    add("")
    add("| panel | block | what it draws |")
    add("| --- | --- | --- |")
    for panel in report["figure"]["panels"]:
        add("| **%s** | %s | %s |" % (panel["id"], panel["block"],
                                      panel["shows"]))
    add("")
    add("## What the plate does not claim")
    add("")
    add("* **No damage verdict.** Every ordinate is a deviation from the case's "
        "*own* reference, in that reference's robust scales — two plates are not "
        "comparable in magnitude, which is why the comparison they support is a "
        "comparison of *shape* and of *where the movement sits*.")
    add("* **No significance on the curve.** Both intervals of the plate — the "
        "contrast's and the held-out class's — are differences of medians, not "
        "bands on the curve: the monthly medians are drawn with their `n` and "
        "without one.")
    add("* **No directional score.** `R_B` and `R_G` keep their frozen "
        "definitions; no directional or persistence score is added, and the "
        "verdict never reads one.")
    add("* **The reference maximum is not a boundary.** The envelope is the "
        "q05/q50/q95 of the monthly reference states; the maximum is printed "
        "only so the reader can see how far q95 sits from it.")
    add("* **The Summer-2024 subset is post-hoc.** It is descriptive, untested "
        "and not an input to the verdict.")
    add("* **The environment is beside the OCV numbers, never inside them.** "
        "No weather correction is applied and no causal attribution is made.")
    add("* **No interpolation.** A month without a chip, the split event month "
        "and the event itself break the line: nothing is drawn across a hole.")
    add("* **A month is a bin, not a date.** Its median can mix both orbits; "
        "panel c2 says how many chips of each the bin holds, and on the sites "
        "with two orbits the two directions alternate in blocks of days.")
    add("* **No pooling across sites.** The plate is per site by construction; "
        "`figures/cases/fig_ocv_trajectory.sh` runs the six and nothing else.")
    if not report["post"]:
        add("* **No contrast at all on this site.** The event-free control has "
            "no tail, so the plate is the trajectory alone — which is the point "
            "of a control.")
    if report["event"] is None:
        add("* **No event line.** There is no event in this record to draw, so "
            "the plate carries none; the reference window is the whole record.")
    add("")
    add("## Pins")
    add("")
    add("`--verify` re-derives every number of this report from `%s` and fails "
        "if the JSON, the report and the table disagree; the table's sha256 is "
        "`%s`."
        % (report["table"], report["sha256"]))
    add("")
    add("The assessment's pins are constants, asserted by `--selftest`: seed "
        "**%d**, block **%d**, draws **%d**, decision quantiles **%s**, the "
        "magnitude boundary **q%d**, and a persistence rule of **%d** consecutive "
        "months. The bootstrap is seeded and no wall-clock value enters any "
        "output, so the JSON, the report and the plate reproduce deterministically."
        % (SEED, BLOCK, BOOT, "/".join("q%02d" % q for q in QUANTILES),
           DECISION_QUANTILE, PERSIST_MIN))
    add("")
    return "\n".join(out)


#: The documented windows of ``code/cases/ocv_trajectory_spec.md``, as numbers.
#: ``--selftest`` asserts every one of them against the table, so the spec and the
#: code cannot drift apart without a failing check.
#:
#: ``reference`` is the *chosen* reference (R2), ``state`` the whole reference
#: state it was chosen from, ``outside`` the healthy chips the choice leaves out
#: (drawn, standardising nothing). Only a site with a ``chosen`` window can have
#: chips outside, and for the rest ``state == reference`` and ``outside == 0``.
EXPECTED = {
    "carola": {"reference": 121, "ascending": 81, "descending": 40,
               "post": 252, "reference_window": ("2022-09-03", "2023-12-27"),
               "post_window": ("2024-09-13", "2026-08-23"), "excluded": {},
               "chosen": ("2022-09", "2023-12"), "state": 184,
               "state_window": ("2022-09-03", "2024-09-09"),
               "outside": 63, "outside_window": ("2024-01-01", "2024-09-09")},
    "cts": {"reference": 136, "ascending": 136, "descending": 0, "post": 42,
            "reference_window": ("2015-09-21", "2021-06-21"),
            "post_window": ("2021-07-03", "2022-12-25"), "excluded": {},
            "chosen": None, "state": 136,
            "state_window": ("2015-09-21", "2021-06-21"), "outside": 0},
    "morandi": {"reference": 198, "ascending": 91, "descending": 107,
                "post": 9, "reference_window": ("2016-08-24", "2018-08-08"),
                "post_window": ("2018-08-14", "2018-09-13"), "excluded": {},
                "chosen": None, "state": 198,
                "state_window": ("2016-08-24", "2018-08-08"), "outside": 0},
    "kdlo": {"reference": 16, "ascending": 0, "descending": 16, "post": 43,
             "reference_window": ("2022-06-02", "2022-12-11"),
             "post_window": ("2023-05-04", "2024-09-19"),
             "excluded": {"single": 1}, "chosen": None, "state": 16,
             "state_window": ("2022-06-02", "2022-12-11"), "outside": 0},
    "haltern": {"reference": 96, "ascending": 44, "descending": 52, "post": 11,
                "reference_window": ("2021-01-01", "2021-09-17"),
                "post_window": ("2021-10-01", "2021-10-17"), "excluded": {},
                "chosen": None, "state": 96,
                "state_window": ("2021-01-01", "2021-09-17"), "outside": 0},
    "espoo": {"reference": 75, "ascending": 31, "descending": 44, "post": 0,
              "reference_window": ("2024-09-05", "2026-09-07"),
              "post_window": None, "excluded": {}, "chosen": None, "state": 75,
              "state_window": ("2024-09-05", "2026-09-07"), "outside": 0},
}


def first_difference(left, right, path="root"):
    """The first path at which two JSON structures disagree (``None`` if equal)."""
    if type(left) is not type(right) and not (
            isinstance(left, (int, float)) and isinstance(right, (int, float))):
        return "%s: %r vs %r" % (path, left, right)
    if isinstance(left, dict):
        for key in sorted(set(left) | set(right)):
            if key not in left or key not in right:
                return "%s/%s: present on one side only" % (path, key)
            found = first_difference(left[key], right[key], "%s/%s" % (path, key))
            if found:
                return found
        return None
    if isinstance(left, list):
        if len(left) != len(right):
            return "%s: %d vs %d items" % (path, len(left), len(right))
        for index, (a, b) in enumerate(zip(left, right)):
            found = first_difference(a, b, "%s[%d]" % (path, index))
            if found:
                return found
        return None
    if isinstance(left, float) and isinstance(right, float):
        if abs(left - right) > 1e-9:
            return "%s: %r vs %r" % (path, left, right)
        return None
    if left != right:
        return "%s: %r vs %r" % (path, left, right)
    return None


def selftest():
    """The rules of the template, checked without writing a single file."""
    checks = []

    def check(name, ok, detail=""):
        checks.append((name, bool(ok), detail))

    check("registry: exactly the six case sites",
          sorted(SITES) == sorted(EXPECTED), "%s" % sorted(SITES))
    problems = []
    for site, spec in SITES.items():
        states = set(spec["states"])
        post = set(spec["post"])
        excluded = set(spec["excluded"])
        if spec["reference"] is not None and spec["reference"] not in states:
            problems.append("%s: reference %r" % (site, spec["reference"]))
        if not post <= states:
            problems.append("%s: post %s" % (site, sorted(post - states)))
        if excluded & post:
            problems.append("%s: excluded %s is a post state"
                            % (site, sorted(excluded & post)))
        if not excluded <= states:
            problems.append("%s: excluded %s" % (site, sorted(excluded - states)))
        if spec["event"] is None and post:
            problems.append("%s: an event-free site with a tail" % site)
        if spec["event"] is not None and not post:
            problems.append("%s: an event without a tail" % site)
    check("registry: reference, post and excluded states are consistent",
          not problems, "; ".join(problems))
    check("blocks: every dimension sits in exactly one block",
          set(BLOCK_OF) == set(DIMS) and len(BLOCK_OF) == len(DIMS),
          str(BLOCK_OF))
    check("blocks: the plane dimensions belong to their block",
          set(GEOMETRY_PLANE) <= set(BLOCKS["geometry"])
          and set(BLOCKS["behaviour"]).isdisjoint(BLOCKS["geometry"]))

    try:
        sys.path.insert(0, os.path.join(REPO, "code", "sds"))
        import ocv_reconfigure as reconfig          # noqa: PLC0415
    except Exception as error:                       # noqa: BLE001
        reconfig = None
        check("ocv_reconfigure: importable", False, "%s: %s"
              % (type(error).__name__, error))
    if reconfig is not None:
        check("ocv_reconfigure.BLOCKS matches this template",
              {k: tuple(v) for k, v in reconfig.BLOCKS.items()} == BLOCKS,
              str(reconfig.BLOCKS))
        check("ocv_reconfigure.GEOMETRY_PLANE matches this template",
              tuple(reconfig.GEOMETRY_PLANE) == GEOMETRY_PLANE,
              str(reconfig.GEOMETRY_PLANE))
        # Bautzen's gamma2 is floor-corrected and signed, so it stays on identity
        # there; what this template must agree with is the rule on the sites both
        # registries know (Carola, Morandi, Haltern, KDLO).
        shared = [site for site in sorted(SITES) if site in reconfig.TRANSFORM]
        drift = [(site, dim, reconfig.TRANSFORM[site][dim], TRANSFORM[dim])
                 for site in shared for dim in DIMS
                 if reconfig.TRANSFORM[site][dim] != TRANSFORM[dim]]
        check("ocv_reconfigure.TRANSFORM agrees on the four shared sites",
              not drift and len(shared) == 4,
              "%s %s" % (shared, drift))

    loaded = {}
    for site in sorted(SITES):
        path = table_path(site)
        check("%s: the support table exists" % site, os.path.exists(path),
              os.path.relpath(path, REPO))
        if not os.path.exists(path):
            continue
        chips, audit = load_chips(site)
        loaded[site] = chips
        want = EXPECTED[site]
        ref = [c for c in chips if in_reference(site, c)]
        post = [c for c in chips if c["state"] in SITES[site]["post"]]
        check("%s: reference window holds %d chips" % (site, want["reference"]),
              len(ref) == want["reference"], "got %d" % len(ref))
        check("%s: reference split ASC %d / DESC %d"
              % (site, want["ascending"], want["descending"]),
              (sum(1 for c in ref if c["orbit"] == "ASCENDING"),
               sum(1 for c in ref if c["orbit"] == "DESCENDING"))
              == (want["ascending"], want["descending"]))
        check("%s: reference window dates" % site,
              window(ref) == list(want["reference_window"]), str(window(ref)))
        chosen = SITES[site]["reference_window"]
        check("%s: the chosen reference window is the pinned one" % site,
              chosen == want["chosen"], "%r" % (chosen,))
        ref_chips, outside, state = reference_split(site, chips)
        check("%s: the reference state holds %d chips" % (site, want["state"]),
              len(state) == want["state"], "got %d" % len(state))
        check("%s: the reference state keeps its own span" % site,
              window(state) == list(want["state_window"]), str(window(state)))
        check("%s: the choice leaves %d healthy chips outside" % (site,
                                                                want["outside"]),
              len(outside) == want["outside"], "got %d" % len(outside))
        check("%s: reference and outside partition the reference state" % site,
              len(ref_chips) == len(ref)
              and len(ref_chips) + len(outside) == len(state),
              "%d + %d vs %d" % (len(ref_chips), len(outside), len(state)))
        if want["outside"]:
            check("%s: the outside window dates" % site,
                  window(outside) == list(want["outside_window"]),
                  str(window(outside)))
            check("%s: every outside chip is outside the chosen months" % site,
                  all(not within(c["day"], chosen) for c in outside))
            check("%s: the realised reference lies inside the chosen months"
                  % site, chosen[0] <= window(ref_chips)[0][:7]
                  and window(ref_chips)[1][:7] <= chosen[1],
                  "%s .. %s in %s" % (window(ref_chips)[0], window(ref_chips)[1],
                                      chosen))
        else:
            check("%s: with no chosen window nothing is left outside" % site,
                  chosen is None and len(state) == len(ref))
        check("%s: the reference is big enough to standardise (%d chips, %d "
              "months)" % (site, len(ref_chips),
                           len({month_key(c["day"]) for c in ref_chips})),
              len(ref_chips) >= 12
              and len({month_key(c["day"]) for c in ref_chips}) >= 3)
        check("%s: the tail holds %d chips" % (site, want["post"]),
              len(post) == want["post"], "got %d" % len(post))
        if want["post_window"]:
            check("%s: tail window dates" % site,
                  window(post) == list(want["post_window"]), str(window(post)))
        check("%s: the registry's excluded states are counted" % site,
              audit["excluded"] == want["excluded"], str(audit["excluded"]))
        check("%s: every chip of the cut carries an echo (A >= %g)"
              % (site, ECHO_MIN), all(c["echo"] for c in chips),
              "%d of %d" % (sum(1 for c in chips if c["echo"]), len(chips)))
        bad = [(c["day"].isoformat(), dim) for c in chips for dim in DIMS
               if TRANSFORM[dim] == "log"
               and (c["values"][dim] is None or c["values"][dim] <= 0.0)]
        check("%s: every log channel is strictly positive" % site, not bad,
              str(bad[:3]))
        if SITES[site]["event"]:
            event = dt.date.fromisoformat(SITES[site]["event"])
            inside = (chips[0]["day"] <= event <= chips[-1]["day"])
            check("%s: the event lies inside the record" % site, inside,
                  "%s in %s .. %s" % (SITES[site]["event"],
                                      chips[0]["day"].isoformat(),
                                      chips[-1]["day"].isoformat()))
        stats, weather_stats, clipped = standardise(site, chips)
        zero = []
        for (dim, label), _cell in stats.items():
            group = [c for c in chips
                     if in_reference(site, c) and group_label(site, c) == label]
            values = [c["z"][dim] for c in group if c["z"][dim] is not None]
            zero.append(abs(median(values)))
        check("%s: the reference median standardises to zero" % site,
              max(zero) < 1e-9, "max |median z| = %.3g" % max(zero))
        check("%s: every reference scale is positive" % site,
              all(cell["scale"] > 0.0 for cell in stats.values()))
        sizes = {(dim, label): cell["n"] for (dim, label), cell in stats.items()}
        behaviour_ok = all(
            sizes[(dim, label)]
            == sum(1 for c in chips if in_reference(site, c)
                   and group_label(site, c) == label and c["echo"])
            for dim in BLOCKS["behaviour"] for label in
            {lab for _d, lab in stats})
        check("%s: R3 gives the behaviour block the echo-bearing reference"
              % site, behaviour_ok)
        geometry_ok = all(
            sizes[(dim, label)]
            == sum(1 for c in chips if in_reference(site, c)
                   and group_label(site, c) == label)
            for dim in BLOCKS["geometry"] for label in
            {lab for _d, lab in stats})
        check("%s: R3 gives the geometry block every reference chip"
              % site, geometry_ok)
        extremes = [abs(c["z"][dim]) for c in chips for dim in DIMS
                    if c["z"][dim] is not None]
        check("%s: no standardised value exceeds the clip" % site,
              max(extremes) <= Z_CLIP + 1e-12, "max |z| = %.4f" % max(extremes))
        months = monthly(site, chips)
        check("%s: the monthly bins partition the chips" % site,
              sum(m["n"] for m in months) == len(chips),
              "%d vs %d" % (sum(m["n"] for m in months), len(chips)))
        bad = []
        event = event_day(site)
        for month in months:
            key = (int(month["month"][:4]), int(month["month"][5:]))
            group = [c for c in chips if month_key(c["day"]) == key
                     and (month["part"] is None
                          or (c["day"] < event) == (month["part"] == "pre"))]
            for dim in DIMS:
                values = [c["z"][dim] for c in group if c["z"][dim] is not None]
                if values and abs(median(values) - month["z"][dim]) > 1e-9:
                    bad.append((month["label"], dim))
        check("%s: every monthly median is the median of its own bin" % site,
              not bad, str(bad[:3]))
        check("%s: every month is classified (and only a tail is 'post')" % site,
              all(month_class(m) in ("ref", "outside", "post", "mixed")
                  for m in months)
              and (bool(SITES[site]["post"])
                   or all(month_class(m) != "post" for m in months)))
        check("%s: no month sits on both sides of the reference window" % site,
              all(m["n_reference"] == 0 or m["n_outside"] == 0
                  for m in months))
        check("%s: outside chips only land in outside or mixed months" % site,
              all(month_class(m) in ("outside", "mixed")
                  for m in months if m["n_outside"]))
        check("%s: an outside month holds nothing but outside chips" % site,
              all(m["n_outside"] == m["n"] for m in months
                  if month_class(m) == "outside"))
        check("%s: the monthly bins account for every outside chip" % site,
              sum(m["n_outside"] for m in months) == len(outside))
        check("%s: no reference scale is degenerate" % site,
              all(cell["source"] != "degenerate" for cell in stats.values()),
              str(sorted({key[0] for key, cell in stats.items()
                          if cell["source"] == "degenerate"})))

    # R6 -- R2 declares a third class, so the plate reads it too: the healthy
    # chips the chosen window leaves out are tested against the reference with
    # the contrast's own interval, and the report prints that table.
    applying = [site for site in sorted(loaded)
                if outside_contrast(site, loaded[site]).get("applies")]
    choosing = [site for site in sorted(loaded) if EXPECTED[site]["outside"]]
    check("R6: the outside class is tested exactly where a window is chosen",
          applying == choosing, "%s vs %s" % (applying, choosing))
    for site in sorted(loaded):
        block = outside_contrast(site, loaded[site])
        want = EXPECTED[site]["outside"]
        held = sum(1 for c in loaded[site]
                   if reference_state(site, c) and not in_reference(site, c))
        check("%s: the outside block holds the %d chips the window leaves out"
              % (site, want),
              block["n_outside"] == want == held
              and block["n_reference"] == EXPECTED[site]["reference"],
              "%d / %d / %d reference %d"
              % (block["n_outside"], want, held, block["n_reference"]))
    carola_block = outside_contrast("carola", loaded["carola"])
    carola_ref = [c for c in loaded["carola"] if in_reference("carola", c)]
    carola_out = [c for c in loaded["carola"]
                  if reference_state("carola", c)
                  and not in_reference("carola", c)]

    def channel(chips, key):
        """The per-chip values of one row of the outside table."""
        return [c["z"][key] if key in DIMS else c[key] for c in chips
                if (c["z"][key] if key in DIMS else c[key]) is not None]

    drift = [(key, carola_block["dims"][key]["delta"],
              median(channel(carola_out, key)) - median(channel(carola_ref, key)))
             for key in carola_block["dims"]
             if abs(carola_block["dims"][key]["delta"]
                    - (median(channel(carola_out, key))
                       - median(channel(carola_ref, key)))) > 1e-9]
    check("R6: the outside delta is the difference of the two medians",
          not drift, str(drift[:2]))
    check("R6: the outside bootstrap is seeded and reproducible",
          outside_contrast("carola", loaded["carola"]) == carola_block)
    strays = [(key, cell["delta"], cell["ci"])
              for key, cell in carola_block["dims"].items()
              if not cell["ci"][0] <= cell["delta"] <= cell["ci"][1]]
    check("R6: every outside interval brackets its own delta", not strays,
          str(strays[:2]))

    check("R5: the block magnitude is the RMS of the available dimensions",
          rms([1.0, None, -1.0, 1.0]) == math.sqrt(1.0)
          and rms([2.0, 2.0]) == 2.0 and rms([None]) is None,
          "%r %r" % (rms([1.0, None, -1.0, 1.0]), rms([2.0, 2.0])))
    series = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    first = resampled_medians(series, 200, BLOCK, random.Random(SEED))
    second = resampled_medians(series, 200, BLOCK, random.Random(SEED))
    check("R5: the bootstrap is seeded and reproducible", first == second)
    check("R5: a constant series keeps its median exactly",
          all(v == 2.5 for v in resampled_medians([2.5] * 12, 100, BLOCK,
                                                  random.Random(SEED))))
    check("R5: the bootstrap median stays inside the observed range",
          min(first) >= min(series) and max(first) <= max(series),
          "%.3f .. %.3f" % (min(first), max(first)))

    # --- the frozen pins of the held-out assessment -----------------------
    check("assessment: seed, block and draws are the plate's frozen pins",
          (SEED, BLOCK, BOOT) == (7, 3, 2000), "%r %r %r" % (SEED, BLOCK, BOOT))
    check("assessment: the decision quantiles are q05/q50/q95",
          QUANTILES == (5, 50, 95) and DECISION_QUANTILE == 95)
    check("assessment: a sustained departure needs %d consecutive months"
          % PERSIST_MIN, PERSIST_MIN == 2)
    check("assessment: the Summer-2024 window is fixed and pre-event",
          SUMMER_2024 == ("2024-06-01", "2024-08-31")
          and SUMMER_2024[1] < SITES["carola"]["event"])
    check("assessment: the verdict vocabulary is the five allowed words",
          len(ALLOWED_VERDICTS) == len(set(ALLOWED_VERDICTS)) == 5)

    # --- the three frozen classes ------------------------------------------
    for site in sorted(loaded):
        chips = loaded[site]
        membership = class_members(site, chips)
        identifiers = {name: {id(c) for c in group}
                       for name, group in membership.items()}
        disjoint = all(not (identifiers[a] & identifiers[b])
                       for a in identifiers for b in identifiers if a < b)
        check("%s: REFERENCE, OUTSIDE, POST and unassigned are disjoint" % site,
              disjoint)
        check("%s: the classes reconcile with the cut" % site,
              sum(len(group) for group in membership.values()) == len(chips),
              str({name: len(group) for name, group in membership.items()}))

    carola_chips = loaded["carola"]
    carola_members = class_members("carola", carola_chips)
    check("carola: REFERENCE + OUTSIDE + POST reconcile with the 436 kept rows",
          (len(carola_members["reference"]), len(carola_members["outside"]),
           len(carola_members["post"]), len(carola_members["unassigned"]))
          == (121, 63, 252, 0),
          str({name: len(group) for name, group in carola_members.items()}))
    check("carola: OUTSIDE is after the reference window and strictly pre-event",
          all(c["day"] > dt.date(2023, 12, 31)
              and c["day"] < dt.date.fromisoformat(SITES["carola"]["event"])
              for c in carola_members["outside"]))

    mixed_bins = []
    for site in sorted(loaded):
        chips = loaded[site]
        split = split_event_month(site, chips)
        event = event_day(site)
        groups = {}
        for chip in chips:
            key = (month_key(chip["day"]), None)
            if split is not None and month_key(chip["day"]) == split:
                key = (month_key(chip["day"]),
                       "pre" if chip["day"] < event else "post")
            groups.setdefault(key, []).append(chip)
        for key, group in groups.items():
            if len({chip_class(site, chip) for chip in group}) > 1:
                mixed_bins.append((site, key))
    check("every monthly bin holds one class: none mixes pre- and post-event chips",
          not mixed_bins, str(mixed_bins[:3]))

    carola_months = monthly("carola", carola_chips)
    check("carola: 2024-09 is split at the event into a pre and a post bin",
          [m["label"] for m in carola_months if m["month"] == "2024-09"]
          == ["2024-09 (pre)", "2024-09 (post)"]
          and [m["n"] for m in carola_months if m["month"] == "2024-09"]
          == [2, 5])

    # The held-out class may not move the reference it is measured against.
    baseline = standardise("carola", [dict(c) for c in carola_chips])[0]
    perturbed = []
    for chip in carola_chips:
        copy = dict(chip)
        if chip_class("carola", chip) == "outside":
            copy["values"] = {key: (value * 1.5 if value is not None else None)
                              for key, value in chip["values"].items()}
        perturbed.append(copy)
    shifted = standardise("carola", perturbed)[0]
    check("carola: held-out chips move no reference median or scale",
          baseline == shifted,
          str([key for key in baseline if baseline[key] != shifted.get(key)][:2]))

    carola_envelope = reference_envelope("carola", carola_chips, carola_months)
    reference_labels = {m["label"] for m in carola_months
                        if m["n"] == m["n_reference"]}
    check("carola: the monthly envelope is built from the reference months only",
          {pair["month"] for pair in carola_envelope["monthly"]["R_B"]["values"]}
          == reference_labels == set(carola_envelope["months"]),
          "%d vs %d" % (len(reference_labels),
                        carola_envelope["n_reference_months"]))
    check("carola: the envelope holds the 16 reference monthly states",
          carola_envelope["n_reference_months"] == 16)
    recomputed = quantile_block([pair["value"] for pair in
                                 carola_envelope["monthly"]["R_B"]["values"]])
    check("carola: q05/q50/q95 are the empirical quantiles of those states",
          all(abs(recomputed[q] - carola_envelope["monthly"]["R_B"][q]) < 1e-12
              for q in ("q05", "q50", "q95")))
    check("carola: q95 is below the reference maximum in both block magnitudes",
          carola_envelope["monthly"]["R_B"]["q95"]
          < carola_envelope["monthly"]["R_B"]["max"]
          and carola_envelope["monthly"]["R_G"]["q95"]
          < carola_envelope["monthly"]["R_G"]["max"])

    check("carola: R_B is the RMS of gamma2 and P only",
          all(abs(c["R_B"] - rms([c["z"]["gamma2"], c["z"]["P"]])) < 1e-12
              for c in carola_chips)
          and tuple(BLOCKS["behaviour"]) == ("gamma2", "P"))
    check("carola: R_G is the RMS of A, D, F and S only",
          all(abs(c["R_G"] - rms([c["z"][dim] for dim in ("A", "D", "F", "S")]))
              < 1e-12 for c in carola_chips)
          and tuple(BLOCKS["geometry"]) == ("A", "D", "F", "S"))

    carola_outside = outside_contrast("carola", carola_chips, carola_months,
                                      carola_envelope, baseline)
    check("carola: the outside block holds 63 chips over 9 monthly intervals",
          carola_outside["n_outside"] == 63
          and carola_outside["n_intervals"] == 9
          and carola_outside["intervals"][-1] == "2024-09 (pre)",
          "%d / %d" % (carola_outside["n_outside"],
                       carola_outside["n_intervals"]))
    check("carola: every outside cell reports r, z, its interval and a position",
          all(set(("r", "z", "ci", "reference_comparison")) <= set(cell)
              for dim, cell in carola_outside["dims"].items() if dim in DIMS)
          and all(set(("value", "delta", "ci", "reference_comparison"))
                  <= set(cell) for dim, cell in carola_outside["dims"].items()
                  if dim in ("R_B", "R_G")))
    check("carola: r over the reference scale reproduces every z",
          all(abs(cell["z"] - cell["r"] / cell["reference"]["scale"]) < 1e-9
              for dim, cell in carola_outside["dims"].items() if dim in DIMS))
    check("carola: every persistence threshold is its monthly reference q95",
          all(abs(carola_outside["persistence"][name]["threshold"]
                  - carola_envelope["monthly"][name]["q95"]) < 1e-12
              for name in ("R_B", "R_G")))
    check("carola: 'sustained' is exactly a run of the fixed minimum length",
          all(carola_outside["persistence"][name]["sustained"]
              == bool(carola_outside["persistence"][name]["longest_run"]
                      and carola_outside["persistence"][name]["longest_run"]["length"]
                      >= PERSIST_MIN) for name in ("R_B", "R_G")))

    carola_summer = summer_block("carola", carola_chips, carola_months,
                                 carola_envelope, carola_outside)
    check("carola: the Summer-2024 subset is post-hoc, descriptive and untested",
          carola_summer["applies"] and carola_summer["post_hoc"]
          and carola_summer["descriptive_only"]
          and not carola_summer["inferential"]
          and 0 < carola_summer["n_chips"] < len(carola_members["outside"])
          and carola_summer["n_intervals"] == 3
          and not any(key in carola_summer
                      for key in ("p_value", "significance", "verdict")))
    check("carola: the summer subset lands inside the held-out months",
          all("2024-06" <= label[:7] <= "2024-08"
              for label in carola_summer["intervals"]))

    carola_environment = environment_context("carola", carola_months,
                                             carola_outside)
    carola_comparison = class_comparison("carola", carola_chips, carola_months,
                                         carola_envelope)
    carola_assessment = assessment("carola", carola_chips, carola_months,
                                   carola_envelope, carola_outside,
                                   carola_comparison, carola_summer,
                                   carola_environment)
    check("carola: the verdict is one of the allowed words",
          carola_assessment["verdict"] in ALLOWED_VERDICTS,
          carola_assessment["verdict"])
    check("carola: the verdict reads R_B and R_G only, never the summer subset",
          carola_assessment["verdict_inputs"] == ["R_B", "R_G"]
          and carola_assessment["summer_is_descriptive_only"]
          and carola_assessment["verdict"]
          == verdict_word(carola_outside["persistence"]["R_B"],
                          carola_outside["persistence"]["R_G"]))
    check("carola: the class comparison covers all eight quantities",
          carola_comparison["order"] == list(COMPARE_ORDER)
          and all(set(carola_comparison["classes"][name]["values"])
                  == set(COMPARE_ORDER)
                  for name in ("reference", "outside", "post")))

    first_report = build("carola")
    second_report = build("carola")
    check("carola: two builds give the same JSON",
          json.dumps(round_floats(first_report), sort_keys=True)
          == json.dumps(round_floats(second_report), sort_keys=True))
    check("carola: two builds give the same report",
          report_md(first_report) == report_md(second_report))

    # --- a window selected by rule, and the audit's fourth bucket -----------
    registry_window = SITES["morandi"]["reference_window"]
    check("window: the accessor reads the registry when nothing is overridden",
          reference_window("morandi") == registry_window
          and window_rule("morandi") is None)
    morandi_chips, _cut = load_chips("morandi")
    probe = ("2017-01", "2017-12")
    with window_override("morandi", probe, "selftest: first full year"):
        inside = [c for c in morandi_chips if in_reference("morandi", c)]
        held = outside_chips("morandi", morandi_chips)
        before = [c for c in morandi_chips
                  if chip_class("morandi", c) == "prefix"]
        check("window: an overridden window moves the set that standardises",
              bool(inside) and all(within(c["day"], probe) for c in inside),
              "%d inside" % len(inside))
        check("outside: no held-out chip sits before the frozen window",
              bool(held)
              and all(month_text(month_key(c["day"])) >= probe[0]
                      for c in held), "%d held out" % len(held))
        check("prefix: what precedes the window is its own bucket",
              bool(before)
              and all(month_text(month_key(c["day"])) < probe[0]
                      for c in before), "%d before it" % len(before))
        check("classes: the four buckets partition the kept chips",
              len(inside) + len(held) + len(before)
              + sum(1 for c in morandi_chips
                    if chip_class("morandi", c) in ("post", "unassigned"))
              == len(morandi_chips), "%d chips" % len(morandi_chips))
        check("window: the rule travels with the window as a sentence",
              window_rule("morandi") == "selftest: first full year")
        probe_report = build("morandi")
        audit = probe_report["classes"]["reconciliation"]
        check("build: a rule-selected window keeps every class in the audit",
              audit["kept"] == len(morandi_chips)
              and audit["prefix"] == len(before)
              and probe_report["reference"]["choice_rule"]
              == "selftest: first full year"
              and probe_report["reference"]["choice"].startswith(
                  "the window the stated rule selected"),
              str(audit))
    check("window: leaving the block restores the registry",
          reference_window("morandi") == registry_window
          and window_rule("morandi") is None
          and outside_chips("morandi", morandi_chips)
          == [c for c in morandi_chips
              if chip_class("morandi", c) == "outside"])
    plain = build("morandi")
    without_tail = [c for c in morandi_chips
                    if chip_class("morandi", c) != "post"]
    trimmed = build("morandi", chips=without_tail)
    # Deleting the tail relabels the divided event month (nothing is left on its
    # post side), so the comparison is the numbers: the frozen census and every
    # held-out departure must be unchanged, byte for byte.
    frozen_keys = ("dims", "persistence", "directional")
    check("build: deleting the tail moves no frozen or held-out number",
          trimmed["post"] is None
          and trimmed["reference"]["n_chips"] == plain["reference"]["n_chips"]
          and trimmed["reference"]["window"] == plain["reference"]["window"]
          and trimmed["reference"]["span"] == plain["reference"]["span"]
          and trimmed["cut"]["rows"] == len(without_tail)
          and json.dumps(round_floats({k: trimmed["outside"][k]
                                       for k in frozen_keys}),
                         sort_keys=True)
          == json.dumps(round_floats({k: plain["outside"][k]
                                      for k in frozen_keys}),
                        sort_keys=True))

    width = max(len(name) for name, _ok, _detail in checks)
    failed = 0
    for name, ok, detail in checks:
        if not ok:
            failed += 1
        line = "%s  %-*s" % ("PASS" if ok else "FAIL", width, name)
        if detail and not ok:
            line += "  " + detail
        print(line)
    print("\n%d checks, %d failed" % (len(checks), failed))
    return 1 if failed else 0


def publish(site, paths=None, draw_figures=True):
    """Write the JSON, the report and (unless asked not to) the plate."""
    report = build(site)
    paths = dict(paths or out_paths(site))
    for key in ("json", "md"):
        os.makedirs(os.path.dirname(paths[key]), exist_ok=True)
    with open(paths["json"], "w") as handle:
        json.dump(round_floats(report), handle, indent=1, sort_keys=True)
        handle.write("\n")
    with open(paths["md"], "w") as handle:
        handle.write(report_md(report))
    if draw_figures:
        os.makedirs(os.path.dirname(paths["jpg"]), exist_ok=True)
        draw(report, paths["jpg"])
    return report, paths


def verify(site):
    """Re-derive everything and compare it with the committed JSON and report."""
    report = build(site)
    paths = out_paths(site)
    problems = []
    if not os.path.exists(paths["json"]):
        problems.append("%s: missing — run the script first"
                        % os.path.relpath(paths["json"], REPO))
    else:
        with open(paths["json"]) as handle:
            committed = json.load(handle)
        difference = first_difference(committed, round_floats(report))
        if difference:
            problems.append("%s: %s" % (os.path.relpath(paths["json"], REPO),
                                        difference))
    if not os.path.exists(paths["md"]):
        problems.append("%s: missing — run the script first"
                        % os.path.relpath(paths["md"], REPO))
    else:
        with open(paths["md"]) as handle:
            committed_text = handle.read()
        if committed_text != report_md(report):
            problems.append("%s: the report does not match the table"
                            % os.path.relpath(paths["md"], REPO))
    # The publication plate must reproduce deterministically too: two draws of
    # the same report are byte-identical, and the committed file is that draw.
    if not os.path.exists(paths["jpg"]):
        problems.append("%s: missing — run the script first"
                        % os.path.relpath(paths["jpg"], REPO))
    else:
        digests = []
        for _ in range(2):
            handle, temporary = tempfile.mkstemp(suffix=".jpg")
            os.close(handle)
            try:
                draw(report, temporary)
                digests.append(sha256_of(temporary))
            finally:
                os.unlink(temporary)
        if len(set(digests)) != 1:
            problems.append("%s: the plate is not deterministic (two draws differ)"
                            % os.path.relpath(paths["jpg"], REPO))
        else:
            committed = sha256_of(paths["jpg"])
            if committed != digests[0]:
                problems.append("%s: the committed plate is not the one this "
                                "report draws (%s vs %s)"
                                % (os.path.relpath(paths["jpg"], REPO),
                                   committed[:16], digests[0][:16]))
    return problems


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="ocv_trajectory.py",
        description="The standard OCV trajectory plate of one case site "
                    "(support-sized chip, per-block standardisation, the "
                    "post-event tail drawn).")
    parser.add_argument("--site", choices=sorted(SITES),
                        help="the site to draw")
    parser.add_argument("--all", action="store_true",
                        help="every registered site, in registry order")
    parser.add_argument("--json", default=None, help="override the JSON path")
    parser.add_argument("--md", default=None, help="override the report path")
    parser.add_argument("--jpg", default=None, help="override the figure path")
    parser.add_argument("--no-figures", action="store_true",
                        help="numbers and report only")
    parser.add_argument("--selftest", action="store_true",
                        help="check the template's rules, write nothing")
    parser.add_argument("--verify", action="store_true",
                        help="re-derive and compare with the committed files")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)
    if args.selftest:
        return selftest()
    sites = sorted(SITES) if args.all else ([args.site] if args.site else [])
    if not sites:
        parser.error("give --site <name>, --all or --selftest")
    failures = 0
    for site in sites:
        paths = out_paths(site)
        if args.json:
            paths["json"] = args.json
        if args.md:
            paths["md"] = args.md
        if args.jpg:
            paths["jpg"] = args.jpg
        if args.verify:
            problems = verify(site)
            if problems:
                failures += 1
                for problem in problems:
                    print("FAIL  %s" % problem)
            elif not args.quiet:
                print("OK    %s: the JSON and the report match the table (%s)"
                      % (site, report_brief(build(site))))
            continue
        report, written = publish(site, paths,
                                  draw_figures=not args.no_figures)
        if not args.quiet:
            print("wrote %s (+ %s%s)"
                  % (os.path.relpath(written["json"], REPO),
                     os.path.relpath(written["md"], REPO),
                     ", " + os.path.relpath(written["jpg"], REPO)
                     if not args.no_figures else ""))
            print("      %s" % report_brief(report))
    return 1 if failures else 0


def report_brief(report):
    """One line: the window, the tail and the contrast verdicts, for the terminal."""
    ref = report["reference"]
    post = report["post"]
    text = "%d chips · reference %s %d (%s → %s)" % (
        report["cut"]["kept"], ref["state"], ref["n_chips"],
        ref["window"][0], ref["window"][1])
    if post:
        text += " · tail %s %d" % ("/".join(post["states"]), post["n_chips"])
    con = report["contrast"]
    if con.get("applies"):
        resolved = [key for key, cell in con["dims"].items()
                    if cell["excludes_zero"]]
        text += " · resolves %s" % (", ".join(resolved) if resolved else "nothing")
    else:
        text += " · no contrast (event-free)"
    assessment = report.get("assessment", {})
    if assessment.get("applies"):
        text += " · held-out: %s" % assessment["verdict"]
    return text


if __name__ == "__main__":
    sys.exit(main())
