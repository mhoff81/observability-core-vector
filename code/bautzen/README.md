# `code/bautzen/` — Bautzen OpenLabs bridge package: analysis and plotting scripts

The bridge counterpart of [`../lumo/`](../lumo/README.md): same file names
(`…_stats.py`, `…_masks.py`, `…_core.py`, `…_channels_csv.py`,
plus the figure script `fig_bautzen_ocv_orbits.py`), same
JSON/Markdown/CSV shape, same determinism rules —
**only the mask layer is swapped**. The site is the bridge of the Bautzen
OpenLabs test field; the intervention is a retrofitted structural health
monitoring (SHM) installation, documented as the states RS / DS1 / DS2.

Paths below are relative to the **repository root**; the scripts are started from
there (the wrapper `../../figures/bautzen/fig_bautzen_ocv_orbits.sh` does that):
inputs and outputs in `data/bautzen/`, figures in `figures/bautzen/`.

### `bautzen_ocv_stats.py` — statistics layer (port of the LUMO primitives)

`_stats`, `welch_mw_test`, `cliffs_delta`, `cliffs_delta_bruteforce`,
`_u_from_ranks`, `bootstrap_delta_ci`, `delta_block`, `sign_test`, `sds_from`,
`by_state`, `vals`, `pairs_block` — the same functions, with the site's own
constants: `STATE_ORDER = [RS, DS1, DS2]`, `SEVERITY = {RS: 0, DS1: 1, DS2: 2}`,
`REFERENCE_STATE = "RS"` (the pre-intervention state, LUMO's "healthy"),
`DAMAGED = [DS1, DS2]`, seeds `RNG_SEED = 7`, `BOOT_SEED = 42`, `N_PERM =
20000`, `N_BOOT = 5000`, `BOOT_MAIN = 10000`, `BOOT_STRATA = 2000`,
`PLACEBO_WEEKS = (2, 4, 6, 8, 10, 12)` and the `WINTER_MONTHS` / `SUMMER_MONTHS`
strata.

### `bautzen_ocv_masks.py` — mask layer (the deck edge, the site's only deviation)

The bridge deck is a bright, geometrically known line, so the mask is **anchored
geodetically** instead of by threshold: the deck band sits at the rect centre
(`CENTER_ROW = CENTER_COL = 40`, from the 28.90 m deck axis / 5.90 m deck width
and the 80 x 80 px rect), the search ROI is the band ± `ANCHOR_HALF = 4` rows and
± `COL_HALF = 4` columns, the band itself is the anchor row ± `BAND_HALF = 1`.
A pixel is masked when `|z|^2 >= K_LOCAL = 1.5 x median(neighbour azimuth rows)`;
the house gates are `MIN_N_MASKED = 2` and a visible deck edge from
`CONTRAST_DETECT = 1.5` (band/neighbour contrast). A *global* echo rule of the
LUMO kind is not usable here — the scene clutter pulls the median up and n
collapses to ~0.

`load_intensity`, `load_complex`, `largest_component`, `gamma2` (on the selected
pixels, floor-corrected), `house_mask`, `date_mask`, `analyse_series` (takes
**intensities**, one date per row — see the note in the module),
`add_persistence` (frequency map built **per series**, because the deck mask
differs between bridge and control and between ASC and DESC),
`state_metrics`, `discrimination`, `counts_of`, `unions`.

### `bautzen_ocv_core.py` — data layer and constants

Paths (`CSV_PATH`, `CSV_META_PATH`, `METEO_PATH`, `REFERENCE_FILES`), the four
series (`ASC`/`DESC` bridge and control 500 m east of the bridge, committed
German keys `SERIES` + English `SERIES_EN`), the state rule `state_of()` — assigned
by **overpass time**, not by date, because both state boundaries fall in the
middle of a day (`OVERPASS_LOCAL = {ASC: 18:51, DESC: 07:08}`, `CUT_DS1 =
2025-05-12 10:00`, `CUT_DS2 = 2025-09-29 10:00`, windows in `STATE_WINDOWS`) —
the channel/feature order (`FEATURES = [gamma2, P, D, A, F, S]`), the model sets
of Figure D (`MODELS`, `HEADLINE_PAIR = "gamma2PD"`), the seeds
(`N_PERM_PAIR`, `RNG_SEED`, …) and the loaders: `load_csv()`,
`load_csv_meta()`, `load_reference()`, `load_meteo()`, `pixels_to_mask()` (rebuilds
a per-date mask from the committed `mask_pixels` column), `geom_of()`,
`is_bridge()`, `month_of()`, `by_series()`, `n_by_state()`, `feature_dataset()`,
`sha256()`.

### `bautzen_ocv_channels_csv.py` — generator of the channel table

Enumerates the four committed rect text files (`data/bautzen/bautzen_rects_*.txt`)
into `data/bautzen/bautzen_ocv_channels.csv`: one row per date and series with
`gamma2` and the deck-edge channels `A, D, F, S, P` plus the auditable
ingredients (`anchor_row, peak_row, band_contrast, mask_pixels, …`) and the state
of the date. No burst cache, no network — a fresh clone rebuilds the table
offline. It verifies against the two committed reference JSONs
(`data/bautzen/reference/bautzen_deck_edge_mask.json`,
`…_state.json`) date by date, and cross-checks the two independent traversals
(channel path vs. analysis path) against each other:

```bash
python3 code/bautzen/bautzen_ocv_channels_csv.py                # build + verify
python3 code/bautzen/bautzen_ocv_channels_csv.py --verify-only  # verify only
```

### `fig_bautzen_ocv_orbits.py` — OCV time series per dimension and orbit

A second **standalone companion** to the paper figure: one panel per OCV
dimension and orbit, i.e. a **6 x 2 grid** — rows in the formal OCV order
`gamma2`, `P`, `D`, `F`, `A`, `S`; columns ASC (IW2, ~18:51 local) and DESC
(IW3, ~07:08 local). Each panel is
a time series over the **actual local overpass instant** (`core.acq_dt`), not
midnight, because both state cuts fall in the middle of a day; a date-only axis
would put the boundary-day overpasses in the wrong state. The three state
windows RS / DS1 / DS2 are background bands in different colours, the two cuts
are dashed lines, the line connects the bridge values in acquisition order and
the overpasses are drawn as state-coloured markers (`o` RS, `s` DS1, `^` DS2).
The bridge series are complete (all six channels on all 36 ASC and 35 DESC
overpasses); the plotting path would leave a missing value as a gap, and nothing
is interpolated.

The plate carries **keys, not titles**: each panel is tagged with its position in
the grid by `figure_style.tag_axes(ax, figure_style.panel_id(row, col))` — `a1`
… `f2`, the letter the dimension row and the number the orbit column — there is
no figure title and no column title, and the columns are named once, in the key
line of the legend band under the grid (`columns — 1: ASC, IW2, ~18:51 local ·
2: DESC, IW3, ~07:08 local`, built from `core.OVERPASS_LOCAL`, the same constant
the state rule reads). The height stays 245 mm: the band the removed titles kept
free goes to the panels. What the plate shows is the report's `## Caption`, the
text that goes under the figure.

The figure is **bridge only** (ASC/DESC `Brücke`). The control rectangle is the
counter-check of figures A–E, but its ASC `gamma2` range (-0.94 .. +0.93) is much
wider than the bridge range (-0.14 .. +0.24); a shared panel would compress the
bridge signal this figure is meant to show.

It recomputes every point from `data/bautzen/bautzen_ocv_channels.csv` and pins
it against the committed metadata — 142 rows, the state rule
(`core.state_of(date, orbit) == state`), per-series date/state counts and date
lists, unique dates, `A == n_masked`, channel presence on both orbits, and the
sha256 of every rect file — **41 checks, 0 failures**, exit code != 0 on any
failure. The plate is drawn on the shared print frame of `code/figure_style.py`
— 170.0 mm wide (`COL_WIDTH`), the shared type ladder and pens — and written by
`figure_style.save`, which refuses a layout that overflows that width and
measures all three containers afterwards. Outputs are the 630 dpi JPEG figure
`figures/bautzen/fig_bautzen_ocv_orbits.jpg` with its `.pdf` and `.eps` under the
same stem, the result JSON `data/bautzen/fig_bautzen_ocv_orbits.json` and the
report `figures/bautzen/fig_bautzen_ocv_orbits.md`:

```bash
python3 code/bautzen/fig_bautzen_ocv_orbits.py            # ~3 s
PYTHON=/usr/bin/python3 bash figures/bautzen/fig_bautzen_ocv_orbits.sh
```

Status of the mirror: `bautzen_ocv_stats.py`, `bautzen_ocv_masks.py`,
`bautzen_ocv_core.py`, `bautzen_ocv_channels_csv.py` and
`fig_bautzen_ocv_orbits.py` are complete. The
generator verifies **every**
committed per-date mask/state value and all effect sizes of the two reference
JSONs (`VERIFICATION OK`, zero problems, two independent traversals
cross-checked), and the orbit plate pins every value it draws (8656 checks,
0 failures) over the 56 German doc strings of the reference layer.


