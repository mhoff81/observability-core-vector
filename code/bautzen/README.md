# `code/bautzen/` — Bautzen OpenLabs bridge package: analysis and plotting scripts

The bridge counterpart of [`../lumo/`](../lumo/README.md): same file names
(`…_stats.py`, `…_masks.py`, `…_core.py`, `…_channels_csv.py`,
`fig_…_ocv_paper.py`, plus the companion `fig_…_ocv_amp_phase.py`), same
JSON/Markdown/CSV shape, same determinism rules —
**only the mask layer is swapped**. The site is the bridge of the Bautzen
OpenLabs test field; the intervention is a retrofitted structural health
monitoring (SHM) installation, documented as the states RS / DS1 / DS2.

Paths below are relative to the **repository root**; the scripts are started from
there (the wrapper `../../figures/bautzen/fig_bautzen_ocv_paper.sh` does that):
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

### `fig_bautzen_ocv_paper.py` — the paper figure, JSON result, report, pin block

The plotting counterpart of `../lumo/fig_lumo_ocv_paper.py`, sharing its
structure — recompute every plotted number from
`data/bautzen/bautzen_ocv_channels.csv`, pin it against the two committed
reference JSONs, one figure set, one result JSON, one Markdown report, no
timestamps, exit code != 0 on any pin failure. The two reference JSONs, the
channel table + meta and `core.LABEL_EN` (for the German doc strings) are the
only inputs. Five figures, one per question:

| file | panel content |
| --- | --- |
| `figures/bautzen/fig_bautzen_ocv_paper_A.png` | the six channels one panel each on one shared y axis, per-overpass dots and the series/state median bar, one column group per series |
| `figures/bautzen/fig_bautzen_ocv_paper_B.png` | the deck-edge mask layer: band contrast per series and state against the detect gate, mask visibility (Wilson 95% CI), deck-edge pixel count, and the in-series gamma^2 effect size per state |
| `figures/bautzen/fig_bautzen_ocv_paper_C.png` | the two cut dates against cuts shifted by +/-2..12 weeks (placebo sweeps), for band contrast and for gamma^2 |
| `figures/bautzen/fig_bautzen_ocv_paper_D.png` | bridge minus control, date by date: the paired per-date difference with the RS / DS1 / DS2 bands and the step/placebo reading of that difference |
| `figures/bautzen/fig_bautzen_ocv_paper_E.png` | seasonality and environment controls: monthly medians per series, and gamma^2 against air temperature / relative humidity at the overpass instant |

It writes the figures to `figures/bautzen/`, the result JSON to
`data/bautzen/fig_bautzen_ocv_paper.json` and the report to
`figures/bautzen/fig_bautzen_ocv_paper.md`.

```bash
python3 code/bautzen/fig_bautzen_ocv_paper.py            # full run (~6 min)
python3 code/bautzen/fig_bautzen_ocv_paper.py --quick    # placebo sweeps copied
bash figures/bautzen/fig_bautzen_ocv_paper.sh            # the same, via the wrapper
```

Where the site deviates from the LUMO definitions (the deck band is only reached
through the local geometry-anchored gate, and the state and the season are
confounded), the deviation is documented in the module docstring and in the
generated report — the numbers are not silently made to look like LUMO's.

### `fig_bautzen_ocv_amp_phase.py` — amplitude/phase time series of the deck edge

The **standalone companion** of the paper figure: instead of the six OCV channels
it reports the polar decomposition of the coherent deck-edge sum, `A(t)` and
`phi(t)`, **over the deck-edge pixels only** (not the whole chip), plus the
**paired differential phase** `dphi(t)` between the bridge and its control
rectangle. With the complex chip `z` of one overpass and `z_sel = z[mask]` the
deck-edge pixels of `masks.date_mask` (`N` pixels, energy `S2 = sum |z_sel|^2`):

```
A(t)    = |sum z_sel| / sqrt(N * S2)                in [0, 1]    coherence amplitude
phi(t)  = arg(sum z_sel)                            in (-pi, pi] interferometric phase
dphi(t) = wrap(phi_bridge(t) - phi_control(t))      in (-pi, pi] paired differential phase
```

`A(t)^2` is *exactly* the site's coherence channel, so the script cross-checks
`A^2 == gamma2_band_raw` (plus the pixel count `n_masked` and the band contrast)
against every committed row of `data/bautzen/bautzen_ocv_channels.csv` — **515
checks** (412 from the table, 67 for the pairing and its wrap identity, 36 for the
state deltas), exit code != 0 on any failure. The four series (bridge/control ×
ASC/DESC) are carried through the amplitude and phase panels, the states RS/DS1/DS2
are the background bands and the two cut dates the dashed lines; the phase is
summarised **circularly** (circular mean and resultant length `R`), never
arithmetically.

The differential panel works on the **pair**: bridge and control share every
overpass time, so `dphi` is formed date by date without interpolation — ASC 28
pairs (8 control overpasses have an empty deck-edge mask and are left out), DESC
35 pairs. It removes what both rectangles see in common (season, atmosphere,
repeat-pass baseline); the panel is drawn once per geometry, under its bridge
column.

Section 4 closes with the **state-to-state** shift
`delta = wrap(circmean(state) - circmean(RS))` (RS → DS1, RS → DS2), for the paired
`dphi` and for the absolute `phi` of each of the four series as counter-check. It
is deliberately *not* the difference of the series' own state phases — the
circular mean is not linear — and each delta carries a seeded percentile bootstrap
interval (20 000 resamples, fixed seed, so the run stays deterministic). Every one
of those intervals covers nearly the whole circle, so all twelve deltas are
reported as **not resolvable** at 5–23 overpasses per state.

It adds **one** figure and does **not** touch the committed reference JSONs or
the figures A–E pin block (the reference files hold no amplitude or phase). It
writes the figure to `figures/bautzen/fig_bautzen_ocv_amp_phase.png`, the result
JSON to `data/bautzen/fig_bautzen_ocv_amp_phase.json` and the report to
`figures/bautzen/fig_bautzen_ocv_amp_phase.md`:

```bash
python3 code/bautzen/fig_bautzen_ocv_amp_phase.py            # ~2 s
bash figures/bautzen/fig_bautzen_ocv_amp_phase.sh            # the same, via the wrapper
```

Status of the mirror: `bautzen_ocv_stats.py`, `bautzen_ocv_masks.py`,
`bautzen_ocv_core.py`, `bautzen_ocv_channels_csv.py`, `fig_bautzen_ocv_paper.py`
and `fig_bautzen_ocv_amp_phase.py` are complete. The generator verifies **every**
committed per-date mask/state value and all effect sizes of the two reference
JSONs (`VERIFICATION OK`, zero problems, two independent traversals
cross-checked); the figure script reproduces the two reference JSONs pin by pin
(8656 checks, 0 failures) and translates the 56 German doc strings.

