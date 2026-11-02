# `code/espoo/` — Espoo site package: analysis and plotting scripts

Five Python modules, no package, no dependencies beyond the standard library plus
numpy / matplotlib / scipy (scipy is only used for the p-values; without it the
scripts still run, just with fewer p-values). This is the Espoo Kurttila mast
(Espoo, Finland) half of the repository; the LUMO tower lives in `../lumo/`,
Bautzen in `../bautzen/`, KDLO in `../kdlo/`, Carola in `../carola/`, Morandi in
`../morandi/`, CTS in `../cts/` and Haltern in `../haltern/`.

The plotting half needs **matplotlib** (`pip install matplotlib`, 3.11.1 here),
which draws the pdf/eps/jpg of every `fig_*.py`; its own dependencies arrive with
it, and its `numpy>=1.25` floor is already met by the numpy the committed `*.json`
was computed with (2.5.4), so a plain install upgrades nothing. matplotlib is
imported inside the drawing functions, never at module level, so the analysis
half of this package also runs on an interpreter without it.

## What this site is, and what it is for

The asset is the **Espoo Kurttila communications mast** (Kurttila, Espoo; the
committed `espoo_mast_observability.json` carries the asset id
`c4ac1eb8-63e6-481d-997e-4bfe90405d69` and OSM node 2233589804). Its record is
**event-free**:

```
first acquisition   2024-09-05
last acquisition    2026-09-07      150 acquisitions in 25 months, 6 per month
event               none — the window contains no damage of any kind
labels              damage_label / structural_state / condition_label are empty
                    on all 150 rows, and the committed upstream verdict says so
                    itself: "This site has no ground-truth state label"
```

There is therefore no severity axis to regress against and no event to anchor a
window on. This package is the family's **event-free control**: the LUMO package
asks whether an OCV separates four *damage* states, Morandi, CTS, Carola and Haltern
ask the same across a collapse, and Espoo asks what the same machinery reports
when there is demonstrably nothing to find. The only contrast the record supports
is the one between its two **orbit geometries** — and `espoo_ocv_core.state_of()`
is that whole definition:

```
ASCENDING   ascending orbit, afternoon pass    n = 70
DESCENDING  descending orbit, morning pass     n = 80
```

## The observability core vector of this site

`OCV = [gamma2, A]`; both are measured quantities and both are properties of the
*geometry*:

* `gamma2` — the pipeline's interferometric coherence (`coherence_gamma2`);
* `A` = `coherence_masked_pixels` — the size in px of the coherent echo mask the
  pipeline found around the mast.

The record delivers exactly **one** mask-geometry quantity and **no mask raster**:
`peak_intensity`, `sub_aperture_brightness`, `dwell_s`, `scatterers`, `amplitude`,
`phase_rad`, `mast_peak_row/col` and every `peak_*` column are empty on all 150
rows — 56 of the export's 78 columns are. So `DIMS = ["A"]`, `D` (density), `F`
(fragmentation), `S` (centroid ↔ peak shift) and `P` (persistence) are
deliberately **absent**: a generator that produced them would be inventing an
echo-mask layer this record does not have. `A` is **copied, not recomputed**, and
the generator proves that copy the identity row by row against the committed
column *and* against the committed reference statistics.

The phase ladder (`phase_coherence`, `phase_rms_rad`) exists for **80 of 150**
acquisitions (53.3 %, 2024-09 … 2025-10, nothing after) and `phase_snr_db` for
**0**. It is therefore reported as *availability* in panel E4 and §A4 of the
report, never as a third observability channel: a channel missing on half the
record cannot carry a claim about the record.

All paths below are relative to the **repository root** and the scripts are always
started from there (the wrappers `../../figures/espoo/*.sh` do that): inputs in
`data/espoo/`, figures in `figures/espoo/`.

### `espoo_ocv_stats.py` — statistics layer (verbatim port) plus this site's two layers

The first block is a **byte-for-byte copy of `../lumo/lumo_ocv_stats.py`**, itself
a verbatim port of the LUMO project analysis scripts: `_stats`, `welch_mw_test`,
`cliffs_delta` (+ `cliffs_delta_bruteforce`), `bootstrap_delta_ci`, `delta_block`,
`sign_test`, `sds_from`, `by_state`, `vals`, `pairs_block`, Jonckheere–Terpstra
(`_jt_statistic`, `jonckheere`, `monotonicity`), `season_of`, `strata`,
`fit_lda`, `predict_lda`, `loo_cv`, `permutation_test`, `fit_linear_score`,
`score_row`, `oof_scores`, `delta_sds`, `permutation_p`, `paired_bootstrap_diff`.
So a difference in the Espoo numbers cannot come from a difference in the
statistics.

The Espoo difference is naming and scope only — and it is a large one:

* **the two states are orbits, not damage labels.** `STATE_ORDER` is
  `ASCENDING` / `DESCENDING` and the field is `state`; no row of the record
  carries a `damage_label`. `DAMAGED = ["DESCENDING"]` and `SEVERITY` keep their
  LUMO names because the ported helpers expect them, but here they are the *group
  side* and an *index* of a two-level geometry cut, not a claim: the figure script
  reports `monotonicity` as not applicable rather than reading it.
* **one measured mask dimension exists.** `COHERENCE_PEAK_FRAC = 0.30`,
  `COHERENCE_MEDIAN_MULT = 5.0`, `COHERENCE_MIN_MASKED = 2` and
  `PRED_MIN_COHERENCE = 0.15` are carried because the upstream classifier uses
  them, but no mask raster is recomputed here — `A` is the pipeline's own
  `coherence_masked_pixels`.
* **two layers the event sites do not need**, ported line by line so that the pin
  block can reproduce them:
  * `group_stats`, `classify`, `CLASS_TEXT`, `reference_for`, `per_group`,
    `compare_orbit`, `orbit_block` port `espoo_mast_observability.py` — the
    per-pass / per-orbit observability block, its reference-relative class
    (`compact_echo_strong` against the LUMO medians `LUMO_ASC_G2_MEDIAN = 0.0617`
    / `LUMO_DESC_G2_MEDIAN = 0.0261` / `LUMO_HEALTHY_NMASKED_MEDIAN = 19` /
    `LUMO_CLUTTER_FLOOR = 0.04`) and the orbit test;
  * `mann_kendall`, `theil_sen`, `block_bootstrap_slope`, `lag1_autocorr`,
    `runs_test`, `ljung_box`, `annual_harmonic`, `split_half_series`,
    `monthly_means_series`, `stationarity_flags`, `stationarity_of`,
    `phase_series_block` (memoised), `stationarity_summary`, `flagged_groups`,
    `primary_series` port `espoo_phase_stationarity.py` — the Theil–Sen /
    Mann–Kendall stationarity layer. The phase scan reads a *different* input file
    than the channel table, but `espoo_phase_stationarity.json` commits the full
    per-acquisition series of all 14 of its groups under `series_by_group`, so
    every statistical number of that file — median, IQR, Theil–Sen slope and its
    moving-block bootstrap CI, Mann–Kendall p, lag-1 autocorrelation (raw and
    detrended), runs test, Ljung–Box Q, split-half Cliff's delta with its CI, the
    annual harmonic and all 25 monthly means — is recomputed here and pinned
    field by field, together with the flag rule that turns those primitives into
    the committed `verdict`. Only the keys the file does not carry (`dwell_s`,
    `n_sub`, `nyquist_hz`, `n_observable`, `dwell_control`) are not re-derived;
    the report names them.

Conventions are the origin's: `SDS = 10 * |Cliff's delta|`,
`delta(ref, grp) > 0` means the *group* is larger than the reference, bootstrap
seeds 7 (in-sample) resp. 11/12/13 (out-of-fold) and permutation seed 7/8.


All paths below are relative to the **repository root** and the scripts are always
started from there (the wrappers `../../figures/espoo/*.sh` do that): inputs in
`data/espoo/`, figures in `figures/espoo/`.

### `espoo_ocv_core.py` — data layer, constants and loaders

Paths (`DATA`, `REF`, `FIGDIR`, `RAW_CSV_PATH`, `CSV_PATH`, `CSV_META_PATH`,
`REFERENCE_FILES`), the vector definition (`DIMS = ["A"]`, `DIM_LABEL`,
`FEATURES = ["gamma2", "A"]`, `COVARIATES = ["wind_speed_ms", "temperature_c"]`),
the model registry (`MODELS` / `MODEL_LABEL` with `HEADLINE_PAIR = "gamma2A"`,
`BASELINE = "gamma2"`, `SECONDARY_PAIR = "A"`), the phase block (`PHASE_DIMS`), the
identity map `COHERENCE_COLS`, the seeds and budgets (`LAMBDA_GRID`,
`HEADLINE_LAMBDA = 0.5`, `N_PERM_2CLASS = N_PERM_PAIR = 1000`, `N_BOOT = 10000`,
`N_BOOT_STRATA = 2000`, `RNG_SEED = 7`, `OOT_FRACTION = 0.5`) and the loaders
(`load_csv`, `load_reference`, `load_csv_meta`).

`state_of()` is the whole state definition of the package (the orbit), and
`SITE` / `ASSET_ID` / `ASSET_NAME` / `OSM_NODE` / `GROUP_LABEL` carry the site's
identity verbatim from the committed provenance, with `EVENT = None` on purpose —
every figure says so. `row_brief`, `n_by_state`, `feature_dataset` and `sha256`
are the small shared helpers; `python3 code/espoo/espoo_ocv_core.py` prints a
read-only brief of the layer (rows, span, per-orbit split, the generator's checks,
the committed orbit test and verdict, the stationarity summary and the phase
ladder count).

### `espoo_ocv_channels_csv.py` — build `data/espoo/espoo_ocv_channels.csv`

Reads the committed site record table `data/espoo/espoo_channels.csv` (the site's
own export: the same 150 acquisitions that `espoo_mast_observability.json`
commits under `acquisitions`, in the same order, with the same 78 fields) and
writes the extended table plus `espoo_ocv_channels_meta.json`. It does the honest
amount of work:

* **appends** `A`, `doy`, `phase_available` (`A` = `coherence_masked_pixels`, the
  OCV name of the one measured mask quantity) and **fills in place** `month` and
  `acquisition_date`, which the export carries but leaves empty on all 150 rows —
  3 appended + 2 filled, nothing dropped, 150 rows in committed order because the
  seeded bootstrap CIs and permutation tests depend on that order;
* **verifies** the copy against the committed JSON: 4050 copied fields, 235
  reference channel statistics, `A_equals_committed_masked_pixels`, 0 rows left
  empty, 25 months;
* writes the meta file (source and output sha256, `mask_definition`,
  `derived_checks`, `filled_columns`, `new_columns`, the 25-month composition and
  the per-orbit counts) that `espoo_ocv_core` and both figure scripts read back.

```bash
python3 code/espoo/espoo_ocv_channels_csv.py                # build + verify
python3 code/espoo/espoo_ocv_channels_csv.py --verify-only  # verify only
```

The build is byte-idempotent: two runs produce the same CSV and the same meta
JSON.


### `espoo_base_width.py` — the one geometric number the record does not carry

`asset_evidence.json` holds `base_width_m` as `unconfirmed_default` (12.0 m, the
registry default for a lattice tower) with the reason "nobody read a width". The
mast is 36 m tall and documented, the footprint width is not — and it is the only
geometric input missing for a support rectangle (the chip itself is sized by the
200 m reference ring, not by the mast).

```bash
python3 code/espoo/espoo_base_width.py --panos                 # the pano table
python3 code/espoo/espoo_base_width.py --selftest              # the geometry and the guard
python3 code/espoo/espoo_base_width.py \
    --pano AEYr7sSKMbOVLnwPVDRY3g --px 116 --fov 60 \
    --pano SbCRupSUGyADKLwb7veNGQ --px 148 --fov 60 --adopt
```

Measured quantity: **outer edge to outer edge of the two outer legs at ground
level**, at 60 deg horizontal over 640 px (10.67 px/deg, the scale
`asset_evidence.json` already uses for its "12 m = 110 px at 64.5 m" note). The
rule is deliberately unforgiving: **two independent readings within 1.5 m** or
nothing is adopted, and `--adopt` writes the value plus both readings and the
formula into `asset_evidence.json` (keeping `.bak-width`) and refuses to
overwrite a value that is already confirmed.

Why 1.5 m is enough: the width only sets the support rectangle inside the chip,
and the 30 m guard means a 1.5 m error cannot clip the mast. Why it is still
worth doing: it makes `A` admissible instead of defaulted.

Pano table (from `espoo_site_geometry.md`; 12 m ≈ px at 60 deg / 640 px):
`AEYr7sSKMbOVLnwPVDRY3g` 63.1 m → 116 px (best view), `SbCRupSUGyADKLwb7veNGQ`
49.0 m → 149 px, `GveCM2wS-gN-zWrNFmII_w` 49.9 m → 146 px, `lSdEAOCz4_FlRTNG352ziQ`
53.0 m → 138 px, `W9HFZy8cxmXKeMXtRohb3A` 54.1 m → 135 px, `vvhBpOg7YlLxmWTeFPMLbQ`
68.2 m → 107 px, `eggBJwZBZH1iEGBuWAXE3w` 50.2 m → 145 px,
`f2MnW9O3tbDOVvzhfNaBLg` 60.1 m → 122 px, `9c6Glb7byX0fkTMNHBHhqA` 74.0 m → 99 px.

### The echo mask: what the first S3 probe found (2026-10-22)

The support-chip rule **can** be applied at Espoo, but the site is at the edge of
the swath, and that limits which overpasses are cuttable:

```
python3 tools/ocv_migrate.py --probe 3 --site espoo \
    --env-file secrets/espoo_slc.env --auth static

2024-09-05 ASCE   35 x 129  A=65                     cut
2024-09-10 ASCE  REFUSED: window ends at column 23344, past the raster's 23280 samples
2024-09-05 DESC  REFUSED: window ends at column 23062, past the raster's 22998 samples
```

* **IW1 is the only sub-swath that images the mast.** IW2 and IW3 return `NaN`
  for the anchor (the point is outside their grids), so `subswath` in
  `tools/ocv_migrate.py` is pinned to `IW1`.
* **In IW1 the mast sits on the far range edge** of the burst raster: the chip
  needs `half_bins = 64` (200 m ring + margin at 3.70 m ground range), and on
  two of the three probed dates only ~0 of those 64 columns exist beyond the
  anchor. The refusal is the tool doing its job — it will not crop the ring.
* The offset **varies per acquisition** by a few hundred metres (the burst
  window is placed per date), so the eligibility is a property of the date, not
  of the site: 1 of 3 probed overpasses was cuttable.

What this means for a cut: the mask would exist only for the subset of dates
whose burst window leaves the ring margin. That subset is deterministic but not
random — it selects on burst-window placement, which co-varies with orbit
geometry and season — so a cut must report the excluded dates and the acceptance
rate, and the stability figure must say which series it is actually looking at.

Alternatives, in order of rule fidelity: (1) cut the eligible subset and report
it, (2) a **half-ring** (azimuth ±200 m, range 30-200 m on the near side only)
if every date is needed — a documented deviation from the eight-site rule, (3)
no mask, which is where this stood before the probe.

The eight-date probe (2026-10-22) makes the pattern quantitative:

| date | orbit | result |
|---|---|---|
| 2024-09-05 | ASC | **cut**, chip 35 x 129, `A = 65` |
| 2024-09-10 | ASC | refused (column 23344 > raster 23280) |
| 2024-09-17 | ASC | **cut**, `A = 61` |
| 2024-10-04 | ASC | refused (23344 > 23280) |
| 2024-09-05 | DESC | refused (23062 > 22998) |
| 2024-09-10 | DESC | **cut**, `A = 392` |
| 2024-09-17 | DESC | refused (23061 > 22997) |
| 2024-10-04 | DESC | **cut**, `A = 708` |

**4 of 8 cut, and the refusals alternate by acquisition** (same column offsets to the
sample): the burst-window placement flips, so a cut of the eligible subset would
keep every other overpass, not a random half. Two consequences for any later
figure: the series is halved *and* selected on burst geometry, and `A` at the
mast differs by an **order of magnitude between the orbits** (ASC 61-65,
DESC 392-708) - a property of the observation geometry, not of the tower.

### `fig_espoo_ocv_chips.py` — the first OCV vector of the chip, and its three tests

The 75 support chips carry what the pipeline export never had: `A, D, F, S, P`
and `gamma2_mask`. This module runs the same three tests as the pipeline figure — retired with the
twelve-plate scope, its numbers frozen at `data/espoo/fig_espoo_ocv_orbits.json` —
on
that vector — with one structural change: **the pairing is gone**. The cut is
restricted to the two tracks that place the mast inside the swath (ASC =
relative orbit **160**, DESC = **51**), and those tracks have only **three shared
dates**, so the orbit contrast is a two-sample comparison.

```bash
python3 code/espoo/fig_espoo_ocv_chips.py             # ~25 s, JPG + PDF + EPS + MD + JSON
python3 code/espoo/fig_espoo_ocv_chips.py --selftest  # 24 checks, the last twelve draw a plate
python3 code/espoo/fig_espoo_ocv_chips.py --verify    # 42 fields + 133 figure pins
```

The figure is one **170 mm plate**: `6.6929 × 7.5500 in`, six panels of
74.6 × 45.7 mm in reading order `A D / F S / P γ²`, written as JPG (630 dpi),
PDF (Type-42 vector) and EPS (alpha pre-blended over white). The canvas *is*
the print width — `bbox_inches="tight"` is never used, so the file is not
re-cut after the fact and no glyph shrinks; the content has to stay inside the
canvas instead, which four guards assert on every build: the statistics line
inside its axes, the panel tag clear of the DESC row by 1 mm, the legend clear of
the first statistics band, and the drawn content at least 0.3 mm from every edge.
`--verify` measures all three containers again and re-derives the 133 pinned
numbers from the CSV.

The plate carries **keys, not titles**: no headline and no panel titles. Each
cell is tagged with its grid position, `a1` … `c2` (`style.panel_id`), names its
dimension in a one-line y label (`A  [cells]`, `D  [1]`, …), and the two track
counts sit in one legend in the band the headline used to start. The tag takes
`TAG_PAD_PT = 24` instead of the grammar's 2 pt, because this plate's statistics
band already owns the space above its axes — an *inside* tag was measured to cover
a drawn point in panel `c1`, so it is drawn above the band and the guard re-checks
the clearance on every build.

Typography follows `figure_style`: panel tag 8.8, axis label 7.6, ticks and
legend 6.8, statistics 6.6 pt — no title, nothing below 6 pt. The dimensions are
labelled canonically in the y labels and in the pinned `label` of each panel
(`A — observable area (cells)` … `γ² — coherence`); the column names of the table
are not drawn, and the two statistics rows sit in one fixed band above every
panel: DESC at +10.0 pt, ASC at +1.0 pt, each the slope per year, its bootstrap
interval, the Mann–Kendall p and the verdict. The prose definitions — "observable
area", "centroid to peak" — and the definition of the metric are in the caption.

The caption and a machine-readable audit of the plate — target and measured
width, `figsize`, dpi, final point sizes, bounds, EPS alpha policy, one pinned
row per panel and track — are in `figures/espoo/fig_espoo_ocv_chips.md`. The
scientific JSON is untouched by the layout: it stays byte-identical, which is
what `--verify` re-derives field by field.

Result, in three lines:

* **No drift anywhere.** All six dimensions on both tracks get `trend = no`; the
  closest case is `P` on ASC (slope −0.022/a, interval [−0.044, +0.0007], MK
  p = 0.013 — the interval covers zero, so the flag rule says no).
* **One year-on-year difference**: `P` on ASC drops by −0.018 [−0.032, −0.002]
  with 2 of 9 months positive. Everything else repeats within its interval.
* **The two tracks are not two samples of one distribution.** The orbit contrast
  is a complete separation: `A` 76 vs 528, `D` 0.021 vs 0.117, `F` 19 vs 133
  (δ = −1.00), `S` 25.7 vs 14.4, `P` 0.234 vs 0.147, `gamma2_mask` 0.0164 vs
  0.0007 (δ = +0.84 … +0.94), every p < 0.001.

So "Espoo is stable" holds **within a track** and not **between** tracks: the
series is flat over two years, and the two acquisition geometries sit at
completely different levels. And the tautology travels with every number:
`D = A / bbox_area` exactly, `gamma2_mask` has `n = A` in its denominator, and
`A` itself differs sevenfold between the tracks.

### Width sensitivity — run at 11, 12 and 13 m, and it does **not** reach the cells

The 12 m base width is still `unconfirmed_default`, so the chip was re-cut at
**11 m and 13 m** on the same 75 dates (`--base-m` in `tools/ocv_migrate.py`,
outputs in `/tmp/espoo_w11`, `/tmp/espoo_w13`):

| width | support band | chip | support cells | support_m | margin | A median |
|---|---|---|---|---|---|---|
| 11 m | (5, 33) | (35, 129) | 165 | (139.113, 69.540) | (1, 1) | 329 |
| 12 m | (5, 33) | (35, 129) | 165 | (139.113, 69.540) | (1, 1) | 329 |
| 13 m | (5, 33) | (35, 129) | 165 | (139.113, 69.540) | (1, 1) | 329 |

The `(date, orbit, A, gamma2)` series are **byte-identical** at all three widths
(row hash `41fc35be…`), and so is the support geometry — which means every other
dimension is identical too, because `D, F, S, P` are functions of the same chip
and the same masks.

**Why:** the band is quantised in cells, `band = 2·ceil(half) + 1`. One metre of
base width is **0.07 lines** (13.9 m/line) and **0.27 bins** (3.7 m/bin) — below
the rounding threshold. The support is also dominated by the **height layover**
(36 m → 12.88 bins), so a ±1 m change never approaches a cell boundary.

So the pre-registered width test passes **by construction**, and the honest
reading is narrower than "the widths agree": the trend decision cannot differ
between 11 and 13 m at this site because the mask is the same. A confirmed width
inside ±2 m would change nothing in the numbers; it would only make the statement
"the support is the mast base" literally true instead of defaulted. A width would
have to move by about **4 m in range** (one bin) or **14 m in azimuth** (one
line) before it could touch a single cell.
