# `code/espoo/` — Espoo site package: analysis and plotting scripts

Five Python modules, no package, no dependencies beyond the standard library plus
numpy / matplotlib / scipy (scipy is only used for the p-values; without it the
scripts still run, just with fewer p-values). This is the Espoo Kurttila mast
(Espoo, Finland) half of the repository; the LUMO tower lives in `../lumo/`,
Bautzen in `../bautzen/`, KDLO in `../kdlo/`, Carola in `../carola/`, Morandi in
`../morandi/`, CTS in `../cts/` and YWF in `../ywf/`.

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
asks whether an OCV separates four *damage* states, Morandi, CTS, Carola and YWF
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


### `fig_espoo_ocv_paper.py` — figures A–E, JSON result, the report and the pins

```bash
python3 code/espoo/fig_espoo_ocv_paper.py                   # recompute everything (~4.5 min)
python3 code/espoo/fig_espoo_ocv_paper.py --quick           # no permutation nulls
python3 code/espoo/fig_espoo_ocv_paper.py --write-reference # (once) define the lock
python3 code/espoo/fig_espoo_ocv_paper.py --no-figures      # JSON + report + pins only
bash   figures/espoo/fig_espoo_ocv_paper.sh                 # the same, from the root
```

Panels: **A** raw distributions of the two measured channels per orbit (+ the
phase ladder), **B** in-sample effect sizes (Cliff's delta with a bootstrap CI,
SDS, the season strata and the two weather controls), **C** the 2-class LDA with
leave-one-out and the label-permutation null (chance 0.5, the two-state
counterpart of LUMO's four-state 0.25), **D** the headline pair across the four
feature sets plus the out-of-time split, **E** the no-event controls — E1 the 14
committed stationarity series, E2 the month/season composition, E3 the same effect
size recomputed *inside* each wind stratum, E4 the phase-ladder availability.

Every number it draws is recomputed from the channel table and then **pinned
field by field** against the three committed reference JSONs in
`data/espoo/reference/` — `espoo_mast_observability.json`,
`espoo_phase_stationarity.json`, `fig_espoo_channels.json` — plus the generator's
own `output_sha256`. `--write-reference` defines the package's **own** lock
`espoo_ocv_findings.json` once; every later run re-derives it and pins it too. Any
deviation aborts with exit code ≠ 0, and the claim of the package is not "we
assert" but **"we reproduce and read off"**.

### `fig_espoo_ocv_months.py` — the record month by month (the companion)

The **standalone companion** of the A–E figures, in the shape of
`../carola/fig_carola_ocv_months.py` and `../cts/fig_cts_ocv_months.py`. Where A–E
pool the whole record, this one keeps the *time axis*: the record is event-free,
so nothing is cut, anchored or windowed and all 150 acquisitions (25 months at
exactly 6 each) appear. Its six panels are (a) `gamma^2` per month, (b) `A` per
month — both with per-acquisition dots and the monthly median per orbit —
(c) the per-month Cliff's delta with a 2000-draw bootstrap CI, quoted only for the
months where both geometries carry ≥ 3 acquisitions (20 of 25), (d) the two
primary stationarity series of `3-8 Hz · strip400_pipeline`, (e) coverage and the
phase ladder, (f) the monthly weather medians.

It reads only `data/espoo/espoo_ocv_channels.csv` and
`data/espoo/reference/espoo_phase_stationarity.json`, adds **no** number to the
A–E pin contract and touches none of it; its own pin block re-derives the two
primary series (`ASC` 0.2474 / `DESC` 0.2165 as literal claims to 4 decimals, the
25 committed monthly means of both, the trend signs), the per-orbit pooled
medians, the 25 × 6 cadence and the phase ladder, and exits non-zero on any
deviation.

```bash
python3 code/espoo/fig_espoo_ocv_months.py            # ~4 s
bash figures/espoo/fig_espoo_ocv_months.sh            # the same, via the wrapper
python3 code/espoo/fig_espoo_ocv_months.py --figdir /tmp/x --json /tmp/y.json --md /tmp/z.md
```

