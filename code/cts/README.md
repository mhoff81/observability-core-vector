# `code/cts/` — Champlain Towers South site package: analysis and plotting scripts

Five Python modules, no package, no dependencies beyond the standard library plus
numpy / (scipy only for the p-values; without it the scripts still run, just with
fewer p-values). This is the **Champlain Towers South** (Surfside, Miami-Dade,
FL, 25.87306 N 80.12083 W) half of the repository; the LUMO tower lives in
`../lumo/`, Bautzen in `../bautzen/`, KDLO in `../kdlo/`, Carola in `../carola/`
and Morandi in `../morandi/`.

The site is the 12-storey reinforced-concrete condominium that partially
collapsed in the early hours of `2021-06-24` (01:22 EDT = 05:22 UTC); the
standing remnant was brought down on `2021-07-05`. The export observes the site
along **four Sentinel-1 footprints of the same orbit** (ASC rel-48 IW3), one
80 x 80 complex chip per acquisition per track — see `data/cts/cts_tracks.txt`:

```
cts_asc_48_iw3     the tower itself            role target
ctn_asc_48_iw3     Champlain Towers North      role control_ctn   (165 m N)
cte_asc_48_iw3     Champlain Towers East       role control_cte   (~90 m N)
beach_asc_48_iw3   the beach east of the site  role reference_beach
```

The package analyses the two states of the *same* structure around the collapse
date:

```
pre-collapse (healthy)   <   2021-06-24   <=   post-collapse
```

The four footprints share **identical strata, dates and counts** (178 dates
each, 712 raw rows, 543 pre / 168 post echo masks), so this is the one site in
the repository that can do what a single-structure site cannot: a
**difference-in-differences** (DiD) that contrasts the *target* tower against a
*control* tower which sees the same atmosphere and the same orbit. That layer is
first class here (`cts_ocv_did.py`) and is pinned against the three committed
`data/cts/reference/cts_reference_did_*.json` files. The post state is **bare
sand** (the cleared lot): its intensity is moisture-driven, so a change in the
cleared-lot window is *not* evidence about the collapse.

The observability core vector is `OCV = [gamma2, P, D]`; the full mask vector is
`x = [gamma2, P, D, A, F, S]`, carried in **two layers**:

* **echo mask** — the whole 80 x 80 chip (the layer of the 6D vector), the
  project-wide deck echo rule `intensity >= 0.30 x peak` **and**
  `intensity >= 5 x np.median(intensity)`, `A >= 2`;
* **per-track deck mask** — the very same rule on the *track's own* chip, i.e.
  the same six dimensions restricted to `cts` / `ctn` / `cte` / `beach`.

The mask layer rests on a **committed cache**
(`data/cts/cts_windows_mask_cache.txt`): the per-date `cts_windows_*.npz` cubes
do not belong in the repository, so the cache carries the echo mask of every
chip plus the sufficient statistics, and the build re-verifies the rule on it.

Neither layer's `A` is the site's own `coherence_masked_pixels` column: that
column is a *fixed 640-px (10 %) quantile mask* while this package recomputes
the echo mask from the payloads. The CSV keeps the site's fixed 640 as
provenance (`cts_tracks.txt` / the pipeline meta), beside the recomputed
`A` / `gamma2_mask`.

All paths below are relative to the **repository root** and the scripts are
always started from there: inputs in `data/cts/`, figures in `figures/cts/`.

### `cts_ocv_stats.py` — statistics layer (verbatim port)

`_stats`, `welch_mw_test`, `cliffs_delta` (+ `cliffs_delta_bruteforce`),
`bootstrap_delta_ci`, `delta_block`, `sign_test`, `sds_from`, `by_state`, `vals`,
`season_of`, `strata_pair`, `fit_lda`, `predict_lda`, `loo_cv`,
`permutation_test`, `fit_linear_score`, `score_row`, `oof_scores`, `delta_sds`,
`permutation_p`, `paired_bootstrap_diff`, plus the port of the site's own
correlation/regression battery (`spearman`, `mannwhitney`, `pair_block`,
`ols_fit`, `f_test`). Copied line by line from the LUMO project analysis scripts
(via `../carola/`, `../kdlo/`, `../morandi/`), so a difference in the CTS numbers
cannot come from a difference in the statistics.

### `cts_ocv_core.py` — data layer, constants and the reference join

Paths, the analysis channels, the site constants (`SITE`, `COLLAPSE_DATE`,
`EVENT`), the four tracks and their roles, the DiD contrasts, the loaders of the
committed tables (`load_measurements`, `load_tracks`, `load_cache`,
`load_manifest`, `load_csv`, `load_csv_meta`, `load_reference`) and the small
aggregate helpers the generators share (`row_brief`, `summ`, `by_girder`,
`by_track`, `summary_by_state`, `state_of`, `echo_mode_of`, `track_index_of`).

### `cts_ocv_masks.py` — the two mask layers from the committed window cache

Owns the echo mask (port of `analyze_carola_coherence.coherence_gamma2` and
`carola_segment_c_mask_gate`): `largest_component`, `echo_mask`,
`gamma2_of`, `vector_from_row` (re-verifies the rule on every committed mask),
`add_persistence` (`P` = mean pixel frequency inside the own mask), the cache
I/O (`load_cache` / `load_manifest`) and the one-time `--extract` that writes the
cache from the `cts_windows_*.npz` cubes (records every source file's sha256 in
the manifest).

```bash
python3 code/cts/cts_ocv_masks.py --extract   # once, needs the npz cubes
python3 code/cts/cts_ocv_masks.py --check     # acceptance test, no payload
python3 code/cts/cts_ocv_masks.py --json      # print the derived 6D vectors
```

### `cts_ocv_channels_csv.py` — the channel table

Builds `data/cts/cts_ocv_channels.csv` (one row per chip) and
`cts_ocv_channels_meta.json` from the committed cache plus the measurement
extract and the track table, recomputing the six mask dimensions with
`vector_from_row` and deriving `segment`, `track_index`, `state`,
`pre_collapse`, `day`, `month`, `year`, `season`, `echo_mode`.

```bash
python3 code/cts/cts_ocv_channels_csv.py               # build + verify
python3 code/cts/cts_ocv_channels_csv.py --verify-only # verify only
```

### `cts_ocv_did.py` — the difference-in-differences layer

The layer the CTS site makes possible and a single-structure site cannot have.
A target and a control tower see the *same* atmosphere, orbit drift and any
processing effect, so modelling

```
y ~ 1 + T + P + T:P        T = is_target, P = is_post
DiD = coefficient of T:P
```

removes exactly those common terms. The module is a verbatim port of the site's
own generator (`asset_analysis.reference_did` — same `cells` /
`analyse_channel` / `run` / `report`) onto the CTS measurement table, so the
three committed references can be reproduced and pinned:

```
data/cts/reference/cts_reference_did_ctn.json               target vs. CTN
data/cts/reference/cts_reference_did_cte.json               target vs. CTE
data/cts/reference/cts_reference_did_placebo_ctn_vs_cte.json  placebo CTN vs. CTE
```

The placebo (CTN vs. CTE, both controls) is the *specificity* check: neither
control saw an event, so a significant placebo DiD would mean the contrast is
not specific to the collapse. Three things are reported beside the DiD, because
a DiD without them is not interpretable: **pre-period comparability** (Welch,
target vs. reference before the event), **the reference's own step** (a p-value
for whether the reference *also* steps at the event date) and **parallel
pre-trends** (Spearman rho of each series against time). ASCENDING and
DESCENDING are never pooled (EQS-9); here every footprint is ASCENDING, so the
DESCENDING rows are empty.

Nine channels are reported in the committed order — `phase_coherence`,
`phase_snr_db`, `phase_scatterer_count`, `displacement_los_m`,
`coherence_gamma2`, `coherence`, `coherence_masked_pixels`, `intensity`,
`brightness_ratio`. One of them is **degenerate**: `coherence_masked_pixels` is
the site's *fixed* 640-px quantile mask, constant for every row, so its OLS
interaction is pure floating-point noise (`-5.07e-14` in every contrast). It is
pinned structurally (`did` / `n` / `cells` / deltas / `r2`), **never** on
`did_t` / `did_p`.

```bash
python3 code/cts/cts_ocv_did.py --verify            # pin the three refs (series + reports)
python3 code/cts/cts_ocv_did.py --verify --json     # ... machine-readable
python3 code/cts/cts_ocv_did.py --run did_ctn       # print one DiD document (JSON)
python3 code/cts/cts_ocv_did.py --run did_ctn --report  # ... its markdown report
```

`--verify` reproduces every numeric cell of the three committed reference JSONs
to 1e-9 **and** re-renders the three committed markdown reports and compares them
cell by cell (skipping the degenerate channel's t/p), then prints
`PIN PASS`. A single deviation exits non-zero.

### `fig_cts_ocv_paper.py` — the paper figures A–F

The plotting end of the package (run by `figures/cts/fig_cts_ocv_paper.sh` from
the repository root, so its defaults resolve). It reads only
`data/cts/cts_ocv_channels.csv` (plus the cache manifest and the three reference
JSONs) and writes the six PNGs to `figures/cts/`, the English report
`figures/cts/fig_cts_ocv_paper.md` and the JSON artifact
`data/cts/fig_cts_ocv_paper.json`. Panels **A–E** mirror the other sites
(per-state distributions, in-sample effect sizes, learnability, the OCV plane,
the confounder panel); panel **F** is what only a DiD site has — the interaction
t per channel for both controls and the CTN-vs-CTE placebo, the differencing
assumption and the pre-trend parallelism. The script re-derives the three
committed references and its `pin_block` fails the whole run (non-zero exit) on
any deviation, so a successful run certifies the DiD layer too.

```bash
python3 code/cts/fig_cts_ocv_paper.py            # figures A–F + report (~15 min, permutation null)
python3 code/cts/fig_cts_ocv_paper.py --quick    # ~50 s, no permutation null
python3 code/cts/fig_cts_ocv_paper.py --figdir figures/cts \
        --json data/cts/fig_cts_ocv_paper.json   # explicit outputs
```

### `fig_cts_ocv_months.py` — the echo-mask dimensions from the first images to the month after the collapse

The time-axis **companion** of `fig_cts_ocv_paper.py` (like
`../carola/fig_carola_ocv_months.py` and `../morandi/fig_morandi_ocv_months.py`
are companions of their paper scripts), run by
`figures/cts/fig_cts_ocv_months.sh` from the repository root. It reads only
`data/cts/cts_ocv_channels.csv` and keeps the time axis instead of the two pooled
states: the **first Sentinel-1 images available for the footprint** (`F0`, the
2015-09 record start), the three one-month windows before the collapse and the
**first month after** it (`2021-06-24` opens `M+1`). It adds no number to the
A–F pin contract and touches none of it.

Why the baseline is 2015 and not 2014: Sentinel-1 data for **this** footprint do
not exist in 2014 — Sentinel-1A launched on 2014-04-03, and the suitable orbit
(ASC rel-48 IW3) has its **first acquisition on 2015-09-21**
(`data/cts/cts_acquisition_census.md`; the Step-0 census already queried from
2014-10-01 and found nothing usable earlier). `F0` is therefore the record-start
epoch that carries the first images that actually exist, and the report documents
that rather than inventing a 2014 window.

The figure is 2 x 3: (a)–(d) one jittered strip per window for `gamma2`, `P`,
`S`, `D` with the window median, the IQR and the collapse cut as a dashed line
(plus the dotted record gap between `F0` and `M-3`); (e) the `(P, D)` plane, one
dot per chip coloured by window (symlog y); (f) the 6D fingerprint — each window's
median divided by the `F0` (first-images) baseline median on a log axis. It writes
the PNG to `figures/cts/`, the JSON to `data/cts/fig_cts_ocv_months.json` and the
report to `figures/cts/fig_cts_ocv_months.md`. It is deterministic (fixed seeds,
no timestamps), so two runs are byte-identical.

```bash
python3 code/cts/fig_cts_ocv_months.py            # ~10 s
bash figures/cts/fig_cts_ocv_months.sh            # the same, via the wrapper
python3 code/cts/fig_cts_ocv_months.py --figdir /tmp/x --json /tmp/y.json --md /tmp/z.md
```

## Reading the contrast

None of the reported DiD effects is significant at the 5 % level, and the
channels whose *reference also steps* (e.g. `phase_coherence`,
`coherence_gamma2`, `intensity`, `brightness_ratio`) have small reference-step
p-values, so their DiD is not a DiD result — the differencing assumption failed
for those channels. `intensity`'s `pre comparable p = 0.000` says the target and
reference were already on different scales before the event, which is exactly the
post-state bare-sand caveat the site carries. The package therefore reads the
CTS record as an **honest null** on the collapse signal in this window; the DiD
values stand as reported contrasts with their comparability and pre-trend
co-variates attached, not as evidence of damage.

## Verification and determinism

Every module is offline and deterministic (no timestamps in any output, no
network access). The three build/verify commands above are self-checking:

```
python3 code/cts/cts_ocv_masks.py --check           # mask rule on every cache row
python3 code/cts/cts_ocv_channels_csv.py            # re-derive every CSV cell
python3 code/cts/cts_ocv_did.py --verify            # pin the three DiD references
```

The module and site names keep their `cts_` prefix although the folder already
says `cts/`: the import statements (`import cts_ocv_core as core`), the reference
file names and every command line in the reports stay parallel to the other
sites, which is what makes the packages readable as one diff.

