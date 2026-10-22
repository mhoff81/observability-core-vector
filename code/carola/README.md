# `code/carola/` — Carola site package: analysis and plotting scripts

Eight Python modules, no package, no dependencies beyond the standard library plus
numpy / matplotlib / scipy (scipy is only used for the p-values; without it the
scripts still run, just with fewer p-values). This is the Carola Brücke
(Dresden) half of the repository; the LUMO tower lives in `../lumo/`, Bautzen in
`../bautzen/`, KDLO in `../kdlo/` and Morandi in `../morandi/`.

The site is the collapsed **Carolabrücke** (bridge `way` B 170, Wikidata
Q1044279). Its record is a **per-girder chip export**: per acquisition the
pipeline stored one 80 x 80 complex chip per girder request (segment 0/1/2 →
Girder A/B/C) plus 7 x 7 chips for the wrongly typed tower asset. The package
analyses the two states of the *same* structure around the collapse date
`2024-09-11`:

```
pre-collapse (healthy)   <   2024-09-11   <=   post-collapse
```

so season, weather, traffic and pipeline generation all move with the state. The
observability core vector is `OCV = [gamma2, P, D]`; the full mask vector is
`x = [gamma2, P, D, A, F, S]`, carried in **two layers**:

* **echo mask** — the whole 80 x 80 chip (the layer of the 6D vector), the
  project-wide bridge-deck echo rule `intensity >= 0.30 x peak` **and**
  `intensity >= 5 x np.median(intensity)`, `A >= 2`;
* **deck mask** — the very same rule on the *girder's own* chip, i.e. the same
  six dimensions restricted to Girder A / B / C (the layer this package adds).

The mask layer rests on a **committed cache**
(`data/carola/carola_windows_mask_cache.txt`): the 182 MB chip payloads do not
belong in the repository, so the cache carries the echo mask of every chip plus
the sufficient statistics, and the build re-verifies the rule on it.

Neither layer's `A` is the pipeline's own `coherence_masked_pixels` column: that
column was computed on the pre-fix *asset-point* chip while the committed export
holds the per-girder re-extraction, so the two agree only by coincidence
(**60 of 1,713** rows). The site's analysis project therefore recomputed the mask
from the payloads — and this package re-derives it from the committed cache.

All paths below are relative to the **repository root** and the scripts are
always started from there (the wrapper `../../figures/carola/fig_carola_ocv_paper.sh`
does that): inputs in `data/carola/`, figures in `figures/carola/`.

### `carola_ocv_stats.py` — statistics layer (verbatim port)

`_stats`, `welch_mw_test`, `cliffs_delta` (+ `cliffs_delta_bruteforce`),
`bootstrap_delta_ci`, `delta_block`, `sign_test`, `sds_from`, `by_state`, `vals`,
`season_of`, `strata_pair`, `fit_lda`, `predict_lda`, `loo_cv`,
`permutation_test`, `fit_linear_score`, `score_row`, `oof_scores`, `delta_sds`,
`permutation_p`, `paired_bootstrap_diff`, plus the port of the site's own
correlation/regression battery (`spearman`, `mannwhitney`, `pair_block`,
`ols_fit`, `f_test`). Copied line by line from `../kdlo/kdlo_ocv_stats.py`, which
is itself a verbatim port of the LUMO project analysis scripts, so a difference
in the Carola numbers cannot come from a difference in the statistics. The only
Carola edits are naming: the state field is `state` and the two states are the
site's collapse alphabet `pre-collapse (healthy)` / `post-collapse`
(`analyze_carola_coherence.STATE_LABELS`).

### `carola_ocv_core.py` — data layer, constants and loaders

Paths, the analysis channels (`OCV = [gamma2, P, D]`), the 6D feature vector
(`FEATURES = [gamma2, P, D, A, F, S]`, `DIMS`), the axis/title labels, the
weather and control columns (`WEATHER`, `CONTROLS`), the three girders
(`GIRDERS`), the model/seed constants (`MODELS`, `HEADLINE_PAIR = x_6d`,
`LAMBDA_GRID`, `HEADLINE_LAMBDA = 0.5`, `RNG_SEED = 7`, `N_BOOT = 10000`,
`N_PERM = 1000`) and the loaders (`load_csv`, `load_csv_meta`, `load_cache`,
`load_manifest`, `load_measurements`, `load_segments`, `load_reference`).
`MASK_COLUMNS` lists the echo-mask columns the generator appends, so that
`CSV_COLUMNS` (and the CSV header) carry the whole 6D vector.

It also carries the small shared aggregates (`row_brief`, `summ`, `by_girder`,
`summary_by_state`, `day_of`, `state_of`, `echo_mode_of`).

### `carola_ocv_masks.py` — the two mask layers (port of `../kdlo/kdlo_ocv_masks.py`)

The module owns the mask rule of both layers: `echo_mask` (`0.30 x peak`,
`5 x np.median`, `A >= 2`), `largest_component` / `connected_components_8`
(port of the LUMO `_largest`), `gamma2_of` (`|sum z|^2 / (A sum |z|^2)`,
clamped to [0, 1]), `vector_from_row` (`A`, `D = A / bbox_area`, `F`, `S`, the
six dimensions, the bbox/centroid/peak geometry and the masked `gamma2`),
`add_persistence` (`P` = mean pixel frequency of the own mask over all chips of
the majority shape) and the cache I/O (`load_cache`, `load_manifest`, `extract`).

`extract()` is the one-time offline step: it reads the 182 MB
`carola_windows_full.txt`, applies the mask rule to every chip, resolves the
per-girder segment from the committed `carola_windows_index.txt` +
`carola_segments.txt`, drops byte-identical duplicate payloads, and writes the
cache plus its manifest (the source file's sha256 included). After that the
cache is the only mask input. `python3 code/carola/carola_ocv_masks.py --check`
is the acceptance test of the committed cache: every masked pixel must satisfy
both thresholds, the brightest masked pixel must *be* the peak
(`max_masked == peak`), `A >= 2`, and the recomputed `gamma2` must agree with the
stored sufficient statistics (472 duplicate payloads deduplicated, 1,717
80 x 80 rows, no chip without an echo mask).

### `carola_ocv_channels_csv.py` — build `data/carola/carola_ocv_channels.csv`

Reads the committed cache and the committed measurement extract
(`carola_measurements_full.txt`), recomputes the six mask dimensions of every
chip, derives `segment`, `girder_index`, `state`, `pre_collapse`, `day`, `month`,
`year`, `season` and `echo_mode`, and writes the channel table (1,717 rows,
57 columns) plus `carola_ocv_channels_meta.json` (file sha256s, the mask rule,
the guard counts, the derived-column definitions, the `coherence_masked_pixels`
cross-check and the per-state summary of every channel). The build **refuses to
write** if a chip has no echo mask or if the mask layer cannot be verified.

```bash
python3 code/carola/carola_ocv_channels_csv.py                 # build + verify
python3 code/carola/carola_ocv_channels_csv.py --verify-only   # verify only
# one-time re-extraction of the mask cache (needs the 182 MB payload file):
python3 code/carola/carola_ocv_masks.py --extract
python3 code/carola/carola_ocv_masks.py --check
```

### `fig_carola_ocv_paper.py` — figures A–E, JSON result, report, pin block

The main script. `compute()` recalculates every number of the five figures from
`data/carola/carola_ocv_channels.csv`; `pin_all()` then compares each number with
the three committed reference JSONs in `data/carola/reference/`
(**474 checks**: 414 exact, 60 within `1e-12` — the closed-form OLS report of the
site's analysis, whose sums ran in a different order). Recomputed and pinned
number by number: the state table, the per-girder and per-request (8-char
prefix) tables, the SDS block, the nested OLS / inverse-N battery of
`carola_echo_mask_gamma2.json`; the per-window `gamma2` and `n_masked` of
`carola_coherence_states.json` (**925** of its **1,454** stored windows are still
in the current export, and all 925 match exactly); and the committed map extract
`carola_bridge_osm.json`, pinned structurally (Carolabrücke, B 170, Q1044279,
two 12-node bridge ways, two razed DVB ways). A single deviation ends the run
with exit code 1.

`main()` writes the five figures to `figures/carola/`, the result JSON to
`data/carola/fig_carola_ocv_paper.json` and the English report to
`figures/carola/fig_carola_ocv_paper.md`.

```bash
python3 code/carola/fig_carola_ocv_paper.py            # full run, ~21 min
python3 code/carola/fig_carola_ocv_paper.py --quick    # ~6 min, without the permutation null
python3 code/carola/fig_carola_ocv_paper.py --figdir /tmp/x --json /tmp/y.json --md /tmp/z.md
```

`CAROLA_RES_CACHE=<path>` pickles the `compute()` result and reuses it on the
next run (the second run then takes seconds) — a developer aid, never used by the
wrapper.

Wall-clock times are only printed to stdout, never written into an artifact: with
the fixed seeds (RNG seed 7, 10000 bootstrap draws, 2000 for the strata, 1000
label permutations) two runs produce identical JSON, Markdown and PNG files (the
JSON's `quick` flag is the only difference between a `--quick` and a full run).

One string of the result JSON stays in the reference's own wording: the
`scenario_reading` of the `sds` block is a *pinned* value of
`carola_echo_mask_gamma2.json` and is therefore copied verbatim; the report
renders that reading in English instead.

### `fig_carola_ocv_months.py` — the echo-mask dimensions over the 24 months before the collapse

The **standalone companion** of the paper figure (the Carola analogue of
`../bautzen/fig_bautzen_ocv_amp_phase.py`). Where the A–E figures pool the two
states (`pre-collapse (healthy)` vs. `post-collapse`), this script keeps the
*time axis* and splits the committed chip export into the **24 rolling one-month
windows anchored at the event** (M-24 … M-1) plus the month the collapse opens
(M+1), so the reader sees **when** the mask moved, not only **that** it moved:

```
M-24  [2022-09-11, 2022-10-11)     pre-collapse (healthy)     16 chips   first record month
…
M-1   [2024-08-11, 2024-09-11)     pre-collapse (healthy)     56 chips   last
M+1   [2024-09-11, 2024-10-11)     post-collapse              60 chips
```

The collapse date opens M+1, so M-1 is pure pre and M+1 pure post and the two
are the only adjacent-month contrast the script tests (per-chip Cliff's delta
with a seeded percentile bootstrap CI, plus the same comparison against the
24 pre windows pooled). The fingerprint baseline is the **record-start** month
`M-24` (2022-09), and `compute()` reads only
`data/carola/carola_ocv_channels.csv` and recomputes each dimension's per-window
block with `carola_ocv_stats`; the figure is a 2 x 3 panel — (a)-(d) one dot per
chip per month for `gamma2`, `P`, `D`, `A` with the month median and the IQR,
(e) the `(P, D)` plane coloured by month on a symlog y axis, (f) the 6D
fingerprint, each month's median over the `M-24` baseline median on a log axis.
`main()` writes the figure to `figures/carola/fig_carola_ocv_months.png`, the
result JSON to `data/carola/fig_carola_ocv_months.json` and the English report
to `figures/carola/fig_carola_ocv_months.md`:

```bash
python3 code/carola/fig_carola_ocv_months.py            # ~35 s
bash figures/carola/fig_carola_ocv_months.sh            # the same, via the wrapper
python3 code/carola/fig_carola_ocv_months.py --figdir /tmp/x --json /tmp/y.json --md /tmp/z.md
```

It adds **no** number to the A–E pin contract and touches none of it (it never
reads the committed reference JSONs), and it is deterministic like the rest of
the package (fixed seed for the jitter and the bootstrap, no timestamps), so
repeated runs produce byte-identical PNG, JSON and Markdown. Its reading is
deliberately narrower than the pooled one: the 24 pre months carry 892 chips
(median 40 per month) and are **not** stationary (`gamma2`'s monthly median
ranges 0.0133 … 0.2975 over M-24..M-1, and the record start M-24 sits at 0.0687),
and against that month-to-month noise only `D` clears the bootstrap CI in the
adjacent-month contrast M+1 vs. M-1 (median 0.00729 → 0.00579, delta −0.355); `P`
— the sharpest pooled difference of the A–E figures — moves by only −0.192 here,
with a CI that spans zero, and γ2 runs *the other way* (0.0365 → 0.0624, delta
+0.078, CI spanning zero), so the month axis narrows the pooled verdict rather
than reproducing it.

The module and site names keep their `carola_` prefix although the folder already
says `carola/`: the import statements (`import carola_ocv_core as core`), the
reference file names and every command line in the reports stay parallel to the
other seven packages, which is what makes the eight packages readable as one diff.

### `carola_ocv_precollapse.py` — trend and late level on the pre-collapse record

Look at the orbit figure and the deck seems to move before the collapse, and its
ASC channel visibly alternates between two acquisition geometries. A look is not
a test, so this module asks the two questions that look raises, and only those:

1. **trend** — inside the pre-collapse observations, does `D`, `F` or `A` change
   systematically with time? Theil-Sen slope (per year, and as percent of the
   baseline median per year), a **moving-block bootstrap** CI (block 5 overpasses,
   2000 resamples) and Mann-Kendall with tie correction;
2. **late level** — is the last **three months** before the collapse different
   from the older pre-collapse baseline? Mann-Whitney U, Cliff's delta (positive
   = late-pre larger, so a sinking dimension is negative), moving-block CIs
   (block 3, 5000 resamples) for delta **and** the median difference, and **Holm**
   over the three dimensions per channel.

It invents no statistic: `theil_sen`, `mann_kendall`, `block_bootstrap_slope`,
`lag1_autocorr` and `runs_test` come from `../espoo/espoo_ocv_stats.py`,
`cliffs_delta` and `welch_mw_test` from `../lumo/lumo_ocv_stats.py`. The two
pieces the house did not have are `holm()` and the moving-block delta/median CI
(the house `bootstrap_delta_ci` resamples independently and would understate the
serial dependence).

```bash
python3 code/carola/carola_ocv_precollapse.py            # ~60 s, writes JSON + MD
python3 code/carola/carola_ocv_precollapse.py --selftest  # offline checks
python3 code/carola/carola_ocv_precollapse.py --verify    # recompute, compare
```

Outputs: `data/carola/support/carola_ocv_precollapse.json` and `.md` (event
`2024-09-11`, cut `2024-06-11`). `--verify` recomputes 126 fields from the
committed channel table and refuses a JSON whose `sha256` belongs to another
table; `--selftest` checks Holm by hand, the bootstrap's determinism and block
behaviour, that a planted trend/drop is found while noise is not, the calendar
arithmetic, and — as a guard — that **no state word reaches the report**.

What it finds, and what it does not:

* **No trend clears the interval.** ASC (n = 124, 2022-09-03 … 2024-09-09,
  revisit 5 d) gives Theil-Sen slopes of −0.0028 m/a, −8.8 Hz/a and
  −27 counts/a with CIs spanning zero (Holm p 0.74); DESC (n = 60,
  2022-09-12 … 2024-08-20) is positive and also spans zero (Holm p 0.68). The
  ASC series is strongly serially dependent (lag-1 −0.54, runs p < 0.001), which
  is exactly why the block interval, not the Mann-Kendall p, is the claim.
* **No late level shift is detectable in `D`, `F`, `A`.** All six cells are null
  (delta −0.04 … −0.09, CIs ±0.2 … 0.4, Holm p 1.00). A *large* sinking in these
  three dimensions before the collapse is therefore excluded; a modest one
  (|delta| < 0.2) is not, at 16 late ASC and 6 late DESC overpasses.
* **The window does not rescue it.** The 2/3/4/6-month ladder with the event
  moved ±12 d (24 cells) has no significant cell; the deltas lean negative
  (worst −0.22 at 4 months ASC), so the direction is consistent but unresolved.
* DESC has **no overpass in the last 22 d** before the collapse, so its late
  window ends before the event does — the report prints that gap per channel
  rather than letting a short window read as "unchanged until the end".


### `carola_ocv_decompose.py` — what is left of the post-collapse signature

*This is now a thin wrapper: the method lives in `../orbits/ocv_decompose.py`, whose registry describes every site (table, event, geometry rule, targets, model tree). The Carola numbers are its regression pin — 6096 leaves of the site JSON reproduce field by field.*

The orbit figure shows a large pre/post difference; this module asks how much of
it survives once the four candidate contributions are in the model: the observed
area `A` (explanatory, never a target), the acquisition geometry, temperature and
season, and a residual post effect.

Two **parallel ladders** are fitted per stratum and target:

```
A-adjusted   M0 (state) -> M1 (+log A) -> M2 (+T) -> M3 (+T, season)
A-free       M0 (state) -> M2- (+T)   -> M3- (+T, season)
```

`g(Y) = log Y` for `gamma2, P, D, F` and the identity for `S` (pixels). `f(log A)`
is linear with a 3-df restricted cubic spline as sensitivity, `h(T)` is linear and
`season` is the two-term annual harmonic in day-of-year (a twelve-month dummy
variant is the sensitivity). Every post coefficient carries a **moving-block
bootstrap** CI (block 5, 3 and 7 as stability check) in which blocks are drawn
inside each state, so none crosses the pre/post boundary or a stratum.

Strata are fitted separately — `ASC 36.9` (62/85), `ASC 44.8` (62/81), `DESC`
(60/86) — because the two ASC tracks are a strict 2-cycle and their post moves
point in opposite directions for `A`, `D` and `F`. Pooled ASC medians read
`A +10 %`, `D +8 %`, `F +10 %`: the strata cancel, so there is no pooled ASC
headline. The joint ASC model (`post`, `angle`, `post x angle`) is family B and
confirms it: the interaction is significant for `gamma2` (F 8.75, p 0.003) and
`P` (F 6.29, p 0.013).

```bash
python3 code/carola/carola_ocv_decompose.py            # ~90 s, writes JSON + MD
python3 code/carola/carola_ocv_decompose.py --selftest  # four planted worlds
python3 code/carola/carola_ocv_decompose.py --verify    # 228 fields recomputed
python3 code/orbits/ocv_decompose.py --site all        # Carola, CTS, Haltern, Morandi
                                                       # + data/decompose_summary.md
```

Outputs `data/carola/support/carola_ocv_decompose.{json,md}`. `--selftest` builds
four synthetic worlds where the answer is known before the fit — the effect made
by the area alone (A-adjusted CI covers no effect, A-free does not), by the
season (M0 sees it, M3 removes it), a genuine step (survives at −0.40) and nothing
(all CIs cover no effect) — plus the bootstrap's determinism, the no-crossing
rule, the refusal of `T24`/`T72`, and a wording guard that the report never
contains `area-driven`, `caused by` or `explained by`.

What it finds: after area, geometry, temperature and season, and after Holm over
the five targets in each stratum, **one cell survives** — `gamma2` in ASC 36.9,
ratio 2.10, CI [1.50, 2.85], Holm < 0.001. Two more have intervals that exclude
1 without surviving Holm (`ASC 44.8 F` 0.965 [0.927, 0.992], Holm 0.110;
`DESC P` 0.797 [0.678, 0.956], Holm 0.100). Most of the rest **co-moves with the
observed area**: `DESC D` falls to 0.32 raw and to 0.99 after `log A`, `DESC F`
from 0.31 to 1.05. `ASC 44.8 F` even changes sign between the ladders, which the
verdict column reports as a sign flip rather than a result.

Run across the four sites (`--site all`, summary in
`data/decompose_summary.md`), the answer is that the Carola pattern **does**
generalise — the post-collapse signature is largely co-movement with the
observed area — with exactly one exception after the pre-registered robustness
gate (`data/analysis_freeze.md`): **Morandi `D`** (ratio **0.31**, CI
[0.22, 0.52], Holm < 0.001, placebo rank 1 of 8, n_post = 9). At Carola itself
every dimension collapses onto the area, including `gamma2` (A-free +0.540 →
A-adjusted +0.274, CI [−0.20, 0.61]); at CTS `P` is direction-consistent but never
Holm-firm (0.93, Holm 0.21); Haltern's `gamma2` shrinks to 1.40 with a ladder
that disagrees, and Morandi's `gamma2`, `P` and `S` show no independent effect
(`gamma2` placebo rank 9 of 8). Morandi and Haltern are `thin` (9 and 11 post
overpasses) and carry that label wherever they appear.

Two warnings the report states, because they are easy to get backwards:

* **`A` is not a neutral confounder.** The collapse moves the area (DESC −77 %),
  so conditioning on it removes part of the change being measured — and for a
  mediator it can also *create* association where the raw data show none
  (`ASC 36.9 P`: raw ratio 1.03, A-adjusted 1.05 with an interval that grazes 1).
  The A-adjusted column measures co-movement; the A-free column is the structural
  signature. Both are reported for that reason.
* **`T24` and `T72` do not exist in this record** — one temperature per overpass,
  5 to 12 days apart. `T0` is primary and `T_smooth30` (a 30-day centred mean of
  the overpass temperatures) the sensitivity; neither is a deck temperature.

The multivariate check (standardised `M3-` residuals, circular-shift null) finds
no global shift beyond the covariates: p = 0.21 (ASC 36.9), 0.07 (ASC 44.8),
0.20 (DESC). The placebo panel ranks the real signature against eight pseudo
events inside the pre period: rank 1 of 8 for DESC (raw and size-matched) and for
ASC 44.8 raw (5 of 8 when size-matched), but 9 of 8 for ASC 36.9 — its large
`gamma2` move is *not* exceptional against the pre period's own variability.
DESC has 56 % of its post area below the pre 5th percentile, so its A-adjusted number
is an extrapolation and the common-support subset is the one to read.
