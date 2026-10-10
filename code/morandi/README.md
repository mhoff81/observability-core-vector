# `code/morandi/` — Morandi site package: analysis and plotting scripts

Six Python modules, no package, no dependencies beyond the standard library plus
numpy / matplotlib / scipy (scipy is only used for the p-values; without it the
scripts still run, just with fewer p-values). This is the Ponte Morandi
(Polcevera, Genova) half of the repository; the LUMO tower lives in `../lumo/`,
Bautzen in `../bautzen/`, KDLO in `../kdlo/` and Carola in `../carola/`.

The site is the collapsed **Ponte Morandi / Viadotto Polcevera** (Genova,
44.4243 N 8.8906 E), a three-span balanced-cantilever viaduct. Its record is a
**per-track chip export**: per acquisition the pipeline stored one 80 x 80
complex chip per Sentinel-1 observation geometry of the deck — `A_asc`
(ASC rel-15 IW1 VV) and `A_des` (DESC rel-66 IW1 VV), see
`data/morandi/morandi_tracks.txt`. The package analyses the two states of the
*same* structure around the collapse date `2018-08-14`:

```
pre-collapse (healthy)   <   2018-08-14   <=   post-collapse
```

so season, weather and pipeline generation all move with the state. The
observability core vector is `OCV = [gamma2, P, D]`; the full mask vector is
`x = [gamma2, P, D, A, F, S]`, carried in **two layers**:

* **echo mask** — the whole 80 x 80 chip (the layer of the 6D vector), the
  project-wide bridge-deck echo rule `intensity >= 0.30 x peak` **and**
  `intensity >= 5 x np.median(intensity)`, `A >= 2`;
* **deck mask** — the very same rule on the *track's own* chip, i.e. the same
  six dimensions restricted to `A_asc` / `A_des` (the layer this package adds).

The mask layer rests on a **committed cache**
(`data/morandi/morandi_windows_mask_cache.txt`): the per-date 400 x 400 payloads
do not belong in the repository, so the cache carries the echo mask of every
chip plus the sufficient statistics, and the build re-verifies the rule on it.

Neither layer's `A` is the site's own `coherence_masked_pixels` column: that
column is a *fixed 640-px (10 %) quantile mask* (`morandi_deck_channels.csv` /
`morandi_registration.json`) while this package recomputes the echo mask from
the payloads — the two agree on **0 of 232** rows, so the CSV keeps the export
as provenance only (exported as `A_export` in the figures, panel E).

Two structural caveats run through every panel: the analysed 80 x 80 window is
**not on the collapsed span** and the archive holds no persistent scatterer at
this footprint (see `morandi_analysis/MORANDI_SCATTERER_LIMITS.md`); and only
**11** chips (5 `A_asc` + 6 `A_des`) survive the event against **221** before
it, so every post-state interval is wide. The figure is an **honest null**.

All paths below are relative to the **repository root** and the scripts are
always started from there (the wrappers `../../figures/morandi/*.sh` do that):
inputs in `data/morandi/`, figures in `figures/morandi/`.

### `morandi_ocv_stats.py` — statistics layer (verbatim port)

`_stats`, `welch_mw_test`, `cliffs_delta` (+ `cliffs_delta_bruteforce`),
`bootstrap_delta_ci`, `delta_block`, `sign_test`, `sds_from`, `by_state`, `vals`,
`season_of`, `strata_pair`, `fit_lda`, `predict_lda`, `loo_cv`,
`permutation_test`, `fit_linear_score`, `score_row`, `oof_scores`, `delta_sds`,
`permutation_p`, `paired_bootstrap_diff`, plus the port of the site's own
correlation/regression battery (`spearman`, `mannwhitney`, `pair_block`,
`ols_fit`, `f_test`). Copied line by line from `../carola/carola_ocv_stats.py`,
which is itself a verbatim port of the LUMO project analysis scripts, so a
difference in the Morandi numbers cannot come from a difference in the
statistics. The only Morandi edit is naming: the two states are the site's
collapse alphabet `pre-collapse (healthy)` / `post-collapse`
(`STATE_ORDER`, `STATE_SHORT`).

### `morandi_ocv_core.py` — data layer, constants and loaders

Paths, the analysis channels (`OCV = [gamma2, P, D]`), the 6D feature vector
(`FEATURES = [gamma2, P, D, A, F, S]`, `DIMS`), the axis/title labels, the
weather and control columns (`WEATHER`, `CONTROLS`), the two tracks
(`GIRDERS = [A_asc, A_des]`), the model/seed constants (`MODELS`,
`HEADLINE_PAIR = x_6d`, `LAMBDA_GRID`, `HEADLINE_LAMBDA = 0.5`, `RNG_SEED = 7`,
`N_BOOT = 10000`, `N_PERM = 1000`) and the loaders (`load_csv`, `load_csv_meta`,
`load_cache`, `load_manifest`, `load_measurements`, `load_segments`,
`load_reference`). `MASK_COLUMNS` lists the echo-mask columns the generator
appends, so that `CSV_COLUMNS` (and the CSV header) carry the whole 6D vector;
the loader also injects `coherence_masked_pixels` from the committed
registration (`REF_STATES`) so the export column is provenance in the CSV too.

It also carries the small shared aggregates (`row_brief`, `summ`, `by_girder` /
`by_track`, `summary_by_state`, `day_of`, `state_of`, `echo_mode_of`,
`track_index_of`). The track is called a "girder" only through the alias
`by_girder` and the constant `GIRDERS`, so the figure scripts read identically
to Carola's.

### `morandi_ocv_masks.py` — the two mask layers (port of `../carola/carola_ocv_masks.py`)

The module owns the mask rule of both layers: `echo_mask` (`0.30 x peak`,
`5 x np.median`, `A >= 2`), `largest_component` / `connected_components_8`
(port of the LUMO `_largest`), `gamma2_of` (`|sum z|^2 / (A sum |z|^2)`,
clamped to [0, 1]), `vector_from_row` (`A`, `D = A / bbox_area`, `F`, `S`, the
six dimensions, the bbox/centroid/peak geometry and the masked `gamma2`),
`add_persistence` (`P` = mean pixel frequency of the own mask over all chips of
the majority shape) and the cache I/O (`load_cache`, `load_manifest`, `extract`,
`check`).

`extract()` is the one-time offline step: it reads the per-date 400 x 400 rect
payloads of the site project (`/home/projects/morandi_analysis`,
`ref_rects_A_asc_400.txt` / `ref_rects_A_des_400.txt`), keeps the central
`[160:240, 160:240]` 80 x 80 crop, applies the mask rule to every chip, drops
byte-identical duplicate payloads, and writes the cache plus its manifest (the
source files' sha256 included, and the counts of the 234 `400x400` payloads and
the 5 clipped/odd-shaped ones excluded). After that the cache is the only mask
input. `python3 code/morandi/morandi_ocv_masks.py --check` is the acceptance
test of the committed cache: every masked pixel must satisfy both thresholds,
`A >= 2`, and the recomputed `gamma2` must agree with the stored sufficient
statistics (232 `80 x 80` rows; 5 clipped/odd payloads excluded and 2 of the 234
`400x400` payloads carried no echo mask, so both are dropped, leaving 232).

### `morandi_ocv_channels_csv.py` — build `data/morandi/morandi_ocv_channels.csv`

Reads the committed cache and the committed measurement extract
(`morandi_weather_table.csv`), recomputes the six mask dimensions of every chip,
derives `segment`, `track_index`, `state`, `pre_collapse`, `day`, `month`,
`year`, `season` and `echo_mode`, injects `coherence_masked_pixels` from the
committed registration as provenance, and writes the channel table (**232** rows,
**59** columns) plus `morandi_ocv_channels_meta.json` (file sha256s, the mask
rule, the guard counts, the derived-column definitions, the
`coherence_masked_pixels` cross-check, and the per-state summary of every
channel). The build **refuses to write** if a chip has no echo mask or if the
mask layer cannot be verified.

```bash
python3 code/morandi/morandi_ocv_channels_csv.py                 # build + verify
python3 code/morandi/morandi_ocv_channels_csv.py --verify-only   # verify only
# one-time re-extraction of the mask cache (needs the site's 400 x 400 rects):
python3 code/morandi/morandi_ocv_masks.py --extract
python3 code/morandi/morandi_ocv_masks.py --check
```


### `fig_morandi_ocv_paper.py` — figures A–E, JSON result, report, pin block

The main script. `compute()` recalculates every number of the five figures from
`data/morandi/morandi_ocv_channels.csv`; `pin_block()` then compares each number
with the three committed references in `data/morandi/reference/` (**47 checks**).
Recomputed and pinned number by number: the generated table's meta block (site,
event, generator, 232 rows, 221 pre / 11 post, 59 columns, `features` / `ocv` /
`girders`, the `coherence_masked_pixels` cross-check), the committed mask-cache
manifest (`n_rows`, the `[160,240,160,240]` crop, the window rule, the two
source sha256s, the 234 `400x400` and 2 no-echo counts), the per-track counts, the
site's fixed 640-px deck-mask table (`morandi_deck_channels.json`) and the
viaduct geometry (`morandi_geometry.json`, pinned structurally). A single
deviation ends the run with exit code 2.

```
A  raw distributions of the six channels per state (all 232 chips)
B  in-sample effect sizes: Cliff's delta + bootstrap CI and SDS, pooled and per track
C  2-class LDA (pre vs. post) leave-one-out accuracy + permutation null
D  the OCV plane (gamma2, P, D): does the core vector separate the states?
E  the evidence against reading the contrast as damage: weather controls,
   orbit/track composition, echo-mode mix, and the export-vs-recomputed mask
```

`main()` writes the five figures to `figures/morandi/`, the result JSON to
`data/morandi/fig_morandi_ocv_paper.json` and the English report to
`figures/morandi/fig_morandi_ocv_paper.md`.

```bash
python3 code/morandi/fig_morandi_ocv_paper.py            # full run, ~3.5 min (permutation null)
python3 code/morandi/fig_morandi_ocv_paper.py --quick    # ~23 s, without the permutation null
python3 code/morandi/fig_morandi_ocv_paper.py --figdir /tmp/x --json /tmp/y.json --md /tmp/z.md
```

Its headline is an **honest null**. Of the six echo-mask channels, **none**
clears the bootstrap CI in the pooled contrast (`gamma2` +0.218, CI
[-0.164, 0.582]; `P` -0.130, CI [-0.392, 0.221]); the only channel that does is
`temperature_c` (+0.635, CI [0.454, 0.806]) — a co-variate, which is exactly the
control the panels are there to show. The 6D LDA reaches 0.746 raw accuracy but
only **0.521 balanced accuracy** (chance 0.5) on the 221 / 11 imbalance, and its
label-permutation p is 0.077 (1000 draws), so it does not clear the null; the
weather-only competitor reaches 0.685 / 0.662 (balanced-accuracy p 0.040), and
the out-of-fold OCV scores carry SDS 0.06 against 1.64 for the 6D vector. The
states are not separable in the analysed window. The per-track deck-mask table
*does* separate within one track — `A_asc` moves (`D` -0.58, `S` +0.71, CI≠0)
while `A_des` does not — but with only 5-6 post chips per track the report
carries that split as one more co-variate to read, not a verdict.

### `fig_morandi_ocv_months.py` — the echo-mask dimensions month by month around the collapse

The **standalone companion** of the paper figure. Where the A–E figures pool the
two states (`pre-collapse (healthy)` vs. `post-collapse`), this script keeps the
*time axis* and splits the committed chip export into the four rolling one-month
windows anchored at the event, so the reader sees **when** the mask moved, not
only **that** it moved:

```
M-3   [2018-05-14, 2018-06-14)     pre-collapse (healthy)     10 chips
M-2   [2018-06-14, 2018-07-14)     pre-collapse (healthy)      8 chips
M-1   [2018-07-14, 2018-08-14)     pre-collapse (healthy)     10 chips
M+1   [2018-08-14, 2018-09-14)     post-collapse              11 chips
```

The collapse date opens M+1, so M-1 is pure pre and M+1 pure post and the two
are the only adjacent-month contrast the script tests (per-chip Cliff's delta
with a seeded percentile bootstrap CI, plus the same comparison against the
three pre windows pooled). `compute()` reads only
`data/morandi/morandi_ocv_channels.csv` and recomputes each dimension's
per-window block with `morandi_ocv_stats`; the figure is a 2 x 3 panel — (a)-(d)
jittered strips per window for `gamma2`, `P`, `S`, `D` with the window median,
the IQR and the collapse cut, (e) the `(P, D)` plane coloured by window on a
symlog y axis, (f) the 6D fingerprint, each window's median over the
pre-baseline median on a log axis. `main()` writes the figure to
`figures/morandi/fig_morandi_ocv_months.png`, the result JSON to
`data/morandi/fig_morandi_ocv_months.json` and the English report to
`figures/morandi/fig_morandi_ocv_months.md`:

```bash
python3 code/morandi/fig_morandi_ocv_months.py            # ~7 s
bash figures/morandi/fig_morandi_ocv_months.sh            # the same, via the wrapper
python3 code/morandi/fig_morandi_ocv_months.py --figdir /tmp/x --json /tmp/y.json --md /tmp/z.md
```

It adds **no** number to the A–E pin contract and touches none of it (it never
reads the committed reference JSONs), and it is deterministic like the rest of
the package (fixed seed for the jitter and the bootstrap, no timestamps), so
repeated runs produce byte-identical PNG, JSON and Markdown. Its reading is
narrower still than the pooled one: with only 8-11 chips per window **no**
dimension clears the bootstrap CI in either contrast (M+1 vs. M-1 or M+1 vs. the
pooled pre side), so the month axis shows the co-variate noise the pooled figure
already flags rather than a resolvable step.

Wall-clock times are only printed to stdout, never written into an artifact: with
the fixed seeds (RNG seed 7, 10000 bootstrap draws, 1000 label permutations) two
runs produce identical JSON, Markdown and PNG files (the JSON's `quick` flag is
the only difference between a `--quick` and a full run).

The module and site names keep their `morandi_` prefix although the folder
already says `morandi/`: the import statements (`import morandi_ocv_core as
core`), the reference file names and every command line in the reports stay
parallel to the other sites, which is what makes the packages readable as one
diff.

