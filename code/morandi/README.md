# `code/morandi/` — Morandi site package: analysis and plotting scripts

Four Python modules, no package, no dependencies beyond the standard library plus
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

## The shipped plate

This package's figure in `figures/` is
[`figures/morandi/fig_morandi_ocv_trajectory.md`](../../figures/morandi/fig_morandi_ocv_trajectory.md)
— the standard OCV trajectory plate on the **support-sized** chip, drawn by
`code/cases/ocv_trajectory.py --site morandi` from this package's registry entry
(reference `pre-collapse (healthy)`, 198 chips; tail `post-collapse`, 9 chips). It
reads `data/morandi/support/morandi_ocv_channels.csv` and pins every number it
draws; the wide interval of its contrast is a property of the 9-chip tail and is
printed as such.

The A–E paper plate and the 24-month companion that this package once fed were
retired with the twelve-plate scope. Their frozen JSONs
(`data/morandi/fig_morandi_ocv_paper.json`, `..._months.json`) stay committed as
provenance, but no script recomputes them any more.


