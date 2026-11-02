# `code/carola/` — Carola site package: analysis and plotting scripts

Four Python modules, no package, no dependencies beyond the standard library plus
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
always started from there (the wrapper `../../figures/carola/fig_carola_ocv_trajectory.sh`
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

## The shipped plate

This package's figure in `figures/` is
[`figures/carola/fig_carola_ocv_trajectory.md`](../../figures/carola/fig_carola_ocv_trajectory.md)
— the standard OCV trajectory plate on the **support-sized** chip, drawn by
`code/cases/ocv_trajectory.py --site carola` from this package's registry entry
(reference `pre-collapse (healthy)`, restricted to 2022-09 … 2023-12; tail
`post-collapse`; the 63 healthy chips outside the window drawn as the held-out
class). It reads `data/carola/support/carola_ocv_channels.csv` and pins every
number it draws.

The A–E paper plate and the 24-month companion that this package once fed were
retired with the twelve-plate scope. Their frozen JSONs
(`data/carola/fig_carola_ocv_paper.json`, `..._months.json`) stay committed as
provenance, but no script recomputes them any more; the channel table below is
what every surviving reader builds on.

