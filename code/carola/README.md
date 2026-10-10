# `code/carola/` — Carola site package: analysis and plotting scripts

Five Python modules, no package, no dependencies beyond the standard library plus
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

### `fig_carola_ocv_months.py` — the echo-mask dimensions month by month around the collapse

The **standalone companion** of the paper figure (the Carola analogue of
`../bautzen/fig_bautzen_ocv_amp_phase.py`). Where the A–E figures pool the two
states (`pre-collapse (healthy)` vs. `post-collapse`), this script keeps the
*time axis* and splits the committed chip export into the four rolling
one-month windows anchored at the event, so the reader sees **when** the mask
moved, not only **that** it moved:

```
M-3   [2024-06-11, 2024-07-11)     pre-collapse (healthy)     42 chips
M-2   [2024-07-11, 2024-08-11)     pre-collapse (healthy)     52 chips
M-1   [2024-08-11, 2024-09-11)     pre-collapse (healthy)     56 chips
M+1   [2024-09-11, 2024-10-11)     post-collapse              60 chips
```

The collapse date opens M+1, so M-1 is pure pre and M+1 pure post and the two
are the only adjacent-month contrast the script tests (per-chip Cliff's delta
with a seeded percentile bootstrap CI, plus the same comparison against the
three pre windows pooled). `compute()` reads only
`data/carola/carola_ocv_channels.csv` and recomputes each dimension's per-window
block with `carola_ocv_stats`; the figure is a 2 x 3 panel — (a)-(d) jittered
strips per window for `gamma2`, `P`, `S`, `D` with the window median, the IQR and
the collapse cut, (e) the `(P, D)` plane coloured by window on a symlog y axis,
(f) the 6D fingerprint, each window's median over the pre-baseline median on a
log axis. `main()` writes the figure to `figures/carola/fig_carola_ocv_months.png`,
the result JSON to `data/carola/fig_carola_ocv_months.json` and the English
report to `figures/carola/fig_carola_ocv_months.md`:

```bash
python3 code/carola/fig_carola_ocv_months.py            # ~14 s
bash figures/carola/fig_carola_ocv_months.sh            # the same, via the wrapper
python3 code/carola/fig_carola_ocv_months.py --figdir /tmp/x --json /tmp/y.json --md /tmp/z.md
```

It adds **no** number to the A–E pin contract and touches none of it (it never
reads the committed reference JSONs), and it is deterministic like the rest of
the package (fixed seed for the jitter and the bootstrap, no timestamps), so
repeated runs produce byte-identical PNG, JSON and Markdown. Its reading is
deliberately narrower than the pooled one: the three pre windows are **not**
stationary (`gamma2`'s median roughly halves twice over M-3..M-1, 0.158 → 0.072
→ 0.037), and against that month-to-month noise only `D` clears the bootstrap CI
in the adjacent-month contrast M+1 vs. M-1 (median 0.0073 → 0.0058, delta −0.36);
`P` — the sharpest pooled difference of the A–E figures — moves by only −0.19
here, with a CI that spans zero, so the month axis narrows the pooled verdict
rather than reproducing it.

The module and site names keep their `carola_` prefix although the folder already
says `carola/`: the import statements (`import carola_ocv_core as core`), the
reference file names and every command line in the reports stay parallel to the
other seven packages, which is what makes the eight packages readable as one diff.
