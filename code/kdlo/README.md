# `code/kdlo/` — KDLO site package: analysis and plotting scripts

Five Python modules, no package, no dependencies beyond the standard library plus
numpy / matplotlib / scipy (scipy is only used for the p-values of the control
tests; without it the scripts still run, just with fewer p-values). This is the
KDLO half of the repository; the LUMO tower site lives in `../lumo/` and Bautzen
in `../bautzen/`.

The mask layer follows the LUMO rule (`0.30 x peak`, `5 x median`) and carries the
same 6D echo-mask vector `x_KDLO = [gamma2, P, D, A, F, S]`, but the KDLO backend
stored only the mask *size* (`A = coherence_masked_pixels`). The masks are
therefore recomputed from a **committed** cache of the 64 full-dwell strips the
pipeline decoded (`data/kdlo/strips/`, 11 x 400 complex samples each, ~70 KB per
file), produced by `tower/backend/src/cli/analysis/kdlo_strip_capture.rs` and read
here by `kdlo_ocv_masks.py`. The OCV proper is `[gamma2, P, D]`; `A` and the
masked coherence are pinned to the pipeline's own `coherence_masked_pixels` and
`coherence_gamma2` on every acquisition.

All paths below are relative to the **repository root** and the scripts are always
started from there (the wrapper `../../figures/kdlo/fig_kdlo_ocv_paper.sh` does
that): inputs and outputs in `data/kdlo/`, figures in `figures/kdlo/`.

### `kdlo_ocv_stats.py` — statistics layer (verbatim port)

`_stats`, `welch_mw_test`, `cliffs_delta` (+ `cliffs_delta_bruteforce`),
`bootstrap_delta_ci`, `delta_block`, `sign_test`, `sds_from`, `by_state`, `vals`,
`season_of`, `strata_pair`, `fit_lda`, `predict_lda`, `loo_cv`,
`permutation_test`, `fit_linear_score`, `score_row`, `oof_scores`, `delta_sds`,
`permutation_p`, `paired_bootstrap_diff`. Copied line by line from
`../lumo/lumo_ocv_stats.py`, which is itself a verbatim port of the LUMO project
analysis scripts, so a difference in the numbers cannot come from a difference in
the statistics. The only KDLO edits are naming: the state field is `state`
(LUMO: `damage_label`) and the two states are **epochs** (`pre`, `rebuild`)
instead of damage labels.

### `kdlo_ocv_core.py` — data layer, constants and the committed aggregates

Paths, the analysis channels (`OCV = [gamma2, P, D]`), the 6D feature vector
(`FEATURES = [gamma2, P, D, A, F, S]`, `DIMS`, `DIM_LABEL`) and the six controls
(`CONTROLS`, `WEATHER`), the epoch windows (`STATE_WINDOWS`, `SINGLETON_DAY`), the
pipeline thresholds (`COMPACT_MASK_PX = 10`, `MASK_CLUTTER_PX = 25`), the
model/seed constants (`MODELS`, `HEADLINE_PAIR`, `LAMBDA_GRID`,
`HEADLINE_LAMBDA = 0.5`, `RNG_SEED = 7`, `N_BOOT = 10000`) and the loaders
(`load_csv`, `load_csv_meta`, `load_raw_parts`, `load_db_extra`, `load_reference`).
`MASK_COLUMNS` lists the echo-mask columns the generator appends, so that
`CSV_COLUMNS` (and the CSV header) carry the whole 6D vector.

It also carries the two aggregate blocks the committed channel JSONs were written
with — `summ`, `row_brief`, `echo_modes_of`, `phase_observable_of`,
`part_stats` — and the two pin-only mirrors of the files that cannot be
recomputed from the record: `availability_checks` (128 structural checks of
`kdlo_s1_availability.json`) and `fem_checks` (67 checks of
`kdlo_fem_frequencies.json`).

### `kdlo_ocv_masks.py` — the echo-mask layer (port of `../lumo/lumo_ocv_masks.py`)

`read_complex_bin`, `connected_components_8`, `vector_from_strip` (returns the mask
plus `A`, `D = A / bbox_area`, `F`, `S`, the bbox/centroid/peak geometry and the
masked `gamma2`), `add_persistence` (`P` = mean pixel frequency of the own mask
over all dates of the majority strip shape), `load_vectors` and `check`.

The one deliberate difference from the LUMO port is the background median of the
mask rule: it is the **upper** median `sorted(v)[len // 2]`, the value the Rust
pipeline uses (`sarsolve::coherence::coherence_mast_echo`), because `A` has to
reproduce `coherence_masked_pixels` bit for bit. `np.median` is kept as the
diagnostic `A_np_median` / `gamma2_np_median` pair and never reported.

```bash
python3 code/kdlo/kdlo_ocv_masks.py           # acceptance test on all strips
python3 code/kdlo/kdlo_ocv_masks.py --json    # same, plus the derived vectors
```

The check recomputes the mask of every committed strip and requires
`A == coherence_masked_pixels`, `abs(gamma2_mask - coherence_gamma2) <= 1e-12` and
a 400 x 11 strip shape. It is a self-contained acceptance test of the mask layer
that needs neither database nor network — the per-strip sha256 values of
`strips/manifest.json` are the pins behind it.

### `kdlo_ocv_channels_csv.py` — generator of the acquisition table

Merges the two committed raw API dumps (`data/kdlo/fig_kdlo_1_raw.json`,
`fig_kdlo_2_raw.json`) with the **one-time, read-only** Postgres extraction of
nine columns that the API rows do not contain (`kdlo_db_extra.json`: phase SNR,
phase observability, snow, precipitation, gusts, intensity, incidence angle),
derives `part`, `state`, `month`, `year`, `season` and `echo_mode`, attaches the
echo-mask columns of the 6D vector (`A`, `D`, `F`, `S`, `P`, the mask geometry and
the recomputed `gamma2_mask`) from the committed strip cache via
`kdlo_ocv_masks.attach_masks`, and writes `data/kdlo/kdlo_ocv_channels.csv`
(60 rows, one per acquisition, 58 columns) plus `kdlo_ocv_channels_meta.json`
(provenance: file sha256s, the windows, the definitions, the counts, and the
`strip_cache` block with the mask rule and the manifest sha256).

The build **fails** if a row has no echo mask or if any `A`/`gamma2_mask`
deviates from the record, so the CSV cannot be written with a lookalike mask. The
verification repeats that comparison and additionally re-hashes every committed
strip (`--no-strip-verify` skips only the hashing).

```bash
python3 code/kdlo/kdlo_ocv_channels_csv.py                 # build + verify
python3 code/kdlo/kdlo_ocv_channels_csv.py --verify-only   # verify only
python3 code/kdlo/kdlo_ocv_channels_csv.py --verify-only --no-strip-verify
# one-time re-extraction of the nine DB columns (needs the `tower-postgres`
# container; the result is vendored, so this is not part of a normal run):
python3 code/kdlo/kdlo_ocv_channels_csv.py --fetch-db-extra
```

### `kdlo_capture_requests.py` — the request list of the strip cache

Reads the 64 measurement rows of `onboarder.insar_measurements` for the KDLO asset
(read-only, `docker exec psql`) and writes `data/kdlo/strips/capture_requests.json`:
site point, strip spec and, per acquisition, the CDSE burst Id, orbit direction,
relative orbit, sub-swath, polarisation, product name and content start. The Rust
binary turns that list into the strips; this script is what makes the capture
reproducible from the record.

### `fig_kdlo_ocv_paper.py` — figures A–E, JSON result, report, pin block

The main script. `compute()` recalculates every number of the five figures from
`data/kdlo/kdlo_ocv_channels.csv`; `pin_all()` then compares each number with the
four committed reference JSONs in `data/kdlo/reference/`, row by row (2985 checks
in a full run, 2975 in `--quick`; 128 + 67 of them are the structural mirrors of
the two pin-only files), verifies the echo-mask layer (the `strip_cache` block of
the metadata, every committed strip's sha256 through the manifest sha256, and per
row `A == coherence_masked_pixels`, `abs(gamma2_mask - coherence_gamma2) <= 1e-12`,
`D == A / bbox_area`, `P` in [0, 1]). A single deviation ends the run with exit
code 1. `main()` writes the five figures to `figures/kdlo/`, the result JSON to
`data/kdlo/fig_kdlo_ocv_paper.json` and the English report to
`figures/kdlo/fig_kdlo_ocv_paper.md`.

```bash
python3 code/kdlo/fig_kdlo_ocv_paper.py            # full run, ~1 min
python3 code/kdlo/fig_kdlo_ocv_paper.py --quick    # without the permutation null
python3 code/kdlo/fig_kdlo_ocv_paper.py --figdir /tmp/x --json /tmp/y.json --md /tmp/z.md
```

Wall-clock times are only printed to stdout, never written into an artifact: with
the fixed seeds (RNG seed 7, 10000 bootstrap draws, 2000 for the strata, 1000
label permutations) two runs produce identical JSON, Markdown and PNG files.

The module and site names keep their `kdlo_` prefix although the folder already
says `kdlo/`: the import statements (`import kdlo_ocv_core as core`), the
reference file names and every command line in the reports stay parallel to the
other two sites, which is what makes the three packages readable as one diff.
