# `code/` — analysis and plotting scripts

Five Python modules, no package, no dependencies beyond the standard library plus
numpy / matplotlib / scipy (scipy is only used for the p-values of the control
tests; without it the scripts still run, just with fewer p-values). This folder
is the repository root of every figure: all `data/` paths below are relative to
the repository root, and the scripts are always started from there (the wrapper
in `figures/` does that).

### `lumo_ocv_stats.py` — statistics layer (verbatim port)

`_stats`, `welch_mw_test`, `cliffs_delta`, `delta_sds`, `bootstrap_delta_ci`,
`monotonicity` (Jonckheere-Terpstra + Spearman), `strata`, `loo_cv`,
`permutation_test`, `oof_scores`, `permutation_p`, `paired_bootstrap_diff`,
`feature_dataset`, `mask_definition`, `by_state`, `vals`. Ported line by line
from the LUMO project's analysis scripts so that the numbers stay comparable with
the committed results — no reimplementation, no "improved" variant.

### `lumo_ocv_masks.py` — mask layer (verbatim port)

The echo-mask rule of `analyze_lumo_echo_mask_vector.py`: threshold at
`0.30 x peak` and `5 x median`, at least 2 pixels, then
`A = n_masked`, `D = A / bbox_area`, `F = 8-neighbourhood components`,
`S = ||mask centroid - peak||`, plus the persistence `P` (mean pixel frequency
across the 177 majority-shape masks).

### `lumo_ocv_core.py` — data layer and constants

Paths, the channel/feature order (`FEATURES = [gamma2, P, D, A, F, S]`), the
three feature sets of Figure D (`MODELS`), the shrinkage/seed/permutation
constants, the label maps (`SHORT`, `SHORT_ALC`, `LABEL_EN`) and the loaders:
`load_csv()` (the extended overpass table), `load_csv_meta()`,
`load_reference()` (the five committed JSONs), `feature_dataset()`,
`n_by_state()`, `sha256()`.

### `lumo_ocv_channels_csv.py` — generator of the extended overpass table

Reads `data/lumo_channels.csv` (the copy of the committed overpass table) and
the machine-local burst cache `../../lumo_dam6_analysis/monthly_bursts/`
(a pipeline artefact, ~100 MB, not committed), computes the mask geometry of the
178 coherence overpasses with `lumo_ocv_masks.py` and appends the columns
`A, D, F, S, P` (plus the auditable ingredients `mask_rows, mask_cols,
bbox_area, centroid_r/c, peak_r/c`) to `data/lumo_ocv_channels.csv`. It verifies
while writing: every overpass of the committed table must be found in the cache,
and `A` must equal the committed `n_masked`.

```bash
python3 code/lumo_ocv_channels_csv.py                  # build + verify
python3 code/lumo_ocv_channels_csv.py --verify-only    # verify only
```

### `fig_lumo_ocv_paper.py` — figures A–E, JSON result, report, pin block

The main script. `compute()` recalculates every number of the five figures from
`data/lumo_ocv_channels.csv`; `pin_all()` then compares each number with the
committed reference JSONs in `data/reference/`, key by key and row by row (2982
checks in a full run, 2892 in `--quick`; the German documentation strings of the
references are translated
through `core.LABEL_EN`). A single deviation ends the run with exit code 1.
`main()` writes the five figures to `figures/`, the result JSON to
`data/fig_lumo_ocv_paper.json` and the English report to
`figures/fig_lumo_ocv_paper.md`.

```bash
python3 code/fig_lumo_ocv_paper.py            # full run, ~5 min
python3 code/fig_lumo_ocv_paper.py --quick    # ~2 min: the two permutation
                                              # tests are copied from the refs
python3 code/fig_lumo_ocv_paper.py --figdir /tmp/x --json /tmp/y.json --md /tmp/z.md
```

Wall-clock times are only printed to stdout, never written into an artifact: with
the fixed seeds (`RNG_SEED = 11`, 2000 bootstrap draws, 1000 / 500 permutation
draws) two runs produce identical JSON, Markdown and PNG files.
