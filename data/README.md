# `data/` — one subfolder per site: input tables and reference results

`lumo/` and `bautzen/`. Everything in here is either a **copy** of a committed
artifact of the site's analysis project or a file **generated** by a script in
`../code/<site>/`; nothing is edited by hand. Every script reads and writes only
inside its own `data/<site>/`, so the two sites can never overwrite each other.

## `data/lumo/` — LUMO lattice tower (Hannover)

### Inputs (committed copies)

| file | what it is | read by |
| --- | --- | --- |
| `lumo_channels.csv` | the LUMO overpass table as committed in the LUMO project: one row per Sentinel-1 overpass of the tower (date, orbit, acquisition mode, `coherence_gamma2`, `coherence_masked_pixels`, `peak_intensity`, wind speed, …) | `code/lumo/lumo_ocv_channels_csv.py` |
| `reference/lumo_tower_coherence_states.json` | committed LUMO result: the 178-overpass coherence table (`bursts`), `per_state`, `per_orbit`, Welch/Mann-Whitney tests, wind control, mask parameters | `code/lumo/lumo_ocv_core.py` (`load_reference`), verified in `pin_all()` |
| `reference/lumo_gamma2_pairs.json` | committed LUMO result: gamma^2 pair comparisons (Cliff's delta, SDS, bootstrap CI), monotonicity/trend test, season/orbit controls, mask controls | same |
| `reference/lumo_echo_mask_vector.json` | committed LUMO result: the per-dimension EMV aggregates (A/D/F/S/P + gamma^2), the ranking by pooled SDS and the meta block | same |
| `reference/lumo_emv_multivariate_cv.json` | committed LUMO result: 4-class LDA, leave-one-out over the shrinkage grid, confusion matrix and label-permutation test | same |
| `reference/lumo_emv_pairwise_triplets_cv.json` | committed LUMO result: out-of-fold scores and SDS/CI/permutation values of the six binary state pairs for the four feature sets | same |

The five JSONs are the *reference* of this repository: the whole point of the
figure script is to recompute their numbers and to fail loudly if any of them
does not come out identically.

### Generated files

| file | written by | content |
| --- | --- | --- |
| `lumo_ocv_channels.csv` | `code/lumo/lumo_ocv_channels_csv.py` | `lumo_channels.csv` plus the mask-geometry columns `A, D, F, S, P` and their ingredients (`mask_rows, mask_cols, bbox_area, centroid_r/c, peak_r/c`) for the 178 coherence overpasses. Only input of the figure pipeline — the burst cache is not needed afterwards. |
| `lumo_ocv_channels_meta.json` | `code/lumo/lumo_ocv_channels_csv.py` | provenance: source CSV + sha256, burst-cache path, number of cache directories matched / skipped, the mask thresholds, the resulting `n_masked` statistics and the cross-check against the committed values |
| `fig_lumo_ocv_paper.json` | `code/lumo/fig_lumo_ocv_paper.py` | the complete result of figures A–E in machine-readable form (descriptives, effect sizes, CV grid, per-pair scores, control tests) plus the summary of the pin block and the versions/seeds of the run |

To rebuild the extended table the burst cache must be present (the generator
looks four levels above `code/lumo/`, i.e. `../../lumo_dam6_analysis/monthly_bursts/`,
overridable with `--burst-dir`). The figures never read that directory.

## `data/bautzen/` — Bautzen OpenLabs bridge

### Inputs (committed copies)

| file | what it is | read by |
| --- | --- | --- |
| `bautzen_rects_asc_iw2.txt`, `bautzen_rects_desc_iw3.txt` | the bridge rect series (ASC / DESC), one line per date: `date\|ORBIT\|rows\|cols\|[[re,im], …]` — the 80 x 80 complex chip of the rect, row-major. ~2 MB each; committed so that a fresh clone rebuilds the whole site offline (the LUMO burst cache is far too large for that). | `code/bautzen/bautzen_ocv_channels_csv.py`, `bautzen_ocv_masks.py` |
| `bautzen_rects_asc_ctrl_500m_east.txt`, `bautzen_rects_ctrl_500m_east.txt` | the same for the control area 500 m east of the bridge (rules out that the deck band is a processing artefact) | same |
| `*.manifest.json` (one per rect file) | per date: `burst_id`, rows/cols, acquisition timestamp, orbit — the time base of the state rule | `code/bautzen/bautzen_ocv_channels_csv.py` |
| `openmeteo_bautzen_2025.json` | hourly air temperature / wind speed / gusts for the site: the environment covariate of the report | `code/bautzen/bautzen_ocv_core.py` (`load_meteo`), figure script |
| `reference/bautzen_deck_edge_mask.json` | committed result of the deck-edge mask analysis: per series `file`, `per_date` (anchor row, peak row/offset, `n_masked`, band contrast, components, centroid), `summary`, the `union50`/`union30` pixel lists and the `config` block of the mask constants | `code/bautzen/bautzen_ocv_core.py` (`load_reference`), verified in `verify()` |
| `reference/bautzen_deck_edge_state.json` | committed result: per series `per_date` (state, band contrast, `n_masked`, the four gamma^2 variants + pixel counts), `n_by_state`, `dates_by_state`, per-state `metrics`, `jaccard` overlaps, plus the site blocks `step_tests`, `placebo`, `delta_vs_ctrl`, `seasonality`, `gamma2_*`, `gamma2_environment` | same |

These two JSONs were produced by the bridge project's own analysis scripts
(`bautzen_deck_edge_*`); this package recomputes their numbers from the rect
files and pins them, so the mirror is checked, not asserted.

### Generated files

| file | written by | content |
| --- | --- | --- |
| `bautzen_ocv_channels.csv` | `code/bautzen/bautzen_ocv_channels_csv.py` | one row per date and series: `gamma2` (deck band, raw and floor-corrected, plus the house variant), the deck-edge channels `A, D, F, S, P`, their ingredients (`anchor_row, peak_row, peak_row_offset, band_contrast, n_masked, n_components, centroid_row/col, mask_rows, mask_cols, bbox_area, row_span, col_span`), the state and the `mask_pixels` column (the masked pixels, row-major — enough to rebuild each per-date mask without the rect files) |
| `bautzen_ocv_channels_meta.json` | same | provenance: the four rect files + sha256, the mask constants (`mask_definition`), the per-series detection summary (`n_dates`, `n_by_state`, `dates_by_state`, `persistence`, `summary`, `union50`/`union30`), the reference sources and the verdict of the verify/cross-check blocks |
| `fig_bautzen_ocv_paper.json` | `code/bautzen/fig_bautzen_ocv_paper.py` (pending) | the complete Bautzen result in machine-readable form plus the pin block summary and the versions/seeds of the run |

The CSV and the meta JSON are both committed: together they are the whole input
of the figure pipeline, so the figure can be regenerated without the rect files
and without any cache.