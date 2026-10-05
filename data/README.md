# `data/` — input tables and reference results

Everything in here is either a **copy** of a committed LUMO artifact or a file
**generated** by a script in `../code/`. No file in this folder is edited by
hand.

### Inputs (committed copies)

| file | what it is | read by |
| --- | --- | --- |
| `lumo_channels.csv` | the LUMO overpass table as committed in the LUMO project: one row per Sentinel-1 overpass of the tower (date, orbit, acquisition mode, `gamma2`, `n_masked`, `peak_intensity`, wind speed, …) | `code/lumo_ocv_channels_csv.py` |
| `reference/lumo_tower_coherence_states.json` | committed LUMO result: the 178-overpass coherence table (`bursts`), `per_state`, `per_orbit`, Welch/Mann-Whitney tests, wind control, mask parameters | `code/lumo_ocv_core.py` (`load_reference`), verified in `pin_all()` |
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
| `lumo_ocv_channels.csv` | `code/lumo_ocv_channels_csv.py` | `lumo_channels.csv` plus the mask-geometry columns `A, D, F, S, P` and their ingredients (`mask_rows, mask_cols, bbox_area, centroid_r/c, peak_r/c`) for the 178 coherence overpasses. Only input of the figure pipeline — the burst cache is not needed afterwards. |
| `lumo_ocv_channels_meta.json` | `code/lumo_ocv_channels_csv.py` | provenance: source CSV + sha256, burst-cache path, number of cache directories matched / skipped, the mask thresholds, the resulting `n_masked` statistics and the cross-check against the committed values |
| `fig_lumo_ocv_paper.json` | `code/fig_lumo_ocv_paper.py` | the complete result of figures A–E in machine-readable form (descriptives, effect sizes, CV grid, per-pair scores, control tests) plus the summary of the pin block and the versions/seeds of the run |

To rebuild the extended table the burst cache must be present (the generator
looks for `../../lumo_dam6_analysis/monthly_bursts/`, overridable with
`--burst-dir`). The figures never read that directory.
