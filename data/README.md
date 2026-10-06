# `data/` — one subfolder per site: input tables and reference results

`lumo/`, `bautzen/`, `kdlo/` and `carola/`. Everything in here is either a **copy**
of a committed artifact of the site's analysis project or a file **generated** by
a script in `../code/<site>/`; nothing is edited by hand. Every script reads and
writes only inside its own `data/<site>/`, so the sites can never overwrite each
other.

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

## `data/kdlo/` — KDLO media tower (Garden City, SD)

This site is the exception in the repository: **the measurement rows are not
recomputed here** — the two raw files are the committed API dumps of an *external*
analysis record, and that record delivered only the mask *size*
(`coherence_masked_pixels`). Everything else of the echo mask is recomputed from a
**committed strip cache**: `strips/` holds the 64 full-dwell `11 x 400` complex
strips the pipeline decoded (one per measurement of the asset, ~70 KB each,
`u32 width, u32 height` + interleaved little-endian f64), written by
`tower/backend/src/cli/analysis/kdlo_strip_capture.rs`, and the mask geometry
`A, D, F, S, P` plus the masked coherence are derived from them. That is what makes
this package comparable to the LUMO one without inventing a mask.

### Inputs (committed copies)

| file | what it is | read by |
| --- | --- | --- |
| `fig_kdlo_1_raw.json`, `fig_kdlo_2_raw.json` | the two committed raw API dumps: "part 1" (window 2022-06-01..2022-12-23, 17 acquisitions) and "part 2" (2023-05-01..2024-09-30, 43 acquisitions), one row per acquisition with the 27 source columns and the request-id provenance | `code/kdlo/kdlo_ocv_channels_csv.py`, `code/kdlo/kdlo_ocv_core.py` (`load_raw_parts`, verified in `pin_all()`) |
| `kdlo_db_extra.json` | **one-time, read-only** Postgres extraction of the nine columns the API rows do not carry (`phase_snr_db`, `phase_observable`, `phase_unobservable_reason`, `snow_depth_m`, `snowfall_cm`, `precipitation_mm`, `wind_gust_ms`, `intensity`, `incidence_angle_deg`) for the same 60 acquisitions, keyed by row id | same (`load_db_extra`) |
| `reference/fig_kdlo_1_channels.json` | committed result (part 1): window, `n_measurements`, `passes`, `alternate_pass_days`, `orbits`, the `summ` blocks of `coherence_gamma2`, `coherence_masked_pixels`, `phase_coherence`, `phase_snr_db`, `phase_rms_rad` and the four weather columns, `phase_observable`, `echo_modes`, `collapse_pair`, `request_ids` | `code/kdlo/fig_kdlo_ocv_paper.py` (`pin_all`), recomputed exactly |
| `reference/fig_kdlo_2_channels.json` | committed result (part 2): the same blocks plus `intact_reference` (part 1 as its reference) | same |
| `reference/kdlo_s1_availability.json` | committed Sentinel-1 **catalogue query** result: site coordinates, the repeated-orbit cycle, the collapse bracket (`last_before` = 2022-12-11, `first_after` = 2022-12-23), the nullcheck window and the per-part product/burst accounting | `code/kdlo/kdlo_ocv_core.py` (`availability_checks`), **pinned structurally** (128 checks) |
| `reference/kdlo_fem_frequencies.json` | committed **eigenmode-solver** output (structural frequencies of the tower) | `code/kdlo/kdlo_ocv_core.py` (`fem_checks`), **pinned structurally** (67 checks) |
| `strips/capture_requests.json` | the resolved request list of the strip cache: site point, strip spec and per acquisition the CDSE burst Id, orbit direction, relative orbit, sub-swath, polarisation, product name and content start | `code/kdlo/kdlo_capture_requests.py` (writes it), `tower/backend/src/cli/analysis/kdlo_strip_capture.rs` (consumes it) |
| `strips/strip_<date>.bin` | the 64 committed full-dwell strips: `u32 width, u32 height` then `2*w*h` little-endian f64 (re, im), row-major, decoded by `decode_slc_rect_at_location(annotation, measurement, lat, lon, 11, 400)` — the exact window the KDLO pipeline saw | `code/kdlo/kdlo_ocv_masks.py` (mask + geometry), `code/kdlo/kdlo_ocv_channels_csv.py` (attaches the columns) |
| `strips/manifest.json` | provenance of the capture: per acquisition the burst/product ids, sub-swath, polarisation, rows/cols, byte count, sha256, mean intensity and the Rust-side `gamma2`/masked-pixel cross-check | the figure's mask pins and `kdlo_ocv_channels_csv.py --verify` |

The last two files are *not* rederivable from the record — one is a catalogue
query, the other a solver output — so the package pins their internal structure
and the parts of them that follow from the acquisition dates, and says so in the
report instead of pretending to recompute them.

### Generated files

| file | written by | content |
| --- | --- | --- |
| `kdlo_ocv_channels.csv` | `code/kdlo/kdlo_ocv_channels_csv.py` | one row per acquisition: the derived columns `part, state, month, year, season, echo_mode`, the 27 source columns, the 9 vendored DB columns and the 16 echo-mask columns of the 6D vector (`A, D, F, S, P`, the mask geometry, `peak_intensity`, the recomputed `gamma2_mask`). 60 rows: 16 `pre` (2022-06-02..2022-12-11), 43 `rebuild` (2023-05-04..2024-09-19) and the singleton `2022-12-23`, which belongs to no state |
| `kdlo_ocv_channels_meta.json` | same | provenance: the raw files + sha256, the DB-extra file + sha256, the reference files + sha256, the windows, the derived-column definitions, the mask columns, the `strip_cache` block (paths, manifest sha256, per-column mask rule, acceptance criteria) and `n_by_state` / `n_by_part` / `n_by_echo_mode` |
| `fig_kdlo_ocv_paper.json` | `code/kdlo/fig_kdlo_ocv_paper.py` | the complete result of figures A–E in machine-readable form (descriptives, effect sizes, both epoch interpretations, the LOO grid, the permutation verdict, the controls, the rederived channel JSONs) plus the pin summary |

The CSV and the meta JSON are both committed: they are the whole input of the
figure pipeline, so the figures can be regenerated offline from a fresh clone.
Re-extracting the nine DB columns needs the `tower-postgres` container
(`--fetch-db-extra`); that is a maintenance step, never part of a normal run.

## `data/carola/` — Carolabrücke (Dresden), collapsed on 2024-09-11

The collapsed bridge. Like KDLO, this site cannot be rebuilt from a small table:
the analysis project's `carola_windows_full.txt` holds the 80 x 80 complex chip
of every girder request (182 MB, not committed), which is why the mask layer is
carried as a **committed window-payload cache** —
`carola_windows_mask_cache.txt`: the echo mask (pixel coordinates) of every chip
plus the sufficient statistics the 6D vector needs. Everything else is derived
from it here. Unlike KDLO, the record *does* have a per-girder resolution: three
girders (Girder A/B/C, 1,189 + 264 + 264 chips) and two mask layers (the echo
mask of the whole chip and the same rule on the girder's own chip).

### Inputs (committed copies)

| file | what it is | read by |
| --- | --- | --- |
| `carola_windows_mask_cache.txt` | the committed mask cache: one line per chip (`mid`, asset, segment, date, orbit, `w h`, `peak`, `median`, the mask statistics and the masked-pixel coordinates), written once by `carola_ocv_masks.py --extract` from the 182 MB payload file and re-verified by `--check` | `code/carola/carola_ocv_masks.py` (`load_cache`), `code/carola/carola_ocv_channels_csv.py` |
| `carola_windows_mask_cache.manifest.json` | provenance of that extraction: source file + sha256, window sizes (1,717 `80x80`, 265 `7x7` excluded), the mask rule, the 472 de-duplicated byte-identical payloads, the 264 asset/date groups with several girder segments and the resolution flag | the figure's mask pins and `carola_ocv_masks.py --check` |
| `carola_measurements_full.txt` | the committed measurement extract, one line per acquisition (`id`, asset, `request_id`, `segment_index`, timestamp, orbit, pass, wind, temperature, precipitation, gust, humidity, incidence angle, intensity, the pipeline's `coherence_masked_pixels` / `coherence_gamma2`, …) | `code/carola/carola_ocv_channels_csv.py`, `code/carola/carola_ocv_core.py` (`load_measurements`) |
| `carola_segments.txt` | the girder labels + FEM baselines: `(asset_id, segment_index) -> {label, fundamental_hz}` — Girder A / B / C and the `Carola Brücke` tower label | `code/carola/carola_ocv_masks.py` (segment resolution), `code/carola/carola_ocv_core.py` (`load_segments`) |
| `carola_windows_index.txt` | the export index: per chip the asset, `request_id`, `segment_index`, content start and asset-point coordinates — the per-girder re-extraction that makes the segment resolvable | `code/carola/carola_ocv_masks.py` (`extract`) |
| `reference/carola_coherence_states.json` | committed result (first analysis): the 1,454 stored windows (`records` with `date`, orbit, wind, temperature, `gamma2`, `n_masked`), `per_state`, `per_orbit`, `per_girder`, `per_pass`, `bursts`, the monthly table and the wind control | `code/carola/fig_carola_ocv_paper.py` (`pin_all`), recomputed per window |
| `reference/carola_echo_mask_gamma2.json` | committed result (echo-mask analysis): the state table (`healthy` / `damaged` / `pooled` / `by_orbit` / `committed_reference`), `models` (nested OLS + inverse-N + the two F-tests), `sds`, `by_segment`, `by_request`, `by_asset` and the guard/meta block | same, **recomputed number by number** |
| `reference/carola_bridge_osm.json` | committed Overpass API extract of the bridge: two `bridge=yes` ways (Carolabrücke, B 170, 12 nodes each) and the two razed DVB tram ways of the collapse | same, **pinned structurally** (it is a map extract, not rederivable) |

The `coherence_masked_pixels` column of the measurement extract is **not** the
mask of these chips: the pipeline computed it on the pre-fix *asset-point* chip
while the committed export holds the per-girder re-extraction, so it coincides
with the recomputed `A` on only **60 of 1,713** rows. The generator keeps it as
`coherence_masked_pixels`/`gamma2_export` provenance and reports the difference.

### Generated files

| file | written by | content |
| --- | --- | --- |
| `carola_ocv_channels.csv` | `code/carola/carola_ocv_channels_csv.py` | one row per chip (1,717 rows, 57 columns): the source columns, the derived columns (`segment`, `girder_index`, `state`, `pre_collapse`, `day`, `month`, `year`, `season`, `echo_mode`) and the 18 mask columns of the 6D vector (`A, D, F, S, P, gamma2_mask`, the bbox/centroid/peak geometry, `n_components`, the largest component, the fragmentation) |
| `carola_ocv_channels_meta.json` | same | provenance: the cache, the manifest, the measurements and the segments file + sha256, the mask rule and window, the girder counts, the state definition, `n_echo_masks`, `n_rows_without_echo_mask`, the `coherence_masked_pixels` cross-check and the per-state summary of every channel |
| `fig_carola_ocv_paper.json` | `code/carola/fig_carola_ocv_paper.py` | the complete result of figures A–E in machine-readable form (channel and girder descriptives/effect sizes, the strata, the nested OLS, the LDA grid and OOF scores, the mask-size scatter block, the state summaries) plus the pin summary and the seeds of the run |

The CSV and the meta JSON are both committed: they are the whole input of the
figure pipeline, so the figures and the pin block can be regenerated offline from
a fresh clone. `carola_ocv_masks.py --extract` (needs the 182 MB payload file) is
the only step that must not be part of a normal run.

