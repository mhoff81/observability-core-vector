# `data/` — one subfolder per site: input tables and reference results

`lumo/`, `bautzen/`, `kdlo/`, `carola/`, `morandi/`, `cts/`, `ywf/` and `espoo/`. Everything in here is either a **copy**
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


## `data/morandi/` — Ponte Morandi / Polcevera (Genova), collapsed on 2018-08-14

The collapsed viaduct. Like Carola, this site cannot be rebuilt from a small
table: the site project's per-date **400 x 400 rect payloads** (in
`/home/projects/morandi_analysis`, not committed) hold the complex chip of every
acquisition, which is why the mask layer is carried as a **committed
window-payload cache** — `morandi_windows_mask_cache.txt`: the echo mask (pixel
coordinates) of every chip plus the sufficient statistics the 6D vector needs.
The payloads are cropped to the central `[160:240, 160:240]` 80 x 80 window
before the rule is applied. Where Carola has per-girder chips, this site has the
two **Sentinel-1 tracks** of the deck (`A_asc` ASC rel-15 IW1 VV, `A_des` DESC
rel-66 IW1 VV, `morandi_tracks.txt`) and two mask layers (the echo mask of the
whole chip and the same rule on the track's own chip).

The record is tiny after the event: **232** rows, **221** pre / **11** post.

### Inputs (committed copies)

| file | what it is | read by |
| --- | --- | --- |
| `morandi_windows_mask_cache.txt` | the committed mask cache: one line per chip (`mid`, track, date, orbit, `w h`, `peak`, `median`, the mask statistics and the masked-pixel coordinates), written once by `morandi_ocv_masks.py --extract` from the site's 400 x 400 rects and re-verified by `--check` | `code/morandi/morandi_ocv_masks.py` (`load_cache`), `code/morandi/morandi_ocv_channels_csv.py` |
| `morandi_windows_mask_cache.manifest.json` | provenance of that extraction: the two source rect files + sha256 and track, the `400` source size, the `[160,240,160,240]` crop, the window sizes (234 `400x400`, 5 clipped/odd `200x400`/`268x400` excluded), the mask rule and the per-track counts (119 / 113) | the figure's mask pins and `morandi_ocv_masks.py --check` |
| `morandi_weather_table.csv` | the committed measurement + weather extract, one pipe-delimited line per acquisition (`track`, `orbit`, `date`, `platform`, the wind/temperature/precipitation/humidity radiation columns, `intensity`, `gamma2`, the phases, `coh_masked`, `is_master`, …) | `code/morandi/morandi_ocv_masks.py` (`extract`), `code/morandi/morandi_ocv_channels_csv.py`, `code/morandi/morandi_ocv_core.py` (`load_measurements`) |
| `morandi_tracks.txt` | the track table: `track -> {label, orbit, rel_orbit, platform, rect/npz file, first, last, master, n_dates, n_pre, n_post, n_masked_pixels}` — `A_asc` 119 dates, `A_des` 115, shipped `n_masked_pixels` 640 | `code/morandi/morandi_ocv_core.py` (`load_tracks`) |
| `morandi_deck_channels.csv` | the site's own per-track deck channels: a **fixed 640-px (10 %) quantile mask**, *not* the echo mask this package recomputes | the figure's deck-mask pins and the export cross-check |
| `reference/morandi_registration.json` | committed registration/coregistration result: master dates, coregistration and the `n_masked_pixels = 640` deck mask | `code/morandi/fig_morandi_ocv_paper.py` (`pin_block`), pinned structurally |
| `reference/morandi_deck_channels.json` | committed deck-channel result: the fixed 640-px mask table | same |
| `reference/morandi_geometry.json` | committed viaduct geometry (Ponte Morandi / Polcevera, Genova) | same, **pinned structurally** (a geometry extract, not rederivable) |

The `coh_masked` / registration `n_masked_pixels` column is **not** the mask of
these chips: it is a *fixed 640-px (10 %) quantile mask* while the committed
cache holds the recomputed echo mask, so it coincides with the recomputed `A` on
**0 of 232** rows. The generator keeps it as provenance and reports the
difference (panel E).

### Generated files

| file | written by | content |
| --- | --- | --- |
| `morandi_ocv_channels.csv` | `code/morandi/morandi_ocv_channels_csv.py` | one row per chip (**232 rows, 59 columns**): the source columns, the derived columns (`segment`, `track_index`, `state`, `pre_collapse`, `day`, `month`, `year`, `season`, `echo_mode`) and the mask columns of the 6D vector (`A, D, F, S, P, gamma2_mask`, the bbox/centroid/peak geometry, `n_components`, the largest component, the fragmentation) |
| `morandi_ocv_channels_meta.json` | same | provenance: the cache, the manifest, the weather table and the tracks file + sha256, the mask rule and window, the track counts, the state definition, `n_echo_masks`, `n_rows_without_echo_mask`, the `coherence_masked_pixels` cross-check (232 compared, 0 matched) and the per-state summary of every channel |
| `fig_morandi_ocv_paper.json` | `code/morandi/fig_morandi_ocv_paper.py` | the complete result of figures A–E in machine-readable form (channel and track descriptives/effect sizes, the strata, the weather controls, the LDA grid and OOF scores, the state summaries) plus the pin summary and the seeds of the run |
| `fig_morandi_ocv_months.json` | `code/morandi/fig_morandi_ocv_months.py` | the companion's per-window blocks of the four one-month windows (counts, medians, Cliff's deltas) and the 6D fingerprint |

The CSV and the meta JSON are both committed: they are the whole input of the
figure pipeline, so the figures and the pin block can be regenerated offline from
a fresh clone. `morandi_ocv_masks.py --extract` (needs the site's 400 x 400 rect
files) is the only step that must not be part of a normal run.

## `data/cts/` — Champlain Towers South (Surfside), collapsed on 2021-06-24

The collapsed condominium. The record is a **four-footprint chip export of one
orbit**: per acquisition the analysis project stored one 80 x 80 complex chip per
Sentinel-1 footprint of the *same* geometry (ASC rel-48 IW3) — the tower
(`cts`), two standing control towers (`ctn`, 165 m N; `cte`, ~90 m N) and the
beach (`beach`) — see `cts_tracks.txt`. The payloads (the site project's
`cts_windows_*.npz` cubes, not committed) are carried as a **committed window
cache** — `cts_windows_mask_cache.txt`: the echo mask (pixel coordinates) of
every chip plus the sufficient statistics the 6D vector needs. Where Morandi has
two tracks, this site has four, and two mask layers (the echo mask of the whole
chip and the same rule on the track's own chip).

The four footprints share **identical strata, dates and counts**, so this is the
one site that can do a **difference-in-differences** (DiD): the collapsed tower
against a control tower whose common atmosphere cancels. The record is **711**
rows, **543** pre / **168** post (one `cte` date has no echo).

### Inputs (committed copies)

| file | what it is | read by |
| --- | --- | --- |
| `cts_windows_mask_cache.txt` | the committed mask cache: one line per chip (`mid`, track, segment, date, orbit, `w h`, `peak`, `median`, the mask statistics and the masked-pixel coordinates), written once by `cts_ocv_masks.py --extract` from the `cts_windows_*.npz` cubes and re-verified by `--check` | `code/cts/cts_ocv_masks.py` (`load_cache`), `code/cts/cts_ocv_channels_csv.py` |
| `cts_windows_mask_cache.manifest.json` | provenance of that extraction: the four source npz files + sha256 and track, the `80` window, the mask rule, the per-track counts (178/178/177/178), the one row without echo and the site's `quantile_mask_pixels` (640) | `cts_ocv_masks.py --check`, the channel meta | 
| `cts_measurements_full.txt` | the committed measurement extract, one pipe-delimited line per acquisition per track (712 rows, 28 columns: `asset_id`…`track`, `intensity`, `coherence_gamma2`, `coherence`, `phase_coherence`, `coherence_masked_pixels`, `phase_scatterer_count`, `brightness_ratio`, `phase_snr_db`, `displacement_los_m`, the weather columns, the registration shifts) | `code/cts/cts_ocv_core.py` (`load_measurements`), `cts_ocv_channels_csv.py`, `cts_ocv_did.py` |
| `cts_tracks.txt` | the track table: `track -> {label, role, orbit_direction, rel_orbit, subswath, platform, npz_file, first, last, n_dates, n_pre, n_post, n_masked_pixels}` (178 dates each, role `target` / `control_ctn` / `control_cte` / `reference_beach`) | `code/cts/cts_ocv_core.py` (`load_tracks`) |
| `cts_segments.txt` | provenance only: a "segment" is a track here, so `SEGMENTS_PATH` aliases `cts_tracks.txt` | `code/cts/cts_ocv_core.py` (`load_segments`) |
| `cts_asset.json`, `cts_events.json`, `cts_coherence_states.json`, `cts_acquisition_census.{json,md}`, `cts_phase0_audit.{json,md}` | the site's Phase 0 provenance chain (asset, the collapse/demolition events, the coherence-state histogram, the acquisition census and the Step-0 audit) | `code/cts/cts_ocv_core.py` (`load_reference`), the reports |
| `reference/cts_reference_did_ctn.json` (`.md`) | committed DiD reference: `target` vs `control_ctn` | `code/cts/cts_ocv_did.py` (`--verify`), **pinned series + markdown** |
| `reference/cts_reference_did_cte.json` (`.md`) | committed DiD reference: `target` vs `control_cte` | same |
| `reference/cts_reference_did_placebo_ctn_vs_cte.json` (`.md`) | committed DiD reference (placebo): `control_ctn` vs `control_cte` | same |

The site's own `coherence_masked_pixels` column is **not** the mask of these
chips: it is a *fixed 640-px (10 %) quantile mask* while the committed cache
holds the recomputed echo mask. The generator keeps it as provenance, and in the
DiD it is the **degenerate** channel (constant per row, so its OLS interaction is
pure floating-point noise `-5.07e-14`, pinned structurally only).

### Generated files

| file | written by | content |
| --- | --- | --- |
| `cts_ocv_channels.csv` | `code/cts/cts_ocv_channels_csv.py` | one row per chip (**711 rows**; the header has 54 entries — `coherence_masked_pixels` is carried once in the source block and once in the mask block, 53 distinct names): the source columns, the derived columns (`segment`, `track_index`, `state`, `pre_collapse`, `day`, `month`, `year`, `season`, `echo_mode`) and the mask columns of the 6D vector (`A, D, F, S, P, gamma2_mask`, the bbox/centroid/peak geometry, `n_components`, the largest component, the fragmentation) |
| `cts_ocv_channels_meta.json` | same | provenance: the cache, the manifest, the measurements and the tracks file + sha256, the mask rule and window, the track/role counts, the state definition, `n_echo_masks`, `n_rows_without_echo_mask`, the `coherence_masked_pixels` cross-check and the per-state summary of every channel |
| `fig_cts_ocv_paper.json` | `code/cts/fig_cts_ocv_paper.py` | the OCV-paper figure result: the pooled and per-track Cliff's-delta blocks, the fitted 2-class LDA models, the pinned ASCENDING `did_asc` rows of the three contrasts and the pin verdict (the figure script re-derives the three committed references and exits non-zero on any deviation) |
| `fig_cts_ocv_months.json` | `code/cts/fig_cts_ocv_months.py` | the companion's per-window blocks (the five windows `F0` + `M-3..M+1` with counts, medians and dates; the Cliff's deltas and bootstrap CIs of `M+1` vs. `F0` and vs. `M-1`, plus every dimension vs. the pooled pre windows; the 6D fingerprint relative to the `F0` first-images baseline) |

The CSV and the meta JSON are both committed: they are the whole input of the
DiD layer, so the three reference contrasts can be regenerated offline from a
fresh clone. `cts_ocv_masks.py --extract` (needs the site's `cts_windows_*.npz`
cubes) is the only step that must not be part of a normal run.

## `data/ywf/` — Yeongdeok Wind Farm, Unit 21 (Samgye-ri), collapsed on 2026-02-02

The steel monopole tower (100 m, 4.6 m base diameter) that collapsed onto a public
road. The record is a **7 x 7 complex window** per acquisition per *mast-section
request* — and it is **small** (57 windows x 49 complex samples = 38 kB), so unlike
the bridge sites the payload is **committed in full** and the mask layer is
recomputed from it directly. It is also **not segment-resolved**: 24 of the 57
windows are byte-identical duplicates that differ only in the requested section
(2 vs. 3) on one date, so both requests decoded the same chip at the asset point
— the package therefore carries **one** mask layer (the echo mask of the whole
window; a per-mast-section layer would be an artefact) and only mast section 3
appears in the record at all.

The record holds **33 unique chips** on **32 dates**; **17** of them carry an echo
(**16** of 27 pre-event chips, **1** of 6 post-event chips). The site's collapsed
tower therefore supports an **echo-coverage** statement with a power caveat
rather than a contrast. The same coverage is also stated per **orbit direction**
(ascending / descending, the record's only second cross-cut beside the state):
**10 of 20** ascending chips carry an echo against **7 of 13** descending ones
(Fisher exact p = 1.0000), the one post-event echo is an ascending pass, and the
block is pinned as a *description* — both orbits straddle the collapse date.

Two facts the export does not carry are committed as references resolved from
public catalogues (`reference/ywf_bursts.json`, `reference/ywf_osm_road.json` and
`reference/ywf_geometry.json`, each with its own offline `--check`): the record's
33 chips are **two sub-swaths** — ascending `IW3` / relative orbit `54` at 09:2x UT
(20 chips) and descending `IW2` / relative orbit `61` at 21:2x UT (13 chips), with
**VV and VH interleaved in both** — and the collapsed tower stood **on** OSM way
`577440814` (`highway=tertiary`, `surface=asphalt`). The road is an anchor only:
a ~7 m carriageway is 2.0 px across range but 0.5 px across azimuth, and the
window is 7 x 7 px with no pixel-spacing spec of its own, so this record supports
no edge/band layer (see `code/ywf/README.md`, "the four structural blockers").

### Inputs (committed copies)

| file | what it is | read by |
| --- | --- | --- |
| `ywf_windows_full.txt` | the committed window payload: one line per window, `mid\|w\|h\|[[re, im], …]` (57 lines, 38 kB) | `code/ywf/ywf_ocv_masks.py`, `ywf_ocv_core.py` (`chip_table`), `ywf_ocv_channels_csv.py` |
| `ywf_windows_index.txt` | the window index: `mid\|asset_id\|asset_name\|segment_index\|acquisition_ts\|burst_id\|w\|h` — the 6th field is the **CDSE burst id** of the SLC burst the window was cut from (**not** the measurement extract's own `request_id`, a different column with a different value); the index labels a window's mast section and is the time base of the state rule | same |
| `ywf_segments.txt` | the five mast sections of the asset with their FEM fundamental: `asset_id\|asset_name\|segment_index\|label\|fundamental_hz\|…` | `ywf_ocv_core.py` (`load_segments`) |
| `ywf_measurements_full.txt` | the committed measurement extract: one pipe-delimited line per measurement (58 lines, ~50 columns: `intensity`, `displacement_los_m`, `brightness_ratio`, the weather columns, `coherence`, `coherence_gamma2`, `coherence_masked_pixels`, `structural_frequency_hz`, …) | `ywf_ocv_core.py` (`load_measurements`), `ywf_ocv_channels_csv.py`, the figure script |
| `ywf_windows_mask_cache.txt` | the committed mask cache — **only the 17 windows that carry an echo** (one line per chip: `mid`, asset, segment, date, orbit, `w h`, `peak`, `median`, the mask statistics and the masked-pixel coordinates), written by `ywf_ocv_masks.py --extract` and re-verified by `--check` | `ywf_ocv_masks.py` (`load_cache`), `ywf_ocv_core.py`, the figure script |
| `ywf_windows_mask_cache.manifest.json` | provenance of that extraction: the three source files + sha256, the mask rule, the guard counts (`n_rows` 17, `n_without_echo` 16, `n_deduplicated_same_payload` 24, `segment_resolved` false, the section groups), written by the same run | `ywf_ocv_masks.py --check`, `ywf_ocv_core.py`, the channel meta |
| `reference/ywf_ocv_findings.json` | the committed **analysis-layer reference** — the only reference of the repository that has *no* upstream project behind it: it is **defined** by `fig_ywf_ocv_paper.py --write-reference` and re-derived field by field (1179 checks) on every later run | `code/ywf/fig_ywf_ocv_paper.py` (`pin_all`), `fig_ywf_ocv_weeks.py` |
| `reference/ywf_bursts.json` | the **resolved pass geometry** of the record (49 kB): all 33 burst ids of the window index queried against the public CDSE OData burst catalogue — sub-swath, relative orbit, polarisation, azimuth time, parent product, footprint, S3 path, plus the sha256 of the index it was resolved from. A **mutable external input**, so it is re-verified *structurally* offline (`code/ywf/ywf_bursts_resolve.py --check`, the default) and pinned by sha256 in the figures | `ywf_bursts_resolve.py`, `code/ywf/ywf_ocv_core.py` (`burst_provenance`) |
| `reference/ywf_osm_road.json` | the committed **OSM map extract** of the road the collapsed tower fell on (91 kB, way `577440814`, `highway=tertiary`, © OpenStreetMap contributors, ODbL 1.0) | `code/ywf/ywf_osm_anchor.py` |
| `reference/ywf_geometry.json` | the derived **anchor** (its `status` is `anchor reference — no edge layer claimed`): the nearest point of that way (36.85 m from the asset point), the carriageway centreline and bearing, the **estimated burst grid** (`chip_grid`: 3.35 / 3.45 m per pixel across range, 15.3 / 14.4 across azimuth), `road_vs_pixel_grid` — a 7 m carriageway is 2.0-2.1 px across range but 0.46-0.49 px across azimuth — and `road_vs_window`, the **containment proof**: the near carriageway edge is 8.10-9.33 px across range from the window centre while a 7 x 7 frame reaches 6 px, so no committed pixel crosses the road (containment would need 19-21 px per side, 27-29 px with Bautzen's ±4-row ROI) | `code/ywf/ywf_osm_anchor.py` (`--check`, the default) |

### Generated files

| file | written by | content |
| --- | --- | --- |
| `ywf_ocv_channels.csv` | `code/ywf/ywf_ocv_channels_csv.py` | one row per **unique chip** (**33 rows**, 81 columns): the 34 source columns of the measurement extract, the derived columns (`segment`, `section_index`, `state`, `pre_collapse`, `day`, `month`, `year`, `season`, `echo_mode`, `echo_mask_present`, `payload_group_size`, `payload_shared_with_segments`) and the recomputed mask columns (`A, D, F, S, P, gamma2_mask`, the bbox/centroid/peak geometry, `n_components`, the largest component, the fragmentation). `echo_mask_present` marks the 17 echo-bearing rows; the mask columns are empty on the other 16 |
| `ywf_ocv_channels_meta.json` | same | provenance: the payload, index, segments, measurements, cache and manifest + sha256, the mask rule and window, the guard counts (`n_windows` 57, `n_unique_chips` 33, `n_duplicate_windows` 24, `n_echo_masks` 17, `n_rows_without_echo_mask` 16, `n_by_orbit` ASC 20 / DESC 13), the echo-coverage table with its Fisher test, **the same coverage split by orbit direction** (`coverage_by_orbit`: per-orbit chips / echoes / pre-post cross-tab, the ASC-vs-DESC Fisher test and the two within-orbit pre/post tests), the echo-mode cuts (cut on this site's 49-px window), the `coherence_masked_pixels` cross-check and the per-state summary of every channel |
| `fig_ywf_ocv_paper.json` | `code/ywf/fig_ywf_ocv_paper.py` | the OCV-paper figure result: the coverage/power block, the per-orbit coverage block, the pooled Cliff's-delta blocks of the mask channels and of the co-variate controls, the echo-mode mix, the month/section composition, the 32-row timeline, the payload-digest groups and the pin verdict |

Neither the CSV nor the meta JSON contains anything that is not derivable from
the committed payload, so the figures and the pin block run offline from a fresh
clone; `ywf_ocv_masks.py --extract` is the only step that rewrites the cache (and
it is deterministic, byte for byte).

## `data/espoo/` — Espoo Kurttila mast (Espoo, Finland), event-free

The communications mast at Kurttila, whose record is the family's **event-free
control**. It is the only site in this repository with **no damage axis at all**:
the committed export behind it (`espoo_channels.csv`) carries 150 acquisitions between
`2024-09-05` and `2026-09-07` (25 months at exactly 6 each) and `damage_label`,
`structural_state` and `condition_label` are empty on **every** row — the upstream
verdict says so itself ("This site has no ground-truth state label"). There is no
event to anchor a window on and no severity to regress against, so the only contrast
the record supports is the one between its two **orbit geometries** (`ASCENDING` /
afternoon pass n = 70, `DESCENDING` / morning pass n = 80).

The export also has **no mask raster**: `peak_intensity`, `sub_aperture_brightness`,
`scatterers`, `amplitude`, `phase_rad`, `mast_peak_row/col` and every other geometry
column are empty on all 150 rows (56 of the 78 columns). The observability core
vector of the site is therefore `OCV = [gamma2, A]` with `A` =
`coherence_masked_pixels`, the mask *size* the pipeline delivered: `D`, `F`, `S` and
`P` are deliberately absent, and `A` is **copied and proved**, never recomputed. The
phase ladder (`phase_coherence`, `phase_rms_rad`) exists for 80 of 150 acquisitions
(53.3 %, nothing after 2025-10) and `phase_snr_db` for 0, so it is stored and
reported as *availability* only.

### Inputs (committed copies)

| file | what it is | read by |
| --- | --- | --- |
| `espoo_channels.csv` | the committed Espoo mast export as the site's own generator wrote it: one row per acquisition (**150 rows, 78 columns** — `acquisition_ts`, `pass_label`, `orbit_direction`, `coherence_gamma2`, `coherence_masked_pixels` (= `A`), the phase-ladder columns, `wind_speed_ms`, `temperature_c`, `ndvi`, `status` and the empty geometry columns) | `code/espoo/espoo_ocv_channels_csv.py`, `espoo_ocv_core.py` | 
| `reference/espoo_mast_observability.json` | the committed **analysis-layer reference** of the site: `per_pass`, `per_orbit`, `orbit_test` (Welch / Mann-Whitney), `verdict.classes` ("no ground-truth state label"), `constants`, the LUMO reference scale and the provenance / catalogue blocks | `code/espoo/espoo_ocv_channels_csv.py`, `fig_espoo_ocv_paper.py` (`pin_all`) |
| `reference/espoo_phase_stationarity.json` | the committed **stationarity table**: all 14 phase-coherence series (`1-3 / 3-8 / 8-18 Hz` x `full_burst / strip400_annotation / strip400_pipeline` x ASCENDING / DESCENDING) with median, IQR, Theil-Sen slope + moving-block bootstrap CI, Mann-Kendall, lag-1, runs test, Ljung-Box, half-split shift, annual harmonic, the 25 monthly means, the flags and the verdict — plus every number the file's own prose quotes | `code/espoo/fig_espoo_ocv_paper.py`, `fig_espoo_ocv_months.py` |
| `reference/fig_espoo_channels.json` | the committed **channel statistics**: the 15 pooled / per-orbit `{n, median}` blocks, the 25-month composition and the LUMO reference scale | `code/espoo/espoo_ocv_channels_csv.py`, `fig_espoo_ocv_paper.py` |
| `reference/espoo_ocv_findings.json` | this package's **own derivation lock** — *defined* by `code/espoo/fig_espoo_ocv_paper.py --write-reference` once (there is no upstream file behind it) and re-derived field by field on every later run | `code/espoo/fig_espoo_ocv_paper.py` (`pin_all`) |

The first four files are the *reference* of this package: the whole point of the
generator and the figure scripts is to recompute their numbers and to fail loudly if
any of them does not come out identically.

### Generated files

| file | written by | content |
| --- | --- | --- |
| `espoo_ocv_channels.csv` | `code/espoo/espoo_ocv_channels_csv.py` | `espoo_channels.csv` plus the derived bookkeeping columns (state = the orbit geometry, month, season, the acquisition order). **`A` is copied, not recomputed**: the generator proves the copy against the committed column *and* against the committed reference statistics row by row (4050 copied fields, 235 reference channel statistics, `A_equals_committed_masked_pixels`, 0 rows left empty, 25 months) |
| `espoo_ocv_channels_meta.json` | same | provenance: source and output sha256, `mask_definition` (why there is no mask layer), `derived_checks`, `filled_columns`, `new_columns`, `n_by_orbit` (ASC 70 / DESC 80), `date_range` and the `monthly` 25-month composition |
| `fig_espoo_ocv_paper.json` | `code/espoo/fig_espoo_ocv_paper.py` | the OCV-paper figure result: the pooled / per-orbit channel blocks, `effects`, `models` (LDA + shrinkage grid + permutation null), `oof`, `oot`, the four `controls` blocks (stationarity, season, wind strata, phase ladder), `headline`, `shrinkage` and the pin verdict |
| `fig_espoo_ocv_months.json` | `code/espoo/fig_espoo_ocv_months.py` | the month companion's result: the 25 `months` / `per_month` composition, `deltas` (the per-month Cliff's delta + CI), `pooled`, `primary` (the two `3-8 Hz · strip400_pipeline` series), `availability`, `signs`, `standing` and its own pin verdict |

Neither the channel table nor the meta JSON contains anything that is not derivable
from the committed `espoo_channels.csv` and the committed reference JSONs, so the
figures and every pin block run offline from a fresh clone. Two runs of the
generator produce the same CSV and the same meta JSON, byte for byte.

