# `code/ywf/` — Yeongdeok Wind Farm, the collapsed Unit 21 monopole

**Site** Yeongdeok Wind Farm / Changpo Wind Power Complex, Samgye-ri
(36.4236 N 129.4203 E) · **asset** `Yeongdeok Wind Turbine (Samgye-ri)`
(`af683c7e-416a-4395-86cc-094ee83d9497`) · **event** the tower **collapsed onto a
public road on 2026-02-02** · **states** `pre-collapse (healthy)` (Aug 2025 -
Jan 2026) and `post-collapse` (Feb - Apr 2026) of the *same* structure.

This is the seventh package of the family (after `lumo/`, `bautzen/`, `kdlo/`,
`carola/`, `morandi/`, `cts/`) and it keeps their shape: a stats layer, a data
layer, a mask layer, a channel-table generator, a figure/report script with a
pin block, and a committed reference under `data/ywf/reference/`. What differs is
the *size* of the record and, with it, what the package is allowed to claim.

```bash
python3 code/ywf/ywf_ocv_masks.py --check     # re-verify the mask rule on the cache
python3 code/ywf/ywf_ocv_masks.py --extract   # rewrite data/ywf/ywf_windows_mask_cache.txt
python3 code/ywf/ywf_ocv_channels_csv.py      # build + verify the channel table
python3 code/ywf/ywf_ocv_core.py              # read-only brief of the committed layer
python3 code/ywf/ywf_ocv_core.py --chips      # one line per unique chip, with its 6D vector
python3 code/ywf/ywf_bursts_resolve.py        # --check the pass geometry of the record (default)
python3 code/ywf/ywf_osm_anchor.py            # --check the road anchor + the estimated grid (default)
python3 code/ywf/ywf_mask_sensitivity.py      # the offline mask-rule sensitivity (writes nothing)
python3 code/ywf/ywf_bursts_resolve.py --write   # re-resolve against the live CDSE catalogue
python3 code/ywf/ywf_osm_anchor.py --fetch --write   # re-fetch the OSM way and re-derive the anchor
bash    figures/ywf/fig_ywf_ocv_paper.sh      # figures A-E + report + pins (~8 s)
bash    figures/ywf/fig_ywf_ocv_paper.sh --write-reference   # (re)define the reference
```

## The record

| fact | value |
| --- | --- |
| committed windows | **57** on **32** acquisition dates (`data/ywf/ywf_windows_full.txt`, 38 kB, committed **in full**) |
| unique chips | **33** — 24 of the 57 windows are **byte-identical duplicates** that differ only in the requested mast section (index 2 vs. 3, same date, same payload) |
| `segment_resolved` | **false** — the two requests decoded the *same* chip at the asset point |
| mast sections present | `turbine-mast-section-3` only (of the five in `data/ywf/ywf_segments.txt`) |
| chips with an echo | **17** of 33 — **16** of 27 pre-event chips, **1** of 6 post-event chips |
| window size | 7 x 7 px (49 px), one complex window per request |

Because the record is not segment-resolved, this package carries **one** mask
layer — the echo mask of the whole 7 x 7 window. A per-mast-section (per-girder)
layer, the second layer of Carola/Morandi/CTS, would be an artefact of a
duplication and is deliberately absent; figure E d shows the duplication.

The one cross-cut the record *does* carry beside the state is the **sub-swath** of
the pass geometry (`orbit_direction`: **ASCENDING** / **DESCENDING**), and it is
not a plain pass split: the ascending pass is **IW3 on relative orbit 54**
(09:2x UT, 20 of the 33 chips), the descending one **IW2 on relative orbit 61**
(21:2x UT, 13 chips) — the reverse of Bautzen, whose ascending pass is IW2 — and
the record interleaves **VV (14 chips) and VH (19)** in both
(`data/ywf/reference/ywf_bursts.json`). The export's own `pass_label` follows the
*UTC* clock of `acquisition_ts` (`morning` = 09:2x UT, `afternoon` = 21:2x UT);
at 129.4 E the ascending pass is the local **evening** one (18:00 local solar
time) and the descending one the local morning (06:00), so `pass_label` is an
orbit clue and reads *inverted* as a statement about local time or illumination.

It is reported as a second, explicitly descriptive coverage block —
`coverage_by_orbit` in `ywf_ocv_core.py`, the `echo_coverage_by_orbit` block of
the reference and §A2 of the report — never as an orbit contrast: both orbits
straddle the collapse date, the two echo rates differ by 3.8 percentage points,
and the single post-event echo is an ascending pass.

## The mask layer

`data/ywf/ywf_windows_full.txt` is the committed payload (this site is small
enough to commit it), and `data/ywf/ywf_windows_mask_cache.txt` is the committed
mask cache the other packages have, written by `ywf_ocv_masks.py --extract`:

```
intensity = |z|^2
mask      = intensity >= 0.30 * max(intensity)  AND  intensity >= 5.0 * np.median(intensity)
A         = count(mask)                    (no echo when A < 2)
gamma2    = |sum(z[mask])|^2 / (A * sum(|z[mask]|^2)),  clamped to [0, 1]
D         = A / bbox_area          F = 8-connected components      S = ||centroid - peak||
P         = mean pixel frequency inside the own mask
```

The cache holds **only the 17 windows that carry an echo** (8 kB); its manifest
records the 16 dropped ones (`n_without_echo`), the 24 de-duplicated payloads and
the section counts. `--check` re-derives every dimension of every cache row from
the committed payload and fails on the first mismatch.

## The channel table

`data/ywf/ywf_ocv_channels.csv` — **one row per unique chip** (33 rows, 81
columns, built by `ywf_ocv_channels_csv.py`, which **refuses to write** unless the
mask layer verifies). `echo_mask_present` marks the 17 echo-bearing rows; on the
other 16 the six mask dimensions are empty, because a mask vector without a mask
does not exist:

| group | columns |
| --- | --- |
| provenance | `id`, `asset_*`, `request_id`, `segment_index`, `acquisition_ts`, `pass_label`, `orbit_direction`, `status`, `…` (34 source columns) |
| derived | `segment`, `section_index`, `state`, `pre_collapse`, `day`, `month`, `year`, `season`, `echo_mode`, `echo_mask_present`, `payload_group_size`, `payload_shared_with_segments` |
| recomputed mask | `A`, `D`, `F`, `S`, `P`, `gamma2_mask`, `mask_bbox_area`, `mask_row_span`, `mask_col_span`, `mask_centroid_row/col`, `mask_peak_row/col`, `mask_n_components`, `mask_largest_component`, `mask_fragmentation` |
| pipeline column | `coherence_masked_pixels`, `coherence_gamma2`, `intensity`, `brightness_ratio`, … (kept as provenance, not as the mask) |

`echo_mode` is cut on this site's own scale (`A <= 2` compact, `3..5`
intermediate, `>= 6` distributed of 49 px) — the project-wide cuts belong to the
80 x 80 bridge chips and would put every YWF echo in one bucket. The cut used is
recorded in the metadata (`echo_mode_cuts`).

### `request_id` is not the burst id

`request_id` in the provenance group is the **measurement extract's own** request
column (field 5 of `data/ywf/ywf_measurements_full.txt`, one single distinct value
in this record, and no value of it occurs in the window index). It is *not* the id
of the window's burst: the committed window index carries a **CDSE burst id** as
its 6th field, and this package's earlier label collided the two names. The index
reader (`ywf_ocv_masks.py::read_index`, `ywf_ocv_core.chip_table`) now calls that
field **`burst_id`**, and its value is resolved per chip by
`ywf_bursts_resolve.py` (never taken for the extract's `request_id`, which the CSV
reads out of the extract itself).

### `A` is not `coherence_masked_pixels`

The pipeline's own column was computed on the **pre-fix asset-point chip of the
whole-scene window** (75, 133, … px) while this package recomputes the echo mask
of the committed 7 x 7 window (2..7 px). The two agree on **0 of 17** rows
(`n_coherence_masked_pixels_matched` in the metadata, panel E a), so the CSV keeps
both and reports the coincidence instead of merging them.

## The reference file — defined here, then pinned

The other sites pin against JSONs their upstream analysis project committed.
**This site has no upstream analysis project** (the record *is* the site's own
export), so the reference layer is defined by this package:

* `fig_ywf_ocv_paper.py --write-reference` writes
  `data/ywf/reference/ywf_ocv_findings.json` and proves the round trip (the file is
  re-read and re-compared field by field);
* every later run **without** the flag re-derives the whole file from the
  committed tables and aborts (exit code != 0) on the first deviation.

The reference carries no timestamp and is written with sorted keys, so two runs
are byte-identical; the mask cache, the window payload, the measurements, the
segments, the manifest and the channel CSV are pinned by **sha256** as well. A
`nan` — a test that cannot be computed, here the Welch/Mann-Whitney pair of a
control whose two states are nearly identical (`coherence` = 0.999 throughout) —
is stored as JSON `null`, never as the non-standard `NaN` literal.

## The two resolved references — the pass and the road

Two facts the export does **not** carry are resolved from public catalogues and
committed under `data/ywf/reference/`, each with its own module, its own writer
and its own offline `--check`:

| reference | written by | what it fixes |
| --- | --- | --- |
| `ywf_bursts.json` (49 kB) | `ywf_bursts_resolve.py --write` | the **CDSE burst** of every committed window: the index's 6th field *is* the burst id, and each of the 33 ids resolves to the SLC burst it was cut from — sub-swath, relative orbit, polarisation, azimuth time, parent product, S3 object and footprint |
| `ywf_osm_road.json` (91 kB) | `ywf_osm_anchor.py --fetch` | the committed **map extract**: `OSM way 577440814` (v8, `highway=tertiary`, `surface=asphalt`, 346 nodes), ODbL |
| `ywf_geometry.json` | `ywf_osm_anchor.py --write` | the derived **anchor**: the nearest point of that way (36.4234168 N 129.4206815 E, 36.85 m from the asset point), the carriageway centreline and its bearing, the **estimated grid** — metres per pixel per axis, from the committed burst footprints — and the two readings of it: `road_vs_pixel_grid` (the carriageway as a number of pixels) and `road_vs_window` (the **containment proof**: the road lies *outside* the 7 x 7 frame) |

`ywf_bursts.json` and `ywf_osm_road.json` are committed **catalogue/map extracts**
(each carries its own `source` note); `ywf_geometry.json` is derived from them and
is the file the rest of the package reads. The road anchor is the **carriageway
centreline**, never an area or a polygon: the family's Bergsøysund analysis
measured a −81..+7 m shape error (rms 56 m) for a polygon anchor.

The catalogue is a **mutable external input**, so `ywf_bursts.json` is pinned
*structurally* rather than re-derived: `--check` re-verifies every id against the
committed window index (whose sha256 the file records) and every id's
`OrbitDirection` against the measurement extract's own `orbit_direction` column.
Both references are pinned by **sha256** in the figures' pin block, through the
`sources` block that `fig_ywf_ocv_paper.py` writes into
`reference/ywf_ocv_findings.json`.

The estimate cross-checks in two directions: 3.35 m/px (IW3) and 3.45 m/px (IW2)
**across range** against Bautzen's *committed* rect spec of 3.37 m, and against
the pipeline's own conversion constant `MIN_GROUND_RANGE_PIXEL_M = 3.4` m
(`tower/backend/src/asset_onboarder/insar_monitor.rs`); 14.4-15.3 m/px **across
azimuth** against Bautzen's 12.38 m. The scale is therefore consistent — and it is
still an estimate, because this record carries no rect spec of its own.

## What the record cannot do — the four structural blockers

The echo mask is the **only** layer this record can carry, and the reason is
structural. None of the four facts below is about the state of the tower:

| # | structural fact of the export | what it costs |
| --- | --- | --- |
| 1 | **one 7 x 7 px window** per acquisition — 49 complex samples. The size is the backend's tower constant (`TOWER_WINDOW_SIZE = 7`, `tower/backend/src/asset_onboarder/insar_monitor.rs`): "a tower resolves as a point scatterer in a single IW SLC scene, so a small window around the location is the correct measure" | every number in this package is computed on 49 samples; a band, edge or profile test needs a *resolvable ground extent per pixel* along the band axis |
| 2 | **no rect spec** — a window carries no pixel-spacing fields | the ground extent is only **estimated**, from the committed footprints (`chip_grid`); Bautzen, which has a committed rect spec, is the one site of the family that gets a deck-edge layer |
| 3 | **no per-section coordinates** — the export's `segment_latitude` / `segment_longitude` are empty for all five mast sections | no second anchor and no per-girder re-extraction, so the record is **not segment-resolved** (`segment_resolved = false`, 24 duplicate windows); a per-mast-section layer would be an artefact of that duplication |
| 4 | the target is **outside the frame and sub-pixel** — two independent facts, both committed in `road_vs_window`: (a) **containment** — the road anchor is 36.85 m from the asset point (`road.distance_from_asset_point_m`), so the carriageway's **near edge** lies **8.10 px (ASC/IW3) / 9.33 px (DESC/IW2) across range** from the window centre, while a 7 x 7 frame reaches 6 px at most (3 px from its centre); (b) **resolution** — a carriageway taken as 7.0 m (`road.carriageway_width_m`, itself an *assumption*: the OSM way carries no `width` tag) is **2.03-2.09 px across range but 0.46-0.49 px across azimuth** (`road_vs_pixel_grid`) | **not one committed pixel crosses the carriageway**, so no deck-like edge exists in this record, at any state of the tower — and a frame that *did* contain the road would still have to resolve it |

Blocker 1 is a property of the **export** (an upstream constant), not of this
package: nothing here can widen a committed window, which is why the fourth
column of the anchor file is `"status": "anchor reference — no edge layer
claimed"` instead of a layer. Blocker 4 is committed in two halves there
(`road_vs_window`), and they price the layer from both sides: **containment**
needs 19-21 px per side (27-29 px with Bautzen's ±4-row ROI) — 361-841 samples
against the committed 49, a factor of 7.4 to 17.2, quadratic in the side — while
**resolution** needs a *finer* grid rather than a wider frame (the same ~7 m
cross-section is 2.03-2.09 px across range at 3.35-3.45 m/px and would be 4.7 px
at 1.5 m/px). Neither is a re-reading of the committed record: both are
properties of the acquisition, so the layer stays gated on `TOWER_WINDOW_SIZE`
(the staged step C of this package) and on the product level of the export. The
levers, and who owns each one, are listed in the section below.

### The two no-echo modes

The 16 echo-free chips are echo-free in **two** different ways, and neither is an
empty payload (every payload is non-zero). The rule is `intensity >= 0.30 * peak
AND intensity >= 5.0 * median(intensity)` with `A = count(mask) >= 2`:

| mode | count | how it arises |
| --- | --- | --- |
| `A == 1` | **11** | exactly **one** pixel clears both thresholds, so the mask is *refused* although a bright pixel exists |
| `A == 0` | **5** | **no** pixel clears both: the peak does not stand 5x above the median — a flat window |

```bash
/usr/bin/python3 - <<'PY'    # the two modes over the 16 echo-free chips
import sys, collections; sys.path.insert(0, "code/ywf")
import numpy as np; import ywf_ocv_core as c
mode = collections.Counter()
for r in c.chip_table()[1]:
    if r["mask"] is None:
        I = np.abs(r["z"]) ** 2
        mode[int(((I >= 0.30 * I.max()) & (I >= 5.0 * np.median(I))).sum())] += 1
print(dict(sorted(mode.items())))          # {0: 5, 1: 11}
PY
```

Neither mode is predicted by a structural condition: both occur in **both**
orbits/sub-swaths (ASC-IW3: 4 x `A=0`, 6 x `A=1`; DESC-IW2: 1 x `A=0`, 5 x
`A=1`) and with both polarisations. The split is the *rule boundary* itself — one
further pixel either clears the second threshold or none does — and
`echo_mask_present` collapses both into the same `0`. The echo-free side of the
coverage table is therefore **not** one physical regime, and the CSV, the JSON and
the figure do not pretend that it is.

## What could improve this record — and what cannot

Every limit above is a property of the **export**, so the useful question is not
"which analysis would change the margin" but "which lever exists, and who owns
it". Two of the levers are inside this repository, four are not, and
`ywf_mask_sensitivity.py` measures the first one offline:

| lever | what it would buy | can this repository pull it? |
| --- | --- | --- |
| the rule of the mask layer (`peak_frac`, `median_mult`, `min_n_masked`) | a different census: at `min_n_masked` 2 the unique chips carrying an echo fall **33 / 17 / 5 / 2 of 33** as `median_mult` goes 2.5 / 5 / 7.5 / 10, and the committed 16/27 pre against 1/6 post is reproduced by **5 of the 20** (`peak_frac`, `median_mult`) pairs at the committed `median_mult` 5; `min_n_masked` 3 gives 13 of 33 (12/27 pre). Each knob has a knife edge several steps wide — `peak_frac` 0.3 moves over [0.02, 0.42] and `median_mult` 5 over [5, 5.25] without changing a single count — so the committed census is not balanced on one crossing | **yes — and that is why it is not done.** The rule is *family-wide* (`sarsolve::coherence::coherence_mast_echo`, Carola's `coherence_gamma2`) and every published number of the eight packages hangs on its committed version, so a site-local variant would be a new layer, not a better one. `ywf_mask_sensitivity.py` measures it and never switches it on — the neighbour-referenced variant it also evaluates, for instance, moves the census *without* reproducing the headline (32 of 33 chips at k = 1.5, 23 at k = 3) |
| a **wider window** (`TOWER_WINDOW_SIZE = 7`, backend constant) | containment: the near carriageway edge is 8.10 px (ASC/IW3) / 9.33 px (DESC/IW2) across range from the window centre, so **19-21 px per side** would hold it (27-29 px with Bautzen's ±4-row ROI) — 361-841 samples against the committed 49 (`road_vs_window`) | **no** — the samples were never exported, so no re-reading of this record can produce them; the window width is the backend's step C, frozen at 7 |
| a **finer ground sampling** (say 1.5 m/px instead of 3.35-3.45 m/px across range) | resolution: the same ~7 m cross-section becomes 4.7 px across range instead of 2.03-2.09 — a cross-section a band test could measure instead of a sub-pixel one | **no** — a window carries no rect spec, the committed grid is only *estimated* from the burst footprints, and there is nothing in the record to resample |
| **per-section coordinates** in the export | a second anchor and per-girder re-extraction, i.e. a segment-resolved layer instead of the 24 duplicated windows | **no** — `segment_latitude` / `segment_longitude` are empty for all five mast sections (`segment_resolved = false`) |
| **more post-event acquisitions** | power: with the post side echo-free, 16/27 would already be significant — 4 post-event chips at zero echo suffice at alpha = 0.05 (`echo_coverage.power.post_chips_needed_at_zero_echo_for_fisher_p_lt_alpha` in the reference) — against the 6 committed, one of which carries the surviving echo | **no** — the record is what was acquired; the honest output is the coverage table with its power note, not a contrast |
| the **orbit / sub-swath split** (ASCENDING IW3 / DESCENDING IW2) | a second, independent cross-cut of the same coverage (10 of 20 against 7 of 13 chips, Fisher p = 1.0000, within-orbit pre/post tests p = 0.5820 / 0.1923) | **already in the layer** — `echo_coverage_by_orbit` is pinned like the state block, and it cannot repair the power limit, which is set by the 6 post-event chips |

The ranking is therefore not "what would make the number significant" but "what
is a property of the acquisition": four of the six levers lie outside this
analysis, one is already spent (the orbit split), and the only one the package
*does* control — the rule — is the one it must not touch, because five sibling
packages share it. That asymmetry is why the deliverable here is a *limit*
(`status`: "anchor reference — no edge layer is claimed") plus the coverage
table, and not a better-looking verdict.

## What the figure shows (and what it does not)

| panel | content |
| --- | --- |
| `A` | the six echo-mask channels per state, over the **17 echo-bearing chips** — a box where a state has ≥ 3 of them, a single marker for the one post-event chip |
| `B` | the **headline**: (a) echo coverage per state with its Fisher exact test; (b) the coverage timeline across the collapse (marker area = unique chips of the day) |
| `C` | in-sample effect sizes: (a) Cliff's delta + bootstrap CI95 and SDS of the six mask channels; (b) the same for the co-variate controls |
| `D` | the OCV plane (`gamma2` vs `P`, `gamma2` vs `D`, `P` vs `D`), 16 pre chips and 1 post chip |
| `E` | why the contrast is not a damage signal: (a) the export `A` vs. the recomputed `A`; (b) the echo-mode mix; (c) the season composition; (d) the 24 duplicated windows |

**No classifier is fitted.** Panel D is a scatter, not a decision boundary: one
post-event echo chip is not a class, so a leave-one-out accuracy computed from it
would be a property of the sample size, not of the tower. The report says so
instead of printing a number.

**The site's claim is an echo-coverage statement.** 16 of 27 pre-event chips
carry an echo, 1 of 6 post-event chips does; the Fisher exact test gives
p = 0.085 — not significant, and it cannot be, because that single surviving echo
is in the post side. Had all six post-event chips been echo-free, four chips
would already have sufficed at alpha = 0.05
(`echo_coverage.power.post_chips_needed_at_zero_echo_for_fisher_p_lt_alpha` in
the reference). What the record supports is the *observability* statement: the
echo mask is the quantity that can be observed on this collapsed tower, and its
coverage is the number to extend before any state contrast is attempted — never a
damage verdict.

**The orbit split does not change that.** Both orbits carry the same coverage
(10 of 20 ascending chips, 7 of 13 descending; Fisher exact p = 1.0000), and the
one surviving post-event echo is an ascending pass, so the descending side holds
none at all. The block is pinned like the state block
(`echo_coverage_by_orbit`, with the within-orbit pre/post Fisher tests
p = 0.5820 / 0.1923), and it exists so that the record's *second* cross-cut is
stated rather than assumed — it cannot repair the power limit above, which is set
by the 6 post-event chips.

## The modules

| file | role |
| --- | --- |
| `ywf_ocv_stats.py` | the statistics of the family (Cliff's delta + bootstrap CI, SDS, Welch/Mann-Whitney, sign test, Spearman/Mann-Whitney, shrinkage LDA + LOO + label permutation, block-permutation correlation, OLS) — ported verbatim, unused parts included so the eight packages read as one diff (the ported `pair_block(..., orbit=...)` parameter stays unused here: this package reports the orbit *coverage* of the mask, and the 17 echo chips carry no contrast to correlate) |
| `ywf_ocv_masks.py` | the mask rule, the payload reader, the cache writer (`--extract`), the acceptance test (`--check`), the 6D vector of a chip |
| `ywf_ocv_core.py` | paths, the site/event/state/orbit constants, the state, echo-mode and orbit definitions, the loaders of the committed tables, `chip_table`, `coverage_block`, `coverage_by_orbit`, `fisher_p`, and a read-only CLI brief |
| `ywf_ocv_channels_csv.py` | the channel-table generator (build + verify) |
| `ywf_bursts_resolve.py` | the pass geometry of the record: resolves the 33 committed window ids against the public **CDSE** burst catalogue (`--write`, network) and re-verifies the committed file **offline** against the committed index and the measurement extract (`--check`, the default) |
| `ywf_osm_anchor.py` | the road anchor: fetches the **OSM** way (`--fetch`), re-reads the committed map extract, derives the nearest point, the centreline and the estimated burst grid, commits the **containment proof** of the edge layer's absence (`road_vs_window`) and states the limit of the record (`--check`, the default) |
| `ywf_mask_sensitivity.py` | the **offline** sensitivity study of the mask rule: the coarse (`peak_frac` x `median_mult` x `min_n_masked`) grid, the knife-edge scan of every knob, the one neighbour-referenced variant of the rule and a `--json` dump — it re-derives the committed coverage first, writes no cache/CSV/reference, and never switches a rule version on |
| `fig_ywf_ocv_paper.py` | panels A–E, the result JSON (`data/ywf/fig_ywf_ocv_paper.json`), the report (`figures/ywf/fig_ywf_ocv_paper.md`), the reference writer and the pin block |
| `fig_ywf_ocv_weeks.py` | the pre-window companion figure: the same channels over the **27 pre-event chips / 23 ISO weeks**, pinned against this package's result JSON and reference (`figures/ywf/fig_ywf_ocv_weeks.md`) |

## Pins

Every run prints its pin summary and exits non-zero on any deviation (1179
checks on the current record; 920 when `--write-reference` defines the reference
itself, 422 in `fig_ywf_ocv_weeks.py`):

| pinned against | what |
| --- | --- |
| `data/ywf/reference/ywf_ocv_findings.json` | the echo-coverage table and its Fisher test, the same coverage split by orbit direction (with the within-orbit pre/post tests), the mask dimensions per state, the pipeline column, the co-variate controls, the month/section composition, the 32-row timeline and the payload-digest groups — field by field |
| `data/ywf/ywf_ocv_channels_meta.json` | the table's own definition, guard counts, per-state summaries **and the per-orbit coverage block** against this package's recomputation (the generator itself refuses to write unless the per-orbit cross-check passes) |
| `data/ywf/ywf_windows_mask_cache.txt` | each cached chip's `n_masked` / `n_components` against the channel table's `A` / `F` (17 rows) |
| the committed inputs | their sha256 digests (payload, cache, manifest, measurements, segments) |
| the two resolved references | the sha256 digests of `ywf_bursts.json`, `ywf_osm_road.json` and `ywf_geometry.json`, recorded in the `sources` block of `reference/ywf_ocv_findings.json` and re-derived from the files on every run — by `fig_ywf_ocv_paper.py` **and** by `fig_ywf_ocv_weeks.py` |

Two runs of `fig_ywf_ocv_paper.py` with the fixed seed are byte-identical — the
result JSON, the report and all five PNGs (checked with `sha256sum`, including a
run into a scratch directory); no wall-clock value is written and the runtime
stays on stdout. `--write-reference` is idempotent: re-writing the reference
leaves the committed `data/ywf/reference/ywf_ocv_findings.json` byte-for-byte
unchanged.

## Reading

Before the collapse **16 of 27** unique chips carry an echo, after it **1 of 6**.
With six post-event acquisitions and one surviving echo the site cannot separate
the states, and the single post-event echo chip lies inside the pre-event cloud of
`gamma2` / `P` / `D` (panel D) — so the effect sizes of panel C are reported to
*document* the co-variation of season, weather and acquisition generation with the
state, not as a mechanism. The honest output of this package is the coverage table
of panel B and the power statement that goes with it.

