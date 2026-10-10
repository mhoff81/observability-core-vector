# Yeongdeok Wind Farm (Samgye-ri) — observability core vector, pre/post collapse

**Site** Yeongdeok Wind Farm (Samgye-ri) · **asset** Yeongdeok Wind Turbine (Samgye-ri) · **event** 2026-02-02 · **rows** 33 unique 7x7 chips (57 committed windows on 32 dates) · **states** `pre` / `post`

```
OCV  = [gamma2, P, D]              observability core vector
x    = [gamma2, P, D, A, F, S]      echo-mask vector
mask = intensity >= 0.30*peak AND >= 5.0*np.median(intensity), A >= 2
layer: echo mask (whole 7x7 window) — the record is not segment-resolved, so there is no per-mast-section layer
cross-cuts: state (pre/post) · orbit (ASC/DESC, descriptive — both orbits straddle the event)
```

## What is recomputed, and what is pinned

| input | role |
| --- | --- |
| `data/ywf/ywf_windows_full.txt` | the committed window payloads (38 kB, committed in full) |
| `data/ywf/ywf_windows_mask_cache.txt` | the echo mask of every committed chip — the mask layer's committed input |
| `data/ywf/ywf_segments.txt` | the mast-section labels + FEM baselines |
| `data/ywf/ywf_ocv_channels.csv` | generated: one row per unique chip, 81 columns |
| `data/ywf/reference/ywf_ocv_findings.json` | the committed analysis-layer reference this script defines and pins |
| `data/ywf/reference/ywf_bursts.json` | resolved: which CDSE burst, sub-swath, relative orbit and polarisation each of the 33 chips is (`ywf_bursts_resolve.py --check`) |
| `data/ywf/reference/ywf_osm_road.json` + `ywf_geometry.json` | resolved: the OSM way the tower fell on, and the burst grid estimated from the committed footprints (`ywf_osm_anchor.py --check`) |

The record holds **57** committed windows on **32** acquisition dates, **33** of them unique chips after dropping **24** byte-identical duplicates; `segment_resolved = False`. All 33 unique chips are the same mast section (`turbine-mast-section-3`): the 24 duplicated windows differ only in the requested section index (2 vs. 3) and carry the identical payload, so **one** mask layer is all the record supports (fig. E d).

`A` is *not* the pipeline's own `coherence_masked_pixels`: that column was computed on the asset-point chip of the whole-scene window while this package recomputes the echo mask of the committed 7x7 window — the two are equal on **0** of **17** rows compared (fig. E a).

## A — the headline: echo coverage

The 6D vector exists only where the mask rule holds, and the collapse contrast this site can carry is the *coverage* of that echo:

| state | unique chips | with echo | without echo | echo rate |
| --- | --- | --- | --- | --- |
| pre | 27 | 16 | 11 | 59.3 % |
| post | 6 | 1 | 5 | 16.7 % |

Fisher exact on the 2x2 table `['pre', 'post'] x ['echo', 'no_echo']` = [[16, 11], [1, 5]]: **p = 0.0854** — not significant, and it cannot be: the post side holds 6 chips of which 1 still carries an echo. Had all 6 been echo-free, the same test would need only **4** post-event chips at alpha = 0.05 — so this record's coverage contrast is limited by the *surviving echo*, not only by the length of the record.

### A2 — the same coverage, split by orbit direction

Every unique chip of this record is the same mast section, so beside the state the only cross-cut the record carries is the pass geometry (ascending / descending). Both orbits pass over the site *and* both straddle the collapse date:

| orbit | unique chips | with echo | without echo | echo rate | pre | post |
| --- | --- | --- | --- | --- | --- | --- |
| ASC | 20 | 10 | 10 | 50.0 % | 9/16 | 1/4 |
| DESC | 13 | 7 | 6 | 53.8 % | 7/11 | 0/2 |

The two orbits are **indistinguishable** in this record: 10 of 20 ascending chips carry an echo (0.500) against 7 of 13 descending ones (0.538). Fisher exact on `['ASCENDING', 'DESCENDING'] x ['echo', 'no_echo']` = [[10, 10], [7, 6]]: **p = 1.0000**. `pass_label` (morning / afternoon) is the identical cut — ASC = morning, DESC = afternoon.

The one surviving post-event echo (`2026-03-16`) is an **ascending** pass, so the descending side holds no post-event echo at all (0 of 2 chips). The block is therefore read as a *description* of where the 17 echoes sit over the two pass geometries, never as an orbit contrast: with 6 post-event chips it cannot repair the power limit of §A, and the pre/post tests within an orbit stay descriptive too (`ASC` p = 0.5820, `DESC` p = 0.1923).

## B — the channels (where the mask holds)

| channel | n pre/post | pre median | post median | delta | SDS | CI95 excludes 0 |
| --- | --- | --- | --- | --- | --- | --- |
| gamma2 | 16/1 | 0.4345 | 0.1886 | -0.375 | 3.750 | no CI |
| P | 16/1 | 0.1348 | 0.1647 | 0.250 | 2.500 | no CI |
| D | 16/1 | 0.2778 | 0.1190 | -0.750 | 7.500 | no CI |
| A | 16/1 | 3.0000 | 5.0000 | 0.812 | 8.125 | no CI |
| F | 16/1 | 2.0000 | 3.0000 | 0.938 | 9.375 | no CI |
| S | 16/1 | 1.5926 | 1.4142 | -0.125 | 1.250 | no CI |

The post column is **one** chip (`1` of the 6 post-event chips carry an echo), so the bootstrap CI is undefined there (`st.delta_block` declines to resample a single value) and `SDS` documents only the spread of the pre side. Fig. A draws those distributions, fig. C their effect sizes beside the co-variate controls.

## C — the OCV plane

Fig. D plots `gamma2` against `P` and `D` for the 17 echo-bearing chips. No classifier is fitted and no leave-one-out accuracy is reported: one post-event echo chip is not a class, and an accuracy computed from it would be an artefact of the sample size rather than a property of the site.

## D — what else moves with the collapse date?

The pre side is Aug-Jan and the post side Feb-Apr, so season, weather and acquisition generation move with the state. The co-variate controls beside the mask channels (fig. C b, fig. E):

| control | n pre/post | pre median | post median | delta | CI95 |
| --- | --- | --- | --- | --- | --- |
| wind_speed_ms | 27/6 | 3.130 | 4.860 | 0.167 | [-0.407, 0.704] |
| temperature_c | 27/6 | 17.200 | 11.100 | -0.302 | [-0.679, 0.111] |
| brightness_ratio | 27/6 | 0.317 | 0.818 | 0.259 | [-0.259, 0.704] |
| intensity | 27/6 | 2551.960 | 6585.013 | 0.259 | [-0.259, 0.704] |
| coherence | 27/6 | 0.999 | 0.999 | 0.000 | [0.000, 0.000] |
| coherence_gamma2 | 27/6 | 0.010 | 0.034 | 0.580 | [0.222, 0.864] |
| cumulative_mm | 27/6 | 0.270 | 14.167 | 0.222 | [-0.309, 0.704] |
| structural_frequency_hz | 6/2 | 2.493 | 2.738 | 0.167 | [-0.667, 1.000] |

Echo-mode mix (on this site's own scale, `A <= 2`, `3 <= A <= 5`, `A >= 6` of the 49-px window): pre {'compact': 4, 'intermediate': 11, 'distributed': 1} over 16 echo chips, post {'compact': 0, 'intermediate': 1, 'distributed': 0} over 1.

## What this record cannot do

The echo mask is the **only** layer this record can carry. Four structural facts of the export set that limit, and none of them is about the state of the tower:

| structural fact | consequence |
| --- | --- |
| one **7 x 7 px window** per acquisition (49 px, `TOWER_WINDOW_SIZE = 7`) | every number above is computed on 49 complex samples; a band or an edge test needs a resolvable ground extent per pixel *along* the band axis |
| **no rect spec** — a window carries no pixel spacing | the ground extent is only *estimated*, from the committed burst footprints (`data/ywf/reference/ywf_geometry.json`); Bautzen, which carries a rect spec, is the one site of the family that gets a deck-edge layer |
| **no per-section coordinates** — the export's `segment_latitude` / `segment_longitude` are empty for all five mast sections | there is no second anchor and no per-girder re-extraction, so the record is not segment-resolved (`segment_resolved = False`, 24 duplicate windows) |
| the target is **outside the frame and sub-pixel** — the road anchor stands ~36.85 m from the asset point, so the carriageway's near edge is 8.10 px (ascending/IW3) / 9.33 px (descending/IW2) across range from the window centre against the 6 px a 7 x 7 frame reaches (`road_vs_window`), and a ~7 m carriageway is a couple of pixels across the range axis but under one across the azimuth axis (`road_vs_pixel_grid`) | not one committed pixel crosses the road, so no deck-like edge exists in this record, at any state of the tower |

The 16 echo-free chips are echo-free in **two** different ways, and neither is an empty payload: some have exactly one pixel passing both thresholds (the rule requires `A >= 2`, so the mask is refused), the others have *none*, because the peak does not stand 5x above the median — a single bright pixel and a flat window both collapse into the same `echo_mask_present = 0` here.

The record's second cross-cut is a **sub-swath** split, not just a pass split: ascending is IW3 on relative orbit 54 (09:2x UT), descending is IW2 on relative orbit 61 (21:2x UT). At 129.4 E the ascending pass is the local **evening** one (18:00 local solar time) and the descending one the local morning (06:00), so the export's own `pass_label` (ascending = `morning`) follows the UTC clock and is not a statement about local time or illumination. The record also interleaves **VV and VH** acquisitions; the mask rule is a per-window relative threshold, so the mask channels are invariant to that gain, but the export's absolute columns (`intensity`, `brightness_ratio`) are not — they stay provenance (`data/ywf/reference/ywf_bursts.json` resolves the polarisation of every chip).

## Pins

Every number above was re-derived from the committed tables and compared field by field against `data/ywf/reference/ywf_ocv_findings.json`.

1179 checks (1179 exact, 0 within 1e-12), **0 deviations**.

| committed layer | what is pinned |
| --- | --- |
| `data/ywf/reference/ywf_ocv_findings.json` | the echo-coverage table, its Fisher test, the mask dimensions per state, the pipeline column, the co-variate controls, the timeline and the record's provenance — field by field |
| `data/ywf/ywf_ocv_channels_meta.json` | the table's own definition and summary blocks against this script's recomputation |
| `data/ywf/ywf_windows_mask_cache.txt` | each cached chip's `n_masked` / `n_components` against the CSV's `A` / `F` (17 rows) |
| the committed inputs | their sha256 digests (`ywf_windows_full.txt`, cache, manifest, measurements, segments) |
| the two resolved references | the sha256 digests of `ywf_bursts.json` (the pass geometry of the chips), `ywf_osm_road.json` and `ywf_geometry.json` (the road and the estimated grid), each re-derived and checked by its own module |

### Reading (not a claim of damage)

Before the collapse **16 of 27** unique chips carry an echo; after it **1 of 6**. With 6 post-event acquisitions the site cannot reach significance, and the single post-event echo chip lies inside the pre-event cloud of `gamma2` / `P` / `D` (fig. D). What the record supports is the *observability* statement: the echo mask is the quantity that can be observed on this collapsed tower, and its coverage is the number to extend before any state contrast is attempted.

```bash
python3 code/ywf/ywf_ocv_masks.py --check                  # the mask rule
python3 code/ywf/ywf_ocv_channels_csv.py                 # build + verify
python3 code/ywf/fig_ywf_ocv_paper.py                    # figures + pins
python3 code/ywf/fig_ywf_ocv_paper.py --write-reference  # re-define
```
