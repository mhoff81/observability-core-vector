# Haltern AV9 Wind Park (Kreis Recklinghausen, NRW) — observability core vector, pre/post collapse

**Site** Haltern AV9 Wind Park (Kreis Recklinghausen, NRW) · **asset** Haltern AV9 Wind Turbine (Nordex N149/4.5) · **event** 2021-09-29 · **rows** 120 unique 7x7 chips (168 committed windows on 120 dates) · **states** `pre` / `post`

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
| `data/haltern/haltern_windows_full.txt` | the committed window payloads (38 kB, committed in full) |
| `data/haltern/haltern_windows_mask_cache.txt` | the echo mask of every committed chip — the mask layer's committed input |
| `data/haltern/haltern_segments.txt` | the mast-section labels + FEM baselines |
| `data/haltern/haltern_ocv_channels.csv` | generated: one row per unique chip, 81 columns |
| `data/haltern/reference/haltern_ocv_findings.json` | the committed analysis-layer reference this script defines and pins |
| `data/haltern/reference/haltern_bursts.json` | resolved: which CDSE burst, sub-swath, relative orbit and polarisation each of the 120 chips is (`haltern_bursts_resolve.py --check`) |
| `data/haltern/reference/haltern_osm_node.json` + `haltern_geometry.json` | resolved: the turbine's **own OSM node** (not a carriageway - this tower fell inside the wind farm), the anchor derived from it and the burst grid estimated from the committed footprints (`haltern_osm_anchor.py --check`) |

The record holds **168** committed windows on **120** acquisition dates, **120** of them unique chips after dropping **48** byte-identical duplicates across **48** payload groups; `segment_resolved = False`. All 120 unique chips are the same mast section (`turbine-mast-section-0`, `turbine-mast-section-1`): the duplicated windows differ only in the requested section index (0 vs. 1) and carry the identical payload, so **one** mask layer is all the record supports (fig. E d).

`A` is *not* the pipeline's own `coherence_masked_pixels`: that column was computed on the asset-point chip of the whole-scene window while this package recomputes the echo mask of the committed 7x7 window — the two are equal on **6** of **66** rows compared (fig. E a).

## A — the headline: echo coverage

The 6D vector exists only where the mask rule holds, and the collapse contrast this site can carry is the *coverage* of that echo:

| state | unique chips | with echo | without echo | echo rate |
| --- | --- | --- | --- | --- |
| pre | 108 | 62 | 46 | 57.4 % |
| post | 12 | 4 | 8 | 33.3 % |

Fisher exact on the 2x2 table `['pre', 'post'] x ['echo', 'no_echo']` = [[62, 46], [4, 8]]: **p = 0.1343** — not significant, and it cannot be: the post side holds 12 chips of which 4 still carries an echo. Had all 12 been echo-free, the same test would need only **4** post-event chips at alpha = 0.05 — so this record's coverage contrast is limited by the *surviving echo*, not only by the length of the record.

### A2 — the same coverage, split by orbit direction

Every unique chip of this record is the same mast section, so beside the state the only cross-cut the record carries is the pass geometry (ascending / descending). Both orbits pass over the site *and* both straddle the collapse date:

| orbit | unique chips | with echo | without echo | echo rate | pre | post |
| --- | --- | --- | --- | --- | --- | --- |
| ASC | 62 | 33 | 29 | 53.2 % | 31/56 | 2/6 |
| DESC | 58 | 33 | 25 | 56.9 % | 31/52 | 2/6 |

The two orbits are **indistinguishable** in this record: 33 of 62 ascending chips carry an echo (0.532) against 33 of 58 descending ones (0.569). Fisher exact on `['ASCENDING', 'DESCENDING'] x ['echo', 'no_echo']` = [[33, 29], [33, 25]]: **p = 0.7166**. `pass_label` (morning / afternoon) is the identical cut — ASC = morning, DESC = afternoon.

The 4 surviving post-event echo chips (`2021-10-02`, `2021-10-08`, `2021-10-10`, `2021-10-17`) carry ASC/DESC pass geometry, and the  side holds 0 post-event echo chips. The block is therefore read as a *description* of where the 66 echoes sit over the two pass geometries, never as an orbit contrast: with 12 post-event chips it cannot repair the power limit of §A, and the pre/post tests within an orbit stay descriptive too (`ASC` p = 0.4054, `DESC` p = 0.3867).

## B — the channels (where the mask holds)

| channel | n pre/post | pre median | post median | delta | SDS | CI95 excludes 0 |
| --- | --- | --- | --- | --- | --- | --- |
| gamma2 | 62/4 | 0.5032 | 0.2795 | -0.315 | 3.145 | no |
| P | 62/4 | 0.0901 | 0.0997 | 0.367 | 3.669 | no |
| D | 62/4 | 0.4750 | 0.3929 | -0.097 | 0.968 | no |
| A | 62/4 | 3.0000 | 4.0000 | 0.089 | 0.887 | no |
| F | 62/4 | 2.0000 | 2.5000 | 0.355 | 3.548 | no |
| S | 62/4 | 0.8667 | 1.6817 | 0.363 | 3.629 | no |

The post column is **one** chip (`4` of the 12 post-event chips carry an echo), so the bootstrap CI is undefined there (`st.delta_block` declines to resample a single value) and `SDS` documents only the spread of the pre side. Fig. A draws those distributions, fig. C their effect sizes beside the co-variate controls.

## C — the OCV plane

Fig. D plots `gamma2` against `P` and `D` for the 17 echo-bearing chips. No classifier is fitted and no leave-one-out accuracy is reported: one post-event echo chip is not a class, and an accuracy computed from it would be an artefact of the sample size rather than a property of the site.

## D — what else moves with the collapse date?

The pre side is Aug-Jan and the post side Feb-Apr, so season, weather and acquisition generation move with the state. The co-variate controls beside the mask channels (fig. C b, fig. E):

| control | n pre/post | pre median | post median | delta | CI95 |
| --- | --- | --- | --- | --- | --- |
| wind_speed_ms | 108/12 | 3.065 | 2.975 | -0.191 | [-0.526, 0.156] |
| temperature_c | 108/12 | 11.450 | 11.400 | 0.008 | [-0.201, 0.209] |
| brightness_ratio | 60/12 | 1.000 | 1.000 | -0.163 | [-0.504, 0.164] |
| intensity | 108/12 | 4586.387 | 3596.364 | -0.230 | [-0.590, 0.150] |
| coherence | 60/12 | 0.999 | 0.999 | 0.000 | [0.000, 0.000] |
| coherence_gamma2 | 108/12 | 0.026 | 0.258 | 0.369 | [0.002, 0.705] |

Carried in the CSV but not contrasted here, because one side holds no value: `cumulative_mm` (0/0), `structural_frequency_hz` (0/0).

Echo-mode mix (on this site's own scale, `A <= 2`, `3 <= A <= 5`, `A >= 6` of the 49-px window): pre {'compact': 22, 'intermediate': 31, 'distributed': 9} over 62 echo chips, post {'compact': 2, 'intermediate': 0, 'distributed': 2} over 4.

## What this record cannot do

The echo mask is the **only** layer this record can carry. Four structural facts of the export set that limit, and none of them is about the state of the tower:

| structural fact | consequence |
| --- | --- |
| one **7 x 7 px window** per acquisition (49 px, `TOWER_WINDOW_SIZE = 7`) | every number above is computed on 49 complex samples; a band or an edge test needs a resolvable ground extent per pixel *along* the band axis |
| **no rect spec** — a window carries no pixel spacing | the ground extent is only *estimated*, from the committed burst footprints (`data/haltern/reference/haltern_geometry.json`); Bautzen, which carries a rect spec, is the one site of the family that gets a deck-edge layer |
| **no per-section coordinates** — the export's `segment_latitude` / `segment_longitude` are empty for all five mast sections | there is no second anchor and no per-girder re-extraction, so the record is not segment-resolved (`segment_resolved = False`, 48 duplicate windows) |
| the target is **outside the frame** — the tower base alone is 2.4-2.6 px across range and 0.6-0.6 px across azimuth, against the 6 px a 7 x 7 frame reaches (3 px from its centre), so **not one** committed pixel spans the base on any orbit (`structure_vs_window`); the N149 rotor (149 m = 39.6-43.4 px across range) *would* span a wider frame, but a rotor is not a band/edge object — its ground projection is a line whose length and direction are set by the yaw, which neither the export nor the catalogue carries (`structure_vs_pixel_grid`) | no deck-like edge exists in this record, at any state of the tower |

The 54 echo-free chips are echo-free in **two** different ways, and neither is an empty payload: some have exactly one pixel passing both thresholds (the rule requires `A >= 2`, so the mask is refused), the others have *none*, because the peak does not stand 5x above the median — a single bright pixel and a flat window both collapse into the same `echo_mask_present = 0` here.

The record's second cross-cut is a **sub-swath** split, not just a pass split, and the sub-swath is not fixed per direction: `swath`=IW1/IW3 `rel. orbit`=15/88 `UTC`=15:16/15:25/17:16/17:17/17:24/17:25 `pol`=VH/VV over 62 days; `swath`=IW1/IW3 `rel. orbit`=37/139 `UTC`=05:41/05:42/05:49/05:50 `pol`=VH/VV over 58 days (the daily pick is the strongest window among the offered bursts, so it may alternate sub-swaths). Haltern sits at 7.09 E, where local solar time runs 28 min ahead of UT, and the export's own `pass_label` (afternoon, morning) follows the *platform's* naming convention, not local time or illumination. The record also interleaves VH and VV acquisitions; the mask rule is a per-window relative threshold, so the mask channels are invariant to that gain, but the export's absolute columns (`intensity`, `brightness_ratio`) are not — they stay provenance (`data/haltern/reference/haltern_bursts.json` resolves the polarisation of every chip).

## Pins

Every number above was re-derived from the committed tables and compared field by field against `data/haltern/reference/haltern_ocv_findings.json`.

2117 checks (2115 exact, 2 within 1e-12), **2 deviations**.

| committed layer | what is pinned |
| --- | --- |
| `data/haltern/reference/haltern_ocv_findings.json` | the echo-coverage table, its Fisher test, the mask dimensions per state, the pipeline column, the co-variate controls, the timeline and the record's provenance — field by field |
| `data/haltern/haltern_ocv_channels_meta.json` | the table's own definition and summary blocks against this script's recomputation |
| `data/haltern/haltern_windows_mask_cache.txt` | each cached chip's `n_masked` / `n_components` against the CSV's `A` / `F` (66 rows) |
| the committed inputs | their sha256 digests (`haltern_windows_full.txt`, cache, manifest, measurements, segments) |
| the two resolved references | the sha256 digests of `haltern_bursts.json` (the pass geometry of the chips), `haltern_osm_node.json` and `haltern_geometry.json` (the turbine's own OSM node, the anchor derived from it and the estimated grid), each re-derived and checked by its own module |

### Reading (not a claim of damage)

Before the collapse **62 of 108** unique chips carry an echo; after it **4 of 12**. With 12 post-event acquisitions the site cannot reach significance, and the single post-event echo chip lies inside the pre-event cloud of `gamma2` / `P` / `D` (fig. D). What the record supports is the *observability* statement: the echo mask is the quantity that can be observed on this collapsed tower, and its coverage is the number to extend before any state contrast is attempted.

```bash
python3 code/haltern/haltern_ocv_masks.py --check                  # the mask rule
python3 code/haltern/haltern_ocv_channels_csv.py                 # build + verify
python3 code/haltern/fig_haltern_ocv_paper.py                    # figures + pins
python3 code/haltern/fig_haltern_ocv_paper.py --write-reference  # re-define
```
