# Yeongdeok Wind Farm (Samgye-ri) — observability core vector, pre/post collapse

**Site** Yeongdeok Wind Farm (Samgye-ri) · **asset** Yeongdeok Wind Turbine (Samgye-ri) · **event** 2026-02-02 · **rows** 33 unique 7x7 chips (57 committed windows on 32 dates) · **states** `pre` / `post`

```
OCV  = [gamma2, P, D]              observability core vector
x    = [gamma2, P, D, A, F, S]      echo-mask vector
mask = intensity >= 0.30*peak AND >= 5.0*np.median(intensity), A >= 2
layer: echo mask (whole 7x7 window) — the record is not segment-resolved, so there is no per-mast-section layer
```

## What is recomputed, and what is pinned

| input | role |
| --- | --- |
| `data/ywf/ywf_windows_full.txt` | the committed window payloads (38 kB, committed in full) |
| `data/ywf/ywf_windows_mask_cache.txt` | the echo mask of every committed chip — the mask layer's committed input |
| `data/ywf/ywf_segments.txt` | the mast-section labels + FEM baselines |
| `data/ywf/ywf_ocv_channels.csv` | generated: one row per unique chip, 81 columns |
| `data/ywf/reference/ywf_ocv_findings.json` | the committed analysis-layer reference this script defines and pins |

The record holds **57** committed windows on **32** acquisition dates, **33** of them unique chips after dropping **24** byte-identical duplicates; `segment_resolved = False`. All 33 unique chips are the same mast section (`turbine-mast-section-3`): the 24 duplicated windows differ only in the requested section index (2 vs. 3) and carry the identical payload, so **one** mask layer is all the record supports (fig. E d).

`A` is *not* the pipeline's own `coherence_masked_pixels`: that column was computed on the asset-point chip of the whole-scene window while this package recomputes the echo mask of the committed 7x7 window — the two are equal on **0** of **17** rows compared (fig. E a).

## A — the headline: echo coverage

The 6D vector exists only where the mask rule holds, and the collapse contrast this site can carry is the *coverage* of that echo:

| state | unique chips | with echo | without echo | echo rate |
| --- | --- | --- | --- | --- |
| pre | 27 | 16 | 11 | 59.3 % |
| post | 6 | 1 | 5 | 16.7 % |

Fisher exact on the 2x2 table `['pre', 'post'] x ['echo', 'no_echo']` = [[16, 11], [1, 5]]: **p = 0.0854** — not significant, and it cannot be: the post side holds 6 chips of which 1 still carries an echo. Had all 6 been echo-free, the same test would need only **4** post-event chips at alpha = 0.05 — so this record's coverage contrast is limited by the *surviving echo*, not only by the length of the record.

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

## Pins

Every number above was re-derived from the committed tables and compared field by field against `data/ywf/reference/ywf_ocv_findings.json`.

1080 checks (1080 exact, 0 within 1e-12), **0 deviations**.

| committed layer | what is pinned |
| --- | --- |
| `data/ywf/reference/ywf_ocv_findings.json` | the echo-coverage table, its Fisher test, the mask dimensions per state, the pipeline column, the co-variate controls, the timeline and the record's provenance — field by field |
| `data/ywf/ywf_ocv_channels_meta.json` | the table's own definition and summary blocks against this script's recomputation |
| `data/ywf/ywf_windows_mask_cache.txt` | each cached chip's `n_masked` / `n_components` against the CSV's `A` / `F` (17 rows) |
| the committed inputs | their sha256 digests (`ywf_windows_full.txt`, cache, manifest, measurements, segments) |

### Reading (not a claim of damage)

Before the collapse **16 of 27** unique chips carry an echo; after it **1 of 6**. With 6 post-event acquisitions the site cannot reach significance, and the single post-event echo chip lies inside the pre-event cloud of `gamma2` / `P` / `D` (fig. D). What the record supports is the *observability* statement: the echo mask is the quantity that can be observed on this collapsed tower, and its coverage is the number to extend before any state contrast is attempted.

```bash
python3 code/ywf/ywf_ocv_masks.py --check                  # the mask rule
python3 code/ywf/ywf_ocv_channels_csv.py                 # build + verify
python3 code/ywf/fig_ywf_ocv_paper.py                    # figures + pins
python3 code/ywf/fig_ywf_ocv_paper.py --write-reference  # re-define
```
