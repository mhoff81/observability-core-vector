# Champlain Towers South, Surfside (FL) — the echo-mask dimensions from the first images to the month after the collapse

**Site** Champlain Towers South, Surfside (FL) · **event** 2021-06-24 · **chips** 48 (windows F0 + M-3..M+1 only)

```
x     = [gamma2, P, D, A, F, S]      echo-mask vector
mask  = intensity >= 0.30*peak AND >= 5.0*np.median(intensity), A >= 2
F0    [2015-09-21, 2015-11-01)       pre-collapse (first images available)
M-3   [2021-03-24, 2021-04-24)       pre-collapse (healthy)
M-2   [2021-04-24, 2021-05-24)       pre-collapse (healthy)
M-1   [2021-05-24, 2021-06-24)       pre-collapse (healthy)
M+1   [2021-06-24, 2021-07-24)       post-collapse (month after)
```

## Why the baseline is 2015 and not 2014

Sentinel-1 data for **this** footprint do not exist in 2014. Sentinel-1A was launched on 2014-04-03, and the suitable orbit (ASC rel-48 IW3) has its **first acquisition on 2015-09-21** (`data/cts/cts_acquisition_census.json`; the Step-0 census already queried from 2014-10-01 and found nothing usable earlier — the descending stratum carries only 6 acquisitions in six years). `F0` is therefore the record-start epoch that carries the first images that actually exist, and is read here as the *first-images* baseline. No 2014 window is invented.

## The five windows

`F0` opens the record where the first Sentinel-1 images of this footprint exist; the collapse date opens `M+1`, so the pre side is the three `M-*` windows (plus the `F0` baseline) and the post side is `M+1` alone. This is a companion of `fig_cts_ocv_paper.py` (figures A–F): it adds one figure and does not touch the A–F pin contract.

| window | range | state | chips | first day | last day | days |
| --- | --- | --- | --- | --- | --- | --- |
| F0 | [2015-09-21, 2015-11-01) | pre | 8 | 2015-09-21 | 2015-10-15 | 2 |
| M-3 | [2021-03-24, 2021-04-24) | pre | 12 | 2021-03-29 | 2021-04-22 | 3 |
| M-2 | [2021-04-24, 2021-05-24) | pre | 8 | 2021-05-04 | 2021-05-16 | 2 |
| M-1 | [2021-05-24, 2021-06-24) | pre | 12 | 2021-05-28 | 2021-06-21 | 3 |
| M+1 | [2021-06-24, 2021-07-24) | post | 8 | 2021-07-03 | 2021-07-15 | 2 |

## Per-window medians of the six dimensions

| dim | F0 (first images) | F0 | M-3 | M-2 | M-1 | M+1 |
| --- | --- | --- | --- | --- | --- | --- |
| gamma2 | 0.02987 | 0.02987 | 0.0281 | 0.0486 | 0.06106 | 0.03362 |
| P | 0.01578 | 0.01578 | 0.02464 | 0.0187 | 0.01963 | 0.01979 |
| D | 0.01835 | 0.01835 | 0.01136 | 0.009964 | 0.01342 | 0.01264 |
| A | 49 | 49 | 31 | 39.5 | 29 | 28.5 |
| F | 21 | 21 | 16 | 16.5 | 14.5 | 16 |
| S | 21.82 | 21.82 | 17.07 | 17.84 | 17.99 | 17.26 |

`F0 (first images)` = the record-start baseline (the 1.0 line of panel f).

## The long baseline — M+1 vs F0 (the requested contrast)

| dim | median F0 | median M+1 | Cliff's delta | SDS | CI95 excludes 0 | CI95 |
| --- | --- | --- | --- | --- | --- | --- |
| gamma2 | 0.02987 | 0.03362 | 0.03125 | 0.3125 | no | [-0.5625, 0.625] |
| P | 0.01578 | 0.01979 | 0.2188 | 2.188 | no | [-0.4062, 0.75] |
| D | 0.01835 | 0.01264 | -0.1562 | 1.562 | no | [-0.75, 0.4375] |
| A | 49 | 28.5 | -0.3125 | 3.125 | no | [-0.8125, 0.3125] |
| F | 21 | 16 | -0.25 | 2.5 | no | [-0.7812, 0.375] |
| S | 21.82 | 17.26 | -0.2812 | 2.812 | no | [-0.8438, 0.3438] |

## The adjacent window — M+1 vs M-1

| dim | median M-1 | median M+1 | Cliff's delta | SDS | CI95 excludes 0 | CI95 |
| --- | --- | --- | --- | --- | --- | --- |
| gamma2 | 0.06106 | 0.03362 | -0.1875 | 1.875 | no | [-0.6875, 0.3333] |
| P | 0.01963 | 0.01979 | 0.02083 | 0.2083 | no | [-0.6042, 0.6458] |
| D | 0.01342 | 0.01264 | -0.04167 | 0.4167 | no | [-0.625, 0.5417] |
| A | 29 | 28.5 | -0.02083 | 0.2083 | no | [-0.6042, 0.5625] |
| F | 14.5 | 16 | 0.2083 | 2.083 | no | [-0.3542, 0.75] |
| S | 17.99 | 17.26 | -0.1458 | 1.458 | no | [-0.6667, 0.3958] |

Cliff's delta > 0 means M+1 larger. The same comparison against the three pre windows pooled is in the JSON (`deltas[dim].vs_pre_pool`).

## Reading (descriptive — not a claim of damage)

Against the **first-images** epoch `F0` the long baseline clears the bootstrap CI in **0** of the six dimensions; against the adjacent `M-1` it is **0** of six. The `F0` baseline carries only 8 chips, seven years and a different scattering state before the event, so its intervals are wide and the `F0 -> M+1` step conflates long-term drift, season and the collapse. The month axis therefore does **not** by itself reproduce a damage signal; it shows the full spread of the record and where the collapse sits inside it.

Both sides live on the *same* structure, and the post state is **bare sand** whose intensity is moisture-driven, so season, weather, pipeline generation and the cleared lot move with the collapse date as well; the panels are read as that co-variation, **never** as damage.

```bash
python3 code/cts/cts_ocv_channels_csv.py    # build + verify the table
python3 code/cts/fig_cts_ocv_months.py      # this figure
```
