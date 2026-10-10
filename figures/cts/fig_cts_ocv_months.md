# Champlain Towers South, Surfside (FL) — the echo-mask dimensions from the first images through the 24 months before the collapse to the month after it

**Site** Champlain Towers South, Surfside (FL) · **event** 2021-06-24 · **chips** 252 (windows F0 + M-24..M+1)

```
x     = [gamma2, P, D, A, F, S]      echo-mask vector
mask  = intensity >= 0.30*peak AND >= 5.0*np.median(intensity), A >= 2
F0    [2015-09-21, 2015-11-01)       pre-collapse (first images available)
M-24  [2019-06-24, 2019-07-24)       pre-collapse (healthy)  first record month
M-1   [2021-05-24, 2021-06-24)       pre-collapse (healthy)  last
M+1   [2021-06-24, 2021-07-24)       post-collapse (month after)
```

## Why the baseline is 2015 and not 2014

Sentinel-1 data for **this** footprint do not exist in 2014. Sentinel-1A was launched on 2014-04-03, and the suitable orbit (ASC rel-48 IW3) has its **first acquisition on 2015-09-21** (`data/cts/cts_acquisition_census.json`; the Step-0 census already queried from 2014-10-01 and found nothing usable earlier — the descending stratum carries only 6 acquisitions in six years). `F0` is therefore the record-start epoch that carries the first images that actually exist, and is read here as the *first-images* baseline. No 2014 window is invented.

## F0 and the 24 pre-collapse months

`F0` opens the record where the first Sentinel-1 images of this footprint exist; the 24 rolling one-month buckets anchor at the collapse date and cover the two years up to the event, and the collapse date itself opens `M+1`. The pre side is therefore the 24 `M-*` windows (plus the `F0` baseline) and the post side is `M+1` alone. This is a companion of `fig_cts_ocv_paper.py` (figures A–F): it adds one figure and does not touch the A–F pin contract.

| window | range | state | chips | first day | last day | days |
| --- | --- | --- | --- | --- | --- | --- |
| F0 | [2015-09-21, 2015-11-01) | pre | 8 | 2015-09-21 | 2015-10-15 | 2 |
| M-24 | [2019-06-24, 2019-07-24) | pre | 8 | 2019-07-02 | 2019-07-14 | 2 |
| M-23 | [2019-07-24, 2019-08-24) | pre | 12 | 2019-07-26 | 2019-08-19 | 3 |
| M-22 | [2019-08-24, 2019-09-24) | pre | 12 | 2019-08-31 | 2019-09-12 | 3 |
| M-21 | [2019-09-24, 2019-10-24) | pre | 12 | 2019-09-24 | 2019-10-18 | 3 |
| M-20 | [2019-10-24, 2019-11-24) | pre | 12 | 2019-10-30 | 2019-11-23 | 3 |
| M-19 | [2019-11-24, 2019-12-24) | pre | 8 | 2019-12-05 | 2019-12-17 | 2 |
| M-18 | [2019-12-24, 2020-01-24) | pre | 12 | 2019-12-29 | 2020-01-22 | 3 |
| M-17 | [2020-01-24, 2020-02-24) | pre | 8 | 2020-02-03 | 2020-02-15 | 2 |
| M-16 | [2020-02-24, 2020-03-24) | pre | 12 | 2020-02-27 | 2020-03-22 | 3 |
| M-15 | [2020-03-24, 2020-04-24) | pre | 8 | 2020-04-03 | 2020-04-15 | 2 |
| M-14 | [2020-04-24, 2020-05-24) | pre | 12 | 2020-04-27 | 2020-05-21 | 3 |
| M-13 | [2020-05-24, 2020-06-24) | pre | 8 | 2020-06-02 | 2020-06-14 | 2 |
| M-12 | [2020-06-24, 2020-07-24) | pre | 12 | 2020-06-26 | 2020-07-20 | 3 |
| M-11 | [2020-07-24, 2020-08-24) | pre | 8 | 2020-08-01 | 2020-08-13 | 2 |
| M-10 | [2020-08-24, 2020-09-24) | pre | 8 | 2020-08-25 | 2020-09-06 | 2 |
| M-9 | [2020-09-24, 2020-10-24) | pre | 8 | 2020-09-30 | 2020-10-12 | 2 |
| M-8 | [2020-10-24, 2020-11-24) | pre | 8 | 2020-10-24 | 2020-11-05 | 2 |
| M-7 | [2020-11-24, 2020-12-24) | pre | 12 | 2020-11-29 | 2020-12-23 | 3 |
| M-6 | [2020-12-24, 2021-01-24) | pre | 4 | 2021-01-04 | 2021-01-04 | 1 |
| M-5 | [2021-01-24, 2021-02-24) | pre | 12 | 2021-01-28 | 2021-02-21 | 3 |
| M-4 | [2021-02-24, 2021-03-24) | pre | 8 | 2021-03-05 | 2021-03-17 | 2 |
| M-3 | [2021-03-24, 2021-04-24) | pre | 12 | 2021-03-29 | 2021-04-22 | 3 |
| M-2 | [2021-04-24, 2021-05-24) | pre | 8 | 2021-05-04 | 2021-05-16 | 2 |
| M-1 | [2021-05-24, 2021-06-24) | pre | 12 | 2021-05-28 | 2021-06-21 | 3 |
| M+1 | [2021-06-24, 2021-07-24) | post | 8 | 2021-07-03 | 2021-07-15 | 2 |

## Median of each dimension per month

| month | n | gamma2 | P | D | A | F | S |
| --- | --- | --- | --- | --- | --- | --- | --- |
| F0 | 8 | 0.02987 | 0.01578 | 0.01835 | 49 | 21 | 21.82 |
| M-24 | 8 | 0.07587 | 0.03217 | 0.01727 | 36.5 | 15.5 | 17.34 |
| M-23 | 12 | 0.03668 | 0.01679 | 0.008986 | 23.5 | 10 | 18.84 |
| M-22 | 12 | 0.02951 | 0.03214 | 0.01321 | 32.5 | 14.5 | 13.99 |
| M-21 | 12 | 0.07371 | 0.02007 | 0.01013 | 35.5 | 19 | 16.75 |
| M-20 | 12 | 0.06617 | 0.02846 | 0.0119 | 25 | 12.5 | 23.79 |
| M-19 | 8 | 0.0933 | 0.02302 | 0.0123 | 36 | 16.5 | 19.05 |
| M-18 | 12 | 0.01773 | 0.03115 | 0.007901 | 26.5 | 13 | 26.73 |
| M-17 | 8 | 0.1115 | 0.02432 | 0.01284 | 23.5 | 11 | 15.78 |
| M-16 | 12 | 0.08674 | 0.02569 | 0.009125 | 22.5 | 13 | 15.15 |
| M-15 | 8 | 0.09278 | 0.02591 | 0.008328 | 17 | 9.5 | 14.4 |
| M-14 | 12 | 0.0505 | 0.02136 | 0.007506 | 14 | 7 | 18.86 |
| M-13 | 8 | 0.05186 | 0.03145 | 0.009848 | 19.5 | 11 | 21.21 |
| M-12 | 12 | 0.0247 | 0.02095 | 0.01425 | 32 | 14 | 22 |
| M-11 | 8 | 0.02732 | 0.03089 | 0.01106 | 33 | 16.5 | 14.67 |
| M-10 | 8 | 0.09453 | 0.02043 | 0.007012 | 23.5 | 9.5 | 18.01 |
| M-9 | 8 | 0.04844 | 0.0303 | 0.01071 | 37 | 19 | 19.08 |
| M-8 | 8 | 0.03769 | 0.02977 | 0.00812 | 31 | 13.5 | 20.36 |
| M-7 | 12 | 0.03524 | 0.02663 | 0.0111 | 24.5 | 11 | 18.67 |
| M-6 | 4 | 0.1536 | 0.04165 | 0.00551 | 15.5 | 8 | 18.8 |
| M-5 | 12 | 0.07324 | 0.03338 | 0.007436 | 28.5 | 13.5 | 17.71 |
| M-4 | 8 | 0.08128 | 0.0285 | 0.0121 | 34 | 13 | 25.39 |
| M-3 | 12 | 0.0281 | 0.02464 | 0.01136 | 31 | 16 | 17.07 |
| M-2 | 8 | 0.0486 | 0.0187 | 0.009964 | 39.5 | 16.5 | 17.84 |
| M-1 | 12 | 0.06106 | 0.01963 | 0.01342 | 29 | 14.5 | 17.99 |
| M+1 | 8 | 0.03362 | 0.01979 | 0.01264 | 28.5 | 16 | 17.26 |

`baseline` = the `F0` (first-images) medians — the 1.0 line of panel f: `gamma2` 0.02987, `P` 0.01578, `D` 0.01835, `A` 49, `F` 21, `S` 21.82.

## The 6D fingerprint at the cut (month / baseline)

| month | gamma2 | P | D | A | F | S |
| --- | --- | --- | --- | --- | --- | --- |
| F0 | 1 | 1 | 1 | 1 | 1 | 1 |
| M-24 | 2.54 | 2.039 | 0.9414 | 0.7449 | 0.7381 | 0.7946 |
| M-12 | 0.8269 | 1.327 | 0.7766 | 0.6531 | 0.6667 | 1.008 |
| M-1 | 2.044 | 1.244 | 0.7317 | 0.5918 | 0.6905 | 0.8245 |
| M+1 | 1.126 | 1.254 | 0.6887 | 0.5816 | 0.7619 | 0.7909 |

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

Cliff's delta > 0 means M+1 larger. The same comparison against all 24 pre months pooled is in the JSON (`deltas[dim].vs_pre_pool`).

## Reading (descriptive — not a claim of damage)

The 24 pre-collapse months carry **236** chips (median 10 per month) and the first post month **8**, so the monthly intervals are wide and this is a descriptive companion, not a test.

Relative to the **first-images** epoch `F0` (= 1.0) the **last** pre month M-1 sits at `gamma2` 2.044, `P` 1.244, `D` 0.7317, `A` 0.5918, `F` 0.6905, `S` 0.8245; the first post month M+1 sits at `gamma2` 1.126, `P` 1.254, `D` 0.6887, `A` 0.5816, `F` 0.7619, `S` 0.7909. The pre-collapse record is **not** a flat line — season, weather, pipeline generation and the bare-sand post state move the six dimensions month to month well before the event, and the `F0 -> M+1` step spans 6 years and conflates that drift with the collapse. Against the **first-images** baseline `F0` the long contrast clears the bootstrap CI in **0** of the six dimensions; against the adjacent `M-1` it is **0** of six. The month axis therefore does **not** by itself reproduce a damage signal; it shows the full spread of the record and where the collapse sits inside it.

Both sides live on the *same* structure, and the post state is **bare sand** whose intensity is moisture-driven, so season, weather, pipeline generation and the cleared lot move with the collapse date as well; the panels are read as that co-variation, **never** as damage.

```bash
python3 code/cts/cts_ocv_channels_csv.py    # build + verify the table
python3 code/cts/fig_cts_ocv_months.py      # this figure
```
