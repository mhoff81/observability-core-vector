# Haltern AV9 Wind Park (Kreis Recklinghausen, NRW) — the orbit split of the echo-mask vector

**Site** Haltern AV9 Wind Park (Kreis Recklinghausen, NRW) · **asset** Haltern AV9 Wind Turbine (Nordex N149/4.5) · **event** 2021-09-29 (Wednesday) · **record** 2021-01-01 -> 2021-10-17 (120 acquisition days, 120 unique chips, 168 committed windows) · **figure** `figures/haltern/fig_haltern_ocv_orbits.png` (600 dpi) · **result** `data/haltern/fig_haltern_ocv_orbits.json`

**What this figure is.** The record's four Sentinel-1 pass geometries fall into two directions that **alternate in blocks of two days**, so no acquisition day carries both a usable ascending and a usable descending echo — all 120 of the 120 days carry a single direction — and a *per-day* orbit difference therefore does not exist. What does exist is the **trailing-window difference**

```
delta_orbit(t) = median(DESC in window(t)) - median(ASC in window(t))
```

with `window(t)` the 8-week window ending at the echo-bearing day `t` (at least 2 chips of each orbit, no interpolation), drawn on the twin axis of each panel beside the two orbits' own chips; the dashed line is the **pooled** difference over the whole record. A line joins two points only when they are one 3-day cadence step apart, so the nine holes of the record are gaps in the series and are never bridged. Panels (g) and (h) carry the pooled/pre-only contrast, the composition of the record and the two Fisher tests.

**What it is not.** Median(DESC) - median(ASC) over the echo-bearing chips: a descriptive pass-geometry difference, never a damage or orbit verdict. The two orbits are disjoint sets of acquisition days, so the chips are not paired and no causal reading is available. The record's four post-event echo chips are two ascending and two descending passes, so the pre-only column drops a balanced pair: the pooled difference and the before-the-collapse difference differ by those four chips alone, and printing both shows how little of the pooled number is state. The split adds no number to the A-E pin contract of `fig_haltern_ocv_paper.py`.

## The four pass geometries behind the labels

The record's passes are ASC is IW1/IW3 on relative orbit 15 + 88 (15:2x + 17:2x UT, VV x 33 + VH x 29); DESC is IW1/IW3 on relative orbit 37 + 139 (05:2x UT, VH x 31 + VV x 27).

| orbit | chips | swath | relative orbit | UT hour | VV / VH chips | label in the export |
| --- | --- | --- | --- | --- | --- | --- |
| ASC | 62 | IW1/IW3 | 15 + 88 | 15:2x + 17:2x | 33 / 29 | afternoon |
| DESC | 58 | IW1/IW3 | 37 + 139 | 05:2x | 27 / 31 | morning |

The chain is checked on every chip: the channel table's `orbit_direction`/`pass_label` agree with the committed measurement extract (120 of 120 chips present, **0 orbit mismatches**, **0 label mismatches**) and with the resolved burst record (`reference/haltern_bursts.json`, 120 resolved bursts, **0 orbit mismatches**); the 120 chips map to 120 distinct bursts. Which pass a chip is, is re-derived here and not assumed.

## The pooled orbit split

| channel | what it is | n ASC | n DESC | median ASC | median DESC | Δ (DESC − ASC) | CI 2.5% | CI 97.5% | Cliff's δ | CI excludes 0 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gamma2 | masked coherence | 33 | 33 | 0.407 | 0.567 | 0.160 | -0.162 | 0.478 | 0.179 | no |
| P | mask persistence | 33 | 33 | 0.091 | 0.091 | 0.000 | -0.014 | 0.010 | -0.071 | no |
| D | mask density A / bbox | 33 | 33 | 0.500 | 0.450 | -0.050 | -0.167 | 0.267 | 0.019 | no |
| F | 8-connected components | 33 | 33 | 2.000 | 2.000 | 0.000 | -1.000 | 1.000 | 0.049 | no |
| S | \|\|centroid - peak\|\|, px | 33 | 33 | 0.943 | 0.949 | 0.007 | -0.559 | 0.584 | 0.062 | no |
| A | echo pixels n_masked | 33 | 33 | 3.000 | 4.000 | 1.000 | -1.000 | 2.000 | 0.224 | no |

The pooled difference is **gamma2 0.1595** with a seeded 95% CI of [-0.162, 0.478] over 33 ascending and 33 descending echo chips, and **P 0.0000** with [-0.014, 0.010]. Every channel's CI spans zero — the widest is gamma2, whose upper end sits at 0.478 against only 33 descending chips — so on this record the mask vector is **orbit-invariant**, and the three mask-geometry channels (A 1.0000 (CI -1.000 .. 2.000), F 0.0000 (CI -1.000 .. 1.000), S 0.0065 (CI -0.559 .. 0.584)) plus D -0.0500 point the same way. Cliff's delta is called with the ascending chips as reference, so it carries the same sign as the bars and as the difference.

No channel has a pooled CI that excludes zero; the 10 000-resample bootstrap is seeded (7) and resamples each orbit independently, because the two orbits are disjoint sets of acquisition days and the chips are not paired.

## The same split before the collapse

| channel | what it is | n ASC | n DESC | median ASC | median DESC | Δ (DESC − ASC) | CI 2.5% | CI 97.5% | Cliff's δ | CI excludes 0 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gamma2 | masked coherence | 31 | 31 | 0.407 | 0.567 | 0.160 | -0.272 | 0.520 | 0.176 | no |
| P | mask persistence | 31 | 31 | 0.091 | 0.089 | -0.002 | -0.015 | 0.011 | -0.080 | no |
| D | mask density A / bbox | 31 | 31 | 0.500 | 0.450 | -0.050 | -0.233 | 0.267 | -0.004 | no |
| F | 8-connected components | 31 | 31 | 2.000 | 2.000 | 0.000 | -1.000 | 1.000 | 0.077 | no |
| S | \|\|centroid - peak\|\|, px | 31 | 31 | 0.791 | 0.949 | 0.159 | -0.495 | 0.601 | 0.092 | no |
| A | echo pixels n_masked | 31 | 31 | 3.000 | 4.000 | 1.000 | 0.000 | 2.000 | 0.242 | no |

The record's four post-event echo chips are two ascending and two descending passes, so the pre-only column drops a balanced pair: the pooled difference and the before-the-collapse difference differ by those four chips alone, and printing both shows how little of the pooled number is state. 4 of the record's 66 echo chips are post-event — 2 on the ascending and 2 on the descending geometry — so the pre-only column is the same contrast with a balanced pair dropped. gamma2 moves from 0.1595 to 0.1595 and P from 0.0000 to -0.0017, and gamma2's CI from [-0.162, 0.478] to [-0.272, 0.520] — both span zero, and the whole mask geometry is unchanged. Neither column is a damage reading.

## The trailing-window series of every channel

The series is drawn at the 61 echo-bearing days whose 8-week trailing window holds at least 2 chips of each orbit (``2021-02-01``, ``2021-02-03``, ``2021-02-04``, ``2021-02-06``, ``2021-02-07``, ``2021-02-10``, ``2021-02-12``, ``2021-02-15``, ``2021-02-16``, ``2021-02-19``, ``2021-03-02``, ``2021-03-06``, ``2021-03-12``, ``2021-03-14``, ``2021-03-15``, ``2021-03-17``, ``2021-04-01``, ``2021-04-05``, ``2021-04-08``, ``2021-04-11``, ``2021-04-13``, ``2021-04-16``, ``2021-04-17``, ``2021-05-01``, ``2021-05-07``, ``2021-05-10``, ``2021-05-13``, ``2021-05-14``, ``2021-06-01``, ``2021-06-04``, ``2021-06-09``, ``2021-06-10``, ``2021-06-12``, ``2021-06-15``, ``2021-06-16``, ``2021-06-18``, ``2021-06-19``, ``2021-07-01``, ``2021-07-07``, ``2021-07-09``, ``2021-07-12``, ``2021-07-13``, ``2021-07-16``, ``2021-07-18``, ``2021-07-19``, ``2021-08-03``, ``2021-08-06``, ``2021-08-08``, ``2021-08-14``, ``2021-08-15``, ``2021-08-18``, ``2021-09-01``, ``2021-09-04``, ``2021-09-13``, ``2021-09-14``, ``2021-09-16``, ``2021-09-17``, ``2021-10-02``, ``2021-10-08``, ``2021-10-10``, ``2021-10-17``), one point per day and channel: 366 points in all, 24 of them at the four post-event echo days.

| day | n ASC | n DESC | state | gamma2 | P | D | F | S | A |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2021-02-01 | 2 | 4 | pre | 0.648 | 0.017 | -0.181 | 0.500 | 0.109 | 1.000 |
| 2021-02-03 | 2 | 5 | pre | 0.670 | 0.019 | -0.095 | 0.500 | -0.120 | 1.000 |
| 2021-02-04 | 2 | 6 | pre | 0.648 | 0.017 | -0.081 | 0.500 | 0.018 | 1.500 |
| 2021-02-06 | 3 | 6 | pre | 0.795 | 0.011 | 0.086 | 0.000 | -0.445 | 2.500 |
| 2021-02-07 | 4 | 6 | pre | 0.648 | 0.011 | 0.086 | 0.500 | -0.112 | 2.500 |
| 2021-02-10 | 4 | 7 | pre | 0.625 | 0.010 | 0.100 | 0.500 | -0.250 | 2.000 |
| 2021-02-12 | 5 | 7 | pre | 0.748 | 0.010 | 0.100 | 0.000 | -0.264 | 2.000 |
| 2021-02-15 | 5 | 8 | pre | 0.770 | 0.011 | 0.100 | -0.500 | -0.376 | 1.500 |
| 2021-02-16 | 5 | 9 | pre | 0.748 | 0.010 | 0.100 | -1.000 | -0.488 | 1.000 |
| 2021-02-19 | 6 | 9 | pre | 0.613 | 0.014 | 0.100 | -0.500 | -0.315 | 1.000 |
| 2021-03-02 | 7 | 9 | pre | 0.478 | 0.010 | 0.100 | 0.000 | -0.141 | 1.000 |
| 2021-03-06 | 6 | 9 | pre | 0.261 | 0.013 | 0.125 | -0.500 | -0.381 | 1.500 |
| 2021-03-12 | 6 | 9 | pre | 0.195 | 0.013 | 0.125 | -0.500 | -0.381 | 1.500 |
| 2021-03-14 | 7 | 8 | pre | 0.444 | 0.011 | 0.312 | -1.000 | -0.525 | 1.000 |
| 2021-03-15 | 8 | 7 | pre | 0.195 | 0.015 | 0.500 | -0.500 | -0.381 | 0.000 |
| 2021-03-17 | 8 | 8 | pre | 0.271 | 0.011 | 0.312 | -0.500 | -0.381 | 0.500 |
| 2021-04-01 | 8 | 7 | pre | 0.143 | 0.000 | 0.417 | -0.500 | -0.354 | 0.000 |
| 2021-04-05 | 6 | 7 | pre | 0.143 | 0.006 | 0.167 | -0.500 | -0.250 | 0.000 |
| 2021-04-08 | 7 | 6 | pre | 0.417 | -0.008 | 0.083 | -1.000 | -0.443 | 0.500 |
| 2021-04-11 | 6 | 7 | pre | 0.117 | -0.011 | -0.333 | -0.500 | -0.221 | 1.000 |
| 2021-04-13 | 7 | 6 | pre | 0.307 | -0.016 | -0.192 | 0.500 | -0.178 | 1.500 |
| 2021-04-16 | 7 | 6 | pre | 0.396 | -0.016 | -0.192 | 0.500 | -0.148 | 2.000 |
| 2021-04-17 | 6 | 7 | pre | 0.514 | -0.015 | -0.133 | 0.500 | -0.266 | 2.000 |
| 2021-05-01 | 6 | 7 | pre | 0.514 | -0.015 | -0.133 | 0.500 | -0.266 | 2.000 |
| 2021-05-07 | 7 | 6 | pre | 0.474 | -0.016 | -0.075 | 0.000 | -0.036 | 2.000 |
| 2021-05-10 | 6 | 6 | pre | 0.425 | -0.016 | -0.225 | 0.500 | 0.476 | 2.000 |
| 2021-05-13 | 6 | 5 | pre | 0.335 | -0.011 | -0.317 | 0.500 | 0.751 | 2.000 |
| 2021-05-14 | 7 | 5 | pre | 0.285 | -0.008 | -0.400 | 1.000 | 0.793 | 2.000 |
| 2021-06-01 | 7 | 4 | pre | -0.096 | -0.008 | -0.307 | 1.500 | 1.746 | 2.000 |
| 2021-06-04 | 6 | 5 | pre | -0.229 | -0.004 | -0.483 | 1.000 | 0.896 | 2.000 |
| 2021-06-09 | 5 | 5 | pre | -0.703 | 0.000 | -0.733 | 1.000 | 1.000 | 2.000 |
| 2021-06-10 | 5 | 6 | pre | -0.704 | -0.003 | -0.717 | 1.000 | 0.937 | 2.000 |
| 2021-06-12 | 6 | 5 | pre | -0.662 | 0.002 | -0.567 | 1.000 | 1.000 | 1.000 |
| 2021-06-15 | 6 | 5 | pre | -0.662 | 0.002 | -0.533 | 1.000 | 0.874 | 0.000 |
| 2021-06-16 | 6 | 6 | pre | -0.701 | -0.001 | -0.505 | 1.000 | 1.248 | 1.000 |
| 2021-06-18 | 7 | 6 | pre | -0.658 | -0.003 | -0.421 | 1.000 | 1.248 | 0.000 |
| 2021-06-19 | 8 | 6 | pre | -0.375 | -0.001 | -0.380 | 1.000 | 1.197 | 0.000 |
| 2021-07-01 | 8 | 6 | pre | -0.267 | -0.001 | -0.171 | 0.500 | 1.052 | 0.000 |
| 2021-07-07 | 8 | 5 | pre | -0.229 | -0.010 | -0.143 | 0.500 | 0.515 | 0.000 |
| 2021-07-09 | 7 | 6 | pre | -0.053 | -0.018 | 0.158 | -0.500 | -0.077 | -0.500 |
| 2021-07-12 | 7 | 6 | pre | -0.053 | -0.035 | 0.158 | -0.500 | -0.077 | -0.500 |
| 2021-07-13 | 8 | 6 | pre | -0.130 | -0.027 | 0.158 | -0.500 | -0.077 | -0.500 |
| 2021-07-16 | 8 | 7 | pre | -0.129 | -0.023 | 0.024 | 0.000 | 0.256 | 0.000 |
| 2021-07-18 | 9 | 7 | pre | -0.052 | -0.015 | 0.024 | 0.000 | 0.256 | 0.000 |
| 2021-07-19 | 10 | 7 | pre | -0.045 | -0.018 | 0.024 | 0.000 | 0.094 | 0.000 |
| 2021-08-03 | 9 | 7 | pre | -0.013 | -0.015 | 0.000 | 0.000 | -0.068 | -1.000 |
| 2021-08-06 | 10 | 5 | pre | 0.251 | -0.020 | 0.024 | 0.000 | -0.248 | 0.000 |
| 2021-08-08 | 9 | 6 | pre | 0.229 | -0.012 | 0.012 | 0.000 | -0.298 | -0.500 |
| 2021-08-14 | 8 | 5 | pre | 0.115 | -0.023 | 0.017 | 0.000 | -0.272 | 0.500 |
| 2021-08-15 | 7 | 6 | pre | 0.020 | -0.017 | 0.000 | 0.000 | -0.298 | 0.500 |
| 2021-08-18 | 8 | 6 | pre | -0.077 | -0.013 | 0.000 | 0.000 | -0.290 | 0.500 |
| 2021-09-01 | 7 | 7 | pre | 0.031 | -0.019 | 0.000 | 0.000 | -0.340 | 0.000 |
| 2021-09-04 | 7 | 6 | pre | 0.020 | -0.008 | 0.000 | 0.000 | -0.282 | -1.500 |
| 2021-09-13 | 4 | 6 | pre | -0.164 | 0.006 | 0.017 | 0.000 | 0.078 | -1.000 |
| 2021-09-14 | 3 | 7 | pre | -0.506 | 0.023 | 0.000 | 0.000 | 0.495 | -1.000 |
| 2021-09-16 | 4 | 7 | pre | -0.302 | 0.016 | 0.017 | 0.000 | 0.120 | -1.500 |
| 2021-09-17 | 5 | 7 | pre | -0.171 | 0.019 | 0.000 | 0.000 | 0.411 | -1.000 |
| 2021-10-02 | 4 | 7 | post | -0.134 | 0.028 | -0.031 | 1.000 | 0.772 | -1.500 |
| 2021-10-08 | 4 | 7 | post | -0.134 | 0.028 | -0.031 | 1.000 | 0.772 | -1.500 |
| 2021-10-10 | 5 | 6 | post | -0.092 | 0.027 | 0.343 | 0.000 | 0.571 | -2.500 |
| 2021-10-17 | 5 | 5 | post | -0.226 | 0.018 | -0.014 | 1.000 | 0.768 | -2.000 |

| channel | points | min Δ | max Δ | max \|Δ\| | sign changes | distinct values |
| --- | --- | --- | --- | --- | --- | --- |
| gamma2 | 61 | -0.704 | 0.795 | 0.795 | 5 | 51 |
| P | 61 | -0.035 | 0.028 | 0.035 | 4 | 39 |
| D | 61 | -0.733 | 0.500 | 0.733 | 6 | 37 |
| F | 61 | -1.000 | 1.500 | 1.500 | 4 | 6 |
| S | 61 | -0.525 | 1.746 | 1.746 | 8 | 50 |
| A | 61 | -2.500 | 2.500 | 2.500 | 3 | 11 |

**The series moves where the pooled number cannot say so.** D swings between -0.7333 and 0.5000 (6 sign changes), S between -0.5246 and 1.7458 (8 sign changes), F between -1.0000 and 1.5000, and A down to -2.5000 — all of them channels whose pooled CI spans zero. A pooled zero is an average over orbit-dependent windows, not the absence of an orbit effect; and with the two directions never sharing a day, a per-day reading of any of it does not exist.

**The series keeps its points past the event.** The record's four post-event echoes fall on two ascending and two descending days, so the window of the first of them (`2021-10-02`, ASC) already holds 4 ascending and 7 descending echo chips and a window difference exists there: the drawn line runs into the post-event stretch instead of stopping at the collapse.

The 9 steps of the record longer than the 3-day cadence (11, 12, 14, 15 days) are holes: the drawn line breaks at each of them and no value is interpolated across one.

| after | before | days | crosses the event |
| --- | --- | --- | --- |
| 2021-01-17 | 2021-02-01 | 15 | no |
| 2021-02-19 | 2021-03-02 | 11 | no |
| 2021-03-18 | 2021-04-01 | 14 | no |
| 2021-04-17 | 2021-05-01 | 14 | no |
| 2021-05-17 | 2021-06-01 | 15 | no |
| 2021-06-19 | 2021-07-01 | 12 | no |
| 2021-07-19 | 2021-08-02 | 14 | no |
| 2021-08-18 | 2021-09-01 | 14 | no |
| 2021-09-17 | 2021-10-01 | 14 | yes |

## Where the echoes sit over the two pass directions

| orbit | state | chips | with echo | without echo | echo rate |
| --- | --- | --- | --- | --- | --- |
| ASC | pre | 56 | 31 | 25 | 0.554 |
| ASC | post | 6 | 2 | 4 | 0.333 |
| ASC | all | 62 | 33 | 29 | 0.532 |
| DESC | pre | 52 | 31 | 21 | 0.596 |
| DESC | post | 6 | 2 | 4 | 0.333 |
| DESC | all | 58 | 33 | 25 | 0.569 |

The coverage of the two directions is a **distribution, not a contrast**: pooled Fisher p = 0.717 (33+29 against 33+25), and within an orbit pre-vs-post p = 0.405 (ASC) and 0.387 (DESC). Both directions straddle the collapse date and the four post-event echo chips are two ascending and two descending passes: this block describes where the echoes sit over the four pass geometries, it is not an orbit contrast of the states.

Panel (h) draws the same block: the record's 120 acquisition days and 120 chips, 66 of them echoes (62 pre, 4 post), 62 chips on the ascending geometry and 58 on the descending one.

## What else moves with the pass geometry

What else moves with the pass geometry: the per-orbit medians of the export's co-variates. Ascending passes are the 15:2x and 17:2x UT ones, descending the 05:2x UT ones, so temperature, wind and intensity move with the orbit by construction — they would have to be held fixed before an orbit split of the mask channels meant anything physical. The medians below are recomputed twice — from the channel table and from the committed measurement extract — and pinned against each other.

All 120 chips:

| co-variate | label | median ASC | median DESC | Δ (DESC − ASC) |
| --- | --- | --- | --- | --- |
| temperature_c | temperature (C) | 13.950000 | 8.850000 | -5.100000 |
| wind_speed_ms | wind (m/s) | 3.750000 | 2.615000 | -1.135000 |
| wind_gust_ms | gust (m/s) | 8.400000 | 5.100000 | -3.300000 |
| intensity | intensity | 6566.392880 | 4300.653936 | -2265.738944 |
| coherence_gamma2 | gamma2 (export) | 0.008287 | 0.131980 | 0.123693 |
| precipitation_mm | precip (mm) | 0.000000 | 0.000000 | 0.000000 |

The 66 echo-bearing chips:

| co-variate | label | median ASC | median DESC | Δ (DESC − ASC) |
| --- | --- | --- | --- | --- |
| temperature_c | temperature (C) | 14.700000 | 8.600000 | -6.100000 |
| wind_speed_ms | wind (m/s) | 4.000000 | 2.610000 | -1.390000 |
| wind_gust_ms | gust (m/s) | 9.400000 | 5.000000 | -4.400000 |
| intensity | intensity | 5852.264499 | 3748.981934 | -2103.282564 |
| coherence_gamma2 | gamma2 (export) | 0.008902 | 0.163193 | 0.154292 |
| precipitation_mm | precip (mm) | 0.000000 | 0.000000 | 0.000000 |

The ascending passes are the warmer ones (-5.1000 C on the record, -6.1000 C among the echo chips) and the darker ones in SAR intensity (-2103.3); no chip of the record carries rain (`precipitation_mm` is 0.0 on every chip). Whatever an orbit split of the mask channels would mean physically, these are the terms that move with it.

## The record: one direction per day, and its nine holes

The record is 120 acquisition days and 168 committed windows (48 of them a repeat of a payload already in the table), which compile to 120 unique chips: 66 carry an echo mask (62 pre, 4 post) and 54 do not. The two directions alternate in blocks of two days and never share one: all 120 of the 120 days carry a single direction, so the two orbits of the record never meet inside one day and this figure's ordinate never has two orbits to subtract on it.

| day | state | orbits | windows | chips | repeats | with echo | without echo |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 2021-01-01 | pre | ASC | 2 | 1 | 1 | 0 | 1 |
| 2021-01-02 | pre | ASC | 2 | 1 | 1 | 0 | 1 |
| 2021-01-04 | pre | DESC | 2 | 1 | 1 | 0 | 1 |
| 2021-01-05 | pre | DESC | 2 | 1 | 1 | 1 | 0 |
| 2021-01-07 | pre | ASC | 2 | 1 | 1 | 1 | 0 |
| 2021-01-08 | pre | ASC | 2 | 1 | 1 | 0 | 1 |
| 2021-01-10 | pre | DESC | 2 | 1 | 1 | 1 | 0 |
| 2021-01-11 | pre | DESC | 2 | 1 | 1 | 0 | 1 |
| 2021-01-13 | pre | ASC | 2 | 1 | 1 | 0 | 1 |
| 2021-01-14 | pre | ASC | 2 | 1 | 1 | 0 | 1 |
| 2021-01-16 | pre | DESC | 2 | 1 | 1 | 1 | 0 |
| 2021-01-17 | pre | DESC | 2 | 1 | 1 | 1 | 0 |
| 2021-02-01 | pre | ASC | 2 | 1 | 1 | 1 | 0 |
| 2021-02-03 | pre | DESC | 2 | 1 | 1 | 1 | 0 |
| 2021-02-04 | pre | DESC | 2 | 1 | 1 | 1 | 0 |
| 2021-02-06 | pre | ASC | 2 | 1 | 1 | 1 | 0 |
| 2021-02-07 | pre | ASC | 2 | 1 | 1 | 1 | 0 |
| 2021-02-09 | pre | DESC | 2 | 1 | 1 | 0 | 1 |
| 2021-02-10 | pre | DESC | 2 | 1 | 1 | 1 | 0 |
| 2021-02-12 | pre | ASC | 2 | 1 | 1 | 1 | 0 |
| 2021-02-15 | pre | DESC | 2 | 1 | 1 | 1 | 0 |
| 2021-02-16 | pre | DESC | 2 | 1 | 1 | 1 | 0 |
| 2021-02-18 | pre | ASC | 2 | 1 | 1 | 0 | 1 |
| 2021-02-19 | pre | ASC | 2 | 1 | 1 | 1 | 0 |
| 2021-03-02 | pre | ASC | 2 | 1 | 1 | 1 | 0 |
| 2021-03-03 | pre | ASC | 2 | 1 | 1 | 0 | 1 |
| 2021-03-05 | pre | DESC | 2 | 1 | 1 | 0 | 1 |
| 2021-03-06 | pre | DESC | 2 | 1 | 1 | 1 | 0 |
| 2021-03-08 | pre | ASC | 2 | 1 | 1 | 0 | 1 |
| 2021-03-09 | pre | ASC | 2 | 1 | 1 | 0 | 1 |
| 2021-03-11 | pre | DESC | 2 | 1 | 1 | 0 | 1 |
| 2021-03-12 | pre | DESC | 2 | 1 | 1 | 1 | 0 |
| 2021-03-14 | pre | ASC | 2 | 1 | 1 | 1 | 0 |
| 2021-03-15 | pre | ASC | 2 | 1 | 1 | 1 | 0 |
| 2021-03-17 | pre | DESC | 2 | 1 | 1 | 1 | 0 |
| 2021-03-18 | pre | DESC | 2 | 1 | 1 | 0 | 1 |
| 2021-04-01 | pre | ASC | 2 | 1 | 1 | 1 | 0 |
| 2021-04-02 | pre | ASC | 2 | 1 | 1 | 0 | 1 |
| 2021-04-04 | pre | DESC | 2 | 1 | 1 | 0 | 1 |
| 2021-04-05 | pre | DESC | 2 | 1 | 1 | 1 | 0 |
| 2021-04-07 | pre | ASC | 2 | 1 | 1 | 0 | 1 |
| 2021-04-08 | pre | ASC | 2 | 1 | 1 | 1 | 0 |
| 2021-04-10 | pre | DESC | 2 | 1 | 1 | 0 | 1 |
| 2021-04-11 | pre | DESC | 2 | 1 | 1 | 1 | 0 |
| 2021-04-13 | pre | ASC | 2 | 1 | 1 | 1 | 0 |
| 2021-04-14 | pre | ASC | 2 | 1 | 1 | 0 | 1 |
| 2021-04-16 | pre | DESC | 2 | 1 | 1 | 1 | 0 |
| 2021-04-17 | pre | DESC | 2 | 1 | 1 | 1 | 0 |
| 2021-05-01 | pre | ASC | 1 | 1 | 0 | 1 | 0 |
| 2021-05-02 | pre | ASC | 1 | 1 | 0 | 0 | 1 |
| 2021-05-04 | pre | DESC | 1 | 1 | 0 | 0 | 1 |
| 2021-05-05 | pre | DESC | 1 | 1 | 0 | 0 | 1 |
| 2021-05-07 | pre | ASC | 1 | 1 | 0 | 1 | 0 |
| 2021-05-08 | pre | ASC | 1 | 1 | 0 | 0 | 1 |
| 2021-05-10 | pre | DESC | 1 | 1 | 0 | 1 | 0 |
| 2021-05-11 | pre | DESC | 1 | 1 | 0 | 0 | 1 |
| 2021-05-13 | pre | ASC | 1 | 1 | 0 | 1 | 0 |
| 2021-05-14 | pre | ASC | 1 | 1 | 0 | 1 | 0 |
| 2021-05-16 | pre | DESC | 1 | 1 | 0 | 0 | 1 |
| 2021-05-17 | pre | DESC | 1 | 1 | 0 | 0 | 1 |
| 2021-06-01 | pre | ASC | 1 | 1 | 0 | 1 | 0 |
| 2021-06-04 | pre | DESC | 1 | 1 | 0 | 1 | 0 |
| 2021-06-06 | pre | ASC | 1 | 1 | 0 | 0 | 1 |
| 2021-06-07 | pre | ASC | 1 | 1 | 0 | 0 | 1 |
| 2021-06-09 | pre | DESC | 1 | 1 | 0 | 1 | 0 |
| 2021-06-10 | pre | DESC | 1 | 1 | 0 | 1 | 0 |
| 2021-06-12 | pre | ASC | 1 | 1 | 0 | 1 | 0 |
| 2021-06-13 | pre | ASC | 1 | 1 | 0 | 0 | 1 |
| 2021-06-15 | pre | DESC | 1 | 1 | 0 | 1 | 0 |
| 2021-06-16 | pre | DESC | 1 | 1 | 0 | 1 | 0 |
| 2021-06-18 | pre | ASC | 1 | 1 | 0 | 1 | 0 |
| 2021-06-19 | pre | ASC | 1 | 1 | 0 | 1 | 0 |
| 2021-07-01 | pre | ASC | 1 | 1 | 0 | 1 | 0 |
| 2021-07-03 | pre | DESC | 1 | 1 | 0 | 0 | 1 |
| 2021-07-06 | pre | ASC | 1 | 1 | 0 | 0 | 1 |
| 2021-07-07 | pre | ASC | 1 | 1 | 0 | 1 | 0 |
| 2021-07-09 | pre | DESC | 1 | 1 | 0 | 1 | 0 |
| 2021-07-10 | pre | DESC | 1 | 1 | 0 | 0 | 1 |
| 2021-07-12 | pre | ASC | 1 | 1 | 0 | 1 | 0 |
| 2021-07-13 | pre | ASC | 1 | 1 | 0 | 1 | 0 |
| 2021-07-15 | pre | DESC | 1 | 1 | 0 | 0 | 1 |
| 2021-07-16 | pre | DESC | 1 | 1 | 0 | 1 | 0 |
| 2021-07-18 | pre | ASC | 1 | 1 | 0 | 1 | 0 |
| 2021-07-19 | pre | ASC | 1 | 1 | 0 | 1 | 0 |
| 2021-08-02 | pre | DESC | 1 | 1 | 0 | 0 | 1 |
| 2021-08-03 | pre | DESC | 1 | 1 | 0 | 1 | 0 |
| 2021-08-05 | pre | ASC | 1 | 1 | 0 | 0 | 1 |
| 2021-08-06 | pre | ASC | 1 | 1 | 0 | 1 | 0 |
| 2021-08-08 | pre | DESC | 1 | 1 | 0 | 1 | 0 |
| 2021-08-09 | pre | DESC | 1 | 1 | 0 | 0 | 1 |
| 2021-08-11 | pre | ASC | 1 | 1 | 0 | 0 | 1 |
| 2021-08-12 | pre | ASC | 1 | 1 | 0 | 0 | 1 |
| 2021-08-14 | pre | DESC | 1 | 1 | 0 | 1 | 0 |
| 2021-08-15 | pre | DESC | 1 | 1 | 0 | 1 | 0 |
| 2021-08-17 | pre | ASC | 1 | 1 | 0 | 0 | 1 |
| 2021-08-18 | pre | ASC | 1 | 1 | 0 | 1 | 0 |
| 2021-09-01 | pre | DESC | 1 | 1 | 0 | 1 | 0 |
| 2021-09-02 | pre | DESC | 1 | 1 | 0 | 0 | 1 |
| 2021-09-04 | pre | ASC | 1 | 1 | 0 | 1 | 0 |
| 2021-09-05 | pre | ASC | 1 | 1 | 0 | 0 | 1 |
| 2021-09-07 | pre | DESC | 1 | 1 | 0 | 0 | 1 |
| 2021-09-08 | pre | DESC | 1 | 1 | 0 | 0 | 1 |
| 2021-09-10 | pre | ASC | 1 | 1 | 0 | 0 | 1 |
| 2021-09-11 | pre | ASC | 1 | 1 | 0 | 0 | 1 |
| 2021-09-13 | pre | DESC | 1 | 1 | 0 | 1 | 0 |
| 2021-09-14 | pre | DESC | 1 | 1 | 0 | 1 | 0 |
| 2021-09-16 | pre | ASC | 1 | 1 | 0 | 1 | 0 |
| 2021-09-17 | pre | ASC | 1 | 1 | 0 | 1 | 0 |
| 2021-10-01 | post | DESC | 1 | 1 | 0 | 0 | 1 |
| 2021-10-02 | post | DESC | 1 | 1 | 0 | 1 | 0 |
| 2021-10-04 | post | ASC | 1 | 1 | 0 | 0 | 1 |
| 2021-10-05 | post | ASC | 1 | 1 | 0 | 0 | 1 |
| 2021-10-07 | post | DESC | 1 | 1 | 0 | 0 | 1 |
| 2021-10-08 | post | DESC | 1 | 1 | 0 | 1 | 0 |
| 2021-10-10 | post | ASC | 1 | 1 | 0 | 1 | 0 |
| 2021-10-11 | post | ASC | 1 | 1 | 0 | 0 | 1 |
| 2021-10-13 | post | DESC | 1 | 1 | 0 | 0 | 1 |
| 2021-10-14 | post | DESC | 1 | 1 | 0 | 0 | 1 |
| 2021-10-16 | post | ASC | 1 | 1 | 0 | 0 | 1 |
| 2021-10-17 | post | ASC | 1 | 1 | 0 | 1 | 0 |

## The pins and the provenance

3540 pin comparisons against five committed layers (the channel table's generated metadata, the mask-cache manifest, the committed mask cache, the committed paper result and the committed analysis reference, plus the raw measurement extract and the three resolved reference records): **3477 exact**, **63 within 1e-12** (the seeded numbers and the six-decimal record facts are pinned at 1e-06), **0 failed**. The recomputed number is in `recomputed`, the committed one in `committed`, and no artifact is written before they agree.

| committed input | layer | sha256 |
| --- | --- | --- |
| `data/haltern/haltern_ocv_channels.csv` | the channel table (this figure's rows) | `a82d716056ef7fa4…` |
| `data/haltern/haltern_windows_full.txt` | the window file (the chips) | `e837309d7e7545b0…` |
| `data/haltern/haltern_measurements_full.txt` | the raw measurement extract | `d179ff230d42ff4d…` |
| `data/haltern/haltern_windows_mask_cache.txt` | the committed mask cache (the vector) | `70031d3fe516a5cf…` |
| `data/haltern/haltern_windows_mask_cache.manifest.json` | the cache's manifest | `1a679ead9ba43ab5…` |
| `data/haltern/haltern_segments.txt` | the segment table | `9052f8abb8610c69…` |
| `data/haltern/reference/haltern_bursts.json` | the resolved bursts (the pass) | `74a9e1fd530cbe8f…` |
| `data/haltern/reference/haltern_geometry.json` | the geometry anchor | `fd71b4b3d9b3ea31…` |
| `data/haltern/reference/haltern_osm_node.json` | the OSM node extract | `a37502af41bca70d…` |

Figures: `figures/haltern/fig_haltern_ocv_orbits.png` (600 dpi). JSON: `data/haltern/fig_haltern_ocv_orbits.json`.
