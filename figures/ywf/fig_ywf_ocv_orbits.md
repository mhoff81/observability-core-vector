# Yeongdeok Wind Farm (Samgye-ri) — the orbit split of the echo-mask vector

**Site** Yeongdeok Wind Farm (Samgye-ri) · **asset** Yeongdeok Wind Turbine (Samgye-ri) · **event** 2026-02-02 (Monday) · **record** 2025-08-06 -> 2026-04-27 (32 acquisition days, 33 unique chips, 57 committed windows) · **figure** `figures/ywf/fig_ywf_ocv_orbits.png` (600 dpi) · **result** `data/ywf/fig_ywf_ocv_orbits.json`

**What this figure is.** The record's two Sentinel-1 pass geometries **alternate**: one orbit per 6-day acquisition day, so no acquisition day carries both a usable ascending and a usable descending echo — 31 of the 32 days carry a single orbit, and the one exception (`2025-08-12`) carries two echo-free chips. A *per-day* orbit difference therefore does not exist. What does exist is the **trailing-window difference**

```
delta_orbit(t) = median(DESC in window(t)) - median(ASC in window(t))
```

with `window(t)` the 8-week window ending at the echo-bearing day `t` (at least 2 chips of each orbit, no interpolation), drawn on the twin axis of each panel beside the two orbits' own chips; the dashed line is the **pooled** difference over the whole record. A line joins two points only when they are one 6-day cadence step apart, so the nine holes of the record are gaps in the series and are never bridged. Panels (g) and (h) carry the pooled/pre-only contrast, the composition of the record and the two Fisher tests.

**What it is not.** Median(DESC) - median(ASC) over the echo-bearing chips: a descriptive pass-geometry difference, never a damage or orbit verdict. The two orbits are disjoint sets of acquisition days, so the chips are not paired and no causal reading is available. The record's only post-event echo chip is an ascending pass, so the pooled difference and the before-the-collapse difference are two different contrasts: the pre-only column drops that chip and shows how much of the pooled number is state rather than orbit. The split adds no number to the A-E pin contract of `fig_ywf_ocv_paper.py`.

## The two pass geometries behind the labels

The record's passes are ASC is IW3 on relative orbit 54 (09:2x UT, VV x 11 + VH x 9); DESC is IW2 on relative orbit 61 (21:2x UT, VH x 10 + VV x 3).

| orbit | chips | swath | relative orbit | UT hour | VV / VH chips | label in the export |
| --- | --- | --- | --- | --- | --- | --- |
| ASC | 20 | IW3 | 54 | 09:2x | 11 / 9 | morning |
| DESC | 13 | IW2 | 61 | 21:2x | 3 / 10 | afternoon |

The chain is checked on every chip: the channel table's `orbit_direction`/`pass_label` agree with the committed measurement extract (33 of 33 chips present, **0 orbit mismatches**, **0 label mismatches**) and with the resolved burst record (`reference/ywf_bursts.json`, 33 resolved bursts, **0 orbit mismatches**); the 33 chips map to 33 distinct bursts. Which pass a chip is, is re-derived here and not assumed.

## The pooled orbit split

| channel | what it is | n ASC | n DESC | median ASC | median DESC | Δ (DESC − ASC) | CI 2.5% | CI 97.5% | Cliff's δ | CI excludes 0 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gamma2 | masked coherence | 10 | 7 | 0.227 | 0.815 | 0.588 | 0.020 | 0.730 | 0.543 | yes |
| P | mask persistence | 10 | 7 | 0.118 | 0.206 | 0.088 | 0.015 | 0.141 | 0.743 | yes |
| D | mask density A / bbox | 10 | 7 | 0.267 | 0.222 | -0.044 | -0.233 | 0.567 | 0.029 | no |
| F | 8-connected components | 10 | 7 | 2.000 | 2.000 | 0.000 | -1.000 | 0.000 | -0.129 | no |
| S | \|\|centroid - peak\|\|, px | 10 | 7 | 1.540 | 1.519 | -0.022 | -1.057 | 1.168 | 0.057 | no |
| A | echo pixels n_masked | 10 | 7 | 3.000 | 3.000 | 0.000 | -2.000 | 1.000 | -0.114 | no |

The pooled difference is **gamma2 0.5878** with a seeded 95% CI of [0.020, 0.730] and **P 0.0882** with [0.015, 0.141]: both exclude zero, but the gamma2 lower end sits at 0.020 against 7 descending and 10 ascending echo chips, so the headline channel is resolved only marginally. The three mask-geometry channels are **orbit-invariant on the pooled record**: A 0.0000 (CI -2.000 .. 1.000), F 0.0000 (CI -1.000 .. 0.000), S -0.0218 (CI -1.057 .. 1.168), and D is -0.0444. Cliff's delta is called with the ascending chips as reference, so it carries the same sign as the bars and as the difference.

Only gamma2, P have a pooled CI that excludes zero; the 10 000-resample bootstrap is seeded (7) and resamples each orbit independently, because the two orbits are disjoint sets of acquisition days and the chips are not paired.

## The same split before the collapse

| channel | what it is | n ASC | n DESC | median ASC | median DESC | Δ (DESC − ASC) | CI 2.5% | CI 97.5% | Cliff's δ | CI excludes 0 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gamma2 | masked coherence | 9 | 7 | 0.266 | 0.815 | 0.549 | -0.022 | 0.751 | 0.524 | no |
| P | mask persistence | 9 | 7 | 0.118 | 0.206 | 0.088 | 0.024 | 0.141 | 0.778 | yes |
| D | mask density A / bbox | 9 | 7 | 0.333 | 0.222 | -0.111 | -0.248 | 0.550 | -0.016 | no |
| F | 8-connected components | 9 | 7 | 2.000 | 2.000 | 0.000 | -1.000 | 0.000 | -0.048 | no |
| S | \|\|centroid - peak\|\|, px | 9 | 7 | 1.667 | 1.519 | -0.148 | -1.057 | 1.250 | 0.048 | no |
| A | echo pixels n_masked | 9 | 7 | 3.000 | 3.000 | 0.000 | -1.000 | 1.000 | -0.048 | no |

The record's only post-event echo chip is an ascending pass, so the pooled difference and the before-the-collapse difference are two different contrasts: the pre-only column drops that chip and shows how much of the pooled number is state rather than orbit. 1 of the record's 17 echo chips is post-event and it is the ascending (ASC) pass, so the pre-only column is the same contrast with that one chip dropped. gamma2 moves from 0.5878 to 0.5493 and its CI from [0.020, 0.730] to [-0.022, 0.751] — the lower end now spans zero — while P stays at 0.0882 with [0.024, 0.141], still clear of zero, and the whole mask geometry is unchanged. The pooled number is therefore a mixture of orbit and state; both columns are printed and neither is a damage reading.

## The trailing-window series of every channel

The series is drawn at the 10 echo-bearing days whose 8-week trailing window holds at least 2 chips of each orbit (``2025-10-17``, ``2025-10-23``, ``2025-10-29``, ``2025-11-10``, ``2025-12-04``, ``2025-12-10``, ``2025-12-16``, ``2026-01-03``, ``2026-01-15``, ``2026-01-27``), one point per day and channel: 60 points in all, every one of them pre-event.

| day | n ASC | n DESC | state | gamma2 | P | D | F | S | A |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2025-10-17 | 3 | 2 | pre | 0.133 | 0.126 | 0.408 | -0.500 | -0.567 | -0.500 |
| 2025-10-23 | 4 | 2 | pre | 0.110 | 0.115 | 0.325 | -0.500 | -0.291 | -1.000 |
| 2025-10-29 | 4 | 3 | pre | -0.003 | 0.129 | 0.150 | 0.000 | -0.252 | -0.500 |
| 2025-11-10 | 3 | 3 | pre | 0.427 | 0.141 | -0.111 | 0.000 | 0.025 | -1.000 |
| 2025-12-04 | 2 | 4 | pre | 0.071 | 0.115 | 0.236 | -0.500 | -0.516 | 0.000 |
| 2025-12-10 | 2 | 4 | pre | 0.203 | 0.115 | 0.219 | -0.500 | -0.646 | 0.000 |
| 2025-12-16 | 2 | 4 | pre | 0.203 | 0.071 | -0.106 | 0.000 | 0.232 | 0.000 |
| 2026-01-03 | 2 | 4 | pre | 0.430 | 0.077 | -0.189 | 0.500 | 0.604 | 0.000 |
| 2026-01-15 | 3 | 3 | pre | 0.630 | 0.080 | -0.208 | 0.000 | 0.144 | 0.000 |
| 2026-01-27 | 4 | 3 | pre | 0.589 | 0.065 | -0.121 | 0.000 | -0.199 | 0.000 |

| channel | points | min Δ | max Δ | max \|Δ\| | sign changes | distinct values |
| --- | --- | --- | --- | --- | --- | --- |
| gamma2 | 10 | -0.003 | 0.630 | 0.630 | 2 | 9 |
| P | 10 | 0.065 | 0.141 | 0.141 | 0 | 9 |
| D | 10 | -0.208 | 0.408 | 0.408 | 3 | 10 |
| F | 10 | -0.500 | 0.500 | 0.500 | 1 | 3 |
| S | 10 | -0.646 | 0.604 | 0.646 | 4 | 10 |
| A | 10 | -1.000 | 0.000 | 1.000 | 0 | 3 |

**The series moves where the pooled number cannot say so.** D swings between -0.2083 and 0.4083 (3 sign changes), S between -0.6457 and 0.6035 (4 sign changes), F between -0.5000 and 0.5000, and A down to -1.0000 — all of them channels whose pooled CI spans zero. A pooled zero is an average over orbit-dependent weeks, not the absence of an orbit effect; and with a single day of the record carrying two orbits, a per-day reading of any of it does not exist.

**No point at the post-event chip.** The record's single post-event echo is `2026-03-16` (ASC); its own 56-day window holds 2 ascending and 0 descending echo chips, so no descending median exists there and the series keeps no point: the omission is a number of this figure and not a silent absence.

The 9 steps of the record longer than the 6-day cadence (12, 24 days) are holes: the drawn line breaks at each of them and no value is interpolated across one.

| after | before | days | crosses the event |
| --- | --- | --- | --- |
| 2025-08-12 | 2025-08-24 | 12 | no |
| 2025-08-24 | 2025-09-05 | 12 | no |
| 2026-01-03 | 2026-01-15 | 12 | no |
| 2026-01-15 | 2026-01-27 | 12 | no |
| 2026-01-27 | 2026-02-08 | 12 | yes |
| 2026-02-08 | 2026-02-20 | 12 | no |
| 2026-02-20 | 2026-03-16 | 24 | no |
| 2026-03-16 | 2026-04-09 | 24 | no |
| 2026-04-15 | 2026-04-27 | 12 | no |

## Where the echoes sit over the two pass geometries

| orbit | state | chips | with echo | without echo | echo rate |
| --- | --- | --- | --- | --- | --- |
| ASC | pre | 16 | 9 | 7 | 0.562 |
| ASC | post | 4 | 1 | 3 | 0.250 |
| ASC | all | 20 | 10 | 10 | 0.500 |
| DESC | pre | 11 | 7 | 4 | 0.636 |
| DESC | post | 2 | 0 | 2 | 0.000 |
| DESC | all | 13 | 7 | 6 | 0.538 |

The coverage of the two geometries is a **distribution, not a contrast**: pooled Fisher p = 1.000 (10+10 against 7+6), and within an orbit pre-vs-post p = 0.582 (ASC) and 0.192 (DESC). Both orbits straddle the collapse date and the only post-event echo chip is an ascending pass: this block describes where the echoes sit over the two pass geometries, it is not an orbit contrast of the states.

Panel (h) draws the same block: the record's 32 acquisition days and 33 chips, 17 of them echoes (16 pre, 1 post), 20 chips on the ascending geometry and 13 on the descending one.

## What else moves with the pass geometry

What else moves with the pass geometry: the per-orbit medians of the export's co-variates. Ascending passes are the 09:2x UT ones, descending the 21:2x UT ones, so temperature and intensity move with the orbit by construction — they would have to be held fixed before an orbit split of the mask channels meant anything physical. The medians below are recomputed twice — from the channel table and from the committed measurement extract — and pinned against each other.

All 33 chips:

| co-variate | label | median ASC | median DESC | Δ (DESC − ASC) |
| --- | --- | --- | --- | --- |
| temperature_c | temperature (C) | 14.800000 | 10.600000 | -4.200000 |
| wind_speed_ms | wind (m/s) | 3.730000 | 2.380000 | -1.350000 |
| wind_gust_ms | gust (m/s) | 8.700000 | 5.000000 | -3.700000 |
| intensity | intensity | 6585.013262 | 1321.508762 | -5263.504501 |
| coherence_gamma2 | gamma2 (export) | 0.010965 | 0.017830 | 0.006865 |
| precipitation_mm | precip (mm) | 0.000000 | 0.000000 | 0.000000 |

The 17 echo-bearing chips:

| co-variate | label | median ASC | median DESC | Δ (DESC − ASC) |
| --- | --- | --- | --- | --- |
| temperature_c | temperature (C) | 14.750000 | 6.400000 | -8.350000 |
| wind_speed_ms | wind (m/s) | 3.580000 | 2.890000 | -0.690000 |
| wind_gust_ms | gust (m/s) | 9.600000 | 5.900000 | -3.700000 |
| intensity | intensity | 5233.029175 | 1301.816010 | -3931.213165 |
| coherence_gamma2 | gamma2 (export) | 0.010871 | 0.017830 | 0.006959 |
| precipitation_mm | precip (mm) | 0.000000 | 0.000000 | 0.000000 |

The ascending passes are the warmer ones (-4.2000 C on the record, -8.3500 C among the echo chips) and the darker ones in SAR intensity (-3931.2); no chip of the record carries rain (`precipitation_mm` is 0.0 on every chip). Whatever an orbit split of the mask channels would mean physically, these are the terms that move with it.

## The record: one orbit per day, and its nine holes

The record is 32 acquisition days and 57 committed windows (24 of them a repeat of a payload already in the table), which compile to 33 unique chips: 17 carry an echo mask (16 pre, 1 post) and 16 do not. 31 of the 32 days carry a single orbit; `2025-08-12` carries both, and both of its chips are echo-free, which is why the two orbits of a day never meet inside this figure's ordinate.

| day | state | orbits | windows | chips | repeats | with echo | without echo |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 2025-08-06 | pre | ASC | 1 | 1 | 0 | 1 | 0 |
| 2025-08-12 | pre | ASC, DESC | 2 | 2 | 0 | 0 | 2 |
| 2025-08-24 | pre | ASC | 1 | 1 | 0 | 0 | 1 |
| 2025-09-05 | pre | ASC | 1 | 1 | 0 | 1 | 0 |
| 2025-09-11 | pre | DESC | 1 | 1 | 0 | 1 | 0 |
| 2025-09-17 | pre | ASC | 1 | 1 | 0 | 1 | 0 |
| 2025-09-23 | pre | DESC | 1 | 1 | 0 | 0 | 1 |
| 2025-09-29 | pre | ASC | 1 | 1 | 0 | 0 | 1 |
| 2025-10-05 | pre | DESC | 2 | 1 | 1 | 0 | 1 |
| 2025-10-11 | pre | ASC | 2 | 1 | 1 | 1 | 0 |
| 2025-10-17 | pre | DESC | 2 | 1 | 1 | 1 | 0 |
| 2025-10-23 | pre | ASC | 2 | 1 | 1 | 1 | 0 |
| 2025-10-29 | pre | DESC | 2 | 1 | 1 | 1 | 0 |
| 2025-11-04 | pre | ASC | 2 | 1 | 1 | 0 | 1 |
| 2025-11-10 | pre | DESC | 2 | 1 | 1 | 1 | 0 |
| 2025-11-16 | pre | ASC | 2 | 1 | 1 | 0 | 1 |
| 2025-11-22 | pre | DESC | 2 | 1 | 1 | 0 | 1 |
| 2025-11-28 | pre | ASC | 2 | 1 | 1 | 0 | 1 |
| 2025-12-04 | pre | DESC | 2 | 1 | 1 | 1 | 0 |
| 2025-12-10 | pre | ASC | 2 | 1 | 1 | 1 | 0 |
| 2025-12-16 | pre | DESC | 2 | 1 | 1 | 1 | 0 |
| 2025-12-22 | pre | ASC | 2 | 1 | 1 | 0 | 1 |
| 2025-12-28 | pre | DESC | 2 | 1 | 1 | 1 | 0 |
| 2026-01-03 | pre | ASC | 2 | 1 | 1 | 1 | 0 |
| 2026-01-15 | pre | ASC | 2 | 1 | 1 | 1 | 0 |
| 2026-01-27 | pre | ASC | 2 | 1 | 1 | 1 | 0 |
| 2026-02-08 | post | ASC | 2 | 1 | 1 | 0 | 1 |
| 2026-02-20 | post | ASC | 2 | 1 | 1 | 0 | 1 |
| 2026-03-16 | post | ASC | 2 | 1 | 1 | 1 | 0 |
| 2026-04-09 | post | ASC | 2 | 1 | 1 | 0 | 1 |
| 2026-04-15 | post | DESC | 2 | 1 | 1 | 0 | 1 |
| 2026-04-27 | post | DESC | 2 | 1 | 1 | 0 | 1 |

## The pins and the provenance

1341 pin comparisons against five committed layers (the channel table's generated metadata, the mask-cache manifest, the committed mask cache, the committed paper result and the committed analysis reference, plus the raw measurement extract and the three resolved reference records): **1283 exact**, **58 within 1e-12** (the seeded numbers and the six-decimal record facts are pinned at 1e-06), **0 failed**. The recomputed number is in `recomputed`, the committed one in `committed`, and no artifact is written before they agree.

| committed input | layer | sha256 |
| --- | --- | --- |
| `data/ywf/ywf_ocv_channels.csv` | the channel table (this figure's rows) | `b5b17240246f6133…` |
| `data/ywf/ywf_windows_full.txt` | the window file (the chips) | `41b8838e80aebb21…` |
| `data/ywf/ywf_measurements_full.txt` | the raw measurement extract | `652a6c1d57766c31…` |
| `data/ywf/ywf_windows_mask_cache.txt` | the committed mask cache (the vector) | `0c29524d8569cd28…` |
| `data/ywf/ywf_windows_mask_cache.manifest.json` | the cache's manifest | `a14a32d23620c01b…` |
| `data/ywf/ywf_segments.txt` | the segment table | `5d48b2b1118d952e…` |
| `data/ywf/reference/ywf_bursts.json` | the resolved bursts (the pass) | `29e48c491585c4e7…` |
| `data/ywf/reference/ywf_geometry.json` | the geometry anchor | `57ac7a9815857df2…` |
| `data/ywf/reference/ywf_osm_road.json` | the road extract | `f88c4d6f2cd58c85…` |

Figures: `figures/ywf/fig_ywf_ocv_orbits.png` (600 dpi). JSON: `data/ywf/fig_ywf_ocv_orbits.json`.
