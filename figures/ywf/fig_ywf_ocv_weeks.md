# Yeongdeok Wind Farm (Samgye-ri) — the OCV time series in the weeks before the collapse

**Site** Yeongdeok Wind Farm (Samgye-ri) · **asset** Yeongdeok Wind Turbine (Samgye-ri) · **event** 2026-02-02 (Monday) · **window** 2025-08-06 -> 2026-01-27 (25.7 -> 0.9 weeks before the collapse) · **chips** 27 on 26 acquisition days, 16 with an echo

```
OCV  = [gamma2, P, D]              observability core vector
x    = -(2026-02-02 - day) / 7      weeks before the collapse (0 = event)
tick = the 1st of every month, its weeks before the collapse under it
mask = intensity >= 0.30*peak AND >= 5.0*np.median(intensity), A >= 2
step = 6 d (the two orbits alternate), 12 d = a gap: drawn, never interpolated over
post event chips: none in this figure — the window is the subject
```

## The window

| fact | value |
| --- | --- |
| acquisition days | **26**, 2025-08-06 -> 2026-01-27, a 6-day step with 4 gaps of 12 days |
| chips | **27** unique 7x7 chips (57 committed windows, 33 unique chips in the whole record) |
| echo coverage | **16 of 27** = 59.3 % (post-event, for reference: 1 of 6 = 16.7 %, Fisher exact p = 0.0854) |
| longest echo-free run | **3** consecutive chips without an echo |
| ISO weeks | **23** Monday-anchored buckets, 16 of them hold an echo |
| weather | temperature median 17.2 C (p5..p95 1.8..31.0 C), wind 3.13 m/s, gust 7.80 m/s |
| the pipeline `coherence` column | 0.998917748918 on all 33 rows (`n_unique` = 1) — a constant, so it cannot be a co-variate of the window |
| the `P` reference | the **17** echo-bearing chips of the whole record; `P` is a window-wide reference, not a per-chip quantity |

## The channels over the weeks

Dots are one chip each (no pooling); the rolling median is the trailing 4-week median of the echo-bearing chips with the gap honoured — a chip after a 12-day step sees the chips of the previous weeks, never an interpolated point.

| channel | n | median | p5 | p95 | rolling first | rolling last | rho vs. time | p | perm p |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gamma2 | 16 | 0.435 | 0.114 | 0.988 | 0.135 | 0.185 | 0.135 | 0.617 | 1.000 |
| P | 16 | 0.135 | 0.085 | 0.243 | 0.162 | 0.118 | -0.172 | 0.523 | 1.000 |
| D | 16 | 0.278 | 0.096 | 0.750 | 0.283 | 0.375 | -0.250 | 0.351 | 0.497 |
| A | 16 | 3.000 | 2.000 | 5.500 | - | - | 0.120 | 0.659 | 1.000 |

`A` (the echo of the window) is the coverage flag: 16 of 27 chips carry one (fig. a), and the three OCV channels are defined exactly there. The rolling medians of `gamma2`, `P` and `D` stay inside ±605 % of their first value over the whole window (figs. b-d), while the coverage moves with the halves of the window and the temperature falls (figs. a, f).

## The trend tests, and what they can carry

Spearman rho against **acquisition time forward** (days since 2025-08-06), with the moving-block permutation p of `ywf_ocv_stats.block_perm_p` beside it.

| series | n | rho | p | perm p |
| --- | --- | --- | --- | --- |
| gamma2 (echo chips) | 16 | 0.135 | 0.617 | 1.000 |
| P (echo chips) | 16 | -0.172 | 0.523 | 1.000 |
| D (echo chips) | 16 | -0.250 | 0.351 | 0.497 |
| A (echo chips) | 16 | 0.120 | 0.659 | 1.000 |
| echo flag (all chips) | 27 | 0.252 | 0.205 | 0.674 |
| temperature (C) (all chips) | 27 | -0.870 | 0.000 | 0.179 |
| wind (m/s) (all chips) | 27 | 0.387 | 0.046 | 0.346 |
| gust (m/s) (all chips) | 27 | 0.319 | 0.104 | 0.505 |
| intensity (all chips) | 27 | -0.181 | 0.365 | 0.332 |

**The permutation cannot resolve this window.** `block_perm_p` cuts blocks of `max(8, n // 10)` along acquisition order, so the 27-point series of the trend table is cut into 3 blocks of 8 chips and the null has only 6 distinct arrangements: its p can only take 6 values, and the number in that column is one of them. The Spearman p beside it is the informative one — and on this window what it reports is the *absence* of a channel trend.

## The coverage drift, the halves and the gaps

| half | days | chips | with echo | echo rate |
| --- | --- | --- | --- | --- |
| first half | 2025-08-06 -> 2025-10-29 | 14 | 8 | 0.571 |
| second half | 2025-11-04 -> 2026-01-27 | 13 | 8 | 0.615 |

The two halves are the window cut at its median acquisition index: the echo is not equally present along the axis, so a *pooled* pre-event median is a median of the chips the weather happened to make observable, and the figure says so instead of smoothing it.

| after | before | days | weeks before (after -> before) |
| --- | --- | --- | --- |
| 2025-08-12 | 2025-08-24 | 12 | 24.9 -> 23.1 |
| 2025-08-24 | 2025-09-05 | 12 | 23.1 -> 21.4 |
| 2026-01-03 | 2026-01-15 | 12 | 4.3 -> 2.6 |
| 2026-01-15 | 2026-01-27 | 12 | 2.6 -> 0.9 |

Every other step of the window is the 6-day repeat step of the two alternating orbits (ascending and descending passes); the four gaps above are drawn as bands in all six panels and are never interpolated across.

## The acquisition table

| day | W- | orbit | echo | mode | A | F | S | gamma2 | P | D | intensity | T (C) | wind | gust | status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2025-08-06 | 25.714 | ASCE | yes | intermediate | 3.000 | 2.000 | 1.667 | 0.165 | 0.196 | 0.600 | 8054.055 | 32.300 | 4.010 | 12.300 | processed |
| 2025-08-12 | 24.857 | ASCE | no | - | - | - | - | - | - | - | 10058.999 | 24.200 | 2.570 | 7.100 | processed |
| 2025-08-12 | 24.857 | DESC | no | - | - | - | - | - | - | - | 1409.322 | 22.400 | 0.520 | 1.200 | decorrelation_exceeded |
| 2025-08-24 | 23.143 | ASCE | no | - | - | - | - | - | - | - | 2463.409 | 33.200 | 3.820 | 9.000 | decorrelation_exceeded |
| 2025-09-05 | 21.429 | ASCE | yes | intermediate | 3.000 | 2.000 | 1.944 | 0.116 | 0.118 | 0.167 | 2792.128 | 28.100 | 4.450 | 11.200 | decorrelation_exceeded |
| 2025-09-11 | 20.571 | DESC | yes | compact | 2.000 | 2.000 | 2.000 | 0.154 | 0.206 | 0.400 | 1301.816 | 21.600 | 3.870 | 8.000 | decorrelation_exceeded |
| 2025-09-17 | 19.714 | ASCE | yes | intermediate | 5.000 | 2.000 | 0.200 | 0.401 | 0.094 | 0.333 | 2774.058 | 26.400 | 1.500 | 7.200 | decorrelation_exceeded |
| 2025-09-23 | 18.857 | DESC | no | - | - | - | - | - | - | - | 5190.896 | 20.900 | 1.910 | 3.500 | decorrelation_exceeded |
| 2025-09-29 | 18.000 | ASCE | no | - | - | - | - | - | - | - | 1870.342 | 22.000 | 3.640 | 8.400 | decorrelation_exceeded |
| 2025-10-05 | 17.143 | DESC | no | - | - | - | - | - | - | - | 1735.594 | 19.600 | 2.220 | 4.300 | decorrelation_exceeded |
| 2025-10-11 | 16.286 | ASCE | yes | compact | 2.000 | 2.000 | 1.803 | 0.986 | 0.088 | 0.167 | 9371.393 | 25.100 | 1.330 | 4.900 | processed |
| 2025-10-17 | 15.429 | DESC | yes | intermediate | 3.000 | 1.000 | 0.471 | 0.915 | 0.235 | 0.750 | 1430.639 | 18.300 | 0.430 | 1.200 | decorrelation_exceeded |
| 2025-10-23 | 14.571 | ASCE | yes | intermediate | 4.000 | 2.000 | 1.250 | 0.448 | 0.132 | 0.333 | 2736.837 | 16.500 | 2.070 | 5.000 | processed |
| 2025-10-29 | 13.714 | DESC | yes | intermediate | 4.000 | 2.000 | 1.275 | 0.421 | 0.265 | 0.222 | 1321.509 | 4.900 | 2.380 | 5.000 | decorrelation_exceeded |
| 2025-11-04 | 12.857 | ASCE | no | - | - | - | - | - | - | - | 2551.960 | 17.200 | 2.830 | 6.900 | decorrelation_exceeded |
| 2025-11-10 | 12.000 | DESC | yes | intermediate | 3.000 | 2.000 | 4.631 | 0.875 | 0.118 | 0.086 | 642.544 | 6.400 | 8.200 | 16.200 | decorrelation_exceeded |
| 2025-11-16 | 11.143 | ASCE | no | - | - | - | - | - | - | - | 9607.982 | 17.700 | 3.130 | 7.000 | decorrelation_exceeded |
| 2025-11-22 | 10.286 | DESC | no | - | - | - | - | - | - | - | 973.945 | 7.800 | 4.850 | 10.100 | decorrelation_exceeded |
| 2025-11-28 | 9.429 | ASCE | no | - | - | - | - | - | - | - | 7573.048 | 7.800 | 5.280 | 12.100 | decorrelation_exceeded |
| 2025-12-04 | 8.571 | DESC | yes | intermediate | 3.000 | 1.000 | 0.745 | 0.701 | 0.216 | 0.750 | 932.864 | 0.600 | 6.750 | 13.600 | decorrelation_exceeded |
| 2025-12-10 | 7.714 | ASCE | yes | compact | 2.000 | 2.000 | 2.062 | 0.723 | 0.088 | 0.200 | 12732.160 | 10.600 | 3.150 | 8.000 | decorrelation_exceeded |
| 2025-12-16 | 6.857 | DESC | yes | compact | 2.000 | 2.000 | 2.500 | 0.996 | 0.147 | 0.100 | 836.231 | 6.900 | 1.590 | 3.600 | decorrelation_exceeded |
| 2025-12-22 | 6.000 | ASCE | no | - | - | - | - | - | - | - | 1509.594 | 7.000 | 2.620 | 7.800 | decorrelation_exceeded |
| 2025-12-28 | 5.143 | DESC | yes | distributed | 7.000 | 3.000 | 1.519 | 0.815 | 0.168 | 0.167 | 2580.424 | 3.400 | 2.890 | 5.900 | processed |
| 2026-01-03 | 4.286 | ASCE | yes | intermediate | 4.000 | 1.000 | 0.750 | 0.107 | 0.074 | 0.444 | 1283.737 | 1.500 | 7.530 | 16.500 | processed |
| 2026-01-15 | 2.571 | ASCE | yes | intermediate | 3.000 | 2.000 | 1.374 | 0.185 | 0.118 | 0.375 | 8023.634 | 13.000 | 7.950 | 17.000 | decorrelation_exceeded |
| 2026-01-27 | 0.857 | ASCE | yes | intermediate | 3.000 | 2.000 | 3.145 | 0.266 | 0.137 | 0.150 | 2562.549 | 2.500 | 6.310 | 14.600 | decorrelation_exceeded |

`A`, `F`, `S`, `gamma2`, `P` and `D` are empty where the chip carries no echo: they are undefined there, not zero.

## The ISO-week buckets (23 weeks)

| week (Mon) | days | W- | chips | echo | A med | gamma2 med | P med | D med |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2025-08-04 | 08-06 | 25.714 | 1 | 1 | 3.000 | 0.165 | 0.196 | 0.600 |
| 2025-08-11 | 08-12 | 24.857 | 2 | 0 | - | - | - | - |
| 2025-08-18 | 08-24 | 23.143 | 1 | 0 | - | - | - | - |
| 2025-09-01 | 09-05 | 21.429 | 1 | 1 | 3.000 | 0.116 | 0.118 | 0.167 |
| 2025-09-08 | 09-11 | 20.571 | 1 | 1 | 2.000 | 0.154 | 0.206 | 0.400 |
| 2025-09-15 | 09-17 | 19.714 | 1 | 1 | 5.000 | 0.401 | 0.094 | 0.333 |
| 2025-09-22 | 09-23 | 18.857 | 1 | 0 | - | - | - | - |
| 2025-09-29 | 09-29, 10-05 | 18.000 | 2 | 0 | - | - | - | - |
| 2025-10-06 | 10-11 | 16.286 | 1 | 1 | 2.000 | 0.986 | 0.088 | 0.167 |
| 2025-10-13 | 10-17 | 15.429 | 1 | 1 | 3.000 | 0.915 | 0.235 | 0.750 |
| 2025-10-20 | 10-23 | 14.571 | 1 | 1 | 4.000 | 0.448 | 0.132 | 0.333 |
| 2025-10-27 | 10-29 | 13.714 | 1 | 1 | 4.000 | 0.421 | 0.265 | 0.222 |
| 2025-11-03 | 11-04 | 12.857 | 1 | 0 | - | - | - | - |
| 2025-11-10 | 11-10, 11-16 | 12.000 | 2 | 1 | 3.000 | 0.875 | 0.118 | 0.086 |
| 2025-11-17 | 11-22 | 10.286 | 1 | 0 | - | - | - | - |
| 2025-11-24 | 11-28 | 9.429 | 1 | 0 | - | - | - | - |
| 2025-12-01 | 12-04 | 8.571 | 1 | 1 | 3.000 | 0.701 | 0.216 | 0.750 |
| 2025-12-08 | 12-10 | 7.714 | 1 | 1 | 2.000 | 0.723 | 0.088 | 0.200 |
| 2025-12-15 | 12-16 | 6.857 | 1 | 1 | 2.000 | 0.996 | 0.147 | 0.100 |
| 2025-12-22 | 12-22, 12-28 | 6.000 | 2 | 1 | 7.000 | 0.815 | 0.168 | 0.167 |
| 2025-12-29 | 01-03 | 4.286 | 1 | 1 | 4.000 | 0.107 | 0.074 | 0.444 |
| 2026-01-12 | 01-15 | 2.571 | 1 | 1 | 3.000 | 0.185 | 0.118 | 0.375 |
| 2026-01-26 | 01-27 | 0.857 | 1 | 1 | 3.000 | 0.266 | 0.137 | 0.150 |

## Reading (a statement about coverage, not about damage)

Inside this window the echo is present on **16 of 27** chips and the longest run without one is **3** chips; the OCV values themselves show no monotone drift (gamma2 rho +0.135, p 0.62; P -0.172, p 0.52; D -0.250, p 0.35), while the temperature falls from late summer to deep winter (rho -0.870, p 0.0000) and the coverage moves between the halves. The pooled pre/post contrast of the paper figure rests on a window whose coordinate, coverage and weather are all moving: that is the reason this companion prints the time axis instead of only its median.

Consequence for the headline of [`fig_ywf_ocv_paper.md`](fig_ywf_ocv_paper.md): the 16/27 against 1/6 coverage table (Fisher p = 0.0854) must be read with the drift above attached — the pre side is not a stationary baseline, and the post side holds six chips. The figure claims neither a trend before the event nor a change at it.

## What this figure does not claim

* it draws **no** post-event chip: it is not a contrast, and it adds no number to the A-E pin contract of `fig_ywf_ocv_paper.py` (27 of the record's 33 chips are drawn, the 6 post-event ones are not);
* `A` here is the recomputed echo mask of the committed 7x7 window — **not** the pipeline's `coherence_masked_pixels`, which coincides with `A` on 0 of 17 echo-bearing rows;
* `P` is computed over the 17 echo-bearing chips of the whole record, so the chips of this window are not independent draws of it and its rolling median must be read as a shape statistic, not as a per-chip measurement;
* the moving-block permutation of this window cannot resolve p-values finer than the number of block arrangements, so no significance claim is made from it.

## Cross-checks

Each run re-derives the window from the committed channel table and compares it, field by field, with the committed mask cache (all 17 echo rows, the mask rule re-verified on each of them), the generated metadata, the paper result JSON and the analysis reference. **422 checks (421 exact, 1 within 1e-12), 0 deviations.**

| committed layer | what is pinned |
| --- | --- |
| `data/ywf/ywf_windows_mask_cache.txt` | the mask vector `A, D, F, S, gamma2` re-derived from the committed coordinates and `P` from the cache's own reference set, against the CSV — 17 rows x 6 channels, the mask rule re-verified row by row |
| `data/ywf/ywf_ocv_channels_meta.json` | the guard counts (`n_rows`, `n_unique_chips`, `n_echo_masks`, `n_P_set`, `n_pre_collapse`, `n_post_collapse`), the mask rule, the `coherence_masked_pixels` cross-check and the sha256 of all five committed inputs |
| `data/ywf/fig_ywf_ocv_paper.json` | the record counts, the coverage table and its Fisher test, the per-day timeline of the 26 pre-event days, the pre-side co-variate summary, the mask rule |
| `data/ywf/reference/ywf_ocv_findings.json` | the same coverage table, its power caveat and the per-day timeline — the two committed timelines must agree day by day |
| `data/ywf/reference/ywf_bursts.json` | the CDSE pass geometry of the record's chips (sub-swath, relative orbit, polarisation, azimuth time) and `ywf_osm_road.json` / `ywf_geometry.json` (the road the tower fell on, the estimated grid) — by the sha256 digests the reference records |

```bash
python3 code/ywf/ywf_ocv_masks.py --check              # the mask rule
python3 code/ywf/ywf_ocv_channels_csv.py             # build + verify
bash figures/ywf/fig_ywf_ocv_weeks.sh                 # this figure + pins
bash figures/ywf/fig_ywf_ocv_paper.sh                 # the pooled A-E figure
```
