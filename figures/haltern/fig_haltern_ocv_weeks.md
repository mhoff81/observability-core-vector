# Haltern AV9 Wind Park (Kreis Recklinghausen, NRW) — the OCV time series in the weeks before the collapse

**Site** Haltern AV9 Wind Park (Kreis Recklinghausen, NRW) · **asset** Haltern AV9 Wind Turbine (Nordex N149/4.5) · **event** 2021-09-29 (Wednesday) · **window** 2021-01-01 -> 2021-09-17 (38.7 -> 1.7 weeks before the collapse) · **chips** 108 on 108 acquisition days, 62 with an echo

```
OCV  = [gamma2, P, D]              observability core vector
x    = -(2021-09-29 - day) / 7      weeks before the collapse (0 = event)
tick = the 1st of every month, its weeks before the collapse under it
mask = intensity >= 0.30*peak AND >= 5.0*np.median(intensity), A >= 2
step = 6 d (the two orbits alternate), 12 d = a gap: drawn, never interpolated over
post event chips: none in this figure — the window is the subject
```

## The window

| fact | value |
| --- | --- |
| acquisition days | **108**, 2021-01-01 -> 2021-09-17, a 6-day step with 8 gaps of 12 days |
| chips | **108** unique 7x7 chips (168 committed windows, 120 unique chips in the whole record) |
| echo coverage | **62 of 108** = 57.4 % (post-event, for reference: 4 of 12 = 33.3 %, Fisher exact p = 0.1343) |
| longest echo-free run | **5** consecutive chips without an echo |
| ISO weeks | **29** Monday-anchored buckets, 26 of them hold an echo |
| weather | temperature median 11.4 C (p5..p95 -1.3..23.1 C), wind 3.06 m/s, gust 7.20 m/s |
| the pipeline `coherence` column | 0.999285714286 on all 120 rows (`n_unique` = 1) — a constant, so it cannot be a co-variate of the window |
| the `P` reference | the **66** echo-bearing chips of the whole record; `P` is a window-wide reference, not a per-chip quantity |

## The channels over the weeks

Dots are one chip each (no pooling); the rolling median is the trailing 4-week median of the echo-bearing chips with the gap honoured — a chip after a 12-day step sees the chips of the previous weeks, never an interpolated point.

| channel | n | median | p5 | p95 | rolling first | rolling last | rho vs. time | p | perm p |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gamma2 | 62 | 0.503 | 0.019 | 0.990 | 0.396 | 0.321 | -0.145 | 0.262 | 0.038 |
| P | 62 | 0.090 | 0.061 | 0.117 | 0.075 | 0.077 | 0.090 | 0.488 | 0.617 |
| D | 62 | 0.475 | 0.126 | 1.000 | 0.717 | 0.275 | -0.336 | 0.008 | 0.017 |
| A | 62 | 3.000 | 2.000 | 6.000 | - | - | 0.047 | 0.715 | 0.823 |

`A` (the echo of the window) is the coverage flag: 62 of 108 chips carry one (fig. a), and the three OCV channels are defined exactly there. The rolling medians of `gamma2`, `P` and `D` stay inside ±131 % of their first value over the whole window (figs. b-d), while the coverage moves with the halves of the window and the temperature falls (figs. a, f).

## The trend tests, and what they can carry

Spearman rho against **acquisition time forward** (days since 2025-08-06), with the moving-block permutation p of `haltern_ocv_stats.block_perm_p` beside it.

| series | n | rho | p | perm p |
| --- | --- | --- | --- | --- |
| gamma2 (echo chips) | 62 | -0.145 | 0.262 | 0.038 |
| P (echo chips) | 62 | 0.090 | 0.488 | 0.617 |
| D (echo chips) | 62 | -0.336 | 0.008 | 0.017 |
| A (echo chips) | 62 | 0.047 | 0.715 | 0.823 |
| echo flag (all chips) | 108 | 0.005 | 0.961 | 0.980 |
| temperature (C) (all chips) | 108 | 0.782 | 0.000 | 0.001 |
| wind (m/s) (all chips) | 108 | -0.177 | 0.068 | 0.134 |
| gust (m/s) (all chips) | 107 | -0.146 | 0.134 | 0.227 |
| intensity (all chips) | 108 | -0.005 | 0.963 | 0.963 |

**The permutation cannot resolve this window.** `block_perm_p` cuts blocks of `max(8, n // 10)` along acquisition order, so the 108-point series of the trend table is cut into 10 blocks of 10 chips and the null has only 3628800 distinct arrangements: its p can only take 3628800 values, and the number in that column is one of them. The Spearman p beside it is the informative one — and on this window what it reports is the *absence* of a channel trend.

## The coverage drift, the halves and the gaps

| half | days | chips | with echo | echo rate |
| --- | --- | --- | --- | --- |
| first half | 2021-01-01 -> 2021-05-08 | 54 | 30 | 0.556 |
| second half | 2021-05-10 -> 2021-09-17 | 54 | 32 | 0.593 |

The two halves are the window cut at its median acquisition index: the echo is not equally present along the axis, so a *pooled* pre-event median is a median of the chips the weather happened to make observable, and the figure says so instead of smoothing it.

| after | before | days | weeks before (after -> before) |
| --- | --- | --- | --- |
| 2021-01-17 | 2021-02-01 | 15 | 36.4 -> 34.3 |
| 2021-02-19 | 2021-03-02 | 11 | 31.7 -> 30.1 |
| 2021-03-18 | 2021-04-01 | 14 | 27.9 -> 25.9 |
| 2021-04-17 | 2021-05-01 | 14 | 23.6 -> 21.6 |
| 2021-05-17 | 2021-06-01 | 15 | 19.3 -> 17.1 |
| 2021-06-19 | 2021-07-01 | 12 | 14.6 -> 12.9 |
| 2021-07-19 | 2021-08-02 | 14 | 10.3 -> 8.3 |
| 2021-08-18 | 2021-09-01 | 14 | 6.0 -> 4.0 |

Every other step of the window is the 6-day repeat step of the two alternating orbits (ascending and descending passes); the four gaps above are drawn as bands in all six panels and are never interpolated across.

## The acquisition table

| day | W- | orbit | echo | mode | A | F | S | gamma2 | P | D | intensity | T (C) | wind | gust | status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2021-01-01 | 38.714 | ASCE | no | - | - | - | - | - | - | - | 2477.117 | 2.900 | 2.940 | 5.700 | acquired |
| 2021-01-02 | 38.571 | ASCE | no | - | - | - | - | - | - | - | 3302.563 | 2.500 | 0.810 | 2.000 | acquired |
| 2021-01-04 | 38.286 | DESC | no | - | - | - | - | - | - | - | 11555.021 | 1.400 | 4.470 | 9.200 | acquired |
| 2021-01-05 | 38.143 | DESC | yes | distributed | 6.000 | 2.000 | 2.062 | 0.783 | 0.086 | 0.600 | 1199.625 | 1.200 | 3.860 | 7.700 | acquired |
| 2021-01-07 | 37.857 | ASCE | yes | intermediate | 5.000 | 1.000 | 0.447 | 0.010 | 0.064 | 0.833 | 86361.853 | 1.600 | 3.090 | 6.200 | acquired |
| 2021-01-08 | 37.714 | ASCE | no | - | - | - | - | - | - | - | 4285.850 | 0.900 | 1.770 | 3.500 | acquired |
| 2021-01-10 | 37.429 | DESC | yes | intermediate | 4.000 | 2.000 | 0.791 | 0.946 | 0.091 | 0.333 | 20641.581 | -0.100 | 1.060 | 2.300 | acquired |
| 2021-01-11 | 37.286 | DESC | no | - | - | - | - | - | - | - | 4275.398 | 0.800 | 4.690 | 9.500 | acquired |
| 2021-01-13 | 37.000 | ASCE | no | - | - | - | - | - | - | - | 9849.153 | 3.200 | 5.920 | 12.200 | acquired |
| 2021-01-14 | 36.857 | ASCE | no | - | - | - | - | - | - | - | 9582.758 | 1.200 | 3.220 | 6.700 | acquired |
| 2021-01-16 | 36.571 | DESC | yes | distributed | 6.000 | 2.000 | 0.500 | 0.827 | 0.088 | 0.400 | 15395.108 | -0.800 | 1.840 | 3.300 | acquired |
| 2021-01-17 | 36.429 | DESC | yes | intermediate | 4.000 | 2.000 | 1.250 | 0.567 | 0.072 | 0.571 | 1047.538 | -0.600 | 2.570 | 5.100 | acquired |
| 2021-02-01 | 34.286 | ASCE | yes | intermediate | 3.000 | 2.000 | 1.374 | 0.305 | 0.076 | 0.500 | 5617.699 | 2.400 | 1.910 | 3.600 | acquired |
| 2021-02-03 | 34.000 | DESC | yes | intermediate | 5.000 | 1.000 | 0.566 | 0.978 | 0.106 | 0.625 | 23499.492 | 8.600 | 5.170 | 10.300 | acquired |
| 2021-02-04 | 33.857 | DESC | yes | distributed | 6.000 | 1.000 | 1.067 | 0.035 | 0.073 | 0.600 | 2807.959 | 6.500 | 5.570 | 11.900 | acquired |
| 2021-02-06 | 33.571 | ASCE | yes | intermediate | 3.000 | 3.000 | 2.749 | 0.006 | 0.091 | 0.250 | 7382.683 | 1.000 | 7.000 | 13.800 | acquired |
| 2021-02-07 | 33.429 | ASCE | yes | compact | 2.000 | 1.000 | 0.707 | 0.977 | 0.076 | 0.500 | 1887.853 | -5.900 | 6.360 | 15.600 | acquired |
| 2021-02-09 | 33.143 | DESC | no | - | - | - | - | - | - | - | 11918.950 | -9.900 | 1.660 | 3.600 | acquired |
| 2021-02-10 | 33.000 | DESC | yes | compact | 2.000 | 1.000 | 0.500 | 0.762 | 0.083 | 1.000 | 2237.968 | -11.400 | 2.020 | 3.800 | acquired |
| 2021-02-12 | 32.714 | ASCE | yes | intermediate | 3.000 | 2.000 | 1.054 | 0.035 | 0.035 | 0.375 | 1929.113 | -5.300 | 5.000 | 11.300 | acquired |
| 2021-02-15 | 32.286 | DESC | yes | compact | 2.000 | 1.000 | 0.500 | 0.992 | 0.091 | 1.000 | 17268.983 | -0.400 | 4.810 | 9.600 | acquired |
| 2021-02-16 | 32.143 | DESC | yes | compact | 2.000 | 1.000 | 0.500 | 0.735 | 0.083 | 1.000 | 1784.077 | 3.700 | 3.570 | 7.200 | acquired |
| 2021-02-18 | 31.857 | ASCE | no | - | - | - | - | - | - | - | 12535.173 | 9.700 | 6.200 | 14.200 | acquired |
| 2021-02-19 | 31.714 | ASCE | yes | compact | 2.000 | 1.000 | 0.500 | 0.829 | 0.068 | 1.000 | 4341.082 | 9.800 | 3.800 | 8.000 | acquired |
| 2021-03-02 | 30.143 | ASCE | yes | compact | 2.000 | 1.000 | 0.500 | 0.990 | 0.098 | 1.000 | 9993.991 | 14.000 | 1.360 | 4.700 | acquired |
| 2021-03-03 | 30.000 | ASCE | no | - | - | - | - | - | - | - | 15077.937 | 13.500 | 2.510 | 7.500 | acquired |
| 2021-03-05 | 29.714 | DESC | no | - | - | - | - | - | - | - | 3618.813 | -1.600 | 4.020 | 7.800 | acquired |
| 2021-03-06 | 29.571 | DESC | yes | compact | 2.000 | 1.000 | 0.500 | 0.999 | 0.106 | 1.000 | 2116.051 | -2.700 | 1.550 | 2.500 | acquired |
| 2021-03-08 | 29.286 | ASCE | no | - | - | - | - | - | - | - | 2276.091 | 6.400 | 1.210 | 4.100 | acquired |
| 2021-03-09 | 29.143 | ASCE | no | - | - | - | - | - | - | - | 4601.974 | 6.400 | 2.560 | 6.300 | acquired |
| 2021-03-11 | 28.857 | DESC | no | - | - | - | - | - | - | - | 7170.692 | 7.300 | 8.390 | 18.000 | acquired |
| 2021-03-12 | 28.714 | DESC | yes | intermediate | 4.000 | 2.000 | 0.559 | 0.311 | 0.095 | 0.400 | 4491.419 | 4.300 | 5.990 | 12.800 | acquired |
| 2021-03-14 | 28.429 | ASCE | yes | compact | 2.000 | 2.000 | 1.118 | 0.019 | 0.068 | 0.333 | 9518.566 | 7.000 | 5.590 | 11.500 | acquired |
| 2021-03-15 | 28.286 | ASCE | yes | compact | 2.000 | 1.000 | 0.500 | 0.951 | 0.114 | 1.000 | 11378.124 | 7.400 | 6.660 | 14.300 | acquired |
| 2021-03-17 | 28.000 | DESC | yes | intermediate | 3.000 | 1.000 | 0.333 | 0.914 | 0.061 | 0.500 | 12646.016 | 3.700 | 3.810 | 7.700 | acquired |
| 2021-03-18 | 27.857 | DESC | no | - | - | - | - | - | - | - | 4428.029 | 0.500 | 1.270 | 3.200 | acquired |
| 2021-04-01 | 25.857 | ASCE | yes | compact | 2.000 | 2.000 | 1.000 | 0.407 | 0.114 | 0.667 | 3054.305 | 12.300 | 5.010 | 11.100 | acquired |
| 2021-04-02 | 25.714 | ASCE | no | - | - | - | - | - | - | - | 10879.247 | 7.000 | 5.090 | 11.200 | acquired |
| 2021-04-04 | 25.429 | DESC | no | - | - | - | - | - | - | - | 23843.911 | 4.100 | 1.840 | 4.100 | acquired |
| 2021-04-05 | 25.286 | DESC | yes | distributed | 9.000 | 2.000 | 0.949 | 0.550 | 0.089 | 0.450 | 3009.615 | 5.900 | 6.140 | 12.900 | acquired |
| 2021-04-07 | 25.000 | ASCE | no | - | - | - | - | - | - | - | 3008.148 | 4.000 | 6.390 | 15.800 | acquired |
| 2021-04-08 | 24.857 | ASCE | yes | intermediate | 3.000 | 2.000 | 0.943 | 0.308 | 0.106 | 0.500 | 9804.077 | 7.100 | 4.890 | 10.900 | acquired |
| 2021-04-10 | 24.571 | DESC | no | - | - | - | - | - | - | - | 4325.909 | 4.900 | 3.550 | 6.900 | acquired |
| 2021-04-11 | 24.429 | DESC | yes | intermediate | 4.000 | 2.000 | 3.606 | 0.693 | 0.091 | 0.082 | 1306.535 | 4.100 | 3.450 | 7.100 | acquired |
| 2021-04-13 | 24.143 | ASCE | yes | compact | 2.000 | 1.000 | 0.707 | 0.169 | 0.106 | 0.500 | 5852.264 | 7.800 | 4.310 | 10.700 | acquired |
| 2021-04-14 | 24.000 | ASCE | no | - | - | - | - | - | - | - | 11470.345 | 6.300 | 4.900 | 10.000 | acquired |
| 2021-04-16 | 23.714 | DESC | yes | intermediate | 4.000 | 1.000 | 0.559 | 0.928 | 0.068 | 0.667 | 4445.832 | 2.200 | 1.770 | 4.400 | acquired |
| 2021-04-17 | 23.571 | DESC | yes | intermediate | 4.000 | 3.000 | 1.500 | 0.872 | 0.102 | 0.267 | 570.669 | 1.300 | 2.610 | 4.700 | acquired |
| 2021-05-01 | 21.571 | ASCE | yes | compact | 2.000 | 1.000 | 0.500 | 0.977 | 0.098 | 1.000 | 10141.718 | 11.600 | 3.040 | 7.700 | processed |
| 2021-05-02 | 21.429 | ASCE | no | - | - | - | - | - | - | - | 11434.476 | 10.300 | 4.830 | 11.100 | processed |
| 2021-05-04 | 21.143 | DESC | no | - | - | - | - | - | - | - | 4959.526 | 8.500 | 8.720 | 17.700 | processed |
| 2021-05-05 | 21.000 | DESC | no | - | - | - | - | - | - | - | 688.965 | 5.200 | 5.690 | 11.600 | processed |
| 2021-05-07 | 20.714 | ASCE | yes | intermediate | 4.000 | 2.000 | 0.791 | 0.218 | 0.057 | 0.143 | 6342.750 | 11.300 | 5.710 | 12.400 | processed |
| 2021-05-08 | 20.571 | ASCE | no | - | - | - | - | - | - | - | 3290.428 | 11.800 | 4.880 | 9.600 | processed |
| 2021-05-10 | 20.286 | DESC | yes | intermediate | 5.000 | 3.000 | 3.406 | 0.019 | 0.091 | 0.119 | 5829.522 | 19.400 | 7.100 | 14.700 | processed |
| 2021-05-11 | 20.143 | DESC | no | - | - | - | - | - | - | - | 547.394 | 12.800 | 2.340 | 4.300 | processed |
| 2021-05-13 | 19.857 | ASCE | yes | compact | 2.000 | 1.000 | 0.500 | 0.878 | 0.091 | 1.000 | 2107.936 | 16.400 | 1.680 | 7.600 | processed |
| 2021-05-14 | 19.714 | ASCE | yes | compact | 2.000 | 1.000 | 0.500 | 0.987 | 0.083 | 1.000 | 11016.013 | 15.800 | 2.660 | 6.900 | processed |
| 2021-05-16 | 19.429 | DESC | no | - | - | - | - | - | - | - | 3844.506 | 9.200 | 4.280 | 8.600 | processed |
| 2021-05-17 | 19.286 | DESC | no | - | - | - | - | - | - | - | 702.259 | 9.100 | 3.440 | 6.600 | processed |
| 2021-06-01 | 17.143 | ASCE | yes | intermediate | 4.000 | 2.000 | 1.118 | 0.965 | 0.117 | 0.333 | 14136.588 | 23.500 | 4.020 | 9.100 | processed |
| 2021-06-04 | 16.714 | DESC | yes | distributed | 10.000 | 1.000 | 0.412 | 0.262 | 0.079 | 0.625 | 2486.865 | 16.600 | 1.000 | 2.900 | processed |
| 2021-06-06 | 16.429 | ASCE | no | - | - | - | - | - | - | - | 8748.777 | 18.700 | 3.850 | 8.500 | processed |
| 2021-06-07 | 16.286 | ASCE | no | - | - | - | - | - | - | - | 12830.446 | 19.700 | 2.280 | 6.600 | processed |
| 2021-06-09 | 16.000 | DESC | yes | compact | 2.000 | 2.000 | 2.121 | 0.051 | 0.091 | 0.125 | 14295.311 | 16.700 | 1.040 | 2.400 | processed |
| 2021-06-10 | 15.857 | DESC | yes | intermediate | 3.000 | 2.000 | 1.374 | 0.259 | 0.086 | 0.300 | 2604.556 | 16.000 | 1.780 | 3.200 | processed |
| 2021-06-12 | 15.571 | ASCE | yes | intermediate | 4.000 | 1.000 | 0.250 | 0.313 | 0.087 | 0.667 | 2623.424 | 20.400 | 5.450 | 12.600 | processed |
| 2021-06-13 | 15.429 | ASCE | no | - | - | - | - | - | - | - | 19883.151 | 22.000 | 1.430 | 5.300 | processed |
| 2021-06-15 | 15.143 | DESC | yes | compact | 2.000 | 1.000 | 0.500 | 0.964 | 0.098 | 1.000 | 13116.104 | 15.600 | 2.620 | 5.900 | processed |
| 2021-06-16 | 15.000 | DESC | yes | intermediate | 5.000 | 3.000 | 2.828 | 0.181 | 0.076 | 0.357 | 4570.800 | 15.800 | 2.820 | 5.600 | processed |
| 2021-06-18 | 14.714 | ASCE | yes | distributed | 6.000 | 1.000 | 0.601 | 0.073 | 0.101 | 0.750 | 6790.035 | 30.000 | 2.960 | 9.500 | processed |
| 2021-06-19 | 14.571 | ASCE | yes | intermediate | 4.000 | 2.000 | 1.768 | 0.003 | 0.083 | 0.267 | 5270.608 | 26.400 | 4.000 | 9.400 | processed |
| 2021-07-01 | 12.857 | ASCE | yes | compact | 2.000 | 2.000 | 1.118 | 0.662 | 0.129 | 0.333 | 4028.717 | 14.700 | 4.920 | 10.800 | processed |
| 2021-07-03 | 12.571 | DESC | no | - | - | - | - | - | - | - | 12049.490 | 13.300 | 1.210 | 3.000 | processed |
| 2021-07-06 | 12.143 | ASCE | no | - | - | - | - | - | - | - | 12465.094 | 20.800 | 7.290 | 16.400 | processed |
| 2021-07-07 | 12.000 | ASCE | yes | compact | 2.000 | 2.000 | 2.121 | 0.131 | 0.136 | 0.125 | 2437.162 | 22.300 | 2.970 | 7.200 | processed |
| 2021-07-09 | 11.714 | DESC | yes | intermediate | 4.000 | 1.000 | 0.707 | 0.961 | 0.080 | 1.000 | 4471.375 | 15.800 | 2.200 | 4.400 | processed |
| 2021-07-10 | 11.571 | DESC | no | - | - | - | - | - | - | - | 12749.692 | 15.100 | 1.940 | 3.700 | processed |
| 2021-07-12 | 11.286 | ASCE | yes | intermediate | 3.000 | 2.000 | 1.491 | 0.468 | 0.126 | 0.300 | 6199.669 | 23.200 | 1.660 | 5.200 | processed |
| 2021-07-13 | 11.143 | ASCE | yes | intermediate | 4.000 | 1.000 | 0.500 | 0.895 | 0.080 | 1.000 | 14274.946 | 20.100 | 2.730 | 6.000 | processed |
| 2021-07-15 | 10.857 | DESC | no | - | - | - | - | - | - | - | 4456.678 | 17.100 | 0.510 | 1.600 | processed |
| 2021-07-16 | 10.714 | DESC | yes | distributed | 6.000 | 2.000 | 1.740 | 0.557 | 0.109 | 0.240 | 1348.262 | 15.300 | 2.590 | 6.400 | processed |
| 2021-07-18 | 10.429 | ASCE | yes | intermediate | 5.000 | 2.000 | 1.442 | 0.084 | 0.097 | 0.333 | 2519.650 | 24.800 | 4.700 | 10.400 | processed |
| 2021-07-19 | 10.286 | ASCE | yes | compact | 2.000 | 2.000 | 3.041 | 0.299 | 0.106 | 0.143 | 4563.544 | 20.000 | 3.700 | 8.300 | processed |
| 2021-08-02 | 8.286 | DESC | no | - | - | - | - | - | - | - | 2415.144 | 13.300 | 2.010 | 3.800 | processed |
| 2021-08-03 | 8.143 | DESC | yes | intermediate | 3.000 | 2.000 | 1.202 | 0.286 | 0.066 | 0.333 | 2509.921 | 12.900 | 1.300 | 2.300 | processed |
| 2021-08-05 | 7.857 | ASCE | no | - | - | - | - | - | - | - | 12705.219 | 21.600 | 2.020 | 5.300 | processed |
| 2021-08-06 | 7.714 | ASCE | yes | intermediate | 4.000 | 2.000 | 1.458 | 0.843 | 0.098 | 0.333 | 3820.977 | 19.200 | 5.000 | 13.600 | processed |
| 2021-08-08 | 7.429 | DESC | yes | compact | 2.000 | 2.000 | 1.118 | 0.499 | 0.098 | 0.333 | 2958.166 | 14.300 | 6.150 | 13.300 | processed |
| 2021-08-09 | 7.286 | DESC | no | - | - | - | - | - | - | - | 1489.297 | 13.900 | 4.030 | 7.900 | processed |
| 2021-08-11 | 7.000 | ASCE | no | - | - | - | - | - | - | - | 6884.513 | 21.900 | 0.540 | 3.200 | processed |
| 2021-08-12 | 6.857 | ASCE | no | - | - | - | - | - | - | - | 15266.726 | 24.500 | 1.430 | 4.900 | processed |
| 2021-08-14 | 6.571 | DESC | yes | intermediate | 4.000 | 3.000 | 1.521 | 0.337 | 0.080 | 0.160 | 3748.982 | 14.000 | 2.270 | 4.000 | processed |
| 2021-08-15 | 6.429 | DESC | yes | compact | 2.000 | 1.000 | 0.500 | 0.477 | 0.098 | 1.000 | 7052.828 | 15.100 | 2.690 | 4.800 | processed |
| 2021-08-17 | 6.143 | ASCE | no | - | - | - | - | - | - | - | 13395.522 | 14.000 | 4.490 | 10.300 | processed |
| 2021-08-18 | 6.000 | ASCE | yes | intermediate | 3.000 | 1.000 | 0.333 | 0.938 | 0.051 | 0.500 | 14488.888 | 18.800 | 5.380 | 12.600 | processed |
| 2021-09-01 | 4.000 | DESC | yes | compact | 2.000 | 1.000 | 0.500 | 0.990 | 0.061 | 1.000 | 13661.511 | 14.000 | 1.730 | 3.200 | processed |
| 2021-09-02 | 3.857 | DESC | no | - | - | - | - | - | - | - | 2569.575 | 10.000 | 1.330 | 1.500 | processed |
| 2021-09-04 | 3.571 | ASCE | yes | distributed | 6.000 | 2.000 | 0.707 | 0.042 | 0.068 | 0.300 | 3423.416 | 20.200 | 2.500 | 6.200 | processed |
| 2021-09-05 | 3.429 | ASCE | no | - | - | - | - | - | - | - | 3197.007 | 23.000 | 1.750 | 4.700 | processed |
| 2021-09-07 | 3.143 | DESC | no | - | - | - | - | - | - | - | 4402.313 | 14.300 | 1.300 | 2.700 | processed |
| 2021-09-08 | 3.000 | DESC | no | - | - | - | - | - | - | - | 9251.444 | 10.000 | 5.000 | - | processed |
| 2021-09-10 | 2.714 | ASCE | no | - | - | - | - | - | - | - | 9188.065 | 19.600 | 2.860 | 7.900 | processed |
| 2021-09-11 | 2.571 | ASCE | no | - | - | - | - | - | - | - | 3223.153 | 19.700 | 4.310 | 10.100 | processed |
| 2021-09-13 | 2.286 | DESC | yes | distributed | 6.000 | 4.000 | 2.224 | 0.208 | 0.098 | 0.250 | 3773.146 | 10.000 | 2.020 | 2.800 | processed |
| 2021-09-14 | 2.143 | DESC | yes | intermediate | 3.000 | 3.000 | 2.749 | 0.137 | 0.091 | 0.167 | 9685.211 | 11.600 | 2.650 | 5.000 | processed |
| 2021-09-16 | 1.857 | ASCE | yes | intermediate | 5.000 | 3.000 | 1.456 | 0.435 | 0.082 | 0.167 | 1968.878 | 18.800 | 4.260 | 9.800 | processed |
| 2021-09-17 | 1.714 | ASCE | yes | intermediate | 4.000 | 2.000 | 0.791 | 0.507 | 0.072 | 0.333 | 16157.073 | 19.700 | 2.330 | 6.000 | processed |

`A`, `F`, `S`, `gamma2`, `P` and `D` are empty where the chip carries no echo: they are undefined there, not zero.

## The ISO-week buckets (29 weeks)

| week (Mon) | days | W- | chips | echo | A med | gamma2 med | P med | D med |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2020-12-28 | 01-01, 01-02 | 38.714 | 2 | 0 | - | - | - | - |
| 2021-01-04 | 01-04, 01-05, 01-07, 01-08, 01-10 | 38.286 | 5 | 3 | 5.000 | 0.783 | 0.086 | 0.600 |
| 2021-01-11 | 01-11, 01-13, 01-14, 01-16, 01-17 | 37.286 | 5 | 2 | 5.000 | 0.697 | 0.080 | 0.486 |
| 2021-02-01 | 02-01, 02-03, 02-04, 02-06, 02-07 | 34.286 | 5 | 5 | 3.000 | 0.305 | 0.076 | 0.500 |
| 2021-02-08 | 02-09, 02-10, 02-12 | 33.143 | 3 | 2 | 2.500 | 0.398 | 0.059 | 0.688 |
| 2021-02-15 | 02-15, 02-16, 02-18, 02-19 | 32.286 | 4 | 3 | 2.000 | 0.829 | 0.083 | 1.000 |
| 2021-03-01 | 03-02, 03-03, 03-05, 03-06 | 30.143 | 4 | 2 | 2.000 | 0.995 | 0.102 | 1.000 |
| 2021-03-08 | 03-08, 03-09, 03-11, 03-12, 03-14 | 29.286 | 5 | 2 | 3.000 | 0.165 | 0.081 | 0.367 |
| 2021-03-15 | 03-15, 03-17, 03-18 | 28.286 | 3 | 2 | 2.500 | 0.932 | 0.087 | 0.750 |
| 2021-03-29 | 04-01, 04-02, 04-04 | 25.857 | 3 | 1 | 2.000 | 0.407 | 0.114 | 0.667 |
| 2021-04-05 | 04-05, 04-07, 04-08, 04-10, 04-11 | 25.286 | 5 | 3 | 4.000 | 0.550 | 0.091 | 0.450 |
| 2021-04-12 | 04-13, 04-14, 04-16, 04-17 | 24.143 | 4 | 3 | 4.000 | 0.872 | 0.102 | 0.500 |
| 2021-04-26 | 05-01, 05-02 | 21.571 | 2 | 1 | 2.000 | 0.977 | 0.098 | 1.000 |
| 2021-05-03 | 05-04, 05-05, 05-07, 05-08 | 21.143 | 4 | 1 | 4.000 | 0.218 | 0.057 | 0.143 |
| 2021-05-10 | 05-10, 05-11, 05-13, 05-14, 05-16 | 20.286 | 5 | 3 | 2.000 | 0.878 | 0.091 | 1.000 |
| 2021-05-17 | 05-17 | 19.286 | 1 | 0 | - | - | - | - |
| 2021-05-31 | 06-01, 06-04, 06-06 | 17.143 | 3 | 2 | 7.000 | 0.613 | 0.098 | 0.479 |
| 2021-06-07 | 06-07, 06-09, 06-10, 06-12, 06-13 | 16.286 | 5 | 3 | 3.000 | 0.259 | 0.087 | 0.300 |
| 2021-06-14 | 06-15, 06-16, 06-18, 06-19 | 15.143 | 4 | 4 | 4.500 | 0.127 | 0.091 | 0.554 |
| 2021-06-28 | 07-01, 07-03 | 12.857 | 2 | 1 | 2.000 | 0.662 | 0.129 | 0.333 |
| 2021-07-05 | 07-06, 07-07, 07-09, 07-10 | 12.143 | 4 | 2 | 3.000 | 0.546 | 0.108 | 0.562 |
| 2021-07-12 | 07-12, 07-13, 07-15, 07-16, 07-18 | 11.286 | 5 | 4 | 4.500 | 0.513 | 0.103 | 0.317 |
| 2021-07-19 | 07-19 | 10.286 | 1 | 1 | 2.000 | 0.299 | 0.106 | 0.143 |
| 2021-08-02 | 08-02, 08-03, 08-05, 08-06, 08-08 | 8.286 | 5 | 3 | 3.000 | 0.499 | 0.098 | 0.333 |
| 2021-08-09 | 08-09, 08-11, 08-12, 08-14, 08-15 | 7.286 | 5 | 2 | 3.000 | 0.407 | 0.089 | 0.580 |
| 2021-08-16 | 08-17, 08-18 | 6.143 | 2 | 1 | 3.000 | 0.938 | 0.051 | 0.500 |
| 2021-08-30 | 09-01, 09-02, 09-04, 09-05 | 4.000 | 4 | 2 | 4.000 | 0.516 | 0.064 | 0.650 |
| 2021-09-06 | 09-07, 09-08, 09-10, 09-11 | 3.143 | 4 | 0 | - | - | - | - |
| 2021-09-13 | 09-13, 09-14, 09-16, 09-17 | 2.286 | 4 | 4 | 4.500 | 0.321 | 0.086 | 0.208 |

## Reading (a statement about coverage, not about damage)

Inside this window the echo is present on **62 of 108** chips and the longest run without one is **5** chips; the OCV values themselves show no monotone drift (gamma2 rho -0.145, p 0.26; P +0.090, p 0.49; D -0.336, p 0.01), while the temperature falls from late summer to deep winter (rho +0.782, p 0.0000) and the coverage moves between the halves. The pooled pre/post contrast of the paper figure rests on a window whose coordinate, coverage and weather are all moving: that is the reason this companion prints the time axis instead of only its median.

Consequence for the headline of [`fig_haltern_ocv_paper.md`](fig_haltern_ocv_paper.md): the 62/108 against 4/12 coverage table (Fisher p = 0.1343) must be read with the drift above attached — the pre side is not a stationary baseline, and the post side holds six chips. The figure claims neither a trend before the event nor a change at it.

## What this figure does not claim

* it draws **no** post-event chip: it is not a contrast, and it adds no number to the A-E pin contract of `fig_haltern_ocv_paper.py` (108 of the record's 120 chips are drawn, the 12 post-event ones are not);
* `A` here is the recomputed echo mask of the committed 7x7 window — **not** the pipeline's `coherence_masked_pixels`, which coincides with `A` on 6 of 66 echo-bearing rows;
* `P` is computed over the 66 echo-bearing chips of the whole record, so the chips of this window are not independent draws of it and its rolling median must be read as a shape statistic, not as a per-chip measurement;
* the moving-block permutation of this window cannot resolve p-values finer than the number of block arrangements, so no significance claim is made from it.

## Cross-checks

Each run re-derives the window from the committed channel table and compares it, field by field, with the committed mask cache (all 66 echo rows, the mask rule re-verified on each of them), the generated metadata, the paper result JSON and the analysis reference. **1097 checks (1081 exact, 16 within 1e-12), 16 deviations.**

| committed layer | what is pinned |
| --- | --- |
| `data/haltern/haltern_windows_mask_cache.txt` | the mask vector `A, D, F, S, gamma2` re-derived from the committed coordinates and `P` from the cache's own reference set, against the CSV — 66 rows x 6 channels, the mask rule re-verified row by row |
| `data/haltern/haltern_ocv_channels_meta.json` | the guard counts (`n_rows`, `n_unique_chips`, `n_echo_masks`, `n_P_set`, `n_pre_collapse`, `n_post_collapse`), the mask rule, the `coherence_masked_pixels` cross-check and the sha256 of all five committed inputs |
| `data/haltern/fig_haltern_ocv_paper.json` | the record counts, the coverage table and its Fisher test, the per-day timeline of the 108 pre-event days, the pre-side co-variate summary, the mask rule |
| `data/haltern/reference/haltern_ocv_findings.json` | the same coverage table, its power caveat and the per-day timeline — the two committed timelines must agree day by day |
| `data/haltern/reference/haltern_bursts.json` | the CDSE pass geometry of the record's chips (sub-swath, relative orbit, polarisation, azimuth time) and `haltern_osm_node.json` / `haltern_geometry.json` (the turbine's own OSM node, the anchor derived from it, the estimated grid) — by the sha256 digests the reference records |

| failed check | recomputed | committed |
| --- | --- | --- |
| `window.event_weekday` | `'Wednesday'` | `'Monday'` |
| `gaps[0].days` | `15` | `12` |
| `gaps[1].days` | `11` | `12` |
| `gaps[2].days` | `14` | `12` |
| `gaps[3].days` | `14` | `12` |
| `gaps[4].days` | `15` | `12` |
| `gaps[6].days` | `14` | `12` |
| `gaps[7].days` | `14` | `12` |
| `gaps.steps_are_cadence_or_gap` | `'<len 4>'` | `'<len 1>'` |
| `weeks.n` | `29` | `23` |
| `weeks.echo_per_day_max` | `5` | `1` |
| `echo_free_run` | `5` | `3` |
| `weeks.first_wb` | `38.714` | `25.714` |
| `weeks.last_wb` | `1.714` | `0.857` |
| `months.labels` | `'<len 9>'` | `'<len 6>'` |
| `coherence.constant` | `0.999285714` | `0.998917749` |

```bash
python3 code/haltern/haltern_ocv_masks.py --check              # the mask rule
python3 code/haltern/haltern_ocv_channels_csv.py             # build + verify
bash figures/haltern/fig_haltern_ocv_weeks.sh                 # this figure + pins
bash figures/haltern/fig_haltern_ocv_paper.sh                 # the pooled A-E figure
```
