# Yeongdeok Wind Farm (Samgye-ri) — the echo-mask vector over the whole record

**Site** Yeongdeok Wind Farm (Samgye-ri) · **asset** Yeongdeok Wind Turbine (Samgye-ri) · **event** 2026-02-02 (Monday, inside the axis, inside a step) · **record** 2025-08-06 -> 2026-04-27 (265 days) · **chips** 33 unique 7x7 chips on 32 acquisition days, 17 carrying an echo

```
OCV  = [gamma2, P, D]               observability core vector
x    = the acquisition day, one dot per chip (2025-08-06 -> 2026-04-27)
y    = the value over that dimension's own pre-window median, taken
       over the 16 pre-event echo chips; 1.0 = that median
mask = intensity >= 0.30*peak AND >= 5.0*np.median(intensity), A >= 2 (window 7x7)
step = 6 d apart (the two orbits alternate); 9 of 31 steps are longer and
       are drawn as bands — a hole is never interpolated over
post = 1 of 6 post-event chips carry an echo, so the post side carries
       one point and no line, no median and no trend
```

## The record

| fact | value |
| --- | --- |
| chips | **33** unique 7x7 chips behind **57** committed windows (24 byte-identical duplicates of the windows file) |
| acquisition days | **32**, 2025-08-06 -> 2026-04-27 (265 days): a 6-day cadence with 9 longer steps and 22 cadence steps |
| pre-event | **27 chips on 26 days**, 2025-08-06 -> 2026-01-27 (175 d), 16 with an echo (59.3 %) |
| post-event | **6 chips on 6 days**, 2026-02-08 -> 2026-04-27 (79 d), 1 with an echo (16.7 %) |
| the event falls inside a step | 2026-02-02 sits in the **12-day** step 2026-01-27 -> 2026-02-08; no chip on the event day, so the two stretches meet at a band and not at a point |
| the only two-orbit day | 2025-08-12 carries both passes (Ascending, Descending), and **neither** echoes |
| the `P` reference | `P` is a window-wide reference, not a per-chip quantity: the committed cache carries `P_reference_chips` = 17 — the same 17 echo-bearing chips of the whole record of which the mean pixel frequency inside the own mask is taken |
| the pipeline `coherence_masked_pixels` column | it runs 5.000 .. 167.000 masked pixels over the 27 pre-event rows (median 75.000) and agrees with the recomputed `A` on **0 of 17** echo chips — see the caveat below |
| the constant control `coherence` | 0.998917748918 on every one of the 33 chips (1 distinct value): reported as a control, never drawn |

## The six panels

Every panel draws the same axis: one dot per echo-bearing chip, y = that dimension's value over **its own** pre-window median. The committed quartiles of the pre window are divided by the same median to make the box, so every drawn ordinate is a pinned number.

| panel | dimension | y scale | base (pre median, raw) | pre n | pre min (ratio) | pre p25 (ratio) | pre p75 (ratio) | pre max (ratio) | post n | the one post chip (ratio) | post median (raw) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| a | gamma2 | log | 0.435 | 16 | 0.247 | 0.415 | 1.910 | 2.292 | 1 | 0.434 | 0.189 |
| b | A | log | 3.000 | 16 | 0.667 | 0.917 | 1.333 | 2.333 | 1 | 1.667 | 5.000 |
| c | P | linear | 0.135 | 16 | 0.545 | 0.829 | 1.473 | 1.964 | 1 | 1.222 | 0.165 |
| d | D | log | 0.278 | 16 | 0.309 | 0.600 | 1.480 | 2.700 | 1 | 0.429 | 0.119 |
| e | F | log | 2.000 | 16 | 0.500 | 1.000 | 1.000 | 1.500 | 1 | 1.500 | 3.000 |
| f | S | log | 1.593 | 16 | 0.126 | 0.706 | 1.265 | 2.908 | 1 | 0.888 | 1.414 |

`P` is the only linear panel: it is the mean pixel frequency inside the mask over the 17 reference chips — a bounded frequency, so its ratio is read as a difference. The other five are multiplicative magnitudes (49 is the window, the denominator of `A` and `D`); all six pre windows stay inside 0.12..2.91 of 1.0, so no log axis hides a level. The centre line of every box is 1.0 by construction: it *is* the pre-window median the ordinate is divided by.

## The steps of the record

Acquisitions come in pairs: 9 of the 31 steps between two consecutive acquisition days are longer than the 6-day cadence. Each of them is drawn as a grey band and **no line crosses it** — a line from one chip to the next would be an interpolation over a hole this record does not have data for.

| step | how many | what it is | the steps |
| --- | --- | --- | --- |
| 6 d | 22 | the cadence | 2025-08-06 -> 2025-08-12, 2025-09-05 -> 2025-09-11, 2025-09-11 -> 2025-09-17, 2025-09-17 -> 2025-09-23, 2025-09-23 -> 2025-09-29, 2025-09-29 -> 2025-10-05, 2025-10-05 -> 2025-10-11, 2025-10-11 -> 2025-10-17, 2025-10-17 -> 2025-10-23, 2025-10-23 -> 2025-10-29, 2025-10-29 -> 2025-11-04, 2025-11-04 -> 2025-11-10, 2025-11-10 -> 2025-11-16, 2025-11-16 -> 2025-11-22, 2025-11-22 -> 2025-11-28, 2025-11-28 -> 2025-12-04, 2025-12-04 -> 2025-12-10, 2025-12-10 -> 2025-12-16, 2025-12-16 -> 2025-12-22, 2025-12-22 -> 2025-12-28, 2025-12-28 -> 2026-01-03, 2026-04-09 -> 2026-04-15 |
| 12 d | 7 | a hole, drawn as a band | 2025-08-12 -> 2025-08-24, 2025-08-24 -> 2025-09-05, 2026-01-03 -> 2026-01-15, 2026-01-15 -> 2026-01-27, 2026-01-27 -> 2026-02-08, 2026-02-08 -> 2026-02-20, 2026-04-15 -> 2026-04-27 |
| 24 d | 2 | a hole, drawn as a band | 2026-02-20 -> 2026-03-16, 2026-03-16 -> 2026-04-09 |

## The two stretches

| stretch | first day | last day | calendar days | chips | acquisition days | committed windows | duplicates | with echo | without echo | echo rate | orbits | steps | longer steps | step lengths |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| pre-collapse (healthy) | 2025-08-06 | 2026-01-27 | 175 | 27 | 26 | 45 | 18 | 16 | 11 | 0.593 | Ascending, Descending | 25 | 4 | 6 d, 12 d |
| post-collapse | 2026-02-08 | 2026-04-27 | 79 | 6 | 6 | 12 | 6 | 1 | 5 | 0.167 | Ascending, Descending | 5 | 4 | 6 d, 12 d, 24 d |

The two stretches are unequal in every respect: **175 days and 27 chips** before the collapse against **79 days and 6 chips** after it. The collapse step (2026-01-27 -> 2026-02-08) is the only step of the record that holds the event date, and the dashed cut of every panel falls inside it: the two stretches of the figure meet at a grey band, never at a point, and the nearest chips sit 6 days either side of the event date.

## The coverage: what the collapse left behind

| state | echo | no echo |
| --- | --- | --- |
| pre | 16 | 11 |
| post | 1 | 5 |

Fisher's exact test on that 2x2 table: **p = 0.0854** (two-sided, alpha = 0.05). With 16 of 27 pre-event chips carrying an echo and 1 of 6 post-event chips doing so, the record does not reject "the coverage is unchanged" — and it cannot, with 6 post chips.

**The power question of the record.** If all post-event chips came back echo-free, it would take **4 post-event chips** for the same test to fall below alpha = 0.05 against the pre-event rate of 59.3 %. The record has 6 chips there, one of which echoes.

1 of 6 post-event chips carries an echo, so there is no post-state distribution: the collapse contrast of this site is an echo-coverage statement with this caveat attached, not a damage verdict.

## Where the echoes sit over the two pass geometries

| orbit | chips | with echo | without echo | echo rate | pre chips | pre echo rate | post chips | post echo rate | Fisher within the orbit |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ASCENDING | 20 | 10 | 10 | 0.500 | 16 | 0.562 | 4 | 0.250 | 0.582 |
| DESCENDING | 13 | 7 | 6 | 0.538 | 11 | 0.636 | 2 | 0.000 | 0.192 |

Across the whole record the two geometries are indistinguishable on coverage (Fisher p = 1.000). Both orbits straddle the collapse date and the only post-event echo chip is an ascending pass: this block describes where the echoes sit over the two pass geometries, it is not an orbit contrast of the states.

## The composition of the record

| month | pre chips | post chips |
| --- | --- | --- |
| Jan | 3 | 0 |
| Feb | 0 | 2 |
| Mar | 0 | 1 |
| Apr | 0 | 3 |
| Aug | 4 | 0 |
| Sep | 5 | 0 |
| Oct | 5 | 0 |
| Nov | 5 | 0 |
| Dec | 5 | 0 |

27 of the 33 chips are pre-event and 6 post-event. The asset carries 1 labelled section (`turbine-mast-section-3`) — a single turbine-mast section, so no panel of this figure is a comparison across sections.

## The echo geometries

| state | echo chips | compact | distributed | intermediate |
| --- | --- | --- | --- | --- |
| pre | 16 | 4 | 1 | 11 |
| post | 1 | 0 | 0 | 1 |

## The record, day by day

Every acquisition day of the record, with what was requested, what the mask rule found and which pass geometry carried it. This is the committed timeline (`n_windows`, `n_unique_chips`, `n_duplicated`, `n_with_echo`, `state`) as the channel table recomputes it.

| day | state | committed windows | unique chips | duplicates | with echo | orbits | echo geometry |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 2025-08-06 | pre | 1 | 1 | 0 | 1 | ASCE | intermediate |
| 2025-08-12 | pre | 2 | 2 | 0 | 0 | ASCE, DESC | - |
| 2025-08-24 | pre | 1 | 1 | 0 | 0 | ASCE | - |
| 2025-09-05 | pre | 1 | 1 | 0 | 1 | ASCE | intermediate |
| 2025-09-11 | pre | 1 | 1 | 0 | 1 | DESC | compact |
| 2025-09-17 | pre | 1 | 1 | 0 | 1 | ASCE | intermediate |
| 2025-09-23 | pre | 1 | 1 | 0 | 0 | DESC | - |
| 2025-09-29 | pre | 1 | 1 | 0 | 0 | ASCE | - |
| 2025-10-05 | pre | 2 | 1 | 1 | 0 | DESC | - |
| 2025-10-11 | pre | 2 | 1 | 1 | 1 | ASCE | compact |
| 2025-10-17 | pre | 2 | 1 | 1 | 1 | DESC | intermediate |
| 2025-10-23 | pre | 2 | 1 | 1 | 1 | ASCE | intermediate |
| 2025-10-29 | pre | 2 | 1 | 1 | 1 | DESC | intermediate |
| 2025-11-04 | pre | 2 | 1 | 1 | 0 | ASCE | - |
| 2025-11-10 | pre | 2 | 1 | 1 | 1 | DESC | intermediate |
| 2025-11-16 | pre | 2 | 1 | 1 | 0 | ASCE | - |
| 2025-11-22 | pre | 2 | 1 | 1 | 0 | DESC | - |
| 2025-11-28 | pre | 2 | 1 | 1 | 0 | ASCE | - |
| 2025-12-04 | pre | 2 | 1 | 1 | 1 | DESC | intermediate |
| 2025-12-10 | pre | 2 | 1 | 1 | 1 | ASCE | compact |
| 2025-12-16 | pre | 2 | 1 | 1 | 1 | DESC | compact |
| 2025-12-22 | pre | 2 | 1 | 1 | 0 | ASCE | - |
| 2025-12-28 | pre | 2 | 1 | 1 | 1 | DESC | distributed |
| 2026-01-03 | pre | 2 | 1 | 1 | 1 | ASCE | intermediate |
| 2026-01-15 | pre | 2 | 1 | 1 | 1 | ASCE | intermediate |
| 2026-01-27 | pre | 2 | 1 | 1 | 1 | ASCE | intermediate |
| 2026-02-08 | post | 2 | 1 | 1 | 0 | ASCE | - |
| 2026-02-20 | post | 2 | 1 | 1 | 0 | ASCE | - |
| 2026-03-16 | post | 2 | 1 | 1 | 1 | ASCE | intermediate |
| 2026-04-09 | post | 2 | 1 | 1 | 0 | ASCE | - |
| 2026-04-15 | post | 2 | 1 | 1 | 0 | DESC | - |
| 2026-04-27 | post | 2 | 1 | 1 | 0 | DESC | - |

## The payload groups of the committed window file

The committed window file holds **57 windows** for **33 unique chips**: `payload_group_size` counts the byte-identical duplicates of one measurement window. 17 of the chips carry an echo and 16 do not, so the 24 duplicates split 13 onto echo-bearing chips and 11 onto echo-free ones. The cache's own tally of the same duplication is 24 rows against 57 windows (17 + 24); the reference calls the same quantity `n_windows` = 41.

| window digest | day | windows behind it | echo | sections |
| --- | --- | --- | --- | --- |
| 6c05e03c | 2025-08-06 | 1 | yes | section-3 |
| 34ff57ef | 2025-08-12 | 1 | no | section-3 |
| b7475bd2 | 2025-08-12 | 1 | no | section-3 |
| f4cf9b6f | 2025-08-24 | 1 | no | section-3 |
| d8c52627 | 2025-09-05 | 1 | yes | section-3 |
| 7615e5de | 2025-09-11 | 1 | yes | section-3 |
| 226a2a50 | 2025-09-17 | 1 | yes | section-3 |
| 06e39942 | 2025-09-23 | 1 | no | section-3 |
| fd65890a | 2025-09-29 | 1 | no | section-3 |
| 3c9dfd0b | 2025-10-05 | 2 | no | section-2, section-3 |
| a40f3f2c | 2025-10-11 | 2 | yes | section-2, section-3 |
| 1e494779 | 2025-10-17 | 2 | yes | section-2, section-3 |
| deb437c1 | 2025-10-23 | 2 | yes | section-2, section-3 |
| 6c361c78 | 2025-10-29 | 2 | yes | section-2, section-3 |
| 4a047373 | 2025-11-04 | 2 | no | section-2, section-3 |
| b9579f83 | 2025-11-10 | 2 | yes | section-2, section-3 |
| c38ab6fa | 2025-11-16 | 2 | no | section-2, section-3 |
| d119a4ae | 2025-11-22 | 2 | no | section-2, section-3 |
| a3767e49 | 2025-11-28 | 2 | no | section-2, section-3 |
| 4ab34cb4 | 2025-12-04 | 2 | yes | section-2, section-3 |
| 9f6afc59 | 2025-12-10 | 2 | yes | section-2, section-3 |
| d5cad0f9 | 2025-12-16 | 2 | yes | section-2, section-3 |
| cdab9cd4 | 2025-12-22 | 2 | no | section-2, section-3 |
| 6c2164f0 | 2025-12-28 | 2 | yes | section-2, section-3 |
| 51a954be | 2026-01-03 | 2 | yes | section-2, section-3 |
| bfb8bdb8 | 2026-01-15 | 2 | yes | section-2, section-3 |
| 05994a29 | 2026-01-27 | 2 | yes | section-2, section-3 |
| b1c1e50d | 2026-02-08 | 2 | no | section-2, section-3 |
| b9880a19 | 2026-02-20 | 2 | no | section-2, section-3 |
| eaf450d9 | 2026-03-16 | 2 | yes | section-2, section-3 |
| cb0157e3 | 2026-04-09 | 2 | no | section-2, section-3 |
| a9b60ed2 | 2026-04-15 | 2 | no | section-2, section-3 |
| 2e98fd96 | 2026-04-27 | 2 | no | section-2, section-3 |

## The caveats of the record

**The pipeline column is a different generation.** The column was computed on the asset-point chip of the whole-scene window, while this package recomputes the echo mask of the committed 7x7 window. The two therefore disagree, and they do: the committed column matches the recomputed `A` on **0 of 17** echo chips. Its pre-event median is 75.000 masked pixels against 3.000 recomputed; the figure draws the recomputed `A`, because that is the quantity this record's mask rule defines, and pins the committed column beside it rather than plotting it.

**A control with no variance.** `coherence` is 0.998917748918 on all 33 chips — 1 distinct value — so it separates nothing: it is reported, never drawn. All 8 co-variate controls of the channel table (`wind_speed_ms`, `temperature_c`, `brightness_ratio`, `intensity`, `coherence`, `coherence_gamma2`, `cumulative_mm`, `structural_frequency_hz`) are pinned in the JSON, and one of them (`structural_frequency_hz`) is committed for only 8 of the 33 chips, so it is not drawn either.

## The pins and the provenance

2637 pin comparisons against five committed layers (the channel table's meta, the mask-cache manifest, the committed mask cache, the paper result JSON and the analysis reference): **2633 exact**, **4 within 1e-12**, **0 failed**. The recomputed number is in `recomputed`, the committed one in `committed`, and no artifact is written before they agree.

| committed input | sha256 |
| --- | --- |
| `data/ywf/reference/ywf_bursts.json` | `29e48c491585c4e7127e088d4000c051577f28abad137e900ae8b4b30be9500a` |
| `data/ywf/ywf_ocv_channels.csv` | `b5b17240246f613357814068ec3df62fceedead11b8426e01e3ebd6b006c5f33` |
| `data/ywf/ywf_windows_mask_cache.txt` | `0c29524d8569cd2861bbae187501d3ac270c49ad94dbde2d597109f5e70bfa6b` |
| `data/ywf/ywf_windows_mask_cache.manifest.json` | `a14a32d23620c01bf2d051c2e912d9e8d8d0be9cfb49245a4b4fc4e8afb6d0c8` |
| `data/ywf/ywf_measurements_full.txt` | `652a6c1d57766c319b09fa23f73865b9fed034bd6b58b3dd104497a571e34950` |
| `data/ywf/reference/ywf_geometry.json` | `57ac7a9815857df20ad55e010e0111b537ceabaa11c671cbf1d08316c47237f7` |
| `data/ywf/reference/ywf_osm_road.json` | `f88c4d6f2cd58c85bc932d3016e5455b3e1cb46339fda039005a46c045dd1699` |
| `data/ywf/ywf_segments.txt` | `5d48b2b1118d952edb1c3a212136e34047c8d6be7e2100204c3b0efd37682f09` |
| `data/ywf/ywf_windows_full.txt` | `41b8838e80aebb2184a5da8bf9941ff55afc0a280a24742cfa620caa053b21bf` |

Figures: `figures/ywf/fig_ywf_ocv_series.png` (600 dpi). JSON: `data/ywf/fig_ywf_ocv_series.json`.
