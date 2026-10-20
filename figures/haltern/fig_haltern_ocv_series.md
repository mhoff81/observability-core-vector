# Haltern AV9 Wind Park (Kreis Recklinghausen, NRW) — the echo-mask vector over the whole record

**Site** Haltern AV9 Wind Park (Kreis Recklinghausen, NRW) · **asset** Haltern AV9 Wind Turbine (Nordex N149/4.5) · **event** 2021-09-29 (Wednesday, inside the axis, inside a step) · **record** 2021-01-01 -> 2021-10-17 (290 days) · **chips** 120 unique 7x7 chips on 120 acquisition days, 66 carrying an echo

```
OCV  = [gamma2, P, D]               observability core vector
x    = the acquisition day, one dot per chip (2021-01-01 -> 2021-10-17)
y    = the value over that dimension's own pre-window median, taken
       over the 62 pre-event echo chips; 1.0 = that median
mask = intensity >= 0.30*peak AND >= 5.0*np.median(intensity), A >= 2 (window 7x7)
step = 3 d or less (the record alternates in two-day blocks); 9 of 119 steps are longer
       and are drawn as bands — a hole is never interpolated over
post = 4 of 12 post-event chips carry an echo, so the post side carries
       four points and no line, no median and no trend
```

## The record

| fact | value |
| --- | --- |
| chips | **120** unique 7x7 chips behind **168** committed windows (48 byte-identical duplicates of the windows file) |
| acquisition days | **120**, 2021-01-01 -> 2021-10-17 (290 days): steps of 1..3 days (the two-day alternation blocks) with 9 longer steps and 110 compliant steps |
| pre-event | **108 chips on 108 days**, 2021-01-01 -> 2021-09-17 (260 d), 62 with an echo (57.4 %) |
| post-event | **12 chips on 12 days**, 2021-10-01 -> 2021-10-17 (17 d), 4 with an echo (33.3 %) |
| the event falls inside a step | 2021-09-29 sits in the **14-day** step 2021-09-17 -> 2021-10-01; no chip on the event day, so the two stretches meet at a band and not at a point |
| the pass geometries | no acquisition day carries two chips: the record's 120 chips sit on the two committed pass labels (`afternoon` 62, `morning` 58), so the orbit offset the x axis adds is cosmetic and never separates two chips of one day |
| the `P` reference | `P` is a window-wide reference, not a per-chip quantity: the committed cache carries `P_reference_chips` = 66 — the same 66 echo-bearing chips of the whole record of which the mean pixel frequency inside the own mask is taken |
| the pipeline `coherence_masked_pixels` column | it runs 3.000 .. 257.650 masked pixels over the 108 pre-event rows (median 40.000) and agrees with the recomputed `A` on **6 of 66** echo chips — see the caveat below |
| the constant control `coherence` | 0.999285714286 on every one of the 120 chips (1 distinct value): reported as a control, never drawn |

## The six panels

Every panel draws the same axis: one dot per echo-bearing chip, y = that dimension's value over **its own** pre-window median. The committed quartiles of the pre window are divided by the same median to make the box, so every drawn ordinate is a pinned number.

| panel | dimension | y scale | base (pre median, raw) | pre n | pre min (ratio) | pre p25 (ratio) | pre p75 (ratio) | pre max (ratio) | post n | post readings (ratio) | post median (raw) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| a | gamma2 | log | 0.503 | 62 | 0.006 | 0.419 | 1.807 | 1.986 | 4 | 0.192, 1.261, 0.864, 0.246 | 0.279 |
| b | A | log | 3.000 | 62 | 0.667 | 0.667 | 1.333 | 3.333 | 4 | 2.000, 0.667, 2.000, 0.667 | 4.000 |
| c | P | linear | 0.090 | 62 | 0.393 | 0.841 | 1.093 | 1.514 | 4 | 1.121, 1.093, 0.897, 1.262 | 0.100 |
| d | D | log | 0.475 | 62 | 0.172 | 0.632 | 1.711 | 2.105 | 4 | 0.602, 2.105, 0.505, 1.053 | 0.393 |
| e | F | log | 2.000 | 62 | 0.500 | 0.500 | 1.000 | 2.000 | 4 | 1.500, 0.500, 1.500, 1.000 | 2.500 |
| f | S | log | 0.867 | 62 | 0.288 | 0.577 | 1.681 | 4.160 | 4 | 3.041, 0.577, 2.150, 1.731 | 1.682 |

`P` is the only linear panel: it is the mean pixel frequency inside the mask over the 66 reference chips — a bounded frequency, so its ratio is read as a difference. The other five are multiplicative magnitudes (49 is the window, the denominator of `A` and `D`); all six pre windows stay inside 0.005..4.2 of 1.0, so no log axis hides a level. The centre line of every box is 1.0 by construction: it *is* the pre-window median the ordinate is divided by.

## The steps of the record

Acquisitions alternate in blocks of two days: 9 of the 119 steps between two consecutive acquisition days are longer than the 3-day cadence. Each of them is drawn as a grey band and **no line crosses it** — a line from one chip to the next would be an interpolation over a hole this record does not have data for.

| step | how many | what it is | the steps |
| --- | --- | --- | --- |
| 1 d | 57 | the cadence | 2021-01-01 -> 2021-01-02, 2021-01-04 -> 2021-01-05, 2021-01-07 -> 2021-01-08, 2021-01-10 -> 2021-01-11, 2021-01-13 -> 2021-01-14, 2021-01-16 -> 2021-01-17, 2021-02-03 -> 2021-02-04, 2021-02-06 -> 2021-02-07, 2021-02-09 -> 2021-02-10, 2021-02-15 -> 2021-02-16, 2021-02-18 -> 2021-02-19, 2021-03-02 -> 2021-03-03, 2021-03-05 -> 2021-03-06, 2021-03-08 -> 2021-03-09, 2021-03-11 -> 2021-03-12, 2021-03-14 -> 2021-03-15, 2021-03-17 -> 2021-03-18, 2021-04-01 -> 2021-04-02, 2021-04-04 -> 2021-04-05, 2021-04-07 -> 2021-04-08, 2021-04-10 -> 2021-04-11, 2021-04-13 -> 2021-04-14, 2021-04-16 -> 2021-04-17, 2021-05-01 -> 2021-05-02, 2021-05-04 -> 2021-05-05, 2021-05-07 -> 2021-05-08, 2021-05-10 -> 2021-05-11, 2021-05-13 -> 2021-05-14, 2021-05-16 -> 2021-05-17, 2021-06-06 -> 2021-06-07, 2021-06-09 -> 2021-06-10, 2021-06-12 -> 2021-06-13, 2021-06-15 -> 2021-06-16, 2021-06-18 -> 2021-06-19, 2021-07-06 -> 2021-07-07, 2021-07-09 -> 2021-07-10, 2021-07-12 -> 2021-07-13, 2021-07-15 -> 2021-07-16, 2021-07-18 -> 2021-07-19, 2021-08-02 -> 2021-08-03, 2021-08-05 -> 2021-08-06, 2021-08-08 -> 2021-08-09, 2021-08-11 -> 2021-08-12, 2021-08-14 -> 2021-08-15, 2021-08-17 -> 2021-08-18, 2021-09-01 -> 2021-09-02, 2021-09-04 -> 2021-09-05, 2021-09-07 -> 2021-09-08, 2021-09-10 -> 2021-09-11, 2021-09-13 -> 2021-09-14, 2021-09-16 -> 2021-09-17, 2021-10-01 -> 2021-10-02, 2021-10-04 -> 2021-10-05, 2021-10-07 -> 2021-10-08, 2021-10-10 -> 2021-10-11, 2021-10-13 -> 2021-10-14, 2021-10-16 -> 2021-10-17 |
| 2 d | 50 | the cadence | 2021-01-02 -> 2021-01-04, 2021-01-05 -> 2021-01-07, 2021-01-08 -> 2021-01-10, 2021-01-11 -> 2021-01-13, 2021-01-14 -> 2021-01-16, 2021-02-01 -> 2021-02-03, 2021-02-04 -> 2021-02-06, 2021-02-07 -> 2021-02-09, 2021-02-10 -> 2021-02-12, 2021-02-16 -> 2021-02-18, 2021-03-03 -> 2021-03-05, 2021-03-06 -> 2021-03-08, 2021-03-09 -> 2021-03-11, 2021-03-12 -> 2021-03-14, 2021-03-15 -> 2021-03-17, 2021-04-02 -> 2021-04-04, 2021-04-05 -> 2021-04-07, 2021-04-08 -> 2021-04-10, 2021-04-11 -> 2021-04-13, 2021-04-14 -> 2021-04-16, 2021-05-02 -> 2021-05-04, 2021-05-05 -> 2021-05-07, 2021-05-08 -> 2021-05-10, 2021-05-11 -> 2021-05-13, 2021-05-14 -> 2021-05-16, 2021-06-04 -> 2021-06-06, 2021-06-07 -> 2021-06-09, 2021-06-10 -> 2021-06-12, 2021-06-13 -> 2021-06-15, 2021-06-16 -> 2021-06-18, 2021-07-01 -> 2021-07-03, 2021-07-07 -> 2021-07-09, 2021-07-10 -> 2021-07-12, 2021-07-13 -> 2021-07-15, 2021-07-16 -> 2021-07-18, 2021-08-03 -> 2021-08-05, 2021-08-06 -> 2021-08-08, 2021-08-09 -> 2021-08-11, 2021-08-12 -> 2021-08-14, 2021-08-15 -> 2021-08-17, 2021-09-02 -> 2021-09-04, 2021-09-05 -> 2021-09-07, 2021-09-08 -> 2021-09-10, 2021-09-11 -> 2021-09-13, 2021-09-14 -> 2021-09-16, 2021-10-02 -> 2021-10-04, 2021-10-05 -> 2021-10-07, 2021-10-08 -> 2021-10-10, 2021-10-11 -> 2021-10-13, 2021-10-14 -> 2021-10-16 |
| 3 d | 3 | the cadence | 2021-02-12 -> 2021-02-15, 2021-06-01 -> 2021-06-04, 2021-07-03 -> 2021-07-06 |
| 11 d | 1 | a hole, drawn as a band | 2021-02-19 -> 2021-03-02 |
| 12 d | 1 | a hole, drawn as a band | 2021-06-19 -> 2021-07-01 |
| 14 d | 5 | a hole, drawn as a band | 2021-03-18 -> 2021-04-01, 2021-04-17 -> 2021-05-01, 2021-07-19 -> 2021-08-02, 2021-08-18 -> 2021-09-01, 2021-09-17 -> 2021-10-01 |
| 15 d | 2 | a hole, drawn as a band | 2021-01-17 -> 2021-02-01, 2021-05-17 -> 2021-06-01 |

## The two stretches

| stretch | first day | last day | calendar days | chips | acquisition days | committed windows | duplicates | with echo | without echo | echo rate | orbits | steps | longer steps | step lengths |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| pre-collapse (healthy) | 2021-01-01 | 2021-09-17 | 260 | 108 | 108 | 156 | 48 | 62 | 46 | 0.574 | Ascending, Descending | 107 | 8 | 1 d, 2 d, 3 d, 11 d, 12 d, 14 d, 15 d |
| post-collapse | 2021-10-01 | 2021-10-17 | 17 | 12 | 12 | 12 | 0 | 4 | 8 | 0.333 | Ascending, Descending | 11 | 0 | 1 d, 2 d |

The two stretches are unequal in every respect: **260 days and 108 chips** before the collapse against **17 days and 12 chips** after it. The collapse step (2021-09-17 -> 2021-10-01) is the only step of the record that holds the event date, and the dashed cut of every panel falls inside it: the two stretches of the figure meet at a grey band, never at a point, and the nearest chips sit 12 days before and 2 days after the event date.

## The coverage: what the collapse left behind

| state | echo | no echo |
| --- | --- | --- |
| pre | 62 | 46 |
| post | 4 | 8 |

Fisher's exact test on that 2x2 table: **p = 0.1343** (two-sided, alpha = 0.05). With 62 of 108 pre-event chips carrying an echo and 4 of 12 post-event chips doing so, the record does not reject "the coverage is unchanged" — and it cannot, with 12 post chips.

**The power question of the record.** If all post-event chips came back echo-free, it would take **4 post-event chips** for the same test to fall below alpha = 0.05 against the pre-event rate of 57.4 %. The record has 12 chips there, 4 of which echo.

4 of 12 post-event chips carry an echo, so the post-state distribution is drawn from 4 echoes: the collapse contrast of this site is an echo-coverage statement with this caveat attached, not a damage verdict.

## Where the echoes sit over the two pass geometries

| orbit | chips | with echo | without echo | echo rate | pre chips | pre echo rate | post chips | post echo rate | Fisher within the orbit |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ASCENDING | 62 | 33 | 29 | 0.532 | 56 | 0.554 | 6 | 0.333 | 0.405 |
| DESCENDING | 58 | 33 | 25 | 0.569 | 52 | 0.596 | 6 | 0.333 | 0.387 |

Across the whole record the two geometries are indistinguishable on coverage (Fisher p = 0.717). Both directions straddle the collapse date and the four post-event echo chips are two ascending and two descending passes: this block describes where the echoes sit over the four pass geometries, it is not an orbit contrast of the states.

## The composition of the record

| month | pre chips | post chips |
| --- | --- | --- |
| Jan | 12 | 0 |
| Feb | 12 | 0 |
| Mar | 12 | 0 |
| Apr | 12 | 0 |
| May | 12 | 0 |
| Jun | 12 | 0 |
| Jul | 12 | 0 |
| Aug | 12 | 0 |
| Sep | 12 | 0 |
| Oct | 0 | 12 |

108 of the 120 chips are pre-event and 12 post-event. The asset carries 2 labelled sections (`turbine-mast-section-1`, `turbine-mast-section-0`) — two sections of one turbine mast, not two structures — so no panel of this figure is a comparison across sections.

## The echo geometries

| state | echo chips | compact | distributed | intermediate |
| --- | --- | --- | --- | --- |
| pre | 62 | 22 | 9 | 31 |
| post | 4 | 2 | 2 | 0 |

## The record, day by day

Every acquisition day of the record, with what was requested, what the mask rule found and which pass geometry carried it. This is the committed timeline (`n_windows`, `n_unique_chips`, `n_duplicated`, `n_with_echo`, `state`) as the channel table recomputes it.

| day | state | committed windows | unique chips | duplicates | with echo | orbits | echo geometry |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 2021-01-01 | pre | 2 | 1 | 1 | 0 | ASCE | - |
| 2021-01-02 | pre | 2 | 1 | 1 | 0 | ASCE | - |
| 2021-01-04 | pre | 2 | 1 | 1 | 0 | DESC | - |
| 2021-01-05 | pre | 2 | 1 | 1 | 1 | DESC | distributed |
| 2021-01-07 | pre | 2 | 1 | 1 | 1 | ASCE | intermediate |
| 2021-01-08 | pre | 2 | 1 | 1 | 0 | ASCE | - |
| 2021-01-10 | pre | 2 | 1 | 1 | 1 | DESC | intermediate |
| 2021-01-11 | pre | 2 | 1 | 1 | 0 | DESC | - |
| 2021-01-13 | pre | 2 | 1 | 1 | 0 | ASCE | - |
| 2021-01-14 | pre | 2 | 1 | 1 | 0 | ASCE | - |
| 2021-01-16 | pre | 2 | 1 | 1 | 1 | DESC | distributed |
| 2021-01-17 | pre | 2 | 1 | 1 | 1 | DESC | intermediate |
| 2021-02-01 | pre | 2 | 1 | 1 | 1 | ASCE | intermediate |
| 2021-02-03 | pre | 2 | 1 | 1 | 1 | DESC | intermediate |
| 2021-02-04 | pre | 2 | 1 | 1 | 1 | DESC | distributed |
| 2021-02-06 | pre | 2 | 1 | 1 | 1 | ASCE | intermediate |
| 2021-02-07 | pre | 2 | 1 | 1 | 1 | ASCE | compact |
| 2021-02-09 | pre | 2 | 1 | 1 | 0 | DESC | - |
| 2021-02-10 | pre | 2 | 1 | 1 | 1 | DESC | compact |
| 2021-02-12 | pre | 2 | 1 | 1 | 1 | ASCE | intermediate |
| 2021-02-15 | pre | 2 | 1 | 1 | 1 | DESC | compact |
| 2021-02-16 | pre | 2 | 1 | 1 | 1 | DESC | compact |
| 2021-02-18 | pre | 2 | 1 | 1 | 0 | ASCE | - |
| 2021-02-19 | pre | 2 | 1 | 1 | 1 | ASCE | compact |
| 2021-03-02 | pre | 2 | 1 | 1 | 1 | ASCE | compact |
| 2021-03-03 | pre | 2 | 1 | 1 | 0 | ASCE | - |
| 2021-03-05 | pre | 2 | 1 | 1 | 0 | DESC | - |
| 2021-03-06 | pre | 2 | 1 | 1 | 1 | DESC | compact |
| 2021-03-08 | pre | 2 | 1 | 1 | 0 | ASCE | - |
| 2021-03-09 | pre | 2 | 1 | 1 | 0 | ASCE | - |
| 2021-03-11 | pre | 2 | 1 | 1 | 0 | DESC | - |
| 2021-03-12 | pre | 2 | 1 | 1 | 1 | DESC | intermediate |
| 2021-03-14 | pre | 2 | 1 | 1 | 1 | ASCE | compact |
| 2021-03-15 | pre | 2 | 1 | 1 | 1 | ASCE | compact |
| 2021-03-17 | pre | 2 | 1 | 1 | 1 | DESC | intermediate |
| 2021-03-18 | pre | 2 | 1 | 1 | 0 | DESC | - |
| 2021-04-01 | pre | 2 | 1 | 1 | 1 | ASCE | compact |
| 2021-04-02 | pre | 2 | 1 | 1 | 0 | ASCE | - |
| 2021-04-04 | pre | 2 | 1 | 1 | 0 | DESC | - |
| 2021-04-05 | pre | 2 | 1 | 1 | 1 | DESC | distributed |
| 2021-04-07 | pre | 2 | 1 | 1 | 0 | ASCE | - |
| 2021-04-08 | pre | 2 | 1 | 1 | 1 | ASCE | intermediate |
| 2021-04-10 | pre | 2 | 1 | 1 | 0 | DESC | - |
| 2021-04-11 | pre | 2 | 1 | 1 | 1 | DESC | intermediate |
| 2021-04-13 | pre | 2 | 1 | 1 | 1 | ASCE | compact |
| 2021-04-14 | pre | 2 | 1 | 1 | 0 | ASCE | - |
| 2021-04-16 | pre | 2 | 1 | 1 | 1 | DESC | intermediate |
| 2021-04-17 | pre | 2 | 1 | 1 | 1 | DESC | intermediate |
| 2021-05-01 | pre | 1 | 1 | 0 | 1 | ASCE | compact |
| 2021-05-02 | pre | 1 | 1 | 0 | 0 | ASCE | - |
| 2021-05-04 | pre | 1 | 1 | 0 | 0 | DESC | - |
| 2021-05-05 | pre | 1 | 1 | 0 | 0 | DESC | - |
| 2021-05-07 | pre | 1 | 1 | 0 | 1 | ASCE | intermediate |
| 2021-05-08 | pre | 1 | 1 | 0 | 0 | ASCE | - |
| 2021-05-10 | pre | 1 | 1 | 0 | 1 | DESC | intermediate |
| 2021-05-11 | pre | 1 | 1 | 0 | 0 | DESC | - |
| 2021-05-13 | pre | 1 | 1 | 0 | 1 | ASCE | compact |
| 2021-05-14 | pre | 1 | 1 | 0 | 1 | ASCE | compact |
| 2021-05-16 | pre | 1 | 1 | 0 | 0 | DESC | - |
| 2021-05-17 | pre | 1 | 1 | 0 | 0 | DESC | - |
| 2021-06-01 | pre | 1 | 1 | 0 | 1 | ASCE | intermediate |
| 2021-06-04 | pre | 1 | 1 | 0 | 1 | DESC | distributed |
| 2021-06-06 | pre | 1 | 1 | 0 | 0 | ASCE | - |
| 2021-06-07 | pre | 1 | 1 | 0 | 0 | ASCE | - |
| 2021-06-09 | pre | 1 | 1 | 0 | 1 | DESC | compact |
| 2021-06-10 | pre | 1 | 1 | 0 | 1 | DESC | intermediate |
| 2021-06-12 | pre | 1 | 1 | 0 | 1 | ASCE | intermediate |
| 2021-06-13 | pre | 1 | 1 | 0 | 0 | ASCE | - |
| 2021-06-15 | pre | 1 | 1 | 0 | 1 | DESC | compact |
| 2021-06-16 | pre | 1 | 1 | 0 | 1 | DESC | intermediate |
| 2021-06-18 | pre | 1 | 1 | 0 | 1 | ASCE | distributed |
| 2021-06-19 | pre | 1 | 1 | 0 | 1 | ASCE | intermediate |
| 2021-07-01 | pre | 1 | 1 | 0 | 1 | ASCE | compact |
| 2021-07-03 | pre | 1 | 1 | 0 | 0 | DESC | - |
| 2021-07-06 | pre | 1 | 1 | 0 | 0 | ASCE | - |
| 2021-07-07 | pre | 1 | 1 | 0 | 1 | ASCE | compact |
| 2021-07-09 | pre | 1 | 1 | 0 | 1 | DESC | intermediate |
| 2021-07-10 | pre | 1 | 1 | 0 | 0 | DESC | - |
| 2021-07-12 | pre | 1 | 1 | 0 | 1 | ASCE | intermediate |
| 2021-07-13 | pre | 1 | 1 | 0 | 1 | ASCE | intermediate |
| 2021-07-15 | pre | 1 | 1 | 0 | 0 | DESC | - |
| 2021-07-16 | pre | 1 | 1 | 0 | 1 | DESC | distributed |
| 2021-07-18 | pre | 1 | 1 | 0 | 1 | ASCE | intermediate |
| 2021-07-19 | pre | 1 | 1 | 0 | 1 | ASCE | compact |
| 2021-08-02 | pre | 1 | 1 | 0 | 0 | DESC | - |
| 2021-08-03 | pre | 1 | 1 | 0 | 1 | DESC | intermediate |
| 2021-08-05 | pre | 1 | 1 | 0 | 0 | ASCE | - |
| 2021-08-06 | pre | 1 | 1 | 0 | 1 | ASCE | intermediate |
| 2021-08-08 | pre | 1 | 1 | 0 | 1 | DESC | compact |
| 2021-08-09 | pre | 1 | 1 | 0 | 0 | DESC | - |
| 2021-08-11 | pre | 1 | 1 | 0 | 0 | ASCE | - |
| 2021-08-12 | pre | 1 | 1 | 0 | 0 | ASCE | - |
| 2021-08-14 | pre | 1 | 1 | 0 | 1 | DESC | intermediate |
| 2021-08-15 | pre | 1 | 1 | 0 | 1 | DESC | compact |
| 2021-08-17 | pre | 1 | 1 | 0 | 0 | ASCE | - |
| 2021-08-18 | pre | 1 | 1 | 0 | 1 | ASCE | intermediate |
| 2021-09-01 | pre | 1 | 1 | 0 | 1 | DESC | compact |
| 2021-09-02 | pre | 1 | 1 | 0 | 0 | DESC | - |
| 2021-09-04 | pre | 1 | 1 | 0 | 1 | ASCE | distributed |
| 2021-09-05 | pre | 1 | 1 | 0 | 0 | ASCE | - |
| 2021-09-07 | pre | 1 | 1 | 0 | 0 | DESC | - |
| 2021-09-08 | pre | 1 | 1 | 0 | 0 | DESC | - |
| 2021-09-10 | pre | 1 | 1 | 0 | 0 | ASCE | - |
| 2021-09-11 | pre | 1 | 1 | 0 | 0 | ASCE | - |
| 2021-09-13 | pre | 1 | 1 | 0 | 1 | DESC | distributed |
| 2021-09-14 | pre | 1 | 1 | 0 | 1 | DESC | intermediate |
| 2021-09-16 | pre | 1 | 1 | 0 | 1 | ASCE | intermediate |
| 2021-09-17 | pre | 1 | 1 | 0 | 1 | ASCE | intermediate |
| 2021-10-01 | post | 1 | 1 | 0 | 0 | DESC | - |
| 2021-10-02 | post | 1 | 1 | 0 | 1 | DESC | distributed |
| 2021-10-04 | post | 1 | 1 | 0 | 0 | ASCE | - |
| 2021-10-05 | post | 1 | 1 | 0 | 0 | ASCE | - |
| 2021-10-07 | post | 1 | 1 | 0 | 0 | DESC | - |
| 2021-10-08 | post | 1 | 1 | 0 | 1 | DESC | compact |
| 2021-10-10 | post | 1 | 1 | 0 | 1 | ASCE | distributed |
| 2021-10-11 | post | 1 | 1 | 0 | 0 | ASCE | - |
| 2021-10-13 | post | 1 | 1 | 0 | 0 | DESC | - |
| 2021-10-14 | post | 1 | 1 | 0 | 0 | DESC | - |
| 2021-10-16 | post | 1 | 1 | 0 | 0 | ASCE | - |
| 2021-10-17 | post | 1 | 1 | 0 | 1 | ASCE | compact |

## The payload groups of the committed window file

The committed window file holds **168 windows** for **120 unique chips**: `payload_group_size` counts the byte-identical duplicates of one measurement window. 66 of the chips carry an echo and 54 do not, so the 48 duplicates split 28 onto echo-bearing chips and 20 onto echo-free ones. The cache's own tally of the same duplication is 48 rows against 168 windows (66 + 48); the reference calls the same quantity `n_windows` = 114.

| window digest | day | windows behind it | echo | sections |
| --- | --- | --- | --- | --- |
| 87f3deb6 | 2021-01-01 | 2 | no | section-0, section-1 |
| 94db0e13 | 2021-01-02 | 2 | no | section-0, section-1 |
| 4f68e5b6 | 2021-01-04 | 2 | no | section-0, section-1 |
| bee0654b | 2021-01-05 | 2 | yes | section-0, section-1 |
| bb01065c | 2021-01-07 | 2 | yes | section-0, section-1 |
| 6870acd6 | 2021-01-08 | 2 | no | section-0, section-1 |
| 5e69cd4a | 2021-01-10 | 2 | yes | section-0, section-1 |
| c3d706d8 | 2021-01-11 | 2 | no | section-0, section-1 |
| 96f5e4bf | 2021-01-13 | 2 | no | section-0, section-1 |
| e736422c | 2021-01-14 | 2 | no | section-0, section-1 |
| 8a32731b | 2021-01-16 | 2 | yes | section-0, section-1 |
| 50c3638e | 2021-01-17 | 2 | yes | section-0, section-1 |
| 56184a8a | 2021-02-01 | 2 | yes | section-0, section-1 |
| be2472b1 | 2021-02-03 | 2 | yes | section-0, section-1 |
| 3b6bc7a6 | 2021-02-04 | 2 | yes | section-0, section-1 |
| 058d6440 | 2021-02-06 | 2 | yes | section-0, section-1 |
| f25cdad9 | 2021-02-07 | 2 | yes | section-0, section-1 |
| deedec21 | 2021-02-09 | 2 | no | section-0, section-1 |
| 2828183c | 2021-02-10 | 2 | yes | section-0, section-1 |
| ceaa449c | 2021-02-12 | 2 | yes | section-0, section-1 |
| 3a1a6fd5 | 2021-02-15 | 2 | yes | section-0, section-1 |
| d3d17b69 | 2021-02-16 | 2 | yes | section-0, section-1 |
| 9af25d0b | 2021-02-18 | 2 | no | section-0, section-1 |
| cca2b921 | 2021-02-19 | 2 | yes | section-0, section-1 |
| 42fa10fe | 2021-03-02 | 2 | yes | section-0, section-1 |
| 3760fb60 | 2021-03-03 | 2 | no | section-0, section-1 |
| e5642685 | 2021-03-05 | 2 | no | section-0, section-1 |
| 6889eae1 | 2021-03-06 | 2 | yes | section-0, section-1 |
| f5cfe68e | 2021-03-08 | 2 | no | section-0, section-1 |
| cd921283 | 2021-03-09 | 2 | no | section-0, section-1 |
| ff76d2a5 | 2021-03-11 | 2 | no | section-0, section-1 |
| 1ed793db | 2021-03-12 | 2 | yes | section-0, section-1 |
| 0c3a022c | 2021-03-14 | 2 | yes | section-0, section-1 |
| b72b3422 | 2021-03-15 | 2 | yes | section-0, section-1 |
| f1a096d3 | 2021-03-17 | 2 | yes | section-0, section-1 |
| 39e14df3 | 2021-03-18 | 2 | no | section-0, section-1 |
| 2bec5966 | 2021-04-01 | 2 | yes | section-0, section-1 |
| 42ceed76 | 2021-04-02 | 2 | no | section-0, section-1 |
| f7573b2f | 2021-04-04 | 2 | no | section-0, section-1 |
| d67e6842 | 2021-04-05 | 2 | yes | section-0, section-1 |
| 5d702839 | 2021-04-07 | 2 | no | section-0, section-1 |
| d5867cc9 | 2021-04-08 | 2 | yes | section-0, section-1 |
| ee2af69b | 2021-04-10 | 2 | no | section-0, section-1 |
| ef0e7dab | 2021-04-11 | 2 | yes | section-0, section-1 |
| c6217f82 | 2021-04-13 | 2 | yes | section-0, section-1 |
| 8d5cdec7 | 2021-04-14 | 2 | no | section-0, section-1 |
| 984e1d33 | 2021-04-16 | 2 | yes | section-0, section-1 |
| bde33c49 | 2021-04-17 | 2 | yes | section-0, section-1 |
| 9ff7301b | 2021-05-01 | 1 | yes | section-0 |
| fe65808c | 2021-05-02 | 1 | no | section-0 |
| 34f0630a | 2021-05-04 | 1 | no | section-0 |
| 23d26f0e | 2021-05-05 | 1 | no | section-0 |
| ed02e6ad | 2021-05-07 | 1 | yes | section-0 |
| 1d1e8054 | 2021-05-08 | 1 | no | section-0 |
| 0e0eac92 | 2021-05-10 | 1 | yes | section-0 |
| d84b53ba | 2021-05-11 | 1 | no | section-0 |
| 97f48586 | 2021-05-13 | 1 | yes | section-0 |
| 7ac827d9 | 2021-05-14 | 1 | yes | section-0 |
| 05b66376 | 2021-05-16 | 1 | no | section-0 |
| 687c0ca5 | 2021-05-17 | 1 | no | section-0 |
| 5370cb71 | 2021-06-01 | 1 | yes | section-0 |
| c95c7e25 | 2021-06-04 | 1 | yes | section-0 |
| b1f7edc0 | 2021-06-06 | 1 | no | section-0 |
| 0507e478 | 2021-06-07 | 1 | no | section-0 |
| f4e4f626 | 2021-06-09 | 1 | yes | section-0 |
| 80faa10e | 2021-06-10 | 1 | yes | section-0 |
| 881e10d2 | 2021-06-12 | 1 | yes | section-0 |
| 11ebe1dd | 2021-06-13 | 1 | no | section-0 |
| 8570da6d | 2021-06-15 | 1 | yes | section-0 |
| 91f4b0b7 | 2021-06-16 | 1 | yes | section-0 |
| 5bb6757d | 2021-06-18 | 1 | yes | section-0 |
| e2e52f0c | 2021-06-19 | 1 | yes | section-0 |
| ccb3fc31 | 2021-07-01 | 1 | yes | section-0 |
| 47286ae9 | 2021-07-03 | 1 | no | section-0 |
| 3beb92cf | 2021-07-06 | 1 | no | section-0 |
| 437c7330 | 2021-07-07 | 1 | yes | section-0 |
| 23ae1b91 | 2021-07-09 | 1 | yes | section-0 |
| ebd4ae10 | 2021-07-10 | 1 | no | section-0 |
| eaf6e071 | 2021-07-12 | 1 | yes | section-0 |
| e8241c18 | 2021-07-13 | 1 | yes | section-0 |
| 351c642d | 2021-07-15 | 1 | no | section-0 |
| ade5ec45 | 2021-07-16 | 1 | yes | section-0 |
| c9678fa7 | 2021-07-18 | 1 | yes | section-0 |
| 4628445d | 2021-07-19 | 1 | yes | section-0 |
| 2954aae1 | 2021-08-02 | 1 | no | section-0 |
| 470efdcf | 2021-08-03 | 1 | yes | section-0 |
| ed5301b9 | 2021-08-05 | 1 | no | section-0 |
| 637cc087 | 2021-08-06 | 1 | yes | section-0 |
| 6b1e0915 | 2021-08-08 | 1 | yes | section-0 |
| 5fa80d61 | 2021-08-09 | 1 | no | section-0 |
| 2a9e109c | 2021-08-11 | 1 | no | section-0 |
| 121be4d9 | 2021-08-12 | 1 | no | section-0 |
| 9523606a | 2021-08-14 | 1 | yes | section-0 |
| 8fb17461 | 2021-08-15 | 1 | yes | section-0 |
| 113b00bb | 2021-08-17 | 1 | no | section-0 |
| 3e7b2e51 | 2021-08-18 | 1 | yes | section-0 |
| 43fa08a4 | 2021-09-01 | 1 | yes | section-0 |
| 279a80cf | 2021-09-02 | 1 | no | section-0 |
| 943343c8 | 2021-09-04 | 1 | yes | section-0 |
| 5c7a5f2b | 2021-09-05 | 1 | no | section-0 |
| 95ed617d | 2021-09-07 | 1 | no | section-0 |
| 8b536979 | 2021-09-08 | 1 | no | section-0 |
| f9795c4e | 2021-09-10 | 1 | no | section-0 |
| 67581581 | 2021-09-11 | 1 | no | section-0 |
| b5c419dd | 2021-09-13 | 1 | yes | section-0 |
| af9418ec | 2021-09-14 | 1 | yes | section-0 |
| 4e19eea9 | 2021-09-16 | 1 | yes | section-0 |
| 6be2cf60 | 2021-09-17 | 1 | yes | section-0 |
| 7d93cc79 | 2021-10-01 | 1 | no | section-0 |
| 081dfed5 | 2021-10-02 | 1 | yes | section-0 |
| 0e576409 | 2021-10-04 | 1 | no | section-0 |
| 19e8a155 | 2021-10-05 | 1 | no | section-0 |
| c9cf3011 | 2021-10-07 | 1 | no | section-0 |
| 2bd66a14 | 2021-10-08 | 1 | yes | section-0 |
| 6456e108 | 2021-10-10 | 1 | yes | section-0 |
| 7aeaeb48 | 2021-10-11 | 1 | no | section-0 |
| 0d940063 | 2021-10-13 | 1 | no | section-0 |
| 1f205781 | 2021-10-14 | 1 | no | section-0 |
| 8adcea05 | 2021-10-16 | 1 | no | section-0 |
| 5ea0b8ef | 2021-10-17 | 1 | yes | section-0 |

## The caveats of the record

**The pipeline column is a different generation.** The column was computed on the asset-point chip of the whole-scene window, while this package recomputes the echo mask of the committed 7x7 window. The two therefore disagree, and they do: the committed column matches the recomputed `A` on **6 of 66** echo chips. Its pre-event median is 40.000 masked pixels against 3.000 recomputed; the figure draws the recomputed `A`, because that is the quantity this record's mask rule defines, and pins the committed column beside it rather than plotting it.

**A control with no variance.** `coherence` is 0.999285714286 on all 120 chips — 1 distinct value — so it separates nothing: it is reported, never drawn. All 8 co-variate controls of the channel table (`wind_speed_ms`, `temperature_c`, `brightness_ratio`, `intensity`, `coherence`, `coherence_gamma2`, `cumulative_mm`, `structural_frequency_hz`) are pinned in the JSON, and one of them (`structural_frequency_hz`) is committed for only 0 of the 120 chips, so it is not drawn either.

## The pins and the provenance

6541 pin comparisons against five committed layers (the channel table's meta, the mask-cache manifest, the committed mask cache, the paper result JSON and the analysis reference): **6533 exact**, **8 within 1e-12**, **0 failed**. The recomputed number is in `recomputed`, the committed one in `committed`, and no artifact is written before they agree.

| committed input | sha256 |
| --- | --- |
| `data/haltern/reference/haltern_bursts.json` | `74a9e1fd530cbe8f111fd6af0118d0337e35cad69ece82fa16ef858c8cee5c94` |
| `data/haltern/haltern_ocv_channels.csv` | `a82d716056ef7fa468f8b65d30302f9a03db6b1faba518de1c9b0296d52fc883` |
| `data/haltern/haltern_windows_mask_cache.txt` | `70031d3fe516a5cf982cc968e5ef1a18a4399a6b496b7b6e5d6a795d25ad08da` |
| `data/haltern/haltern_windows_mask_cache.manifest.json` | `1a679ead9ba43ab57074fcb88859f74310aec881425ff7e6f28f71a8cf78b2ba` |
| `data/haltern/haltern_measurements_full.txt` | `d179ff230d42ff4db7bd9fe7493353108e048501d37dde923dc747da0ccb115f` |
| `data/haltern/reference/haltern_geometry.json` | `fd71b4b3d9b3ea3162ed1c30312a876094a4fc20b4e22ca4b743e4a39eff9638` |
| `data/haltern/reference/haltern_osm_node.json` | `a37502af41bca70d5bf8e4bf3be93be9bc3a32ef9ccc046b7ae4033d28da617e` |
| `data/haltern/haltern_segments.txt` | `9052f8abb8610c69ae4784673839ca39805f087354700e845acffe32355f275e` |
| `data/haltern/haltern_windows_full.txt` | `e837309d7e7545b04829aad62d28c60cfd49bdbdb45a01aa196d85cd32a31dcd` |

Figures: `figures/haltern/fig_haltern_ocv_series.png` (600 dpi). JSON: `data/haltern/fig_haltern_ocv_series.json`.
