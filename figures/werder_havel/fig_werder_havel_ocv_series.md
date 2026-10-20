# fig_werder_havel_ocv_series

**Site** Havelbrücke Werder (A 10) (`40f420c6-c38e-5eab-8f2b-ac00d66707d0`) - Werder/Havel, the A 10 crossing of the Havel. **Axis** 2024-04-01 -> 2025-10-05 (553 d, 241 acquisitions, one per day). **Stretches** pre 31 (2024-04-01 .. 2024-06-29, all `W2_bestand`, all S1A), the documented repair 2024-07-01 .. 2025-09-05 (432 d drawn, 189 acquisitions, all `W3_sanierung`), post 21 (2025-09-05 .. 2025-10-05, `W4_post` truncated). **Plate** six channels x 241 acquisitions, 1438 drawn values, 183 mm, 600 dpi.

```
channels = [coherence_gamma2, phase_scatterer_count, phase_coherence,
            intensity, brightness_ratio, coherence]      the export's own columns
x        = the acquisition day (2024-04-01 -> 2025-10-05), one dot per
           acquisition; the four pass geometries alternate, no day carries two
y        = value / that channel's own pre-window median (1.0 = the median of
           the 31 pre acquisitions); all six panels logarithmic
base     = coherence_gamma2 0.000128534 / phase_scatterer_count 324 / phase_coherence 0.0124072
           intensity 12237.6 / brightness_ratio 0.94795 / coherence 0.0546966
band     = 0.002 .. 22: every drawn ordinate and every
           drawn level stays inside it (pinned, not a rule of thumb)
no line  = the mean step is 2.30 d over four interleaved geometries, so a
           line would draw the interleave and not the record
event    = the documented repair, drawn as a band and never as a level
```

## The axis

| stretch | first | last | n | days | committed stratum | that stratum's committed span | stratum truncated? | platforms | geometries |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| the three months before the documented start | 2024-04-01 | 2024-06-29 | 31 | 90 | W2_bestand | 2022-07-02 .. 2024-06-29 | yes | S1A 31 | asc_146_iw1 8, asc_44_iw3 8, desc_168_iw1 8, desc_95_iw3 7 |
| the band (the documented repair), 2024-07-01 .. 2025-09-05 | 2024-07-03 | 2025-09-04 | 189 | 429 | W3_sanierung | 2024-07-03 .. 2025-09-04 | no | S1A 143, S1C 46 | asc_146_iw1 49, asc_44_iw3 42, desc_168_iw1 49, desc_95_iw3 49 |
| the month after the documented finish | 2025-09-05 | 2025-10-05 | 21 | 31 | W4_post | 2025-09-05 .. 2025-12-04 | yes | S1A 10, S1C 11 | asc_146_iw1 6, asc_44_iw3 5, desc_168_iw1 5, desc_95_iw3 5 |

The record itself is 991 rows x 28 columns over 491 acquisition dates (2022-07-02 .. 2025-12-04), of which this plate draws the 241 between 2024-04-01 and 2025-10-05 (210 earlier, 40 later are not drawn). The chain inside the axis is strictly interleaved - the four pass geometries cycle, 0 days carry two acquisitions - and the three stretches cut on the record's own strata, not on the calendar: the documented start 2024-07-01 falls between `W2_bestand` and `W3_sanierung` and the documented finish 2025-09-05 on the first `W4_post` acquisition (Autobahn GmbH via havel_werder_analysis/havel_werder_asset_hierarchy.md (Step-0 A0 document) and havel_werder_asset.json). Only the band covers its stratum end to end; the pre and post stretches are marked truncated in the JSON because both strata run past the axis this figure documents.

## The six panels

Each panel draws one column of the export, and every ordinate is that column over the **pre window's own median** of it - re-derived here from the 31 pre acquisitions, never copied. Six bases: `coherence_gamma2` 0.000128534, `phase_scatterer_count` 324, `phase_coherence` 0.0124072, `intensity` 12237.6, `brightness_ratio` 0.94795, `coherence` 0.0546966. The pre column of the table below is therefore a distribution of ratios about 1.0, and it is the only window in which 1.0 means anything.

| panel | column | pre median (raw, the base) | n drawn | n empty | pre min | pre median | pre max | repair median | post median | post/pre (pooled) | pre-median spread over the four geometries |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| a | coherence_gamma2 | 0.000128534 | 241 | 0 | 0.085 | 1.000 | 8.429 | 1.312 | 1.182 | 1.182 | 10.882 |
| b | phase_scatterer_count | 324.000 | 241 | 0 | 0.645 | 1.000 | 1.133 | 0.935 | 0.926 | 0.926 | 1.216 |
| c | phase_coherence | 0.012 | 241 | 0 | 0.167 | 1.000 | 2.224 | 1.092 | 1.045 | 1.045 | 1.855 |
| d | intensity | 12237.622 | 241 | 0 | 0.364 | 1.000 | 1.758 | 0.994 | 0.964 | 0.964 | 3.279 |
| e | brightness_ratio | 0.948 | 241 | 0 | 0.723 | 1.000 | 1.299 | 0.991 | 0.993 | 0.993 | 1.069 |
| f | coherence | 0.055 | 233 | 8 | 0.197 | 1.000 | 6.544 | 0.783 | 0.838 | 0.838 | 2.453 |

All six channels are strictly positive on every drawn row (a squared coherence, a count, a ratio, a power), so every panel is logarithmic and no ordinate of this plate can reach 0; the drawn extent of the 1438 values is pinned as 0.002 .. 22 of their own base and the ordinate ticks are 0.003, 0.01, 0.03, 0.1, 0.3, 0.5, 1, 3, 10, 30.

The **pooled post/pre ratio** in the last-but-one column is the one number this plate refuses to read as a level. On `coherence_gamma2` it is 1.182, and the four geometries disagree about it in **sign**: 0.94 (asc_146_iw1), 1.97 (desc_168_iw1), 0.58 (asc_44_iw3), 0.54 (desc_95_iw3). Their own pre medians already span 10.9x inside the pre window alone, so the pooled number is a **mixture of geometries**, not a level of the bridge: a ratio of a median of four populations is not the median of their four ratios, and with sign disagreement it is not even the direction. `phase_scatterer_count` and `brightness_ratio` are the two channels whose four geometries agree (their pre medians span 1.069x and 1.216x respectively), and even there the direction is per geometry: see the table in the next section.

No trend line, no fit and no test is drawn: on a sawtooth that changes its level by 10.9x with the pass geometry, any line through the dots would measure the orbits. The honest per-geometry numbers are in the JSON (`geometries`), and the reference-chain reading is its own section below.

One caveat on the panel pair. `brightness_ratio` is `intensity` divided by a **per-geometry constant** - 20286.607 (asc_146_iw1), 14510.520 (desc_168_iw1), 6495.294 (asc_44_iw3), 9935.879 (desc_95_iw3) - so panels (d) and (e) are one channel rescaled and not two measurements (the constant holds to 1e-9 relative on all 241 rows). Both are drawn because the export carries both columns and this figure does not choose for the reader; the JSON says so, and so does the self-check.

## The four pass geometries

| geometry | orbit | sub-swath | platforms | n pre | n repair | n post | first | last | own pre median of gamma2 (ratio to the pooled base) | own post/pre of gamma2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| asc_146_iw1 | ASCENDING | IW1 | S1A 47, S1C 16 | 8 | 49 | 6 | 2024-04-01 | 2025-10-05 | 1.348 | 0.938 |
| desc_168_iw1 | DESCENDING | IW1 | S1A 46, S1C 16 | 8 | 49 | 5 | 2024-04-03 | 2025-10-01 | 0.600 | 1.969 |
| asc_44_iw3 | ASCENDING | IW3 | S1A 46, S1C 9 | 8 | 42 | 5 | 2024-04-06 | 2025-10-04 | 0.405 | 0.575 |
| desc_95_iw3 | DESCENDING | IW3 | S1A 45, S1C 16 | 7 | 49 | 5 | 2024-04-10 | 2025-10-02 | 4.411 | 0.541 |

The plate draws, inside the pre window only, **each geometry's own pre median** of that channel as a solid rule in that geometry's own colour, and the post window's pooled median as one dashed rule; the two are not the same quantity and the dashed rule is the pooled one, because a per-geometry post median over 5-6 acquisitions per geometry would be a number this figure cannot defend. A dot is read against its own geometry's rule, never across the sawtooth.

| channel | asc_146_iw1 | desc_168_iw1 | asc_44_iw3 | desc_95_iw3 |
| --- | --- | --- | --- | --- |
| coherence_gamma2 | 0.938 | 1.969 | 0.575 | 0.541 |
| phase_scatterer_count | 0.995 | 0.993 | 0.294 | 0.831 |
| phase_coherence | 1.167 | 0.772 | 0.858 | 1.337 |
| intensity | 0.921 | 1.091 | 1.859 | 0.964 |
| brightness_ratio | 0.921 | 1.091 | 1.859 | 0.964 |
| coherence | 1.985 | 2.093 | 0.119 | 0.528 |

That table is the honest reading of the plate, and it is why the pooled ratio in the previous section is quarantined: the same four post weeks move `coherence_gamma2` **up** on desc_168_iw1 and **down** on asc_146_iw1, asc_44_iw3, desc_95_iw3, and the same is true of `phase_coherence`, `intensity` and `coherence`. Only `phase_scatterer_count` agrees in sign across all four geometries, and there the drop is not uniform either (1.00 (asc_146_iw1), 0.99 (desc_168_iw1), 0.29 (asc_44_iw3), 0.83 (desc_95_iw3)).

## The coverage change inside the band

The axis is not uniform, and this is drawn as a **dotted vertical and never as a level**: S1C joins the constellation on 2025-03-30 (11 acquisitions that month, then 18 and 19 in the two months after), so 57 of the 241 drawn acquisitions are S1C and the axis' monthly count goes from 9-11 to 18-20. The mechanism is the revisit: asc_146_iw1 [12] -> [6], desc_168_iw1 [12] -> [6], asc_44_iw3 [12] -> [6, 12], desc_95_iw3 [12, 24] -> [6] days. Every geometry gains a 6-day step; three of the four drop the 12-day one, `asc_44_iw3` keeps it alongside and `desc_95_iw3` had one 24-day miss before. Two acquisitions more per month is not two more observations of the bridge - it is the same bridge sampled by two satellites, and the plate marks the change instead of smoothing it.

| month | acquisitions | platforms | geometries | first | last |
| --- | --- | --- | --- | --- | --- |
| 2024-04 | 11 | S1A 11 | asc_146_iw1 3, asc_44_iw3 3, desc_168_iw1 3, desc_95_iw3 2 | 2024-04-01 | 2024-04-30 |
| 2024-05 | 10 | S1A 10 | asc_146_iw1 3, asc_44_iw3 2, desc_168_iw1 2, desc_95_iw3 3 | 2024-05-04 | 2024-05-31 |
| 2024-06 | 10 | S1A 10 | asc_146_iw1 2, asc_44_iw3 3, desc_168_iw1 3, desc_95_iw3 2 | 2024-06-02 | 2024-06-29 |
| 2024-07 | 10 | S1A 10 | asc_146_iw1 3, asc_44_iw3 2, desc_168_iw1 2, desc_95_iw3 3 | 2024-07-03 | 2024-07-30 |
| 2024-08 | 10 | S1A 10 | asc_146_iw1 2, asc_44_iw3 3, desc_168_iw1 3, desc_95_iw3 2 | 2024-08-01 | 2024-08-28 |
| 2024-09 | 10 | S1A 10 | asc_146_iw1 3, asc_44_iw3 2, desc_168_iw1 3, desc_95_iw3 2 | 2024-09-04 | 2024-09-30 |
| 2024-10 | 10 | S1A 10 | asc_146_iw1 2, asc_44_iw3 3, desc_168_iw1 2, desc_95_iw3 3 | 2024-10-03 | 2024-10-31 |
| 2024-11 | 10 | S1A 10 | asc_146_iw1 3, asc_44_iw3 2, desc_168_iw1 3, desc_95_iw3 2 | 2024-11-03 | 2024-11-29 |
| 2024-12 | 10 | S1A 10 | asc_146_iw1 2, asc_44_iw3 3, desc_168_iw1 2, desc_95_iw3 3 | 2024-12-02 | 2024-12-30 |
| 2025-01 | 11 | S1A 11 | asc_146_iw1 3, asc_44_iw3 3, desc_168_iw1 3, desc_95_iw3 2 | 2025-01-02 | 2025-01-31 |
| 2025-02 | 9 | S1A 9 | asc_146_iw1 2, asc_44_iw3 2, desc_168_iw1 2, desc_95_iw3 3 | 2025-02-04 | 2025-02-28 |
| 2025-03 | 11 | S1A 10, S1C 1 | asc_146_iw1 3, asc_44_iw3 2, desc_168_iw1 3, desc_95_iw3 3 | 2025-03-03 | 2025-03-30 |
| 2025-04 | 18 | S1A 10, S1C 8 | asc_146_iw1 5, asc_44_iw3 3, desc_168_iw1 5, desc_95_iw3 5 | 2025-04-01 | 2025-04-29 |
| 2025-05 | 19 | S1A 11, S1C 8 | asc_146_iw1 5, asc_44_iw3 4, desc_168_iw1 5, desc_95_iw3 5 | 2025-05-02 | 2025-05-31 |
| 2025-06 | 19 | S1A 9, S1C 10 | asc_146_iw1 5, asc_44_iw3 4, desc_168_iw1 5, desc_95_iw3 5 | 2025-06-01 | 2025-06-28 |
| 2025-07 | 19 | S1A 11, S1C 8 | asc_146_iw1 6, asc_44_iw3 3, desc_168_iw1 5, desc_95_iw3 5 | 2025-07-01 | 2025-07-31 |
| 2025-08 | 20 | S1A 10, S1C 10 | asc_146_iw1 5, asc_44_iw3 5, desc_168_iw1 5, desc_95_iw3 5 | 2025-08-02 | 2025-08-30 |
| 2025-09 | 20 | S1A 10, S1C 10 | asc_146_iw1 5, asc_44_iw3 5, desc_168_iw1 5, desc_95_iw3 5 | 2025-09-01 | 2025-09-29 |
| 2025-10 | 4 | S1A 2, S1C 2 | asc_146_iw1 1, asc_44_iw3 1, desc_168_iw1 1, desc_95_iw3 1 | 2025-10-01 | 2025-10-05 |

## The steps of the chain

| step | how many | what it is | the steps |
| --- | --- | --- | --- |
| 1 | 57 | the cadence | 2025-03-29 -> 2025-03-30, 2025-04-01 -> 2025-04-02, 2025-04-04 -> 2025-04-05, 2025-04-10 -> 2025-04-11, 2025-04-13 -> 2025-04-14, 2025-04-16 -> 2025-04-17, ... and 51 more - all 57 are in the JSON (`steps`) |
| 2 | 87 | the cadence | 2024-04-01 -> 2024-04-03, 2024-04-13 -> 2024-04-15, 2024-04-25 -> 2024-04-27, 2024-05-07 -> 2024-05-09, 2024-05-19 -> 2024-05-21, 2024-05-31 -> 2024-06-02, ... and 81 more - all 87 are in the JSON (`steps`) |
| 3 | 66 | the cadence | 2024-04-03 -> 2024-04-06, 2024-04-10 -> 2024-04-13, 2024-04-15 -> 2024-04-18, 2024-04-22 -> 2024-04-25, 2024-04-27 -> 2024-04-30, 2024-05-04 -> 2024-05-07, ... and 60 more - all 66 are in the JSON (`steps`) |
| 4 | 29 | the cadence | 2024-04-06 -> 2024-04-10, 2024-04-18 -> 2024-04-22, 2024-04-30 -> 2024-05-04, 2024-05-12 -> 2024-05-16, 2024-05-24 -> 2024-05-28, 2024-06-05 -> 2024-06-09, ... and 23 more - all 29 are in the JSON (`steps`) |
| 7 | 1 | a hole, drawn as a band | 2024-08-28 -> 2024-09-04 |

No line crosses the plate and no line crosses a hole: the mean step is 2.30 days over four interleaved geometries, and the one step longer than the 4-day cadence (2024-08-28 -> 2024-09-04, 7 days) is drawn as a grey band and never interpolated over. The table names the first 6 pairs of each step length; all 240 steps of the chain are in the JSON (`steps`).

## The reference chains of the two IW1 geometries

Two of the four geometries carry a target window **and** its own reference window on the very same dates - 125 of the axis' 241 acquisitions ({'asc_146_iw1': 63, 'desc_168_iw1': 62}; pre 16, repair 98, post 11). That is the one comparison of this record that does **not** mix pass geometries, because both windows are measured on the same pass at different places. It is not drawn - it would need a second ordinate - and it is not a test; it is reported:

| stretch | n dates | target gamma2 (median) | reference gamma2 | outside-reference gamma2 | target / reference | target / outside-reference | target / reference (scatterers) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| pre | 16 | 0.00014977 | 0.000107938 | 0.000123909 | 1.323 | 1.397 | 0.899 |
| repair | 98 | 0.00013998 | 0.000100074 | 0.000131438 | 1.226 | 1.254 | 0.975 |
| post | 11 | 0.000151941 | 0.00019889 | 4.54517e-05 | 0.666 | 4.309 | 0.974 |

In the post window that reading is **0.666** - the target window's gamma2 over its own same-pass reference window, i.e. *below* the reference - while the pooled post/pre ratio of the same column is 1.182. In the pre window the same reading is 1.323. Both numbers are in the JSON and neither is drawn: the figure's plate carries the six channels, and a reader who wants the reference reading has it here rather than reconstructed from a plot that does not show it.

## The holes and the artefacts of the record

An empty field in the export is a **missing value** and this figure never turns one into a 0. Inside the axis:

| channel | empty rows of 241 | which dates |
| --- | --- | --- |
| coherence_gamma2 | 0 | - |
| phase_scatterer_count | 0 | - |
| phase_coherence | 0 | - |
| intensity | 0 | - |
| brightness_ratio | 0 | - |
| coherence | 8 | 2024-07-03, 2024-07-06, 2024-07-08, 2024-07-11, 2025-09-05, 2025-09-07, 2025-09-08, 2025-09-10 |

The 8 empty `coherence` rows are the first four acquisitions of the repair stratum and the first four of the post stratum - the two places where the chain changes state - and they are drawn as ticks on the panel floor, never as zeros. Across the whole record (991 rows) two columns are empty on **every** row: `displacement_los_m`, `phase_snr_db` - they are not drawn at all. `coherence` is empty on 24 of the record's 991 rows. `coherence_masked_pixels` is 640 = 10 % of the 6400-px window in every row of the record: a definition (the brightest decile), not a measurement, and never drawn. The record's window note: one 80x80 complex window per acquisition at the site point; the export committed no latitude/longitude, no segment_index and no burst id (the anchor is the OSM bridge outline, see werder_havel_aoi.py), and a window carries no pixel-spacing spec.

The four co-variate controls of the record are carried on every drawn row and none of them is plotted (a second ordinate the plate does not have):

| control | median pre | median repair | median post | axis min | axis max |
| --- | --- | --- | --- | --- | --- |
| temperature [degC] | 15.600 | 14.300 | 14.600 | -5.000 | 25.800 |
| precipitation [mm] | 0.400 | 0.100 | 0.300 | 0.000 | 32.800 |
| shift_x [m] | 0.032 | 0.002 | -0.077 | -0.776 | 0.820 |
| shift_y [m] | -0.072 | -0.080 | -0.199 | -174.150 | 2.106 |

And the record carries 14 rows whose `shift_y` is about -174 m where every other row of the record stays within 2.21 m - all of them `asc_44_iw3`, all S1C, 9 of them inside this axis. They are pinned in the JSON, never drawn and never removed:

| date | geometry | platform | shift_x [m] | shift_y [m] |
| --- | --- | --- | --- | --- |
| 2025-05-13 | asc_44_iw3 | S1C | 0.527 | -174.128 |
| 2025-06-06 | asc_44_iw3 | S1C | 0.335 | -173.652 |
| 2025-06-18 | asc_44_iw3 | S1C | 0.366 | -174.087 |
| 2025-08-05 | asc_44_iw3 | S1C | 0.081 | -173.821 |
| 2025-08-17 | asc_44_iw3 | S1C | 0.099 | -173.462 |
| 2025-08-29 | asc_44_iw3 | S1C | -0.093 | -173.939 |
| 2025-09-10 | asc_44_iw3 | S1C | 0.184 | -173.819 |
| 2025-09-22 | asc_44_iw3 | S1C | -0.096 | -174.150 |
| 2025-10-04 | asc_44_iw3 | S1C | -0.030 | -173.988 |
| 2025-10-16 | asc_44_iw3 | S1C | -0.122 | -173.216 |
| 2025-10-28 | asc_44_iw3 | S1C | 0.033 | -173.814 |
| 2025-11-09 | asc_44_iw3 | S1C | 0.124 | -173.556 |
| 2025-11-21 | asc_44_iw3 | S1C | 0.037 | -173.567 |
| 2025-12-03 | asc_44_iw3 | S1C | 0.654 | -174.014 |
## What this figure asserts

19 of 19 self-checks pass, each pinned against the layer it belongs to, plus 66 numeric pins (41 exact, 25 within tolerance, 0 failed). The list is in the JSON under `self_check`; a failed check aborts the run before anything is written.

- **ok** - the committed measurement export re-verifies field by field: ocv_masks.check(): no problem
- **ok** - the axis is 241 acquisitions over 553 days, one acquisition per day: 241 rows on 241 distinct days, 2024-04-01 .. 2025-10-05 (553 d, mean step 2.30 d)
- **ok** - the three stretches partition the axis, in order and without overlap: pre 31 + repair 189 + post 21 = 241
- **ok** - every drawn acquisition is inside its own committed stratum: pre = 2024-04-01 .. 2024-06-29; repair = 2024-07-03 .. 2025-09-04; post = 2025-09-05 .. 2025-10-05 (only the repair covers its stratum end to end, so the other two are marked truncated in the JSON)
- **ok** - the repair band is drawn on the documented span and the other two on the axis: band 2024-07-01 .. 2025-09-05 = 432 d
- **ok** - the base of every channel is the pre window's own median, over 31 values: six bases: coherence_gamma2 0.000128534, phase_scatterer_count 324, phase_coherence 0.0124072, intensity 12237.6, brightness_ratio 0.94795, coherence 0.0546966
- **ok** - every drawn ordinate is positive and inside the pinned band: 1438 drawn values, 0.002072 .. 21.72 of their own base, band 0.002 .. 22
- **ok** - the four geometries disagree on gamma2's post/pre sign and on its pre level: post/pre 0.938, 1.969, 0.575, 0.541; pre medians span 10.9x
- **ok** - brightness_ratio is the most geometry-stable channel and phase_scatterer_count the second: pre-median spread: brightness_ratio 1.069x, phase_scatterer_count 1.216x, coherence_gamma2 10.9x
- **ok** - brightness_ratio is intensity over a per-geometry constant: constants 20286.607 / 14510.521 / 6495.294 / 9935.879, worst relative spread 1.1e-09: panels d and e are one channel rescaled, not two measurements
- **ok** - coherence is empty on the first four rows of the repair and of the post stretch, and on nothing else: 8 empty coherence rows: 2024-07-03, 2024-07-06, 2024-07-08, 2024-07-11, 2025-09-05, 2025-09-07, 2025-09-08, 2025-09-10; the other five channels are empty on 0 of 241
- **ok** - the two columns empty on every row of the record are not drawn: empty on all 991 rows: displacement_los_m, phase_snr_db; coherence is empty on 24 of the record's 991 rows (8 of them inside this axis); coherence_masked_pixels = 640 = 10 % of the 6400 px window, one value
- **ok** - the chain has exactly one step longer than the cadence, and it is drawn: steps 2024-08-28 -> 2024-09-04 (7 d); every other step is 1..4 days
- **ok** - the coverage change (S1C joining) is pinned, drawn and never read as a level: S1C from 2025-03-30 (11 acquisitions that month), 57 of the 241 drawn rows; per-geometry revisit steps asc_146_iw1 [12] -> [6], desc_168_iw1 [12] -> [6], asc_44_iw3 [12] -> [6, 12], desc_95_iw3 [12, 24] -> [6]; no day carries two acquisitions
- **ok** - all four co-variate controls are carried on every drawn row and none is drawn: temperature_c, precipitation_mm, shift_x, shift_y: 964 of 964 cells filled, reported in the JSON, not plotted
- **ok** - the shift_y artefacts are reported, not removed: 14 of the record's 991 rows carry shift_y ~ -174 m (-174.2 .. -173.2), all asc_44_iw3, all S1C, 9 of them inside this axis; the other 477 rows stay within 2.21 m
- **ok** - the two IW1 geometries carry their own reference chain on the same dates: 125 of the axis' 241 acquisitions carry a target and a reference window ({'asc_146_iw1': 63, 'desc_168_iw1': 62}; pre 16, repair 98, post 11), the other 143 do not
- **ok** - the same-pass reference reading disagrees with the pooled post/pre ratio: target/reference gamma2 in the post window 0.666 (pre 1.323, outside-reference in the post window 4.31), against the pooled post/pre 1.182 of panel a - neither is drawn and neither is a test, both are in the JSON
- **ok** - every pinned number agrees with the layer it is pinned against: 66 comparisons (41 exact, 25 within tolerance): 0 failed, max |deviation| 0

## Provenance

- the committed export: `data/werder_havel/werder_havel_measurements_full.txt`, sha256 `b1f9fcc02e4645e3416e4acc5e695ef8faff9ce7d78199249337240b99b25d98`, 991 rows x 28 columns, 491 acquisition dates (2022-07-02 .. 2025-12-04), roles reference 250, reference_outside 250, target 491; platforms S1A, S1C; state windows W2_bestand 241, W3_sanierung 189, W4_post 61.
- the committed segment list: `data/werder_havel/werder_havel_segments.txt`, sha256 `962a974c6651799aed7a89fad3b29804e6ff4ca8ebdfaf7569b63431996ddacf`.
- the committed AOI reference: `data/werder_havel/reference/werder_havel_aoi.json`, re-derived field by field and pinned; the documented change is `2024-07-01` .. `2025-09-05` (Autobahn GmbH via havel_werder_analysis/havel_werder_asset_hierarchy.md (Step-0 A0 document) and havel_werder_asset.json).
- the record's own window note: one 80x80 complex window per acquisition at the site point; the export committed no latitude/longitude, no segment_index and no burst id (the anchor is the OSM bridge outline, see werder_havel_aoi.py), and a window carries no pixel-spacing spec
- the plate: `figures/werder_havel/fig_werder_havel_ocv_series.png` (183 mm, 600 dpi, 3 x 2 panels); the numbers: `data/werder_havel/fig_werder_havel_ocv_series.json`.
- this run: 0.0 s, stdlib + numpy + matplotlib only, offline, no cache read and none written.
