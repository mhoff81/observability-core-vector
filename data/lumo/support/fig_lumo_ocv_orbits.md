# LUMO OCV time series by dimension and orbit

## 1. Question and layout

How do the six OCV components of the LUMO lattice tower move through the documented damage states, and do the two Sentinel-1 pass geometries tell the same story? The figure draws one panel per **dimension x orbit**: six rows in the formal OCV order (`gamma2`, `P`, `D`, `F`, `A`, `S`) and two columns (ASC / DESC), i.e. twelve panels.

The x coordinate is the **actual local overpass instant**, not a calendar date. Coherence rows carry no timestamp, so the acquisition time is joined from the `record_class="acquisition"` rows by `burst_id` and converted from UTC to `Europe/Berlin` with `zoneinfo`. The background bands follow the **13 consecutive `damage_label` runs** of the time-sorted cohort, not fixed windows; the twelve dashed lines mark the run transitions.

## 2. Data, join and time base

- source: `data/lumo/support/lumo_ocv_channels.csv` (178 coherence rows, 482 acquisition rows)
- committed metadata: `data/lumo/support/lumo_ocv_channels_meta.json`
- source channel table: `lumo_ocv_channels.csv`
- reference: `data/lumo/reference/lumo_tower_coherence_states.json`
- join: **178 / 178** coherence rows matched one-to-one by `burst_id`; 304 acquisition rows have no coherence row and are not plotted
- dates: **178** unique; bursts: **178** unique
- orbits: **57 ASC**, **121 DESC**
- states: **healthy 121**, **DAM3 25**, **DAM4 16**, **DAM6 16**

Local overpass clock and UTC offsets (DST-aware):

| orbit | earliest | latest | UTC offsets |
| --- | --- | --- | --- |
| ASC | 2020-12-04 17:09:00+01:00 | 2021-07-20 19:09:04+02:00 | +01:00: 23, +02:00: 34 |
| DESC | 2020-08-15 05:33:56+02:00 | 2021-07-22 07:42:12+02:00 | +01:00: 52, +02:00: 69 |

Sub-swath provenance of the plotted rows:

| orbit | IW1 | IW2 |
| --- | --- | --- |
| ASC | 1 | 56 |
| DESC | 108 | 13 |

## 3. The 13 consecutive state runs

A run is a maximal block of equal `damage_label` after sorting by acquisition instant. Bands span the midpoint between the last overpass of one run and the first of the next; the first band starts 10 days before the first overpass and the last ends 10 days after the last one.

| run | state | n | ASC | DESC | first | last |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | healthy | 36 | 12 | 24 | 2020-08-02 | 2020-10-11 |
| 2 | DAM6 | 8 | 2 | 6 | 2020-10-13 | 2020-10-26 |
| 3 | healthy | 6 | 2 | 4 | 2020-10-29 | 2020-11-07 |
| 4 | DAM4 | 8 | 3 | 5 | 2020-11-10 | 2020-11-24 |
| 5 | healthy | 54 | 17 | 37 | 2020-11-25 | 2021-03-13 |
| 6 | DAM3 | 17 | 5 | 12 | 2021-03-18 | 2021-04-18 |
| 7 | healthy | 7 | 3 | 4 | 2021-04-21 | 2021-05-03 |
| 8 | DAM6 | 8 | 2 | 6 | 2021-05-05 | 2021-05-18 |
| 9 | healthy | 4 | 2 | 2 | 2021-05-21 | 2021-05-27 |
| 10 | DAM4 | 8 | 2 | 6 | 2021-05-29 | 2021-06-14 |
| 11 | healthy | 5 | 1 | 4 | 2021-06-16 | 2021-06-23 |
| 12 | DAM3 | 8 | 3 | 5 | 2021-06-26 | 2021-07-11 |
| 13 | healthy | 9 | 3 | 6 | 2021-07-14 | 2021-07-29 |

## 4. Per-state medians

Median of each channel inside each state, with the number of finite overpasses in parentheses. The log-scaled channels are strictly positive in every state x orbit cell. The single missing `P` (2021-07-14 ASC) is drawn as a gap and is not interpolated.

| orbit | dimension | healthy | DAM3 | DAM4 | DAM6 |
| --- | --- | --- | --- | --- | --- |
| ASC | gamma2 | 0.1241 (n=40) | 0.1289 (n=8) | 0.0683 (n=5) | 0.3738 (n=4) |
| ASC | P | 0.0935 (n=40) | 0.0622 (n=8) | 0.0791 (n=5) | 0.1001 (n=4) |
| ASC | D | 0.0042 (n=40) | 0.0045 (n=8) | 0.0050 (n=5) | 0.5019 (n=4) |
| ASC | F | 4.0000 (n=40) | 3.0000 (n=8) | 4.0000 (n=5) | 1.5000 (n=4) |
| ASC | A | 10.0000 (n=40) | 12.5000 (n=8) | 11.0000 (n=5) | 5.5000 (n=4) |
| ASC | S | 110.9221 (n=40) | 138.7067 (n=8) | 108.2076 (n=5) | 73.0865 (n=4) |
| DESC | gamma2 | 0.6509 (n=81) | 0.6086 (n=17) | 0.7066 (n=11) | 0.4556 (n=12) |
| DESC | P | 0.1917 (n=81) | 0.1333 (n=17) | 0.2400 (n=11) | 0.1869 (n=12) |
| DESC | D | 0.0102 (n=81) | 0.0107 (n=17) | 1.0000 (n=11) | 0.0071 (n=12) |
| DESC | F | 2.0000 (n=81) | 2.0000 (n=17) | 1.0000 (n=11) | 2.0000 (n=12) |
| DESC | A | 4.0000 (n=81) | 4.0000 (n=17) | 4.0000 (n=11) | 4.0000 (n=12) |
| DESC | S | 14.2215 (n=81) | 13.2877 (n=17) | 0.7071 (n=11) | 28.4463 (n=12) |

Means, minima and maxima are in the result JSON under `summary.<orbit>.<dimension>.<state>`.

## 5. Cross-checks

The script recomputes every plotted point from the table it reads (`data/lumo/support/lumo_ocv_channels.csv`) and checks it against the metadata, the paper figure JSON and the committed tower reference: **123 checks, 0 failures**.

**Failures: none.** The checks cover row and date counts, the 482-row acquisition join, the 13-run sequence and lengths, sub-swath counts, per-state/per-orbit counts, row-by-row `burst_id` / date / orbit / label parity with the tower reference, `A == n_masked`, `A <= chip`, strict positivity on the log axes, local-time ranges and UTC offsets, and the ring definition pins.

| check group | checks |
| --- | --- |
| A_equals_n_masked | 1 |
| data | 11 |
| join | 13 |
| local | 5 |
| log_positive | 32 |
| meta | 15 |
| n_masked_positive | 1 |
| paper | 8 |
| reference | 12 |
| ring | 6 |
| runs | 10 |
| series | 8 |
| subswath | 1 |

**Superseded by the reference rule (24 pins, not failures).** This table measures with a reference ring that is disjoint from the support, so the pins that describe the committed rule -- the median inside the search window, the 400 x 11 chip, channel completeness and the `n_masked` / `gamma2` values the paper reported -- cannot hold here. They are kept with both sides in `superseded` in the JSON:

| check | ring value | committed value |
| --- | --- | --- |
| channels.finite | {'gamma2': 178, 'P': 178, 'D': 178, 'F': 178, 'A': 178, 'S': 178} | {'gamma2': 178, 'P': 177, 'D': 178, 'F': 178, 'A': 178, 'S': 178} |
| channels.n_rows_with_missing | 0 | 1 |
| channels.missing_dims | [] | ['P'] |
| channels.complete_6d | 178 | 177 |
| channels.complete_6d_by_state | {'healthy': 121, 'DAM3': 25, 'DAM4': 16, 'DAM6': 16} | {'DAM3': 25, 'DAM4': 16, 'DAM6': 16, 'healthy': 120} |
| P.missing.n | 0 | 1 |
| P.missing.burst_id | [] | ['3801dfa4-5872-44fe-8538-84c178cac83f'] |
| P.missing.date | [] | ['2021-07-14'] |
| P.missing.orbit | [] | ['ASCENDING'] |
| P.missing.label | [] | ['healthy'] |
| P.missing.subswath | [] | ['iw1'] |
| paper.n_by_state_complete_6d | {'DAM3': 25, 'DAM4': 16, 'DAM6': 16, 'healthy': 120} | {'healthy': 121, 'DAM3': 25, 'DAM4': 16, 'DAM6': 16} |
| paper.n_masked_values.healthy | [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 | [2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 7, 7, 7 |
| paper.n_masked_median.healthy | 4.0 | 15.0 |
| paper.n_masked_values.DAM3 | [2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 6, 6, 7, 8, 8, 9, 11, 12, 13, 13, 13, 15 | [2, 4, 5, 5, 7, 7, 9, 9, 10, 10, 12, 12, 13, 13, 14, 14, 14, 16, 17, 17, 27, 36, |
| paper.n_masked_median.DAM3 | 5.0 | 13.0 |
| paper.n_masked_values.DAM4 | [2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 8, 11, 11, 14, 16] | [2, 5, 5, 6, 8, 10, 10, 13, 15, 16, 20, 23, 28, 46, 81, 193] |
| paper.n_masked_median.DAM4 | 4.0 | 14.0 |
| paper.n_masked_values.DAM6 | [2, 2, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 6, 7, 8, 11] | [4, 4, 5, 7, 8, 11, 12, 12, 12, 12, 15, 19, 21, 117, 124, 218] |
| paper.n_masked_median.DAM6 | 4.0 | 12.0 |
| reference.gamma2 | [0.43104165509818504, 0.0016436554898093363, 0.04282334886686867, 0.011730188657 | [0.8624535129217228, 0.31696055541415974, 0.03550748548152263, 0.602804267273875 |
| reference.n_masked | [7, 78, 15, 9, 19, 20, 3, 117, 16, 35, 117, 9, 11, 10, 10, 2, 149, 18, 16, 13, 5 | [5, 10, 10, 3, 4, 14, 2, 4, 15, 3, 2, 10, 7, 12, 15, 25, 4, 4, 2, 7, 9, 2, 2, 14 |
| reference.per_state.n_masked.median | {'DAM0(healthy)': 4.0, 'DAM3': 5.0, 'DAM4': 4.0, 'DAM6': 4.0} | {'DAM0(healthy)': 15.0, 'DAM3': 13.0, 'DAM4': 14.0, 'DAM6': 12.0} |
| data.n_nonfinite | {'gamma2': 0, 'P': 0, 'D': 0, 'F': 0, 'A': 0, 'S': 0} | {'gamma2': 0, 'P': 1, 'D': 0, 'F': 0, 'A': 0, 'S': 0} |

## 6. Reproduce

```sh
cd <repo root>
python3 code/lumo/fig_lumo_ocv_orbits.py \
    --csv data/lumo/support/lumo_ocv_channels.csv \
    --meta data/lumo/support/lumo_ocv_channels_meta.json
```

Needs numpy and matplotlib. This report was generated with python 3.13.5, numpy 2.5.4, matplotlib 3.11.2. No output carries a timestamp, so repeated runs are identical.

## 7. Caveats

- **The 13 runs are one deployment, not 13 independent experiments.** States repeat (healthy seven times, each damage state twice), but the repeats are separated in time and season. The bands make the sequence visible; they do not turn it into a randomized treatment.
- **State and season are confounded.** The first healthy block covers late summer and autumn 2020, the long healthy block spans winter 2020/21, and the damage states are short visits in spring and summer 2021.
- **Orbit and sub-swath are not balanced.** ASC is almost entirely IW2 (56 of 57 rows), DESC mostly IW1 (108 of 121); the columns therefore differ in both viewing geometry and sub-swath.
- **The x axis is local wall time.** ASC overpasses are around 17–19 local, DESC around 05–08 local. UTC storage and the DST offset change are handled explicitly; a date-only axis would hide the two clocks.
- **Missing values are gaps.** Exactly one of 178 rows lacks `P` (2021-07-14 ASC, the one non-majority chip shape excluded from the persistence map). Nothing is interpolated.
- **`P`, `D`, `F` and `A` use log axes** because their ranges span orders of magnitude; `gamma2` and `S` are linear.

---

Generated by `code/lumo/fig_lumo_ocv_orbits.py` (123 cross-checks, 0 failures).
