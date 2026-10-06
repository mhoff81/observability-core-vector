# Carola Bruecke, Dresden — the echo-mask dimensions month by month

**Site** Carola Bruecke, Dresden · **event** 2024-09-11 · **chips** 210 (windows M-3..M+1 only)

```
x     = [gamma2, P, D, A, F, S]      echo-mask vector
mask  = intensity >= 0.30*peak AND >= 5.0*np.median(intensity), A >= 2
M-3   [2024-06-11, 2024-07-11)       pre-collapse (healthy)
M-2   [2024-07-11, 2024-08-11)       pre-collapse (healthy)
M-1   [2024-08-11, 2024-09-11)       pre-collapse (healthy)
M+1   [2024-09-11, 2024-10-11)       post-collapse
```

## The four windows

Rolling one-month buckets anchored at the collapse date (the event opens M+1, so the pre side is the three M-* windows and the post side is M+1 alone). This is a companion of `fig_carola_ocv_paper.py` (figures A–E): it adds one figure and does not touch the A–E pin contract.

| window | range | state | chips | first day | last day | days |
| --- | --- | --- | --- | --- | --- | --- |
| M-3 | [2024-06-11, 2024-07-11) | pre | 42 | 2024-06-12 | 2024-07-06 | 7 |
| M-2 | [2024-07-11, 2024-08-11) | pre | 52 | 2024-07-11 | 2024-08-08 | 8 |
| M-1 | [2024-08-11, 2024-09-11) | pre | 56 | 2024-08-11 | 2024-09-09 | 7 |
| M+1 | [2024-09-11, 2024-10-11) | post | 60 | 2024-09-13 | 2024-10-10 | 8 |

## Per-window medians of the six dimensions

| dim | pre-baseline | M-3 | M-2 | M-1 | M+1 |
| --- | --- | --- | --- | --- | --- |
| gamma2 | 0.08882 | 0.1578 | 0.07215 | 0.0365 | 0.06237 |
| P | 0.02677 | 0.02199 | 0.02563 | 0.0327 | 0.03182 |
| D | 0.005724 | 0.004717 | 0.005169 | 0.007286 | 0.005789 |
| A | 15.33 | 13 | 13 | 20 | 17 |
| F | 7.333 | 7 | 7 | 8 | 9 |
| S | 19.74 | 25 | 15.11 | 19.09 | 10.65 |

`pre-baseline` = mean of the three pre-window medians (the 1.0 line of panel f).

## What moved across the cut (M+1 vs M-1)

| dim | median M-1 | median M+1 | Cliff's delta | SDS | CI95 excludes 0 | CI95 |
| --- | --- | --- | --- | --- | --- | --- |
| gamma2 | 0.0365 | 0.06237 | 0.07798 | 0.7798 | no | [-0.15, 0.3036] |
| P | 0.0327 | 0.03182 | -0.1917 | 1.917 | no | [-0.4024, 0.02292] |
| D | 0.007286 | 0.005789 | -0.3554 | 3.554 | yes | [-0.5503, -0.1571] |
| A | 20 | 17 | 0.1821 | 1.821 | no | [-0.02827, 0.3839] |
| F | 8 | 9 | 0.1664 | 1.664 | no | [-0.04673, 0.3723] |
| S | 19.09 | 10.65 | 0.01786 | 0.1786 | no | [-0.2021, 0.2327] |

Cliff's delta > 0 means M+1 larger than M-1. The same comparison against the three pre windows pooled is in the JSON (`deltas[dim].vs_pre_pool`).

## Reading (descriptive — not a claim of damage)

The three pre windows are **not** stationary: `gamma2`'s median roughly halves twice over M-3..M-1 (0.158 → 0.072 → 0.037) and `S` swings 25 → 15 → 19, so the month-to-month spread is large even before the event. Against that noise the single adjacent-month contrast (M+1 vs M-1) is mostly weak: only `D` clears the bootstrap CI (median 0.0073 → 0.0058, delta −0.36), while `P` — the sharpest difference in the pooled A–E figures (delta −0.661) — moves here by only −0.19 with a CI that spans zero. Switching the baseline to the three pre windows pooled puts `D` back at −0.03 (CI includes 0) and lifts `A` (+0.29) and `F` (+0.24), so the sizes and signs depend on which pre baseline is chosen. The month axis therefore does **not** reproduce the pooled headline by itself; it narrows it and shows how much of the pre-collapse record already moves on its own. Both sides live on the *same* bridge, so season, weather, traffic and pipeline generation move with the collapse date as well.

```bash
python3 code/carola/carola_ocv_channels_csv.py       # build + verify the table
python3 code/carola/fig_carola_ocv_months.py          # this figure
```
