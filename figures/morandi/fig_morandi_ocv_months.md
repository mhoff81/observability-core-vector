# Ponte Morandi (Polcevera), Genova — the echo-mask dimensions month by month around the collapse

**Site** Ponte Morandi (Polcevera), Genova · **event** 2018-08-14 · **chips** 39 in the four one-month windows M-3 … M+1

```
M-3  [2018-05-14, 2018-06-14)   pre-collapse (healthy)
M-2  [2018-06-14, 2018-07-14)   pre-collapse (healthy)
M-1  [2018-07-14, 2018-08-14)   pre-collapse (healthy)
M+1  [2018-08-14, 2018-09-14)   post-collapse
```

## Windows

| window | interval | state | n | days |
| --- | --- | --- | --- | --- |
| M-3 | [2018-05-14, 2018-06-14) | pre | 10 | 2018-05-16 … 2018-06-11 |
| M-2 | [2018-06-14, 2018-07-14) | pre | 8 | 2018-06-17 … 2018-07-11 |
| M-1 | [2018-07-14, 2018-08-14) | pre | 10 | 2018-07-15 … 2018-08-10 |
| M+1 | [2018-08-14, 2018-09-14) | post | 11 | 2018-08-14 … 2018-09-13 |

## Median of each dimension per window

| dimension | M-3 | M-2 | M-1 | M+1 |
| --- | --- | --- | --- | --- |
| gamma2 | 0.08577 | 0.1082 | 0.258 | 0.1016 |
| P | 0.05981 | 0.07085 | 0.04918 | 0.07256 |
| D | 0.004879 | 0.003505 | 0.008443 | 0.005159 |
| A | 11.5 | 9.5 | 8.5 | 13 |
| F | 5 | 5 | 4 | 7 |
| S | 7.767 | 5.466 | 14.47 | 13.78 |

## M+1 against the pre side

| dimension | δ vs M-1 | CI≠0 vs M-1 | δ vs pooled pre | CI≠0 vs pooled pre |
| --- | --- | --- | --- | --- |
| gamma2 | -0.1091 | no | 0.05195 | no |
| P | -0.09091 | no | -0.1299 | no |
| D | -0.1182 | no | 0.06818 | no |
| A | 0.1455 | no | 0.1234 | no |
| F | 0.07273 | no | 0.09091 | no |
| S | -0.07273 | no | 0.01299 | no |

## What the figure shows

The first post month carries **11** chips against 28 before the event, so its intervals are wide and this is a descriptive companion, not a test. Of the six dimensions, **0** show a Cliff's-delta 95 % CI that excludes zero against M-1; the (P, D) plane (panel e) and the 6D fingerprint (panel f) show *when* the mask moved across the four windows.

Both sides of the event live on the *same* bridge, so season, weather and pipeline generation move with the collapse date too; the panels are read as that co-variation, **never** as damage.

```bash
python3 code/morandi/morandi_ocv_channels_csv.py    # build + verify the table
python3 code/morandi/fig_morandi_ocv_months.py     # this figure + report
```
