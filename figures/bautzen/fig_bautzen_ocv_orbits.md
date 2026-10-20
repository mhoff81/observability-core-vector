# Bautzen OCV time series by dimension and orbit

## 1. Question and layout

How do the six OCV components of the Bautzen OpenLabs bridge deck edge move through the three documented states RS, DS1 and DS2, and do the two Sentinel-1 pass geometries tell the same story? The figure draws one panel per **dimension x orbit**: six rows in the formal OCV order (`gamma2`, `P`, `D`, `F`, `A`, `S`) and two columns (ASC IW2, DESC IW3), i.e. twelve panels.

The x coordinate is the **actual local overpass instant**, not midnight. The state boundaries fall in the middle of a day (DESC 07:08, ASC 18:51 local), so a date-only axis would place the boundary-day points in the wrong state band. The three state windows are background bands with different colours and the two cut dates are dashed vertical lines.

Only the bridge series is drawn. The control rectangle 500 m east is the counter-check of the paper figure, but its ASC `gamma2` range (-0.94 .. +0.93) is much wider than the bridge range (-0.14 .. +0.24); a shared panel would compress the bridge signal this figure is meant to show. The control remains in `data/bautzen/bautzen_ocv_channels.csv` and in figures A-E.

## 2. Data and state model

- source: `data/bautzen/bautzen_ocv_channels.csv`
- committed metadata: `data/bautzen/bautzen_ocv_channels_meta.json`
- overpasses: **142** in four rect series (ASC Brücke: 36, DESC Brücke: 35, ASC Kontrolle 500 m östlich: 36, DESC Kontrolle 500 m östlich: 35)
- states: **28 RS**, **92 DS1**, **22 DS2**
- cuts: DS1 `2025-05-12T10:00:00`, DS2 `2025-09-29T10:00:00` local

State assignment is by overpass time, not date:

- **RS**: until 2025-05-12 10:00 local (before DS1 was installed)
- **DS1**: 2025-05-12 10:00 .. 2025-09-29 10:00 local
- **DS2**: from 2025-09-29 10:00 local

The bridge series used by the figure:

| orbit | series | dates | RS | DS1 | DS2 |
| --- | --- | --- | --- | --- | --- |
| ASC | ASC Brücke | 36 | 7 | 23 | 6 |
| DESC | DESC Brücke | 35 | 7 | 23 | 5 |

## 3. Per-state medians

Median of each channel inside each state, with the number of finite overpasses in parentheses. The bridge panels are complete (all 36 ASC and 35 DESC overpasses); the plotting path would leave any missing value as a gap, and nothing is interpolated.

| orbit | dimension | RS median (n) | DS1 median (n) | DS2 median (n) |
| --- | --- | --- | --- | --- |
| ASC | gamma2 | -0.010 (n=7) | 0.022 (n=23) | 0.104 (n=6) |
| ASC | P | 0.704 (n=7) | 0.727 (n=23) | 0.700 (n=6) |
| ASC | D | 0.667 (n=7) | 0.667 (n=23) | 0.722 (n=6) |
| ASC | F | 1.000 (n=7) | 1.000 (n=23) | 1.000 (n=6) |
| ASC | A | 18.000 (n=7) | 18.000 (n=23) | 19.500 (n=6) |
| ASC | S | 2.082 (n=7) | 2.239 (n=23) | 3.059 (n=6) |
| DESC | gamma2 | -0.034 (n=7) | -0.055 (n=23) | -0.081 (n=5) |
| DESC | P | 0.367 (n=7) | 0.357 (n=23) | 0.376 (n=5) |
| DESC | D | 0.381 (n=7) | 0.417 (n=23) | 0.370 (n=5) |
| DESC | F | 2.000 (n=7) | 2.000 (n=23) | 2.000 (n=5) |
| DESC | A | 7.000 (n=7) | 9.000 (n=23) | 10.000 (n=5) |
| DESC | S | 2.508 (n=7) | 2.167 (n=23) | 2.276 (n=5) |

The corresponding means, minima and maxima are in the result JSON under `summary.<orbit>.<dimension>.<state>`.

## 4. Cross-checks

The script recomputes every plotted point from the committed channel table and checks it against the committed metadata: **41 checks, 0 failures**. The checks cover the row count, the state rule, the per-series date/state counts, uniqueness of dates, `A == n_masked`, the presence and completeness of every channel on both bridge orbits, and the sha256 of every rect file.

**Failures: none.** Every state value reproduces `core.state_of`, every per-series count and date list reproduces the committed metadata, and every rect-file sha256 matches.

## 5. Reproduce

```sh
cd <repo root>
python3 code/bautzen/fig_bautzen_ocv_orbits.py
PYTHON=/usr/bin/python3 bash figures/bautzen/fig_bautzen_ocv_orbits.sh
```

Needs numpy and matplotlib. This report was generated with python 3.13.5, numpy 2.4.6, matplotlib 3.11.1. No output carries a timestamp, so repeated runs are identical.

## 6. Caveats

- **States are consecutive time intervals, not a randomized treatment.** RS, DS1 and DS2 are one deployment observed over six months; state and season are confounded. The background bands make that structure visible, they do not remove it.
- **No control series in this figure.** The control rectangle is the counter-check in figures A-E; it is omitted here only because its ASC `gamma2` range would flatten the bridge signal in a shared panel.
- **Missing values are gaps.** The bridge series are complete: all 36 ASC and 35 DESC overpasses carry all six channels. The 14 missing `gamma2` rows and the 8 missing `D`, `S`, `P` rows are in the control series, which this figure does not draw; no value is interpolated.
- **`gamma2` is the floor-corrected deck-band coherence** (`gamma2_band_floor_corrected`), the same channel as paper figure A, not the raw estimator.
- **The x axis is local overpass time.** ASC is drawn at 18:51 and DESC at 07:08 local; this is what makes the state assignment at the two cut days correct.

---

Generated by `code/bautzen/fig_bautzen_ocv_orbits.py` (41 cross-checks, 0 failures).
