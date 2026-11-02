# Espoo mast (OSM 2233589804) — support chips — the OCV vector of the support chips

* table `data/espoo/support/espoo_ocv_channels.csv` (sha256 `82a8cc04208e5105…`), 75 chips
* 2024-09-05 … 2026-09-07; orbits: ASCENDING 31, DESCENDING 44; **3 shared dates**
* **second instrument**: these are chip values, not the pipeline columns of the retired pipeline figure (frozen at `data/espoo/fig_espoo_ocv_orbits.json`: r(log A, log `coherence_masked_pixels`) = −0.24, r(`gamma2_mask`, `coherence_gamma2`) = 0.13)
* **two tracks**: ASC = relative orbit 160, DESC = 51; the other two tracks place the mast on the far range edge, where the 200 m ring runs past the swath
* **no pairing**: the tracks have three shared dates, so the orbit contrast below is a two-sample comparison
* **tautology**: `D = A / bbox_area` exactly, and `gamma2_mask` carries `n = A` in its denominator; `A` differs by a factor 7 between the orbits, so both dimensions are dominated by the mask size

## The figure — one 170 mm plate

* **canvas** `6.6929 × 7.5500 in` = **170.00 × 191.77 mm**; measured from the written files: eps 170.00 × 191.77 mm, jpeg 169.98 × 191.75 mm, pdf 170.00 × 191.77 mm. The plate is authored at the print width — it is not cropped to its content, so nothing here is rescaled.
* **typography** header 9.0 · panel title 8.8 · axis label 7.6 · ticks 6.8 · legend 6.8 · statistics 6.6 pt — nothing below 6 pt, every role inside the `figure_style` band.
* **panels**, reading order: `A — observable area (cells)` · `D — density (fraction of bbox)` · `F — fragmentation (components)` · `S — spatial shift (cells)` · `P — persistence (pixel frequency)` · `$\gamma^2$ — coherence`. The column names of the table are not on the plate: `gamma2_mask` is drawn as the coherence it measures.
* **statistics rows**: one fixed band above every panel, identical in all six — DESC at +10.0 pt, ASC at +1.0 pt, 6.6 pt, flush with the left spine. Each row is the Theil–Sen slope per year, its bootstrap 95 % interval, the Mann–Kendall p and the verdict of the ported flag rule (`trend` when it fires, `flat` when it does not).
* **marks** 3.6 pt markers with a 0.5 pt white edge above the 1.2 pt trend line; solid where that verdict is `trend`, dashed where it is `flat`; the colour is the track's. One figure legend carries the track counts at the right end of the header line.
* **bounds** content 168.0 × 190.1 mm inside the 170.0 × 191.8 mm canvas, 0.77 mm to the nearest edge (floor 0.3 mm, asserted) — a guard, not a crop. LaTeX: the file is 170 mm wide — included with `width=\linewidth`; no trim, clip or resizebox.
* **EPS** 108 translucent artist(s) pre-blended over white (`eps="blend"`), 0 rasterised, text Type 42 (pdf/ps fonttype 42): the EPS is vector, like the PDF.
* **moved into this caption**: the orbit ids (ASC = relative orbit 160, DESC = 51), the three shared dates, the trend semantics and the definitions of the six dimensions.

```json
{
 "target_width_mm": 170.0,
 "target_height_mm": 191.77,
 "figsize_inches": [
  6.692913,
  7.55
 ],
 "native_dpi": 630,
 "final_font_sizes_pt": {
  "suptitle": 9.0,
  "title": 8.8,
  "label": 7.6,
  "tick": 6.8,
  "legend": 6.8,
  "annotation": 6.6
 },
 "native_width_mm": {
  "jpeg": 169.98,
  "pdf": 170.0,
  "eps": 170.0
 },
 "native_height_mm": {
  "jpeg": 191.75,
  "pdf": 191.77,
  "eps": 191.77
 },
 "latex_display_width_mm": 170.0,
 "latex_scaling": "the file is 170 mm wide — included with `width=\\linewidth`; no trim, clip or resizebox",
 "bounds": {
  "export": "canvas-authorised (figure_style.save); bbox_inches is never 'tight'",
  "canvas_mm": [
   170.0,
   191.77
  ],
  "content_mm": [
   168.0,
   190.1
  ],
  "edge_mm": {
   "left": 1.02,
   "right": 1.02,
   "bottom": 0.91,
   "top": 0.77
  },
  "minimum_edge_mm": 0.3
 },
 "eps": {
  "alpha_artists": 108,
  "strategy": "blend",
  "background": "white",
  "rasterised_artists": 0,
  "vector_text": "Type 42 (pdf/ps fonttype 42)"
 },
 "panel_order": [
  "A",
  "D",
  "F",
  "S",
  "P",
  "gamma2_mask"
 ],
 "statistics": {
  "pt": 6.6,
  "row_offset_pt": {
   "ASCENDING": 1.0,
   "DESCENDING": 10.0
  },
  "anchor": "x=0.0 axes, ha=left, va=bottom -- one rule, six panels"
 },
 "panels": {
  "A": {
   "label": "A — observable area (cells)",
   "x_limits": [
    19934.4,
    20739.6
   ],
   "y_limits": [
    5.85,
    823.15
   ],
   "tracks": {
    "ASCENDING": {
     "n": 31,
     "first": "2024-09-05",
     "last": "2026-09-02",
     "median": 76.0,
     "iqr": [
      62.5,
      87.5
     ],
     "theil_sen_per_year": 0.0,
     "slope_ci95_per_year": [
      -25.785600148388333,
      23.67361111111111
     ],
     "mann_kendall_p": 1.0,
     "verdict": "no trend",
     "trend_line": {
      "first_ordinal": 739134.0,
      "last_ordinal": 739861.0,
      "dashed": true
     }
    },
    "DESCENDING": {
     "n": 44,
     "first": "2024-09-10",
     "last": "2026-09-07",
     "median": 527.5,
     "iqr": [
      392.75,
      616.75
     ],
     "theil_sen_per_year": -23.955609836065577,
     "slope_ci95_per_year": [
      -190.01237975543478,
      97.7238782670552
     ],
     "mann_kendall_p": 0.7385268571345529,
     "verdict": "no trend",
     "trend_line": {
      "first_ordinal": 739139.0,
      "last_ordinal": 739866.0,
      "dashed": true
     }
    }
   }
  },
  "D": {
   "label": "D — density (fraction of bbox)",
   "x_limits": [
    19934.4,
    20739.6
   ],
   "y_limits": [
    0.005208,
    0.182128
   ],
   "tracks": {
    "ASCENDING": {
     "n": 31,
     "first": "2024-09-05",
     "last": "2026-09-02",
     "median": 0.021001,
     "iqr": [
      0.018817968961647124,
      0.02368355855855856
     ],
     "theil_sen_per_year": 0.0001596900278540931,
     "slope_ci95_per_year": [
      -0.005993546232033069,
      0.005807223191444992
     ],
     "mann_kendall_p": 0.8650383150743488,
     "verdict": "no trend",
     "trend_line": {
      "first_ordinal": 739134.0,
      "last_ordinal": 739861.0,
      "dashed": true
     }
    },
    "DESCENDING": {
     "n": 44,
     "first": "2024-09-10",
     "last": "2026-09-07",
     "median": 0.116833,
     "iqr": [
      0.08698781838316723,
      0.1366002214839424
     ],
     "theil_sen_per_year": -0.005305782909427579,
     "slope_ci95_per_year": [
      -0.04208469097573305,
      0.021644269826590293
     ],
     "mann_kendall_p": 0.7385268571345529,
     "verdict": "no trend",
     "trend_line": {
      "first_ordinal": 739139.0,
      "last_ordinal": 739866.0,
      "dashed": true
     }
    }
   }
  },
  "F": {
   "label": "F — fragmentation (components)",
   "x_limits": [
    19934.4,
    20739.6
   ],
   "y_limits": [
    6.75,
    166.25
   ],
   "tracks": {
    "ASCENDING": {
     "n": 31,
     "first": "2024-09-05",
     "last": "2026-09-02",
     "median": 19.0,
     "iqr": [
      16.5,
      21.5
     ],
     "theil_sen_per_year": -0.951171875,
     "slope_ci95_per_year": [
      -3.381944444444444,
      2.3436874999999873
     ],
     "mann_kendall_p": 0.5268577144233753,
     "verdict": "no trend",
     "trend_line": {
      "first_ordinal": 739134.0,
      "last_ordinal": 739861.0,
      "dashed": true
     }
    },
    "DESCENDING": {
     "n": 44,
     "first": "2024-09-10",
     "last": "2026-09-07",
     "median": 133.0,
     "iqr": [
      121.75,
      142.0
     ],
     "theil_sen_per_year": -8.689537582069278,
     "slope_ci95_per_year": [
      -19.5762703274548,
      3.3796888197041737
     ],
     "mann_kendall_p": 0.05554339130406016,
     "verdict": "no trend",
     "trend_line": {
      "first_ordinal": 739139.0,
      "last_ordinal": 739866.0,
      "dashed": true
     }
    }
   }
  },
  "S": {
   "label": "S — spatial shift (cells)",
   "x_limits": [
    19934.4,
    20739.6
   ],
   "y_limits": [
    0.634631,
    33.372305
   ],
   "tracks": {
    "ASCENDING": {
     "n": 31,
     "first": "2024-09-05",
     "last": "2026-09-02",
     "median": 25.722471,
     "iqr": [
      23.316768584107486,
      27.99561296659574
     ],
     "theil_sen_per_year": -1.5497492140071458,
     "slope_ci95_per_year": [
      -7.00536849352294,
      3.0556406136789547
     ],
     "mann_kendall_p": 0.19645146482180031,
     "verdict": "no trend",
     "trend_line": {
      "first_ordinal": 739134.0,
      "last_ordinal": 739861.0,
      "dashed": true
     }
    },
    "DESCENDING": {
     "n": 44,
     "first": "2024-09-10",
     "last": "2026-09-07",
     "median": 14.386891,
     "iqr": [
      8.576674114098868,
      19.529021942506553
     ],
     "theil_sen_per_year": 0.78009139817415,
     "slope_ci95_per_year": [
      -4.777544491722997,
      2.9533783213159133
     ],
     "mann_kendall_p": 0.7538688170084293,
     "verdict": "no trend",
     "trend_line": {
      "first_ordinal": 739139.0,
      "last_ordinal": 739866.0,
      "dashed": true
     }
    }
   }
  },
  "P": {
   "label": "P — persistence (pixel frequency)",
   "x_limits": [
    19934.4,
    20739.6
   ],
   "y_limits": [
    0.106182,
    0.285987
   ],
   "tracks": {
    "ASCENDING": {
     "n": 31,
     "first": "2024-09-05",
     "last": "2026-09-02",
     "median": 0.233953,
     "iqr": [
      0.2216096096096096,
      0.24643817204301074
     ],
     "theil_sen_per_year": -0.021599254218647893,
     "slope_ci95_per_year": [
      -0.04399941288860608,
      0.0007312312312311192
     ],
     "mann_kendall_p": 0.013083919895177765,
     "verdict": "no trend",
     "trend_line": {
      "first_ordinal": 739134.0,
      "last_ordinal": 739861.0,
      "dashed": true
     }
    },
    "DESCENDING": {
     "n": 44,
     "first": "2024-09-10",
     "last": "2026-09-07",
     "median": 0.147182,
     "iqr": [
      0.1416604302459698,
      0.15736810405500515
     ],
     "theil_sen_per_year": 0.0030752288516493224,
     "slope_ci95_per_year": [
      -0.0039184953817895385,
      0.012841015259441427
     ],
     "mann_kendall_p": 0.33662512075788875,
     "verdict": "no trend",
     "trend_line": {
      "first_ordinal": 739139.0,
      "last_ordinal": 739866.0,
      "dashed": true
     }
    }
   }
  },
  "gamma2_mask": {
   "label": "$\\gamma^2$ — coherence",
   "x_limits": [
    19934.4,
    20739.6
   ],
   "y_limits": [
    -0.003211,
    0.067517
   ],
   "tracks": {
    "ASCENDING": {
     "n": 31,
     "first": "2024-09-05",
     "last": "2026-09-02",
     "median": 0.016387,
     "iqr": [
      0.005888811407042432,
      0.02941901234266612
     ],
     "theil_sen_per_year": -0.009932670463953035,
     "slope_ci95_per_year": [
      -0.03012555252733329,
      0.005337540220890599
     ],
     "mann_kendall_p": 0.11789500171022915,
     "verdict": "no trend",
     "trend_line": {
      "first_ordinal": 739134.0,
      "last_ordinal": 739861.0,
      "dashed": true
     }
    },
    "DESCENDING": {
     "n": 44,
     "first": "2024-09-10",
     "last": "2026-09-07",
     "median": 0.000704,
     "iqr": [
      0.00029068148385597784,
      0.0013553681476136511
     ],
     "theil_sen_per_year": -0.00017564750604610595,
     "slope_ci95_per_year": [
      -0.0010035929719896457,
      0.0002745260925630681
     ],
     "mann_kendall_p": 0.29751878206452487,
     "verdict": "no trend",
     "trend_line": {
      "first_ordinal": 739139.0,
      "last_ordinal": 739866.0,
      "dashed": true
     }
    }
   }
  }
 }
}
```

| panel | track | n | first | last | median | IQR | slope / a | 95% CI | MK p | verdict | y-limits |
|---|---|---|---|---|---|---|---|---|---|---|---|
| A — observable area (cells) | ASC R160 | 31 | 2024-09-05 | 2026-09-02 | 76 | [62.5, 87.5] | 0 | [-25.7856, 23.6736] | 1.000 | no trend | [5.85, 823.15] |
| A — observable area (cells) | DESC R51 | 44 | 2024-09-10 | 2026-09-07 | 527.5 | [392.75, 616.75] | -23.9556 | [-190.012, 97.7239] | 0.739 | no trend | [5.85, 823.15] |
| D — density (fraction of bbox) | ASC R160 | 31 | 2024-09-05 | 2026-09-02 | 0.021001 | [0.018818, 0.0236836] | 0.00015969 | [-0.00599355, 0.00580722] | 0.865 | no trend | [0.005208, 0.182128] |
| D — density (fraction of bbox) | DESC R51 | 44 | 2024-09-10 | 2026-09-07 | 0.116833 | [0.0869878, 0.1366] | -0.00530578 | [-0.0420847, 0.0216443] | 0.739 | no trend | [0.005208, 0.182128] |
| F — fragmentation (components) | ASC R160 | 31 | 2024-09-05 | 2026-09-02 | 19 | [16.5, 21.5] | -0.951172 | [-3.38194, 2.34369] | 0.527 | no trend | [6.75, 166.25] |
| F — fragmentation (components) | DESC R51 | 44 | 2024-09-10 | 2026-09-07 | 133 | [121.75, 142] | -8.68954 | [-19.5763, 3.37969] | 0.056 | no trend | [6.75, 166.25] |
| S — spatial shift (cells) | ASC R160 | 31 | 2024-09-05 | 2026-09-02 | 25.7225 | [23.3168, 27.9956] | -1.54975 | [-7.00537, 3.05564] | 0.196 | no trend | [0.634631, 33.372305] |
| S — spatial shift (cells) | DESC R51 | 44 | 2024-09-10 | 2026-09-07 | 14.3869 | [8.57667, 19.529] | 0.780091 | [-4.77754, 2.95338] | 0.754 | no trend | [0.634631, 33.372305] |
| P — persistence (pixel frequency) | ASC R160 | 31 | 2024-09-05 | 2026-09-02 | 0.233953 | [0.22161, 0.246438] | -0.0215993 | [-0.0439994, 0.000731231] | 0.013 | no trend | [0.106182, 0.285987] |
| P — persistence (pixel frequency) | DESC R51 | 44 | 2024-09-10 | 2026-09-07 | 0.147182 | [0.14166, 0.157368] | 0.00307523 | [-0.0039185, 0.012841] | 0.337 | no trend | [0.106182, 0.285987] |
| $\gamma^2$ — coherence | ASC R160 | 31 | 2024-09-05 | 2026-09-02 | 0.016387 | [0.00588881, 0.029419] | -0.00993267 | [-0.0301256, 0.00533754] | 0.118 | no trend | [-0.003211, 0.067517] |
| $\gamma^2$ — coherence | DESC R51 | 44 | 2024-09-10 | 2026-09-07 | 0.000704 | [0.000290681, 0.00135537] | -0.000175648 | [-0.00100359, 0.000274526] | 0.298 | no trend | [-0.003211, 0.067517] |


## 1. Drift across the two years

| dimension | orbit | n | slope / year | 95% CI | MK p | lag-1 | split-half delta | trend | serial flags | ties |
|---|---|---|---|---|---|---|---|---|---|---|
| `A` | ASC R160 | 31 | +0 | [-25.8, +23.7] | 1.000 | +0.038 | +0.19 | no | none | 0.06 at 61.0 |
| `A` | DESC R51 | 44 | -24 | [-190, +97.7] | 0.739 | +0.19 | -0.012 | no | none | 0.05 at 392.0 |
| `D` | ASC R160 | 31 | +0.00016 | [-0.00599, +0.00581] | 0.865 | +0.021 | +0.22 | no | none | 0.03 at 0.018889858 |
| `D` | DESC R51 | 44 | -0.00531 | [-0.0421, +0.0216] | 0.739 | +0.19 | -0.012 | no | none | 0.05 at 0.086821705 |
| `F` | ASC R160 | 31 | -0.951 | [-3.38, +2.34] | 0.527 | +0.0075 | +0.0042 | no | none | 0.16 at 18.0 |
| `F` | DESC R51 | 44 | -8.69 | [-19.6, +3.38] | 0.056 | +0.16 | -0.21 | no | none | 0.09 at 133.0 |
| `S` | ASC R160 | 31 | -1.55 | [-7.01, +3.06] | 0.196 | +0.14 | -0.2 | no | none | 0.03 at 28.096039857 |
| `S` | DESC R51 | 44 | +0.78 | [-4.78, +2.95] | 0.754 | +0.11 | -0.079 | no | none | 0.02 at 6.920030672 |
| `P` | ASC R160 | 31 | -0.0216 | [-0.044, +0.000731] | 0.013 | +0.06 | -0.37 | no | none | 0.03 at 0.271179487 |
| `P` | DESC R51 | 44 | +0.00308 | [-0.00392, +0.0128] | 0.337 | -0.068 | +0.033 | no | none | 0.02 at 0.162346939 |
| `gamma2_mask` | ASC R160 | 31 | -0.00993 | [-0.0301, +0.00534] | 0.118 | +0.037 | -0.36 | no | none | 0.03 at 0.028037644 |
| `gamma2_mask` | DESC R51 | 44 | -0.000176 | [-0.001, +0.000275] | 0.298 | +0.056 | -0.21 | no | none | 0.02 at 0.001871815 |

## 2. Year repeat — the same calendar month in both campaign years

| dimension | orbit | months | median Δ (y2 − y1) | 95% CI | months + / − |
|---|---|---|---|---|---|
| `A` | ASC R160 | 11 | +6 | [-27.5, +13] | 7 / 4 |
| `A` | DESC R51 | 12 | +6 | [-55.5, +162] | 7 / 5 |
| `D` | ASC R160 | 11 | +0.00173 | [-0.00738, +0.0045] | 6 / 5 |
| `D` | DESC R51 | 12 | +0.00133 | [-0.0123, +0.0359] | 7 / 5 |
| `F` | ASC R160 | 11 | +0.5 | [-2.5, +3] | 6 / 5 |
| `F` | DESC R51 | 12 | -5.75 | [-12.8, +19] | 6 / 6 |
| `S` | ASC R160 | 11 | -1.76 | [-3.77, +0.929] | 4 / 7 |
| `S` | DESC R51 | 12 | -0.165 | [-5.58, +8.76] | 6 / 6 |
| `P` | ASC R160 | 11 | -0.0184 | [-0.0317, -0.00238] | 2 / 9 |
| `P` | DESC R51 | 12 | -0.00275 | [-0.00858, +0.00796] | 5 / 7 |
| `gamma2_mask` | ASC R160 | 11 | -0.0193 | [-0.0199, +0.00868] | 5 / 6 |
| `gamma2_mask` | DESC R51 | 12 | -9.11e-05 | [-0.0013, +0.000327] | 6 / 6 |

## 3. Orbit contrast, unpaired (ASC R160 vs DESC R51)

| dimension | n ASC/DESC | median ASC | median DESC | Cliff's δ | 95% CI | Mann-Whitney p |
|---|---|---|---|---|---|---|
| `A` | 31/44 | 76 | 528 | -1 | [-1, -0.99] | <0.001 |
| `D` | 31/44 | 0.021 | 0.117 | -1 | [-1, -0.99] | <0.001 |
| `F` | 31/44 | 19 | 133 | -1 | [-1, -1] | <0.001 |
| `S` | 31/44 | 25.7 | 14.4 | +0.84 | [+0.7, +0.95] | <0.001 |
| `P` | 31/44 | 0.234 | 0.147 | +0.94 | [+0.85, +1] | <0.001 |
| `gamma2_mask` | 31/44 | 0.0164 | 0.000704 | +0.89 | [+0.74, +0.99] | <0.001 |

## What this can and cannot say

* **It is a different instrument from the pipeline figure.** The two must not be mixed: they correlate at −0.24 (`A`) and 0.13 (`gamma2`).
* **The series is halved and selected by track.** 75 of 150 overpasses; the excluded tracks (ASC R58, DESC R153) are named in the site README.
* **No paired orbit contrast.** Three shared dates; the contrast is a two-sample comparison and is reported as such.
* **`D` and `gamma2_mask` are functions of `A`**, which itself differs by a factor 7 between the orbits.

