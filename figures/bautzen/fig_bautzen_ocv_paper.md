# Bautzen deck-edge OCV — figures A–E

Recomputation of the five figures of the Bautzen deck-edge observability package from the committed channel table. The table is the only input; every number below is recomputed here and pinned against the two committed reference JSONs (section 5).

|  |  |
|---|---|
| input | `data/bautzen/bautzen_ocv_channels.csv` (sha256 `32b2180cc28551ae…`) |
| rows | 142 overpasses in 4 rect series (RS n=28, DS1 n=92, DS2 n=22) |
| mask layer | `bautzen_deck_edge_mask.json` and `bautzen_ocv_channels_meta.json` (rect-file sha256s verified) |
| mode | full recomputation |
| figures | `figures/bautzen/fig_bautzen_ocv_paper_A..E.png` (600 dpi) |
| result JSON | `data/bautzen/fig_bautzen_ocv_paper.json` |
| pins | 8656 checks, 0 failures, 56 doc-string translations |
| cross-checks | 43 meta checks, 0 failed |

## 1. Question and design

The site is a pedestrian bridge (`bautzen_openbridge`): a bare deck edge crosses the SAR pixel grid of the InSAR OCV (observability core vector), so a bright, geometrically fixed line is expected in the coherence image. The deck edge was modified twice in the deployment (`12 May 2025` and `29 Sep 2025`); the three resulting states are RS (before, until 2025-05-12 10:00 local (before DS1 was installed)), DS1 (2025-05-12 10:00 .. 2025-09-29 10:00 local) and DS2 (from 2025-09-29 10:00 local).

Two design properties decide how the numbers below may be read:

1. **Four series, two of them counter-checks.** The same geometry was analysed for both orbit directions and for a control rectangle 500 m east of the bridge. Only a difference that appears in the bridge series and *not* in the control is evidence about the deck edge; the control fixes orbit, season and acquisition effects.
2. **State and season are confounded.** RS, DS1 and DS2 are time intervals of one deployment, so a state effect and a seasonal effect are the same contrast unless the control series moves with the bridge. The single cut dates are therefore tested against cuts shifted by +/-2, 4, 6, 8, 10, 12 weeks (section C), and the season is probed directly with air temperature / humidity at the overpass instant (figure E).

The vector examined here is the site's coherence channel `gamma^2` of the deck band *plus* the five mask-morphology channels `A`, `D`, `F`, `S`, `P` of `bautzen_ocv_masks.py` (area, density, fragmentation, peak shift, persistence). The mask layer is a documented deviation from the LUMO package (section 6).

## 2. Data

One row per overpass and rect series (the committed channel table). The mask share counts the dates on which the deck-edge mask rule of the origin analysis fires (band contrast and the deck-edge pixel count above the detection gate).

| series | n | n RS/DS1/DS2 | dates with mask | median contrast | share >= thr | abs peak offset med | share within 1 px | n_masked med | detectable |
|---|---|---|---|---|---|---|---|---|---|
| ASC bridge | 36 | 7/23/6 | 27/36 | 3.13 | 0.75 | 1.0 | 0.64 | 18.0 | yes |
| DESC bridge | 35 | 7/23/5 | 4/35 | 1.00 | 0.11 | 3.0 | 0.20 | 9.0 | no |
| ASC control 500 m east | 36 | 7/23/6 | 1/36 | 0.38 | 0.03 | 3.0 | 0.03 | 3.0 | no |
| DESC control 500 m east | 35 | 7/23/5 | 11/35 | 0.92 | 0.31 | 3.0 | 0.11 | 9.0 | no |

Detection gate of the origin analysis: band contrast >= 1.5 and at least 2 deck-edge pixels. The `gamma^2` channel is the floor-corrected deck-band coherence (`gamma2_band_floor_corrected`); the raw estimator and the house/peak variant stay in the table as audit columns (`bautzen_ocv_channels.py`).

## 3. Result at a glance

**Bridge/control discrimination of the band contrast** — ratio of the bridge and control medians, paired per date, percentile bootstrap 95% CI:

| orbit direction | RS | DS1 | DS2 |
|---|---|---|---|
| ASC | 6.00 [1.16, 12.67] (n=7) | 7.90 [5.92, 14.20] (n=23) | 6.15 [2.38, 13.78] (n=6) |
| DESC | 0.51 [0.37, 0.66] (n=7) | 1.32 [0.81, 1.75] (n=23) | 1.39 [0.60, 1.83] (n=5) |

**Step at the two cut dates** — median after − median before, two-sided permutation p, and the number of shifted-cut placebos (`k/n`) whose |delta| reaches the observed one:

| band contrast step | DS1 (12 May 2025) | DS2 (29 Sep 2025) |
|---|---|---|
| ASC bridge | -0.82  p=0.889  8/8 | +0.54  p=0.704  7/8 |
| DESC bridge | +0.28  p=0.139  1/8 | -0.09  p=0.738  3/7 |
| ASC control 500 m east | -0.23  p=0.069  1/8 | +0.23  p=0.062  1/8 |
| DESC control 500 m east | -0.69  p=0.035  3/8 | -0.33  p=0.310  5/7 |

| gamma^2 step | DS1 | DS2 |
|---|---|---|
| ASC bridge | band +0.053 (p=0.418, 7/8); band, n_band >= 5 +0.053 (p=0.418, 7/8); house -0.313 (p=0.021, 6/8) | band +0.088 (p=0.086, 1/8); band, n_band >= 5 +0.088 (p=0.086, 1/8); house -0.156 (p=0.372, 8/8) |
| DESC bridge | band -0.022 (p=0.729, 5/8); band, n_band >= 5 -0.023 (p=0.698, 4/8); house -0.060 (p=0.559, 5/8) | band -0.032 (p=0.495, 4/7); band, n_band >= 5 -0.026 (p=0.590, 4/7); house -0.277 (p=0.249, 1/7) |
| ASC control 500 m east | band -0.424 (p=0.107, 2/8); band, n_band >= 5 +0.128 (p=0.712, 6/8); house -0.149 (p=0.106, 6/8) | band -0.169 (p=0.524, 6/8); band, n_band >= 5 -0.393 (p=0.184, 1/7); house -0.169 (p=0.077, 5/8) |
| DESC control 500 m east | band +0.030 (p=0.339, 1/8); band, n_band >= 5 +0.031 (p=0.273, 1/8); house +0.027 (p=0.407, 6/8) | band +0.003 (p=0.949, 7/7); band, n_band >= 5 -0.032 (p=0.430, 3/7); house -0.038 (p=0.380, 1/7) |

**Season proxy** — gamma^2 against the air temperature at the overpass instant (OLS, one row per series; `r` and the two-sided p of the slope):

| series | n | r(T) | p(T) | r(RH) | p(RH) |
|---|---|---|---|---|---|
| ASC bridge | 36 | -0.22 | 0.200 | +0.23 | 0.182 |
| DESC bridge | 35 | +0.04 | 0.827 | -0.14 | 0.408 |
| ASC control 500 m east | 22 | +0.04 | 0.854 | -0.01 | 0.965 |
| DESC control 500 m east | 35 | -0.09 | 0.591 | +0.11 | 0.521 |

## 4. Figures

### Figure A — the six OCV components per series and state

One dot per overpass (state-coloured), the bar is the median of the series and state. The bridge series carries a fixed geometric edge in every state; the control series is the same rectangle shifted 500 m east and only serves as the counter-check.

Median `gamma^2` (deck band, floor-corrected) per series and state:

| series | RS (n=7) | DS1 (n=23) | DS2 (n=6) |
|---|---|---|---|
| ASC bridge | -0.010 | 0.022 | 0.104 |
| DESC bridge | -0.034 | -0.055 | -0.081 |
| ASC control 500 m east | 0.372 | -0.017 | -0.052 |
| DESC control 500 m east | -0.059 | -0.028 | -0.029 |

The descriptives of all six channels (n, median, mean, standard deviation, the 5/25/75/95 percentiles) are in `fig_a.descriptives` of the result JSON.

### Figure B — the deck-edge mask layer

Panels: (a) band contrast per series and state against the detection gate, (b) share of dates with a mask and its Wilson 95% CI, (c) the number of deck-edge pixels, (d) the in-series gamma^2 effect size (Cliff's delta with the RS state as reference).

| series | RS: contrast med | DS1: contrast med | DS2: contrast med | share with mask | detectable |
|---|---|---|---|---|---|
| ASC bridge | 3.58 | 2.76 | 3.67 | 27/36 | yes |
| DESC bridge | 0.77 | 1.15 | 0.93 | 4/35 | no |
| ASC control 500 m east | 0.60 | 0.35 | 0.60 | 1/36 | no |
| DESC control 500 m east | 1.53 | 0.87 | 0.67 | 11/35 | no |

### Figure C — intervention response against shifted cuts

The observed cut is special only if it is the most extreme cut of the series. `k/n` counts the shifted cuts whose |delta| reaches the observed one, so `k = n` means 'no placebo cut is as large' and `real is most extreme` says whether the observed cut leads the ranking.

| series | cut | band contrast delta | p | placebos >= observed | max abs placebo delta | real most extreme |
|---|---|---|---|---|---|---|
| ASC bridge | DS1 | -0.816 | 0.889 | 8/8 | 2.06 | no |
| ASC bridge | DS2 | +0.540 | 0.704 | 7/8 | 2.36 | no |
| DESC bridge | DS1 | +0.280 | 0.139 | 1/8 | 0.28 | yes |
| DESC bridge | DS2 | -0.087 | 0.738 | 3/7 | 0.16 | no |
| ASC control 500 m east | DS1 | -0.229 | 0.069 | 1/8 | 0.23 | yes |
| ASC control 500 m east | DS2 | +0.234 | 0.062 | 1/8 | 0.37 | no |
| DESC control 500 m east | DS1 | -0.690 | 0.035 | 3/8 | 1.24 | no |
| DESC control 500 m east | DS2 | -0.333 | 0.310 | 5/7 | 0.48 | no |

### Figure D — bridge minus control, date by date

The paired difference removes everything the two rectangles have in common (orbit, date, weather). A bridge-specific intervention effect must show up here, in the difference of the two bridge rectangles, and *not* in the single-series step tests of the two control rectangles.

| channel | geom | RS | DS1 | DS2 |
|---|---|---|---|---|
| band contrast | ASC | +2.984 (n=7) | +2.198 (n=23) | +3.118 (n=6) |
| band contrast | DESC | -1.041 (n=7) | +0.086 (n=23) | +0.214 (n=5) |
| gamma2_band_floor_corrected | ASC | -0.407 (n=6) | +0.084 (n=10) | +0.188 (n=6) |
| gamma2_band_floor_corrected | DESC | +0.002 (n=7) | -0.041 (n=23) | -0.038 (n=5) |

| channel | geom | n pairs | step DS1 (p) | step DS2 (p) | placebo at DS1 |
|---|---|---|---|---|---|
| band contrast | ASC | 36 | -0.786 (p=1.000) | +0.527 (p=0.715) | not extreme |
| band contrast | DESC | 35 | +1.191 (p=0.004) | +0.187 (p=0.676) | not extreme |
| gamma2_band_floor_corrected | ASC | 22 | +0.578 (p=0.087) | +0.325 (p=0.360) | not extreme |
| gamma2_band_floor_corrected | DESC | 35 | -0.041 (p=0.792) | -0.018 (p=0.993) | not extreme |

### Figure E — seasonality and environment

Monthly medians (top) and gamma^2 against temperature / humidity at the overpass instant (bottom). Because the three states are consecutive time intervals, any state difference is *also* a seasonal difference; this figure makes that confound visible instead of hiding it.

Monthly median band contrast (empty where the series has no overpass):

| series | 2025-04 | 2025-05 | 2025-06 | 2025-07 | 2025-08 | 2025-09 | 2025-10 |
|---|---|---|---|---|---|---|---|
| ASC bridge | 3.58 | 1.44 | 2.19 | 4.37 | 3.50 | 3.74 | 2.60 |
| DESC bridge | 0.89 | 0.74 | 1.15 | 1.16 | 0.81 | 1.27 | 0.93 |
| ASC control 500 m east | 0.60 | 0.37 | 0.29 | 0.33 | 0.39 | 0.43 | 0.74 |
| DESC control 500 m east | 1.53 | 1.01 | 0.87 | 1.51 | 1.08 | 0.58 | 0.67 |

Monthly median gamma^2 (empty where the series has no overpass):

| series | 2025-04 | 2025-05 | 2025-06 | 2025-07 | 2025-08 | 2025-09 | 2025-10 |
|---|---|---|---|---|---|---|---|
| ASC bridge | -0.01 | -0.06 | -0.02 | 0.04 | 0.07 | 0.03 | 0.11 |
| DESC bridge | 0.03 | -0.06 | -0.04 | -0.03 | -0.02 | -0.08 | -0.08 |
| ASC control 500 m east | 0.58 | -0.18 | 0.25 | 0.26 | 0.07 | 0.72 | -0.10 |
| DESC control 500 m east | -0.04 | -0.03 | 0.19 | -0.03 | -0.03 | -0.10 | -0.03 |

## 5. Verification against the committed Bautzen reference JSONs

Every block of the two committed reference files is recomputed from the channel table and compared leaf by leaf:

- `bautzen_deck_edge_mask.json` — the `config` block (10 constants of the mask module) and the four series blocks (`file`, `per_date` with its exact per-date key set, `summary`, `union50`, `union30`).
- `bautzen_deck_edge_state.json` — the `config` block (22 keys incl. the two German doc strings), the four series blocks (`per_date`, `metrics` with the `gamma2` sub-block, `jaccard`, `n_by_state`, `dates_by_state`) and the ten cross-series blocks (`discrimination`, `step_tests`, `placebo`, `delta_vs_ctrl`, `seasonality`, `gamma2_step_tests`, `gamma2_placebo`, `gamma2_delta_vs_ctrl`, `gamma2_seasonality`, `gamma2_environment`).

| pin result | value |
|---|---|
| checks | 8656 |
| exact (byte for byte) | 8656 |
| within tolerance | 0 (<= 1e-06) |
| failures | 0 |
| max |deviation| | 0.0 |
| doc strings translated | 56 |

**Comparison rules.** Deterministic numbers (medians, effect sizes, SDS, counters, mask sets, unions, gamma^2 values) must match exactly. Permutation and scipy p-values — `mannwhitney_p`, `p`, `p_perm`, `perm_p`, `rho`, `welch_p` — are compared with a tolerance of 1e-06, because they depend on the library's RNG and normalisation internals. Committed German doc strings are compared through `core.LABEL_EN` (an unmapped German string mismatches and becomes a failure). Committed JSON lists and recomputed tuples are normalised before the comparison.

**Meta cross-checks.** 43 independent checks against `bautzen_ocv_channels_meta.json` (0 failed): the file names and geometries of the four series, the per-date row counts, the mask summary, both persistence unions, the mask definition constants, the sha256 of the four rect files and the reference-source names.

**Keys beyond the reference.** The recomputed blocks carry 8 additional keys (the figure payloads and the per-variant gamma^2 gates); they are listed in `pins.recomputed_keys_not_in_reference` and are not failures — the reference files do not contain these values, they are added here.

**Failures: none.** Every pinned number of the two committed reference files is reproduced from the channel table alone.

## 6. Documented deviations from the LUMO figure package

The LUMO package (`code/lumo/fig_lumo_ocv_paper.py`) is the template of this script: same layout, same pin philosophy, same figure style. Four things differ by design, all of them documented in the header of this script and in `bautzen_ocv_masks.py`:

1. **The mask layer is site-specific.** The LUMO mask rules are echo based (the deck echo decides the mask). At this site that rule collapses — the deck-edge echo is not a separated peak — so the mask is *structurally* attributed: a deck band around the known deck-edge geometry (`bautzen_ocv_masks.py`). The house/peak variant stays in the table as a signal-based audit channel and is explicitly *not* a deck-edge statement.
2. **Four series instead of one state series.** Here the four rect series (bridge/control × ASC/DESC) are the counter-check, so every effect is read per series and the paired bridge-minus-control difference is the headline comparison (figure D).
3. **State and season are confounded.** RS, DS1 and DS2 are consecutive time intervals of a single deployment, so state and season cannot be separated by design. The script therefore reports the placebo sweep (figure C) and the environment proxy (figure E) instead of treating the state contrast as causal.
4. **The committed references are German, this report is English.** Doc strings are compared through `core.LABEL_EN`; the number of translated strings is reported in the pin summary (56 in this run). Series keys stay German — they are the committed file names.

`--quick` copies the placebo sweeps of the reference (`placebo`, `gamma2_placebo`, `delta_vs_ctrl[*].placebo`, `gamma2_delta_vs_ctrl[*].placebo`) instead of recomputing them; the copied blocks are recorded in `meta.copied_blocks`. This report was produced in **full** mode.

## 7. Reproduce

```sh
cd <repo root>
python3 code/bautzen/fig_bautzen_ocv_paper.py            # full run
python3 code/bautzen/fig_bautzen_ocv_paper.py --quick    # placebo copied
bash figures/bautzen/fig_bautzen_ocv_paper.sh          # the same, via the shell wrapper
```

The figures need numpy, matplotlib (600 dpi PNG) and scipy. This report was generated with python 3.13.5, numpy 2.4.6, matplotlib 3.11.1, scipy 1.18.0. All random draws use fixed seeds and no output carries a timestamp or a run time, so repeated runs produce identical JSON, Markdown and PNG files.

## 8. Caveats

- **Small samples.** 28 overpasses are in RS in total, 22 in DS2; the per-series, per-state cells go down to 5 overpasses. Every interval in this report (Wilson, bootstrap) is wide for that reason, and the permutation p-values are quantised by the number of permutations and — in figure C — by the number of placebo cuts.
- **The vector is not independent.** `gamma^2` and the mask dimensions describe the same band (`gamma^2` is the coherence *inside* the mask); the DIMS are complementary, not orthogonal, evidence.
- **Thresholds are inherited.** The detectability gate (band contrast, deck-edge pixel count) and the placebo offsets (+/-2..12 weeks) are taken over from the origin analysis, not re-tuned here.
- **State assignment.** States are assigned by overpass time (`core.state_of`), so a date inside the DS1 window is labelled DS1 even if the intervention left no trace in the image. Figures B and C test exactly that, they do not assume it.
- **Season is not removed.** The temperature/humidity regression of figure E is a probe, not a correction: the states are months apart and the covariate is collinear with them by construction.

---

Generated by `code/bautzen/fig_bautzen_ocv_paper.py` (full mode, 8656 pins, 0 failures).
