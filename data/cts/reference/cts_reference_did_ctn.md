# Reference vs target (DiD) — Champlain Towers South

- event: **2021-06-24**
- comparison: `target` vs `control_ctn`
- roles present: control_cte, control_ctn, reference_beach, target
- rows: 712 raw -> 712 observations

## ASCENDING

| channel | n | target Δ | reference Δ | **DiD** | t | p | pre comparable p | reference steps p | pre-trend ρ target/ref |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `phase_coherence` | 356 | -0.001187 | -0.002953 | **+0.002878** | +1.45 | 0.148 | 0.957 | 0.004 | +0.16/-0.07 |
| `phase_snr_db` | 0 | — | — | — | — | — | — | — | eine der vier Zellen hat < 5 Werte |
| `phase_scatterer_count` | 356 | -12.5 | -16 | **-9.69** | -0.71 | 0.480 | 1.000 | 0.574 | +0.05/-0.05 |
| `displacement_los_m` | 0 | — | — | — | — | — | — | — | eine der vier Zellen hat < 5 Werte |
| `coherence_gamma2` | 356 | -3.625e-05 | -4.79e-05 | **+4.164e-05** | +0.53 | 0.596 | 0.191 | 0.017 | +0.10/-0.00 |
| `coherence` | 350 | -0.02372 | -0.004949 | **-0.02116** | -0.65 | 0.518 | 0.667 | 0.090 | +0.25/-0.00 |
| `coherence_masked_pixels` | 356 | +0 | +0 | **-5.07e-14** | -0.23 | 0.821 | n/a | n/a | n/a/n/a |
| `intensity` | 356 | +643.3 | +3499 | **-2830** | -1.88 | 0.060 | 0.000 | 0.004 | -0.08/+0.07 |
| `brightness_ratio` | 356 | +0.016 | +0.07695 | **-0.06023** | -1.71 | 0.088 | 0.755 | 0.004 | -0.08/+0.07 |

## DESCENDING

| channel | n | target Δ | reference Δ | **DiD** | t | p | pre comparable p | reference steps p | pre-trend ρ target/ref |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `phase_coherence` | 0 | — | — | — | — | — | — | — | eine der vier Zellen hat < 5 Werte |
| `phase_snr_db` | 0 | — | — | — | — | — | — | — | eine der vier Zellen hat < 5 Werte |
| `phase_scatterer_count` | 0 | — | — | — | — | — | — | — | eine der vier Zellen hat < 5 Werte |
| `displacement_los_m` | 0 | — | — | — | — | — | — | — | eine der vier Zellen hat < 5 Werte |
| `coherence_gamma2` | 0 | — | — | — | — | — | — | — | eine der vier Zellen hat < 5 Werte |
| `coherence` | 0 | — | — | — | — | — | — | — | eine der vier Zellen hat < 5 Werte |
| `coherence_masked_pixels` | 0 | — | — | — | — | — | — | — | eine der vier Zellen hat < 5 Werte |
| `intensity` | 0 | — | — | — | — | — | — | — | eine der vier Zellen hat < 5 Werte |
| `brightness_ratio` | 0 | — | — | — | — | — | — | — | eine der vier Zellen hat < 5 Werte |

## Reading

- A significant DiD in a channel whose **reference also steps** (`reference steps p` small) is not a DiD result — the differencing assumption failed.
- A significant DiD with a diverging pre-trend or a small `pre comparable p` is reported as such, not as evidence.
- The DiD removes *common* confounders (atmosphere, soil moisture, orbits). It cannot remove anything that affects only the bridge — that is what the temperature and season terms of the model ladder are for.
