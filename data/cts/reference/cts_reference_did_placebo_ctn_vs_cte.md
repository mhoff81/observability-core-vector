# Reference vs target (DiD) — Champlain Towers South

- event: **2021-06-24**
- comparison: `control_ctn` vs `control_cte`
- roles present: control_cte, control_ctn, reference_beach, target
- rows: 712 raw -> 712 observations

## ASCENDING

| channel | n | target Δ | reference Δ | **DiD** | t | p | pre comparable p | reference steps p | pre-trend ρ target/ref |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `phase_coherence` | 356 | -0.002953 | -0.003351 | **-0.001059** | -0.54 | 0.587 | 0.828 | 0.073 | -0.07/-0.05 |
| `phase_snr_db` | 0 | — | — | — | — | — | — | — | eine der vier Zellen hat < 5 Werte |
| `phase_scatterer_count` | 356 | -16 | -29 | **+7.19** | +0.52 | 0.603 | 1.000 | 0.167 | -0.05/-0.09 |
| `displacement_los_m` | 0 | — | — | — | — | — | — | — | eine der vier Zellen hat < 5 Werte |
| `coherence_gamma2` | 356 | -4.79e-05 | -2.836e-05 | **-5.304e-05** | -0.71 | 0.478 | 0.352 | 0.394 | -0.00/-0.00 |
| `coherence` | 350 | -0.004949 | -0.009867 | **+0.01294** | +0.38 | 0.703 | 0.591 | 0.029 | -0.00/+0.09 |
| `coherence_masked_pixels` | 356 | +0 | +0 | **-5.07e-14** | -0.23 | 0.821 | n/a | n/a | n/a/n/a |
| `intensity` | 356 | +3499 | +145.2 | **+2712** | +1.84 | 0.067 | 0.000 | 0.469 | +0.07/-0.13 |
| `brightness_ratio` | 356 | +0.07695 | +0.003615 | **+0.05726** | +1.66 | 0.097 | 0.866 | 0.469 | +0.07/-0.13 |

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
