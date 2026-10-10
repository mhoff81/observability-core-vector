# Reference vs target (DiD) — Champlain Towers South

- event: **2021-06-24**
- comparison: `target` vs `control_cte`
- roles present: control_cte, control_ctn, reference_beach, target
- rows: 712 raw -> 712 observations

## ASCENDING

| channel | n | target Δ | reference Δ | **DiD** | t | p | pre comparable p | reference steps p | pre-trend ρ target/ref |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `phase_coherence` | 356 | -0.001187 | -0.003351 | **+0.001819** | +0.90 | 0.369 | 0.875 | 0.073 | +0.16/-0.05 |
| `phase_snr_db` | 0 | — | — | — | — | — | — | — | eine der vier Zellen hat < 5 Werte |
| `phase_scatterer_count` | 356 | -12.5 | -29 | **-2.5** | -0.19 | 0.850 | 1.000 | 0.167 | +0.05/-0.09 |
| `displacement_los_m` | 0 | — | — | — | — | — | — | — | eine der vier Zellen hat < 5 Werte |
| `coherence_gamma2` | 356 | -3.625e-05 | -2.836e-05 | **-1.14e-05** | -0.14 | 0.890 | 0.671 | 0.394 | +0.10/-0.00 |
| `coherence` | 350 | -0.02372 | -0.009867 | **-0.008225** | -0.24 | 0.807 | 0.904 | 0.029 | +0.25/+0.09 |
| `coherence_masked_pixels` | 356 | +0 | +0 | **-5.07e-14** | -0.23 | 0.821 | n/a | n/a | n/a/n/a |
| `intensity` | 356 | +643.3 | +145.2 | **-118.6** | -0.08 | 0.934 | 0.838 | 0.469 | -0.08/-0.13 |
| `brightness_ratio` | 356 | +0.016 | +0.003615 | **-0.002971** | -0.08 | 0.934 | 0.886 | 0.469 | -0.08/-0.13 |

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
