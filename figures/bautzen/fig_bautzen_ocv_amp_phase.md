# Bautzen deck-edge amplitude and phase — figure, per series and state

## 1. Question and definition

How do the amplitude A(t) and the phase phi(t) of the deck band respond to the two documented interventions at the Bautzen OpenLabs bridge (RS -> DS1 -> DS2)? Both are measured **over the deck-edge pixels only** — the same local, geometry-anchored deck band the channel table uses — and not over the whole bridge chip:

- `A(t)   = |sum z| / sqrt(N * sum |z|^2)`  — coherence amplitude, in [0, 1]
- `phi(t) = arg(sum z)`  — interferometric phase of the deck edge, in (-pi, pi]
- `dphi(t) = wrap(phi_bridge(t) - phi_control(t))`  — paired differential phase of the deck edge, in (-pi, pi]

with `z` the complex chip of the overpass and the sums over the `N` deck-edge pixels `z_sel = z[mask]`. `A(t)^2` is exactly this site's coherence channel, so it is cross-checked against the committed `gamma2_band_raw` column (section 4). The four series (bridge/control × ASC/DESC) are carried through the amplitude and phase panels; the control rectangle 500 m east is the counter-check against a purely seasonal reading and the subtrahend of the differential panel.

The differential panel works on the **pair**. Bridge and control rectangle sit on the same overpasses (identical acquisition times), so their deck-edge phases are differenced pair by pair without any interpolation — wherever both rectangles have a valid deck-edge mask. The difference removes what both rectangles see in common — season, atmosphere, repeat-pass baseline — and keeps the local deck-edge term, which is exactly what the control rectangle is for. `wrap` folds the raw difference back onto the (-pi, pi] branch of `arg`.

Section 4 closes with the **state-to-state** change of both quantities, `delta = wrap(circmean(state) - circmean(RS))` — the shift from the reference state RS to DS1 and to DS2 — for the paired `dphi` and, as counter-check, for the absolute `phi` of each series. It is reported with a seeded bootstrap confidence interval, because a state holds only a handful of overpasses.

The deck-edge mask is rebuilt here with the *same* routine and the same floats as the channel table (`masks.date_mask` on `|z|^2`), so the pixel set is bit-identical and `A(t)^2` reproduces the pinned channel value.

## 2. Inputs and provenance

| input | value |
| --- | --- |
| channel table | `data/bautzen/bautzen_ocv_channels.csv` |
| channel table sha256 | `32b2180cc28551ae85056c8b5994d6dedf02cc724353b2d130880230078c37fa` |
| overpasses in the table | 142 |
| rect files (complex chips) | `bautzen_rects_asc_iw2.txt` (8657db203cfca7e7324b498f515b377963cbab66bc9c8ad3bac6386542fb4ab5), `bautzen_rects_desc_iw3.txt` (ea430cb8338f972fbe8313937ba5724b2f9fa68e5fde587fb2503f293c2ec021), `bautzen_rects_asc_ctrl_500m_east.txt` (54da6d5f1d8d9cbecdfa62d04b46dc9390f7738a45002e5f60febd07700dc2c3), `bautzen_rects_ctrl_500m_east.txt` (1653ba18064c91815910fd8a980e2b50b5c2cdc3941b09e620736e097c294751) |
| states | RS: until 2025-05-12 10:00 local (before DS1 was installed); DS1: 2025-05-12 10:00 .. 2025-09-29 10:00 local; DS2: from 2025-09-29 10:00 local |
| cuts | DS1 = 2025-05-12T10:00:00, DS2 = 2025-09-29T10:00:00 |
| definition | `A = |sum z| / sqrt(N * sum |z|^2), phi = arg(sum z), over the deck-edge mask pixels` |
| differential | `dphi = wrap(phi_bridge - phi_control), paired over identical overpass times, over the deck-edge mask pixels` |

| series | file | dates plotted |
| --- | --- | --- |
| ASC Brücke | `bautzen_rects_asc_iw2.txt` | 36 |
| DESC Brücke | `bautzen_rects_desc_iw3.txt` | 35 |
| ASC Kontrolle 500 m östlich | `bautzen_rects_asc_ctrl_500m_east.txt` | 36 |
| DESC Kontrolle 500 m östlich | `bautzen_rects_ctrl_500m_east.txt` | 35 |

| geometry | bridge | control 500 m east | overpasses (bridge/control) | paired | unpaired |
| --- | --- | --- | --- | --- | --- |
| ASC | ASC Brücke | ASC Kontrolle 500 m östlich | 36/36 | 28 | 8/0 |
| DESC | DESC Brücke | DESC Kontrolle 500 m östlich | 35/35 | 35 | 0/0 |

A pair needs a valid deck-edge mask on *both* overpasses; if one rectangle has none (and therefore no committed `gamma2_band_raw` either), that date is left out of the difference. `unpaired` lists bridge-only / control-only such dates — here all of them are ASC-control dates whose mask is empty.

## 3. Figure

`figures/bautzen/fig_bautzen_ocv_amp_phase.png` — three rows × one column per series. Top row: the coherence amplitude `A(t)` in [0, 1]. Second row: the deck-edge phase `phi(t)` in (-pi, pi]. Third row: the paired differential phase `dphi(t)` in (-pi, pi]. Per-overpass dots are coloured by state (RS / DS1 / DS2), a thin line connects the series in time order, the state windows are the background bands, and the two cut dates are the dashed vertical lines. The bridge columns are the signal, the two control columns the counter-check. The `dphi` panel belongs to the *geometry*, so it is drawn once per geometry under its bridge column (ASC = first, DESC = second); the control columns of that row keep the state bands and a pointer, but carry no series of their own.

## 4. Per-state summary and cross-check

Circular statistics for the phase (arithmetic means would be wrong across the ±pi branch cut): `mean` is the circular mean, `R` its resultant length (0 = uniform angles, 1 = identical angles).

| series | state | n | median A | mean A | phase mean [rad] | R |
| --- | --- | --- | --- | --- | --- | --- |
| ASC bridge | RS (reference) | 7 | 0.2283 | 0.2801 | -1.1664 | 0.1638 |
| ASC bridge | DS1 | 23 | 0.2709 | 0.2641 | 0.8872 | 0.1426 |
| ASC bridge | DS2 | 6 | 0.3830 | 0.4001 | 0.7652 | 0.4186 |
| DESC bridge | RS (reference) | 7 | 0.4063 | 0.4170 | 1.5706 | 0.3719 |
| DESC bridge | DS1 | 23 | 0.2190 | 0.2306 | 1.5849 | 0.1287 |
| DESC bridge | DS2 | 5 | 0.1652 | 0.2660 | -1.2214 | 0.2697 |
| ASC control 500 m east | RS (reference) | 7 | 0.6460 | 0.5629 | 0.7669 | 0.2564 |
| ASC control 500 m east | DS1 | 23 | 0.7857 | 0.7145 | 1.9100 | 0.1172 |
| ASC control 500 m east | DS2 | 6 | 0.4678 | 0.4867 | 0.9709 | 0.3695 |
| DESC control 500 m east | RS (reference) | 7 | 0.1747 | 0.1903 | 1.8881 | 0.1642 |
| DESC control 500 m east | DS1 | 23 | 0.3378 | 0.3653 | -1.7048 | 0.1543 |
| DESC control 500 m east | DS2 | 5 | 0.3611 | 0.4082 | 0.2020 | 0.2108 |

Paired differential phase `dphi(t) = wrap(phi_bridge - phi_control)`, per geometry and state (circular statistics as above; `R` is close to 0 because the differential phase is speckle-dominated and the pair counts per state are small):

| geometry | state | n pairs | dphi mean [rad] | R |
| --- | --- | --- | --- | --- |
| ASC | RS (reference) | 6 | 1.8935 | 0.4449 |
| ASC | DS1 | 16 | -2.6980 | 0.4684 |
| ASC | DS2 | 6 | 2.3308 | 0.1838 |
| DESC | RS (reference) | 7 | 0.0340 | 0.6847 |
| DESC | DS1 | 23 | -1.8984 | 0.1669 |
| DESC | DS2 | 5 | 2.9106 | 0.4415 |

### Δφ between the states (Δφ_DS1−RS, Δφ_DS2−RS)

The state windows are consecutive and disjoint, so there is **no per-date pairing between states**; the shift away from the reference state RS is the difference of the two state circular means, `delta = wrap(circmean(state) - circmean(RS))`. Because the circular mean is not linear this is **not** the same as the difference of the series' own state phases, and both are given: first the figure's paired quantity `dphi`, then the absolute `phi` per series as counter-check. The interval is a seeded percentile bootstrap (20 000 resamples, 2.5/97.5 percentiles); `resolvable` marks a delta whose interval excludes 0.

| geometry | delta | Δφ [rad] | 95% CI [rad] | resolvable |
| --- | --- | --- | --- | --- |
| ASC | Δφ_DS1−RS | 1.6917 | [-2.91, +2.97] | no |
| ASC | Δφ_DS2−RS | 0.4373 | [-2.37, +2.73] | no |
| DESC | Δφ_DS1−RS | -1.9324 | [-2.97, +2.93] | no |
| DESC | Δφ_DS2−RS | 2.8766 | [-3.08, +3.08] | no |

| series | delta | Δφ [rad] | 95% CI [rad] | resolvable |
| --- | --- | --- | --- | --- |
| ASC bridge | Δφ_DS1−RS | 2.0536 | [-3.01, +3.00] | no |
| ASC bridge | Δφ_DS2−RS | 1.9316 | [-2.99, +3.02] | no |
| DESC bridge | Δφ_DS1−RS | 0.0143 | [-2.88, +2.84] | no |
| DESC bridge | Δφ_DS2−RS | -2.7920 | [-3.07, +3.07] | no |
| ASC control 500 m east | Δφ_DS1−RS | 1.1431 | [-2.94, +2.95] | no |
| ASC control 500 m east | Δφ_DS2−RS | 0.2039 | [-2.78, +2.71] | no |
| DESC control 500 m east | Δφ_DS1−RS | 2.6903 | [-3.04, +3.03] | no |
| DESC control 500 m east | Δφ_DS2−RS | -1.6861 | [-2.99, +2.98] | no |

**Verdict: not resolvable.** Every interval above spans almost the whole circle, so no shift from RS to DS1 or DS2 survives the per-state noise (5–23 overpasses per state, resultant length `R` ≈ 0.1–0.7, which puts the uncertainty of a single state mean at order 1 rad). The value also depends on the definition — for ASC DS1 the paired panel gives `+1.69 rad` while `Δφ_bridge − Δφ_control` gives `+0.91 rad` — which is exactly why the definition is fixed above and why neither is reported as an effect.

Cross-checks of the recomputed mask, coherence and pairing against the committed table: **515 checks, 0 failures**. Each date is checked for its deck-edge pixel count (`n_masked`), its band contrast, and the identity `A(t)^2 == gamma2_band_raw` (tolerance 1e-9); each bridge/control pair is checked to cover the same overpasses and to satisfy the wrap identity `wrap_pi(dphi - (phi_bridge - phi_control)) == 0` (tolerance 1e-12); each state delta is checked to equal `wrap_pi(circmean(state) - circmean(RS))` and to lie inside its own bootstrap interval.

**Failures: none.** Every recomputed deck-edge pixel count, band contrast and coherence value reproduces the committed table, every paired date satisfies the wrap identity, and every state delta matches its circular-mean definition and its bootstrap interval.

## 5. Reproduce

```sh
cd <repo root>
python3 code/bautzen/fig_bautzen_ocv_amp_phase.py
PYTHON=/usr/bin/python3 bash figures/bautzen/fig_bautzen_ocv_amp_phase.sh
```

Needs numpy and matplotlib. This report was generated with python 3.13.5, numpy 2.4.6, matplotlib 3.11.1. No output carries a timestamp, so repeated runs are identical.

## 6. Caveats

- **Amplitude = coherence amplitude.** `A(t)` is the square root of the `gamma2_band_raw` channel, so it inherits its limits: the raw estimator is used here (not floor-corrected) and the deck-edge pixel count `N` is small (single digits to a few tens), which widens every per-state reading.
- **Phase without a reference date.** The chips are single SLCs, so `phi(t)` is the *absolute* phase of the coherent deck-edge sum, dominated by the mean backscatter phase of the deck line. It is a time series of that phase, not unwrapped and not referenced to a master date; a real displacement would only appear as a deviation from this baseline.
- **State and season are confounded.** RS / DS1 / DS2 are consecutive time intervals of one deployment, assigned by overpass time (`core.state_of`); the control series and the cut dates are the only counter-checks carried here.
- **The differential phase is not InSAR phase.** `dphi(t)` differences the `arg` of two single-SLC coherent sums; it is defined modulo 2*pi and is *not* an interferometric (master-slave) phase. What it does remove is the component common to the bridge and its control rectangle, so a seasonal or atmospheric swing cancels and a local deck-edge change survives; the residual is still dominated by speckle and by the small pixel count `N`.
- **The state-to-state Δφ is not resolvable.** The shift `wrap(circmean(state) - circmean(RS))` is reported with a seeded bootstrap interval (section 4), but every interval spans almost the whole circle: with 5–23 overpasses per state and `R` ≈ 0.1–0.7 the circular mean of one state carries an uncertainty of order 1 rad. The delta is documented, not claimed, and it is not additive — the circular mean is not linear.
- **No committed pin.** The reference JSONs of the package hold no amplitude or phase, so the only verification is the internal identity `A(t)^2 == gamma2_band_raw` (section 4); the figure does not touch the figures A-E pin block.

---

Generated by `code/bautzen/fig_bautzen_ocv_amp_phase.py` (515 cross-checks, 0 failures).
