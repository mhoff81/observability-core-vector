# Phase 0 audit — Champlain Towers South

Events: kollaps=2021-06-24, abriss=2021-07-05

| track | role | dates | platforms | windows | shift median (y,x) | max \|shift\| | clipped | pairs used | pairs skipped | platform mix |
|---|---|---|---|---|---|---|---|---|---|---|
| `cte_asc_48_iw3` | control_cte | 178 | S1A, S1B | UNKNOWN, W2_prae, W4_kollaps_bis_abriss, W5_geraeumt | (-0.01, -0.03) | 4.60 | 0 | 175 | 2 | S1A->S1A=171, S1A->S1B=2, S1B->S1A=2 |

Reading the audit: a fixed ground footprint shows a median registration
shift of a few pixels at most. A large or bimodal shift means the series is
**not** a stable footprint, and any later temporal signal must be read with
that in mind (the Morandi lesson: a 150 px block offset produced the very
coherence jump that looked like an event). Pairs that straddle a documented
state change are excluded and listed in `pairs_skipped`. Cross-platform pairs
are **kept and labelled** in `platform mix`: an S1A/S1C pair is a valid 6-day
interferogram (the two spacecraft share the orbit plane), and dropping it
would destroy the series in the S1C era — while pooling it unlabelled with
same-platform pairs would hide a platform effect (EQS-10 stratifies on it).
