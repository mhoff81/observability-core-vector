# Acquisition census — Champlain Towers South

Step 0 (metadata only) · candidate `cts` · 25.87306 / -80.12083 · box ±1 km

Query: 2014-10-01T00:00:00Z .. 2026-09-21T11:57:48Z · orbits ASCENDING, DESCENDING · rel-orbit filter any · 300 bursts seen

Geometry source: en.wikipedia 25.87306/-80.12083 (8777 Collins Ave); footprint no longer in OSM (demolished) - take it from county GIS or as-built plans

## Windows

| label | start | end | status | countable | source |
|---|---|---|---|---|---|
| W2_prae | 2014-10-01T00:00:00Z | 2021-06-03T00:00:00Z | given | true | archive start .. before the NIST initiation window |
| W3_initiierung | 2021-06-03T00:00:00Z | 2021-06-24T05:22:00Z | to_pin | false | NIST: 'early June 2021'; exact start day to pin from the report |
| W4_kollaps_bis_abriss | 2021-06-24T05:22:00Z | 2021-07-05T02:00:00Z | given | true | collapse 01:22 EDT (05:22 UTC) .. controlled demolition ~22:00 EDT on 2021-07-04 (hour to pin) |
| W5_geraeumt | 2021-07-05T02:00:00Z | 2022-05-01T00:00:00Z | given | true | debris cleared, site sold in May 2022 |
| W6_neubau | 2022-05-01T00:00:00Z | 2026-09-21T00:00:00Z | to_pin | false | Damac redevelopment; construction start to pin |

## Strata

| stratum (orbit\|platform\|rel\|subswath) | n | first | last | median gap (d) | max gap (d) | pairs | EQS-7 floor |
|---|---|---|---|---|---|---|---|
| ASCENDING|S1A|48|IW3 | 282 | 2015-09-21T23:27:51Z | 2026-06-19T23:28:15Z | 12 | 180 | 271 | yes |
| ASCENDING|S1B|48|IW3 | 2 | 2016-10-15T23:27:16Z | 2019-09-06T23:27:34Z | 1056 | 1056 | 0 | no |
| ASCENDING|S1D|48|IW3 | 7 | 2026-06-26T23:27:34Z | 2026-09-06T23:27:38Z | 12 | 12 | 6 | no |
| DESCENDING|S1A|84|IW2 | 5 | 2016-10-12T11:18:15Z | 2022-11-10T11:18:51Z | 528 | 828 | 0 | no |
| DESCENDING|S1B|84|IW2 | 1 | 2019-09-09T11:17:51Z | 2019-09-09T11:17:51Z | 0 | 0 | 0 | no |

## Window counts per stratum

| stratum | window | n |
|---|---|---|
| ASCENDING|S1A|48|IW3 | W2_prae | 132 |
| ASCENDING|S1A|48|IW3 | W4_kollaps_bis_abriss | 1 |
| ASCENDING|S1A|48|IW3 | W5_geraeumt | 24 |
| ASCENDING|S1B|48|IW3 | W2_prae | 2 |
| ASCENDING|S1B|48|IW3 | W4_kollaps_bis_abriss | 0 |
| ASCENDING|S1B|48|IW3 | W5_geraeumt | 0 |
| ASCENDING|S1D|48|IW3 | W2_prae | 0 |
| ASCENDING|S1D|48|IW3 | W4_kollaps_bis_abriss | 0 |
| ASCENDING|S1D|48|IW3 | W5_geraeumt | 0 |
| DESCENDING|S1A|84|IW2 | W2_prae | 4 |
| DESCENDING|S1A|84|IW2 | W4_kollaps_bis_abriss | 0 |
| DESCENDING|S1A|84|IW2 | W5_geraeumt | 0 |
| DESCENDING|S1B|84|IW2 | W2_prae | 1 |
| DESCENDING|S1B|84|IW2 | W4_kollaps_bis_abriss | 0 |
| DESCENDING|S1B|84|IW2 | W5_geraeumt | 0 |

## Brackets per stratum

| stratum | window | last before | first after | gap before (h) | gap after (h) |
|---|---|---|---|---|---|
| ASCENDING|S1A|48|IW3 | W2_prae | — | 2021-06-09T23:28:23Z | null | 167.47 |
| ASCENDING|S1A|48|IW3 | W4_kollaps_bis_abriss | 2021-06-21T23:28:23Z | 2021-07-15T23:28:25Z | 53.89 | 261.47 |
| ASCENDING|S1A|48|IW3 | W5_geraeumt | 2021-07-03T23:28:24Z | 2022-05-11T23:28:27Z | 26.53 | 263.47 |
| ASCENDING|S1B|48|IW3 | W2_prae | — | — | null | null |
| ASCENDING|S1B|48|IW3 | W4_kollaps_bis_abriss | 2019-09-06T23:27:34Z | — | 15749.91 | null |
| ASCENDING|S1B|48|IW3 | W5_geraeumt | 2019-09-06T23:27:34Z | — | 16010.54 | null |
| ASCENDING|S1D|48|IW3 | W2_prae | — | 2026-06-26T23:27:34Z | null | 44399.46 |
| ASCENDING|S1D|48|IW3 | W4_kollaps_bis_abriss | — | 2026-06-26T23:27:34Z | null | 43629.46 |
| ASCENDING|S1D|48|IW3 | W5_geraeumt | — | 2026-06-26T23:27:34Z | null | 36431.46 |
| DESCENDING|S1A|84|IW2 | W2_prae | — | 2022-11-10T11:18:51Z | null | 12611.31 |
| DESCENDING|S1A|84|IW2 | W4_kollaps_bis_abriss | 2020-08-04T11:18:37Z | 2022-11-10T11:18:51Z | 7770.06 | 11841.31 |
| DESCENDING|S1A|84|IW2 | W5_geraeumt | 2020-08-04T11:18:37Z | 2022-11-10T11:18:51Z | 8030.69 | 4643.31 |
| DESCENDING|S1B|84|IW2 | W2_prae | — | — | null | null |
| DESCENDING|S1B|84|IW2 | W4_kollaps_bis_abriss | 2019-09-09T11:17:51Z | — | 15690.07 | null |
| DESCENDING|S1B|84|IW2 | W5_geraeumt | 2019-09-09T11:17:51Z | — | 15950.7 | null |

## Missing

- W3_initiierung: window boundary to pin (NIST: 'early June 2021'; exact start day to pin from the report)
- W6_neubau: window boundary to pin (Damac redevelopment; construction start to pin)

## Warnings

- discovery mode: 2 relative orbits found — the non-suitable ones must be justified in advance (EQS-9)
