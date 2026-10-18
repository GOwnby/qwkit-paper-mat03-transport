# MAT-03 data verification (15 Oct 2026)

Source: Lee et al., Adv. Electron. Mater. 11, 2400811 (2024), DOI 10.1002/aelm.202400811 (the user's copy).
Figure crops (600 dpi) are made from that copy and are not distributed here.

## Routes
1. `../lee_fig3b.json`: the first, manual digitization of Fig. 3b (used up to 14 Oct 2026).
2. `dig3b.py` -> `dig3b_v2.json`: marker centroids per colour in Fig. 3b (colour distance < 45).
   The 2 V and 3 V colours merge; the 4 V mobility marker at 10 K and the 6 V markers at 2 K are hidden.
3. `dig4.py` -> `dig4.json`: marker centroids in Fig. 4 (mu_H against n_H at 10 K and 150 K). Overlapping
   4 V / 5 V markers at 10 K are split by 2-means on pixel coordinates.
4. `dense5a.py` -> `dense5a.json`: line centre per pixel column in Fig. 5a (MR at 2 K), parabolic fit over 2-7 T.
5. Lee's text: 8 V gives 1999 cm2/Vs and 6.52e13 cm-2 at 2 K, 2102 and 7.04e13 at 10 K, 1262 and 1.08e14 at 150 K.

`build_v2.py` merges them into `../lee_fig3b_v2.json` (each value with its source), which all MAT-03 scripts now read.

## 5 V (the voltage asked about)
| T (K) | n_H manual | n_H centroid | n_H Fig. 4 | mu_H manual | mu_H centroid | mu_H Fig. 4 |
|---|---|---|---|---|---|---|
| 2   | 4.151 | 4.150 | - | 420.3 | 421.9 | - |
| 10  | 4.276 | 4.320 | 4.336 | 409.6 | 412.5 | 406.6 |
| 50  | 4.676 | 4.669 | - | 409.6 | 411.5 | - |
| 100 | 4.818 | 4.831 | - | 379.0 | 382.4 | - |
| 150 | 4.890 | 4.903 | 4.897 | 324.5 | 330.4 | 328.0 |
(n_H in 1e13 cm-2, mu_H in cm2/Vs.) All three routes agree within 1.5%; n_H(2 K)/n_H(150 K) = 0.846-0.849.
MR at 5 V, 2->7 T: -0.07 points (first reading) and -0.117 +- 0.004 (dense fit); the 0.5 / 1.0-point tolerances
of mr.py cover both. MR overlay (digitized points on the figure) confirms the 5 V and 7 V curves.

## Cross-check against Fig. 4, all voltages (centroid or v2 value minus Fig. 4)
4-8 V: n_H within 0.6%, mu_H within 1.5%. 0-1 V: within 0.5% / 5%. 2-3 V at 150 K: the Fig. 4 markers are partly
hidden (small marker areas), differences up to 5% / 11%; the manual values are kept.

## Corrections found
- 4 V mu_H at 10 K: 517 (manual) -> 435.3 (Fig. 4; hidden in Fig. 3b). Changes the 4 V Hall resistance ratio
  R(2K)/R(50K) from 0.579 to 0.612 and the mobility rise 50->2 K from +86% to +76%.
- 6 V: 2 K hidden; 10 K and 150 K now from Fig. 4 (n_H 5.693 and 6.912e13); 50 and 100 K centroids (+-5%).
- 7 V and 8 V densities rise by 0.3-0.9% (centroids); 8 V at 2 K: 6.497e13 and 2010 cm2/Vs (text 6.52e13, 1999).
- MR 7 V 2->7 T: -0.08 (first) vs -0.060 +- 0.005 (dense); within tolerance.
All MAT-03 numbers were regenerated with lee_fig3b_v2.json (outputs before: ../pre_digv2/).
