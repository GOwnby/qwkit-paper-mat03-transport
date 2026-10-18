"""Build the verified data set lee_fig3b_v2.json for MAT-03 from the independent digitizations of Lee et al.,
Adv. Electron. Mater. 11, 2400811 (2024):
  dig3b_v2.json  marker centroids in Fig. 3b (dig3b.py, 600 dpi crop)
  dig4.json      marker centroids in Fig. 4 (mu versus n at 10 K and 150 K; dig4.py)
  ../lee_fig3b.json  the first, manual digitization (used where the centroid colours of 2 V and 3 V merge)
Each value carries its source. The figure crops are made from the user's copy of the paper and are not distributed."""
import json, os
here = os.path.dirname(os.path.abspath(__file__))
J = lambda f: json.load(open(os.path.join(here, f)))
C = J("dig3b_v2.json"); F4 = J("dig4.json"); OLD = J(os.path.join("..", "lee_fig3b.json"))
T = [2, 10, 50, 100, 150]
# Fig. 4 markers by gate voltage (ordered as in dig4.py output; 4 V and 5 V at 10 K split by 2-means)
f4 = {"10": dict(zip(["0V", "1V", "2V", "3V", "4V", "5V", "7V", "6V", "8V"], F4["10K"])),
      "150": dict(zip(["0V", "1V", "2V", "3V", "4V", "5V", "6V", "7V", "8V"], F4["150K"]))}
n, mu, src = {}, {}, {}
for v in [f"{i}V" for i in range(9)]:
    n[v], mu[v], src[v] = [], [], []
    for j, t in enumerate(T):
        c_n = C[v]["n"].get(str(t)); c_m = C[v]["mu"].get(str(t))
        if v in ("2V", "3V"):                       # colours merge in the centroid pass: keep the manual values
            n[v].append(OLD["n_1e13"][v][j]); mu[v].append(OLD["mu"][v][j]); src[v].append("manual"); continue
        if v == "6V":                               # hidden behind 7 V at 2 K; Fig. 4 at 10 and 150 K
            if t == 2: n[v].append(None); mu[v].append(None); src[v].append("hidden")
            elif str(t) in f4:
                p = f4[str(t)][v]; n[v].append(round(p[0] / 1e13, 3)); mu[v].append(p[1]); src[v].append("fig4")
            else: n[v].append(round(c_n[0] / 1e13, 3)); mu[v].append(round(c_m[0], 1)); src[v].append("centroid, +-5%")
            continue
        if v == "4V" and t == 10:                   # mobility marker hidden in Fig. 3b
            p = f4["10"]["4V"]; n[v].append(round(c_n[0] / 1e13, 3)); mu[v].append(p[1]); src[v].append("centroid n, fig4 mu"); continue
        n[v].append(round(c_n[0] / 1e13, 3)); mu[v].append(round(c_m[0], 1)); src[v].append("centroid")
# cross-check against Fig. 4 where both exist
check = {}
for t in ("10", "150"):
    for v, p in f4[t].items():
        j = T.index(int(t))
        if n[v][j] is None or src[v][j].startswith("fig4"): continue
        check[f"{v}@{t}K"] = dict(dn_pct=round(100 * (n[v][j] * 1e13 / p[0] - 1), 2), dmu_pct=round(100 * (mu[v][j] / p[1] - 1), 2))
out = dict(T=T, n_1e13=n, mu=mu, source=src, check_vs_fig4=check,
           note="5 V verified by Fig. 3b centroids, Fig. 4 and the Fig. 5a overlay; Lee's text gives 8 V: 1999 cm2/Vs and "
                "6.52e13 at 2 K, 2102 and 7.04e13 at 10 K, 1262 and 1.08e14 at 150 K.")
json.dump(out, open(os.path.join(here, "..", "lee_fig3b_v2.json"), "w"), indent=1)
for k, v in check.items(): print(k, v)
for v in n: print(v, n[v], mu[v])
