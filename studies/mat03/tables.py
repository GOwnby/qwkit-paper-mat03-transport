"""MAT-03 Table 1 rows (clean BSO) -> ../../paper_mat03/table1.tex"""
import json, os
here = os.path.dirname(os.path.abspath(__file__))
J = lambda f: json.load(open(os.path.join(here, f)))
HA = J("hall_algebra.json"); MR = J("mr.json"); N = J("numbers.json"); L = J("lee_fig3b_v2.json"); CP = J("coupled.json")
fm = lambda v, f: "--" if v is None else format(v, f)
rows = []
for vg in ("4V", "5V", "7V", "8V"):
    nt = L["n_1e13"][vg][-1]; r2 = L["n_1e13"][vg][0] / nt
    f1 = [h for h in HA if h["NA"] == 1e18 and h["vg"] == vg and h["T"] == 2][0]["f1"]
    if vg in MR["sets"]:
        s3 = MR["sets"][vg]["3"]
        slow = f"{s3['0.5']['min_slow']['200.0']:.2f} / {s3['1.0']['min_slow']['200.0']:.2f}"
    else:
        slow = "--"
    rr = [c for c in CP if c["vg"] == vg and c["matched"]]
    if rr:
        v = [c[k]["nH0"] for c in rr for k in ("mult", "add")]
        cou = f"{min(v):.2f}--{max(v):.2f}"
    else:
        cou = "--"
    d1 = [r for r in N["dmin"] if r["NA"] == 1e18 and r["vg"] == vg and r["eps_cap"] == 24.0 and r["Z"] == 1][0]
    d2 = [r for r in N["dmin"] if r["NA"] == 1e18 and r["vg"] == vg and r["eps_cap"] == 24.0 and r["Z"] == 2][0]
    dm = lambda r, k: "$\\leq 1$" if r[k] is not None and r[k] <= 1.0 else fm(r[k], ".1f")
    xm = "50" if d1["xmax150"] >= 49.9 else f"{d1['xmax150']:.1f}"
    rows.append(f"{vg[0]} & {nt:.2f} & {r2:.2f} & {f1:.2f} & {slow} & {cou} & {dm(d1, 'dmin_hall')} / {dm(d1, 'dmin_3a')}"
                f" & {dm(d2, 'dmin_hall')} / {dm(d2, 'dmin_3a')} & {xm} \\\\")
open(os.path.join(here, "..", "..", "paper_mat03", "table1.tex"), "w").write("\n".join(rows) + "\n")
print("\n".join(rows))
