"""MAT-03: every number quoted in the manuscript, from hall_algebra/consistency/mr/coulomb/roughness/conv3(.json) and the
MAT-02 limits (../mat02/gate_curves.json at 300 K, ../mat02/gate_map.json at 150 K). Writes numbers.log, numbers.json."""
import json, os
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
J = lambda f: json.load(open(os.path.join(here, f)))
L = J("lee_fig3b_v2.json"); HA = J("hall_algebra.json"); C = J("consistency.json"); MR = J("mr.json"); CO = J("coulomb.json")
A3 = J("lee_fig3a_Rs.json")
GC = J("../mat02/gate_curves.json"); GM = J("../mat02/gate_map.json")
out = dict(dmin=[], window=[]); lines = []
P = lambda s: (lines.append(s), print(s))
QE = 1.602176634e-19


def xmax(NA, eps, ns, T):
    if T == 300:
        rr = [g for g in GC if g["NA"] == NA and g["eps_cap"] == eps][0]["rows"]
        return float(np.interp(ns, [r["ns"] for r in rr], [r["x_max"] for r in rr]))
    rr = [g for g in GM if g["NA"] == NA and g["eps_cap"] == eps and g["T"] == 150.0][0]["sheet"]
    return float(np.interp(ns, [r["ns"] for r in rr], [r["x_max"] for r in rr]))


def bounds(vg):
    """Largest sigma / (e n_tot) over 2-150 K: from the Hall data, and from Fig. 3a (5-8 V; 50 K minimum of R_s)."""
    nt = L["n_1e13"][vg][-1] * 1e13
    bh = max(n * 1e13 * m for n, m in zip(L["n_1e13"][vg], L["mu"][vg]) if n is not None and m is not None) / nt
    b3 = None
    if vg in A3["Rs"]:
        Rmin = min(R for t, R in zip(A3["T"], A3["Rs"][vg]) if t <= 150)
        b3 = 1 / (QE * Rmin * nt)
    return bh, b3


P("== data (Lee et al. Fig. 3b): n_H(T) relative to 150 K, mu_H, sigma/e")
for vg in [f"{v}V" for v in range(9)]:
    n = L["n_1e13"][vg]; m = L["mu"][vg]
    if n[0] is not None:
        P(f"{vg}: n_H 150 K {n[-1]:.3f}e13, 50 K {n[2]:.3f}, 2 K {n[0]:.3f} | n_H(2K)/n_H(150K) {n[0] / n[-1]:.3f}, "
          f"n_H(50K)/n_H(150K) {n[2] / n[-1]:.3f} | loss 150->2 K {(n[-1] - n[0]):.2f}e13, 50->2 K {(n[2] - n[0]):.2f}e13 | "
          f"mu_H 150/50/2 K {m[-1]:.0f}/{m[2]:.0f}/{m[0]:.0f}")
P("== resistance upturn 50 -> 2 K (consistency.json)")
for u in C["upturn"]:
    s = f"{u['vg']}:"
    for lab in ("fig3a", "hall"):
        if u[lab] is not None:
            s += f" {lab} R(2)/R(50) {u[lab]['ratio']:.3f} A {u[lab]['A']:.0f} kFl {u[lab]['kFl']:.0f};"
    if "dnH" in u:
        s += f" dn_H {100 * u['dnH']:+.1f}% dmu_H {100 * u['dmuH']:+.1f}%"
    if "scale_3a_over_hall" in u:
        s += f" | 3a/Hall {min(u['scale_3a_over_hall']):.2f}-{max(u['scale_3a_over_hall']):.2f}"
    P(s)
P("== magnetoresistance at 2 K (mr.json)")
for vg, w in MR["wl"].items():
    P(f"{vg}: plateau {w['plateau']:.2f}%, change 2->7 T {w['change_2_7']:+.2f}%, R(2K)/R(50K)-1 (Hall) {w['upturn_hall']:.1f}%, "
      f"sigma drop 50->2 K (Hall) {w['sigma_drop_hall']:.1f}%")
for vg, sets in MR["sets"].items():
    for k, rk in sets.items():
        for tol in ("0.5", "1.0"):
            g = rk[tol]
            P(f"{vg} {k} channels (accepted {rk['accepted']}), tol {tol}: flat-MR sets {g['n']}; smallest weight below 25/50/100/200: "
              + "/".join("-" if v is None else f"{v:.3f}" for v in g["min_slow"].values())
              + (f"; slowest mu {g['slowest_mu'][0]:.0f}-{g['slowest_mu'][1]:.0f}, its weight {g['slowest_w'][0]:.3f}-{g['slowest_w'][1]:.3f}" if g["slowest_mu"] else ""))
P("== Hall bend under the two-subband reading (consistency.json, straight-line fits, 2 and 50 K)")
for c in C["linearity"]:
    if c["T"] in (2, 50) and c.get("x") is not None and c["x"] != 1.0:
        P(f"NA {c['NA']:.0e} {c['vg']} {c['T']} K: x {c['x']:.2f}, R_xy off line {100 * c['dev']:.1f}%")
P("== Coulomb scattering by the gate-induced charge (coulomb.json)")
for g in CO:
    for Z in (1, 2):
        rr = [r for r in g["rows"] if r["Z"] == Z]
        dd = np.array([r["d"] for r in rr]); ms = np.array([r["mu_sigma"] for r in rr])
        xs = [r["x"] for r in rr]
        bh, b3 = bounds(g["vg"])
        dmin = lambda b: (float(np.interp(np.log(b), np.log(ms), dd)) if ms[0] < b else float(dd[0]))
        xm = {T: xmax(g["NA"], g["eps_cap"], g["n_tot"], T) for T in (300, 150)}
        row = dict(NA=g["NA"], vg=g["vg"], eps_cap=g["eps_cap"], Z=Z, bound_hall=bh, bound_3a=b3, dmin_hall=dmin(bh),
                   dmin_3a=None if b3 is None else dmin(b3), xmax300=xm[300], xmax150=xm[150])
        mux = float(np.exp(np.interp(xm[150], dd, np.log(ms))))
        row["supp_hall"] = max(bh / mux, 1.0); row["supp_3a"] = None if b3 is None else max(b3 / mux, 1.0)
        if xs[0] is not None:
            # x over the window between d_min (Hall bound) and x_max (150 K), or at x_max if the window is empty
            lo, hi = row["dmin_hall"], xm[150]
            sel = [r["x"] for r in rr if lo <= r["d"] <= hi] or [float(np.interp(hi, dd, xs))]
            row["x_window"] = (min(sel), max(sel)); row["x_d1to6"] = (min(xs[:5]), max(xs[:5]))
            row["nH_window"] = min(r["nH"] / g["n_tot"] for r in rr if lo <= r["d"] <= hi) if any(lo <= r["d"] <= hi for r in rr) else None
        out["dmin"].append(row)
        fm = lambda v, f: "-" if v is None else format(v, f)
        msg = (f"NA {g['NA']:.0e} {g['vg']} eps_cap {g['eps_cap']:.0f} Z {Z}: bound {bh:.0f} (Hall) / {fm(b3, '.0f')} (Fig. 3a)"
               f" -> d_min {row['dmin_hall']:.1f} / {fm(row['dmin_3a'], '.1f')} A | x_max {xm[150]:.1f} (150 K) {xm[300]:.1f} (300 K)"
               f" | suppression needed {row['supp_hall']:.2f} / {fm(row['supp_3a'], '.2f')}")
        if "x_window" in row:
            msg += (f" | x(d 1-6 A) {row['x_d1to6'][0]:.2f}-{row['x_d1to6'][1]:.2f}, x(window) {row['x_window'][0]:.2f}-"
                    f"{row['x_window'][1]:.2f}")
        P(msg)
    if g["dislocations"]:
        xs = [r["x"] for r in g["dislocations"] if r["x"] is not None]
        if xs:
            P(f"   dislocations 1e10-1e11 (d 4 A, eps 24): x {min(xs):.2f}-{max(xs):.2f}")
P("== Coulomb scattering at the measured conductance, coupled magnetotransport (coupled.json; matched rows, d >= d_min)")
CP = J("coupled.json")
for vg in ("5V", "7V", "8V"):
    rr = [c for c in CP if c["vg"] == vg and c["matched"]]
    for lab in ("mult", "add"):
        nh = [c[lab]["nH0"] for c in rr]; nf = [c[lab]["nHfit"] for c in rr]; mr = [c[lab]["mr27"] for c in rr]
        P(f"{vg} {lab}: n_H/n_tot {min(nh):.3f}-{max(nh):.3f} (fit {min(nf):.3f}-{max(nf):.3f}); MR(7T)-MR(2T) {100 * min(mr):.2f}-{100 * max(mr):.2f}%"
          f" over d {min(c['d'] for c in rr):.0f}-{max(c['d'] for c in rr):.0f} A")
if os.path.exists(os.path.join(here, "roughness.json")):
    P("== interface roughness with the Coulomb scattering (roughness.json)")
    for g in J("roughness.json"):
        xs = [r["x"] for r in g["rows"] if r["Lambda"] > 0]
        P(f"{g['vg']}: x {min(xs):.2f}-{max(xs):.2f}; <z> {[round(v, 1) for v in g['zmean']]} A; mu_sigma (Delta, Lambda): "
          + ", ".join(f"({r['Delta']:.0f},{r['Lambda']:.0f}) {r['mu_sigma']:.0f}" for r in g["rows"]))
if os.path.exists(os.path.join(here, "conv3.json")):
    P("== convergence (conv3.json)")
    for c in J("conv3.json"):
        P(f"d {c['d']} dz {c['dz']} nstates {c['nstates']} nth {c['nth']}: mu_sigma {c['mu_sigma']:.1f}")
json.dump(out, open(os.path.join(here, "numbers.json"), "w"), indent=1)
open(os.path.join(here, "numbers.log"), "w").write("\n".join(lines) + "\n")
