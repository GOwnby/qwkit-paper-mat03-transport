"""MAT-03 (9): can trapping at fixed total charge explain the Hall-density drop on cooling?
The free density is taken as the measured Hall density n_H(T) (the trapping reading: trapped electrons do not conduct),
the total charge N_tot >= n_H(150 K) is fixed, and N_tot - n_H(T) electrons sit in traps.

(A) Fermi level. E_F(n, T) of the interface model (common3.solve) for the free density, measured from three
    references a trap level can be tied to: the BSO conduction-band edge at the interface ("int"), the lowest subband
    ("E1"), and the model's fixed frame (far BSO pinned; "abs").
(B) Required trap entropy. A trap fixed relative to a reference fills on cooling between T1 < T2 only if its free
    energy relative to E_F rises with T at least as fast as E_F does: S_req = [E_F(T2) - E_F(T1)] / [k_B (T2 - T1)].
(C) Minimum trap density (linear programme). Traps of any level distribution D(E) >= 0, occupied
    f = 1 / (1 + exp((E - E_F) / kT) / g) with a constant degeneracy (entropy) factor g = exp(S/k_B), S from -40 to
    +5 k_B (negative S favours trapping at low T; the manuscript uses S_t = -S, the entropy gained on emptying). A trap whose barrier freezes it at T_f during cooling keeps the
    occupation it had at T_f until re-warmed past T_f (equal cooling and warming rates; cooling curve = warming curve);
    T_f = 2 K is equilibrium throughout. The LP finds the smallest total number of trap states reproducing n_H(T) at the
    measured temperatures within a tolerance. Variables scaled to 1e13 cm^-2; HiGHS status reported (2 = infeasible;
    numerical failures retried with other HiGHS methods and reported separately); 0.5 meV level grid, +-0.6 eV window.
(D) Release counting (the first version of this argument, kept for the record: it ignores the Fermi-level shift and
    is not a valid test). Writes trapping.json."""
import json
from multiprocessing import Pool
import numpy as np
from scipy.optimize import linprog
from common3 import *

KB = 8.617333262e-5
TF = [2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 50.0, 70.0, 100.0, 150.0]    # freezing temperatures (2 = equilibrium)
VGS = ("4V", "5V", "7V", "8V")
NAS = (5e19, 1e18)
REFS = ("int", "E1", "abs")
SS = (-40.0, -30.0, -25.0, -20.0, -15.0, -10.0, -5.0, -2.0, -0.69, 0.0, 0.69, 2.0, 5.0)
NGRID = np.geomspace(2.8e13, 1.2e14, 24)
SITES = 5.9e14                                   # unit cells per interface plane (a = 4.1 A)
U = 1e13                                         # LP unit (cm^-2)


def ef_job(args):
    NA, T = args
    p = dict(BASE); p["NA"] = NA; p["T"] = T; p["f_int"] = 1.0
    out = []
    for n in NGRID:
        r = solve(float(n), p, extend_to=0.0, nstates=10, dz=0.5)
        z = r["z"]; V0 = float(r["V"][np.nonzero(z > 0)[0][0]])
        out.append(dict(n=float(n), int=float(r["mu"] - V0), E1=float(r["mu"] - r["E"][0]), abs=float(r["mu"])))
    return dict(NA=NA, T=T, rows=out)


def nH_at(vg, T):
    """Measured n_H(T), log-T interpolation between the data temperatures (cooling curve = warming curve)."""
    t = np.log(TEMPS); y = np.array(NH[vg], float)
    ok = np.isfinite(y)
    return float(np.interp(np.log(T), t[ok], y[ok]))


def make_ef(table, NA, ref):
    rows = {row["T"]: row["rows"] for row in table if row["NA"] == NA}
    Ts = sorted(rows)
    def ef(n, T):
        vals = [np.interp(np.log(n), np.log(NGRID), [x[ref] for x in rows[t]]) for t in Ts]
        return float(np.interp(np.log(T), np.log(Ts), vals))
    return ef


def s_required(vg, ef):
    """Entropy (units of k_B) a fixed-level trap needs to fill on cooling across each interval of the data."""
    Tm = [T for T, n in zip(TEMPS, NH[vg]) if n is not None]
    e = [ef(nH_at(vg, T), T) for T in Tm]
    return [dict(T1=Tm[i], T2=Tm[i + 1], dEF_meV=1e3 * (e[i + 1] - e[i]), S_req=(e[i + 1] - e[i]) / (KB * (Tm[i + 1] - Tm[i])))
            for i in range(len(Tm) - 1)]


def _solve(c, A_ub, b_ub, nv):
    for method in ("highs", "highs-ipm", "highs-ds"):
        res = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=[(0, None)] * (nv + 1), method=method)
        if res.status in (0, 2):
            return res, method
    return res, method


def min_traps(vg, ef, g, tol, freezing=True, de=0.0005, window=0.6):
    Tm = [T for T, n in zip(TEMPS, NH[vg]) if n is not None]
    nm = np.array([nH_at(vg, T) for T in Tm]) / U
    efs = [ef(nH_at(vg, T), T) for T in TF]
    E = np.arange(min(efs) - window, max(efs) + window, de)
    classes = TF if freezing else [2.0]
    cols = []
    for Tf in classes:
        occ = np.zeros((len(Tm), len(E)))
        for i, T in enumerate(Tm):
            Tp = max(T, Tf)
            x = (E - ef(nH_at(vg, Tp), Tp)) / (KB * Tp)
            occ[i] = 1.0 / (1.0 + np.exp(np.clip(x, -700, 700)) / g)
        cols.append(occ)
    A = np.hstack(cols)                                            # trapped(T_i) = A @ D
    nv = A.shape[1]
    b = nm[-1] - nm                                                # trapped(T_i) = N_extra + n_H(150) - n_H(T_i)
    Aeq = np.hstack([A, -np.ones((len(Tm), 1))])
    A_ub = np.vstack([Aeq, -Aeq]); b_ub = np.concatenate([b + tol * nm, -(b - tol * nm)])
    c = np.concatenate([np.ones(nv), [0.0]])
    res, method = _solve(c, A_ub, b_ub, nv)
    if res.status == 2:
        return dict(status="infeasible")
    if res.status != 0:
        return dict(status=f"solver_failed({res.status})")
    D = res.x[:nv].reshape(len(classes), len(E)) * U
    tot = D.sum(0)
    w = int(round(0.05 / de)); dens = np.convolve(tot, np.ones(w), "same") / 0.05
    return dict(status="feasible", method=method, total=float(D.sum()), N_extra=float(res.x[-1] * U),
                peak_per_eV=float(dens.max()), E_peak=float(E[np.argmax(dens)]), ef_range=[float(min(efs)), float(max(efs))],
                frozen_share=float(D[1:].sum() / max(D.sum(), 1e-30)) if freezing else 0.0)


def release_count(vg):
    nut = np.log(1e12 * 1e3)
    out = []
    for (T1, n1), (T2, n2) in zip(zip(TEMPS, NH[vg]), list(zip(TEMPS, NH[vg]))[1:]):
        if n1 is None or n2 is None:
            continue
        dE = KB * nut * (T2 - T1)
        out.append(dict(T1=T1, T2=T2, Eb=[KB * nut * T1, KB * nut * T2], dn=n2 - n1, per_eV=(n2 - n1) / dE))
    return out


def lp_job(args):
    NA, ref, vg, table = args
    ef = make_ef(table, NA, ref)
    rows = []
    for S in SS:
        for tol in (0.01, 0.02):
            o = min_traps(vg, ef, np.exp(S), tol, freezing=True)
            o.update(NA=NA, ref=ref, vg=vg, S=S, tol=tol, freezing=True); rows.append(o)
    for tol in (0.01, 0.02):
        o = min_traps(vg, ef, 1.0, tol, freezing=False)
        o.update(NA=NA, ref=ref, vg=vg, S=0.0, tol=tol, freezing=False); rows.append(o)
    return rows


if __name__ == "__main__":
    npr = int(os.environ.get("NPROC", "2"))
    with Pool(npr) as pool:
        table = pool.map(ef_job, [(NA, T) for NA in NAS for T in TEMPS], chunksize=1)
    res = dict(ef_table=table, s_req=[], lp=[], release={vg: release_count(vg) for vg in VGS})
    print("== (A) Fermi level of the free electrons, E_F - reference (meV), 2 K -> 150 K")
    for NA in NAS:
        for ref in REFS:
            ef = make_ef(table, NA, ref)
            print(f"NA {NA:.0e} {ref:>3}: " + " | ".join(
                f"{vg} {1e3 * ef(NH[vg][0], 2.0):.0f} -> {1e3 * ef(NH[vg][-1], 150.0):.0f}" for vg in VGS), flush=True)
    print("== (B) entropy a fixed-level trap needs to fill on cooling, S_req / k_B per interval (2-10, 10-50, 50-100, 100-150 K)")
    for NA in NAS:
        for ref in REFS:
            ef = make_ef(table, NA, ref)
            for vg in VGS:
                s = s_required(vg, ef)
                res["s_req"].append(dict(NA=NA, ref=ref, vg=vg, intervals=s))
                print(f"NA {NA:.0e} {ref:>3} {vg}: " + ", ".join(f"{x['S_req']:.1f}" for x in s)
                      + " | dE_F (meV) " + ", ".join(f"{x['dEF_meV']:.1f}" for x in s), flush=True)
    print("== (C) minimum number of trap states (cm^-2); freezing classes included unless 'eq only'")
    jobs = [(NA, ref, vg, table) for NA in NAS for ref in REFS for vg in VGS]
    with Pool(npr) as pool:
        out = pool.map(lp_job, jobs, chunksize=1)
    for rows in out:
        res["lp"].extend(rows)
        for o in rows:
            tag = "" if o["freezing"] else " (eq only)"
            val = (f"total {o['total']:.2e}, N_extra {o['N_extra']:.1e}, frozen share {o['frozen_share']:.2f}"
                   if o["status"] == "feasible" else o["status"])
            print(f"NA {o['NA']:.0e} {o['ref']:>3} {o['vg']} S/k {o['S']:6.2f} tol {o['tol']:.2f}{tag}: {val}", flush=True)
    print("== (C) summary: smallest total over all S, references, NA (tol 0.02), and the S at which the total first "
          f"drops below the {SITES:.1e} cm^-2 interface sites")
    for vg in VGS:
        f = [o for o in res["lp"] if o["vg"] == vg and o["tol"] == 0.02 and o["status"] == "feasible"]
        fails = [o for o in res["lp"] if o["vg"] == vg and o["status"].startswith("solver")]
        best = min(f, key=lambda o: o["total"]) if f else None
        thr = {}
        for NA in NAS:
            for ref in REFS:
                ok = [o["S"] for o in f if o["NA"] == NA and o["ref"] == ref and o["freezing"] and o["total"] <= SITES]
                thr[f"{NA:.0e}/{ref}"] = max(ok) if ok else None
        print(f"{vg}: min total {best['total']:.2e} at S {best['S']}, {best['ref']}, NA {best['NA']:.0e} | S threshold for "
              f"<= sites: {thr} | solver failures {len(fails)}" if best else f"{vg}: nothing feasible", flush=True)
        # monotonicity in the tolerance
        for o1 in [o for o in res["lp"] if o["vg"] == vg and o["tol"] == 0.01 and o["status"] == "feasible"]:
            o2 = [o for o in res["lp"] if o["vg"] == vg and o["tol"] == 0.02 and o["NA"] == o1["NA"] and o["ref"] == o1["ref"]
                  and o["S"] == o1["S"] and o["freezing"] == o1["freezing"]][0]
            if o2["status"] != "feasible" or o2["total"] > 1.001 * o1["total"]:
                print(f"   NON-MONOTONE: {o1['NA']:.0e} {o1['ref']} S {o1['S']} freeze {o1['freezing']}: tol 0.01 "
                      f"{o1['total']:.2e}, tol 0.02 {o2.get('total', o2['status'])}", flush=True)
    print("== (D) release counting (first version; not a valid test)")
    for vg in VGS:
        for x in res["release"][vg]:
            print(f"{vg} {x['T1']:.0f}-{x['T2']:.0f} K: dn {x['dn']:.2e}, {x['per_eV']:.2e} per eV", flush=True)
    json.dump(res, open(os.path.join(here, "trapping.json"), "w"), indent=1)
