"""MAT-03 (6): what the published magnetoresistance at 2 K (Lee et al. Fig. 5a, 5 and 7 V; lee_fig5a_MR.json) allows.
Coupled subbands in a magnetic field are exactly equivalent to independent Drude channels (eigenmodes of the relaxation
matrix, non-negative weights summing to n_tot; tests/test_transport.py), so the question is which sets of channels
reproduce, at 2 K with n_tot = n_H(150 K): the measured conductance, the measured Hall density as Lee et al. obtained
it (slope of a straight line through R_xy(B) up to 7 T; here a least-squares line through the origin over 0.25-7 T),
and a magnetoresistance that changes between 2 and 7 T by no more than the measured change plus a tolerance (0.5 or
1.0 point, for digitization and a weak-localization or interaction tail).
(a) Weak localization: the negative magnetoresistance saturates within ~1 T; its size bounds the whole
    field-suppressible correction to the 2 K conductance.
(b) Two and three channels with random weights (Dirichlet) and log-uniform mobilities 1-30000 cm^2/Vs, rescaled to the
    measured conductance and accepted if the fitted Hall density is within 1% of the measured one; among the sets with
    a flat magnetoresistance, the smallest total weight of channels slower than 25, 50, 100 and 200 cm^2/Vs, and the
    range of the slowest channel's mobility. Writes mr.json."""
import json
from common3 import *
from qwkit import transport as tr

M = json.load(open(os.path.join(here, "lee_fig5a_MR.json")))
B = M["B_T"]; TOL = (0.5, 1.0); CUTS = (25.0, 50.0, 100.0, 200.0)
BS = np.linspace(0.25, 7.0, 28)
rng = np.random.default_rng(1)


def rho(n, mu, Bv):
    """rho_xx, rho_xy for channel arrays n (k, m), mu (k, m) (cm^-2, cm^2/Vs) at fields Bv (T): arrays (k, len(Bv))."""
    m = mu[..., None] * 1e-4; x = m * Bv[None, None, :]; nn = n[..., None] * 1e4
    sxx = np.sum(nn * m / (1 + x ** 2), axis=1) * QE; sxy = np.sum(nn * m * x / (1 + x ** 2), axis=1) * QE
    d = sxx ** 2 + sxy ** 2
    return sxx / d, sxy / d


QE = 1.602176634e-19
out = dict(wl={}, sets={})
for vg in ("5V", "7V"):
    mr = dict(zip(B, M["MR_percent"][vg]))
    N = NH[vg][-1]; nh = NH[vg][0]; S = NH[vg][0] * MU[vg][0]
    s50 = NH[vg][2] * MU[vg][2]
    meas = mr[7] - mr[2]
    out["wl"][vg] = dict(plateau=mr[2], change_2_7=meas, upturn_hall=100 * (s50 / S - 1), sigma_drop_hall=100 * (1 - S / s50))
    print(f"{vg}: MR plateau (2 T) {mr[2]:.2f}%, change 2->7 T {meas:+.2f}% | R_s(2K)/R_s(50K)-1 (Hall) {100 * (s50 / S - 1):.1f}%,"
          f" sigma(2K)/sigma(50K)-1 {-100 * (1 - S / s50):.1f}% | n_H/n_tot {nh / N:.3f}", flush=True)
    res = {}
    for k in (2, 3):
        n_try = 400000; acc = []
        for chunk in range(n_try // 20000):
            w = rng.dirichlet(np.ones(k), 20000)
            mu = np.exp(rng.uniform(np.log(1.0), np.log(30000.0), (20000, k)))
            mu *= (S / N) / np.sum(w * mu, axis=1)[:, None]
            n = w * N
            rxx, rxy = rho(n, mu, BS)
            sl = np.sum(rxy * BS[None, :], axis=1) / np.sum(BS * BS)
            nfit = 1 / (QE * sl) * 1e-4
            keep = np.abs(nfit / nh - 1) < 0.01
            if not np.any(keep):
                continue
            r0, _ = rho(n[keep], mu[keep], np.array([1e-6]))
            r27, _ = rho(n[keep], mu[keep], np.array([2.0, 7.0]))
            m27 = 100 * (r27[:, 1] - r27[:, 0]) / r0[:, 0]
            for wi, mi, mm in zip(w[keep], mu[keep], m27):
                acc.append((wi, mi, mm))
        rk = dict(accepted=len(acc))
        for tol in TOL:
            good = [(wi, mi) for wi, mi, mm in acc if mm <= meas + tol]
            rk[str(tol)] = dict(n=len(good),
                                min_slow={str(c): (float(min(np.sum(wi[mi < c]) for wi, mi in good)) if good else None) for c in CUTS},
                                slowest_mu=(float(min(mi.min() for _, mi in good)), float(max(mi.min() for _, mi in good))) if good else None,
                                slowest_w=(float(min(wi[np.argmin(mi)] for wi, mi in good)), float(max(wi[np.argmin(mi)] for wi, mi in good))) if good else None)
        res[str(k)] = rk
        for tol in TOL:
            g = rk[str(tol)]
            print(f"   {k} channels: {rk['accepted']} sets reproduce sigma and the fitted n_H; {g['n']} with flat MR (tol {tol}); "
                  f"smallest weight below 25/50/100/200 cm^2/Vs: " + "/".join("-" if v is None else f"{v:.3f}" for v in g["min_slow"].values())
                  + (f"; slowest channel mu {g['slowest_mu'][0]:.0f}-{g['slowest_mu'][1]:.0f}, weight {g['slowest_w'][0]:.3f}-{g['slowest_w'][1]:.3f}" if g["slowest_mu"] else ""),
                  flush=True)
    out["sets"][vg] = res
json.dump(out, open(os.path.join(here, "mr.json"), "w"), indent=1)
