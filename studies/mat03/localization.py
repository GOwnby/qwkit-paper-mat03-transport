"""MAT-03 (10): the alternative to trapping. The electrons stay in the 2DEG (total density and Fermi level fixed), and a
fraction n_s that is immobile at 2 K becomes mobile on warming. Two channels, n_f + n_s = N_tot: at each temperature the
measured conductance sigma = e n_H mu_H and Hall density n_H fix the two mobilities,
  n_f mu_f + n_s mu_s = n_H mu_H,   n_f mu_f^2 + n_s mu_s^2 = n_H mu_H^2.
n_f is set by mu_s(2 K) = 0 (n_f = n_H(2 K)); N_tot = n_H(150 K) (channels equal at 150 K) or 5-10% above it.
This is a reparametrization of the data (two equations, two unknowns at each T), not a test. Reported: mu_f(T),
mu_s(T), the slow channel's sheet conductance in e^2/h (hopping cannot exceed about 1), power-law and Arrhenius fits of
mu_s over 10-100 K (three points: they do not distinguish the forms), and the model's predictions for the Hall
resistance at 7 T (deviation from the low-field line) and the classical magnetoresistance at 7 T at every temperature. Writes localization.json."""
import json, os
import numpy as np
from common3 import NH, MU, TEMPS, here
from qwkit import transport as tr
QE = 1.602176634e-19; G0 = QE ** 2 / 6.62607015e-34


def invert(nH, muH, nf, ns):
    """Both roots (mu_f, mu_s) with mu_s <= mu_f, or None."""
    A = nH * muH; B = nH * muH ** 2
    # n_f a^2 + n_s b^2 = B with a = (A - n_s b)/n_f  ->  quadratic in b
    qa = ns ** 2 / nf + ns; qb = -2 * A * ns / nf; qc = A ** 2 / nf - B
    disc = qb ** 2 - 4 * qa * qc
    if disc < -1e-9 * qb ** 2:
        return None
    disc = max(disc, 0.0)
    b = (-qb - np.sqrt(disc)) / (2 * qa)                  # the smaller root: slow channel
    a = (A - ns * b) / nf
    return float(a), float(b)


out = {}
for vg in ("4V", "5V", "7V", "8V"):
    nH = np.array(NH[vg], float); mu = np.array(MU[vg], float)
    rows = []
    for extra in (1.0, 1.05, 1.10):
        Nt = nH[-1] * extra
        nf = nH[0]; ns = Nt - nf
        # at 2 K mu_s = 0 exactly requires n_H(2K) = n_f; solve every temperature
        sol = [invert(nH[i], mu[i], nf, ns) if (np.isfinite(nH[i]) and np.isfinite(mu[i])) else None for i in range(len(TEMPS))]
        mf = [None if s is None else s[0] for s in sol]; ms = [None if s is None else s[1] for s in sol]
        fit = None
        T = np.array(TEMPS, float)
        ok = [i for i in range(1, 4) if ms[i] is not None and ms[i] > 0]
        if len(ok) >= 2:
            p = np.polyfit(np.log(T[ok]), np.log([ms[i] for i in ok]), 1)
            q = np.polyfit(1 / T[ok], np.log([ms[i] for i in ok]), 1)
            res_pl = float(np.std(np.log([ms[i] for i in ok]) - np.polyval(p, np.log(T[ok]))))
            res_ar = float(np.std(np.log([ms[i] for i in ok]) - np.polyval(q, 1 / T[ok])))
            fit = dict(exponent=float(p[0]), rms_log_powerlaw=res_pl, Ea_meV=float(-q[0] * 8.617e-2), rms_log_arrhenius=res_ar)
        ratio = [None if (a is None or b is None) else b / a for a, b in zip(mf, ms)]
        gs = [None if b is None else float(QE * ns * max(b, 0.0) / G0) for b in ms]          # sigma_s in e^2/h (n in cm^-2, mu in cm^2/Vs)
        bend, mr7 = [], []
        for a, b in zip(mf, ms):
            if a is None or b is None:
                bend.append(None); mr7.append(None); continue
            n2 = [nf, ns]; m2 = [a, max(b, 1e-6)]
            bend.append(float(tr.hall(n2, m2, 0.0)[0] / tr.hall(n2, m2, 7.0)[0] - 1))
            mr7.append(float(tr.magnetoresistance(n2, m2, 7.0)))
        rows.append(dict(N_tot=Nt, n_f=nf, n_s=ns, mu_f=mf, mu_s=ms, ratio=ratio, sigma_s_e2h=gs, Rxy_bend_7T=bend,
                         MR_7T=mr7, fit=fit))
        print(f"{vg} N_tot = {extra:.2f} n_H(150K) = {Nt:.2e}, n_s = {ns:.2e} ({100 * ns / Nt:.0f}%): mu_f "
              f"{[None if x is None else round(x) for x in mf]}, mu_s {[None if x is None else round(x) for x in ms]}"
              + f", mu_s/mu_f {[None if x is None else round(x, 2) for x in ratio]}"
              + f", sigma_s (e^2/h) {[None if x is None else round(x, 1) for x in gs]}"
              + f", R_xy(7T) bend {[None if x is None else round(100 * x, 1) for x in bend]}%"
              + f", MR(7T) {[None if x is None else round(100 * x, 1) for x in mr7]}%"
              + ("" if fit is None else f", 10-100 K: mu_s ~ T^{fit['exponent']:.2f} (rms {fit['rms_log_powerlaw']:.2f}) vs Arrhenius "
                 f"{fit['Ea_meV']:.1f} meV (rms {fit['rms_log_arrhenius']:.2f})"), flush=True)
    out[vg] = rows
json.dump(out, open(os.path.join(here, "localization.json"), "w"), indent=1)
