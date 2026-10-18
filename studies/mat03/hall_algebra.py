"""MAT-03 (2): what the temperature dependence of the Hall density requires, given the model's subband populations.
For each gate voltage the ions are frozen, so the total electron density n_tot is fixed; the Hall density of parallel
Drude channels is n_H = (sum n_i mu_i)^2 / sum n_i mu_i^2 <= n_tot, with equality only for equal mobilities. We take
the minimal assumption n_tot = n_H(150 K) (equal mobilities at 150 K) and, with the subband populations of the model
at n_tot (subbands.json, interpolated in n_s), ask:
  (a) the immobile-upper-subband limit: if only subband 1 conducted at low T, n_H(2 K) = n_1, so n_1 / n_tot must
      equal n_H(2 K) / n_H(150 K); report both;
  (b) a single mobility ratio x = mu_upper / mu_1 (all upper subbands alike) reproducing n_H(T) at each T: both
      branches, x < 1 (upper subbands slower) and x > 1 (upper subbands faster), where they exist;
  (c) the subband mobilities mu_1, mu_upper implied by the measured Hall mobility on each branch.
Also the conductance sigma = e n_H mu_H at each T. Writes hall_algebra.json."""
import json
from scipy.optimize import brentq
from common3 import *

S = json.load(open(os.path.join(here, "subbands.json")))


def pops(NA, T, ntot):
    rows = [s for s in S if s["NA"] == NA and s["T"] == T][0]["rows"]
    ns = np.array([r["ns"] for r in rows])
    f1 = np.array([r["n"][0] / r["ns"] for r in rows])
    return float(np.interp(ntot, ns, f1))


def g(x, f1):
    f2 = 1 - f1
    return (f1 + f2 * x) ** 2 / (f1 + f2 * x ** 2)


def solve_x(target, f1):
    """Mobility ratios x giving n_H / n_tot = target: (x_low in (0,1), x_high > 1); None where absent."""
    lo = hi = None
    if 1 - f1 < 1e-3:                       # a single occupied subband: n_H = n_tot whatever the mobility
        return lo, hi
    if target < 1:
        if g(1e-9, f1) < target < 1:
            lo = brentq(lambda x: g(x, f1) - target, 1e-9, 1.0)
        if g(1e3, f1) < target < 1:          # ratios above 1000 are not considered
            hi = brentq(lambda x: g(x, f1) - target, 1.0, 1e3)
    return lo, hi


if __name__ == "__main__":
    out = []
    for NA in (5e19, 1e18):
        for vg in VG[2:]:
            nt = NH[vg][-1]
            for iT, T in enumerate(TEMPS):
                nh = NH[vg][iT]; muh = MU[vg][iT]
                if nh is None:
                    continue
                f1 = pops(NA, float(T), nt)
                target = nh / nt
                lo, hi = solve_x(target, f1)
                row = dict(NA=NA, vg=vg, T=T, n_tot=nt, nH=nh, muH=muh, ratio=target, f1=f1, x_low=lo, x_high=hi,
                           sigma_e=nh * muh)
                for tag, x in (("low", lo), ("high", hi)):
                    if x is not None:
                        mu1 = muh * (f1 + (1 - f1) * x) / (f1 + (1 - f1) * x ** 2)
                        row[f"mu1_{tag}"] = mu1; row[f"mu2_{tag}"] = mu1 * x
                out.append(row)
                print(f"NA {NA:.0e} {vg} T {T:4.0f}: n_H/n_tot {target:.3f} f1 {f1:.3f} | x_low {lo if lo is None else round(lo, 3)} "
                      f"x_high {hi if hi is None else round(hi, 2)}" + (f" | mu1/mu2 (high) {row['mu1_high']:.0f}/{row['mu2_high']:.0f}"
                                                                        if hi is not None else ""), flush=True)
    json.dump(out, open(os.path.join(here, "hall_algebra.json"), "w"), indent=1)
