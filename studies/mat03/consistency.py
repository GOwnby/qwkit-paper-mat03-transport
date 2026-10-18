"""MAT-03 (5): model-light consistency checks on the published data.
(a) Sheet-resistance upturn below 50 K, from Fig. 3a (digitized R_s, lee_fig3a_Rs.json) and from the Hall data of
    Fig. 3b (R_s = 1 / (e n_H mu_H)); its size in units of the quantum-correction scale: the coefficient
    A = [sigma(50 K) - sigma(2 K)] / [G00 ln(50/2)], G00 = e^2 / (2 pi^2 hbar); weak localization and the
    electron-electron interaction correction give A of order 1 (A <= about 2 together). Also k_F l = (h/e^2) / R_s.
(b) Hall nonlinearity required by the two-subband reading: for each gate voltage and temperature, with the model
    ground-subband fraction f1 at n_tot = n_H(150 K) (subbands.json) and the measured conductance, the mobility ratio x
    for which the Hall density *from a straight-line fit of R_xy(B) over 0-7 T* equals the measured value, and the
    resulting maximum deviation of R_xy(B) from that line (in % of R_xy(7 T)), to set against the statement of Lee et
    al. that R_xy was linear up to 7 T at all gate voltages and temperatures. Writes consistency.json."""
import json
from scipy.optimize import brentq
from common3 import *
from qwkit import transport as tr
from hall_algebra import pops

G00 = 1.602176634e-19 ** 2 / (2 * np.pi ** 2 * 1.054571817e-34)       # S
RK = 25812.807
A3 = json.load(open(os.path.join(here, "lee_fig3a_Rs.json")))
BS = np.linspace(0.25, 7.0, 28)


def fit_n(n, mu):
    R = np.array([B / tr.hall(n, mu, B)[0] for B in BS])
    s = np.sum(R * BS) / np.sum(BS * BS)
    return 1 / s, float(np.max(np.abs(R - s * BS)) / R[-1])


out = dict(upturn=[], linearity=[])
print("== (a) resistance upturn below 50 K")
for vg in ("4V", "5V", "6V", "7V", "8V"):
    T3 = A3["T"]; R3 = A3["Rs"].get(vg)
    r3 = None if R3 is None else {t: R for t, R in zip(T3, R3)}
    rh = {t: (None if n is None or m is None else 1 / (1.602176634e-19 * n * m)) for t, n, m in zip(TEMPS, NH[vg], MU[vg])}
    if rh[2] is None:
        rh = None
    row = dict(vg=vg)
    for lab, R in (("fig3a", r3), ("hall", rh)):
        if R is None:
            row[lab] = None
            continue
        ds = 1 / R[50] - 1 / R[2]
        row[lab] = dict(R2=R[2], R50=R[50], ratio=R[2] / R[50], dsigma=ds, A=ds / (G00 * np.log(25.0)), kFl=RK / R[2])
    if rh is not None:
        row["dnH"] = NH[vg][0] / NH[vg][2] - 1; row["dmuH"] = MU[vg][0] / MU[vg][2] - 1
    if r3 is not None and rh is not None:
        row["scale_3a_over_hall"] = [r3[t] / rh[t] for t in (2, 10, 50, 100, 150)]
    out["upturn"].append(row)
    f = lambda lab, q, fmt: "-" if row[lab] is None else format(row[lab][q], fmt)
    print(f"{vg}: R(2K)/R(50K) = {f('fig3a', 'ratio', '.3f')} (Fig. 3a), {f('hall', 'ratio', '.3f')} (Hall) | A = "
          f"{f('fig3a', 'A', '.0f')} (Fig. 3a), {f('hall', 'A', '.0f')} (Hall) | k_F l = {f('fig3a', 'kFl', '.0f')} / "
          f"{f('hall', 'kFl', '.0f')}" + (f" | 50->2 K: n_H {100 * row['dnH']:+.0f}%, mu_H {100 * row['dmuH']:+.0f}%" if "dnH" in row else "")
          + (f" | R_s(Fig. 3a)/R_s(Hall) at 2-150 K: {np.round(row['scale_3a_over_hall'], 2).tolist()}" if "scale_3a_over_hall" in row else ""),
          flush=True)

print("== (b) two-subband reading with straight-line Hall fits over 0-7 T")
for NA in (1e18, 5e19):
    for vg in ("3V", "4V", "5V", "7V", "8V"):
        N = NH[vg][-1]
        for iT, T in enumerate(TEMPS[:-1]):
            nh, muh = NH[vg][iT], MU[vg][iT]
            if nh is None:
                continue
            f1 = pops(NA, T, N); S1 = nh * muh; n = [N * f1, N * (1 - f1)]
            row = dict(NA=NA, vg=vg, T=T, f1=f1, nH=nh, x=None)
            if nh >= N:
                row.update(x=1.0, dev=0.0)
            elif 1 - f1 > 1e-3:
                def h(lx):
                    x = np.exp(lx); m1 = S1 / (N * (f1 + (1 - f1) * x))
                    return fit_n(n, [m1, m1 * x])[0] - nh
                xs = np.linspace(0.0, np.log(300.0), 120); hv = [h(v) for v in xs]
                for a, b, ha, hb in zip(xs[:-1], xs[1:], hv[:-1], hv[1:]):
                    if ha * hb < 0:
                        lx = brentq(h, a, b, xtol=1e-4); x = float(np.exp(lx)); m1 = S1 / (N * (f1 + (1 - f1) * x))
                        nf, dev = fit_n(n, [m1, m1 * x]); n0 = tr.hall(n, [m1, m1 * x])[0]
                        row.update(x=x, mu1=m1, mu2=m1 * x, nH0=n0, dev=dev)
                        break
            out["linearity"].append(row)
            if row.get("x") == 1.0 and nh >= N:
                print(f"NA {NA:.0e} {vg} {T:3d} K: n_H = n_tot", flush=True)
            elif row["x"] is None:
                print(f"NA {NA:.0e} {vg} {T:3d} K: f1 {f1:.3f} | no two-subband solution", flush=True)
            else:
                print(f"NA {NA:.0e} {vg} {T:3d} K: f1 {f1:.3f} | x {row['x']:.2f} mu1 {row['mu1']:.0f} mu2 {row['mu2']:.0f} | "
                      f"n_H(B->0) {row['nH0']:.2e} vs fit {nh:.2e} | R_xy off the line by {100 * row['dev']:.1f}%", flush=True)
json.dump(out, open(os.path.join(here, "consistency.json"), "w"), indent=1)
