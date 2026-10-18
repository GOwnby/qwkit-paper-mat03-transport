"""Figures for the MAT-03 manuscript from lee_fig3b_v2.json, lee_fig3a_Rs.json, lee_fig5a_MR.json, consistency.json,
mr.json, coulomb.json and numbers.json; Fig. 3 from trapping.json and localization.json (fig3_trapping.py). Writes ../../paper_mat03/figures/*.pdf."""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", ".."))
from qwkit import transport as tr

outdir = os.path.join(here, "..", "..", "paper_mat03", "figures"); os.makedirs(outdir, exist_ok=True)
BLUE, ORANGE, AQUA, PURPLE, RED = "#2a78d6", "#eb6834", "#1baf7a", "#8a4fbf", "#c4314b"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.family": "serif", "font.serif": ["STIX Two Text", "STIXGeneral", "DejaVu Serif"],
                     "mathtext.fontset": "stix", "font.size": 8.5, "axes.labelsize": 8.5, "xtick.labelsize": 7.5,
                     "ytick.labelsize": 7.5, "legend.fontsize": 6.6, "axes.linewidth": 0.6, "axes.edgecolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "lines.linewidth": 1.3, "savefig.bbox": "tight",
                     "savefig.pad_inches": 0.02, "pdf.fonttype": 42})
J = lambda f: json.load(open(os.path.join(here, f)))
L = J("lee_fig3b_v2.json"); A3 = J("lee_fig3a_Rs.json"); M5 = J("lee_fig5a_MR.json")
VC = {"2V": "#c9b458", "3V": "#e39a4b", "4V": RED, "5V": "#2b3a8f", "6V": "#3a6ea5", "7V": "#5aa0e0", "8V": "#2d9a4a"}
G00 = 1.602176634e-19 ** 2 / (2 * np.pi ** 2 * 1.054571817e-34)
QE = 1.602176634e-19


def clean(ax, grid=True):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    if grid:
        ax.grid(color=GRID, lw=0.5); ax.set_axisbelow(True)


def label(ax, s, x=-0.2):
    ax.text(x, 1.03, s, transform=ax.transAxes, fontsize=9.5, fontweight="bold", va="bottom")


def fig1():
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.45), gridspec_kw=dict(wspace=0.42))
    ax = axs[0]; T = np.array(L["T"])
    for vg, c in VC.items():
        n = np.array([np.nan if v is None else v for v in L["n_1e13"][vg]], float)
        ax.plot(T, n / n[-1], color=c, marker="o", ms=2.4, label=vg)
    ax.set_xscale("log"); ax.set_xlabel("$T$ (K)"); ax.set_ylabel(r"$n_H(T)\,/\,n_H(150\ \mathrm{K})$")
    ax.set_ylim(0.55, 1.05); ax.legend(frameon=False, ncol=2, loc="lower right", fontsize=5.8); clean(ax); label(ax, "a")
    ax = axs[1]; T3 = np.array(A3["T"])
    for vg in ("5V", "6V", "7V", "8V"):
        R = np.array(A3["Rs"][vg]); i50 = list(T3).index(50)
        ax.plot(T3[T3 <= 150], (R / R[i50])[T3 <= 150], color=VC[vg], marker="o", ms=2.4, label=f"{vg}")
    Rh = np.array([1 / (QE * n * 1e13 * m) for n, m in zip(L["n_1e13"]["4V"], L["mu"]["4V"])])
    ax.plot(T, Rh / Rh[2], color=VC["4V"], marker="s", ms=2.8, mfc="white", label="4V (Hall)")
    R8 = np.array(A3["Rs"]["8V"]); s50 = 1 / R8[list(T3).index(50)]; Tq = np.geomspace(2, 50, 30)
    ax.plot(Tq, s50 / (s50 - 2 * G00 * np.log(50 / Tq)), color=INK2, lw=0.8, ls="--", label=r"8V, $2G_{00}\ln(50\,\mathrm{K}/T)$")
    ax.set_xscale("log"); ax.set_xlabel("$T$ (K)"); ax.set_ylabel(r"$R_s(T)\,/\,R_s(50\ \mathrm{K})$")
    ax.legend(frameon=False, fontsize=5.2, loc="upper center", ncol=2, columnspacing=0.8); ax.set_ylim(0.5, 1.6); clean(ax); label(ax, "b")
    ax = axs[2]; B = np.array(M5["B_T"])
    MR = J("mr.json")
    for vg in ("3V", "5V", "7V"):
        ax.plot(B, M5["MR_percent"][vg], color=VC[vg], marker="o", ms=2.4, label=f"{vg}, measured")
    HA = J("hall_algebra.json")
    for vg in ("5V", "7V"):
        h = [h for h in HA if h["NA"] == 1e18 and h["vg"] == vg and h["T"] == 2][0]
        f1 = h["f1"]; x = h["x_high"]; N = h["n_tot"]; S = h["nH"] * h["muH"]
        m1 = S / (N * (f1 + (1 - f1) * x)); n = [N * f1, N * (1 - f1)]
        Bq = np.linspace(0, 7, 60); base = M5["MR_percent"][vg][2]
        cl = np.array([100 * tr.magnetoresistance(n, [m1, m1 * x], b) for b in Bq])
        ax.plot(Bq[Bq >= 2], (base + cl - 100 * tr.magnetoresistance(n, [m1, m1 * x], 2.0))[Bq >= 2], color=VC[vg], ls="--", lw=1.0,
                label=f"{vg}, two subbands, $x={x:.1f}$")
    ax.set_xlabel("$B$ (T)"); ax.set_ylabel("magnetoresistance at 2 K (%)"); ax.set_xlim(0, 7.3)
    ax.legend(frameon=False, fontsize=5.4, loc="upper left"); clean(ax); label(ax, "c")
    fig.savefig(os.path.join(outdir, "fig1_data.pdf")); plt.close(fig)


def family(vg, nw=120):
    """Two channels (fraction w at mu_s < mean, 1 - w at mu_f) with the measured 2 K conductance and fitted Hall density:
    for each w the slow mobility from a root search; returns w, mu_s, MR rise 2 -> 7 T (points)."""
    from scipy.optimize import brentq
    N = L["n_1e13"][vg][-1] * 1e13; nh = L["n_1e13"][vg][0] * 1e13; a = nh * L["mu"][vg][0] / N
    BS = np.linspace(0.25, 7.0, 28)

    def nfit(w, ms):
        mf = (a - w * ms) / (1 - w); n = [w * N, (1 - w) * N]
        R = np.array([B / tr.hall(n, [ms, mf], B)[0] for B in BS])
        return 1 / (np.sum(R * BS) / np.sum(BS * BS)), n, [ms, mf]
    out = []
    for w in np.linspace(0.01, 0.6, nw):
        f = lambda ms: nfit(w, ms)[0] / nh - 1
        if f(1e-3) * f(0.999 * a) > 0:
            continue
        ms = brentq(f, 1e-3, 0.999 * a, xtol=1e-3)
        _, n, mu = nfit(w, ms)
        out.append((w, ms, 100 * (tr.magnetoresistance(n, mu, 7.0) - tr.magnetoresistance(n, mu, 2.0))))
    return np.array(out)


def fig2():
    N = J("numbers.json"); CO = J("coulomb.json"); MRJ = J("mr.json"); CP = J("coupled.json")
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.5), gridspec_kw=dict(wspace=0.42))
    ax = axs[0]; V = ["4V", "5V", "7V", "8V"]; xv = [4, 5, 7, 8]
    ax.plot(xv, [L["n_1e13"][v][0] / L["n_1e13"][v][-1] for v in V], color=INK, marker="o", ms=3.2, label="measured, 2 K")
    xs = [5, 7, 8]; lo = []; hi = []
    for v in ("5V", "7V", "8V"):
        rr = [c for c in CP if c["vg"] == v and c["matched"]]
        vals = [c[k]["nH0"] for c in rr for k in ("mult", "add")]
        lo.append(min(vals)); hi.append(max(vals))
    ax.fill_between(xs, lo, hi, color=AQUA, alpha=0.4, lw=0, label="Coulomb scattering,\nmeasured conductance")
    ax.set_xlabel(r"$V_{GS}$ (V)"); ax.set_ylabel(r"$n_H\,/\,n_\mathrm{tot}$ at 2 K"); ax.set_ylim(0.55, 1.02)
    ax.legend(frameon=False, loc="lower left", fontsize=5.8); clean(ax); label(ax, "a")
    ax = axs[1]
    for vg, col in (("5V", VC["5V"]), ("7V", VC["7V"])):
        f = family(vg)
        meas = MRJ["wl"][vg]["change_2_7"]
        ax.plot(f[:, 0], f[:, 2] - meas, color=col, label=f"{vg}")
    ax.axhspan(-0.2, 1.0, color=GRID, lw=0); ax.axhline(0.5, color=INK2, lw=0.6, ls=":")
    ax.set_yscale("symlog", linthresh=1.0); ax.set_xlim(0, 0.6)
    ax.set_xlabel("fraction of electrons in the slow channel"); ax.set_ylabel(r"MR rise 2$\to$7 T above measured (points)")
    ax.legend(frameon=False, loc="upper right"); clean(ax); label(ax, "b")
    ax = axs[2]; vg = "8V"
    for e, ls in ((24.0, "-"), (38.0, "--"), (50.0, ":")):
        for Z, col in ((1, BLUE), (2, ORANGE)):
            g = [g for g in CO if g["NA"] == 1e18 and g["vg"] == vg and g["eps_cap"] == e][0]
            rr = [r for r in g["rows"] if r["Z"] == Z]
            ax.plot([r["d"] for r in rr], [r["mu_sigma"] for r in rr], color=col, ls=ls, lw=1.0, label=(f"$Z={Z}$" if e == 24.0 else None))
    row = [r for r in N["dmin"] if r["NA"] == 1e18 and r["vg"] == vg and r["eps_cap"] == 24.0 and r["Z"] == 1][0]
    ax.axhline(row["bound_hall"], color=INK, lw=0.7, ls=":"); ax.text(12, row["bound_hall"] * 1.1, "bound, Hall", fontsize=6, color=INK2)
    ax.axhline(row["bound_3a"], color=INK2, lw=0.7, ls="-."); ax.text(12, row["bound_3a"] * 0.72, "bound, Fig. 3a", fontsize=6, color=INK2)
    for e, al in ((24.0, 0.28), (38.0, 0.18), (50.0, 0.10)):
        rw = [r for r in N["dmin"] if r["NA"] == 1e18 and r["vg"] == vg and r["eps_cap"] == e and r["Z"] == 1][0]
        ax.axvspan(0.9, rw["xmax150"], color=AQUA, alpha=al, lw=0)
    for e, ls in ((24.0, "-"), (38.0, "--"), (50.0, ":")):
        ax.plot([], [], color=INK2, ls=ls, label=rf"$\varepsilon_\mathrm{{cap}}={e:.0f}$")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.9, 35); ax.set_ylim(150, 6e4)
    ax.set_xlabel(r"distance $d$ of the charge (Å)"); ax.set_ylabel(r"$\mu_\sigma$ (cm$^2$/Vs)")
    ax.set_title("clean BSO, 8 V", fontsize=7.8); ax.legend(frameon=False, fontsize=5.4, loc="upper left"); clean(ax); label(ax, "c")
    fig.savefig(os.path.join(outdir, "fig2_bounds.pdf")); plt.close(fig)


if __name__ == "__main__":
    for f in os.listdir(outdir):
        if f.endswith(".pdf"):
            os.remove(os.path.join(outdir, f))
    fig1(); fig2()
    import runpy
    runpy.run_path(os.path.join(here, "fig3_trapping.py"), run_name="__main__")      # Fig. 3 (trapping, localization)
    print("figures written to", os.path.abspath(outdir))
