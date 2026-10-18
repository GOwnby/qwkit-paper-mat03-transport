"""MAT-03 Fig. 3 from trapping.json and localization.json: (a) the Fermi level of the free electrons in the trapping
reading, E_F - E_1 against temperature (clean and compensated BSO), with the k_B T scale; (b) the two parallel channels
that reproduce the measured conductance and Hall density at fixed total charge. Writes ../../paper_mat03/figures/fig3_trapping.pdf."""
import numpy as np
from paper_figures import plt, J, VC, clean, label, outdir, INK2, os
import trapping as T
from common3 import NH, TEMPS

TR = J("trapping.json"); LO = J("localization.json")
fig, axs = plt.subplots(1, 2, figsize=(5.4, 2.45), gridspec_kw=dict(wspace=0.45))
ax = axs[0]
Tq = np.array(TEMPS, float)
for NA, ls in ((1e18, "-"), (5e19, "--")):
    ef = T.make_ef(TR["ef_table"], NA, "E1")
    for vg in ("5V", "7V", "8V"):
        e = [1e3 * ef(NH[vg][i], Tq[i]) for i in range(len(Tq))]
        ax.plot(Tq, np.array(e) - e[-1], ls=ls, color=VC[vg], marker="o", ms=2.4,
                label=f"{vg}" if NA == 1e18 else None)
ax.fill_between(Tq, -3 * 8.617e-2 * Tq, 3 * 8.617e-2 * Tq, color="#bbbbbb", alpha=0.45, lw=0, label=r"$\pm3k_BT$")
ax.set_xscale("log"); ax.set_xlabel("$T$ (K)"); ax.set_ylabel(r"$E_F-E_1$ rel. to 150 K (meV)")
ax.legend(frameon=False, fontsize=6.2, loc="lower right"); clean(ax); label(ax, "a")
ax = axs[1]
for vg in ("5V", "7V", "8V"):
    r = LO[vg][0]
    ax.plot(Tq, r["mu_f"], color=VC[vg], marker="o", ms=2.4, label=f"{vg}, mobile ({100 * r['n_f'] / r['N_tot']:.0f}%)")
    ax.plot(Tq, [max(x, 1e-9) for x in r["mu_s"]], color=VC[vg], marker="s", ms=2.4, mfc="white", ls="--",
            label=f"{vg}, slow ({100 * r['n_s'] / r['N_tot']:.0f}%)")
ax.set_xscale("log"); ax.set_xlabel("$T$ (K)"); ax.set_ylabel(r"channel mobility (cm$^2$/Vs)"); ax.set_ylim(0, 3600)
ax.legend(frameon=False, fontsize=5.2, loc="upper left", ncol=2, columnspacing=0.8, handlelength=1.8); clean(ax); label(ax, "b")
fig.savefig(os.path.join(outdir, "fig3_trapping.pdf"))
print("fig3 written")
