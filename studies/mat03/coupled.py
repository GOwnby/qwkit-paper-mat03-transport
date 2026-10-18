"""MAT-03 (5): Hall density and magnetoresistance at 2 K produced by the Coulomb scattering of the gate-induced charge
(clean BSO, 5, 7, 8 V, eps_cap = 24, Z = 1, 2, d = 1-24 A), with the coupled-subband magnetotransport
(transport.magnetotransport; exact for isotropic subbands, equivalent to independent Drude eigenmodes). The model
mobilities are larger than the measured ones, so the remaining (unidentified) scattering is added in two limits that
bracket its effect on the mobility differences: (mult) all rates scaled by one factor, which keeps the mobility ratios
and the low-field Hall density of the model; (add) an extra isotropic intrasubband rate, the same for all subbands,
which pushes the ratios towards 1. In both, the conductance is matched to the measured 2 K value (rows where the model
conductance is already below it are flagged: there d < d_min). Reported: low-field and fitted (0.25-7 T) Hall density
over n_tot, and the magnetoresistance rise between 2 and 7 T. Writes coupled.json."""
import json
from multiprocessing import Pool
from common3 import *
from qwkit import transport as tr
from coulomb import setup

BS = np.linspace(0.25, 7.0, 28)


def stats(o, nt):
    r0 = tr.magnetotransport(o, 0.0)[0]
    nH0 = tr.magnetotransport(o, 1e-3)[2]
    R = np.array([tr.magnetotransport(o, B)[1] * B for B in BS])
    nfit = 1 / (np.sum(R * BS) / np.sum(BS * BS)) / tr.QE * 1e-4
    mr = tr.magnetotransport(o, 7.0)[0] / r0 - tr.magnetotransport(o, 2.0)[0] / r0
    return dict(nH0=nH0 / nt, nHfit=nfit / nt, mr27=mr, mu_sigma=float(1 / (tr.QE * r0 * nt * 1e4) * 1e4))


def with_rates(o, a0, a1):
    q = dict(o); q["a0"] = a0; q["a1"] = a1
    q["tau"] = tr.relax(a0, a1, q["k"], np.asarray(q["md"]) * tr.M0)
    q["mu"] = (tr.QE * np.asarray(q["tau"]) / (np.asarray(q["md"]) * tr.M0) * 1e4).tolist()
    return q


def job(args):
    vg, Z, d = args
    r, nt, dn = setup(1e18, vg)
    o = tr.mobilities(r, eps=20.0, sheets=[(-d, dn / Z, Z)], volume=[(0.0, r["w_d"], 1e18, -1)], nth=64, eps_cap=24.0)
    o = {k: (np.asarray(v) if k in ("a0", "a1", "k", "md", "n", "mu", "tau") else v) for k, v in o.items()}
    target = NH[vg][0] * MU[vg][0] / nt                                  # sigma(2 K) / (e n_tot)
    base = stats(o, nt)
    res = dict(vg=vg, Z=Z, d=d, mu=o["mu"].tolist(), model=base, matched=bool(base["mu_sigma"] >= target))
    c = base["mu_sigma"] / target                                         # rates times c
    qm = with_rates(o, o["a0"] * c, o["a1"] * c)
    res["mult"] = dict(stats(qm, nt), mu=qm["mu"])
    lo, hi = 0.0, 1e15
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if stats(with_rates(o, o["a0"] + np.diag(np.full(len(o["k"]), mid)), o["a1"]), nt)["mu_sigma"] > target:
            lo = mid
        else:
            hi = mid
    qa = with_rates(o, o["a0"] + np.diag(np.full(len(o["k"]), 0.5 * (lo + hi))), o["a1"])
    res["add"] = dict(stats(qa, nt), mu=qa["mu"], s=0.5 * (lo + hi))
    m, a = res["mult"], res["add"]
    print(f"{vg} Z {Z} d {d:4.1f} {'matched' if res['matched'] else 'UNMATCHED'}: model mu {np.round(o['mu']).astype(int).tolist()} | "
          f"mult: n_H/n {m['nH0']:.3f} (fit {m['nHfit']:.3f}) MR {100 * m['mr27']:.2f}% | add: mu {np.round(a['mu']).astype(int).tolist()} "
          f"n_H/n {a['nH0']:.3f} (fit {a['nHfit']:.3f}) MR {100 * a['mr27']:.2f}%", flush=True)
    return res


if __name__ == "__main__":
    jobs = [(vg, Z, d) for vg in ("5V", "7V", "8V") for Z in (1, 2) for d in (1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0, 24.0)]
    with Pool(int(os.environ.get("NPROC", "2"))) as pool:
        out = pool.map(job, jobs, chunksize=1)
    json.dump(out, open(os.path.join(here, "coupled.json"), "w"), indent=1)
