"""MAT-03 (7): interface roughness added to the Coulomb scattering by the gate-induced charge (clean BSO, Z = 1 at 6 A,
eps_cap = 24):
abrupt offset 1.6 eV at the interface, Gaussian height correlations with rms Delta = 2 or 4 A (one BaSnO3 unit cell is
4.1 A; AFM of the BaSnO3 surface in Lee et al.: 0.38 nm rms) and correlation length Lambda = 15-200 A. Reports subband
mobilities, the ratio x of the upper subbands to the ground subband, and mu_sigma against the conductance bound.
Also the mean distance <z> of each occupied subband from the interface. Writes roughness.json."""
import json
from multiprocessing import Pool
from common3 import *
from qwkit import transport as tr

LAM = (15.0, 30.0, 60.0, 100.0, 200.0)
DEL = (2.0, 4.0)


def job(vg):
    p = dict(BASE); p["NA"] = 1e18; p["T"] = 2.0; p["f_int"] = 1.0
    r0 = solve(5e12, p, extend_to=0.0, nstates=10, dz=0.5); sigP = 5e12 + r0["N_dep"]
    nt = NH[vg][-1]
    r = solve(nt, p, extend_to=0.0, nstates=10, dz=0.5); dn = nt + r["N_dep"] - sigP
    sb = tr.Subbands(r); z = sb.z * 1e10; dz = z[1] - z[0]
    zm = [float(np.sum(z * sb.psi[:, i] ** 2) * dz * 1e-10) for i in range(sb.ns)]
    base = tr.mobilities(r, eps=20.0, sheets=[(-6.0, dn, 1)], volume=[(0.0, r["w_d"], 1e18, -1)], nth=48, eps_cap=24.0)
    rows = [dict(Delta=0.0, Lambda=0.0, mu=base["mu"].tolist())]
    for D in DEL:
        for Lm in LAM:
            o = tr.mobilities(r, eps=20.0, sheets=[(-6.0, dn, 1)], volume=[(0.0, r["w_d"], 1e18, -1)],
                              roughness=(1.6, D, Lm, 0.0), nth=48, eps_cap=24.0)
            rows.append(dict(Delta=D, Lambda=Lm, mu=o["mu"].tolist()))
    n = base["n"]
    for row in rows:
        mu = np.array(row["mu"])
        row["x"] = float(np.sum(n[1:] * mu[1:]) / np.sum(n[1:]) / mu[0])
        row["mu_sigma"] = float(np.sum(n * mu) / np.sum(n))
        print(f"{vg} Delta {row['Delta']} Lambda {row['Lambda']:5.0f}: mu_i {np.round(mu).astype(int).tolist()} x {row['x']:.2f} "
              f"mu_sigma {row['mu_sigma']:.0f}", flush=True)
    print(f"{vg} <z> of subbands (A): {np.round(zm, 1).tolist()}", flush=True)
    return dict(vg=vg, n=n.tolist(), zmean=zm, rows=rows)


if __name__ == "__main__":
    with Pool(int(os.environ.get("NPROC", "2"))) as pool:
        res = pool.map(job, ["4V", "5V", "7V", "8V"], chunksize=1)
    json.dump(res, open(os.path.join(here, "roughness.json"), "w"), indent=1)
