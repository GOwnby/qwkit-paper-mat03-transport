"""MAT-03 (1): subband populations of the LaInO3/BaSnO3 2DEG versus total density and temperature. Field-free cap (the
positive charge at or next to the interface, as MAT-02 requires; by Gauss's law the BSO side does not depend on where
in the cap it sits while the cap is empty). For n_s = 0.5e13-1.4e14 cm^-2 and T = 2, 10, 50, 100, 150 K, compensated and
clean BSO: subband populations n_i, edges E_i, Fermi level, in-plane masses m_i, the Fermi wavevector of each occupied
subband and the cap fraction. Writes subbands.json."""
import json
from multiprocessing import Pool
from common3 import *

NS = [0.5e13, 1e13, 1.5e13, 2e13, 2.5e13, 3e13, 3.5e13, 4e13, 4.5e13, 5e13, 6e13, 7e13, 8e13, 9e13, 1e14, 1.08e14, 1.2e14, 1.4e14]


def job(args):
    NA, T = args
    p = dict(BASE); p["NA"] = NA; p["T"] = T; p["f_int"] = 1.0
    rows = []
    for ns in NS:
        r = solve(ns, p, extend_to=0.0, nstates=10, dz=0.5)
        z = r["z"]; dzz = z[1] - z[0]
        capw = [float(np.sum(r["psi"][z < 0, j] ** 2) * dzz) for j in range(len(r["E"]))]
        n = r["n_n"]
        occ = [i for i in range(len(n)) if n[i] > 1e-4 * ns]
        kF = [float(np.sqrt(2 * np.pi * max(n[i], 0) * 1e4) * 1e-10) for i in range(len(n))]   # 1/A, spin 2: n = k^2 / 2 pi
        rows.append(dict(ns=ns, E=[float(e) for e in r["E"]], mu=float(r["mu"]), n=[float(x) for x in n], capw=capw,
                         m=[float(x) for x in r["m_n"]], kF=kF, occupied=occ,
                         f1=float(n[0] / ns), cap_frac=float(sum(n[i] for i in range(len(n)) if capw[i] > 0.5) / ns)))
        print(f"NA {NA:.0e} T {T:5.0f} ns {ns:.2e}: n_i/ns {np.round(n[:5] / ns, 3)} | E_F - E_i (meV) "
              f"{np.round((r['mu'] - r['E'][:4]) * 1e3, 0)} | cap {rows[-1]['cap_frac']:.3f}", flush=True)
    return dict(NA=NA, T=T, rows=rows)


if __name__ == "__main__":
    jobs = [(NA, T) for NA in (5e19, 1e18) for T in (2.0, 10.0, 50.0, 100.0, 150.0)]
    with Pool(int(os.environ.get("NPROC", "2"))) as pool:
        res = pool.map(job, jobs, chunksize=1)
    json.dump(res, open(os.path.join(here, "subbands.json"), "w"), indent=1)
