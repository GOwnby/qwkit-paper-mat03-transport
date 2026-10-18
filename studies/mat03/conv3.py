"""MAT-03: numerical convergence of the subband mobilities (8 V, clean BSO, compensating charge Z = 1 at 1 and 6 A,
eps_cap = 24): angular
grid (nth), position grid (dz) and number of solved states. Writes conv3.json."""
import json, time
from common3 import *
from qwkit import transport as tr

p = dict(BASE); p["NA"] = 1e18; p["T"] = 2.0; p["f_int"] = 1.0
nt = NH["8V"][-1]
out = []
for dz, nst in ((0.5, 10), (0.25, 10), (0.5, 14)):
    r0 = solve(5e12, p, extend_to=0.0, nstates=nst, dz=dz); sigP = 5e12 + r0["N_dep"]
    r = solve(nt, p, extend_to=0.0, nstates=nst, dz=dz); dn = nt + r["N_dep"] - sigP
    for nth in ((32, 64, 128) if (dz, nst) == (0.5, 10) else (64,)):
        t = time.time()
        for d in (1.0, 6.0):
            o = tr.mobilities(r, eps=20.0, sheets=[(-d, dn, 1)], volume=[(0.0, r["w_d"], 1e18, -1)], nth=nth, eps_cap=24.0)
            n, mu = o["n"], o["mu"]
            row = dict(d=d, dz=dz, nstates=nst, nth=nth, mu=mu.tolist(), n=n.tolist(),
                       mu_sigma=float(np.sum(n * mu) / np.sum(n)), nH=tr.hall(n, mu)[0], secs=time.time() - t)
            out.append(row)
            print(f"d {d} dz {dz} nstates {nst} nth {nth}: mu_i {np.round(mu, 1).tolist()} mu_sigma {row['mu_sigma']:.1f} "
                  f"n_H {row['nH']:.4e} ({row['secs']:.0f} s)", flush=True)
json.dump(out, open(os.path.join(here, "conv3.json"), "w"), indent=1)
