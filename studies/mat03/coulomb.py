"""MAT-03 (3): low-temperature subband mobilities limited by the positive charge added by gating, for the frozen gate
states of Lee et al. (n_tot = n_H(150 K); 4, 5, 7 and 8 V; compensated and clean BSO). All of the added positive
charge, n_tot + N_dep - sigma_P, as random point charges Ze (Z = 1, 2) in a plane a distance d = 1-32 A inside the
LaInO3 cap, with the cap permittivity eps_cap = 24 (measured), 38 or 50 entering through the dielectric boundary
(transport.mobilities(eps_cap=...)); plus the depleted acceptors of the BSO. Static multisubband screening (eps 20).
Reported per case: subband mobilities, the ratio x of the mean upper-subband mobility to the ground-subband mobility,
the conductivity mobility mu_sigma = sum_i (n_i / n_tot) mu_i, the low-field Hall density and the classical
magnetoresistance between 2 and 7 T. A second set adds charged threading dislocations (1e10-1e11 cm^-2, one charge
per 4.116 A) at d = 4 A, eps_cap = 24. Writes coulomb.json."""
import json
from multiprocessing import Pool
from common3 import *
from qwkit import transport as tr

D = (1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0, 24.0, 32.0)
EPSC = (24.0, 38.0, 50.0)
VGS = ("4V", "5V", "7V", "8V")


def setup(NA, vg):
    p = dict(BASE); p["NA"] = NA; p["T"] = 2.0; p["f_int"] = 1.0
    r0 = solve(5e12, p, extend_to=0.0, nstates=10, dz=0.5)
    sigP = 5e12 + r0["N_dep"]
    nt = NH[vg][-1]
    r = solve(nt, p, extend_to=0.0, nstates=10, dz=0.5)
    return r, nt, nt + r["N_dep"] - sigP


def row(r, nt, n, mu, **kw):
    n = np.asarray(n); mu = np.asarray(mu)
    up = float(np.sum(n[1:] * mu[1:]) / np.sum(n[1:]) / mu[0]) if len(n) > 1 else None
    nH, muH = tr.hall(n, mu)
    return dict(kw, mu=mu.tolist(), n=n.tolist(), x=up, mu_sigma=float(np.sum(n * mu) / nt), nH=nH, muH=muH,
                mr27=tr.magnetoresistance(n, mu, 7.0) - tr.magnetoresistance(n, mu, 2.0))


def job(args):
    NA, vg, epsc = args
    r, nt, dn = setup(NA, vg)
    rows = []
    for Z in (1, 2):
        for d in D:
            o = tr.mobilities(r, eps=20.0, sheets=[(-d, dn / Z, Z)], volume=[(0.0, r["w_d"], NA, -1)], nth=64, eps_cap=epsc)
            rw = row(r, nt, o["n"], o["mu"], Z=Z, d=d)
            rows.append(rw)
            print(f"NA {NA:.0e} {vg} eps_cap {epsc:.0f} Z {Z} d {d:4.1f}: mu_i {np.round(o['mu']).astype(int).tolist()} x "
                  f"{rw['x'] if rw['x'] is None else round(rw['x'], 2)} mu_sigma {rw['mu_sigma']:.0f} n_H/n_tot "
                  f"{rw['nH'] / nt:.3f} MR(7T)-MR(2T) {100 * rw['mr27']:.2f}%", flush=True)
    dis = []
    if epsc == 24.0:
        for Nd in (1e10, 3e10, 1e11):
            for Z in (1, 2):
                o = tr.mobilities(r, eps=20.0, sheets=[(-4.0, dn / Z, Z)], volume=[(0.0, r["w_d"], NA, -1)],
                                  dislocations=(Nd, 1.0, 4.116), nth=64, eps_cap=epsc)
                dis.append(row(r, nt, o["n"], o["mu"], Z=Z, d=4.0, Ndis=Nd))
    return dict(NA=NA, vg=vg, eps_cap=epsc, n_tot=nt, donors=dn, w_d=r["w_d"], rows=rows, dislocations=dis)


if __name__ == "__main__":
    jobs = [(NA, vg, e) for NA in (1e18, 5e19) for vg in VGS for e in EPSC]
    with Pool(int(os.environ.get("NPROC", "2"))) as pool:
        res = pool.map(job, jobs, chunksize=1)
    json.dump(res, open(os.path.join(here, "coulomb.json"), "w"), indent=1)
