"""MAT-03 shared set-up: the MAT-01/MAT-02 interface model, and the digitized Hall data of Lee et al., Adv. Electron.
Mater. 11, 2400811 (2024), Fig. 3b (Hall density and mobility at 2, 10, 50, 100 and 150 K for gate voltages 0-8 V; the
gate was set at 220 K and the ions are frozen below it, so each gate voltage is one fixed positive charge). The data are
the verified set lee_fig3b_v2.json (digitization/build_v2.py: Fig. 3b marker centroids cross-checked against Fig. 4).
The 6 V markers at 2 K are hidden behind the 7 V markers and are not used; 6 V at 10 and 150 K comes from Fig. 4."""
import os, sys, json
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "mat01")); sys.path.insert(0, os.path.join(here, "..", ".."))
from common import *                                     # noqa: F401,F403
here = os.path.dirname(os.path.abspath(__file__))

_L = json.load(open(os.path.join(here, "lee_fig3b_v2.json")))
TEMPS = _L["T"]
VG = [f"{v}V" for v in range(9)]
NH = {k: [None if x is None else x * 1e13 for x in _L["n_1e13"][k]] for k in VG}   # cm^-2
MU = {k: _L["mu"][k] for k in VG}                                                  # cm^2/Vs
NH["6V"][0] = MU["6V"][0] = None                         # hidden marker
