"""Appendix G: 'how slight' candidates and the gravity-boundary table, PDG 2026 edges."""
import sys, math
sys.path.insert(0,'/tmp/claude-0/-home-claude/bbcd54b6-3962-59fe-98a5-d2cef4b7e16c/scratchpad')
from higgs2 import run, beta, MPl
from scipy.optimize import brentq
edgeH, sEH = 128.548, 0.964   # exact-edge Higgs at top 172.60 (higgs_rule2026.py)
edgeT, sET = 170.935, 0.375   # exact-edge top at Higgs 125.13
mH, smH = 125.13, 0.11; mt, smt = 172.60, 0.27
ALS, MW = 0.1180, 80.3625
ang = lambda r: math.degrees(math.atan(r))
print(f"vacuum: Higgs-based {edgeH/mH:.4f} -> {ang(edgeH/mH):.2f} deg; top-based {mt/edgeT:.4f} -> {ang(mt/edgeT):.2f} deg")
obs = [
 ("Matter over antimatter", 1+6.1e-10),
 ("Curvature of the universe", 1+0.0007),
 ("Neutron over proton mass", 939.56542052/938.27208816),
 ("Electron magnet g/2", 1.00115965218059),
 ("Muon magnet g/2", 1.001165920705),
 ("Granite fault, rate and state", 0.73/(0.73-0.005)),
 ("Teflon on steel", 1.0),
 ("Bone with joint fluid", 0.016/0.015),
 ("Dry sand pile", math.tan(math.radians(23.4))/math.tan(math.radians(22.0))),
 ("Spectral tilt 1/n_s", 1/0.9649),
 ("Cube held/action sqrt2", math.sqrt(2)),
 ("alpha/2pi", 1+(1/137.035999)/(2*math.pi)),
 ("lambda/2pi (tree 0.12914)", 1+0.12914/(2*math.pi)),
 ("alpha_s(M_t)/pi (0.1080)", 1+0.10803/math.pi),
 ("y_t^2/16pi^2 (y_t 0.9330)", 1+0.93303**2/(16*math.pi**2)),
]
h1H=h1T=h2H=h2T=0
for n, r in obs:
    pH = edgeH/r; sH = math.hypot(smH, sEH/r); pullH = (mH-pH)/sH
    pT = edgeT*r; sT = math.hypot(smt, sET*r); pullT = (mt-pT)/sT
    h1H += abs(pullH)<1; h1T += abs(pullT)<1; h2H += abs(pullH)<2; h2T += abs(pullT)<2
    print(f"{n:32s} ratio {r:.6f} angle {ang(r):7.3f} | Higgs {pH:7.2f} ({pullH:+5.1f}) | top {pT:7.2f} ({pullT:+5.1f})")
print(f"within 1 sigma: Higgs {h1H}/15, top {h1T}/15; within 2: {h2H}, {h2T}")
# gravity boundary table
Mbar = 2.435e18
def lam_at(Mt, Mh, mu):
    s = run(Mt, MPl, Mh=Mh, als=ALS, MW=MW); return s.sol(math.log(mu/Mt))
y = lam_at(mt, 128.548, Mbar); y0 = list(y); y0[4] = 0.0; B = beta(0, y0)[4]
print(f"B = {B:.3e}")
for A in (1.0, 0.3, 0.1, 0.03):
    target = -B/A
    m = brentq(lambda M: lam_at(mt, M, Mbar)[4]-target, 118, 140, xtol=1e-6)
    print(f"A={A}: lambda(boundary)={target:+.2e} -> m_H {m:.2f}")
lam_meas = lam_at(mt, mH, Mbar)[4]
print(f"measured lambda at Mbar: {lam_meas:.5f} -> needs A = {B/(-lam_meas):.2e}")
