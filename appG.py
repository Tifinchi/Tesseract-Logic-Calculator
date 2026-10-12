"""Appendix G calculations, PDG 2026 inputs (two-loop running + three-loop QCD, Buttazzo et al. 2013 matching; higgs2.py)."""
import sys, math, numpy as np
sys.path.insert(0,'/tmp/claude-0/-home-claude/bbcd54b6-3962-59fe-98a5-d2cef4b7e16c/scratchpad')
from higgs2 import summary, inputs, run, beta, MPl
from scipy.optimize import brentq, fsolve
MH, MT, ALS, MW = 125.13, 172.60, 0.1180, 80.3625
k = 1/(16*math.pi**2)
v = 246.21965
print("== G.1 basics")
print(f"lambda tree = mH^2/2v^2 = {MH**2/(2*v*v):.5f}; v/2 = {v/2:.3f}; v/sqrt2 = {v/math.sqrt(2):.3f}; 2mW = {2*MW:.3f}")
print(f"v/2 vs MH: {100*(MH/(v/2)-1):.3f}% above; pull {(MH-v/2)/0.11:.1f} sigma")
g1,g2,g3,yt,l = inputs(Mt=MT, Mh=MH, als=ALS, MW=MW)
print(f"MSbar at mu=Mt: g1={g1:.5f} g2={g2:.5f} g3={g3:.5f} yt={yt:.5f} lambda={l:.5f}")
s = run(MT, MPl, Mh=MH, als=ALS, MW=MW)
def at(mu): return s.sol(math.log(mu/MT))
print("== G.2 one-loop pieces (units 1/16pi^2) and net per e-fold")
for mu in (1e3, 1e10, 1e18):
    G1,G2,G3,YT,L = at(mu); a2=G2*G2; gp=0.6*G1*G1; T=YT*YT
    top = -6*T*T; carriers = 0.375*(2*a2*a2+(a2+gp)**2); selfmix = 24*L*L + L*(12*T-9*a2-3*gp)
    print(f"mu={mu:.0e}: top {top:+.3f} carriers {carriers:+.3f} self+mix {selfmix:+.3f} net/e-fold {k*(top+carriers+selfmix):+.5f}  (lambda={L:.5f}, yt={YT:.4f})")
print(f"e-folds Mt->MPl: {math.log(MPl/MT):.2f}")
print("== G.3 lambda vs energy")
e_top = brentq(lambda M: summary(Mt=M, Mh=MH, als=ALS, MW=MW)['lmin'], 165, 178, xtol=1e-6)
se = run(e_top, MPl, Mh=MH, als=ALS, MW=MW)
row1=[]; row2=[]
for p in range(3, 20, 2):
    mu = 10.0**p
    row1.append(s.sol(math.log(mu/MT))[4]); row2.append(se.sol(math.log(mu/e_top))[4] if mu<=MPl else float('nan'))
print("energies 1e3..1e19 step 1e2")
print("measured:", " ".join(f"{x:+.4f}" for x in row1))
print(f"edge top {e_top:.3f}:", " ".join(f"{x:+.4f}" for x in row2))
r = summary(Mt=MT, Mh=MH, als=ALS, MW=MW)
print(f"crossing {r['cross']:.3e}, min {r['lmin']:.5f} at {r['mumin']:.3e}, lambda(MPl=1.22e19) {r['lPl']:.5f}")
print("== code check vs Buttazzo et al. 2013 (Mt=173.34, Mh=125.15, als=0.1184, MW=80.384)")
rb = summary(Mt=173.34, Mh=125.15, als=0.1184, MW=80.384)
print(f"lambda(MPl) = {rb['lPl']:.4f} (Buttazzo -0.0143)")
mh_stab = brentq(lambda M: summary(Mt=173.34, Mh=M, als=0.1184, MW=80.384)['lmin'], 120, 140, xtol=1e-5)
print(f"stability bound Mh > {mh_stab:.2f} (Buttazzo 129.6 +- 1.5)")
print("== G.5 gravity boundary lambda(Mbar_Pl)=0, Mbar_Pl = 2.435e18")
Mbar = 2.435e18
def lam_at(Mt, Mh, mu, als=ALS):
    ss = run(Mt, MPl, Mh=Mh, als=als, MW=MW); return ss.sol(math.log(mu/Mt))
mh_grav = brentq(lambda M: lam_at(MT, M, Mbar)[4], 118, 140, xtol=1e-6)
print(f"Higgs from lambda(Mbar)=0 at Mt={MT}: {mh_grav:.3f}")
mh_grav_full = brentq(lambda M: lam_at(MT, M, MPl)[4], 118, 140, xtol=1e-6)
print(f"  (full Planck mass 1.22e19: {mh_grav_full:.3f})")
y = lam_at(MT, mh_grav, Mbar); y0 = list(y); y0[4] = 0.0
B = beta(0, y0)[4]
print(f"B = SM running at lambda=0 at Mbar: {B:.3e} per e-fold")
mh_sw = brentq(lambda M: lam_at(171.3, M, Mbar)[4], 118, 140, xtol=1e-6)
print(f"S&W validation at Mt=171.3: {mh_sw:.2f} (published 126)")
dmh = (brentq(lambda M: lam_at(MT+0.5, M, Mbar)[4], 118, 140, xtol=1e-6) - brentq(lambda M: lam_at(MT-0.5, M, Mbar)[4], 118, 140, xtol=1e-6))
print(f"slope dMh/dMt of the gravity boundary: {dmh:.3f}")
mt_land = brentq(lambda M: lam_at(M, MH, Mbar)[4], 165, 178, xtol=1e-6)
print(f"top that lands the gravity boundary on Mh={MH}: {mt_land:.3f}")
print("== G.5 double root at MPl (lambda = beta_lambda = 0)")
def F(x):
    Mt, Mh = x
    ss = run(Mt, MPl, Mh=Mh, als=ALS, MW=MW); yy = ss.sol(math.log(MPl/Mt))
    return [yy[4]*100, beta(0, yy)[4]*1e4]
sol = fsolve(F, [174.3, 132.1], xtol=1e-10)
print(f"double root: Mt={sol[0]:.3f} Mh={sol[1]:.3f}; residual {F(sol)}")
print("== G.6 bounce (tree level), flat and with gravity")
tauU = 13.8e9*3.156e7/6.582e-25  # in GeV^-1
mus = np.logspace(10, 19.08, 4000)
lam = np.array([s.sol(math.log(m/MT))[4] for m in mus])
neg = lam < 0
S0 = np.where(neg, 8*math.pi**2/(3*np.abs(lam)), np.inf)
marg0 = S0 - 4*np.log(tauU*mus)
i = int(np.argmin(marg0))
print(f"flat: S={S0[i]:.0f} at mu={mus[i]:.2e}; margin S-4ln(tauU mu)={marg0[i]:.0f}")
dS = np.where(neg, 256*math.pi**3*mus**2/(45*Mbar**2*lam**2), np.inf)
S1 = S0 + dS; marg1 = S1 - 4*np.log(tauU*mus); j = int(np.argmin(marg1))
print(f"gravity: S={S1[j]:.0f} at mu={mus[j]:.2e}; margin={marg1[j]:.0f}")
# lifetime estimate: Gamma/V ~ mu^4 e^{-S}; tau ~ (e^{S}/mu^4)^{1/4} ... use p = (tauU mu)^4 e^{-S} -> tau = tauU * e^{margin/4}
for lab, m in (("flat", marg0[i]), ("gravity", marg1[j])):
    print(f"  {lab}: tau ~ tauU*exp(margin/4) = 10^{math.log10(13.8e9)+m/4/math.log(10):.0f} years")
# top mass at which the margin reaches zero (already decayed)
def margin_flat(Mt, grav=False):
    ss = run(Mt, MPl, Mh=MH, als=ALS, MW=MW)
    lam = np.array([ss.sol(math.log(m/Mt))[4] for m in mus]); ng = lam < 0
    if not ng.any(): return 1e9
    S = np.where(ng, 8*math.pi**2/(3*np.abs(lam)), np.inf)
    if grav: S = S + np.where(ng, 256*math.pi**3*mus**2/(45*Mbar**2*lam**2), np.inf)
    return float(np.min(S - 4*np.log(tauU*mus)))
mt_dec = brentq(lambda M: margin_flat(M), 174, 185, xtol=1e-3)
mt_dec_g = brentq(lambda M: margin_flat(M, True), 174, 185, xtol=1e-3)
print(f"decayed already if Mt > {mt_dec:.2f} (flat), {mt_dec_g:.2f} (with gravity)")
print("== G.4 load x bounce identity: (|lambda_min|/8pi^2) * S0 =", abs(r['lmin'])/(8*math.pi**2)*8*math.pi**2/(3*abs(r['lmin'])))
