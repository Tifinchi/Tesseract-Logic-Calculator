# The proton as folded reactions (section 4.6): the moment folded into three reactions,
# and the spin pieces folded to one resolution.
from mpmath import mp, mpf, sqrt, log, zeta
mp.dps = 25

print("=== 1. The moment: unfolded, folded into three reactions, measured")
mup, sup = mpf('2.79284734463'), mpf('0.00000000082')     # CODATA 2022, in nuclear magnetons
mun, sun = mpf('-1.91304276'), mpf('0.00000045')
mpc2 = mpf('938.27208943')                                  # MeV
print("  one reaction spread round the circle 2 pi (4 lam_p) at c: 4 muN, %+.1f%% from the measured %s" % (float((4 - mup)/mup*100), mp.nstr(mup, 12)))
print("  one reaction on its own Compton circle (Dirac, g = 2): 1 muN, %+.1f%%" % float((1 - mup)/mup*100))
# three g = 2 reactions, charges 2/3, 2/3, -1/3 (proton) and 2/3, -1/3, -1/3 (neutron), symmetric spin-flavour state:
# mu_p = (4 mu_u - mu_d)/3, mu_n = (4 mu_d - mu_u)/3, mu_q = (e_q/e)(m_p/m) muN for a reaction of mass m
for label, frac in (("each a third of the proton's mass", mpf(1)/3),):
    k = 1/frac                                              # m_p/m
    mu_u, mu_d = mpf(2)/3*k, -mpf(1)/3*k
    P, N = (4*mu_u - mu_d)/3, (4*mu_d - mu_u)/3
    print("  three folded reactions, %s: mu_p = %s muN (%+.2f%%), mu_n = %s muN (%+.2f%% in size), ratio %s against %s (%+.2f%%)" % (
        label, mp.nstr(P, 6), float((P - mup)/mup*100), mp.nstr(N, 6), float((abs(N) - abs(mun))/abs(mun)*100),
        mp.nstr(N/P, 6), mp.nstr(mun/mup, 9), float((N/P - mun/mup)/(mun/mup)*100)))
print("  measured ratio over the folded -2/3: %+.3f%%" % float((mun/mup/(-mpf(2)/3) - 1)*100))
print("  mass each reaction would need: from mu_p %s MeV; from mu_n %s MeV; a third of m_p is %s MeV" % (
    mp.nstr(mpc2/mup, 6), mp.nstr(mpf(2)/3*mpc2/abs(mun), 6), mp.nstr(mpc2/3, 6)))
print("  the neutron: charge 0, so one reaction has i = q f = 0 and no moment; measured %s muN" % mp.nstr(mun, 9))

print("=== 2. The strong coupling at the two resolutions (4-loop MS-bar running)")
z3 = zeta(3)
def betas(nf):
    b0 = 11 - mpf(2)/3*nf
    b1 = 102 - mpf(38)/3*nf
    b2 = mpf(2857)/2 - mpf(5033)/18*nf + mpf(325)/54*nf**2
    b3 = (mpf(149753)/6 + 3564*z3) - (mpf(1078361)/162 + mpf(6508)/27*z3)*nf + (mpf(50065)/162 + mpf(6472)/81*z3)*nf**2 + mpf(1093)/729*nf**3
    return b0, b1, b2, b3
def run(a_s, mu0, mu1, nf, steps=4000):
    # RK4 in t = ln mu^2 for a = alpha_s/4pi, either direction
    b = betas(nf)
    F = lambda a: -(b[0]*a**2 + b[1]*a**3 + b[2]*a**4 + b[3]*a**5)
    t0, t1 = log(mu0**2), log(mu1**2)
    h = (t1 - t0)/steps
    a = a_s/(4*mp.pi)
    for _ in range(steps):
        k1 = F(a); k2 = F(a + h*k1/2); k3 = F(a + h*k2/2); k4 = F(a + h*k3)
        a += h*(k1 + 2*k2 + 2*k3 + k4)/6
    return a*4*mp.pi
mZ, mb = mpf('91.1876'), mpf('4.18')
out = {}
for key in ('0.1179', '0.1170', '0.1188'):
    aZ = mpf(key)
    ab = run(aZ, mZ, mb, 5)
    a10 = run(ab, mb, sqrt(10), 4)
    a4 = run(ab, mb, mpf(2), 4)
    out[key] = (a10, a4)
    print("  alpha_s(mZ) = %s: alpha_s(sqrt(10) GeV) = %s, alpha_s(2 GeV) = %s" % (key, mp.nstr(a10, 5), mp.nstr(a4, 5)))

print("=== 3. Gluon helicity folded from 10 GeV^2 to 4 GeV^2 (leading order)")
# LO: d(DSigma)/dt = 0; d(DG)/dt = (a_s/2pi)[(3/2)C_F DSigma + (b0/2) DG]  (Ji, Tang, Hoodbhoy, eq. 6)
# with d a_s/dt = -b0 a_s^2/4pi:  a_s(t) DG(t) = a_s(t0) DG(t0) + (4 DSigma/b0)(a_s(t0) - a_s(t))
b0 = betas(4)[0]
DG10, sDG10 = mpf('0.231'), sqrt(mpf('0.017')**2 + mpf('0.044')**2)
hS, shS = mpf('0.190'), mpf('0.022')                         # ETMC 1/2 DSigma at 2 GeV
DS, sDS = 2*hS, 2*shS
def fold(a10, a4, DS=DS):
    return (a10*DG10 + 4*DS/b0*(a10 - a4))/a4
a10, a4 = out['0.1179']
a4_lo = 1/(1/a10 + b0/(4*mp.pi)*log(mpf(4)/10))             # one-loop step from the same a_s(sqrt 10)
for name, A4 in (("4-loop alpha_s", a4), ("one-loop step", a4_lo)):
    DG4 = fold(a10, A4)
    dG_dDG = a10/A4; dG_dDS = 4/b0*(a10 - A4)/A4
    sDG4 = sqrt((dG_dDG*sDG10)**2 + (dG_dDS*sDS)**2)
    L = mpf(1)/2 - hS - DG4
    sL = sqrt((dG_dDG*sDG10)**2 + ((mpf(1)/2 + dG_dDS)*sDS)**2)
    sixth = mpf(1)/6
    print("  %s (alpha_s(2 GeV) = %s): DG(2 GeV) = %s +- %s, %+.2f sigma from 1/6; L = 1/2 - 1/2 DSigma - DG = %s +- %s, %+.2f sigma" % (
        name, mp.nstr(A4, 4), mp.nstr(DG4, 3), mp.nstr(sDG4, 2), float((DG4 - sixth)/sDG4), mp.nstr(L, 3), mp.nstr(sL, 2), float((L - sixth)/sL)))
for aZ in ('0.1170', '0.1188'):
    A10, A4 = out[aZ]
    DG4 = fold(A10, A4)
    print("  alpha_s(mZ) = %s: DG(2 GeV) = %s, L = %s" % (aZ, mp.nstr(DG4, 3), mp.nstr(mpf(1)/2 - hS - DG4, 3)))
print("  with HERMES DSigma = 0.330 for the quark term: DG(2 GeV) = %s" % mp.nstr(fold(a10, a4, mpf('0.330')), 3))
print("  the mixed-resolution remainder, for comparison: L = 0.079 +- 0.052, %+.2f sigma" % float((mpf('0.079') - mpf(1)/6)/mpf('0.052')))

print("=== 4. What folding keeps")
print("  momentum sum at 2 GeV (ETMC): 0.995 +- %s, %+.2f sigma from 1" % (mp.nstr(sqrt(mpf('0.060')**2 + mpf('0.029')**2), 3), float((mpf('0.995') - 1)/sqrt(mpf('0.060')**2 + mpf('0.029')**2))))
print("  spin sum at 2 GeV (ETMC): 0.507 +- %s, %+.2f sigma from 1/2" % (mp.nstr(sqrt(mpf('0.043')**2 + mpf('0.065')**2), 3), float((mpf('0.507') - mpf('0.5'))/sqrt(mpf('0.043')**2 + mpf('0.065')**2))))
for nf in (3, 4, 5, 6):
    print("  infinite resolution, nf = %d: quarks %s, gluons %s of the spin and of the momentum" % (nf, mp.nstr(mpf(3*nf)/(16 + 3*nf), 4), mp.nstr(mpf(16)/(16 + 3*nf), 4)))

print("=== 5. Digits for the text")
for name, A10, A4 in (("4-loop", a10, a4), ("one-loop step", a10, a4_lo), ("0.1170", *out['0.1170']), ("0.1188", *out['0.1188'])):
    DG4 = fold(A10, A4)
    dG_dDG = A10/A4; dG_dDS = 4/b0*(A10 - A4)/A4
    sDG4 = sqrt((dG_dDG*sDG10)**2 + (dG_dDS*sDS)**2)
    L = mpf(1)/2 - hS - DG4
    sL = sqrt((dG_dDG*sDG10)**2 + ((mpf(1)/2 + dG_dDS)*sDS)**2)
    print("  %-14s DG = %s +- %s (%s sigma); L = %s +- %s (%s sigma)" % (name, mp.nstr(DG4, 5), mp.nstr(sDG4, 4), mp.nstr((DG4 - mpf(1)/6)/sDG4, 4), mp.nstr(L, 5), mp.nstr(sL, 4), mp.nstr((L - mpf(1)/6)/sL, 4)))
print("  1/2 DSigma = 0.190 +- 0.022: %s sigma" % mp.nstr((hS - mpf(1)/6)/shS, 4))
