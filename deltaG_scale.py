# Where does the gluon spin equal 1/6?  Leading-order evolution of Delta G (Ji, Tang and Hoodbhoy 1996),
# the same evolution used for the folding in section 5.3, with 4-loop alpha_s running (PDG 2026 alpha_s).
from mpmath import mp, mpf, sqrt, log, zeta, findroot
mp.dps = 30
z3 = zeta(3)
def betas(nf):
    b0 = 11 - mpf(2)/3*nf
    b1 = 102 - mpf(38)/3*nf
    b2 = mpf(2857)/2 - mpf(5033)/18*nf + mpf(325)/54*nf**2
    b3 = (mpf(149753)/6 + 3564*z3) - (mpf(1078361)/162 + mpf(6508)/27*z3)*nf + (mpf(50065)/162 + mpf(6472)/81*z3)*nf**2 + mpf(1093)/729*nf**3
    return b0, b1, b2, b3
def run(a_s, mu0, mu1, nf, steps=3000):
    b = betas(nf); F = lambda a: -(b[0]*a**2 + b[1]*a**3 + b[2]*a**4 + b[3]*a**5)
    t0, t1 = log(mu0**2), log(mu1**2); h = (t1 - t0)/steps; a = a_s/(4*mp.pi)
    for _ in range(steps):
        k1 = F(a); k2 = F(a + h*k1/2); k3 = F(a + h*k2/2); k4 = F(a + h*k3); a += h*(k1 + 2*k2 + 2*k3 + k4)/6
    return a*4*mp.pi
mZ, mb = mpf('91.1876'), mpf('4.18')
def alpha_at(mu, aZ=mpf('0.1180')):
    ab = run(aZ, mZ, mb, 5)
    return run(ab, mb, mu, 4) if mu < mb else run(aZ, mZ, mu, 5)
b0 = betas(4)[0]
sixth = mpf(1)/6
a10 = alpha_at(sqrt(10))
print("alpha_s(sqrt10 GeV) = %s, alpha_s(2 GeV) = %s (alpha_s(mZ) = 0.1180)" % (mp.nstr(a10, 5), mp.nstr(alpha_at(mpf(2)), 5)))
def DG_at(mu, DG10, DS):
    a = alpha_at(mu)
    return (a10*DG10 + 4*DS/b0*(a10 - a))/a
def mu_star(DG10, DS):
    # a(mu*) = a10 (DG10 + 4DS/b0)/(1/6 + 4DS/b0)
    target = a10*(DG10 + 4*DS/b0)/(sixth + 4*DS/b0)
    f = lambda lm: alpha_at(mp.e**lm) - target
    lm = findroot(f, log(2))
    return mp.e**lm
inputs = [("Zhao et al. 2025 (lattice)", mpf('0.231'), sqrt(mpf('0.017')**2 + mpf('0.044')**2)),
          ("Yang et al. 2017 (lattice)", mpf('0.251'), sqrt(mpf('0.047')**2 + mpf('0.016')**2))]
for DS, dslabel in ((mpf('0.380'), "Delta Sigma = 0.380 (ETMC, 2 x 0.190)"), (mpf('0.330'), "Delta Sigma = 0.330 (HERMES)")):
    print("\n--", dslabel)
    for name, DG10, s in inputs:
        pull = (sixth - DG10)/s
        ms = mu_star(DG10, DS); lo = mu_star(DG10 + s, DS); hi = mu_star(DG10 - s, DS)
        print("  %-28s DG(10 GeV^2) = %s +- %s; pull (DG - 1/6)/sigma = %+.2f; DG = 1/6 at mu = %s GeV (mu^2 = %s GeV^2), 1-sigma range %s to %s GeV" % (
            name, DG10, mp.nstr(s, 3), float(-pull), mp.nstr(ms, 3), mp.nstr(ms**2, 3), mp.nstr(lo, 3), mp.nstr(hi, 3)))
print("\n-- Delta G(mu) from Zhao et al. (Delta Sigma = 0.380):")
for mu2 in (mpf(2), mpf(4), mpf("4.5"), mpf(5), mpf(10)):
    mu = sqrt(mu2)
    if mu < mpf('1.3'): continue
    print("  mu^2 = %5s GeV^2: Delta G = %s" % (mp.nstr(mu2, 3), mp.nstr(DG_at(mu, mpf('0.231'), mpf('0.380')), 4)))
# DSSV14 truncated moment, x > 0.05, 90% CL -> 1 sigma
up, dn = mpf('0.06')/mpf('1.6449'), mpf('0.07')/mpf('1.6449')
print("\nDSSV14 truncated (x>0.05, 10 GeV^2): 0.20 +%s/-%s (1 sigma); pull (value - 1/6)/sigma = %+.2f, truncated, not the full moment" % (mp.nstr(up, 2), mp.nstr(dn, 2), float((mpf('0.20') - sixth)/dn)))
print("NNPDF reweighted with STAR 2009 jets (0.05<x<0.2): 0.17 +- 0.06; pull (value - 1/6)/sigma = %+.2f, truncated" % float((mpf('0.17') - sixth)/mpf('0.06')))
print("COMPASS <dg/g> = 0.113 +- %s at x_g ~ 0.10, mu^2 = 3 GeV^2 (a local ratio, not a moment)" % mp.nstr(sqrt(mpf('0.038')**2 + mpf('0.036')**2), 3))
