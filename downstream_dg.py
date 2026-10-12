# The author's selection rule applied to the scale of Delta G: for each candidate alignment (the scale mu_c at
# which Delta G = 1/6 exactly), propagate to the scale where the full moment is measured (10 GeV^2, lattice)
# and see whether it creates a mismatch. Leading-order evolution and 4-loop alpha_s, as in deltaG_scale.py.
import importlib.util, io, contextlib
from mpmath import mp, mpf, sqrt, findroot, log, e
spec = importlib.util.spec_from_file_location("dg", "deltaG_scale.py"); dg = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()): spec.loader.exec_module(dg)
sixth = mpf(1)/6; a10 = dg.a10; b0 = dg.b0
lat = {"Zhao": (mpf('0.231'), sqrt(mpf('0.017')**2 + mpf('0.044')**2)), "Yang": (mpf('0.251'), sqrt(mpf('0.047')**2 + mpf('0.016')**2))}
w = {k: 1/s**2 for k, (x, s) in lat.items()}
comb = (sum(w[k]*lat[k][0] for k in lat)/sum(w.values()), 1/sqrt(sum(w.values())))
def implied10(mu_c, DS):
    a = dg.alpha_at(mu_c)
    return (a*sixth + 4*DS/b0*(a - a10))/a10
print("lattice at 10 GeV^2: Zhao %.3f +- %.3f, Yang %.3f +- %.3f, combined as if independent %.4f +- %.4f" % (
    lat["Zhao"][0], lat["Zhao"][1], lat["Yang"][0], lat["Yang"][1], comb[0], comb[1]))
print("%-26s %12s %12s %10s %10s" % ("alignment: DG = 1/6 at", "DG(10) DS=.38", "DG(10) DS=1/3", "vs comb.", "vs Zhao"))
for label, mu in (("mu = 1.2 GeV", mpf('1.2')), ("mu = 1.5 GeV", mpf('1.5')), ("mu^2 = 3 GeV^2 (COMPASS)", sqrt(3)),
                  ("mu = 2 GeV (lattice)", mpf(2)), ("mu = 2.5 GeV", mpf('2.5')), ("mu^2 = 10 GeV^2 (as written)", sqrt(10)), ("mu = 4 GeV", mpf(4))):
    v38 = implied10(mu, mpf('0.38')); v13 = implied10(mu, mpf(1)/3)
    print("%-26s %12.4f %12.4f %+9.1fs %+9.1fs" % (label, v38, v13, (v38 - comb[0])/comb[1], (v38 - lat["Zhao"][0])/lat["Zhao"][1]))
# window of alignments with no mismatch beyond 1 sigma against the combined lattice value (DS = 0.38)
def mu_for(target, DS=mpf('0.38')):
    f = lambda lm: implied10(e**lm, DS) - target
    return e**findroot(f, log(2))
lo, hi = mu_for(comb[0] + comb[1]), mu_for(comb[0] - comb[1])
print("no-mismatch window (within 1 sigma of the combined lattice): mu = %.2f to %.2f GeV; centre %.2f GeV" % (lo, hi, mu_for(comb[0])))
lo2, hi2 = mu_for(comb[0] + 2*comb[1]), mu_for(comb[0] - 2*comb[1])
print("within 2 sigma: mu = %.2f to %.2f GeV" % (lo2, hi2))
# with the thirds exact at 2 GeV (DS = 1/3), the number the EIC tests at 10 GeV^2
print("thirds exact at 2 GeV (DS = 1/3): DG(10 GeV^2) = %.4f; L(2 GeV) = 1/6 by the sum rule" % implied10(mpf(2), mpf(1)/3))
