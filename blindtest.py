# Numerology or structure? Two tests, rules fixed before running.
#
# Part 1, calibration: does the cube's frozen alphabet (cube-sweep, 2026-10-03; look-elsewhere, 2026-10-05)
#   land on an outside list of measured dimensionless numbers more often than alphabets built by the
#   same recipe from other ingredients? The list is fixed below, before any run, and includes numbers the
#   cube never claimed. Primary statistic: H2 = targets with an alphabet member within 2 sigma.
#   Primary null N3: the same recipe with pi, sqrt2, sqrt3 and phi replaced by random numbers of similar size.
#   Secondary nulls: N1, every ingredient jittered by up to +-25%; N2, only the irrational ingredients jittered.
# Part 2, pre-registration: for each frozen prediction, the chance it passes if the cube is wrong
#   (future result drawn from today's data) and the Bayes factor the future measurement can deliver.
import math, itertools, random, bisect
from statistics import median
random.seed(20261010)
pi = math.pi; phi = (1 + 5**0.5)/2

def build(base):
    vals = sorted(set(round(v, 12) for v in base if v > 0))
    raw = set(vals)
    for a, b in itertools.combinations_with_replacement(vals, 2):
        raw.update((a*b, a/b, b/a))
    xs = sorted(raw); out = [xs[0]]
    for x in xs[1:]:
        if abs(x - out[-1]) > 1e-9*max(1, abs(x)): out.append(x)
    return out

def base_from(P, S2, S3, PHI):
    S5 = 2*PHI - 1
    return [1.0, S2, S3, 2.0, 1/S2, 1/S3, 0.5, 1/3, 0.25, 2/3, 0.75, P/4, P/6, P**2/32, 1 - P/4, 1 - P/6, 1 - P**2/32,
            4/P, 1.5, S3/S2, S3/2, 2*P/3, 1/(2*P), 1/(8*P**2), PHI, 1/PHI, S5, 0.375, 5/3, 3, 4, 6, 8, 12, 16, 24, 27]
CUBE_BASE = base_from(pi, 2**0.5, 3**0.5, phi)
IRRATIONAL = [1, 2, 4, 5, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22, 23, 24, 25, 26]   # indices of the irrational ingredients
cube = build(CUBE_BASE)

# ---- the outside list (value, sigma, source); fixed before running -------------------------------
v = 246.21965                                       # GeV, from G_F
mW, smW, mZ, smZ = 80.3625, 0.0077, 91.1876, 0.0021
deg = pi/180
md_, smd = 4.69, 0.05; mu_, smu = 2.20, 0.07; ms_, sms = 92.74, 0.54
TARGETS = [
    ("mW/mZ", mW/mZ, (mW/mZ)*math.hypot(smW/mW, smZ/mZ), "PDG 2026 W; PDG Z"),
    ("mZ/v", mZ/v, smZ/v, "PDG Z; v from G_F"),
    ("mH/v", 125.13/v, 0.11/v, "PDG 2026"),
    ("mt/v", 172.60/v, 0.27/v, "PDG 2026"),
    ("alpha_s(mZ)", 0.1180, 0.0009, "PDG 2026"),
    ("1/alpha", 137.035999177, 0.000000021, "CODATA 2022"),
    ("m_mu/m_e", 206.7682827, 0.0000046, "CODATA 2022"),
    ("m_tau/m_mu", 1776.93/105.6583755, 0.09/105.6583755, "PDG 2026; CODATA"),
    ("m_u/m_d", mu_/md_, (mu_/md_)*math.hypot(smu/mu_, smd/md_), "PDG 2026 lattice, 2 GeV"),
    ("m_s/m_d", ms_/md_, (ms_/md_)*math.hypot(sms/ms_, smd/md_), "PDG 2026 lattice, 2 GeV"),
    ("|V_us|", 0.22517, 0.00068, "PDG 2026 CKM fit"),
    ("|V_cb|", 0.0407, 0.0013, "PDG 2026"),
    ("|V_ub|", 0.00389, 0.00016, "PDG 2026"),
    ("CKM gamma (rad)", 68*deg, 4*deg, "PDG 2026"),
    ("CKM alpha (rad)", 84.1*deg, 3.35*deg, "PDG 2026, +3.7/-3.0 deg averaged"),
    ("sin2 theta12", 0.3096, 0.0065, "NuFIT 6.1, +0.0057/-0.0073 averaged"),
    ("sin2 theta13", 0.0220, 0.00056, "NuFIT 6.1"),
    ("sin2 theta23", 0.561, 0.014, "NuFIT"),
    ("dm2_21/dm2_31", 7.49e-5/2.513e-3, 0.018*7.49e-5/2.513e-3, "NuFIT 6.1; JUNO (approximate sigma)"),
    ("m_p/m_e", 1836.15267343, 0.00000011, "CODATA 2022"),
    ("mu_p/mu_N", 2.79284734463, 0.00000000082, "CODATA 2022"),
    ("|mu_n/mu_p|", 1.91304276/2.79284734463, 0.00000016, "CODATA 2022"),
    ("g_A", 1.27641, 0.00056, "PERKEO III"),
    ("r_p/lambdabar_p", 3.9977, 0.0030, "CODATA 2022"),
    ("Omega_m", 0.3153, 0.0073, "Planck 2018"),
    ("n_s", 0.9649, 0.0042, "Planck 2018"),
    ("tau_reio", 0.0544, 0.0073, "Planck 2018"),
    ("Omega_c/Omega_b", 5.364, 0.065, "Planck 2018"),
    ("N_eff", 2.99, 0.17, "Planck 2018"),
]

def nearest(alpha, x):
    i = bisect.bisect_left(alpha, x)
    c = [alpha[j] for j in (i - 1, i) if 0 <= j < len(alpha)]
    return min(c, key=lambda a: abs(a - x))

def hits(alpha, targets, k):
    lo, hi = alpha[0], alpha[-1]
    n = 0
    for name, x, s, src in targets:
        if lo <= x <= hi and abs(nearest(alpha, x) - x) <= k*s: n += 1
    return n

print("=== Part 1. Calibration on an outside list")
print("cube alphabet: %d numbers, from %.3g to %.3g" % (len(cube), cube[0], cube[-1]))
inr = [t for t in TARGETS if cube[0] <= t[1] <= cube[-1]]
print("targets: %d listed, %d inside the alphabet's span (the rest are dropped by the rule)" % (len(TARGETS), len(inr)))
print("%-18s %14s %12s %14s %8s %10s" % ("target", "value", "sigma", "nearest", "pull", "expected"))
exp1 = exp2 = 0.0
for name, x, s, src in inr:
    a = nearest(cube, x)
    w = 0.25*x; dens = sum(1 for u in cube if abs(u - x) <= w)/(2*w)
    e2 = 1 - math.exp(-dens*4*s); exp2 += e2; exp1 += 1 - math.exp(-dens*2*s)
    print("%-18s %14.10g %12.3g %14.10g %+8.2f %10.3f" % (name, x, s, a, (x - a)/s, e2))
H1, H2 = hits(cube, inr, 1), hits(cube, inr, 2)
print("cube: H1 = %d (expected from its density %.1f), H2 = %d (expected %.1f)" % (H1, exp1, H2, exp2))

def null_N1():
    return build([b*math.exp(random.uniform(-0.25, 0.25)) for b in CUBE_BASE])
def null_N2():
    out = list(CUBE_BASE)
    for i in IRRATIONAL: out[i] = out[i]*math.exp(random.uniform(-0.25, 0.25))
    return build(out)
def null_N3():
    return build(base_from(random.uniform(2.6, 3.7), random.uniform(1.25, 1.6), random.uniform(1.55, 1.95), random.uniform(1.45, 1.8)))
NR = 2000
for label, gen in (("N3, same recipe, other ingredients (primary)", null_N3), ("N1, every ingredient jittered", null_N1), ("N2, irrational ingredients jittered", null_N2)):
    h1 = []; h2 = []
    for _ in range(NR):
        A = gen(); h1.append(hits(A, inr, 1)); h2.append(hits(A, inr, 2))
    p2 = sum(1 for h in h2 if h >= H2)/NR; p1 = sum(1 for h in h1 if h >= H1)/NR
    print("%-48s mean H2 = %.2f (sd %.2f), P(H2 >= %d) = %.3f; mean H1 = %.2f, P(H1 >= %d) = %.3f" % (
        label, sum(h2)/NR, (sum((h - sum(h2)/NR)**2 for h in h2)/NR)**0.5, H2, p2, sum(h1)/NR, H1, p1))
# decoys: the same alphabet against targets moved by up to +-25%, same relative sigma
d2 = []
for _ in range(NR):
    dec = []
    for name, x, s, src in inr:
        f = math.exp(random.uniform(-0.25, 0.25)); dec.append((name, x*f, s*f, src))
    d2.append(hits(cube, dec, 2))
print("%-48s mean H2 = %.2f, P(H2 >= %d) = %.3f" % ("decoy targets, cube alphabet", sum(d2)/NR, H2, sum(1 for h in d2 if h >= H2)/NR))

print("\n=== Part 2. Pre-registration: what each frozen prediction can decide")
s15 = math.sin(15*deg)**2
PRED = [
    # name, cube value, cube sigma, today, today's sigma below, above, future sigma, decider and date, group
    ("sin2 theta12, 15 deg tilt", (1 - s15)/(3*(1 - s15/3)), 0.0, 0.3096, 0.0073, 0.0057, 0.0016, "JUNO six years, about 2031", "own"),
    ("sin2 theta13, 15 deg tilt", s15/3, 0.0, 0.0220, 0.00056, 0.00056, 0.0004, "reactor and global fits, about 2030 (sigma assumed)", "own"),
    ("CKM alpha = 90 deg", 90.0, 0.0, 84.1, 3.0, 3.7, 2.0, "Belle II and LHCb, about 2030", "own"),
    ("Omega_m = pi^2/32 (flat LCDM)", pi**2/32, 0.0, 0.3045, 0.0032, 0.0032, 0.003, "Euclid and final DESI, 2027-30", "own"),
    ("Delta G = 1/6 at 2 GeV", 1/6, 0.0, 0.158, 0.039, 0.039, 0.016, "EIC, 2030s (about 10%), evolved to 2 GeV", "own"),
    ("Delta Sigma = 1/3", 1/3, 0.0, 0.352, 0.029, 0.029, 0.050, "EIC, 2030s (about 15%)", "own"),
    ("top pole mass 172.57 +- 0.38 GeV", 172.57, 0.38, 172.1, 0.6, 0.6, 0.035, "FCC-ee threshold scan, 2040s", "own"),
    ("r_p = 4 lambdabar_p = 0.841236 fm", 4*0.2103089, 0.0, 0.84075, 0.00064, 0.00064, 0.0002, "hydrogen and muonic hydrogen, about 2030 (sigma assumed)", "prior art"),
    ("m_tau = 1776.969 MeV (Koide)", 1776.969, 0.0, 1776.93, 0.09, 0.09, 0.05, "Belle II, BESIII, STCF (sigma assumed)", "prior art"),
]
def lnN(x, m, s): return -0.5*((x - m)/s)**2 - math.log(s) - 0.5*math.log(2*pi)
def sn_draw(mu, sl, sh):
    z = abs(random.gauss(0, 1))
    return mu + (z*sh if random.random() < sh/(sl + sh) else -z*sl)
def std_pred_logpdf(mu, sl, sh, sf):
    # density of today's split normal convolved with the future measurement's error, on a grid
    lo, hi = mu - 9*max(sl, sf), mu + 9*max(sh, sf)
    n = 6001; h = (hi - lo)/(n - 1)
    xs = [lo + i*h for i in range(n)]
    norm = 2/((2*pi)**0.5*(sl + sh))
    w = [norm*math.exp(-0.5*((x - mu)/(sl if x < mu else sh))**2)*h for x in xs]
    def f(m):
        return math.log(max(1e-300, sum(wi*math.exp(-0.5*((m - x)/sf)**2) for x, wi in zip(xs, w))/((2*pi)**0.5*sf)))
    return f
NS = 6000
own_cube, own_std = [0.0]*NS, [0.0]*NS
rows = []
L = math.log(100)
for name, c, sc, mu, sl, sh, sf, dec, grp in PRED:
    s_c = math.hypot(sc, sf)
    side = sh if c > mu else sl
    pull = (mu - c)/math.hypot(side, sc)
    win = 2*s_c
    lp = std_pred_logpdf(mu, sl, sh, sf)
    lb_c = []; lb_s = []; npass = 0
    for i in range(NS):
        m = random.gauss(c, s_c); b = lnN(m, c, s_c) - lp(m); lb_c.append(b)
        if grp == "own": own_cube[i] += b
        m = sn_draw(mu, sl, sh) + random.gauss(0, sf); b = lnN(m, c, s_c) - lp(m); lb_s.append(b)
        if abs(m - c) <= win: npass += 1
        if grp == "own": own_std[i] += b
    p_std = npass/NS
    rows.append((name, c, pull, p_std, 0.9545/max(p_std, 1e-9), median(lb_c), sum(1 for b in lb_c if b > L)/NS, median(lb_s), sum(1 for b in lb_s if b < -L)/NS, dec, grp))
print("%-34s %10s %7s %8s %8s %9s %9s %9s %9s" % ("prediction", "value", "pull", "P(pass|", "B if", "med lnB", "P(B>100", "med lnB", "P(B<1/100"))
print("%-34s %10s %7s %8s %8s %9s %9s %9s %9s" % ("", "", "today", "wrong)", "pass", "if right", "|right)", "if wrong", "|wrong)"))
for r in rows:
    print("%-34s %10.6g %+7.2f %8.3f %8.1f %9.2f %9.3f %9.2f %9.3f  %s [%s]" % r)
print("cube's own seven together: if right, median ln B = %.1f (B = %.3g), P(B > 100) = %.3f; if wrong, median ln B = %.1f (B = %.3g), P(B < 1/100) = %.3f" % (
    median(own_cube), math.exp(median(own_cube)), sum(1 for b in own_cube if b > L)/NS, median(own_std), math.exp(median(own_std)), sum(1 for b in own_std if b < -L)/NS))
pp = {r[0]: r[3] for r in rows}
j = pp["sin2 theta12, 15 deg tilt"]*pp["CKM alpha = 90 deg"]*pp["Delta G = 1/6 at 2 GeV"]
print("passing theta12, alpha and Delta G together if the cube is wrong: %.4f (1 in %.0f)" % (j, 1/j))
