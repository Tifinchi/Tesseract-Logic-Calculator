# Decide-now tests 3-5 for the Tesseract Logic Calculator (rules fixed before running).
# Test 3, precise constants: of the outside targets known better than 0.1%, how many does the static alphabet hit within 2 sigma?
# Test 4, loads: the main paper's load rule (a static value moved by coupling^2/8pi^2 per loading phase, sign free,
#   anchored on alpha/2pi) applied to every alphabet member with up to three phases of each of three couplings;
#   compare hits on the precise targets with alphabets of the same recipe and other ingredients, given the same loads.
# Test 5, information: the committed numbers' chance probabilities (look-elsewhere test, 2026-10-05), in bits.
import math, itertools, random, bisect
random.seed(20261010 + 5)
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
cube = build(base_from(pi, 2**0.5, 3**0.5, phi))

v = 246.21965
mW, smW, mZ, smZ = 80.3625, 0.0077, 91.1876, 0.0021
PRECISE = [
    ("mW/mZ", mW/mZ, (mW/mZ)*math.hypot(smW/mW, smZ/mZ)),
    ("mZ/v", mZ/v, smZ/v),
    ("mH/v", 125.13/v, 0.11/v),
    ("1/alpha", 137.035999177, 0.000000021),
    ("m_mu/m_e", 206.7682827, 0.0000046),
    ("m_tau/m_mu", 1776.93/105.6583755, 0.09/105.6583755),
    ("m_p/m_e", 1836.15267343, 0.00000011),
    ("mu_p/mu_N", 2.79284734463, 0.00000000082),
    ("|mu_n/mu_p|", 1.91304276/2.79284734463, 0.00000016),
    ("g_A", 1.27641, 0.00056),
    ("r_p/lambdabar_p", 3.9977, 0.0030),
]
def nearest(alpha, x):
    i = bisect.bisect_left(alpha, x)
    c = [alpha[j] for j in (i - 1, i) if 0 <= j < len(alpha)]
    return min(c, key=lambda a: abs(a - x))
def hits(alpha, targets, k=2):
    return [name for name, x, s in targets if abs(nearest(alpha, x) - x) <= k*s]

print("=== Test 3. Precise constants (relative sigma below 0.1%) against the static alphabet")
for name, x, s in PRECISE:
    a = nearest(cube, x)
    print("  %-16s rel sigma %.1e  nearest %.10g  pull %+.3g" % (name, s/x, a, (x - a)/s))
h = hits(cube, PRECISE)
print("  hit within 2 sigma: %d of %d: %s" % (len(h), len(PRECISE), h))

print("\n=== Test 4. With the load rule")
aW = (1/127.951)/0.23122                                 # alpha_2 at mZ = alpha(mZ)/sin^2 theta_W
LOADS = [k*a for a in (1/137.035999177/(2*pi), aW/(2*pi), 0.1180/(2*pi)) for k in (-3, -2, -1, 1, 2, 3)]
print("  loads per phase: alpha/2pi = %.6f, alpha_2/2pi = %.6f, alpha_s/2pi = %.6f; up to 3 phases, either sign" % (
    1/137.035999177/(2*pi), aW/(2*pi), 0.1180/(2*pi)))
def loaded(alpha):
    out = set(alpha)
    for s in alpha:
        for L in LOADS: out.add(s*(1 + L))
    return sorted(out)
cubeL = loaded(cube)
hL = hits(cubeL, PRECISE)
print("  loaded cube alphabet: %d numbers; precise targets hit within 2 sigma: %d: %s" % (len(cubeL), len(hL), hL))
for name, x, s in PRECISE:
    a = nearest(cubeL, x); print("    %-16s nearest loaded %.10g  pull %+.3g" % (name, a, (x - a)/s))
NR = 400
counts = []; counts_static = []
for _ in range(NR):
    A = build(base_from(random.uniform(2.6, 3.7), random.uniform(1.25, 1.6), random.uniform(1.55, 1.95), random.uniform(1.45, 1.8)))
    counts_static.append(len(hits(A, PRECISE)))
    counts.append(len(hits(loaded(A), PRECISE)))
mean = sum(counts)/NR; sd = (sum((c - mean)**2 for c in counts)/NR)**0.5
print("  same recipe, other ingredients, static: mean %.2f, P(>= %d) = %.3f" % (sum(counts_static)/NR, len(h), sum(1 for c in counts_static if c >= len(h))/NR))
print("  same recipe, other ingredients, loaded: mean %.2f (sd %.2f), P(>= %d) = %.3f" % (mean, sd, len(hL), sum(1 for c in counts if c >= len(hL))/NR))

print("\n=== Test 5. Information in the committed numbers (chance probabilities from the look-elsewhere test)")
items = [("15 deg tilt, strict angle set", 0.12), ("15 deg tilt, generous angle set", 1.0), ("the hum, coefficient fixed at 1", 0.052),
         ("the hum, order-one coefficient", 0.57), ("top and Higgs from one load", 0.6), ("r_p = 4 lambdabar_p", 1/6),
         ("Omega_m = pi^2/32 (nine members within 2 sigma)", 1.0), ("Delta Sigma = 1/3", 1.0), ("Delta G = 1/6", 1.0),
         ("Koide Q = 2/3 (Koide 1981)", 1/850), ("lepton phase 2/9 (Brannen 2006)", 1/300)]
for name, p in items:
    print("  %-48s p = %.3g  -> %.1f bits" % (name, p, -math.log2(p)))
def fisher(ps):
    from math import log
    stat = -2*sum(log(p) for p in ps); k = len(ps)
    # chi-square survival with 2k dof (even), closed form
    s = 0.0; term = 1.0
    for i in range(k):
        if i > 0: term *= (stat/2)/i
        s += term
    return math.exp(-stat/2)*s
for label, ps in (("cube's own, strict (tilt, hum fixed, top/Higgs, r_p)", [0.12, 0.052, 0.6, 1/6]),
                  ("cube's own, generous", [1.0, 0.57, 0.6, 1/6]),
                  ("inherited lepton sector (Koide, Brannen)", [1/850, 1/300])):
    p = fisher(ps)
    print("  %-56s combined p = %.3g -> %.1f bits" % (label, p, -math.log2(p)))
print("  thresholds: p = 0.01 is %.1f bits; 3 sigma (p = 0.0027) is %.1f bits; 5 sigma (p = 2.87e-7) is %.1f bits" % (-math.log2(0.01), -math.log2(0.0027), -math.log2(2.87e-7)))
