# Is this how the universe works?  One rule: an outcome stays only if it matches observed behaviour.
# Applied to the Tesseract Logic Calculator and to the standard theories alike (section 22 of the paper).
#   A. the census: which observed behaviours each theory computes, takes as input, misses or leaves silent
#   B. the fit of what each commits to, and what survival does to a fit
#   C. downstream: what the calculator's own numbers force elsewhere
#   D. the calendar: when the frozen predictions can decide, under the frozen scoring rule
#   E. the verdict: the author's four cores at three levels, synthesised by her rule
# Data as of 10 October 2026 (CODATA 2022, PDG 2026, NuFIT 6.0 and 6.1, DESI DR2), the numbers of Tables 63, 74, 78 and 81.
import math
import numpy as np
from scipy.stats import chi2, norm
from scipy.optimize import brentq

print('=== A. The census of observed behaviours')
# classes, calculator: same (restates the standard result), part (first term only), own (its own number), input, silent, fails
# classes, standard theories: computes, input, fails, silent.   x marks an exposure (2 to 5 sigma against a precise data set)
CENSUS = [
 ('gravity', 'orbits, light bending and clock rates near masses', 'same', 'computes', ''),
 ('gravity', 'gyroscopes in orbit', 'same', 'computes', ''),
 ('gravity', 'gravitational waves', 'same', 'computes', ''),
 ('gravity', 'black holes: horizons and shadows', 'same', 'computes', ''),
 ('gravity', 'antimatter falls at g', 'same', 'computes', ''),
 ('quantum', "hydrogen's lines", 'same', 'computes', ''),
 ('quantum', 'the Lamb shift', 'same', 'computes', ''),
 ('quantum', "the electron's magnetic moment", 'part', 'computes', ''),
 ('quantum', 'spin 1/2, its 720 degree turn, antiparticles', 'same', 'computes', ''),
 ('quantum', 'shells and the periodic table', 'same', 'computes', ''),
 ('particles', "one generation's charges, anomaly-free", 'same', 'input', ''),
 ('particles', 'the W and Z masses', 'same', 'computes', ''),
 ('particles', 'the Higgs mass', 'own', 'input', ''),
 ('particles', 'the top mass', 'own', 'input', ''),
 ('particles', 'the tau mass from the electron and muon', 'own', 'input', ''),
 ('particles', 'the neutrino masses', 'own', 'fails', 'C'),
 ('particles', 'neutrino mixing theta12 and theta13', 'own', 'input', ''),
 ('particles', 'theta23 and the neutrino CP phase', 'silent', 'input', ''),
 ('particles', 'quark mixing', 'input', 'fails', ''),
 ('particles', 'quark masses, alpha, the strong phase, three generations', 'input', 'input', ''),
 ('particles', 'lepton universality in B decays', 'silent', 'fails', ''),
 ('strong', "the proton's spin pieces", 'own', 'computes', ''),
 ('strong', "the proton's radius", 'own', 'computes', ''),
 ('strong', "the proton's mass against the electron's", 'silent', 'computes', ''),
 ('cosmos', "the vacuum's energy density", 'own', 'fails', ''),
 ('cosmos', 'the matter share and the acceleration', 'own', 'input', 'C'),
 ('cosmos', 'w = -1 and flat space', 'same', 'input', 'CS'),
 ('cosmos', 'the expansion rate at z = 0.35, E(z)', 'own', 'computes', 'CS'),
 ('cosmos', 'the expansion rate H0', 'fails', 'fails', ''),
 ('cosmos', 'helium and deuterium from the first minutes', 'silent', 'computes', ''),
 ('cosmos', 'lithium-7 from the first minutes', 'silent', 'fails', ''),
 ('cosmos', 'the number of light relics, N_eff', 'silent', 'computes', ''),
 ('cosmos', "the microwave background's acoustic peaks", 'same', 'computes', ''),
 ('cosmos', 'what dark matter is; why matter outnumbers antimatter', 'silent', 'silent', ''),
 ('cosmos', 'the primordial tilt n_s', 'silent', 'input', 'S'),
 ('critical', 'three-dimensional critical exponents', 'silent', 'computes', ''),
]
N = len(CENSUS)
from collections import Counter
cC = Counter(r[2] for r in CENSUS); cS = Counter(r[3] for r in CENSUS)
print('behaviours counted:', N)
print('calculator:', dict(cC), ' exposed:', sum('C' in r[4] for r in CENSUS))
print('standard  :', dict(cS), ' exposed:', sum('S' in r[4] for r in CENSUS))
cross = Counter((r[2], r[3]) for r in CENSUS)
for k in sorted(cross): print('   calculator %-6s / standard %-8s %2d' % (k[0], k[1], cross[k]))
with_outcome_C = cC['same'] + cC['part'] + cC['own']
with_outcome_S = cS['computes'] + cS['fails']
with_outcome_C += cC['fails']
print('behaviours with a computed outcome: calculator %d of %d (own %d, failing %d); standard %d of %d (matching %d, failing %d)' % (
    with_outcome_C, N, cC['own'], cC['fails'], with_outcome_S, N, cS['computes'], cS['fails']))
print('silent or input: calculator %d; standard %d' % (cC['silent'] + cC['input'], cS['silent'] + cS['input']))
print('where the calculator is silent, the standard theories compute %d, fail %d, take %d as input, are silent on %d' % (
    cross[('silent', 'computes')], cross[('silent', 'fails')], cross[('silent', 'input')], cross[('silent', 'silent')]))
print('where the calculator gives its own number, the standard theories take %d as input, fail %d, compute %d' % (
    cross[('own', 'input')], cross[('own', 'fails')], cross[('own', 'computes')]))

print('\n=== B. The fit of what each theory commits to')
def pull(meas, comp, s_lo, s_hi=None, s_c=0.0):
    s_hi = s_lo if s_hi is None else s_hi
    s = s_hi if comp > meas else s_lo          # the side of the measurement's error that faces the computed value
    return (meas - comp)/math.hypot(s, s_c)
s15 = math.sin(math.radians(15))**2
# Table 63, primary measurement of each row; (name, pull, own?). The lepton angle and the lepton phase are computed from
# the same three masses as the tau row and test the same measurement, the tau mass: the lepton sector is counted once.
C_ROWS = [
 ('tau mass', pull(1776.93, 1776.969, 0.09), False),
 ('top pole (direct)', pull(172.60, 172.57, 0.27, s_c=0.38), True),
 ('Higgs', pull(125.13, 125.20, 0.11, s_c=0.92), True),
 ('sin2 theta13', pull(0.02200, s15/3, 0.00056), True),
 ('sin2 theta12', pull(0.3096, (1 - s15)/(3*(1 - s15/3)), 0.0073, 0.0057), True),
 ('CKM alpha', pull(84.1, 90.0, 3.0, 3.7), False),
 ('light bending gamma', 2.1/2.3, False),
 ('geodetic', pull(6601.8, 6604.8, 18.3), False),
 ('frame dragging', pull(37.2, 39.1, 7.2), False),
 ('antihydrogen g', pull(0.75, 1.00, 0.21), False),
 ('Delta Sigma', pull(0.330, 1/3, 0.039), True),
 ('Delta G (lattice)', pull(0.158, 1/6, 0.039), True),
 ('proton radius / lambdabar', pull(3.9977, 4.0, 0.0030), False),
 ('vacuum energy scale', pull(2.24, 2.02, 0.0, s_c=0.20), True),
 ('curvature', pull(0.0007, 0.0, 0.0019), False),
 ('matter share', pull(0.3045, math.pi**2/32, 0.0032), True),
 ('q0', pull(-0.537, 1.5*math.pi**2/32 - 1, 0.026), False),
 ('w', pull(-0.995, -1.0, 0.019), False),
]
for n_, p_, o_ in C_ROWS: print('   %-28s %+6.2f %s' % (n_, p_, 'own' if o_ else ''))
def fit(rows, label):
    x = sum(p*p for p in rows); k = len(rows)
    print('%-58s chi2 = %6.2f for %2d, chi2/k = %.2f, P(chi2 <= this) = %.3f, P(chi2 >= this) = %.3f, max |pull| = %.2f' % (
        label, x, k, x/k, chi2.cdf(x, k), chi2.sf(x, k), max(abs(p) for p in rows)))
    return x, k
fit([p for _, p, _ in C_ROWS], 'calculator, all committed numbers with a measurement')
fit([p for _, p, o in C_ROWS if o], "calculator, its own numbers")
# the standard theories' eleven predictions of independent measurements (Table 78; pulls as std_bits.py computes them)
S_ROWS = [0.0, (137.035999206 - 137.035999166)/math.hypot(15e-9, 11e-9), 2.1/2.3, 5.5/math.hypot(6, 7.7), (0.9983 - 1)/0.0016,
          (6601.8 - 6606.1)/18.3, (0.2449 - 0.24691)/math.hypot(0.00018, 0.0040), (2.55 - 2.57)/math.hypot(0.13, 0.03),
          (2.86 - 3.044)/0.13, (37.2 - 39.2)/7.2, (125.13 - 91)/30]
fit(S_ROWS, "standard theories, eleven predictions of independent measurements")
# what survival does: keep only alignments within 2 sigma of the data, from honest Gaussian errors
rng = np.random.default_rng(20261011)
z = rng.standard_normal((200000, len(C_ROWS))); keep = np.where(np.abs(z) <= 2, z, np.nan)
x_surv = np.nansum(keep**2, axis=1)[np.all(~np.isnan(keep), axis=1)]
print('survivors of a 2 sigma cut, %d numbers: mean chi2 = %.2f, median = %.2f (honest errors: mean %d)' % (len(C_ROWS), x_surv.mean(), np.median(x_surv), len(C_ROWS)))
x_c = sum(p*p for _, p, _ in C_ROWS)
print('P(chi2 <= %.2f) among survivors = %.3f' % (x_c, np.mean(x_surv <= x_c)))

print('\n=== C. Downstream of the calculator\'s own numbers')
D21, D31 = 7.49e-5, 2.513e-3       # eV^2, NuFIT 6.1, normal ordering
def masses(m1): return m1, math.sqrt(m1*m1 + D21), math.sqrt(m1*m1 + D31)
def Q_flip(m1):
    a, b, c = masses(m1)
    return (a + b + c)/(-math.sqrt(a) + math.sqrt(b) + math.sqrt(c))**2 - 2/3
m1 = brentq(Q_flip, 1e-7, 2e-3, xtol=1e-16)
ma, mb, mc = masses(m1)
S = ma + mb + mc
s12, s13 = (1 - s15)/(3*(1 - s15/3)), s15/3
Ue = [(1 - s12)*(1 - s13), s12*(1 - s13), s13]
mbeta = math.sqrt(sum(u*m*m for u, m in zip(Ue, (ma, mb, mc))))
t = [u*m for u, m in zip(Ue, (ma, mb, mc))]
mbb_max = sum(t); mbb_min = max(0.0, max(t) - (sum(t) - max(t)))
print('neutrino masses %.4f, %.4f, %.4f meV; sum %.3f meV; floor %.3f meV' % (ma*1e3, mb*1e3, mc*1e3, S*1e3, (math.sqrt(D21) + math.sqrt(D31))*1e3))
print('|U_e1|^2 = %.6f (2/3 = %.6f); m_beta = %.3f meV; m_bb from %.2f to %.2f meV' % (Ue[0], 2/3, mbeta*1e3, mbb_min*1e3, mbb_max*1e3))
print('sum against the 64 meV cap: %.1f meV under it; against the LCDM effective mass -0.101 +0.047: %+.2f sigma' % ((0.064 - S)*1e3, (-0.101 - S)/0.047))
print('KATRIN cap 0.45 eV is %.0f times m_beta; KamLAND-Zen floor 28 meV is %.1f times the largest m_bb' % (0.45/mbeta, 0.028/mbb_max))
# cosmology
Om = math.pi**2/32; OL = 1 - Om
h_cmb = math.sqrt(0.1430/Om); H_cmb = 100*h_cmb; sH_cmb = H_cmb*0.5*0.0011/0.1430
G = 6.67430e-11; C = 299792458.0; hbarc = 197.3269804e-9; eV = 1.602176634e-19; Mpc = 3.0856775814913673e22
H_hum = lambda meV: math.sqrt(8*math.pi*G*((meV*1e-3)**4/hbarc**3*eV)/(3*C**2*OL))*Mpc/1e3
Hh, Hhp, Hhm = H_hum(2.02), H_hum(2.22) - H_hum(2.02), H_hum(2.02) - H_hum(1.82)
print('H0 from the share and the CMB omega_m h^2 = 0.1430 +- 0.0011: %.2f +- %.2f km/s/Mpc' % (H_cmb, sH_cmb))
print('H0 from the share and the hum: %.2f +%.2f/-%.2f km/s/Mpc' % (Hh, Hhp, Hhm))
print('the two routes differ by %.2f sigma' % ((H_cmb - Hh)/math.hypot(Hhp, sH_cmb)))
for lab, v, s in [('CMB, Planck 2018', 67.4, 0.5), ('CMB, SPT+Planck+ACT', 67.19, 0.38), ('ladder, Riess 2022', 73.04, 1.04),
                  ('ladder, Breuval 2024', 73.17, 0.86), ('ladder, H0DN 2025', 73.50, 0.81), ('JWST ladders, Freedman 2025', 70.4, 1.9)]:
    print('   %-28s share+CMB route %+6.3f sigma; share+hum route %+6.3f sigma' % (
        lab, (v - H_cmb)/math.hypot(s, sH_cmb), (v - Hh)/math.hypot(s, Hhp if v > Hh else Hhm)))
Ez = math.sqrt(Om*1.349**3 + OL)
print('E(z = 0.349) = %.4f against 1.247 +- 0.006: %+.2f sigma; Planck LCDM curve 1.208: %+.2f sigma' % (Ez, (1.247 - Ez)/0.006, (1.247 - 1.208)/0.006))
H0t0 = 2/(3*math.sqrt(OL))*math.asinh(math.sqrt(OL/Om))
print('H0 t0 = %.4f; t0 = %.2f Gyr at the CMB route, %.1f Gyr at the hum route' % (H0t0, H0t0*977.792/H_cmb, H0t0*977.792/Hh))
print('matter = store at z = %.3f; acceleration began at z = %.3f; q0 = %.5f' % ((OL/Om)**(1/3) - 1, (2*OL/Om)**(1/3) - 1, 1.5*Om - 1))
for lab, v, s in [('DES-SN5YR alone', 0.648, 0.017), ('Pantheon+ alone', 0.666, 0.018), ('Union3 alone', 0.644, 0.026)]:
    print('   Omega_Lambda %-18s %+5.2f sigma' % (lab, (v - OL)/s))

print('\n=== D. The calendar of the frozen predictions')
PRED = [  # Table 81: name, value, sigma_c, today, below, above, future sigma, year, group, sigma_c at the decider's input precision
 ('sin2 theta12', (1 - s15)/(3*(1 - s15/3)), 0.0, 0.3096, 0.0073, 0.0057, 0.0016, 2031, 'own', 0.0),
 ('sin2 theta13', s15/3, 0.0, 0.0220, 0.00056, 0.00056, 0.0004, 2030, 'own', 0.0),
 ('CKM alpha', 90.0, 0.0, 84.1, 3.0, 3.7, 2.0, 2030, 'own', 0.0),
 ('Omega_m', math.pi**2/32, 0.0, 0.3045, 0.0032, 0.0032, 0.003, 2030, 'own', 0.0),
 ('Delta G', 1/6, 0.0, 0.158, 0.039, 0.039, 0.016, 2035, 'own', 0.0),
 ('Delta Sigma', 1/3, 0.0, 0.352, 0.029, 0.029, 0.050, 2035, 'own', 0.0),
 ('top pole', 172.57, 0.38, 172.1, 0.6, 0.6, 0.035, 2045, 'own', 0.15),
 ('r_p', 4*0.2103089, 0.0, 0.84075, 0.00064, 0.00064, 0.0002, 2030, 'prior', 0.0),
 ('m_tau', 1776.969, 0.0, 1776.93, 0.09, 0.09, 0.05, 2035, 'prior', 0.0),
]
def log_p0(m, mu, sl, sh, sf):
    """today's split normal convolved with the future measurement's error, in closed form"""
    A = 2/(math.sqrt(2*math.pi)*(sl + sh)); out = np.zeros_like(m)
    for s, upper in ((sl, False), (sh, True)):
        v = s*s + sf*sf; sb = s*sf/math.sqrt(v); xb = (mu*sf*sf + m*s*s)/v
        arg = (xb - mu)/sb if upper else (mu - xb)/sb
        out = out + A*(sb/sf)*np.exp(-(m - mu)**2/(2*v))*norm.cdf(arg)
    return np.log(np.maximum(out, 1e-300))
NS = 200000
L = math.log(100)
rng = np.random.default_rng(20261011)
lnB_right, lnB_wrong, ref_wrong, ref_right = {}, {}, {}, {}
print('%-14s %10s %9s %9s %10s %10s' % ('prediction', 'value', 'pull now', 'pull then', 'refuted?', 'year'))
for name, c, sc, mu, sl, sh, sf, yr, grp, scp in PRED:
    s_c = math.hypot(sc, sf)
    mr = rng.normal(c, s_c, NS)
    side = rng.random(NS) < sh/(sl + sh); zz = np.abs(rng.standard_normal(NS))
    mw = np.where(side, mu + zz*sh, mu - zz*sl) + rng.normal(0, sf, NS)
    lr = -0.5*((mr - c)/s_c)**2 - math.log(s_c) - 0.5*math.log(2*math.pi) - log_p0(mr, mu, sl, sh, sf)
    lw = -0.5*((mw - c)/s_c)**2 - math.log(s_c) - 0.5*math.log(2*math.pi) - log_p0(mw, mu, sl, sh, sf)
    lnB_right[name], lnB_wrong[name] = lr, lw
    y1, y2 = s_c, math.hypot(scp, sf)
    ref_wrong[name] = (np.abs(mw - c)/y1 > 5) & (np.abs(mw - c)/y2 > 5)
    ref_right[name] = (np.abs(mr - c)/y1 > 5) & (np.abs(mr - c)/y2 > 5)
    side_now = sh if c > mu else sl
    p_now = (mu - c)/math.hypot(side_now, sc); p_then = (mu - c)/math.hypot(sc, sf); p_then2 = (mu - c)/math.hypot(scp, sf)
    print('%-14s %10.6g %+9.2f %+9.2f %10s %6d  median B if right %.3g, if wrong %.3g' % (
        name, c, p_now, p_then, 'yes' if (abs(p_then) > 5 and abs(p_then2) > 5) else ('at decider' if abs(p_then2) > 5 else 'no'), yr,
        math.exp(np.median(lr)), math.exp(np.median(lw))))
own = [p for p in PRED if p[8] == 'own']
print('\nthe calculator\'s own seven, cumulative by year (frozen rule: decisive above 100, against below 1/100)')
print('%6s %4s %14s %12s %14s %12s %14s %14s' % ('by', 'n', 'median B|right', 'P(>100|right)', 'median B|wrong', 'P(<1/100|wrong)', 'P(refuted|wrong)', 'P(refuted|right)'))
for Y in (2030, 2031, 2035, 2045):
    sel = [p[0] for p in own if p[7] <= Y]
    R = sum(lnB_right[n] for n in sel); W = sum(lnB_wrong[n] for n in sel)
    rw = np.zeros(NS, bool); rr = np.zeros(NS, bool)
    for n in sel: rw |= ref_wrong[n]; rr |= ref_right[n]
    print('%6d %4d %14.4g %12.4f %14.4g %12.4f %14.4f %14.5f' % (Y, len(sel), math.exp(np.median(R)), np.mean(R > L), math.exp(np.median(W)), np.mean(W < -L), np.mean(rw), np.mean(rr)))
R = sum(lnB_right[p[0]] for p in own)
print('if right, all seven: log2 B median %.1f bits, 10th to 90th percentile %.1f to %.1f bits' % (
    np.median(R)/math.log(2), np.percentile(R, 10)/math.log(2), np.percentile(R, 90)/math.log(2)))
Rp = sum(lnB_right[p[0]] for p in PRED if p[8] == 'prior')
print('prior-art pair if right: median B %.3g' % math.exp(np.median(Rp)))
print('for scale: the standard theories\' eleven predictions carry 102 bits (Table 78)')

print('\n=== E. The verdict: the author\'s four cores at three levels')
NR, RJ, UN = 'NOT_REJECT', 'REJECT', 'UNDETERMINED'
def synth(cores):
    if RJ in cores: return 'REJECTED'
    if UN in cores: return 'UNSTABLE'
    return 'STABILIZED'
CLAIMS = [  # cores in the order outcome, constraint, alternative, consequence (Table 98)
 ('calculator', '1. its laws hold where tested', (NR, NR, NR, NR)),
 ('calculator', '2a. its own numbers already beat chance', (NR, RJ, RJ, UN)),
 ('calculator', '2b. its frozen predictions will hold', (UN, NR, UN, NR)),
 ('calculator', '3. it accounts for every observed behaviour', (RJ, RJ, UN, RJ)),
 ('standard', "1'. their laws hold where tested", (NR, NR, NR, NR)),
 ('standard', "2'. their predictions beat chance", (NR, NR, NR, NR)),
 ('standard', "3'. they account for every observed behaviour", (RJ, RJ, UN, RJ)),
]
for th, cl, cores in CLAIMS: print('   %-10s %-46s %s -> %s' % (th, cl, ' '.join(cores), synth(cores)))
