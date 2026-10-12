# Evidence against chance for the standard theories' advance predictions, on the convention of section 16.8:
# a single number fixed in advance; chance p = width of the window (in ln) / ln 10, the window being
# +-max(|pull|,1) combined sigma around the prediction; bits = -log2 p.  Fisher's method combines.
from mpmath import mp, mpf, log, sqrt, log10
from scipy.stats import chi2
import math
mp.dps = 30
LN10 = math.log(10)
def window_p(pred, sig_lo, sig_hi, pull):
    k = max(abs(pull), 1.0)
    lo = pred - k*sig_lo; hi = pred + k*sig_hi
    return math.log(hi/lo)/LN10
rows = []
def add(name, pred, sig_lo, sig_hi, meas, note=''):
    sig = sig_hi if meas >= pred else sig_lo
    pull = (meas - pred)/sig
    p = window_p(pred, sig_lo, sig_hi, pull)
    rows.append((name, pred, meas, pull, p, -math.log2(p), note))
# light bending: gamma = 1 predicted, measured 1 + (2.1 +- 2.3)e-5
add('light bending, gamma', 1.0, 2.3e-5, 2.3e-5, 1+2.1e-5)
# binary pulsar: ratio to GR 1 predicted, measured 0.9983 +- 0.0016
add('binary pulsar decay / GR', 1.0, 0.0016, 0.0016, 0.9983)
# GW speed: (v-c)/c in [-3e-15, +7e-16]
p_gw = (7e-16 + 3e-15)/LN10
rows.append(('gravitational-wave speed', 1.0, 1.0, 0.0, p_gw, -math.log2(p_gw), 'window = the measured interval'))
# geodetic: GR 6606.1, measured 6601.8 +- 18.3
add('geodetic precession', 6606.1, 18.3, 18.3, 6601.8)
# frame dragging: GR 39.2, measured 37.2 +- 7.2
add('frame dragging', 39.2, 7.2, 7.2, 37.2)
# W mass: SM 80357 +- 6, measured 80362.5 +- 7.7
s = math.hypot(6, 7.7)
add('W mass', 80357, s, s, 80362.5)
# Higgs from EW fit 2011: 91 +30/-23, measured 125.13 +- 0.11
add('Higgs mass, electroweak fit 2011', 91, 23, 30, 125.13)
# electron g-2: alpha^-1 from g-2 137.035999166(15) vs Rb 137.035999206(11)
s = math.hypot(15e-9, 11e-9)
add('electron g-2 against rubidium alpha', 137.035999166, s, s, 137.035999206)
s2 = math.hypot(15e-9, 27e-9)
add('electron g-2 against caesium alpha', 137.035999166, s2, s2, 137.035999046, 'not counted; Rb and Cs disagree')
# N_eff: 3.044 vs 2.86 +- 0.13
add('N_eff', 3.044, 0.13, 0.13, 2.86)
# D/H: predicted 2.57 +- 0.13, observed 2.55 +- 0.03
s = math.hypot(0.13, 0.03)
add('D/H', 2.57, s, s, 2.55)
# Yp: predicted 0.24691 +- 0.00018, observed 0.2449 +- 0.0040
s = math.hypot(0.00018, 0.0040)
add('helium-4 Yp', 0.24691, s, s, 0.2449)
for r in rows:
    print(f"{r[0]:40s} pull {r[3]:+.2f}  p {r[4]:.3g}  bits {r[5]:.1f} {r[6]}")
counted = [r for r in rows if 'not counted' not in r[6]]
X = -2*sum(math.log(r[4]) for r in counted)
k = len(counted)
pc = chi2.sf(X, 2*k)
print('counted', k, 'sum of bits', round(sum(r[5] for r in counted),1), 'Fisher X', round(X,1), 'combined p', pc, 'bits', -math.log2(pc))
# robustness: narrower prior, +-10% (ln range ln(1.1/0.9))
narrow = math.log(1.1/0.9)
print('per-entry change with a +-10% prior:', -math.log2(LN10/narrow), 'bits')
# the calculator, from Table 81 of the earlier paper: combined strict 0.064 (4.0 bits), generous 0.68 (0.6 bits); inherited 5.3e-5 (14.2 bits)
for p in (0.064, 0.68, 5.3e-5):
    print(p, -math.log2(p))
# misses
print('lithium pull', (1.6-4.72)/math.hypot(0.72,0.3))
print('H0 Planck vs Breuval', (73.17-67.4)/math.hypot(0.86,0.5))
