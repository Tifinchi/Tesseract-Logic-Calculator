# Misfits located by their neighbours (section 4.6): the arithmetic.
import math
from mpmath import mp, mpf, pi, exp, mpc
mp.dps = 30
print("=== Electron: one circuit read as the half turn")
# a state turned through 2*pi picks up exp(-i 2 pi J / hbar); landing in the inverse cell (-1) means J = hbar/2
for J in (mpf(1)/2, mpf(1)):
    ph = exp(mpc(0, -2*pi*J))
    print("  J = %s hbar: phase after one circuit = %s" % (J, mp.nstr(ph, 6)))
# g = mu / (J e / 2m) with mu = mu_B = e hbar / 2m and J = hbar/2
print("  g = mu_B / ((hbar/2)(e/2m)) = 2;  with J = hbar: g = 1")
alpha = 1/mpf('137.035999206')
print("  alpha/2pi = %s; measured a_e = 0.00115965218059 (Table 5)" % mp.nstr(alpha/(2*pi), 8))

print("=== Matter share against pi^2/32")
om = math.pi**2/32
for name, v, s in (("DES-SN5YR", 0.352, 0.017), ("DES-Dovekie", 0.330, 0.015), ("Pantheon+", 0.334, 0.018), ("Union3 (carried)", 0.356, 0.027),
                   ("DESI BAO", 0.2975, 0.0086), ("DES-Dovekie + CMB", 0.317, 0.005), ("SN + CMB + BAO, flat LCDM", 0.3045, 0.0032)):
    print("  %-28s %.4f +- %.4f: %+.2f sigma" % (name, v, s, (v-om)/s))

print("=== Acceleration at the cosmos's inverse corner")
q_cube = 1.5*om - 1
print("  cube q0 = 3*Om/2 - 1 = %.5f" % q_cube)
for name, v, s in (("DES-Dovekie+BAO+CMB, flat wCDM", -0.537, 0.026), ("Pantheon+ below z=0.8 (SH0ES)", -0.51, 0.024), ("DES-Dovekie SN only, flat wCDM", -0.423, 0.073)):
    print("  %-34s %.3f +- %.3f: %+.2f sigma" % (name, v, s, (v-q_cube)/s))
Om, w = 0.305, -0.995
print("  q0 implied by Om = 0.305, w = -0.995: %.4f (q0 = Om/2 + (1+3w)(1-Om)/2)" % (Om/2 + (1+3*w)*(1-Om)/2))
Om, w0 = 0.313, -0.803
print("  q0 implied by the w0wa best fit (Om = 0.313, w0 = -0.803): %.4f" % (Om/2 + (1+3*w0)*(1-Om)/2))
print("  constant w = -0.995 +0.019/-0.020: %.2f sigma from -1" % ((-0.995+1)/0.019))

print("=== theta13, the neighbour of theta12 in the 15 degree frame")
print("  0.02200 +- 0.00056 against 0.02233: %+.2f sigma" % ((0.02200-0.02233)/0.00056))
print("=== Curvature: SN+CMB against BAO+CMB")
print("  SN+CMB -0.0063 +- 0.0037: %+.2f sigma; BAO+CMB +0.0025 +- 0.0011: %+.2f sigma; SN+BAO+CMB +0.0026 +- 0.0011: %+.2f sigma" % (-0.0063/0.0037, 0.0025/0.0011, 0.0026/0.0011))
