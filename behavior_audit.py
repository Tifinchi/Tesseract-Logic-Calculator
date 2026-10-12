# The behavior rule across the frozen set (Appendix M.4): down-wind observables for each prediction.
# Inputs are the values the paper already cites; each line names its source.
from math import sqrt, pi
def pull(meas, sm, calc, sc=0.0): return (meas - calc)/sqrt(sm**2 + sc**2)
print("3. alpha = 90: 180 - beta - gamma = 91.0 +- 2.8 (HFLAV gamma): %+.2f; 89.4 +- 4.0 (PDG gamma): %+.2f" % (pull(91.0, 2.8, 90), pull(89.4, 4.0, 90)))
om = pi**2/32; omh2, s_omh2 = 0.1430, 0.0011               # Planck 2018 (Planck Collaboration 2020), TT,TE,EE+lowE+lensing
h = sqrt(omh2/om); sh = 0.5*h*s_omh2/omh2
H0, sH0 = 100*h, 100*sh
print("4. Omega_m = pi^2/32 = %.6f with the CMB's Omega_m h^2 %.4f +- %.4f: H0 = %.2f +- %.2f" % (om, omh2, s_omh2, H0, sH0))
for name, v, s in (("SH0ES 73.04 +- 1.04 (Riess et al. 2022)", 73.04, 1.04), ("H0DN 73.50 +- 0.81 (H0DN 2025)", 73.50, 0.81),
                   ("CCHP 70.39 +- 1.94 (Freedman et al. 2024)", 70.39, 1.94), ("ACT with DESI 68.43 +- 0.27", 68.43, 0.27)):
    print("   against %-42s %+.1f sigma" % (name, pull(v, s, H0, sH0)))
print("   Planck's own H0, 67.4 +- 0.5, against SH0ES: %+.1f sigma (the tension LCDM carries)" % pull(73.04, 1.04, 67.4, 0.5))
print("   the share 1 - pi^2/32 = %.4f against Omega_Lambda 0.6847 +- 0.0073: %+.1f; DESI BAO 0.2975 +- 0.0086 against pi^2/32: %+.1f" % (1 - om, pull(0.6847, 0.0073, 1 - om), pull(0.2975, 0.0086, om)))
print("6. Delta Sigma = 1/3 with Delta G = 1/6 at 2 GeV gives L = 1/6; folded lattice L = 0.152 +- 0.043: %+.2f; HERMES 0.330 +- 0.039: %+.2f; lattice 0.380 +- 0.044: %+.2f" % (
    pull(0.152, 0.043, 1/6), pull(0.330, 0.039, 1/3), pull(0.380, 0.044, 1/3)))
floor = sqrt(7.49e-5) + sqrt(2.513e-3)
print("recorded: normal-ordering floor = %.1f meV; the calculator's 59.2 meV is %.1f meV above it" % (1000*floor, 59.2 - 1000*floor))
print("recorded: vacuum energy 2.01 +- 0.20 against 2.24 meV: %+.2f" % pull(2.24, 0.0, 2.01, 0.20))
