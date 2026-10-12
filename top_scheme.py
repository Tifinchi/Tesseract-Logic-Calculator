# Top-mass determinations by scheme against the balance rule (172.57 +- 0.38 GeV, pole) and the exact edge
# (170.935 +- 0.375 GeV, pole), and the rule converted to the MS-bar mass m_t(m_t).
from math import sqrt
rule, srule = 172.57, 0.38
edge, sedge = 170.935, 0.375
def pull(v, up, dn, c, sc):
    s = dn if v > c else up          # the measurement's error on the side facing the prediction
    return (v - c)/sqrt(s*s + sc*sc)
rows = [
 ("PDG 2026 average of pole determinations", "pole", 172.1, 0.6, 0.6),
 ("PDG 2026 direct average, MC-pole difference 0 +- 0.52 (PDG EW review)", "pole-equivalent", 172.60, sqrt(0.27**2+0.52**2), sqrt(0.27**2+0.52**2)),
 ("direct average with the Butenschoen et al. N2LL calibration (pole = MC - 0.57 +- 0.28)", "pole-equivalent", 172.60-0.57, sqrt(0.27**2+0.28**2), sqrt(0.27**2+0.28**2)),
 ("CMS normalised multi-differential cross sections, 13 TeV", "pole", 170.5, 0.8, 0.8),
 ("CMS tt+jet, 13 TeV", "pole", 172.93, 1.36, 1.36),
 ("ATLAS tt+1 jet, 13 TeV", "pole", 170.7, 1.5, 1.5),
 ("ATLAS+CMS inclusive cross sections, 7 and 8 TeV", "pole", 173.4, 1.8, 2.0),
 ("ATLAS inclusive cross section, 36 fb-1", "pole", 173.1, 2.0, 2.1),
 ("ATLAS Run 2 cross section, 140 fb-1, worked out (G.6)", "pole", 172.72, 1.46, 1.73),
 ("electroweak fit, indirect (PDG 2026)", "pole", 175.2, 1.8, 1.8),
]
for name, kind, v, up, dn in rows:
    print("%-88s %-15s %7.2f +%.2f/-%.2f  rule %+5.1f sigma  edge %+5.1f sigma" % (name, kind, v, up, dn, pull(v, up, dn, rule, srule), pull(v, up, dn, edge, sedge)))
# pole -> MS-bar, 4-loop relation (Marquard et al. 2015), nl = 5 light flavours, alpha_s^(6)(m_t) as in PDG 2026 (0.1088)
c = [0.4244, 0.8345, 2.375, 8.49]
def ratio(a): return 1 + sum(ci*a**(i+1) for i, ci in enumerate(c))
def dratio(a): return sum((i+1)*ci*a**i for i, ci in enumerate(c))
a = 0.1088
print("\nPDG check: pole 172.5 -> m_t(m_t) = %.2f GeV (PDG 2026 quotes 162.69)" % (172.5/ratio(a)))
for name, m, s in (("rule", rule, srule), ("edge", edge, sedge)):
    mbar = m/ratio(a)
    sa = 0.0009*(a/0.1180)**2                       # alpha_s(m_t) error from alpha_s(mZ) +- 0.0009
    s_conv = m*dratio(a)*sa/ratio(a)**2
    print("%s: pole %.3f +- %.3f -> m_t(m_t) = %.2f +- %.2f GeV (conversion alpha_s part %.2f; the 4-loop term is %.2f GeV)" % (
        name, m, s, mbar, sqrt((s/ratio(a))**2 + s_conv**2), s_conv, m*c[3]*a**4/ratio(a)))
print("\nFCC-ee threshold scan (Defranchis et al. 2025): PS mass, stat 4.3, exp 6.8, theory ~35 MeV; MS-bar conversion 7-23 MeV;")
print("pole mass limited by its renormalon ambiguity, ~110 MeV (Beneke et al. 2017); the rule's own error today is %.0f MeV" % (srule*1000))
