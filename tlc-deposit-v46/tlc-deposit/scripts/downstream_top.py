# The author's selection rule applied to the top: which mass the rule's 172.57 +- 0.38 GeV is aligned with,
# and what each alignment implies for the pole-type observables down-wind of it.
from math import sqrt
c = [0.4244, 0.8345, 2.375, 8.49]          # pole/MS-bar at four loops (Marquard et al. 2015), as in top_scheme.py
def ratio(a): return 1 + sum(ci*a**(i+1) for i, ci in enumerate(c))
rule, sr = 172.57, 0.38
obs = {"pole average (PDG 2026)": (172.1, 0.6), "electroweak fit, indirect (PDG 2026)": (175.2, 1.8)}
def show(label, pole, sp):
    print(label)
    print("   implied pole mass %.2f +- %.2f GeV; edge 170.94 -> %s" % (pole, sp, "above (metastable, the rule's side)" if pole > 170.94 else "at or below (stable: breaks the rule)"))
    for k, (v, s) in obs.items():
        print("   %-38s %+.1f sigma" % (k, (v - pole)/sqrt(s**2 + sp**2)))
show("pole = 172.57 (= MS-bar 162.76)", rule, sr)
r = ratio(0.1080)                           # alpha_s at about 172.6 GeV
show("MS-bar m_t(m_t) = 172.57, unconverted", rule*r, sr*r)
show("direct (simulation mass) = 172.57, pole = MC - 0.57 +- 0.28 (Butenschoen et al. 2016)", rule - 0.57, sqrt(sr**2 + 0.28**2))
show("direct = 172.57, MC - pole = 0 +- 0.52 (PDG electroweak review)", rule, sqrt(sr**2 + 0.52**2))
