"""Phase angle read from the masses, and other load candidates, PDG 2026."""
import sys, math
sys.path.insert(0,'/tmp/claude-0/-home-claude/bbcd54b6-3962-59fe-98a5-d2cef4b7e16c/scratchpad')
from higgs2 import summary, inputs
from scipy.optimize import brentq
MH, MT, ALS, MW = 125.13, 172.60, 0.1180, 80.3625
g1, g2, g3, yt, l = inputs(Mt=MT, Mh=MH, als=ALS, MW=MW)
Y, G = yt**2, g2**2
E = 170.93468  # top edge at MH (higgs_rule2026.py: 170.9347)
E = brentq(lambda M: summary(Mt=M, Mh=MH, als=ALS, MW=MW)['lmin'], 165, 178, xtol=1e-7)
print(f"y_t^2={Y:.5f} g^2={G:.5f} edge={E:.5f}")
def angle(top):
    d = top/E-1; z = d*8*math.pi**2
    c = (z*z-Y*Y-G*G)/(2*Y*G)
    return math.degrees(math.acos(max(-1,min(1,c)))), d
sig = math.hypot(0.27, 0.3788/1.00955)  # measured top error and edge error (total incl. theory)
for lab, top in (("best", MT), ("+1sigma", MT+sig), ("-1sigma", MT-sig)):
    a, d = angle(top); print(f"{lab}: top {top:.3f} delta {100*d:.4f}% angle {a:.2f}")
# candidates: y_t^2 alone (perpendicular conductivity), and angles 0, 90, 120, 180
for lab, z in (("y_t^2 alone", Y), ("0 deg", Y+G), ("90 deg", abs(complex(Y, G))), ("120 deg", abs(Y+G*complex(-0.5, math.sqrt(3)/2))), ("180 deg", Y-G)):
    d = z/(8*math.pi**2); top = E*(1+d)
    Mte = MT/(1+d)
    mh = brentq(lambda M: summary(Mt=Mte, Mh=M, als=ALS, MW=MW)['lmin'], 115, 140, xtol=1e-6)
    print(f"{lab}: delta {100*d:.4f}% top {top:.3f} pull {(MT-top)/math.hypot(0.27,0.3788):+.2f} | Higgs from measured top {mh:.3f} pull {(MH-mh)/math.hypot(0.11,0.924):+.2f}")
