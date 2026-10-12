"""Appendix G: top-mass measurements against the rule and the edge; open-data pulls; ADM 2012 check."""
import sys, math
sys.path.insert(0,'/tmp/claude-0/-home-claude/bbcd54b6-3962-59fe-98a5-d2cef4b7e16c/scratchpad')
from higgs2 import summary
from scipy.optimize import brentq
rule_t, s_rule_t = 172.567, 0.379
edge_t, s_edge_t = 170.935, 0.375
rule_h, s_rule_h = 125.196, 0.924
edge_h, s_edge_h = 128.548, 0.964
def pull(m, sm, p, sp): return (m-p)/math.hypot(sm, sp)
def apull(m, up, dn, p, sp):
    # asymmetric measurement: use the side facing the prediction
    s = dn if m > p else up
    return (m-p)/math.hypot(s, sp)
print("open data Higgs 124.20 +- 0.62:")
for lab, p, sp in (("PDG 2026", 125.13, 0.11), ("rule", rule_h, s_rule_h), ("edge", edge_h, s_edge_h)):
    print(f"  vs {lab}: {pull(124.20, 0.62, p, sp):+.2f}")
print("open data top 172.54 +- 0.94:")
for lab, p, sp in (("direct 172.60", 172.60, 0.27), ("pole 172.1", 172.1, 0.6), ("ATLAS 2025 pole 170.7", 170.7, 1.5), ("rule", rule_t, s_rule_t), ("edge", edge_t, s_edge_t)):
    print(f"  vs {lab}: {pull(172.54, 0.94, p, sp):+.2f}")
rows = [("ATLAS Run 1 combination, direct", 172.69, 0.48, 0.48),
        ("ATLAS + CMS Run 1 combination, direct", 172.52, 0.33, 0.33),
        ("Run 2 soft muon 36 fb-1, direct", 174.41, 0.81, 0.81),
        ("Run 2 dilepton 140 fb-1 (prelim), direct", 172.21, 0.80, 0.80),
        ("cross section 36 fb-1, pole", 173.1, 2.0, 2.1),
        ("cross section 140 fb-1 (worked), pole", 172.72, 1.46, 1.73),
        ("ATLAS 2025 tt+1 jet, pole", 170.7, 1.5, 1.5),
        ("PDG 2026 direct", 172.60, 0.27, 0.27),
        ("PDG 2026 pole", 172.1, 0.6, 0.6)]
for lab, m, up, dn in rows:
    print(f"{lab:44s} {m:7.2f}  rule {apull(m, up, dn, rule_t, s_rule_t):+.2f}  edge {apull(m, up, dn, edge_t, s_edge_t):+.2f}")
# worked pole from cross section
th, scale_up, scale_dn, pdf = 834, 21, 30, 21
meas, st, sy, lu = 829.0, 1.0, 13.2, 8.0
slope = 23
m = 172.5+(th-meas)/slope
e_exp = math.sqrt(st**2+sy**2+lu**2)/slope
e_th_up = math.hypot(scale_up, pdf)/slope; e_th_dn = math.hypot(scale_dn, pdf)/slope
print(f"worked pole: {m:.3f} +{math.hypot(e_exp, e_th_up):.2f} -{math.hypot(e_exp, e_th_dn):.2f} (exp {e_exp:.2f}, th +{e_th_up:.2f}/-{e_th_dn:.2f})")
mh = brentq(lambda M: summary(Mt=173.3, Mh=M, als=0.1180, MW=80.3625)['lmin'], 120, 140, xtol=1e-5)
print(f"stability bound at m_t=173.3 (PDG 2026 als, MW): {mh:.2f} (ADM 2012: 129.8 +- 5.6)")
# SM predictions vs PDG 2026
print(f"W: {pull(80.3625, 0.0077, 80.356, 0.005):+.2f}; top fit: {pull(172.60, 0.27, 172.85, 0.55):+.2f}; Higgs 2011: {(125.13-91)/30:+.2f}")
print(f"double root top pull {pull(172.60,0.27,174.289,0.57):+.2f}; Higgs {pull(125.13,0.11,132.133,0.47):+.2f}; FN 1995 top {pull(172.60,0.27,173,5):+.2f} Higgs {pull(125.13,0.11,135,9):+.2f}")
print(f"weak field alone on the Higgs: {128.548/1.00531:.2f}, pull {pull(125.13,0.11,128.548/1.00531,0.964/1.00531):+.2f}")
print(f"(m_t/Mbar)^2 = {(172.60/2.435e18)**2:.2e}; (m_H/Mbar)^2 = {(125.13/2.435e18)**2:.2e}")
print(f"Higgs depth m_H^2 v^2/8 = {125.13**2*246.21965**2/8:.4e} GeV^4")
xi = 4.5+(172.60-173.3)*((11.5-4.5)/(172.1-173.3)); print(f"xi upper bound interpolated at 172.60: {xi:.2f}")
