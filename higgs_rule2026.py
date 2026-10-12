"""Balance rule with the three-phase load, PDG 2026 inputs.
Rule: M_t = M_t,edge(M_H) * (1 + delta3), delta3 = |y_t^2 + g^2 e^{i120}| / (8 pi^2), couplings at mu = M_t.
Edge: min over mu of lambda(mu) = 0 (two-loop running, three-loop QCD; higgs2.py)."""
import sys, math, time
sys.path.insert(0,'/tmp/claude-0/-home-claude/bbcd54b6-3962-59fe-98a5-d2cef4b7e16c/scratchpad')
from higgs2 import summary, inputs
from scipy.optimize import brentq
MH, dMH = 125.13, 0.11
MT, dMT = 172.60, 0.27
ALS, dALS = 0.1180, 0.0009
MW = 80.3625
def lmin(Mt, Mh, als=ALS): return summary(Mt=Mt, Mh=Mh, als=als, MW=MW)['lmin']
def mt_edge(Mh, als=ALS): return brentq(lambda M: lmin(M, Mh, als), 165, 178, xtol=1e-6)
def mh_edge(Mt, als=ALS): return brentq(lambda M: lmin(Mt, M, als), 118, 140, xtol=1e-6)
def delta3(Mt, als=ALS):
    g1, g2, g3, yt, l = inputs(Mt=Mt, Mh=MH, als=als, MW=MW)
    z = yt**2 + g2**2*complex(math.cos(2*math.pi/3), math.sin(2*math.pi/3))
    return abs(z)/(8*math.pi**2), yt, g2
t0 = time.time()
d, yt, g2 = delta3(MT)
print(f"y_t(M_t)={yt:.4f} g(M_t)={g2:.4f} y_t^2={yt**2:.4f} g^2={g2**2:.4f} delta3={100*d:.4f}%")
e = mt_edge(MH); print(f"top edge at M_H={MH}: {e:.4f}")
top = e*(1+d); print(f"rule top: {top:.4f}")
eH = mh_edge(MT); print(f"Higgs edge at M_t={MT}: {eH:.4f}")
# Higgs predicted from the measured top: M_t,edge(M_H) = M_t/(1+delta3)
d_meas = delta3(MT)[0]
mh_pred = mh_edge(MT/(1+d_meas)); print(f"rule Higgs from measured top: {mh_pred:.4f}  (edge top used {MT/(1+d_meas):.4f})")
# sensitivities
h = 0.5
s_top_mh = (mt_edge(MH+h)-mt_edge(MH-h))/(2*h); print(f"dMt_edge/dMh at {MH}: {s_top_mh:.4f} (inverse {1/s_top_mh:.4f})")
s_top_als = (mt_edge(MH, ALS+0.001)-mt_edge(MH, ALS-0.001))/0.002; print(f"dMt_edge/dals: {s_top_als:.2f} GeV per unit -> {s_top_als*dALS:.4f} GeV for dals {dALS}")
Mte = MT/(1+d_meas)
s_mh_mt = (mh_edge(Mte+h)-mh_edge(Mte-h))/(2*h); print(f"dMh_edge/dMt at {Mte:.3f}: {s_mh_mt:.4f}")
s_mh_als = (mh_edge(Mte, ALS+0.001)-mh_edge(Mte, ALS-0.001))/0.002; print(f"dMh_edge/dals at Mte: {s_mh_als:.2f} -> {s_mh_als*dALS:.4f} GeV")
# errors
err_top_mh = s_top_mh*dMH*(1+d); err_top_als = s_top_als*dALS*(1+d)
print(f"rule top error parts: from M_H {err_top_mh:.4f}, from als {err_top_als:.4f}; quad {math.hypot(err_top_mh, err_top_als):.4f}")
err_h_mt = s_mh_mt*dMT/(1+d_meas); err_h_als = s_mh_als*dALS
print(f"rule Higgs error parts: from M_t {err_h_mt:.4f}, from als {err_h_als:.4f}; quad {math.hypot(err_h_mt, err_h_als):.4f}")
pull_top = (MT-top)/math.hypot(dMT, math.hypot(err_top_mh, err_top_als))
print(f"pull top (meas - rule)/sigma = {pull_top:+.3f}")
for tot in (math.hypot(err_h_mt, err_h_als), 0.92):
    print(f"pull Higgs with rule sigma {tot:.3f}: {(MH-mh_pred)/math.hypot(dMH, tot):+.3f}")
# exact-edge rule pulls
eH_err_mt = (mh_edge(MT+h)-mh_edge(MT-h))/(2*h)*dMT; eH_err_als = (mh_edge(MT, ALS+0.001)-mh_edge(MT, ALS-0.001))/0.002*dALS
eHs = math.hypot(eH_err_mt, eH_err_als)
print(f"exact-edge Higgs {eH:.3f} +- {eHs:.3f} (Mt part {eH_err_mt:.3f}, als part {eH_err_als:.3f}); pull {(MH-eH)/math.hypot(dMH, eHs):+.2f}")
e_err = math.hypot(s_top_mh*dMH, s_top_als*dALS)
print(f"exact-edge top {e:.3f} +- {e_err:.3f}; pull vs direct {(MT-e)/math.hypot(dMT, e_err):+.2f}; vs pole 172.1+-0.6 {(172.1-e)/math.hypot(0.6, e_err):+.2f}")
print(f"rule top vs pole 172.1+-0.6: {(172.1-top)/math.hypot(0.6, math.hypot(err_top_mh, err_top_als)):+.2f}")
print(f"{time.time()-t0:.1f} s")
# theory errors of the matching (Buttazzo et al. 2013): lambda(M_t) +-0.00030, y_t(M_t) +-0.00050
def mt_edge_th(Mh, dl=0.0, dy=0.0):
    return brentq(lambda M: summary(Mt=M, Mh=Mh, als=ALS, MW=MW, dl=dl, dy=dy)['lmin'], 165, 178, xtol=1e-6)
def mh_edge_th(Mt, dl=0.0, dy=0.0):
    return brentq(lambda M: summary(Mt=Mt, Mh=M, als=ALS, MW=MW, dl=dl, dy=dy)['lmin'], 118, 140, xtol=1e-6)
tl = (mt_edge_th(MH, dl=0.0003)-mt_edge_th(MH, dl=-0.0003))/2
ty = (mt_edge_th(MH, dy=0.0005)-mt_edge_th(MH, dy=-0.0005))/2
print(f"top edge theory: from lambda {tl:+.4f}, from y_t {ty:+.4f}, quad {math.hypot(tl, ty):.4f}")
top_err = math.sqrt((s_top_mh*dMH)**2+(s_top_als*dALS)**2+tl**2+ty**2)*(1+d)
print(f"rule top total error {top_err:.4f}")
hl = (mh_edge_th(Mte, dl=0.0003)-mh_edge_th(Mte, dl=-0.0003))/2
hy = (mh_edge_th(Mte, dy=0.0005)-mh_edge_th(Mte, dy=-0.0005))/2
print(f"Higgs theory: from lambda {hl:+.4f}, from y_t {hy:+.4f}")
h_err = math.sqrt(err_h_mt**2+err_h_als**2+hl**2+hy**2)
print(f"rule Higgs total error {h_err:.4f}; pull {(MH-mh_pred)/math.hypot(dMH, h_err):+.3f}")
eHl = (mh_edge_th(MT, dl=0.0003)-mh_edge_th(MT, dl=-0.0003))/2
eHy = (mh_edge_th(MT, dy=0.0005)-mh_edge_th(MT, dy=-0.0005))/2
eH_tot = math.sqrt(eH_err_mt**2+eH_err_als**2+eHl**2+eHy**2)
print(f"exact-edge Higgs total error {eH_tot:.4f}; pull {(MH-eH)/math.hypot(dMH, eH_tot):+.3f}")
e_tot = math.sqrt((s_top_mh*dMH)**2+(s_top_als*dALS)**2+tl**2+ty**2)
print(f"exact-edge top total error {e_tot:.4f}; pull vs direct {(MT-e)/math.hypot(dMT, e_tot):+.3f}")
print(f"rule top pull with total error: {(MT-top)/math.hypot(dMT, top_err):+.3f}")
