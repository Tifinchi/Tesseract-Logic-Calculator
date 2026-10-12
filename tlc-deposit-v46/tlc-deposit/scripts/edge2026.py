import sys, math, time
sys.path.insert(0,'/tmp/claude-0/-home-claude/bbcd54b6-3962-59fe-98a5-d2cef4b7e16c/scratchpad')
from higgs2 import *
from scipy.optimize import brentq
t0=time.time()
C=dict(Mh=125.13,Mt=172.60,als=0.1180,MW=80.3625)
r=summary(**C)
print("PDG2026: lambda_min",r['lmin'],"at",r['mumin'],"cross",r['cross'],"lPl",r['lPl'])
mt_c=brentq(lambda M:summary(Mt=M,Mh=125.13,als=0.1180,MW=80.3625)['lmin'],165,175,xtol=1e-5)
mh_c=brentq(lambda M:summary(Mh=M,Mt=172.60,als=0.1180,MW=80.3625)['lmin'],120,140,xtol=1e-5)
print(f"edge top (Mh=125.13): {mt_c:.4f}  edge Higgs (Mt=172.60): {mh_c:.4f}")
for lab,mt in (("delta3",None),):
    pass
yt2,g2=0.8705,0.4194
for ang in (0,90,120,180):
    d=abs(yt2+g2*complex(math.cos(math.radians(ang)),math.sin(math.radians(ang))))/(8*math.pi**2)
    top=mt_c*(1+d)
    print(ang,f"delta={d*100:.4f}% top={top:.3f} pull={(172.60-top)/math.hypot(0.27,0.38):+.2f}")
# slope dMh/dMt along the edge
mh1=brentq(lambda M:summary(Mh=M,Mt=172.60+0.5,als=0.1180,MW=80.3625)['lmin'],120,140,xtol=1e-5)
print("dMh_edge/dMt",(mh1-mh_c)/0.5)
print(time.time()-t0,"s")
