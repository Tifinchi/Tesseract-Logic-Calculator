# Appendix B checks: the framework's V-I-R cube engine.
import itertools, math
import numpy as np
# 24 rotations as signed permutation matrices with det +1
rots=[]
for p in itertools.permutations(range(3)):
    for s in itertools.product([1,-1],repeat=3):
        M=np.zeros((3,3),int)
        for i in range(3): M[i,p[i]]=s[i]
        if round(np.linalg.det(M))==1: rots.append((p,s,M))
print("rotations",len(rots))
var="VIR"
# which variable lands on axis Z (row 2 picks column p[2])
from collections import Counter
byZ=Counter(var[p[2]] for p,s,M in rots)
byorder=Counter("".join(var[p[i]] for i in range(3)) for p,s,M in rots)
print("classes by variable on Z:",dict(byZ))
print("classes by ordered axis assignment:",dict(byorder))
# raw-unit fingerprints
V,I,R=120,0.5,240
print("c_VIR (A):",math.sqrt(V**2+I**2+R**2),"(mA):",math.sqrt(V**2+(I*1000)**2+R**2))
print("sqrt(5774400)=",math.sqrt(5774400),"=sqrt(2400^2+120^2)",math.sqrt(2400**2+120**2))
# bulb held at 240 ohm in the dim layers
for Rl in (2400,14400):
    Il=V/Rl; print("layer",Rl,"current",Il,"bulb power at 240 ohm",Il**2*240)
# log per-unit coordinates
v=math.log(120/120); i=math.log(0.05/0.5); r=math.log(2400/240)
print("2400-ohm layer point",(v,round(i,4),round(r,4)),"on plane v=i+r:",abs(v-(i+r))<1e-12)
# plane through centre perpendicular to (1,-1,-1) cuts cube [-1,1]^3 in a regular hexagon
n=np.array([1,-1,-1])/math.sqrt(3)
pts=[]
edges=[(a,b) for a in itertools.product([-1,1],repeat=3) for b in itertools.product([-1,1],repeat=3) if sum(x!=y for x,y in zip(a,b))==1 and a<b]
for a,b in edges:
    a=np.array(a,float); b=np.array(b,float)
    da=a@n; db=b@n
    if da*db<0:
        t=da/(da-db); pts.append(a+t*(b-a))
pts=np.array(pts); d=np.linalg.norm(pts,axis=1)
print("section vertices",len(pts),"radii",np.round(d,6))
