# Appendix C: converter proofs with PDG 2026 inputs, dimension steps, 4D, T-bar, complex cells, bands and clocks.
import numpy as np, math, sympy as sp
from scipy.integrate import solve_ivp
from scipy.special import ellipk
print("== Scott-T")
V=1.0; VA=1j*math.sqrt(3)/2*V; VB=-V/2; VC=V/2
lines=[VA-VB,VB-VC,VC-VA]
print([round(abs(x),12) for x in lines],[round(math.degrees(np.angle(lines[(i+1)%3]/lines[i])),6) for i in range(3)])
print("== third harmonics at 30%")
t=np.linspace(0,1,200001)
ph=[0,-2*np.pi/3,2*np.pi/3]
va=[np.cos(2*np.pi*t+p)+0.3*np.cos(3*(2*np.pi*t+p)) for p in ph]
K=(2/3)*np.array([[1,-0.5,-0.5],[0,math.sqrt(3)/2,-math.sqrt(3)/2]])
ab=K@np.array(va); r=np.sqrt(ab[0]**2+ab[1]**2)
diag=(va[0]+va[1]+va[2])/3
print("circle radius min/max",r.min().round(6),r.max().round(6),"diag swing",diag.min().round(4),diag.max().round(4),"neutral peak (sum)",(va[0]+va[1]+va[2]).max().round(4))
print("== 4D special case: centre-to-corner sqrt(n)/2 = edge only at n=4:",[n for n in range(1,10) if math.isclose(math.sqrt(n)/2,1)])
print("== Higgs: PDG 2026 inputs")
mW,mZ,mH,GF=80.3625,91.1880,125.13,1.1663788e-5
v=(math.sqrt(2)*GF)**-0.5
leg=math.sqrt(mZ**2-mW**2); thW=math.degrees(math.acos(mW/mZ))
g=2*mW/v; gp=2*leg/v; e=g*gp/math.sqrt(g*g+gp*gp)
lam=mH**2/(2*v*v)
print(f"v={v:.3f} leg g'v/2={leg:.3f} thetaW={thW:.3f} g={g:.5f} g'={gp:.5f} e={e:.5f} lambda={lam:.5f} sin2thW(on-shell)={1-(mW/mZ)**2:.5f}")
for sq in (300,1000,3000,10000):
    s=sq**2; a0no=s/v**2/(16*math.pi); a0h=-(mH**2/v**2)*s/(s-mH**2)/(16*math.pi)
    print(sq,round(a0no,4),round(a0h,5))
print("level -2lambda/16pi",round(-2*lam/(16*math.pi),5),"no-Higgs break sqrt(8pi)v",round(math.sqrt(8*math.pi)*v,1),"LQT bound",round(math.sqrt(16*math.pi/3)*v,1))
for E in (100,1000,10000):
    p=math.sqrt(E*E-mW*mW); print("eps_L",E,round(p/mW,3),round(E/mW,3),round((E/mW)**2-(p/mW)**2,9))
print("== 4D: boost check")
ct,x,b=5,3,0.6; gm=1/math.sqrt(1-b*b); ct2=gm*(ct-b*x); x2=gm*(x-b*ct)
k=math.sqrt((1+b)/(1-b)); u,w=ct+x,ct-x; u2,w2=ct2+x2,ct2-x2
print("s2",ct*ct-x*x,ct2**2-x2**2,"u'/u",u2/u,"1/k",1/k,"v'/v",w2/w,"k",k)
print("proper time 10s at 0.6c:",10*math.sqrt(1-0.36))
print("== T-bar flip rate, box 3x2x1")
X,Y,Z=3.0,2.0,1.0
I1=(Y*Y+Z*Z)/12; I2=(X*X+Z*Z)/12; I3=(X*X+Y*Y)/12
lam_f=math.sqrt((I3-I2)*(I2-I1)/(I1*I3))
d2xy,c2xy,d2yz,c2yz=X*X-Y*Y,X*X+Y*Y,Y*Y-Z*Z,Y*Y+Z*Z
print("lambda formula",lam_f,"rosetta form",math.sqrt(d2xy/c2xy*d2yz/c2yz),"lambda^2",d2xy/c2xy*d2yz/c2yz,"=3/13",3/13)
def euler(t,w):
    w1,w2,w3=w
    return [(I2-I3)*w2*w3/I1,(I3-I1)*w3*w1/I2,(I1-I2)*w1*w2/I3]
sol=solve_ivp(euler,[0,25],[1e-6,1.0,1e-6],rtol=1e-11,atol=1e-14,dense_output=True)
tt=np.linspace(0,12,4000); w1=np.abs(sol.sol(tt)[0])
m=(w1>1e-5)&(w1<1e-2)
slope=np.polyfit(tt[m],np.log(w1[m]),1)[0]
print("simulated growth rate",round(slope,4))
sol=solve_ivp(euler,[0,200],[0.01,1.0,0.01],rtol=1e-10,atol=1e-12,max_step=0.05)
w2=sol.y[1]; flips=np.sum(np.diff(np.sign(w2))!=0)
print("sign changes of omega_2 in 200 units (near-middle-axis start 0.01):",flips)
print("== dimension steps")
print("step angles deg",[round(math.degrees(math.atan(1/math.sqrt(k))),2) for k in (1,2,3,4)])
print("3D->4D triangle (K sqrt3, K, 2K) angles",round(math.degrees(math.atan(1/math.sqrt(3))),2))
for s in (1,math.sqrt(2),math.sqrt(3)):
    print("cuboctahedron of cube side",round(s,4),"radius",round(s*math.sqrt(2)/2,4),"edge",round(s*math.sqrt(2)/2,4),"corner radius",round(s*math.sqrt(3)/2,4))
print("4D simplex centroid distance",np.linalg.norm([0.75]*4))
print("== complex cell facts")
ok=all((a*a+b*b)**2==(a*a-b*b)**2+(2*a*b)**2 for a in range(1,11) for b in range(1,a+1))
print("(a2+b2)^2=(a2-b2)^2+(2ab)^2 all 55:",ok)
R,Xr=4.0,3.0; zs=math.hypot(R,Xr); zp=R*Xr/math.hypot(R,Xr)
print("series",zs,"parallel",zp,"PFs",R/zs,zp/R,"sum squares",(R/zs)**2+(zp/R)**2)
print("== bands and clocks")
L,C,Rr=0.1,100e-6,37.95; w0=1/math.sqrt(L*C); gmm=Rr/(2*L); print("RLC wd",round(math.sqrt(w0*w0-gmm*gmm),2))
for gg in (9.81,9.80665):
    kk=math.sin(math.radians(60)); T=4*math.sqrt(1/gg)*ellipk(kk*kk); print("pendulum 120deg L=1m g",gg,"T",round(T,4),"small-angle",round(2*math.pi*math.sqrt(1/gg),4))
ecc=0.6; M=2*math.pi/6
E=M
for _ in range(100): E=E-(E-ecc*math.sin(E)-M)/(1-ecc*math.cos(E))
nu=2*math.atan(math.sqrt((1+ecc)/(1-ecc))*math.tan(E/2))
print("Kepler e=0.6 t=T/6 true anomaly",round(math.degrees(nu)%360,2))
S=1000; PF=0.8; P=S*PF; Q=math.sqrt(S*S-P*P); print("single-phase P,Q,swing",P,Q,S)
print("three-phase flat power 120V 10A PF0.8:",3*120*10*0.8,"line voltages",[round(math.sqrt(3)*vv,1) for vv in (120,277,347)])
# 1/f from RC relaxations with log-spread time constants
rng=np.random.default_rng(1); taus=10**rng.uniform(-4,2,400)
f=np.logspace(-1,2,200); Sf=np.array([np.sum(taus/(1+(2*np.pi*fi*taus)**2)) for fi in f])
print("1/f slope (0.1-100 Hz, taus 1e-4..1e2 s)",round(np.polyfit(np.log10(f),np.log10(Sf),1)[0],3))
