"""Floquet exponents of y'' + (A - 2q cos 2z) y = 0, A = k^2/m^2 + X/2, q = -3X/4 (Higgs driven by R after inflation)."""
import numpy as np, math
from scipy.integrate import solve_ivp
def mu(A,q):
    f=lambda z,Y:[Y[1],-(A-2*q*math.cos(2*z))*Y[0],Y[3],-(A-2*q*math.cos(2*z))*Y[2]]
    s=solve_ivp(f,[0,math.pi],[1,0,0,1],rtol=1e-11,atol=1e-13)
    M=np.array([[s.y[0,-1],s.y[2,-1]],[s.y[1,-1],s.y[3,-1]]])
    ev=np.linalg.eigvals(M); return max(0.0,math.log(max(abs(ev)))/math.pi)
for X in (0.1,0.3,0.5,1.0,1.33,2,5,10,20):
    m0=mu(X/2,-0.75*X)
    best=max(mu(K+X/2,-0.75*X) for K in np.linspace(0,max(4,2*X),801))
    print(f"X={X}: k=0 {m0:.3f}; max over k {best:.3f}")
# boundaries of the k=0 stable bands
xs=np.linspace(0.01,6,6000); st=[mu(X/2,-0.75*X)<1e-6 for X in xs]
edges=[xs[i] for i in range(1,len(xs)) if st[i]!=st[i-1]]
print("k=0 stability changes at X =", [round(e,3) for e in edges])
