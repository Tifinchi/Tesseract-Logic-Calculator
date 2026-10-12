"""Check the best-case power a 50 Hz line can deliver to the 2012 ball (Cen, Yuan and Xue 2014)."""
import numpy as np, math
eps0=8.8541878128e-12; mu0=4e-7*math.pi
f=49.7; w=2*math.pi*f
# electric: P = 0.5*sigma*|Ein|^2 V, Ein = 3 E0/(3 - i s), s = sigma/(w eps0)
s=np.logspace(-3,3,20001); fac=0.5*9*s/(9+s*s)   # P/(w eps0 E0^2 V)
i=np.argmax(fac); print(f"electric: max P/(w eps0 E0^2 V) = {fac[i]:.4f} at s = {s[i]:.3f} (expect 0.75 at 3)")
def Pe(E0,r): return 0.75*w*eps0*E0**2*(4/3)*math.pi*r**3
for E0 in (133,150):
    print(f"  E0={E0} V/m: core r=0.125 m {Pe(E0,0.125):.2e} W; whole r=2.5 m {Pe(E0,2.5):.2e} W")
# magnetic: chi = -(3/2)[1 - 3/(ka)^2 + 3 cot(ka)/(ka)], k=(1-1j)/delta (e^{iwt}); P = (w/2) Im(-chi)... take |Im|
x=np.logspace(-2,2,40001)  # a/delta
ka=(1-1j)*x
chi=-(1.5)*(1-3/ka**2+3/(ka*np.tan(ka)))
im=np.abs(chi.imag); j=np.argmax(im); print(f"magnetic: max |Im chi| = {im[j]:.4f} at a/delta = {x[j]:.3f}")
for I in (100,1000):
    B0=mu0*I/(2*math.pi*20)*math.sqrt(2)  # peak field from rms current at 20 m (single conductor)
    for r in (0.125,2.5):
        V=(4/3)*math.pi*r**3
        print(f"  I={I} A rms, B0={B0:.2e} T, r={r}: P_max = {(w/2)*im[j]*B0**2*V/mu0:.2e} W")
