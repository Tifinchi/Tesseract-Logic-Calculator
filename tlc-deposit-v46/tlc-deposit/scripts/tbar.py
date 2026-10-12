"""T-bar (Dzhanibekov) flip time against the distance from balance, and time = 2pi dJ/dE."""
import numpy as np
from scipy.integrate import solve_ivp
I = np.array([5, 10, 13])/12.0           # solid box 3 x 2 x 1, unit mass; y is the intermediate axis
lam = np.sqrt((I[1]-I[0])*(I[2]-I[1])/(I[0]*I[2]))   # growth rate at Omega = 1
print('box 3:2:1  I =', I, ' lambda =', lam, ' (sqrt(3/13) =', np.sqrt(3/13), ')')

def rhs(t, w):
    return [(I[1]-I[2])*w[1]*w[2]/I[0], (I[2]-I[0])*w[2]*w[0]/I[1], (I[0]-I[1])*w[0]*w[1]/I[2]]

def flip_times(eps, tmax=400):
    w0 = [eps, 1.0, eps]
    ev = lambda t, w: w[1]
    ev.direction = 0
    sol = solve_ivp(rhs, [0, tmax], w0, events=ev, rtol=1e-12, atol=1e-14, max_step=0.5)
    t = sol.t_events[0]
    return np.diff(t) if len(t) > 2 else None, t

rows = []
for eps in [1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8]:
    d, t = flip_times(eps)
    if d is None: continue
    Tflip = np.mean(d)    # time between successive sign changes of the middle-axis spin (one flip)
    rows.append((eps, Tflip))
    print(f'eps={eps:8.1e}  time between flips = {Tflip:10.5f}   lambda*T = {lam*Tflip:9.4f}   ln(1/eps) = {np.log(1/eps):8.4f}')
x = np.log(1/np.array([r[0] for r in rows])); y = np.array([r[1] for r in rows])
A = np.vstack([x, np.ones_like(x)]).T
slope, icpt = np.linalg.lstsq(A[-6:], y[-6:], rcond=None)[0]
print(f'fit (small eps): T = {slope:.5f} ln(1/eps) + {icpt:.4f};  slope*lambda = {slope*lam:.5f}  (2 expected: two passes of the unstable direction per flip... check)')

# action variable on the angular-momentum sphere: J = (L/2pi) * enclosed solid angle (canonical pair L_z, phi)
# check T = 2pi dJ/dE numerically for orbits circling the major axis z
L = 1.0
def solid_angle_about_z(E):
    # orbit: L1^2/I1 + L2^2/I2 + L3^2/I3 = 2E with |L| = 1, circling z (requires 2E > 1/I2 region near z... z has largest I -> lowest energy)
    phis = np.linspace(0, 2*np.pi, 20001)
    # in spherical coords about z: L3 = cos th, L1 = sin th cos phi, L2 = sin th sin phi
    # sin^2 th (cos^2/I1 + sin^2/I2) + cos^2 th / I3 = 2E -> sin^2 th = (2E - 1/I3)/(cos^2/I1 + sin^2/I2 - 1/I3)
    s2 = (2*E - 1/I[2])/(np.cos(phis)**2/I[0] + np.sin(phis)**2/I[1] - 1/I[2])
    s2 = np.clip(s2, 0, 1)
    costh = np.sqrt(1 - s2)
    return np.trapezoid(1 - costh, phis)      # solid angle of the cap enclosed around z
Emin = 0.5/I[2]; Esep = 0.5/I[1]
for frac in [0.5, 0.9, 0.99, 0.999]:
    E = Emin + frac*(Esep - Emin); h = 1e-7*(Esep-Emin)
    dJdE = (L/(2*np.pi))*(solid_angle_about_z(E+h) - solid_angle_about_z(E-h))/(2*h)
    T_action = 2*np.pi*dJdE
    # simulate: start on the orbit at phi = 0 (L2 = 0): sin^2 th = (2E - 1/I3)/(1/I1 - 1/I3)
    s2 = (2*E - 1/I[2])/(1/I[0] - 1/I[2]); L1 = np.sqrt(s2); L3 = np.sqrt(1 - s2)
    w0 = [L1/I[0], 0.0, L3/I[2]]
    ev = lambda t, w: w[1]; ev.direction = 1
    sol = solve_ivp(rhs, [0, 2000], w0, events=ev, rtol=1e-12, atol=1e-14, max_step=0.2)
    te = sol.t_events[0]; Tsim = np.mean(np.diff(te)) if len(te) > 2 else np.nan
    print(f'E at {frac:5.3f} of the way to the separatrix: 2pi dJ/dE = {T_action:10.5f}   simulated period = {Tsim:10.5f}')
# separatrix geometry: two great circles through the middle axis; lune solid angles
r = np.sqrt((1/I[0]-1/I[1])/(1/I[1]-1/I[2]))
half = np.arctan(1/r)     # half dihedral of the lune around z
print(f'separatrix planes: L3/L1 = +-{r:.5f}; lune about the major axis: dihedral {np.degrees(2*half):.3f} deg, solid angle {4*half:.5f} sr; about the minor axis {np.degrees(np.pi-2*half):.3f} deg, {2*(np.pi-2*half):.5f} sr')
