"""The law of motion, rate = J·grad H, with H the sum table on its two legs:
H = sqrt(A)·sqrt(m²c⁴ + p²c²/B), A and B the time and space factors of the isotropic Schwarzschild metric.
Checks: Mercury's perihelion (three H's), light at the Sun's edge (two H's), the GPS clock rate, hydrogen's Lyman-alpha."""
import numpy as np
from scipy.integrate import solve_ivp

AS = 206264.80624709636                      # arcseconds per radian
GM_SUN = 1.32712440018e20                    # m^3/s^2 (IAU 2015 nominal, TDB-compatible)
C = 299792458.0
RG = GM_SUN / C**2                           # GM/c^2 in metres

def legs(rho):
    u = 1.0 / (2.0 * rho)
    A = ((1 - u) / (1 + u))**2
    B = (1 + u)**4
    dA_du = -4 * (1 - u) / (1 + u)**3
    dB_du = 4 * (1 + u)**3
    du_dr = -2 * u * u
    return A, B, dA_du * du_dr, dB_du * du_dr

def rhs_factory(kind, massless=False):
    def f(t, z):
        x, y, px, py = z
        rho = np.hypot(x, y)
        p2 = px * px + py * py
        A, B, dA, dB = legs(rho)
        if kind == 'time':
            B, dB = 1.0, 0.0
        if kind == 'newton':
            # H = 1 + p^2/2 - 1/rho
            dHdr = 1.0 / rho**2
            return [px, py, -dHdr * x / rho, -dHdr * y / rho]
        if massless:
            p = np.sqrt(p2)
            s = np.sqrt(A / B)                     # H = s·|p|
            ds = 0.5 / s * (dA / B - A * dB / B**2)
            return [s * px / p, s * py / p, -p * ds * x / rho, -p * ds * y / rho]
        S = np.sqrt(1 + p2 / B)
        sA = np.sqrt(A)
        dHdp = sA / (B * S)
        dHdr = dA / (2 * sA) * S + sA * (-p2 * dB / (2 * B * B * S))
        return [dHdp * px, dHdp * py, -dHdr * x / rho, -dHdr * y / rho]
    return f

# ---------------------------------------------------------------- Mercury
a_m = 5.7909050e10 / RG                      # semi-major axis in GM/c^2 (0.387098 au)
e_m = 0.205630
P_days = 87.9691                             # sidereal period
rp, ra = a_m * (1 - e_m), a_m * (1 + e_m)
L2 = (2 / ra - 2 / rp) / ((1 - 2 / ra) / ra**2 - (1 - 2 / rp) / rp**2)
L = np.sqrt(L2)
E = np.sqrt((1 - 2 / rp) * (1 + L2 / rp**2))
rho_p = (rp - 1 + np.sqrt(rp * rp - 2 * rp)) / 2
z0 = [rho_p, 0.0, 0.0, L / rho_p]
A0, B0, _, _ = legs(rho_p)
H0 = np.sqrt(A0) * np.sqrt(1 + (L / rho_p)**2 / B0)
print(f'Mercury: a = {a_m:.6e} GM/c², L = {L:.6e}, E − 1 = {E - 1:.6e}; H at perihelion − E = {H0 - E:.2e}')
T = 2 * np.pi * a_m**1.5
def peri(t, z): return z[0] * z[2] + z[1] * z[3]
peri.direction = 1
res = {}
for kind in ['newton', 'time', 'both']:
    sol = solve_ivp(rhs_factory(kind), [0, 1.6 * T], z0, method='DOP853', rtol=1e-13,
                    atol=[1e-7, 1e-7, 1e-19, 1e-19], events=peri, dense_output=False)
    ts = [t for t in sol.t_events[0] if t > 0.5 * T]
    zs = [z for t, z in zip(sol.t_events[0], sol.y_events[0]) if t > 0.5 * T]
    ang = np.arctan2(zs[0][1], zs[0][0])
    dphi = ang                                # after one full turn the angle wraps to a small positive number
    per_cy = dphi * (36525.0 / P_days) * AS
    res[kind] = per_cy
    print(f'  {kind:7s}: perihelion advance {dphi:.6e} rad per orbit = {per_cy:.3f}″ per century (period {ts[0] / T:.9f} T_Kepler)')
gr = 6 * np.pi / (a_m * (1 - e_m**2)) * (36525.0 / P_days) * AS
print(f'  closed form 6πGM/(c²a(1−e²)): {gr:.3f}″ per century; time leg alone (2+2γ−β)/3 = 1/3 of it: {gr / 3:.3f}″')

# ---------------------------------------------------------------- light at the Sun's edge
R_SUN = 6.957e8 / RG
for kind in ['time', 'both']:
    X = 2e4 * R_SUN
    def done(t, z): return z[0] - X
    done.terminal = True
    f = rhs_factory(kind, massless=True)
    sol = solve_ivp(f, [0, 10 * X], [-X, R_SUN, 1.0, 0.0], method='DOP853', rtol=1e-12, atol=[1e-6, 1e-6, 1e-16, 1e-16], events=done)
    zf = sol.y[:, -1]
    defl = -np.arctan2(zf[3], zf[2])
    # the bending left beyond ±X on both sides, (2·k/b)(1 − X/√(X²+b²)) each, k = 2 for both legs, 1 for the time leg
    k = 2 if kind == 'both' else 1
    tail = 2 * (k / R_SUN) * (1 - X / np.hypot(X, R_SUN))
    print(f'Light at the Sun\'s edge, {kind}: {defl * AS:.5f}″ (+ tail {tail * AS:.2e}″)')

# ---------------------------------------------------------------- GPS
GM_E = 3.986004418e14
a_gps = 26561.75e3
L_G = 6.969290134e-10                        # IAU 2000 B1.9: rate of TT (geoid clocks) against TCG
sat = 1 - np.sqrt(1 - 3 * GM_E / (a_gps * C**2))   # circular geodesic: dτ/dt = sqrt(1 − 3GM/rc²)
frac = L_G - sat
print(f'GPS: satellite clock gains {frac:.5e} (fraction) = {frac * 86400e6:.3f} μs per day on ground clocks')
print(f'     gravity part {(L_G - GM_E / (a_gps * C**2)) * 86400e6:.3f} μs/day (ground L_G minus satellite GM/ac²), speed part {-0.5 * GM_E / (a_gps * C**2) * 86400e6:.3f} μs/day')

# ---------------------------------------------------------------- hydrogen
alpha = 7.2973525643e-3
Rinf = 10973731.568157                       # 1/m (CODATA 2022)
me_mp = 1 / 1836.152673426
red = 1 / (1 + me_mp)
lam_s = 1 / (Rinf * red * 0.75)
print(f'Hydrogen Lyman-alpha, Schrödinger (reduced mass): {lam_s * 1e9:.5f} nm')
def dirac(n, j):
    k = j + 0.5
    d = k - np.sqrt(k * k - alpha**2)
    return (1 + (alpha / (n - d))**2)**-0.5 - 1
mec2_hz = 2 * Rinf * C / alpha**2            # m_e c² / h in Hz
L1S = 8.1727441e9                            # the 1S Lamb shift, 8,172,744.1 kHz (Maisenbacher et al. 2026)
for (j2, lab) in [(0.5, '2P1/2'), (1.5, '2P3/2')]:
    dE = (dirac(2, j2) - dirac(1, 0.5)) * mec2_hz * red    # leading reduced-mass scaling
    print(f'  Dirac {lab} → 1S1/2: {C / dE * 1e9:.5f} nm; with the 1S Lamb shift {C / (dE - L1S) * 1e9:.5f} nm')

# ---------------------------------------------------------------- the share and the hum together fix H0
G = 6.67430e-11; eV = 1.602176634e-19; hbarc = 197.3269804e-9; Mpc = 3.0856775814913673e22
OL = 1 - np.pi**2 / 32 - 9.1e-5
H0 = lambda meV: np.sqrt(8 * np.pi * G * ((meV * 1e-3)**4 / hbarc**3 * eV) / (3 * C**2 * OL)) * Mpc / 1e3
print(f'H0 from the share and the hum (2.02 ± 0.20 meV): {H0(2.02):.1f} +{H0(2.22) - H0(2.02):.1f} / -{H0(2.02) - H0(1.82):.1f} km/s/Mpc')
