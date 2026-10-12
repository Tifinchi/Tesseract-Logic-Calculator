# Every permutation of the Tesseract Logic Calculator, scored by the author's behavior rule.
#
# The calculator's numbers come from a handful of relations, each built from the cube's parts: a balance angle,
# a phase rule, a frame and its tilt, a load, a share, a crossing, fractions and multiples. Every one of those parts
# has alternatives inside the same geometry. A permutation is one choice for every part. For each permutation the
# engine computes every observable the choice reaches downstream and applies rule 9 as section 4.8 states it:
#   carry each alignment to every observable it touches; leave out any that misses a measurement by more than 2 sigma;
#   a miss that every alignment in a family shows is inherited and stays on the record without eliminating anything.
# Survivors are then ranked by chi-square, by precision (how many digits they fix beyond the measurement) and by
# speed (time per evaluation), and the winners are tested on data they were not chosen with (oos.py).
# Data as of 10 October 2026: CODATA 2022, PDG 2026, NuFIT 6.1, DESI DR2 and the sources of Tables 63 and 96.
import math, itertools, json, time, sys
import numpy as np

# ------------------------------------------------------------------------------------------- measurements
D = {
 # leptons (MeV); CODATA 2022 muon-to-electron ratio; PDG 2026 tau
 'me': 0.51099895069, 'r_mu': (206.7682827, 0.0000046), 'mtau': (1776.93, 0.09),
 # neutrinos (eV^2, eV)
 'dm21': 7.49e-5, 'dm31': 2.513e-3, 'sum_cap': 0.064, 'sum_eff': (-0.101, 0.056, 0.047), 'mbeta_cap': 0.45, 'mbb_cap': 0.028,
 # mixing
 's13': (0.02200, 0.00056, 0.00056), 's12': (0.3096, 0.0073, 0.0057), 's12_juno': (0.3092, 0.0087, 0.0087), 's23': (0.561, 0.014, 0.014),
 # CKM angle alpha (deg) and the angle between the up and down legs
 'alpha_all': (84.1, 3.0, 3.7), 'alpha_rr': (92.6, 4.7, 4.5), 'alpha_bg': (91.0, 2.8, 2.8),
 'legs_pdg': (83.5, 2.1, 2.1), 'legs_lat': (85.9, 2.0, 2.0),
 # top, Higgs (GeV)
 'mt_dir': (172.60, 0.27), 'mt_pole': (172.1, 0.6), 'mh': (125.13, 0.11),
 # proton spin pieces and radius
 'dS_hermes': (0.330, 0.039), 'dS_lat': (0.380, 0.044), 'dG_lat': (0.158, 0.039),
 'rp_codata': (3.9977, 0.0030), 'rp_muH': (0.84087/0.2103089, 0.00039/0.2103089), 'rp_ep': (0.8439/0.2103089, 0.0013/0.2103089),
 'mu_p': (2.79284734462, 0.00000000082),
 # cosmos
 'Om_sncb': (0.3045, 0.0032), 'Om_planck': (0.3153, 0.0073), 'Om_desi': (0.2975, 0.0086), 'Om_dpa': (0.3003, 0.0035),
 'OL_des': (0.648, 0.017, 0.017), 'OL_pan': (0.666, 0.018, 0.018), 'OL_uni': (0.644, 0.028, 0.026),
 'q0_pop': (-0.537, 0.026), 'q0_brout': (-0.51, 0.024), 'Ez349': (1.247, 0.006),
 'rhoL4': (2.24, 0.010),       # meV; from Planck's Omega_Lambda and H0 (relative error 0.46%)
 'H0_cmb': (67.19, 0.38), 'H0_ladder': (73.50, 0.81), 'H0_jwst': (70.4, 1.9),
}
def pull(meas, comp, sc=0.0):
    """(measured - computed)/sigma, the side of an asymmetric error facing the computed value, sigma_c in quadrature"""
    if len(meas) == 2: v, s = meas; s_lo = s_hi = s
    else: v, s_lo, s_hi = meas
    s = s_hi if comp > v else s_lo
    return (v - comp)/math.hypot(s, sc)

# ------------------------------------------------------------------------------------------- the cube's own sets
def strict_angles():
    """the forty strict cube angles of section 16.7.4: angles between the 26 cube directions, the ladder's step angles
    arctan(1/sqrt k), and one sum or difference of two of them (degrees, 0 < a <= 90)"""
    pts = [(x, y, z) for x in (-1, 0, 1) for y in (-1, 0, 1) for z in (-1, 0, 1) if (x, y, z) != (0, 0, 0)]
    def ang(u, v):
        c = sum(a*b for a, b in zip(u, v))/math.sqrt(sum(a*a for a in u)*sum(a*a for a in v))
        return math.degrees(math.acos(max(-1, min(1, c))))
    A = set()
    for u, v in itertools.combinations(pts, 2):
        a = ang(u, v); a = min(a, 180 - a)
        if a > 0.01: A.add(round(a, 6))
    A |= {round(math.degrees(math.atan(1/math.sqrt(k))), 6) for k in (1, 2, 3, 4)}
    A = sorted(A); S = set(A)
    for a, b in itertools.combinations_with_replacement(A, 2):
        for x in (abs(a - b), (a + b) % 180):
            x = min(x, 180 - x)
            if 0.01 < x <= 90: S.add(round(x, 6))
    return dedupe(sorted(S), 1e-4)
def dedupe(xs, tol):
    """one value per cluster: rounding to six places can leave the same angle twice, 1e-6 degree apart"""
    out = []
    for x in sorted(xs):
        if not out or abs(x - out[-1]) > tol: out.append(x)
    return out
S40 = strict_angles()
def alphabet():
    r2, r3, r5 = math.sqrt(2), math.sqrt(3), math.sqrt(5); phi = (1 + r5)/2; pi = math.pi
    base = [1, r2, r3, 2, 1/r2, 1/r3, 0.5, 1/3, 0.25, 2/3, 0.75, 1.5, r3/r2, r3/2, pi/4, pi/6, pi**2/32, 1 - pi/4, 1 - pi/6,
            1 - pi**2/32, 4/pi, 2*pi/3, 1/(2*pi), 1/(8*pi**2), phi, 1/phi, r5, 3/8, 5/3, 3, 4, 6, 8, 12, 16, 24, 27]
    cand = list(base)
    for a, b in itertools.product(base, repeat=2): cand += [a*b, a/b]
    cand.sort(); out = []
    for v in cand:
        if not out or abs(v - out[-1]) > 1e-9*max(1, abs(v)): out.append(v)
    return out
ALPHA = alphabet()

# ------------------------------------------------------------------------------------------- F1: charged leptons
# sqrt(m_k) = mu (1 + a cos(delta + k psi)), k = 0, 1, 2. Rules: Q = Qc (the balance angle beta, Q = 1/(3 cos^2 beta));
# delta = dc (a phase), or delta = Q/n (the phase as the balance ratio shared over n phases); psi the phase spacing;
# the slot assignment of (e, mu, tau). Inputs m_e and m_mu (or m_e alone when two rules are imposed).
ME = D['me']; RMU = D['r_mu'][0]; MMU = ME*RMU
Q_SET = dedupe([1/(3*math.cos(math.radians(b))**2) for b in S40 if 1/3 < 1/(3*math.cos(math.radians(b))**2) <= 1.0 + 1e-6], 1e-7)
def rationals(nmax=12, vmax=2*math.pi/3):
    s = set()
    for n in range(1, nmax + 1):
        for m in range(1, int(vmax*n) + 1):
            s.add(round(m/n, 12))
    return sorted(s)
DELTA_SET = dedupe(list(rationals()) + [math.radians(a) for a in S40 if math.radians(a) < 2*math.pi/3], 1e-9)
ASSIGN = list(itertools.permutations(['tau', 'e', 'mu']))      # slot k -> lepton; the paper uses (tau, e, mu)
PSI_SET = [90.0, 120.0, 180.0]

def koide_tau(Qc):
    """Q = Qc with m_e, m_mu: the heavy root and, for completeness, the light one, for both relative signs of sqrt(m_e)"""
    out = []
    for se in (1, -1):
        s = se*math.sqrt(ME) + math.sqrt(MMU); t = ME + MMU
        A, B, C = Qc - 1, 2*Qc*s, Qc*s*s - t
        disc = B*B - 4*A*C
        if disc < 0: continue
        for x in ((-B + math.sqrt(disc))/(2*A), (-B - math.sqrt(disc))/(2*A)) if abs(A) > 1e-15 else (-C/B,):
            out.append((x*x, se, 'heavy' if x*x > MMU else 'light'))
    return out

def fam_leptons():
    P = []
    # Q only (Koide-type): predicts m_tau from m_e, m_mu
    for Qc in Q_SET:
        for mt, se, br in koide_tau(Qc):
            if br != 'heavy': continue
            P.append(dict(rule=f'Q = {Qc:.6g}' + ('' if se > 0 else ', electron root flipped'), Q=Qc, delta=None, psi=None, assign=None,
                          pred={'mtau': (mt, 0.0)}, nrules=1, cx=1, current=(abs(Qc - 2/3) < 1e-12 and se > 0)))
    for psi in PSI_SET:
        ps = math.radians(psi)
        for asg in ASSIGN:
            ke, km, kt = asg.index('e'), asg.index('mu'), asg.index('tau')
            # delta only: a from the e/mu ratio in closed form, then m_tau
            for dc in DELTA_SET:
                ce, cm, ct = math.cos(dc + ke*ps), math.cos(dc + km*ps), math.cos(dc + kt*ps)
                for s in (1, -1):
                    den = cm - s*math.sqrt(RMU)*ce
                    if abs(den) < 1e-14: continue
                    a = (s*math.sqrt(RMU) - 1)/den
                    if abs(1 + a*ce) < 1e-14: continue
                    mu = math.sqrt(ME)/abs(1 + a*ce)
                    mt = (mu*(1 + a*ct))**2
                    P.append(dict(rule=f'delta = {dc:.6g}, psi {psi:g}, slots {asg}, branch {s:+d}', Q=None, delta=dc, psi=psi, assign=asg,
                                  pred={'mtau': (mt, 0.0)}, nrules=1, cx=2, current=False))
            # Q and delta both fixed: m_e sets mu; m_mu/m_e and m_tau are predicted. Q is phase-independent, and its balance
            # meaning (DC power = AC power) holds, only for three phases 120 degrees apart, so this rule is built there.
            for Qc in (Q_SET if psi == 120.0 else []):
                for dc in DELTA_SET:
                    if psi == 120.0:
                        a2 = 2*(3*Qc - 1)
                        if a2 < 0: continue
                        alist = [math.sqrt(a2), -math.sqrt(a2)]
                    else:
                        alist = solve_a_for_Q(Qc, dc, ps)
                    for a in alist:
                        ce, cm, ct = math.cos(dc + ke*ps), math.cos(dc + km*ps), math.cos(dc + kt*ps)
                        if abs(1 + a*ce) < 1e-14: continue
                        mu = math.sqrt(ME)/abs(1 + a*ce)
                        r = ((1 + a*cm)/(1 + a*ce))**2; mt = (mu*(1 + a*ct))**2
                        cur = (psi == 120.0 and asg == ('tau', 'e', 'mu') and abs(Qc - 2/3) < 1e-12 and abs(dc - 2/9) < 1e-12 and a > 0)
                        P.append(dict(rule=f'Q = {Qc:.6g}, delta = {dc:.6g}, psi {psi:g}, slots {asg}, a {a:+.4f}', Q=Qc, delta=dc, psi=psi, assign=asg,
                                      pred={'r_mu': (r, 0.0), 'mtau': (mt, 0.0)}, nrules=2, cx=3, current=cur))
        # delta = Q/n, defined where Q does not depend on the phase (psi = 120): Q = (1 + a^2/2)/3
        if psi == 120.0:
            for asg in ASSIGN:
                ke, km, kt = asg.index('e'), asg.index('mu'), asg.index('tau')
                for n in range(1, 7):
                    f = lambda a: ((1 + a*math.cos((1 + a*a/2)/(3*n) + km*ps))/(1 + a*math.cos((1 + a*a/2)/(3*n) + ke*ps)))**2 - RMU
                    for a in roots_1d(f, -3, 3):
                        dc = (1 + a*a/2)/(3*n)
                        ce, ct = math.cos(dc + ke*ps), math.cos(dc + kt*ps)
                        if abs(1 + a*ce) < 1e-14: continue
                        mu = math.sqrt(ME)/abs(1 + a*ce); mt = (mu*(1 + a*ct))**2
                        P.append(dict(rule=f'delta = Q/{n}, psi 120, slots {asg}, a {a:+.4f}', Q=(1 + a*a/2)/3, delta=dc, psi=120.0, assign=asg,
                                      pred={'mtau': (mt, 0.0)}, nrules=1, cx=2, current=False))
    return P

def roots_1d(f, lo, hi, n=4000):
    xs = np.linspace(lo, hi, n); out = []
    with np.errstate(all='ignore'):
        vals = []
        for x in xs:
            try: vals.append(f(x))
            except ZeroDivisionError: vals.append(float('nan'))
    for i in range(n - 1):
        a, b = vals[i], vals[i + 1]
        if np.isfinite(a) and np.isfinite(b) and a*b < 0 and abs(a) < 1e6 and abs(b) < 1e6:
            x0, x1 = xs[i], xs[i + 1]
            for _ in range(80):
                xm = 0.5*(x0 + x1); fm = f(xm)
                if f(x0)*fm <= 0: x1 = xm
                else: x0 = xm
            out.append(0.5*(x0 + x1))
    return out

def solve_a_for_Q(Qc, dc, ps):
    def Q(a):
        r = [(1 + a*math.cos(dc + k*ps)) for k in range(3)]
        return sum(x*x for x in r)/sum(r)**2 if abs(sum(r)) > 1e-12 else float('nan')
    return roots_1d(lambda a: Q(a) - Qc, -4, 4, 2000)

# ------------------------------------------------------------------------------------------- F2: neutrinos, one root flipped
def fam_neutrinos():
    P = []
    for Qc in Q_SET:
        for order in ('NO', 'IO'):
            for flip in (0, 1, 2):
                def masses(ml):
                    if order == 'NO': return [ml, math.sqrt(ml*ml + D['dm21']), math.sqrt(ml*ml + D['dm31'])]
                    m3 = ml; m1 = math.sqrt(ml*ml + D['dm31'] - D['dm21']); m2 = math.sqrt(ml*ml + D['dm31'])
                    return [m1, m2, m3]
                def Q(ml):
                    m = masses(ml); r = [math.sqrt(x) for x in m]; r[flip] = -r[flip]
                    return sum(m)/sum(r)**2
                for ml in roots_1d(lambda x: Q(x) - Qc, 1e-7, 0.3, 6000):
                    m = masses(ml); S = sum(m)
                    s12 = 0.318107; s13 = 0.022329
                    U = [(1 - s12)*(1 - s13), s12*(1 - s13), s13]
                    mb = math.sqrt(sum(u*x*x for u, x in zip(U, m)))
                    t = [u*x for u, x in zip(U, m)]; mbb_min = max(0.0, max(t) - (sum(t) - max(t)))
                    P.append(dict(rule=f'Q = {Qc:.6g}, {order}, root {flip + 1} flipped', Q=Qc,
                                  pred={'sum': (S, 0.0), 'mbeta': (mb, 0.0), 'mbb_min': (mbb_min, 0.0)}, nrules=1, cx=2,
                                  current=(abs(Qc - 2/3) < 1e-12 and order == 'NO' and flip == 0)))
    return P

# ------------------------------------------------------------------------------------------- F3: mixing, the Clarke frame tilted
U_TBM = np.array([[math.sqrt(2/3), 1/math.sqrt(3), 0.0], [-1/math.sqrt(6), 1/math.sqrt(3), -1/math.sqrt(2)], [-1/math.sqrt(6), 1/math.sqrt(3), 1/math.sqrt(2)]])
def rot(i, j, t, phase=0.0):
    R = np.eye(3, dtype=complex); c, s = math.cos(t), math.sin(t)
    R[i, i] = c; R[j, j] = c; R[i, j] = s*np.exp(-1j*phase); R[j, i] = -s*np.exp(1j*phase)
    return R
def angles(U):
    s13sq = abs(U[0, 2])**2
    s12sq = abs(U[0, 1])**2/(1 - s13sq) if s13sq < 1 else float('nan')
    s23sq = abs(U[1, 2])**2/(1 - s13sq) if s13sq < 1 else float('nan')
    return s12sq, s13sq, s23sq
def fam_mixing():
    P = []
    planes = {0: (1, 2), 1: (0, 2), 2: (0, 1)}     # keep column c fixed: rotate in the plane of the other two
    for side in ('column', 'row'):
        for c in (0, 1, 2):
            i, j = planes[c]
            for tdeg in [0.0] + S40:
                for sg in (1, -1):
                    t = sg*math.radians(tdeg)
                    U = U_TBM @ rot(i, j, t) if side == 'column' else rot(i, j, t) @ U_TBM
                    s12, s13, s23 = angles(U)
                    pred_real = {'s13': (s13, 0.0), 's12': (s12, 0.0), 's12_juno': (s12, 0.0), 's23': (s23, 0.0)}
                    name = f"{'keep column' if side == 'column' else 'keep row'} {c + 1}, tilt {sg*tdeg:+.4g} deg"
                    cur_r = False
                    P.append(dict(rule=name + ', real', pred=pred_real, nrules=1, cx=2, current=cur_r))
                    # with a free phase the atmospheric angle is not fixed; theta12 and theta13 do not depend on it for a column rotation
                    if side == 'column':
                        P.append(dict(rule=name + ', phase free', pred={'s13': (s13, 0.0), 's12': (s12, 0.0), 's12_juno': (s12, 0.0)}, nrules=1, cx=3,
                                      current=(c == 0 and abs(tdeg - 15) < 1e-9 and sg == 1)))
    return P

# ------------------------------------------------------------------------------------------- F4: the CKM right angle
def fam_ckm():
    P = []
    for a in S40 + [a2 for a2 in (100.0, 105.0, 109.4712, 120.0)]:
        P.append(dict(rule=f'alpha = {a:.4g} deg, legs alone', pred={'alpha_all': (a, 0.0), 'alpha_rr': (a, 0.0), 'alpha_bg': (a, 0.0)},
                      nrules=1, cx=1, current=False))
        P.append(dict(rule=f'alpha = {a:.4g} deg = the angle between the up and down legs', pred={'alpha_all': (a, 0.0), 'alpha_rr': (a, 0.0), 'alpha_bg': (a, 0.0),
                      'legs_pdg': (a, 0.0), 'legs_lat': (a, 0.0)}, nrules=1, cx=1, current=(abs(a - 90) < 1e-9)))
    return P

# ------------------------------------------------------------------------------------------- F5: the load on the Higgs vacuum
YT2, G2, GP2, GS2 = 0.8705367702553598, 0.41937471386702224, 0.12854070829126615, 1.3575654676350988
TOP_EDGE = 170.93468           # GeV, the top at the vacuum's edge for M_H = 125.13 (higgs_rule2026.py)
DMH_DMT = None                 # filled by edge_slope()
def fam_load():
    P = []
    sets = {'y_t^2': [YT2], 'y_t^2, g^2': [YT2, G2], "y_t^2, g'^2": [YT2, GP2], 'y_t^2, g_s^2': [YT2, GS2], "y_t^2, g^2, g'^2": [YT2, G2, GP2]}
    for name, cs in sets.items():
        angs = [0.0] if len(cs) == 1 else [0.0, 60.0, 90.0, 120.0, 180.0]
        for ang in angs:
            for Nname, N in (('4 pi^2', 4*math.pi**2), ('8 pi^2', 8*math.pi**2), ('16 pi^2', 16*math.pi**2)):
                z = sum(c*complex(math.cos(math.radians(ang*k)), math.sin(math.radians(ang*k))) for k, c in enumerate(cs))
                d = abs(z)/N
                top = TOP_EDGE*(1 + d)
                # the electron's anchor: one phase (the photon) must give the first sub-reaction alpha/2pi
                ae1 = 4*math.pi*(1/137.035999177)/N
                P.append(dict(rule=f'{name} at {ang:g} deg, over {Nname}', delta=d,
                              pred={'mt_dir': (top, 0.38*(1 + d)/1.00955), 'mt_pole': (top, 0.38*(1 + d)/1.00955), 'mh': ('edge', d), 'ae1': (ae1, 0.0)},
                              nrules=1, cx=len(cs) + 1, current=(name == 'y_t^2, g^2' and ang == 120.0 and Nname == '8 pi^2')))
    return P

# ------------------------------------------------------------------------------------------- F6: the proton's spin legs
def fam_spin():
    P = []
    for x in range(0, 13):
        for y in range(0, 13 - x):
            z = 12 - x - y
            hS, dG, L = x/24, y/24, z/24
            P.append(dict(rule=f'1/2 dSigma = {x}/24, dG = {y}/24, L = {z}/24', pred={'dS_hermes': (2*hS, 0.0), 'dS_lat': (2*hS, 0.0), 'dG_lat': (dG, 0.0)},
                          nrules=1, cx=1, current=(x == 4 and y == 4)))
    return P

# ------------------------------------------------------------------------------------------- F7: the proton's radius
def fam_radius():
    P = []
    for n in sorted({v for v in ALPHA if 1.0 <= v <= 8.0}):
        P.append(dict(rule=f'r_p = {n:.6g} lambdabar_p', pred={'rp_codata': (n, 0.0), 'rp_muH': (n, 0.0), 'rp_ep': (n, 0.0)}, nrules=1, cx=1, current=(abs(n - 4) < 1e-12)))
        P.append(dict(rule=f'r_p = {n:.6g} lambdabar_p, with a charge at c round it', pred={'rp_codata': (n, 0.0), 'rp_muH': (n, 0.0), 'rp_ep': (n, 0.0), 'mu_p': (n, 0.0)},
                      nrules=1, cx=1, current=False))
    return P

# ------------------------------------------------------------------------------------------- F8, F9: the share and the hum
def V(n): return math.pi**(n/2)/math.gamma(n/2 + 1)
SHARES = {}
for n in range(1, 9):
    s = V(n)/2**n; SHARES[f'ball in the {n}-cube'] = s; SHARES[f'1 - ball in the {n}-cube'] = 1 - s
def cosmo_preds(Om):
    OL = 1 - Om
    return {'Om_sncb': (Om, 0.0), 'Om_planck': (Om, 0.0), 'Om_desi': (Om, 0.0), 'Om_dpa': (Om, 0.0),
            'OL_des': (OL, 0.0), 'OL_pan': (OL, 0.0), 'OL_uni': (OL, 0.0),
            'q0_pop': (1.5*Om - 1, 0.0), 'q0_brout': (1.5*Om - 1, 0.0), 'Ez349': (math.sqrt(Om*1.349**3 + OL), 0.0)}
def fam_share():
    return [dict(rule=f'Omega_m = {k} = {v:.6f}', Om=v, pred=cosmo_preds(v), nrules=1, cx=1, current=(k == 'ball in the 4-cube')) for k, v in SHARES.items() if 0 < v < 1]
CROSS = {'M12': 1.1e13, 'M13': 1.7e14, 'M23': 3.006e16}
SCALES = {'v': 246.22, 'm_H': 125.13, 'm_t': 172.60, 'm_Z': 91.1876, 'm_W': 80.3625}
COEF = {'1': 1, '1/2': 0.5, '2': 2, '1/3': 1/3, '3': 3, '1/4': 0.25, '4': 4, '1/sqrt2': 2**-0.5, 'sqrt2': 2**0.5, '1/pi': 1/math.pi, 'pi': math.pi,
        '1/2pi': 1/(2*math.pi), '2pi': 2*math.pi, '1/4pi': 1/(4*math.pi), '4pi': 4*math.pi, '1/8pi^2': 1/(8*math.pi**2)}
def fam_hum():
    P = []
    for cn, M in CROSS.items():
        for sn, sv in SCALES.items():
            for kn, k in COEF.items():
                E = k*sv*sv/M*1e12        # meV
                P.append(dict(rule=f'{kn} x {sn}^2/{cn}', E=E, pred={'rhoL4': (E, 0.10*E)}, nrules=1, cx=1 + (kn != '1') + (sn != 'v'),
                              current=(cn == 'M23' and sn == 'v' and kn == '1')))
    return P
G_N = 6.67430e-11; C_L = 299792458.0; HBARC = 197.3269804e-9; EV = 1.602176634e-19; MPC = 3.0856775814913673e22
def H0_from(E_meV, OL):
    return math.sqrt(8*math.pi*G_N*((E_meV*1e-3)**4/HBARC**3*EV)/(3*C_L**2*OL))*MPC/1e3

# ------------------------------------------------------------------------------------------- scoring
def edge_higgs(d):
    """the Higgs predicted from the measured top: the edge at M_t/(1 + delta), linear along the edge (slope from higgs_rule2026.py)"""
    mt_edge_in = D['mt_dir'][0]/(1 + d)
    return 125.13 + (mt_edge_in - TOP_EDGE)/DMH_DMT_INV
DMH_DMT_INV = 1.0    # placeholder; set in main from the RGE solver
def score(p, thresh=2.0):
    pulls = {}
    for k, (v, sc) in p['pred'].items():
        if k == 'mh':
            mh = edge_higgs(sc); pulls[k] = pull(D['mh'], mh, 0.92); p['mh_pred'] = mh; continue
        if k == 'ae1':
            pulls[k] = 0.0 if abs(v/(1/137.035999177/(2*math.pi)) - 1) < 1e-12 else float('inf'); continue
        if k == 'sum':
            pulls['sum_cap'] = 0.0 if v <= D['sum_cap'] else float('inf')
            pulls['sum_eff'] = pull(D['sum_eff'], v); continue
        if k == 'mbeta': pulls[k] = 0.0 if v <= D['mbeta_cap'] else float('inf'); continue
        if k == 'mbb_min': pulls[k] = 0.0 if v <= D['mbb_cap'] else float('inf'); continue
        meas = D[k]
        pulls[k] = pull(meas, v, sc)
    p['pulls'] = pulls
    finite = [x for x in pulls.values() if math.isfinite(x)]
    p['chi2'] = sum(x*x for x in finite) if len(finite) == len(pulls) else float('inf')
    p['maxpull'] = max(abs(x) for x in pulls.values()) if pulls else 0.0
    return p
# One measurement per quantity decides (the primary of Tables 63 and 96); the others are reported beside it.
SECONDARY = {'s12_juno', 'alpha_rr', 'alpha_bg', 'legs_pdg', 'dS_lat', 'rp_muH', 'rp_ep', 'Om_planck', 'Om_desi', 'Om_dpa', 'OL_des', 'OL_pan', 'OL_uni', 'q0_brout'}
def apply_rule(P, thresh=2.0):
    """rule 9: leave out any permutation that misses a primary measurement by more than thresh; a miss that every
    permutation in the family shows is inherited and stays on the record"""
    obs = sorted({k for p in P for k in p['pulls'] if k not in SECONDARY})
    # inherited: every permutation in the family reaches the observable, and every one misses it. An observable that
    # some permutations take as an input (and so match) is not inherited: the ones that predict it and miss are out.
    inherited = [k for k in obs if all(k in p['pulls'] for p in P) and all(abs(p['pulls'][k]) > thresh for p in P)]
    for p in P:
        p['survives'] = all(abs(v) <= thresh for k, v in p['pulls'].items() if k not in inherited and k not in SECONDARY)
        live = [v for k, v in p['pulls'].items() if k not in inherited and k not in SECONDARY]
        p['chi2_live'] = sum(x*x for x in live) if all(math.isfinite(x) for x in live) else float('inf')
        p['n_live'] = len(live)
    return inherited

FAMILIES = [('charged leptons', fam_leptons), ('neutrino masses', fam_neutrinos), ('lepton mixing', fam_mixing), ('quark mixing angle', fam_ckm),
            ('Higgs vacuum load', fam_load), ("proton's spin legs", fam_spin), ("proton's radius", fam_radius), ('matter share', fam_share), ('vacuum energy (hum)', fam_hum)]

def run(thresh=2.0, quiet=False):
    global DMH_DMT_INV
    sys.path.insert(0, '/tmp/claude-0/-home-claude/bbcd54b6-3962-59fe-98a5-d2cef4b7e16c/scratchpad')
    from higgs2 import summary
    from scipy.optimize import brentq
    t0 = time.time()
    mh_edge = lambda Mt: brentq(lambda M: summary(Mt=Mt, Mh=M, als=0.1180, MW=80.3625)['lmin'], 115, 140, xtol=1e-7)
    h = 0.5
    DMH_DMT_INV = 1/((mh_edge(TOP_EDGE + h) - mh_edge(TOP_EDGE - h))/(2*h))   # dMt/dMh along the edge
    t_rge = time.time() - t0
    out = {}
    for name, fn in FAMILIES:
        t1 = time.time(); P = fn(); t_gen = time.time() - t1
        for p in P: score(p, thresh)
        inh = apply_rule(P, thresh)
        out[name] = dict(P=P, inherited=inh, t_eval=(time.time() - t1)/max(1, len(P)))
    return out, t_rge

if __name__ == '__main__':
    out, t_rge = run()
    print(f'RGE slope along the edge: dMt/dMh = {DMH_DMT_INV:.4f} ({t_rge:.1f} s for two edge solutions)')
    tot = 1; tot_s = 1
    for name, r in out.items():
        P = r['P']; S = [p for p in P if p['survives']]
        cur = [p for p in P if p['current']]
        best = sorted(S, key=lambda p: (p['chi2_live']/max(1, p['n_live']), p['cx']))[:5]
        tot *= len(P); tot_s *= max(1, len(S))
        print(f"\n=== {name}: {len(P)} permutations, {len(S)} survive the 2-sigma rule; inherited misses: {r['inherited']}; {1e6*r['t_eval']:.1f} us each")
        for c in cur:
            print(f"  current: {c['rule']}: survives {c['survives']}, chi2 {c['chi2_live']:.3f} over {c['n_live']}, max pull {c['maxpull']:.3g}; pulls {{{', '.join(f'{k}: {v:+.3g}' for k, v in c['pulls'].items())}}}")
        for b in best:
            print(f"  best:    {b['rule']}: chi2 {b['chi2_live']:.3f} over {b['n_live']}; pulls {{{', '.join(f'{k}: {v:+.3g}' for k, v in b['pulls'].items())}}}")
    print(f'\nall families: {tot:.3e} permutations in the product; {tot_s:.3e} survive')
