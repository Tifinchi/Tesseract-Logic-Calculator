"""The four-layer ring, checked exactly (second part of the corrected mechanics).

Builds on tesseract.py: the cube is the base cell (w = -1) of the tesseract (+-1)^4 and
J: (x, y, z, w) -> (z, -w, -x, y) is the quarter turn fixed by the fabric (q -> i, phi -> v).
"""
import itertools
import numpy as np
from mpmath import mp, mpf, sqrt, pi

unit = {'+': (0, 0, 1), '−': (1, 0, 1), '×': (1, 1, 1), '*': (0, 1, 1),
        '’': (0, 0, 0), '”': (1, 0, 0), ',': (1, 1, 0), '.': (0, 1, 0)}
fabric = {'”': 'q', '.': 'i', '+': 'φ', '×': 'v', '’': 'L/M', '*': 'L/R', ',': 'RC', '−': 'MC'}
faces = {'Top': '+−×*', 'Front': '+’”−', 'Right': '−”,×', 'Left': '*.’+', 'Back': '×*.,', 'Bottom': '’”,.'}
base = {s: np.array([2 * a - 1 for a in p] + [-1]) for s, p in unit.items()}
J = np.array([[0, 0, 1, 0], [0, 0, 0, -1], [-1, 0, 0, 0], [0, 1, 0, 0]])
V = [np.array(v) for v in itertools.product([-1, 1], repeat=4)]

def label(u):
    u = tuple(int(x) for x in u)
    for s, b in base.items():
        if u == tuple(b):
            return s
        if u == tuple(-b):
            return '−' + s
    return str(u)

print('=== 1. J is unique: q, i, phi, v (lifted to w = -1) are an orthogonal basis of R^4')
a, b, c, d = base['”'], base['.'], base['+'], base['×']
Mabcd = np.array([a, b, c, d]).T
print('Gram matrix:\n', (Mabcd.T @ Mabcd).astype(int))
Jsolve = np.array([b, -a, d, -c]).T @ np.linalg.inv(Mabcd)
print('J from Ja=b, Jb=-a, Jc=d, Jd=-c equals the matrix used:', np.allclose(Jsolve, J))

print('\n=== 2. The variable frame: every tesseract vertex in terms of q, i, phi, v')
def coords(u):  # coefficients in the (q, i, phi, v) basis
    return tuple(float(x) for x in np.linalg.solve(Mabcd, u))
for s in ['”', '.', '+', '×', '’', '*', ',', '−']:
    co = coords(base[s])
    print(f'  {s} {fabric[s]:4s} = ' + ' '.join(f'{x:+.1f}{n}' for x, n in zip(co, ['q', 'i', 'φ', 'v'])))
odd = [v for v in V if np.sum(v < 0) % 2 == 1]
print('odd vertices are exactly ±q, ±i, ±φ, ±v:', {tuple(int(x) for x in v) for v in odd} == {tuple(s * u) for u in (a, b, c, d) for s in (1, -1)})
evenv = [v for v in V if np.sum(v < 0) % 2 == 0]
ok = all(sorted(abs(x) for x in coords(v)) == [0.5] * 4 and sum(1 for x in coords(v) if x < 0) % 2 == 1 for v in evenv)
print('even vertices (circuits) are ½(±q±i±φ±v) with an odd number of minus signs:', ok)
axis = [2 * s * np.eye(4, dtype=int)[k] for k in range(4) for s in (1, -1)]
ok = all(sorted(abs(x) for x in coords(p)) == [0.5] * 4 and sum(1 for x in coords(p) if x < 0) % 2 == 0 for p in axis)
print('the 8 axis points (cell centres pushed to radius 2) are ½(±q±i±φ±v) with an even number of minus signs:', ok)

print('\n=== 3. The eight cells as sign patterns of (q, i, φ, v); J carries them in two rings of four')
cells = {}
for k, nm in enumerate('xyzw'):
    for s in (-1, 1):
        cen = s * np.eye(4, dtype=int)[k]
        pat = ''.join('+' if x > 0 else '−' for x in coords(2 * cen))
        cells[(nm, s)] = pat
        members = [label(v) for v in V if v[k] == s]
        print(f'  cell {nm}={s:+d}: pattern (q i φ v) = {pat}; corners {" ".join(members)}')
def jcell(nm, s):
    cen = s * np.eye(4, dtype=int)['xyzw'.index(nm)]
    img = J @ cen
    k = int(np.nonzero(img)[0][0])
    return ('xyzw'[k], int(img[k]))
for start in [('w', -1), ('x', -1)]:
    ring = [start]
    for _ in range(4):
        ring.append(jcell(*ring[-1]))
    print('  ring:', ' -> '.join(f'{n}={s:+d} ({cells[(n, s)]})' for n, s in ring))
print('shared faces around the first ring:')
ringcells = [('w', -1), ('y', 1), ('w', 1), ('y', -1)]
for (n1, s1), (n2, s2) in zip(ringcells, ringcells[1:] + ringcells[:1]):
    k1, k2 = 'xyzw'.index(n1), 'xyzw'.index(n2)
    shared = [label(v) for v in V if v[k1] == s1 and v[k2] == s2]
    print(f'  {n1}={s1:+d} ∩ {n2}={s2:+d}: {" ".join(shared)}')

print('\n=== 4. Inside the base cube, J sends the Front face onto the Back face, letter for letter')
for f in ['Front', 'Back', 'Top', 'Bottom', 'Left', 'Right']:
    img = ''.join(label(J @ base[ch]) if not label(J @ base[ch]).startswith('−') or len(label(J @ base[ch])) == 1 else label(J @ base[ch]) for ch in faces[f])
    imgs = [label(J @ base[ch]) for ch in faces[f]]
    print(f'  J({f} {faces[f]}) = {" ".join(imgs)}')
print('  J(Front) equals the author\'s Back string as written:', [label(J @ base[ch]) for ch in faces['Front']] == list(faces['Back']))
print('  J(Back) equals minus the Front string as written:', [label(J @ base[ch]) for ch in faces['Back']] == ['−' + ch for ch in faces['Front']])

print('\n=== 5. J = a quarter turn about the Front–Back axis (x,z plane) times a quarter turn of that axis into the held axis (y,w plane)')
Rxz = np.array([[0, 0, 1, 0], [0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 0, 1]])
Ryw = np.array([[1, 0, 0, 0], [0, 0, 0, -1], [0, 0, 1, 0], [0, 1, 0, 0]])
print('J = Rxz @ Ryw:', np.array_equal(J, Rxz @ Ryw), '; the two commute:', np.array_equal(Rxz @ Ryw, Ryw @ Rxz))
print('Front faces y = -1 and Back y = +1 (unit y = 0 and 1):', all(unit[ch][1] == 0 for ch in faces['Front']), all(unit[ch][1] == 1 for ch in faces['Back']))

print('\n=== 6. The 3D shadow of the inverse (J^2) is the point inversion; each inverse sits one w-step over the opposite corner')
for s in base:
    inv = -base[s]
    over = [t for t, bb in base.items() if tuple(bb[:3]) == tuple(inv[:3])][0]
    print(f'  −{s} (−{fabric[s]}) sits over {over} ({fabric[over]})', end=';')
print()
print('J^2 = -I (det +1, a rotation in 4D); -I in 3D has det', round(np.linalg.det(-np.eye(3))))
flip = np.diag([-1, 1, 1, -1])
print('one-root mirror x -> -x with w -> -w: det', round(np.linalg.det(flip)), '= a half turn in the (x, w) plane; it sends the base cell to w = +1:', all((flip @ v)[3] == 1 for v in base.values()))

print('\n=== 7. Circles, torus, 24-cell')
fibres = {frozenset(label(np.linalg.matrix_power(J, k) @ v) for k in range(4)) for v in V}
print('vertex circles:', [sorted(f) for f in fibres])
# Clifford parallel: distance from a point on one circle to another circle is constant
def circle(p, t):
    return np.cos(t) * p + np.sin(t) * (J @ p)
ts = np.linspace(0, 2 * np.pi, 721)
dmin = []
for t0 in np.linspace(0, 2 * np.pi, 13):
    p = circle(a.astype(float), t0)
    dmin.append(min(np.linalg.norm(p - circle(c.astype(float), t)) for t in ts))
print('distance from the charge circle to the flux circle, sampled along the charge circle:', np.round(dmin, 3))
mp.dps = 30
print('torus radii |(x,z)| = |(y,w)| = sqrt2; area (2 pi sqrt2)^2 =', mp.nstr((2 * pi * sqrt(2)) ** 2, 12), '= 8 pi^2 =', mp.nstr(8 * pi ** 2, 12))
print('unit-sphere Clifford torus area 2 pi^2 scaled by r^2 = 4:', mp.nstr(2 * pi ** 2 * 4, 12))
pts = [tuple(v) for v in V] + [tuple(p) for p in axis]
d2 = [[float(np.sum((np.array(p) - np.array(q)) ** 2)) for q in pts] for p in pts]
nn = {sum(1 for x in row if abs(x - 4) < 1e-9) for row in d2}
print('24 points at radius 2, nearest neighbours per point:', nn)
sets = {'variables': [tuple(v) for v in odd], 'circuits': [tuple(v) for v in evenv], 'cells': [tuple(p) for p in axis]}
for nm, S in sets.items():
    Sset = set(S)
    print(f'  {nm}: 8 points; J keeps the set: {all(tuple(int(x) for x in J @ np.array(p)) in Sset for p in S)}; a 16-cell (each point has 6 neighbours at distance 2√2 and one antipode): {all(sum(1 for q in S if abs(np.sum((np.array(p) - np.array(q)) ** 2) - 8) < 1e-9) == 6 for p in S)}')

print('\n=== 8. The electron on the fabric corners (exact SI, CODATA 2022)')
h = mpf('6.62607015e-34'); e = mpf('1.602176634e-19'); c0 = mpf(299792458)
hbar = h / (2 * pi); me = mpf('9.1093837139e-31'); mpr = mpf('1.67262192595e-27')
alpha = mpf('7.2973525643e-3'); eps0 = e ** 2 / (2 * alpha * h * c0); mu0 = 1 / (eps0 * c0 ** 2)
lam = hbar / (me * c0); r_e = alpha * lam; a0 = lam / alpha
C = 4 * pi * eps0 * r_e; L = pi * mu0 * a0
w = 1 / sqrt(L * C)
print('omega = 1/sqrt(L C) =', mp.nstr(w, 12), 's^-1; m_e c^2/h =', mp.nstr(me * c0 ** 2 / h, 12))
q = e; i = w * q; phi = L * i; v = q / C
print('q =', mp.nstr(q, 10), 'C;  i = omega q =', mp.nstr(i, 10), 'A;  phi = L i =', mp.nstr(phi, 12), 'Wb (h/e =', mp.nstr(h / e, 12), ');  v = q/C =', mp.nstr(v, 12), 'V (m_e c^2/e =', mp.nstr(me * c0 ** 2 / e, 12), ')')
print('faces: C = q/v =', mp.nstr(q / v, 8), '; L = phi/i =', mp.nstr(phi / i, 8), '; R = v/i =', mp.nstr(v / i, 12), '; M = phi/q =', mp.nstr(phi / q, 12), '; h/e^2 =', mp.nstr(h / e ** 2, 12))
print('rate faces: v/phi =', mp.nstr(v / phi, 10), ' i/q =', mp.nstr(i / q, 10), '(both omega)')
print('energies q v and phi i:', mp.nstr(q * v, 10), mp.nstr(phi * i, 10), ' m_e c^2 =', mp.nstr(me * c0 ** 2, 10), ' (peak store q^2/2C =', mp.nstr(q ** 2 / (2 * C), 8), ')')
print('q phi =', mp.nstr(q * phi, 10), '= h')
muB = e * hbar / (2 * me)
ring_moment = i * pi * lam ** 2
print('one charge per Compton period round the Compton circle: i pi lambda-bar^2 =', mp.nstr(ring_moment, 12), '; mu_B =', mp.nstr(muB, 12), '; ratio', mp.nstr(ring_moment / muB, 15))
print('speed round the circle: 2 pi lambda-bar * (m_e c^2/h) / c =', mp.nstr(2 * pi * lam * (me * c0 ** 2 / h) / c0, 15))
print('angular momentum of that ring m_e c lambda-bar / hbar =', mp.nstr(me * c0 * lam / hbar, 10), '(the electron spin is 1/2)')
print('measured |mu_e|/mu_B = g/2 = 1.00115965218059(13); excess over the ring:', mp.nstr(mpf('1.00115965218059') - 1, 12), '; alpha/2pi =', mp.nstr(alpha / (2 * pi), 10))
lam_p = hbar / (mpr * c0); muN = e * hbar / (2 * mpr)
print('proton: a charge at c round r = 4 lambda-bar_p gives e c r/2 =', mp.nstr(e * c0 * 4 * lam_p / 2 / muN, 10), 'mu_N; measured 2.79284734462(82)')
print('2 pi (4 lambda-bar_p) / (h/(m_p c)) =', mp.nstr(2 * pi * 4 * lam_p / (h / (mpr * c0)), 10), 'Compton wavelengths: one per layer')

print('\n=== 9. Inverse-cell (CPT) pulls, (measured - model)/sigma, model = equality')
def pull(x, s, m=0):
    return (mpf(x) - mpf(m)) / mpf(s)
rows = [
    ('positron g vs electron g, (g+ - g-)/g = (-0.5 ± 2.1)e-12', pull('-0.5', '2.1')),
    ('antimuon g, (g+ - g-)/g = (-0.11 ± 0.12)e-8', pull('-0.11', '0.12')),
    ('muon lifetimes (tau+ - tau-)/tau = (2 ± 8)e-5', pull('2', '8')),
    ('antiproton moment |mu_pbar| - mu_p, in units of combined sigma', (mpf('2.7928473441') - mpf('2.79284734462')) / sqrt(mpf('0.0000000042') ** 2 + mpf('0.00000000082') ** 2)),
    ('top - antitop mass, -0.15 ± sqrt(0.19^2 + 0.09^2) GeV', mpf('-0.15') / sqrt(mpf('0.19') ** 2 + mpf('0.09') ** 2)),
    ('  that sigma (GeV)', sqrt(mpf('0.19') ** 2 + mpf('0.09') ** 2)),
    ('  relative to 172.56 GeV', sqrt(mpf('0.19') ** 2 + mpf('0.09') ** 2) / mpf('172.56')),
    ('antihydrogen fall a_g: (0.75 - 1)/sqrt(0.13^2 + 0.16^2)', (mpf('0.75') - 1) / sqrt(mpf('0.13') ** 2 + mpf('0.16') ** 2)),
    ('  that sigma (g)', sqrt(mpf('0.13') ** 2 + mpf('0.16') ** 2)),
    ('solar nu vs reactor nubar Dm21: (6.10 - 7.50)/sqrt(0.95^2 + 0.12^2) [JUNO]', (mpf('6.10') - mpf('7.50')) / sqrt(mpf('0.95') ** 2 + mpf('0.12') ** 2)),
    ('  vs KamLAND 7.58, sigma sqrt(0.14^2+0.15^2)', (mpf('6.10') - mpf('7.58')) / sqrt(mpf('0.95') ** 2 + mpf('0.14') ** 2 + mpf('0.15') ** 2)),
    ('model Dm21 from 0.364 and 8.66 meV (1e-5 eV^2), a model input', (mpf('8.66') ** 2 - mpf('0.364') ** 2) / 10),
    ('Higgs width CMS 3.2 +2.4/-1.7 vs SM 4.1 MeV (upper error)', pull('3.2', '2.4', '4.1')),
    ('S2 Schwarzschild precession f_SP 1.10 ± 0.19 vs GR 1', pull('1.10', '0.19', 1)),
]
for lab, val in rows:
    print(f'  {lab:78s} {mp.nstr(val, 6)}')

print('\n=== 10. Cosmos ring with Omega_m = pi^2/32 (flat, w = -1)')
Om = pi ** 2 / 32; OL = 1 - Om
print('q0 = 3/2 Om - 1 =', mp.nstr(1.5 * Om - 1, 8), '; with DESI 0.2975:', mp.nstr(1.5 * mpf('0.2975') - 1, 6), '; Planck 0.315:', mp.nstr(1.5 * mpf('0.315') - 1, 6))
print('acceleration starts (q = 0) at z =', mp.nstr((2 * OL / Om) ** (mpf(1) / 3) - 1, 8), '; matter = Lambda at z =', mp.nstr((OL / Om) ** (mpf(1) / 3) - 1, 8))
# jerk j = a''' a^2 / a'^3 = 1 for flat LCDM: check numerically from H(a)
def Hf(x):
    return sqrt(Om * x ** -3 + OL)
x0 = mpf(1); dx = mpf('1e-8')
# q = -1 - a H'/H ; j = q(2q+1) - a dq/da
def qf(x):
    dH = (Hf(x + dx) - Hf(x - dx)) / (2 * dx)
    return -1 - x * dH / Hf(x)
dq = (qf(x0 + dx) - qf(x0 - dx)) / (2 * dx)
jerk = qf(x0) * (2 * qf(x0) + 1) - x0 * dq   # j = q(2q+1) - a dq/da
print('jerk today j0 =', mp.nstr(jerk, 8), '(exactly 1 for flat LCDM)')
print('birefringence 0.264 deg as a fraction of one quarter turn:', mp.nstr(mpf('0.264') / 90, 6))
