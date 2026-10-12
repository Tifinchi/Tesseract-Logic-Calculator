"""Derivatives and inverses of time and scale as tesseract turns, with the cube fixed.

Builds on tesseract.py and ring.py: the cube is the cell w = -1 of the tesseract (+-1)^4,
J: (x, y, z, w) -> (z, -w, -x, y) is the quarter turn fixed by the fabric (q -> i, phi -> v).
"""
import itertools
import cmath
import numpy as np
from mpmath import mp, mpf, sqrt, pi, exp, log

unit = {'+': (0, 0, 1), '−': (1, 0, 1), '×': (1, 1, 1), '*': (0, 1, 1),
        '’': (0, 0, 0), '”': (1, 0, 0), ',': (1, 1, 0), '.': (0, 1, 0)}
fabric = {'”': 'q', '.': 'i', '+': 'φ', '×': 'v', '’': 'L/M', '*': 'L/R', ',': 'RC', '−': 'MC'}
struct = {'”': 'I (scale inverse)', '.': 'D (scale derivative)', '+': 'E (time inverse)', '×': 'T (time derivative)'}
faces = {'Top': '+−×*', 'Front': '+’”−', 'Right': '−”,×', 'Left': '*.’+', 'Back': '×*.,', 'Bottom': '’”,.'}
base = {s: np.array([2 * a - 1 for a in p] + [-1]) for s, p in unit.items()}
V = [np.array(v) for v in itertools.product([-1, 1], repeat=4)]
q, i_, f, v = base['”'], base['.'], base['+'], base['×']
I4 = np.eye(4, dtype=int)

def op(pairs):
    """Linear map from images of the orthogonal basis q, i, phi, v (each |.|^2 = 4)."""
    M = np.zeros((4, 4))
    for src, dst in pairs:
        M += np.outer(dst, src) / 4
    return np.rint(M).astype(int)

J = op([(q, i_), (i_, -q), (f, v), (v, -f)])
K = op([(q, f), (f, -q), (i_, -v), (v, i_)])
L = J @ K
Trev = op([(q, q), (i_, -i_), (f, -f), (v, v)])      # t -> -t: current and flux are odd

def label(u):
    u = tuple(int(x) for x in u)
    for s, b in base.items():
        if u == tuple(b):
            return s
        if u == tuple(-b):
            return '−' + s
    return str(u)

def cell(M):
    img = [M @ b for b in base.values()]
    for k, nm in enumerate('xyzw'):
        col = [u[k] for u in img]
        if len(set(col)) == 1:
            return f'{nm}={col[0]:+d}'
    return '?'

def is_symmetry(M):
    Vs = {tuple(x) for x in V}
    return {tuple(M @ x) for x in V} == Vs

print('=== 1. The second quarter turn K and the quaternion group')
print('K =\n', K)
print('K: (x, y, z, w) ->', tuple(('' if r.sum() > 0 else '-') + 'xyzw'[int(np.nonzero(r)[0][0])] for r in K))
print('L = JK: (x, y, z, w) ->', tuple(('' if r.sum() > 0 else '-') + 'xyzw'[int(np.nonzero(r)[0][0])] for r in L))
print('K maps the tesseract onto itself:', is_symmetry(K), '; L too:', is_symmetry(L))
print('K^2 = -1:', np.array_equal(K @ K, -I4), ' L^2 = -1:', np.array_equal(L @ L, -I4), ' JK = -KJ:', np.array_equal(J @ K, -(K @ J)), ' JKL... JK = L, KL = J, LJ = K:', np.array_equal(K @ L, J), np.array_equal(L @ J, K))
Q8 = {'1': I4, '−1': -I4, 'J': J, '−J': -J, 'K': K, '−K': -K, 'JK': L, '−JK': -L}
cells = {nm: cell(M) for nm, M in Q8.items()}
print('cells reached from the cube:', cells)
print('eight distinct cells (simply transitive):', len(set(cells.values())) == 8)
facename = {'y=+1': 'Back', 'y=-1': 'Front', 'x=-1': 'Left', 'x=+1': 'Right', 'z=+1': 'Top', 'z=-1': 'Bottom'}
print('faces crossed:', {nm: facename.get(c, '—') for nm, c in cells.items() if nm not in ('1', '−1')})

print('\n=== 2. Along which face diagonals each turn moves the tetrahedron corners')
for nm, M in [('J', J), ('K', K), ('JK', L)]:
    moves = {s: label(M @ base[s]) for s in ['”', '.', '+', '×']}
    print(f'  {nm}: ' + ', '.join(f'{s} {struct[s].split()[0]} -> {t}' for s, t in moves.items()))
def diag_face(a, b):
    return [nm for nm, s in faces.items() if a in s and b in s]
for nm, pairs in [('J', [('”', '.'), ('+', '×')]), ('K', [('”', '+'), ('.', '×')]), ('JK', [('”', '×'), ('.', '+')])]:
    print(f'  {nm} runs along the diagonals of', [diag_face(a, b)[0] for a, b in pairs])

print('\n=== 3. Strings carried by K and JK (J carries Front onto Back as written, ring.py)')
def same_cycle(a, b):
    return any(a == b[k:] + b[:k] for k in range(4))
for nm, M in [('J', J), ('K', K), ('JK', L)]:
    for fn, s in faces.items():
        img = [label(M @ base[ch]) for ch in s]
        if all(len(x) == 1 for x in img):
            t = ''.join(img)
            hit = [g for g, s2 in faces.items() if same_cycle(t, s2)]
            hitr = [g for g, s2 in faces.items() if same_cycle(t[::-1], s2)]
            exact = [g for g, s2 in faces.items() if t == s2]
            print(f'  {nm}({fn} {s}) = {t}: ' + (f'{exact[0]} as written, same start' if exact else f'{hit[0]} as written (cyclic shift)' if hit else f'{hitr[0]} reversed' if hitr else 'no face'))

print('\n=== 4. Time reversal t -> -t: current and flux turn over')
print('Trev: (x, y, z, w) ->', tuple(('' if r.sum() > 0 else '-') + 'xyzw'[int(np.nonzero(r)[0][0])] for r in Trev))
print('a symmetry of the tesseract:', is_symmetry(Trev), '; det', round(np.linalg.det(Trev)), '; Trev^2 = 1:', np.array_equal(Trev @ Trev, I4))
print('cell reached:', cell(Trev), '=', facename[cell(Trev)], 'face;', 'corners kept:', [s for s in base if label(Trev @ base[s]) == s], '; swapped:', [(s, label(Trev @ base[s])) for s in base if label(Trev @ base[s]) not in (s,) and len(label(Trev @ base[s])) == 1])
print('Trev J Trev = J^-1 (time reversal turns the derivative into the anti-derivative):', np.array_equal(Trev @ J @ Trev, -J))

print('\n=== 5. The six two-variable turn-overs and the cells they reach')
names = ['q', 'i', 'φ', 'v']
vecs = [q, i_, f, v]
for a_, b_ in itertools.combinations(range(4), 2):
    pairs = [(vecs[k], -vecs[k] if k in (a_, b_) else vecs[k]) for k in range(4)]
    M = op(pairs)
    c = cell(M)
    print(f'  turn over ({names[a_]}, {names[b_]}): symmetry {is_symmetry(M)}, cell {c} ({facename[c]} face kept: {sorted(s for s in base if label(M @ base[s]) in base and unit[s][{"x": 0, "y": 1, "z": 2}[c[0]]] == (1 if c[2] == "+" else 0))})')

print('\n=== 6. Scale and phase in one plane: z = r e^{i theta}, ln z = ln r + i theta')
for z in [2.5 * cmath.exp(0.7j), 0.4 * cmath.exp(-1.2j)]:
    zi = cmath.exp(1j * cmath.log(z))           # z^i, principal branch
    lnr, th = cmath.log(abs(z)).real, cmath.phase(z)
    print(f'  z = {z:.4f}: z^i has log-size {cmath.log(abs(zi)).real:+.6f} (= -theta {-th:+.6f}) and phase {cmath.phase(zi):+.6f} (= ln r {lnr:+.6f})')
    zii = cmath.exp(1j * cmath.log(zi))
    print(f'     (z^i)^i = {zii:.6f}; 1/z = {1/z:.6f}')
# Cauchy-Riemann in polar form: d/dtheta f = i r d/dr f for analytic f
fz = lambda w: w ** 3 + cmath.exp(w) / (1 + w * w)
r0, t0, h = 1.3, 0.8, 1e-6
d_theta = (fz(r0 * cmath.exp(1j * (t0 + h))) - fz(r0 * cmath.exp(1j * (t0 - h)))) / (2 * h)
d_r = (fz((r0 + h) * cmath.exp(1j * t0)) - fz((r0 - h) * cmath.exp(1j * t0))) / (2 * h)
print(f'  Cauchy-Riemann: d f/d theta = {d_theta:.6f}; i r d f/d r = {1j * r0 * d_r:.6f}')

print('\n=== 7. The electron: one charge stepped in time (omega) and in scale (R0)')
mp.dps = 30
h_ = mpf('6.62607015e-34'); e = mpf('1.602176634e-19'); c = mpf(299792458)
hbar = h_ / (2 * pi); me = mpf('9.1093837139e-31'); alpha = mpf('7.2973525643e-3')
eps0 = e ** 2 / (2 * alpha * h_ * c); mu0 = 1 / (eps0 * c ** 2); Z0 = mu0 * c
lam = hbar / (me * c); re = alpha * lam; a0 = lam / alpha
Ce = 4 * pi * eps0 * re; Le = pi * mu0 * a0
w = 1 / sqrt(Le * Ce); R0 = sqrt(Le / Ce)
print('omega =', mp.nstr(w, 12), '; c/(2 pi sqrt(re a0)) =', mp.nstr(c / (2 * pi * sqrt(re * a0)), 12))
print('R0 =', mp.nstr(R0, 12), '; h/e^2 =', mp.nstr(h_ / e ** 2, 12), '; (Z0/2) sqrt(a0/re) =', mp.nstr(Z0 / 2 * sqrt(a0 / re), 12), '; Z0/(2 alpha) =', mp.nstr(Z0 / (2 * alpha), 12))
print('corners e, e w, e R0, e w R0 =', [mp.nstr(x, 10) for x in (e, e * w, e * R0, e * w * R0)])
print('   against q, i = w q, phi = L i, v = q/C:', [mp.nstr(x, 10) for x in (e, w * e, Le * w * e, e / Ce)])
ZC = lambda r: 1 / (w * 4 * pi * eps0 * r)
ZL = lambda r: w * pi * mu0 * r
for r in (re, lam, a0, mpf('1e-14')):
    print(f'  r = {mp.nstr(r, 6)} m: Z_C(r) = {mp.nstr(ZC(r), 10)} ohm; Z_L(lam^2/r) = {mp.nstr(ZL(lam ** 2 / r), 10)} ohm; (Z0/2)(lam/r) = {mp.nstr(Z0 / 2 * lam / r, 10)}')
print('the flip re <-> a0 keeps sqrt(re a0) (the time scale) and inverts a0/re (the scale step):', mp.nstr(sqrt(re * a0) / lam, 15), mp.nstr((a0 / re) * alpha ** 2, 15))

print('\n=== 8. Gravity: the scale flip through the horizon, isotropic coordinates')
M_ = mpf(1)
def areal(rho):
    return rho * (1 + M_ / (2 * rho)) ** 2
def lapse(rho):
    return (1 - M_ / (2 * rho)) / (1 + M_ / (2 * rho))
def psi4(rho):
    return (1 + M_ / (2 * rho)) ** 4
for rho in (mpf('0.2'), mpf('0.5'), mpf('3')):
    rp = M_ ** 2 / (4 * rho)
    drp = M_ ** 2 / (4 * rho ** 2)                   # |d rho'/d rho|
    print(f'  rho = {mp.nstr(rho, 4)} -> {mp.nstr(rp, 6)}: areal r {mp.nstr(areal(rho), 10)} vs {mp.nstr(areal(rp), 10)}; '
          f'radial metric {mp.nstr(psi4(rho), 10)} vs {mp.nstr(psi4(rp) * drp ** 2, 10)}; angular {mp.nstr(psi4(rho) * rho ** 2, 10)} vs {mp.nstr(psi4(rp) * rp ** 2, 10)}; lapse {mp.nstr(lapse(rho), 8)} vs {mp.nstr(lapse(rp), 8)}')
print('  fixed sphere rho = M/2:', M_ ** 2 / (4 * (M_ / 2)) == M_ / 2, '; its areal radius', mp.nstr(areal(M_ / 2), 6), '= 2M (the horizon)')

print('\n=== 9. Efimov: one scale step is a half turn of the log-phase')
s0 = mpf('1.00624')
lam0 = exp(pi / s0)
print('e^(pi/s0) =', mp.nstr(lam0, 8), '; phase per step s0 ln(lambda0) =', mp.nstr(s0 * log(lam0), 12), '= pi')
print('Huang et al. 2014: 21.0 +- 1.3 against', mp.nstr(lam0, 5), '->', mp.nstr((mpf('21.0') - lam0) / mpf('1.3'), 4), 'sigma')
print('energy ratio per step lambda0^2 =', mp.nstr(lam0 ** 2, 6), '(a full turn of phase in energy)')

print('\n=== 10. Electron EDM (time-reversal test)')
de = mpf('-1.3'); sde = sqrt(mpf('2.0') ** 2 + mpf('0.6') ** 2)
print('d_e = (-1.3 +- 2.0 +- 0.6)e-30 e cm: combined sigma', mp.nstr(sde, 4), '; pull from 0:', mp.nstr(de / sde, 4))

print('\n=== 11. A quarter turn of one circle alone (time only or scale only)')
axis = [2 * s * I4[k] for k in range(4) for s in (1, -1)]
odd = [x for x in V if np.sum(x < 0) % 2 == 1]
even = [x for x in V if np.sum(x < 0) % 2 == 0]
setof = lambda S: {tuple(int(t) for t in x) for x in S}
for nm, pairs in [('time circle only (+ -> ×)', [(q, q), (i_, i_), (f, v), (v, -f)]), ('scale circle only (” -> .)', [(q, i_), (i_, -q), (f, f), (v, v)])]:
    M = np.zeros((4, 4))
    for src, dst in pairs:
        M += np.outer(dst, src) / 4
    print(f'  {nm}: tesseract symmetry {setof(M @ x for x in V) == setof(V)}; variables kept {setof(M @ x for x in odd) == setof(odd)}; circuits -> cell centres {setof(M @ x for x in even) == setof(axis)}; cell centres -> circuits {setof(M @ p for p in axis) == setof(even)}')
