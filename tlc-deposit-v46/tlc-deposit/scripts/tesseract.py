"""The corrected cipher mechanics, checked exactly.

The cube never moves: it is the base cell (w = -1) of a tesseract with vertices (+-1)^4.
A derivative is a quarter turn J in 4D (J^2 = -I); an inverse is the half turn J^2 = -I.
J is fixed by the fabric layer: it must send the stores to the rates, charge -> current
and flux -> voltage (i = dq/dt, v = dphi/dt).
"""
import itertools
import numpy as np

unit = {'+': (0, 0, 1), '−': (1, 0, 1), '×': (1, 1, 1), '*': (0, 1, 1),
        '’': (0, 0, 0), '”': (1, 0, 0), ',': (1, 1, 0), '.': (0, 1, 0)}
fabric = {'”': 'q', '.': 'i', '+': 'φ', '×': 'v', '’': 'L/M', '*': 'L/R', ',': 'RC', '−': 'MC'}
base = {s: np.array([2 * a - 1 for a in p] + [-1]) for s, p in unit.items()}
V = [np.array(v) for v in itertools.product([-1, 1], repeat=4)]

def name(v):
    """Name a tesseract vertex by its cell-position delimiter and its cell."""
    v = tuple(int(x) for x in v)
    for s, b in base.items():
        if v == tuple(b):
            return s + ' (base)'
        if v == tuple(-b):
            return 'inverse ' + s
    # other vertices: base position with w = +1
    for s, b in base.items():
        if v[:3] == tuple(b[:3]) and v[3] == 1:
            return s + ' position, w=+1'
    return str(v)

# --- find J from the fabric constraints: q->i, i->-q, phi->v, v->-phi
a, b, c, d = base['”'], base['.'], base['+'], base['×']
print('mutual dot products of q, i, phi, v vectors:', [int(x @ y) for x, y in itertools.combinations([a, b, c, d], 2)])
J = (np.outer(b, a) - np.outer(a, b) + np.outer(d, c) - np.outer(c, d)) / 4
print('J =\n', J.astype(int))
print('J is a signed permutation:', np.all(np.isin(J, [-1, 0, 1])) and np.all(np.abs(J).sum(0) == 1))
print('J^2 = -I:', np.array_equal(J @ J, -np.eye(4)), ' J^4 = I:', np.array_equal(np.linalg.matrix_power(J, 4), np.eye(4)), ' det J =', round(np.linalg.det(J)))
print('J as a map: (x, y, z, w) ->', tuple('xyzw'[int(np.nonzero(row)[0][0])] if row.sum() > 0 else '-' + 'xyzw'[int(np.nonzero(row)[0][0])] for row in J.astype(int)))

# vertices to vertices, edges to edges
Vs = {tuple(v) for v in V}
print('J maps the 16 vertices onto themselves:', {tuple((J @ v).astype(int)) for v in V} == Vs)
edges = {frozenset((tuple(u), tuple(w))) for u in V for w in V if np.sum(u != w) == 1}
print('edges:', len(edges), ' J keeps every edge an edge:', all(frozenset((tuple((J @ np.array(e1)).astype(int)), tuple((J @ np.array(e2)).astype(int)))) in edges for e1, e2 in map(tuple, edges)))

# where each base delimiter goes under J, J^2, J^3
print('\nfour layers for each base corner (instant -> rate -> inverse -> anti-rate):')
for s in ['”', '.', '+', '×', '’', '*', ',', '−']:
    v = base[s]
    orbit = [v, J @ v, J @ J @ v, J @ J @ J @ v]
    labels = []
    for u in orbit:
        u = tuple(int(x) for x in u)
        hit = [t for t, bb in base.items() if tuple(bb) == u]
        hitn = [t for t, bb in base.items() if tuple(-bb) == u]
        labels.append(f'{hit[0]} {fabric[hit[0]]}' if hit else (f'inverse {hitn[0]} (−{fabric[hitn[0]]})' if hitn else str(u)))
    print(f'  {s} {fabric[s]:4s}: ' + '  ->  '.join(labels))

# cells: which cell does J send the base cell to?
def cell_of(vs):
    vs = np.array(vs)
    out = []
    for k in range(4):
        col = vs[:, k]
        if np.all(col == col[0]):
            out.append(f"{'xyzw'[k]}={int(col[0]):+d}")
    return out
img = [J @ v for v in base.values()]
print('\nJ sends the base cell (w=-1) to the cell', cell_of(img))
print('J^2 sends it to', cell_of([J @ J @ v for v in base.values()]), ' J^3 to', cell_of([J @ J @ J @ v for v in base.values()]))
shared = [s for s, v in base.items() if v[1] == 1]
print('base corners on the face shared with the rate layer (y=+1):', shared, '= the Back face' if set(shared) == set('×*.,') else '')

# circles: each vertex's orbit under R(theta) = cos(theta) I + sin(theta) J is a great circle of radius 2
th = np.linspace(0, 2 * np.pi, 9)
r = [np.linalg.norm((np.cos(t) * np.eye(4) + np.sin(t) * J) @ a) for t in th]
print('\nradius along the charge corner\'s circle:', np.round(r, 12))
fibres = {frozenset(tuple(int(x) for x in np.linalg.matrix_power(J, k) @ v) for k in range(4)) for v in V}
print('circles through the 16 vertices:', len(fibres), 'each with', {len(f) for f in fibres}, 'vertices')
# Clifford torus: |(x,z)| = |(y,w)| = sqrt2 for every vertex (J rotates the (x,z) and (y,w) planes)
print('every vertex on the torus |(x,z)| = |(y,w)| = sqrt2:', all(abs(np.hypot(v[0], v[2]) - np.sqrt(2)) < 1e-12 and abs(np.hypot(v[1], v[3]) - np.sqrt(2)) < 1e-12 for v in V))

# parity classes (two 16-cells) and the 24-cell
odd = [v for v in V if np.sum(v < 0) % 2 == 1]
even = [v for v in V if np.sum(v < 0) % 2 == 0]
tet_parity = {s: ('odd' if np.sum(base[s] < 0) % 2 else 'even') for s in base}
print('\nparity of the base corners:', tet_parity)
print('J keeps parity classes:', all(np.sum((J @ v) < 0) % 2 == np.sum(v < 0) % 2 for v in V))
axis = [2 * s * np.eye(4)[k] for k in range(4) for s in (1, -1)]
print('8 axis points at radius', {round(float(np.linalg.norm(p)), 12) for p in axis}, '; 16 corners at radius', {round(float(np.linalg.norm(v)), 12) for v in V})
pts = [tuple(v) for v in V] + [tuple(p) for p in axis]
d2 = sorted({round(float(np.sum((np.array(p) - np.array(q)) ** 2)), 9) for p, q in itertools.combinations(pts, 2)})
nn = [sum(1 for q in pts if round(float(np.sum((np.array(p) - np.array(q)) ** 2)), 9) == d2[0]) for p in pts]
print('24 points: smallest squared distance', d2[0], '; neighbours per point', set(nn), '(24-cell: 8 each)')

# Kruskal: the half turn (T,X)->(-T,-X) keeps X^2 - T^2; the quarter turn (T,X)->(X,-T) negates it
import random
random.seed(1)
ok_half = ok_quarter = True
for _ in range(1000):
    T, X = random.uniform(-3, 3), random.uniform(-3, 3)
    if abs((X * X - T * T) - ((-X) ** 2 - (-T) ** 2)) > 1e-12: ok_half = False
    if abs((X * X - T * T) + ((-T) ** 2 - X ** 2)) > 1e-12: ok_quarter = False
print('\nKruskal: half turn keeps X^2-T^2 (an isometry):', ok_half, '; quarter turn negates it (time and space swap):', ok_quarter)
def region(T, X):
    return 'I' if X > abs(T) else 'II' if T > abs(X) else 'III' if X < -abs(T) else 'IV'
T, X = 0.3, 1.0
seq = [region(T, X)]
for _ in range(4):
    T, X = X, -T
    seq.append(region(T, X))
print('quarter turns carry region', ' -> '.join(seq))
