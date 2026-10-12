"""Checks the group-theory and string claims of section 3 of The Cube Cipher."""
import itertools
import numpy as np

corner = {'+': (0, 0, 1), '−': (1, 0, 1), '×': (1, 1, 1), '*': (0, 1, 1),
          '’': (0, 0, 0), '”': (1, 0, 0), ',': (1, 1, 0), '.': (0, 1, 0)}
pos = {v: k for k, v in corner.items()}
faces = {'Top': '+−×*', 'Front': '+’”−', 'Right': '−”,×', 'Left': '*.’+',
         'Back': '×*.,', 'Bottom': '’”,.'}
outward = {'Top': (0, 0, 1), 'Bottom': (0, 0, -1), 'Front': (0, -1, 0), 'Back': (0, 1, 0),
           'Right': (1, 0, 0), 'Left': (-1, 0, 0)}
tet = set('”.+×'); dual = set('’−*,')

def c(s):  # centred coords
    return np.array([2 * x - 1 for x in corner[s]])

# 1. orientation of each string (normal by right-hand rule)
print('1. string normals')
for f, s in faces.items():
    p = [c(ch) for ch in s]
    n = np.cross(p[1] - p[0], p[2] - p[1])
    n = n / np.linalg.norm(n)
    o = np.array(outward[f])
    print(f'  {f:6s} {s}  normal {n.astype(int)}  {"outward" if n @ o > 0 else "inward"}')

# 2. all 48 matrices of O_h
mats = []
for perm in itertools.permutations(range(3)):
    for signs in itertools.product([1, -1], repeat=3):
        M = np.zeros((3, 3), int)
        for i, j in enumerate(perm):
            M[i, j] = signs[i]
        mats.append(M)
assert len(mats) == 48
def act(M, s):
    v = M @ c(s)
    return pos[tuple(int((x + 1) // 2) for x in v)]
rot = [M for M in mats if round(np.linalg.det(M)) == 1]
print('2. rotations:', len(rot))
keep = [M for M in rot if {act(M, s) for s in tet} == tet]
swap = [M for M in rot if {act(M, s) for s in tet} == dual]
print('   keep tetrahedron:', len(keep), ' swap with dual:', len(swap))
inv = -np.eye(3, dtype=int)
print('   inversion sends tetrahedron to', ''.join(sorted(act(inv, s) for s in tet)), '(dual is', ''.join(sorted(dual)), ')')
# classify rotations by angle
def angle(M):
    return round(np.degrees(np.arccos(max(-1, min(1, (np.trace(M) - 1) / 2)))))
from collections import Counter
print('   rotation angles:', Counter(angle(M) for M in rot))
print('   keepers by angle:', Counter(angle(M) for M in keep))
# face half-turns pair the tetrahedron corners
for name, M in [('z', np.diag([-1, -1, 1])), ('x', np.diag([1, -1, -1])), ('y', np.diag([-1, 1, -1]))]:
    pairs = sorted({tuple(sorted((s, act(M, s)))) for s in tet})
    print(f'   half-turn {name}: tetra pairs {pairs}')

# 3. i-turn cycles vs strings
def cycle_forward(M, s):
    # does M advance each element of string s one step forward?
    return all(act(M, s[k]) == s[(k + 1) % 4] for k in range(4))
def cycle_backward(M, s):
    return all(act(M, s[k]) == s[(k - 1) % 4] for k in range(4))
iz = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])   # (a,b,d) -> (-b,a,d)
ix = np.array([[1, 0, 0], [0, 0, 1], [0, -1, 0]])   # (a,b,d) -> (a,d,-b)
iy = np.array([[0, 0, -1], [0, 1, 0], [1, 0, 0]])   # (a,b,d) -> (-d,b,a)
print('3. strings as i-turn cycles')
for name, M, pair in [('i_z', iz, ('Top', 'Bottom')), ('i_y', iy, ('Front', 'Back')), ('i_x', ix, ('Left', 'Right'))]:
    for f in pair:
        s = faces[f]
        print(f'   {name} on {f:6s}: forward={cycle_forward(M, s)} backward={cycle_backward(M, s)}')

# 4. inversion of strings
print('4. inversion of strings')
for f in ('Top', 'Front', 'Right'):
    s = faces[f]
    t = ''.join(act(inv, ch) for ch in s)
    def same_cycle(a, b):
        return any(a == b[k:] + b[:k] for k in range(4))
    for g, s2 in faces.items():
        if same_cycle(t, s2):
            print(f'   {f} -> {t} = {g} as written')
        elif same_cycle(t[::-1], s2):
            print(f'   {f} -> {t} = {g} reversed')

# 5. omega-turn about the zero axis and the dual-corner circuits
om = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])  # (a,b,d) -> (d,a,b)
print('5. omega-turn:', {s: act(om, s) for s in corner})
neigh = {s: ''.join(t for t in corner if sum(abs(np.array(corner[s]) - np.array(corner[t]))) == 1) for s in corner}
var = {'”': 'q', '.': 'i', '+': 'phi', '×': 'v'}
for d in sorted(dual):
    print(f'   dual corner {d}: neighbours {neigh[d]} -> variables {[var[t] for t in neigh[d]]}; opposite {act(inv, d)} = {var[act(inv, d)]}')
