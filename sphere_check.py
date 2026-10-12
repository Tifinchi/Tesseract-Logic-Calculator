"""The sphere in the cube: levels, shares, the 4D torus, and the layers' fixed spheres."""
import itertools
from mpmath import mp, mpf, pi, sqrt
mp.dps = 30
# levels of the side-2 cube centred on 0
pts = list(itertools.product([-1, 0, 1], repeat=3))
levels = sorted({sum(c*c for c in p) for p in pts})
print("squared levels in 3D:", levels, "-> radii 0, 1, sqrt2, sqrt3")
print("counts:", {k: sum(1 for p in pts if sum(c*c for c in p) == k) for k in levels})
# flip r -> 1/r on the levels
for r in (1, sqrt(2), sqrt(3), 2):
    print(f"  r={mp.nstr(r,6)} -> 1/r={mp.nstr(1/r,6)}  product={mp.nstr(r*(1/r),3)}")
# ball shares
print("3D ball share pi/6 =", mp.nstr(pi/6, 8))
print("4D ball share (pi^2/2)/16 = pi^2/32 =", mp.nstr(pi**2/32, 8), " corners:", mp.nstr(1 - pi**2/32, 8))
print("8 pi^2 =", mp.nstr(8*pi**2, 8), " = 4pi*2pi:", mp.nstr(4*pi*2*pi, 8))
# 4D: tesseract corners on the Clifford torus of two sqrt2 circles
corners = list(itertools.product([-1, 1], repeat=4))
on_torus = all(c[0]**2 + c[1]**2 == 2 and c[2]**2 + c[3]**2 == 2 for c in corners)
print("16 tesseract corners on the torus |z1|=|z2|=sqrt2:", on_torus, " radius:", mp.nstr(sqrt(4), 3))
axis = [p for p in itertools.product([-2, 0, 2], repeat=4) if sum(abs(x) for x in p) == 2]
print("8 axis points:", len(axis), "on core circles:", all((p[0]**2+p[1]**2, p[2]**2+p[3]**2) in [(4, 0), (0, 4)] for p in axis))
# T-bar flip: swap the two planes
flip = lambda c: (c[2], c[3], c[0], c[1])
print("plane swap maps the 16 corners onto themselves:", sorted(map(flip, corners)) == sorted(corners))
# electron's Compton sphere and the horizon
alpha = mpf('7.2973525643e-3'); h = mpf('6.62607015e-34'); me = mpf('9.1093837139e-31'); c = mpf(299792458)
lam = h/(2*pi*me*c); re = alpha*lam; a0 = lam/alpha
print("lambda-bar (fm):", mp.nstr(lam*1e15, 8), " r_e inside:", re < lam, " a0 outside:", a0 > lam, " r_e*a0/lam^2 =", mp.nstr(re*a0/lam**2, 12))
print("flip r -> lam^2/r sends r_e to", mp.nstr(lam**2/re*1e12, 8), "pm = a0", mp.nstr(a0*1e12, 8), "pm")
# face centres = midpoints of the tetrahedron edges (unit cube coords)
tet = {'”': (1,0,0), '.': (0,1,0), '+': (0,0,1), '×': (1,1,1)}
mids = sorted({tuple((a+b)/2 for a, b in zip(tet[p], tet[q])) for p, q in itertools.combinations(tet, 2)})
print("tetrahedron edge midpoints:", mids)
