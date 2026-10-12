"""The full cube inversion overlap: the cube [-1,1]^n and its image under x -> x/|x|^2.
Image of the interior = exterior of the 'flower' F_n = union of 2n balls of radius 1/2 centred at +-e_i/2
(each face plane x_i = +-1 inverts to the sphere |x -+ e_i/2| = 1/2 through the origin).
Overlap O_n = cube minus flower, since F_n lies inside the unit ball inside the cube.
Volumes: V(region r < g(u)) = (1/n) * integral over the unit sphere of g(u)^n dOmega.
cube: g = 1/m(u), m = max|u_i|;  flower: g = m(u).
"""
import mpmath as mp
mp.mp.dps = 30
from math import factorial

# ---------- 2D exact ----------
# flower area = (1/2) int m^2 dtheta = 8 * (1/2) int_0^{pi/4} cos^2 = 4*(pi/8 + 1/4) = pi/2 + 1
F2 = mp.pi/2 + 1
O2 = 4 - F2
print('2D: flower area =', F2, ' overlap area =', O2, ' overlap share of square =', O2/4, ' flower share =', F2/4)

# ---------- 3D: inclusion-exclusion and direct integral ----------
# direct: (1/3) * int m^3 dOmega ; by symmetry 48 * region x>=y>=z>=0 where m = u_x
# parametrise by the projection onto the face x = 1: u = (1, s, t)/sqrt(1+s^2+t^2), 0<=t<=s<=1
# dOmega = ds dt / (1+s^2+t^2)^(3/2); m^3 = (1+s^2+t^2)^(-3/2)
f3 = lambda s, t: (1 + s*s + t*t)**(-3)
I3 = mp.quad(lambda s: mp.quad(lambda t: f3(s, t), [0, s]), [0, 1])
F3 = 48 * I3 / 3
O3 = 8 - F3
print('3D: flower volume =', F3, ' overlap volume =', O3, ' overlap share =', O3/8, ' flower share =', F3/8)
# inclusion-exclusion check: 6 balls of volume pi/6, 12 adjacent lenses, 8 triple overlaps
r = mp.mpf(1)/2; d = 1/mp.sqrt(2)
lens = mp.pi*(4*r + d)*(2*r - d)**2/12
# triple overlap volume of balls at e_x/2, e_y/2, e_z/2: region |x|^2 <= min(x,y,z)
# by symmetry 6 * region x>=y>=z>=0 with |x|^2 <= z ; in spherical-ish: r <= u_z ... use same face param on z?
# do it via the general formula: (1/3) int min(u_x,u_y,u_z)^3 dOmega over the positive octant
g = lambda s, t: (min(1, s, t)) # placeholder not used
def tri():
    # positive octant, u = (a,b,c) unit; min component cubed. Param by face x=1 with s=y/x,t=z/x in [0,inf);
    # split by symmetry: 6 orderings; take x>=y>=z: min = z component = t/sqrt(1+s^2+t^2)
    h = lambda s, t: (t*t)**1.5 * (1 + s*s + t*t)**(-3)
    return 6 * mp.quad(lambda s: mp.quad(lambda t: h(s, t), [0, s]), [0, 1]) / 3
T3 = tri()
F3b = 6*mp.pi/6 - 12*lens + 8*T3
print('   inclusion-exclusion: 6 balls', mp.pi, ' lens', lens, ' triple', T3, ' -> flower', F3b)

# ---------- 4D ----------
# (1/4) int m^4 dOmega_3 ; 4D face param u = (1,s,t,w)/sqrt(1+s^2+t^2+w^2), region 1>=s>=t>=w>=0, 384 = 2^4*4! copies
# dOmega = ds dt dw / (1+q)^2 ; m^4 = (1+q)^-2  -> integrand (1+q)^-4
f4 = lambda s, t, w: (1 + s*s + t*t + w*w)**(-4)
I4 = mp.quad(lambda s: mp.quad(lambda t: mp.quad(lambda w: f4(s, t, w), [0, t]), [0, s]), [0, 1])
F4 = 384 * I4 / 4
O4 = 16 - F4
print('4D: flower volume =', F4, ' overlap volume =', O4, ' overlap share =', O4/16, ' flower share =', F4/16)
# check the cube itself with the same machinery: (1/4) int m^-4 dOmega = 16
I4c = mp.quad(lambda s: mp.quad(lambda t: mp.quad(lambda w: 1.0*(1+s*s+t*t+w*w)**(-2)*(1+s*s+t*t+w*w)**(2), [0, t]), [0, s]), [0, 1])
print('   check cube volume 4D =', 384*I4c/4)
# unit ball volumes for comparison
for n, Fn in [(2, F2), (3, F3), (4, F4)]:
    Vn = mp.pi**(mp.mpf(n)/2)/mp.gamma(mp.mpf(n)/2 + 1)
    print(f'n={n}: unit ball {mp.nstr(Vn,12)}  flower/ball = {mp.nstr(Fn/Vn,12)}  ball share of cube {mp.nstr(Vn/2**n,12)}')
# identification attempts
for name, x in [('F3', F3), ('O3', O3), ('F4', F4), ('O4', O4), ('T3', T3), ('lens', lens)]:
    print(name, mp.nstr(x, 25), ' identify:', mp.identify(x, ['pi', 'sqrt(2)', 'sqrt(3)', 'pi**2', 'pi*sqrt(2)', 'pi*sqrt(3)', 'atan(sqrt(2))', 'log(2)']))
