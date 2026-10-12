# The electron as the fabric's reaction: does it supply the spin's half? (section 4.6)
from mpmath import mp, mpf, pi, sqrt, quad, sin, inf, exp, mpc
mp.dps = 30
h = mpf('6.62607015e-34'); e = mpf('1.602176634e-19'); c = mpf(299792458)
alpha = mpf('7.2973525643e-3'); me = mpf('9.1093837139e-31')      # CODATA 2022
hbar = h/(2*pi); mu0 = 2*alpha*h/(e**2*c); eps0 = 1/(mu0*c**2)
mec2 = me*c**2
f = mec2/h                                   # the ring's circuits per second
lam = hbar/(me*c); re = alpha*lam; a0 = lam/alpha
muB = e*hbar/(2*me)
q, i, phi, v = e, e*f, h/e, mec2/e           # the four corners
Ce, Le = q/v, phi/i
keV = 1000*e

print("=== 1. The two stores")
UC, UL = q*v/2, phi*i/2
print("  held (capacitor)  qv/2   = %s mec2 = %s keV" % (mp.nstr(UC/mec2, 20), mp.nstr(UC/keV, 9)))
print("  moving (inductor) phi*i/2 = %s mec2" % mp.nstr(UL/mec2, 20))
print("  sum = %s mec2" % mp.nstr((UC+UL)/mec2, 20))
print("  charge e on a sphere of radius re: e^2/(8 pi eps0 re) = %s mec2" % mp.nstr(e**2/(8*pi*eps0*re)/mec2, 20))
print("  check Ce = 4 pi eps0 re: %s; Le = pi mu0 a0: %s" % (mp.nstr(Ce/(4*pi*eps0*re), 20), mp.nstr(Le/(pi*mu0*a0), 20)))

print("=== 2. Angular momentum round the Compton circle (radius lam, speed c)")
L_all = (mec2/c)*lam
L_mov = (UL/c)*lam
mu = i*pi*lam**2
print("  all of mec2 going round: %s hbar" % mp.nstr(L_all/hbar, 20))
print("  only the moving half:    %s hbar" % mp.nstr(L_mov/hbar, 20))
print("  moment of the full charge, i*pi*lam^2 = %s muB" % mp.nstr(mu/muB, 20))
g = (mu/muB)/(L_mov/hbar)
print("  g = (mu/muB)/(J/hbar) = %s" % mp.nstr(g, 20))
print("  phase after one circuit for J = hbar/2: %s" % mp.nstr(exp(mpc(0, -2*pi*L_mov/hbar)), 6))
print("  Dirac-field split instead (all the mass at half the charge's rate): m (Omega/2) lam^2 = %s hbar" % mp.nstr(me*(c/lam/2)*lam**2/hbar, 20))

print("=== 3. Classical field angular momentum: charge e on a sphere of radius R, moment muB at its centre")
# l = r x (eps0 E x B); E radial outside R only; B a point dipole along z. z-part integrand:
# eps0 * q mu0 m /(16 pi^2 eps0 r^4) * sin^2(theta), over r^2 dr dOmega
def Lfield(R, m=muB, Q=e):
    ang = quad(lambda t: sin(t)**2*sin(t), [0, pi])*2*pi
    rad = quad(lambda r: 1/r**2, [R, inf])
    return Q*mu0*m/(16*pi**2)*ang*rad
print("  angular integral = %s (8 pi/3 = %s)" % (mp.nstr(quad(lambda t: sin(t)**3, [0, pi])*2*pi, 15), mp.nstr(8*pi/3, 15)))
print("  numerical L(re) = %s hbar; analytic mu0 e muB/(6 pi re) = %s hbar" % (mp.nstr(Lfield(re)/hbar, 15), mp.nstr(mu0*e*muB/(6*pi*re)/hbar, 15)))
print("  radius giving hbar/2: %s re" % mp.nstr(mu0*e*muB/(6*pi*(hbar/2))/re, 15))
M = mpf(1); R = mpf(1); Q = mpf(1)   # Griffiths-type check, uniform magnetisation M: m = 4/3 pi R^3 M -> (2/9) mu0 M Q R^2
print("  check, magnetised charged sphere: L/(mu0 M Q R^2) = %s (2/9 = %s)" % (mp.nstr(Lfield(R, m=4*pi*R**3*M/3, Q=Q)/(mu0*M*Q*R**2), 12), mp.nstr(mpf(2)/9, 12)))

print("=== 4. The ring's rate and the tank's rate")
Tc = h/mec2
print("  sqrt(Le Ce) = %s s; Compton period h/mec2 = %s s; ratio %s" % (mp.nstr(sqrt(Le*Ce), 10), mp.nstr(Tc, 10), mp.nstr(sqrt(Le*Ce)/Tc, 20)))
print("  ring: %s circuits per second; a standard LC tank with Le, Ce: %s cycles per second; ratio %s (2 pi = %s)" % (mp.nstr(f, 8), mp.nstr(1/(2*pi*sqrt(Le*Ce)), 8), mp.nstr(f*2*pi*sqrt(Le*Ce), 15), mp.nstr(2*pi, 15)))
print("  moving energy over the ring's angular rate 2 pi f: %s hbar; over the tank's 1/sqrt(LC) in rad/s: %s hbar" % (mp.nstr(UL/(2*pi*f)/hbar, 15), mp.nstr(UL*sqrt(Le*Ce)/hbar, 15)))
print("  impedances at 1/sqrt(LC): Le/sqrt(LeCe) = %s ohm, sqrt(LeCe)/Ce = %s ohm, h/e^2 = %s ohm" % (mp.nstr(Le/sqrt(Le*Ce), 12), mp.nstr(sqrt(Le*Ce)/Ce, 12), mp.nstr(h/e**2, 12)))

print("=== 5. The charge-flux product over a circle and over a sphere")
print("  q*phi/(2 pi) = %s hbar; q*phi/(4 pi) = %s hbar" % (mp.nstr(q*phi/(2*pi)/hbar, 20), mp.nstr(q*phi/(4*pi)/hbar, 20)))

print("=== 6. The excess over g = 2")
ae = mpf('0.00115965218059')
a1 = alpha/(2*pi); x = alpha/pi
a2 = a1 + mpf('-0.328478965579193')*x**2; a3 = a2 + mpf('1.181241456587')*x**3
print("  measured ae = %s; alpha/2pi = %s; + two loops = %s; + three loops = %s" % (mp.nstr(ae, 12), mp.nstr(a1, 12), mp.nstr(a2, 12), mp.nstr(a3, 12)))
