# The calculator's Version 43 numerics, checked independently: the exact law-of-motion engine (motion_exact.js) against
# 50-digit quadrature of the original integrals, and the off-shell Higgs widths (vv_fast.js) against an adaptive double
# integral of the original integrand. The JavaScript values are produced by Node and compared here.
import json, subprocess, math
import mpmath as mp
from scipy import integrate
D = '/tmp/claude-0/-home-claude/bbcd54b6-3962-59fe-98a5-d2cef4b7e16c/scratchpad/live/v43/'
js = r"""
const X = require('%smotion_exact.js'), V = require('%svv_fast.js');
const r = {
  merc_both: X.advance(0.387098, 0.205630, 'both'), merc_time: X.advance(0.387098, 0.205630, 'time'),
  s200_time: X.advance(200 * 1476.6250 / 1.495978707e11, 0.3, 'time'), s200_both: X.advance(200 * 1476.6250 / 1.495978707e11, 0.3, 'both'),
  e95_time: X.advance(0.1, 0.95, 'time'), e95_both: X.advance(0.1, 0.95, 'both'),
  b1_time: X.bending(1, 1, 'time'), b1_both: X.bending(1, 1, 'both'), b10_both: X.bending(10, 1, 'both'),
  b9e4_both: X.bending(1, 9e4, 'both'), b1e5_time: X.bending(1, 1e5, 'time'),
  WW: V.VVfast(125.13, 80.3692, 2.085, 2, 1.1663788e-5, 4), ZZ: V.VVfast(125.13, 91.1880, 2.4955, 1, 1.1663788e-5, 4) };
console.log(JSON.stringify(Object.fromEntries(Object.entries(r).map(([k, v]) => [k, v.toPrecision(17)]))));
""" % (D, D)
J = json.loads(subprocess.run(['node', '-e', js], capture_output=True, text=True, check=True).stdout)
mp.mp.dps = 60
C = mp.mpf(299792458); GMS = mp.mpf('1.32712440018e20'); AU = mp.mpf('1.495978707e11'); RSUN = mp.mpf('6.957e8'); AS = mp.mpf('206264.80624709636')
def adv_both(aAU, e):
    a = mp.mpf(aAU)*AU/(GMS/C**2); e = mp.mpf(e); ua, up = 1/(a*(1 + e)), 1/(a*(1 - e)); u3 = mp.mpf(1)/2 - ua - up
    m, h = (ua + up)/2, (up - ua)/2      # u = m - h cos t removes the endpoint singularities
    return 2*mp.quad(lambda t: 1/mp.sqrt(2*(u3 - (m - h*mp.cos(t)))), [0, mp.pi]) - 2*mp.pi
def adv_time(aAU, e):
    a = mp.mpf(aAU)*AU/(GMS/C**2); e = mp.mpf(e); rp, ra = a*(1 - e), a*(1 + e)
    L2 = (2/ra - 2/rp)/((1 - 2/ra)/ra**2 - (1 - 2/rp)/rp**2); L = mp.sqrt(L2)
    w0 = 2/(rp - 1 + mp.sqrt(rp*rp - 2*rp)); G = lambda w: ((2 + w)/(2 - w))**2; E2 = (1 + L2*w0**2)/G(w0)
    pr2 = lambda w: E2*G(w) - 1 - L2*w**2
    w1 = mp.findroot(pr2, 2/L2 - w0, tol=mp.mpf(10)**-55); m, h = (w0 + w1)/2, (w0 - w1)/2
    return 2*mp.quad(lambda t: L*h*mp.sin(t)/mp.sqrt(mp.re(pr2(m - h*mp.cos(t)))) if mp.sin(t) != 0 else 0, [0, mp.pi/2, mp.pi]) - 2*mp.pi
def bend(bR, M, kind):
    b = mp.mpf(bR)*RSUN/(GMS*mp.mpf(M)/C**2)
    n = (lambda w: (2 + w)/(2 - w)) if kind == 'time' else (lambda w: (2 + w)**3/(4*(2 - w)))
    w0 = mp.findroot(lambda w: b*w - n(w), 1/b)
    I = mp.quad(lambda w: 1/mp.sqrt(n(w)**2/b**2 - w**2), [0, w0/2, w0])
    return (2*I - mp.pi)*AS
REF = {'merc_both': adv_both('0.387098', '0.205630'), 'merc_time': adv_time('0.387098', '0.205630'),
       's200_time': adv_time(200*mp.mpf('1476.6250')/AU, '0.3'), 's200_both': adv_both(200*mp.mpf('1476.6250')/AU, '0.3'),
       'e95_time': adv_time('0.1', '0.95'), 'e95_both': adv_both('0.1', '0.95'),
       'b1_time': bend(1, 1, 'time'), 'b1_both': bend(1, 1, 'both'), 'b10_both': bend(10, 1, 'both'),
       'b9e4_both': bend(1, '9e4', 'both'), 'b1e5_time': bend(1, '1e5', 'time')}
GF = 1.1663788e-5; mH = 125.13
def VV(mV, GV, delta):
    bw = lambda s: (mV*GV/math.pi)/((s - mV*mV)**2 + (mV*GV)**2)
    def f(q2, q1):
        x = q1*q1/mH**2; y = q2*q2/mH**2; lam = max(0.0, (1 - x - y)**2 - 4*x*y)
        return delta*GF*mH**3/(16*math.sqrt(2)*math.pi)*math.sqrt(lam)*(lam + 12*x*y)*bw(q1*q1)*2*q1*bw(q2*q2)*2*q2
    return integrate.dblquad(f, 0, mH, lambda q1: 0, lambda q1: mH - q1, epsabs=1e-16, epsrel=1e-12)[0]
REF['WW'] = mp.mpf(VV(80.3692, 2.085, 2)); REF['ZZ'] = mp.mpf(VV(91.1880, 2.4955, 1))
worst = 0
for k, v in REF.items():
    r = mp.re(v); d = abs((mp.mpf(J[k]) - r)/r); worst = max(worst, d) if k not in ('WW', 'ZZ') else worst
    print(f'{k:10s} reference {mp.nstr(r, 18):>26s}  calculator {J[k]:>24s}  relative difference {mp.nstr(d, 3)}')
print('motion: largest relative difference', mp.nstr(worst, 3))
print("Mercury per century: both legs %.10f, the time leg %.10f arcsec" % tuple(float(mp.re(REF[k]))*36525/(2*math.pi*math.sqrt((0.387098*1.495978707e11)**3/1.32712440018e20)/86400)*206264.80624709636 for k in ('merc_both', 'merc_time')))
