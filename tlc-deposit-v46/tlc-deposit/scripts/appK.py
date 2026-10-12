"""Appendix K: Koide numbers with PDG 2026 lepton masses; tilt from NuFIT; TM1/TM2 pulls."""
from mpmath import mp, mpf, sqrt, acos, cos, sin, pi, degrees, atan2, asin
mp.dps = 30
me, mmu, mtau, stau = mpf('0.51099895069'), mpf('105.6583755'), mpf('1776.93'), mpf('0.09')
def koide(mt):
    s = sqrt(me)+sqrt(mmu)+sqrt(mt); Q = (me+mmu+mt)/s**2
    a = sqrt(2*(3*Q-1))   # Q = (1 + a^2/2)/3
    ang = degrees(acos((s/sqrt(3))/sqrt(me+mmu+mt)))
    return Q, a, ang
Q, a, ang = koide(mtau)
Qp, ap, angp = koide(mtau+stau)
print(f"Q = {mp.nstr(Q,10)} +- {mp.nstr(abs(Qp-Q),3)}; pull vs 2/3 {mp.nstr((Q-mpf(2)/3)/abs(Qp-Q),3)}")
print(f"a = {mp.nstr(a,10)} +- {mp.nstr(abs(ap-a),3)}; angle = {mp.nstr(ang,10)} +- {mp.nstr(abs(angp-ang),3)}")
# phase delta: sqrt(m_k) = mu[1 + a cos(delta + 2 pi k/3)], k=0 e? fit with a=sqrt2 convention: use the measured a
s = sqrt(me)+sqrt(mmu)+sqrt(mtau); mu = s/3
# tau is the largest: sqrt(m_tau) = mu(1 + a cos delta)
d = acos((sqrt(mtau)/mu-1)/a)
print(f"delta = {mp.nstr(d,10)} rad (2/9 = {mp.nstr(mpf(2)/9,10)})")
# tilt from sin^2 theta13 = sin^2 t / 3
s13, ss13 = mpf('0.02200'), mpf('0.00056')
t = asin(sqrt(3*s13)); dt = 3*ss13/sin(2*t)
print(f"tilt = {mp.nstr(degrees(t),6)} +- {mp.nstr(degrees(dt),3)} deg; 15 deg pull {mp.nstr((degrees(t)-15)/degrees(dt),3)}")
# TM1 and TM2 sin^2 theta12 at measured theta13
tm1 = (1-3*s13)/(3*(1-s13)); tm2 = 1/(3*(1-s13))
print(f"TM1 s12^2 = {mp.nstr(tm1,6)}; TM2 s12^2 = {mp.nstr(tm2,6)}")
print(f"pulls vs 0.3096 (+0.0057/-0.0073): TM1 {mp.nstr((mpf('0.3096')-tm1)/mpf('0.0057'),3)}; TM2 {mp.nstr((mpf('0.3096')-tm2)/mpf('0.0057'),3)}")
s15 = sin(pi/12)**2/3; tm1_15 = (1-3*s15)/(3*(1-s15))
print(f"at 15 deg: s13^2 = {mp.nstr(s15,6)}, s12^2 = {mp.nstr(tm1_15,6)}; pulls {mp.nstr((s13-s15)/ss13,3)}, {mp.nstr((mpf('0.3096')-tm1_15)/mpf('0.0057'),3)}")
print(f"spin-1/2: (716.8-720)/3.8 = {mp.nstr((mpf('716.8')-720)/mpf('3.8'),3)}")
print(f"tau prediction from e, mu: see F.6 (1776.969)")
