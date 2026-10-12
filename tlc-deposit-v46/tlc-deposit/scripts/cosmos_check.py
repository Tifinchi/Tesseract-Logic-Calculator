from mpmath import mp, mpf, pi, sqrt, asinh, log10, cbrt
mp.dps = 30
# constants (CODATA 2022)
c = mpf(299792458); G = mpf('6.67430e-11'); hbarc_GeVfm = mpf('0.1973269804'); eV = mpf('1.602176634e-19')
Mpc = mpf('3.0856775814913673e22'); Gpc = Mpc*1000; Msun = mpf('1.98841e30')
H0 = mpf('67.4')*1000/Mpc
# hum
v = mpf('246.22'); M23 = mpf('3.00577e16')   # the alpha_2 = alpha_3 crossing of higgs2.py's running at PDG 2026 inputs (3.0058e16 GeV)
hum = v**2/M23
print('hum meV', hum*1e12)
OL = mpf('0.685')
rho_c = 3*H0**2/(8*pi*G)
rhoL_J = OL*rho_c*c**2
print('rho_L J/m3', rhoL_J)
GeV_per_m3 = rhoL_J/(eV*1e9)
rhoL_GeV4 = GeV_per_m3*(hbarc_GeVfm*mpf('1e-15'))**3
print('rho_L GeV4', rhoL_GeV4, 'fourth root meV', rhoL_GeV4**mpf('0.25')*1e12)
mH = mpf('125.20')
V0 = mH**2*v**2/8
print('Higgs depth GeV4', V0, 'log10 ratio', log10(V0/rhoL_GeV4))
print('pull hum', (rhoL_GeV4**mpf('0.25')*1e12 - hum*1e12)/mpf('0.20'))
# share
Om = pi**2/32; OLc = 1-Om
print('Om', Om, 'OL', OLc)
for name, m, s in [('SN+CMB+BAO',0.3045,0.0032),('Planck',0.3153,0.0073),('DESI BAO',0.2975,0.0086)]:
    print(name, (mpf(m)-Om)/mpf(s))
zeq = cbrt(OLc/Om)-1; zacc = cbrt(2*OLc/Om)-1
print('z equal', zeq, 'z accel', zacc)
for om in [mpf('0.3153')]:
    print('Planck z equal', cbrt((1-om)/om)-1, 'z acc', cbrt(2*(1-om)/om)-1)
print('q0', 3*Om/2-1, 'Planck q0', 3*mpf('0.3153')/2-1)
print('ratio at z=1', Om*8/OLc)
# E(z) at 0.349
z = mpf('0.349')
E_calc = sqrt(Om*(1+z)**3+OLc); E_pl = sqrt(mpf('0.3153')*(1+z)**3+mpf('0.6847'))
print('E calc', E_calc, 'E Planck', E_pl, 'pulls', (mpf('1.247')-E_calc)/mpf('0.006'), (mpf('1.247')-E_pl)/mpf('0.006'), 'excess %', (mpf('1.247')/E_pl-1)*100)
# H0 t0
def H0t0(om):
    ol = 1-om
    return 2/(3*sqrt(ol))*asinh(sqrt(ol/om))
print('H0t0 calc', H0t0(Om), 'Planck', H0t0(mpf('0.3153')), '2(1-pi/6)', 2*(1-pi/6))
# Hubble sphere
RH = c/H0; MH = c**3/(2*G*H0)
print('c/H0 m', RH, 'Gpc', RH/Gpc, 'mass kg', MH, 'Rs/RH', 2*G*MH/c**2/RH)
print('curv radius Gpc', RH/Gpc/sqrt(mpf('0.0045')), 'x Hubble', 1/sqrt(mpf('0.0045')), 'parent kg', MH/sqrt(mpf('0.0045')))
# OS collapse times
tM = pi*G*Msun/c**3
print('pi GM/c3 per Msun s', tM, '2.6Msun us', tM*mpf('2.6')*1e6, 'SgrA* 4.3e6 s', tM*mpf('4.3e6'), 'M87* 6.5e9 h', tM*mpf('6.5e9')/3600)
# interior volume
Mg = G*Msun/c**2; t = mpf('1e9')*mpf('365.25')*86400*c
V = 3*sqrt(3)*pi*Mg**2*t; Vsun = mpf(4)/3*pi*mpf('6.957e8')**3
print('volume m3', V, 'solar volumes', V/Vsun, '2GM/c2 km', 2*Mg/1000)
# inflation ns
for N in (54,55,57,60):
    print('N', N, '1-2/N', 1-mpf(2)/N)
for name,m,s in [('Planck18',0.9649,0.0042),('+DESI',0.9690,0.0035),('SPT+P+ACT+DESI',0.9737,0.0025),('PDG2026',0.975,0.003)]:
    print(name, 'calc 26/27 pull', (mpf(m)-mpf(26)/27)/mpf(s), 'N=55', (mpf(m)-(1-mpf(2)/55))/mpf(s), 'N=60', (mpf(m)-(1-mpf(2)/60))/mpf(s))
# w
print('w pull', (mpf('-0.995')-(-1))/mpf('0.019'))
# q0 pulls
print('q0 pulls', (mpf('-0.537')-(3*Om/2-1))/mpf('0.026'), (mpf('-0.51')-(3*Om/2-1))/mpf('0.024'))
print('q0 Planck pulls', (mpf('-0.537')-(3*mpf('0.3153')/2-1))/mpf('0.026'), (mpf('-0.51')-(3*mpf('0.3153')/2-1))/mpf('0.024'))
# curvature pulls
for name,m,s in [('Planck+BAO',0.0007,0.0019),('DESI+CMB',0.0023,0.0011),('CMB alone',-0.011,0.006)]:
    print(name, mpf(m)/mpf(s))
# a0/cH0 MOND
print('MOND 1/2pi', 1/(2*pi), 'pull stat', (mpf('0.1834')-1/(2*pi))/mpf('0.0031'), 'with sys', (mpf('0.1834')-1/(2*pi))/sqrt(mpf('0.0031')**2+mpf('0.037')**2))
print('3/16 pull', (mpf('0.1834')-mpf(3)/16)/mpf('0.0031'))
# vacuum energy Planck cutoff
MP = mpf('1.220890e19')
print('Planck cutoff ratio log10', log10(MP**4/rhoL_GeV4), 'reduced', log10((MP/sqrt(8*pi))**4/rhoL_GeV4))
