import math
hbar=1.054571817e-34; kB=1.380649e-23; h=2*math.pi*hbar; e=1.602176634e-19; c=299792458
pi=math.pi
# Ji 2026 IHO: omega_I = 2pi*0.9 kHz -> equivalent temperature hbar*omega/(2pi kB); pure-IHO squeeze at 0.25 ms
wI=2*pi*900.0
print('Ji: T_equiv = hbar wI/(2 pi kB) =', hbar*wI/(2*pi*kB)*1e9, 'nK;  pure IHO squeeze at 0.25 ms:', 20*math.log10(math.exp(wI*0.25e-3)), 'dB; at 0.15 ms', 20*math.log10(math.exp(wI*0.15e-3)))
# Steinhauer: T_H = 0.35 nK -> lambda = 2 pi kB T / hbar
print('Munoz de Nova: lambda = 2 pi kB T/hbar =', 2*pi*kB*0.35e-9/hbar, 's^-1')
# Schmitt: G0 = 1.20 * 4 e^2/h
G0 = 1.20*4*e**2/h; print('Schmitt: 4e^2/h =', 4*e**2/h*1e3, 'mS; 1.20x =', G0*1e3, 'mS; pull', (0.18e-3-G0)/0.02e-3)
# Pardo pulls
print('Pardo pulls: SH0ES', (4.02-4)/0.10, ' Planck', (3.98-4)/0.07)
# Chester vs Lipa
nuE, sE = 0.6709, 0.0001; Ds, sDs = 1.51136, 0.00022
nuB = 1/(3-Ds); snuB = sDs/(3-Ds)**2
print('bootstrap nu =', nuB, '+-', snuB, ' exp nu =', nuE, ' diff/sigma_comb =', (nuB-nuE)/math.hypot(sE, snuB), ' diff/sigma_exp =', (nuB-nuE)/sE)
print('alpha from bootstrap nu: 2-3nu =', 2-3*nuB, ' ; alpha exp -0.0127(3); d_eff = (2-alpha_exp)/nu_B =', (2+0.0127)/nuB)
# LIGO squeezing
for dB in [5.2, 6.1, 10.6]:
    print(f'{dB} dB: variance x{10**(-dB/10):.4f}, amplitude x{10**(-dB/20):.4f}')
# shares vs dark-energy data
for lab, x in [('1-pi^2/32', 1-pi**2/32), ('overlap 3D rho=1', 0.6863271610311136)]:
    rows=[('Planck', .6847, .0073, .0073), ('Planck+DESI DR2', .6973, .0036, .0036), ('Planck+ACT+DESI DR2', .700, .004, .004), ('DESI DR2 BAO+BBN', .7025, .0086, .0086),
          ('DES-SN5YR', .648, .017, .017), ('Pantheon+', .666, .018, .018), ('Union3', .644, .026, .028)]
    out=[]
    for n, m, sp, sm in rows:
        s = sp if m < x else sm   # measured below model -> upper error bar
        out.append(f'{n} {(m-x)/s:+.2f}')
    print(lab, round(x,5), '; '.join(out))
# Hamraoui geometric phase for box 3:2:1 = lune
print('box 3:2:1: (y^4-z^4)/(x^4-z^4) =', (16-1)/(81-1), ' 4 arcsin sqrt =', 4*math.asin(math.sqrt(3/16)), ' dihedral deg', 2*math.degrees(math.asin(math.sqrt(3/16))))
# Higgs bounce with lambda_min from calculator
lm = -0.01034804399435521; S = 8*pi**2/(3*abs(lm)); print('Higgs: lambda_min', lm, ' bounce S =', S, ' load', abs(lm)/(8*pi**2), ' load*S', abs(lm)/(8*pi**2)*S)
# dimension ladder ratios
for n in [2,3,4]: print(n, '2pi/n =', 2*pi/n, '=', math.degrees(2*pi/n), 'deg')
# Unruh / Hawking quick
print('Unruh T at 1 g:', hbar*9.81/(2*pi*c*kB), 'K')
