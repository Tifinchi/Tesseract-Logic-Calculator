"""Rebuild the 998-number alphabet of section 8.1 and count members near each target."""
import math, itertools
r2, r3, r5 = math.sqrt(2), math.sqrt(3), math.sqrt(5); phi = (1+r5)/2; pi = math.pi
base = [1, r2, r3, 2, 1/r2, 1/r3, 0.5, 1/3, 0.25, 2/3, 0.75, 1.5, r3/r2, r3/2,
        pi/4, pi/6, pi**2/32, 1-pi/4, 1-pi/6, 1-pi**2/32, 4/pi, 2*pi/3, 1/(2*pi), 1/(8*pi**2),
        phi, 1/phi, r5, 3/8, 5/3, 3, 4, 6, 8, 12, 16, 24, 27]
assert len(base) == 37
vals = set()
cand = list(base)
for a, b in itertools.product(base, repeat=2):
    cand += [a*b, a/b]
cand.sort()
merged = []
for v in cand:
    if not merged or abs(v-merged[-1]) > 1e-9*max(1, abs(v)):
        merged.append(v)
print("alphabet size:", len(merged), "range", min(merged), max(merged))
targets = [("Omega_Lambda Planck", 0.6847, 0.0073, 0.0073),
           ("Omega_c/Omega_b", 5.364, 0.065, 0.065),
           ("a0/cH0", 0.1834, 0.0031, 0.0031),
           ("n_s Planck", 0.9649, 0.0042, 0.0042),
           ("Cabibbo lambda PDG 2026", 0.22517, 0.00068, 0.00068),
           ("sin2 theta12 NuFIT 6.1", 0.3096, 0.0073, 0.0057),
           ("Delta Sigma HERMES", 0.330, 0.039, 0.039),
           ("Delta G 10 GeV2", 0.231, 0.047, 0.047),
           ("H0 t0 Planck", 0.9510, 0.004, 0.004),
           ("Koide Q", 0.6666645, 0.0000051, 0.0000051),
           ("lepton phase", 0.222225, 0.000006, 0.000006),
           ("r_p/lambdabar_p", 3.9977, 0.0030, 0.0030),
           ("Tully-Fisher slope", 3.85, 0.09, 0.09),
           ("Omega_m Planck", 0.3153, 0.0073, 0.0073)]
for name, m, dn, up in targets:
    within2 = [v for v in merged if m-2*dn <= v <= m+2*up]
    within3 = [v for v in merged if m-3*dn <= v <= m+3*up]
    near = min(merged, key=lambda v: abs(v-m))
    s = up if near > m else dn
    print(f"{name:26s} within 2σ: {len(within2):4d}  within 3σ: {len(within3):4d}  nearest {near:.6f} pull {(m-near)/s:+.2f}")
print("--- other fields")
more = [("Kolmogorov C", 1.62, 0.17), ("von Karman kappa", 0.41, 0.02), ("Strouhal", 0.20, 0.01), ("sphere drag", 0.47, 0.03),
        ("random close packing", 0.64, 0.01), ("Poisson metals", 0.30, 0.04), ("Tully-Fisher 2012", 3.98, 0.06),
        ("Kleiber 0.737", 0.737, 0.026), ("Kleiber 0.686", 0.686, 0.014), ("m_Delta/m_p", 1.3135, 0.002),
        ("spin-1/2 period / 360", 716.8/360, 3.8/360)]
for name, m, s in more:
    w2 = [v for v in merged if abs(v-m) <= 2*s]; w3 = [v for v in merged if abs(v-m) <= 3*s]
    near = min(merged, key=lambda v: abs(v-m))
    print(f"{name:24s} within 2σ: {len(w2):4d} within 3σ: {len(w3):4d} nearest {near:.6f} pull {(m-near)/s:+.2f}")
# density near 2/3
lo, hi = 2/3*0.75, 2/3*1.25
dens = len([v for v in merged if lo <= v <= hi])/(hi-lo)
print(f"density near 2/3: {dens:.0f} per unit")
