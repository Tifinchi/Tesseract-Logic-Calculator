/* The law of motion: rate = J·∇H. On a mass M's field the sum table sits on two legs, the time factor A and the space factor B of
   the isotropic Schwarzschild metric: H = √A·√(m²c⁴ + p²c²/B). Units GM/c² = 1, c = 1. J sends ∇H to (∂H/∂p, −∂H/∂x). */
const MOTION = (() => {
  const C = 299792458, GMS = 1.32712440018e20, AU = 1.495978707e11, RSUN = 6.957e8, AS = 206264.80624709636;
  const legs = r => { const u = 1 / (2 * r), A = ((1 - u) / (1 + u)) ** 2, B = (1 + u) ** 4, du = -2 * u * u;
    return [A, B, -4 * (1 - u) / (1 + u) ** 3 * du, 4 * (1 + u) ** 3 * du]; };
  const rate = (kind, light) => z => {
    const x = z[0], y = z[1], px = z[2], py = z[3], r = Math.hypot(x, y), p2 = px * px + py * py;
    if (kind === "newton") { const g = 1 / (r * r * r); return [px, py, -g * x, -g * y]; }
    let [A, B, dA, dB] = legs(r); if (kind === "time") { B = 1; dB = 0; }
    if (light) { const p = Math.sqrt(p2), s = Math.sqrt(A / B), ds = 0.5 / s * (dA / B - A * dB / (B * B)), k = -p * ds / r;
      return [s * px / p, s * py / p, k * x, k * y]; }
    const S = Math.sqrt(1 + p2 / B), sA = Math.sqrt(A), hp = sA / (B * S), hr = dA / (2 * sA) * S - sA * p2 * dB / (2 * B * B * S);
    return [hp * px, hp * py, -hr * x / r, -hr * y / r];
  };
  const step = (f, z, h) => { const a = f(z), b = f(z.map((v, i) => v + h / 2 * a[i])), c = f(z.map((v, i) => v + h / 2 * b[i])), d = f(z.map((v, i) => v + h * c[i]));
    return z.map((v, i) => v + h / 6 * (a[i] + 2 * b[i] + 2 * c[i] + d[i])); };
  /* perihelion advance per orbit, in radians, for an orbit of semi-major axis aAU and eccentricity e about one solar mass */
  function advance(aAU, e, kind, N = 20000) {
    const a = aAU * AU / (GMS / (C * C)), rp = a * (1 - e), ra = a * (1 + e);
    const L = Math.sqrt((2 / ra - 2 / rp) / ((1 - 2 / ra) / (ra * ra) - (1 - 2 / rp) / (rp * rp))), r0 = (rp - 1 + Math.sqrt(rp * rp - 2 * rp)) / 2;
    const f = rate(kind), T = 2 * Math.PI * a ** 1.5, h = T / N, rad = q => q[0] * q[2] + q[1] * q[3];
    let z = [r0, 0, 0, L / r0];
    for (let i = 0; i < 1.5 * N; i++) {
      const z1 = step(f, z, h);
      if (i > N / 2 && rad(z) < 0 && rad(z1) >= 0) {
        let lo = 0, hi = h;
        for (let k = 0; k < 60; k++) { const m = (lo + hi) / 2; if (rad(step(f, z, m)) < 0) lo = m; else hi = m; }
        const q = step(f, z, (lo + hi) / 2); return Math.atan2(q[1], q[0]);
      }
      z = z1;
    }
    return NaN;
  }
  /* light passing a mass Msun (solar masses) at impact parameter bR (solar radii): total bending in arcseconds */
  function bending(bR, Msun, kind) {
    const b = bR * RSUN / (GMS * Msun / (C * C)), X = 2e4 * b, f = rate(kind, true);
    let z = [-X, b, 1, 0];
    while (z[0] < X) { const r = Math.hypot(z[0], z[1]); z = step(f, z, 0.002 * Math.max(r, b)); }
    return -Math.atan2(z[3], z[2]) * AS;
  }
  /* a clock on a circular orbit of radius rkm about the Earth, against clocks on the geoid (rate L_G below coordinate time) */
  function clock(rkm) {
    const GM = 3.986004418e14, LG = 6.969290134e-10, r = rkm * 1e3, u = GM / (r * C * C);
    return { frac: LG - (1 - Math.sqrt(1 - 3 * u)), grav: LG - u, speed: -u / 2 };
  }
  /* hydrogen Lyman-alpha: Schrödinger (no spin), Dirac (spin), and Dirac with the board's 1S sub-reaction (the Lamb shift) */
  function lymanAlpha(lamb1S = 8.1727441e9) {
    const al = 7.2973525643e-3, Rinf = 10973731.568157, red = 1 / (1 + 1 / 1836.152673426), hz = 2 * Rinf * C / (al * al) * red;
    const dirac = (n, j) => { const k = j + 0.5, d = k - Math.sqrt(k * k - al * al); return Math.pow(1 + (al / (n - d)) ** 2, -0.5) - 1; };
    const f12 = (dirac(2, 0.5) - dirac(1, 0.5)) * hz, f32 = (dirac(2, 1.5) - dirac(1, 0.5)) * hz;
    return { schr: 1e9 / (Rinf * red * 0.75), d12: C / f12 * 1e9, d32: C / f32 * 1e9, l12: C / (f12 - lamb1S) * 1e9, l32: C / (f32 - lamb1S) * 1e9 };
  }
  const perCentury = (rad, aAU) => rad * 36525 / (2 * Math.PI * Math.sqrt((aAU * AU) ** 3 / GMS) / 86400) * AS;
  return { advance, bending, clock, lymanAlpha, perCentury, AS };
})();
if (typeof module !== "undefined") module.exports = MOTION;
