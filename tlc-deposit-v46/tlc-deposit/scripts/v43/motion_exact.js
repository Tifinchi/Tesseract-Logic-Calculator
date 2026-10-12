/* The law of motion, integrated exactly. The orbit and the light path are the same rate = J·∇H as before; instead of
   stepping along them, the angle swept between turning points is written as one integral over a smooth periodic
   variable, with the turning points factored out by divided differences so that no small quantity is formed as the
   difference of two large ones. Units GM/c² = 1, c = 1. */
const MOTIONX = (() => {
  const C = 299792458, GMS = 1.32712440018e20, AU = 1.495978707e11, RSUN = 6.957e8, AS = 206264.80624709636;
  /* trapezoid rule on [0, π] for an even periodic integrand: spectrally accurate; doubles until two estimates agree */
  function periodic(f) {
    let n = 16, prev = NaN;
    for (;;) {
      let s = 0.5 * (f(0) + f(Math.PI));
      for (let k = 1; k < n; k++) s += f(Math.PI * k / n);
      const v = s * Math.PI / n;
      if (Math.abs(v - prev) <= 1e-15 * Math.abs(v) || n >= 1 << 16) return v;
      prev = v; n *= 2;
    }
  }
  /* Gauss–Legendre nodes and weights on [−1, 1] (Newton on the Legendre recurrence), cached by order */
  const GL = {};
  function gauss(n) {
    if (GL[n]) return GL[n];
    const x = new Float64Array(n), w = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      let z = Math.cos(Math.PI * (i + 0.75) / (n + 0.5)), dp = 0;
      for (let it = 0; it < 100; it++) {
        let p0 = 1, p1 = z;
        for (let k = 2; k <= n; k++) { const p2 = ((2 * k - 1) * z * p1 - (k - 1) * p0) / k; p0 = p1; p1 = p2; }
        dp = n * (z * p1 - p0) / (z * z - 1);
        const dz = p1 / dp; z -= dz;
        if (Math.abs(dz) < 1e-16) break;
      }
      x[i] = z; w[i] = 2 / ((1 - z * z) * dp * dp);
    }
    return (GL[n] = { x, w });
  }
  function interval(f, a, b) {                 // Gauss–Legendre, doubling the order until two estimates agree
    let n = 16, prev = NaN;
    for (;;) {
      const { x, w } = gauss(n), m = (a + b) / 2, h = (b - a) / 2;
      let s = 0; for (let i = 0; i < n; i++) s += w[i] * f(m + h * x[i]);
      const v = s * h;
      if (Math.abs(v - prev) <= 1e-15 * Math.abs(v) || n >= 512) return v;
      prev = v; n *= 2;
    }
  }
  const dev = y => y / (Math.sqrt(1 + y) * (1 + Math.sqrt(1 + y)));   // 1 − (1 + y)^(−1/2), formed without cancellation

  /* perihelion advance per orbit (radians) for semi-major axis aAU and eccentricity e about one solar mass; the same
     orbit as the stepped engine: Schwarzschild turning points a(1 ∓ e), started from rest in r at the perihelion */
  function advance(aAU, e, kind) {
    const a = aAU * AU / (GMS / (C * C)), rp = a * (1 - e), ra = a * (1 + e);
    if (kind === "newton") return 0;                                  // Kepler's ellipse closes exactly
    if (kind === "both") {                                            // (du/dφ)² = 2(u − uₐ)(uₚ − u)(u₃ − u), u₃ = ½ − uₐ − uₚ
      const ua = 1 / ra, up = 1 / rp, m = (ua + up) / 2, hw = (up - ua) / 2;
      return 2 * periodic(t => { const x = 2 * (ua + up + m - hw * Math.cos(t)); return -dev(-x); });
    }
    // the time leg alone: A isotropic, B = 1. Γ = 1/A = 1 − 8φ + 16φ², φ = 1/(2 − w), w = 1/r.
    const L2 = (2 / ra - 2 / rp) / ((1 - 2 / ra) / (ra * ra) - (1 - 2 / rp) / (rp * rp));
    const r0 = (rp - 1 + Math.sqrt(rp * rp - 2 * rp)) / 2, w0 = 1 / r0, phi = w => 1 / (2 - w);
    const G0 = (1 + 2 * w0 / (2 - w0)) ** 2;                          // Γ(w₀) = ((2 + w₀)/(2 − w₀))²
    const dd1 = (x, y) => { const p = phi(x), q = phi(y); return p * q * (16 * (p + q) - 8); };
    const dd2 = (x, y, z) => { const p = phi(x), q = phi(y), r = phi(z); return p * q * r * (16 * (p + q + r) - 8); };
    const k0 = 1 + L2 * w0 * w0, Q = w => -dd1(w0, w) * k0 + L2 * G0 * (w0 + w), dQ = w => -k0 * dd2(w0, w, w) + L2 * G0;
    let w1 = 2 / L2 - w0;
    for (let i = 0; i < 8; i++) { const d = Q(w1) / dQ(w1); w1 -= d; if (Math.abs(d) < 1e-17 * w1) break; }
    const m = (w0 + w1) / 2, hw = (w0 - w1) / 2;
    return 2 * periodic(t => { const w = m - hw * Math.cos(t), eta = -k0 * dd2(w0, w1, w) / (L2 * G0); return -dev(eta); });
  }

  /* light from infinity past a mass Msun (solar masses) at impact parameter bR (solar radii): total bending, arcseconds.
     The path bends by the index n = √(B/A): (2 + w)/(2 − w) with the time leg alone, (2 + w)³/(4(2 − w)) with both. */
  function bending(bR, Msun, kind) {
    const b = bR * RSUN / (GMS * Msun / (C * C));
    if (kind === "newton") return 2 * Math.atan(1 / b) * AS;          // a particle at c from infinity on Kepler's hyperbola
    const n = kind === "time" ? w => (2 + w) / (2 - w) : w => (2 + w) ** 3 / (4 * (2 - w));
    // closest approach: the first root of b·w = n(w), below the w where n(w)/w is least (the light ring); none: captured
    const wc = kind === "time" ? 2 * (Math.SQRT2 - 1) : 1 / (1 + Math.sqrt(3) / 2);
    if (b * wc <= n(wc)) return NaN;
    let lo = 0, hi = wc;
    for (let i = 0; i < 200 && hi - lo > 1e-17 * hi; i++) { const m = (lo + hi) / 2; if (b * m < n(m)) lo = m; else hi = m; }
    const w0 = (lo + hi) / 2;
    const a0 = 2 + w0, n0 = n(w0);
    const eps = kind === "time"
      ? s => { const rho = n(w0 * s) / n0; return 4 * w0 * (1 + rho) / ((2 - w0 * s) * a0 * (1 + s)); }
      : s => { const c = 2 + w0 * s, rho = n(w0 * s) / n0; return w0 * (1 + rho) * (4 * (a0 * a0 + a0 * c + c * c) - a0 * c * (a0 + c)) / ((2 - w0 * s) * a0 ** 3 * (1 + s)); };
    return 2 * interval(t => -dev(-eps(Math.sin(t))), 0, Math.PI / 2) * AS;
  }
  const perCentury = (rad, aAU) => rad * 36525 / (2 * Math.PI * Math.sqrt((aAU * AU) ** 3 / GMS) / 86400) * AS;
  return { advance, bending, perCentury };
})();
if (typeof module !== "undefined") module.exports = MOTIONX;
