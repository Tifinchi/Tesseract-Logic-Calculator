const PI = Math.PI;
  const k1 = 1 / (16 * PI * PI), z3 = 1.2020569031595942, MPL = 1.22e19;
  function beta(y) {
    const [g1, g2, g3, yt, l] = y;
    const a1 = g1*g1, a2 = g2*g2, a3 = g3*g3, T = yt*yt, gp = 0.6 * a1;
    const b = [k1*41/10*g1**3, k1*(-19/6)*g2**3, k1*(-7)*g3**3,
      k1*yt*(4.5*T - 8*a3 - 2.25*a2 - 0.85*a1),
      k1*(24*l*l - 6*T*T + 0.375*(2*a2*a2 + (a2 + gp)**2) + l*(12*T - 9*a2 - 3*gp))];
    const K = k1 * k1;
    b[0] += K*g1**3*(199/50*a1 + 27/10*a2 + 44/5*a3 - 17/10*T);
    b[1] += K*g2**3*(9/10*a1 + 35/6*a2 + 12*a3 - 1.5*T);
    b[2] += K*g3**3*(11/10*a1 + 4.5*a2 - 26*a3 - 2*T);
    b[3] += K*yt*(-12*T*T + T*(131/16*gp + 225/16*a2 + 36*a3 - 12*l) + 6*l*l
      + 1187/216*gp*gp - 0.75*gp*a2 + 19/9*gp*a3 - 23/4*a2*a2 + 9*a2*a3 - 108*a3*a3);
    b[4] += K*(-312*l**3 - 144*l*l*T + 36*l*l*(3*a2 + gp) - 3*l*T*T + l*T*(80*a3 + 22.5*a2 + 85/6*gp)
      - 73/8*l*a2*a2 + 39/4*l*a2*gp + 629/24*l*gp*gp
      + 30*T**3 - 32*T*T*a3 - 8/3*T*T*gp - 9/4*T*a2*a2 + 10.5*T*a2*gp - 19/4*T*gp*gp
      + 305/16*a2**3 - 289/48*a2*a2*gp - 559/48*a2*gp*gp - 379/48*gp**3);
    const K3 = k1 ** 3;
    b[2] += K3*g3**3*(65/2*a3*a3 - 40*a3*T + 15*T*T);
    b[3] += K3*yt*(-4166/3 + 640*z3)*a3**3;
    return b;
  }
  function higgsInputs({ Mh = 125.20, Mt = 172.57, als = 0.1180, MW = 80.3692, dl = 0, dy = 0 } = {}) {
    // matching at mu = Mt (Buttazzo et al. 2013, fit formulas)
    const a = (als - 0.1184) / 0.0007, w = (MW - 80.384) / 0.014, d = Mt - 173.34;
    const gY = 0.35830 + 0.00011*d - 0.00020*w, g2 = 0.64779 + 0.00004*d + 0.00011*w;
    const g3 = 1.1666 + 0.00314*a - 0.00046*d, yt = 0.93690 + 0.00556*d - 0.00042*a + dy;
    const l = 0.12604 + 0.00206*(Mh - 125.15) - 0.00004*d + dl;
    return [Math.sqrt(5/3)*gY, g2, g3, yt, l];
  }
  function higgsRun(inp = {}, steps = 4000) {
    const Mt = inp.Mt ?? 172.57, T = Math.log(MPL / Mt), h = T / steps;
    let y = higgsInputs(inp);
    const ts = new Float64Array(steps + 1), lam = new Float64Array(steps + 1);
    lam[0] = y[4];
    const add = (u, v, s) => u.map((x, i) => x + s * v[i]);
    for (let i = 1; i <= steps; i++) {
      const q1 = beta(y), q2 = beta(add(y, q1, h/2)), q3 = beta(add(y, q2, h/2)), q4 = beta(add(y, q3, h));
      y = y.map((x, j) => x + h/6 * (q1[j] + 2*q2[j] + 2*q3[j] + q4[j]));
      ts[i] = i * h; lam[i] = y[4];
    }
    let im = 0; for (let i = 1; i <= steps; i++) if (lam[i] < lam[im]) im = i;
    let lmin = lam[im], tmin = ts[im];
    if (im > 0 && im < steps) {                         // parabola through the three lowest points
      const A = lam[im-1], B = lam[im], D = lam[im+1], den = A - 2*B + D;
      if (den > 0) { const x = 0.5 * (A - D) / den; lmin = B - 0.25 * (A - D) * x; tmin = ts[im] + x * h; }
    }
    let cross = null;
    for (let i = 1; i <= steps; i++) if (lam[i-1] >= 0 && lam[i] < 0) { const f = lam[i-1] / (lam[i-1] - lam[i]); cross = Mt * Math.exp(ts[i-1] + f * h); break; }
    // leading-order tunnelling: bounce action S = 8π²/(3|λ|), rate ~ μ^4 e^{-S}
    let logLife = null;
    for (let i = 0; i <= steps; i++) if (lam[i] < 0) {
      const lt = (8*PI*PI / (3*Math.abs(lam[i]))) / 4 - Math.log(Mt * Math.exp(ts[i]));
      if (logLife === null || lt < logLife) logLife = lt;
    }
    const lifeLog10 = logLife === null ? null : (logLife - Math.log(1.519e24 * 3.156e7)) / Math.LN10;
    return { Mt, ts, lam, lmin, muMin: Mt * Math.exp(tmin), cross, lPl: lam[steps], lMt: lam[0], lifeLog10, yPl: y };
  }
  function bisect(f, lo, hi, tol = 2e-5) {     // Illinois root finder on a bracket (falls back to halving)
    let fa = f(lo), fb = f(hi), side = 0;
    if ((fa > 0) === (fb > 0)) return NaN;
    for (let i = 0; i < 60 && hi - lo > tol; i++) {
      let x = (lo * fb - hi * fa) / (fb - fa);
      if (!(x > lo && x < hi)) x = (lo + hi) / 2;
      const fx = f(x);
      if (Math.abs(fx) < 1e-12) return x;
      if ((fx > 0) === (fb > 0)) { hi = x; fb = fx; if (side === -1) fa /= 2; side = -1; }
      else { lo = x; fa = fx; if (side === 1) fb /= 2; side = 1; }
      if (Math.abs(hi - lo) < tol) break;
    }
    return (lo * fb - hi * fa) / (fb - fa);
  }
  // exact-edge masses: lowest lambda just touches 0
  const edgeTop = (inp = {}) => bisect(Mt => higgsRun({ ...inp, Mt }, 2500).lmin, 160, 185);
  const edgeHiggs = (inp = {}) => bisect(Mh => higgsRun({ ...inp, Mh }, 2500).lmin, 110, 150);
  function edgeErrors(inp, which) {          // spread from the input uncertainties (one at a time, in quadrature)
    const base = which === "top" ? edgeTop(inp) : edgeHiggs(inp);
    const src = which === "top"
      ? [["alpha_s", "als", 0.0009], ["Higgs mass", "Mh", 0.11], ["W mass", "MW", 0.0133], ["matching, λ", "dl", 0.0003], ["matching, y_t", "dy", 0.0005]]
      : [["top mass", "Mt", 0.27], ["alpha_s", "als", 0.0009], ["W mass", "MW", 0.0133], ["matching, λ", "dl", 0.0003], ["matching, y_t", "dy", 0.0005]];
    const parts = src.map(([name, key, s]) => {
      const q = { ...inp }; q[key] = (q[key] ?? { als: 0.1180, Mh: 125.20, MW: 80.3692, dl: 0, dy: 0, Mt: 172.57 }[key]) + s;
      return [name, Math.abs((which === "top" ? edgeTop(q) : edgeHiggs(q)) - base)];
    });
    return { base, parts, total: Math.hypot(...parts.map(p => p[1])) };
  }


module.exports = { beta, higgsInputs, higgsRun, bisect, edgeTop, edgeHiggs, edgeErrors };
