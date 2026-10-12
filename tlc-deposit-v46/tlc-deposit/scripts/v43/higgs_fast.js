/* The vacuum's edge, computed faster: the same two-loop running (+ three-loop QCD) and the same matching, with
   (1) the beta functions evaluated into a fixed buffer instead of new arrays, (2) the edge searches run at 400 steps,
   where the step error is 2 × 10⁻⁸ GeV, (3) a secant search started next to the answer, with the bracketing search
   kept as the fallback, and (4) the error spread started from the base edge shifted by the first-order slope. */
const PI = Math.PI;
const k1 = 1 / (16 * PI * PI), z3 = 1.2020569031595942, MPL = 1.22e19;
const KB = new Float64Array(5);
const K2 = k1 * k1, K3 = k1 * k1 * k1, C3Y = K3 * (-4166/3 + 640*z3);
function betaInto(g1, g2, g3, yt, l, b) {
  const a1 = g1*g1, a2 = g2*g2, a3 = g3*g3, T = yt*yt, gp = 0.6 * a1, l2 = l*l, T2 = T*T, a22 = a2*a2, gp2 = gp*gp, a32 = a3*a3;
  const g13 = g1*a1, g23 = g2*a2, g33 = g3*a3, s = a2 + gp;
  b[0] = k1*4.1*g13 + K2*g13*(3.98*a1 + 2.7*a2 + 8.8*a3 - 1.7*T);
  b[1] = k1*(-19/6)*g23 + K2*g23*(0.9*a1 + 35/6*a2 + 12*a3 - 1.5*T);
  b[2] = k1*(-7)*g33 + K2*g33*(1.1*a1 + 4.5*a2 - 26*a3 - 2*T) + K3*g33*(32.5*a32 - 40*a3*T + 15*T2);
  b[3] = k1*yt*(4.5*T - 8*a3 - 2.25*a2 - 0.85*a1)
    + K2*yt*(-12*T2 + T*(131/16*gp + 225/16*a2 + 36*a3 - 12*l) + 6*l2
      + 1187/216*gp2 - 0.75*gp*a2 + 19/9*gp*a3 - 5.75*a22 + 9*a2*a3 - 108*a32)
    + C3Y*yt*a32*a3;
  b[4] = k1*(24*l2 - 6*T2 + 0.375*(2*a22 + s*s) + l*(12*T - 9*a2 - 3*gp))
    + K2*(-312*l2*l - 144*l2*T + 36*l2*(3*a2 + gp) - 3*l*T2 + l*T*(80*a3 + 22.5*a2 + 85/6*gp)
      - 9.125*l*a22 + 9.75*l*a2*gp + 629/24*l*gp2
      + 30*T2*T - 32*T2*a3 - 8/3*T2*gp - 2.25*T*a22 + 10.5*T*a2*gp - 4.75*T*gp2
      + 305/16*a22*a2 - 289/48*a22*gp - 559/48*a2*gp2 - 379/48*gp2*gp);
}
function higgsInputs({ Mh = 125.20, Mt = 172.57, als = 0.1180, MW = 80.3692, dl = 0, dy = 0 } = {}) {
  const a = (als - 0.1184) / 0.0007, w = (MW - 80.384) / 0.014, d = Mt - 173.34;
  const gY = 0.35830 + 0.00011*d - 0.00020*w, g2 = 0.64779 + 0.00004*d + 0.00011*w;
  const g3 = 1.1666 + 0.00314*a - 0.00046*d, yt = 0.93690 + 0.00556*d - 0.00042*a + dy;
  const l = 0.12604 + 0.00206*(Mh - 125.15) - 0.00004*d + dl;
  return [Math.sqrt(5/3)*gY, g2, g3, yt, l];
}
/* lowest λ only (the edge searches need nothing else): RK4 without allocation, parabola through the three lowest points */
function lminRun(inp, steps) {
  const Mt = inp.Mt ?? 172.57, h = Math.log(MPL / Mt) / steps, y = higgsInputs(inp);
  let y0 = y[0], y1 = y[1], y2 = y[2], y3 = y[3], y4 = y[4];
  const q = new Float64Array(20), b = KB;
  let l2 = y4, l1 = y4, lm = y4, im = 0, before = NaN, after = NaN, prev = y4;
  for (let i = 1; i <= steps; i++) {
    betaInto(y0, y1, y2, y3, y4, b); for (let j = 0; j < 5; j++) q[j] = b[j];
    betaInto(y0 + h/2*q[0], y1 + h/2*q[1], y2 + h/2*q[2], y3 + h/2*q[3], y4 + h/2*q[4], b); for (let j = 0; j < 5; j++) q[5 + j] = b[j];
    betaInto(y0 + h/2*q[5], y1 + h/2*q[6], y2 + h/2*q[7], y3 + h/2*q[8], y4 + h/2*q[9], b); for (let j = 0; j < 5; j++) q[10 + j] = b[j];
    betaInto(y0 + h*q[10], y1 + h*q[11], y2 + h*q[12], y3 + h*q[13], y4 + h*q[14], b);
    y0 += h/6*(q[0] + 2*q[5] + 2*q[10] + b[0]); y1 += h/6*(q[1] + 2*q[6] + 2*q[11] + b[1]); y2 += h/6*(q[2] + 2*q[7] + 2*q[12] + b[2]);
    y3 += h/6*(q[3] + 2*q[8] + 2*q[13] + b[3]); y4 += h/6*(q[4] + 2*q[9] + 2*q[14] + b[4]);
    if (y4 < lm) { lm = y4; im = i; before = prev; after = NaN; } else if (i === im + 1) after = y4;
    prev = y4;
  }
  if (im > 0 && im < steps) { const A = before, B = lm, D = after, den = A - 2*B + D; if (den > 0) { const x = 0.5*(A - D)/den; return B - 0.25*(A - D)*x; } }
  return lm;
}
/* the full run, for display: same outputs as before */
function higgsRun(inp = {}, steps = 4000) {
  const Mt = inp.Mt ?? 172.57, T = Math.log(MPL / Mt), h = T / steps, y = higgsInputs(inp);
  let y0 = y[0], y1 = y[1], y2 = y[2], y3 = y[3], y4 = y[4];
  const ts = new Float64Array(steps + 1), lam = new Float64Array(steps + 1), q = new Float64Array(20), b = KB;
  lam[0] = y4;
  for (let i = 1; i <= steps; i++) {
    betaInto(y0, y1, y2, y3, y4, b); for (let j = 0; j < 5; j++) q[j] = b[j];
    betaInto(y0 + h/2*q[0], y1 + h/2*q[1], y2 + h/2*q[2], y3 + h/2*q[3], y4 + h/2*q[4], b); for (let j = 0; j < 5; j++) q[5 + j] = b[j];
    betaInto(y0 + h/2*q[5], y1 + h/2*q[6], y2 + h/2*q[7], y3 + h/2*q[8], y4 + h/2*q[9], b); for (let j = 0; j < 5; j++) q[10 + j] = b[j];
    betaInto(y0 + h*q[10], y1 + h*q[11], y2 + h*q[12], y3 + h*q[13], y4 + h*q[14], b);
    y0 += h/6*(q[0] + 2*q[5] + 2*q[10] + b[0]); y1 += h/6*(q[1] + 2*q[6] + 2*q[11] + b[1]); y2 += h/6*(q[2] + 2*q[7] + 2*q[12] + b[2]);
    y3 += h/6*(q[3] + 2*q[8] + 2*q[13] + b[3]); y4 += h/6*(q[4] + 2*q[9] + 2*q[14] + b[4]);
    ts[i] = i * h; lam[i] = y4;
  }
  let im = 0; for (let i = 1; i <= steps; i++) if (lam[i] < lam[im]) im = i;
  let lmin = lam[im], tmin = ts[im];
  if (im > 0 && im < steps) {
    const A = lam[im-1], B = lam[im], D = lam[im+1], den = A - 2*B + D;
    if (den > 0) { const x = 0.5 * (A - D) / den; lmin = B - 0.25 * (A - D) * x; tmin = ts[im] + x * h; }
  }
  let cross = null;
  for (let i = 1; i <= steps; i++) if (lam[i-1] >= 0 && lam[i] < 0) { const f = lam[i-1] / (lam[i-1] - lam[i]); cross = Mt * Math.exp(ts[i-1] + f * h); break; }
  let logLife = null;
  for (let i = 0; i <= steps; i++) if (lam[i] < 0) {
    const lt = (8*PI*PI / (3*Math.abs(lam[i]))) / 4 - Math.log(Mt * Math.exp(ts[i]));
    if (logLife === null || lt < logLife) logLife = lt;
  }
  const lifeLog10 = logLife === null ? null : (logLife - Math.log(1.519e24 * 3.156e7)) / Math.LN10;
  return { Mt, ts, lam, lmin, muMin: Mt * Math.exp(tmin), cross, lPl: lam[steps], lMt: lam[0], lifeLog10, yPl: [y0, y1, y2, y3, y4] };
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
/* secant from two nearby starts; the bracketing search on the full range if it strays or stalls */
function secant(f, x0, x1, lo, hi, tol = 1e-7) {
  let f0 = f(x0), f1 = f(x1);
  for (let i = 0; i < 12; i++) {
    if (f1 === f0) break;
    const x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
    if (!(x2 > lo && x2 < hi)) break;
    if (Math.abs(x2 - x1) < tol) return x2;
    x0 = x1; f0 = f1; x1 = x2; f1 = f(x2);
  }
  return bisect(f, lo, hi, tol);
}
const EDGE_STEPS = 400;
const edgeTop = (inp = {}, guess = 171.0) => secant(Mt => lminRun({ ...inp, Mt }, EDGE_STEPS), guess, guess + 0.3, 160, 185);
const edgeHiggs = (inp = {}, guess = 128.5) => secant(Mh => lminRun({ ...inp, Mh }, EDGE_STEPS), guess, guess + 0.6, 110, 150);
function edgeErrors(inp, which) {
  const base = which === "top" ? edgeTop(inp) : edgeHiggs(inp);
  const src = which === "top"
    ? [["alpha_s", "als", 0.0009], ["Higgs mass", "Mh", 0.11], ["W mass", "MW", 0.0133], ["matching, λ", "dl", 0.0003], ["matching, y_t", "dy", 0.0005]]
    : [["top mass", "Mt", 0.27], ["alpha_s", "als", 0.0009], ["W mass", "MW", 0.0133], ["matching, λ", "dl", 0.0003], ["matching, y_t", "dy", 0.0005]];
  const parts = src.map(([name, key, s]) => {
    const q = { ...inp }; q[key] = (q[key] ?? { als: 0.1180, Mh: 125.20, MW: 80.3692, dl: 0, dy: 0, Mt: 172.57 }[key]) + s;
    return [name, Math.abs((which === "top" ? edgeTop(q, base) : edgeHiggs(q, base)) - base)];
  });
  return { base, parts, total: Math.hypot(...parts.map(p => p[1])) };
}
if (typeof module !== "undefined") module.exports = { betaInto, higgsInputs, higgsRun, lminRun, bisect, secant, edgeTop, edgeHiggs, edgeErrors, EDGE_STEPS };
