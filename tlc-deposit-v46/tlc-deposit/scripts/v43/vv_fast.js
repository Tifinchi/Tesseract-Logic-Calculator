/* H -> VV* with both bosons off shell, integrated exactly. By symmetry the lighter boson (q₂ ≤ mH/2) stays far from
   its peak; the heavier one's Breit–Wigner is flattened by s = m² + mΓ sinh u (ds·BW = du/(π cosh u)); the outer
   integral is split where the inner range's end crosses the peak, with a sinh map of width Γ/2 there; double-exponential
   (tanh–sinh) quadrature in both variables. */
const PI = Math.PI, SQ2 = Math.SQRT2;
const TS = {};
function tanhSinh(level) {
  if (TS[level]) return TS[level];
  const h = 2 ** -level, xs = [], ws = [];
  for (let k = -Math.ceil(4.5 / h); k <= Math.ceil(4.5 / h); k++) {
    const t = k * h, u = PI / 2 * Math.sinh(t), c = Math.cosh(u), x = Math.tanh(u), w = h * PI / 2 * Math.cosh(t) / (c * c);
    if (w < 1e-30 || Math.abs(x) >= 1) continue;
    xs.push(x); ws.push(w);
  }
  return (TS[level] = { xs: Float64Array.from(xs), ws: Float64Array.from(ws) });
}
function ts(f, a, b, level) { const { xs, ws } = tanhSinh(level), m = (a + b) / 2, h = (b - a) / 2; let s = 0; for (let i = 0; i < xs.length; i++) s += ws[i] * f(m + h * xs[i]); return s * h; }
function VVfast(mH, mV, GV, delta, GF, level = 4) {
  const m2 = mV * mV, mg = mV * GV, H2 = mH * mH;
  const inner = q2 => {                                   // ∫ BW(s₁) F ds₁ over s₁ from q₂² to (mH − q₂)²
    const y = q2 * q2 / H2, lo = Math.asinh((q2 * q2 - m2) / mg), hi = Math.asinh(((mH - q2) ** 2 - m2) / mg);
    if (!(hi > lo)) return 0;
    return ts(u => { const x = (m2 + mg * Math.sinh(u)) / H2, l = Math.max(0, (1 - x - y) ** 2 - 4 * x * y); return Math.sqrt(l) * (l + 12 * x * y) / Math.cosh(u); }, lo, hi, level) / PI;
  };
  const bw = s => (mg / PI) / ((s - m2) ** 2 + mg * mg);
  const outer = q2 => bw(q2 * q2) * 2 * q2 * inner(q2);
  const qs = mH - mV, w = GV / 2, top = mH / 2;
  let tot;
  if (qs > 0 && qs < top) {
    const g = v => { const q2 = qs + w * Math.sinh(v); return outer(q2) * w * Math.cosh(v); };
    tot = ts(g, Math.asinh(-qs / w), 0, level) + ts(g, 0, Math.asinh((top - qs) / w), level);
  } else tot = ts(outer, 0, top, level);
  return 2 * delta * GF * mH ** 3 / (16 * SQ2 * PI) * tot;
}
if (typeof module !== "undefined") module.exports = { VVfast, tanhSinh };
