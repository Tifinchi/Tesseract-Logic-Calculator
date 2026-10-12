// The calculator's numerics, Version 42 against Version 43: time per call (median of repeated runs after warm-up) and
// error against converged references (50-digit quadrature for motion, a 40,000-step run for the edge, an adaptive
// double integral for the off-shell widths). Run with Node 22 (V8, the engine of Chromium-based browsers).
const R = require('../v41/motion_engine.js'), X = require('./motion_exact.js');
const HO = require('./higgs_old.js'), HN = require('./higgs_fast.js');
const OO = require('./overlaps_old.js'), ON = require('./overlaps_fast.js');
const { VVold } = require('./vv_old.js'), { VVfast } = require('./vv_fast.js');
function time(f, minMs = 300) {
  for (let i = 0; i < 3; i++) f();
  const ts = []; let tot = 0, r;
  while (tot < minMs || ts.length < 5) { const t0 = process.hrtime.bigint(); r = f(); const dt = Number(process.hrtime.bigint() - t0) / 1e6; ts.push(dt); tot += dt; if (ts.length > 2000) break; }
  ts.sort((a, b) => a - b); return [r, ts[ts.length >> 1]];
}
const REF = { merc_both: 5.0186647268741304207e-7, merc_time: 1.672888030253721109397308e-7, e95_both: 1.9082857237397057203e-5, b1_both: 1.7512012730733550276, b1_time: 0.87559808216236236465,
              WW: 0.0008549217818730896, ZZ: 0.0001066297296194355 };
const rows = [];
const rel = (v, r) => Math.abs(v - r) / Math.abs(r);
function row(name, fo, fn, ref, unit = 'rel') {
  const [vo, to] = time(fo), [vn, tn] = time(fn);
  const eo = ref == null ? null : (unit === 'rel' ? rel(vo, ref) : Math.abs(vo - ref)), en = ref == null ? null : (unit === 'rel' ? rel(vn, ref) : Math.abs(vn - ref));
  rows.push({ name, old_ms: to, new_ms: tn, speedup: to / tn, old_err: eo, new_err: en, unit });
}
row("Mercury's perihelion advance, both legs", () => R.advance(0.387098, 0.205630, 'both'), () => X.advance(0.387098, 0.205630, 'both'), REF.merc_both);
row("the same, the time leg alone", () => R.advance(0.387098, 0.205630, 'time'), () => X.advance(0.387098, 0.205630, 'time'), REF.merc_time);
row("an orbit with e = 0.95, both legs", () => R.advance(0.1, 0.95, 'both'), () => X.advance(0.1, 0.95, 'both'), REF.e95_both);
row("light at the Sun's edge, both legs", () => R.bending(1, 1, 'both'), () => X.bending(1, 1, 'both'), REF.b1_both);
row("light at the Sun's edge, the time leg alone", () => R.bending(1, 1, 'time'), () => X.bending(1, 1, 'time'), REF.b1_time);
const inp = { Mh: 125.13, als: 0.1180, MW: 80.3625 }, inpH = { Mt: 172.60, als: 0.1180, MW: 80.3625 };
const eRef = HO.bisect(Mt => HO.higgsRun({ ...inp, Mt }, 40000).lmin, 160, 185, 1e-10), hRef = HO.bisect(Mh => HO.higgsRun({ ...inpH, Mh }, 40000).lmin, 110, 150, 1e-10);
row("the vacuum's edge for the top (GeV)", () => HO.edgeTop(inp), () => HN.edgeTop(inp), eRef, 'abs');
row("the vacuum's edge for the Higgs (GeV)", () => HO.edgeHiggs(inpH), () => HN.edgeHiggs(inpH), hRef, 'abs');
row("the edge's error budget, top", () => HO.edgeErrors(inp, 'top').total, () => HN.edgeErrors(inp, 'top').total, null);
row("the edge's error budget, Higgs", () => HO.edgeErrors(inpH, 'higgs').total, () => HN.edgeErrors(inpH, 'higgs').total, null);
row("the running for display, 3,000 steps", () => HO.higgsRun({ ...inp, Mt: 172.6 }, 3000).lmin, () => HN.higgsRun({ ...inp, Mt: 172.6 }, 3000).lmin, HO.higgsRun({ ...inp, Mt: 172.6 }, 3000).lmin);
row("where the forces cross, 6,000 steps", () => OO.overlaps(inp).M23.mu, () => ON.overlaps(inp).M23.mu, OO.overlaps(inp).M23.mu);
const GF = 1.1663788e-5;
row("H → WW*", () => VVold(125.13, 80.3692, 2.085, 2, GF), () => VVfast(125.13, 80.3692, 2.085, 2, GF, 4), REF.WW);
row("H → ZZ*", () => VVold(125.13, 91.1880, 2.4955, 1, GF), () => VVfast(125.13, 91.1880, 2.4955, 1, GF, 4), REF.ZZ);
const f = x => x == null ? '' : x.toExponential(1);
console.log('computation | V42 ms | V43 ms | speed-up | V42 error | V43 error');
for (const r of rows) console.log(`${r.name} | ${r.old_ms.toPrecision(3)} | ${r.new_ms.toPrecision(3)} | ${r.speedup.toFixed(0)} | ${f(r.old_err)} | ${f(r.new_err)}${r.unit === 'abs' ? ' (GeV)' : ''}`);
console.log(JSON.stringify(rows));
