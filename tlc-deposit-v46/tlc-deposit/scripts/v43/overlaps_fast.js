const F = require('./higgs_fast.js'); const PI = Math.PI, MPL = 1.22e19;
function overlaps(inp = {}, steps = 6000) {
  const Mt = inp.Mt ?? 172.57, h = Math.log(MPL / Mt) / steps, y = F.higgsInputs(inp), b = new Float64Array(5), q = new Float64Array(15);
  let y0 = y[0], y1 = y[1], y2 = y[2], y3 = y[3], y4 = y[4];
  const out = { M12: null, M13: null, M23: null };
  let p0 = 4 * PI / (y0 * y0), p1 = 4 * PI / (y1 * y1), p2 = 4 * PI / (y2 * y2);
  for (let i = 1; i <= steps; i++) {
    F.betaInto(y0, y1, y2, y3, y4, b); for (let j = 0; j < 5; j++) q[j] = b[j];
    F.betaInto(y0 + h/2*q[0], y1 + h/2*q[1], y2 + h/2*q[2], y3 + h/2*q[3], y4 + h/2*q[4], b); for (let j = 0; j < 5; j++) q[5 + j] = b[j];
    F.betaInto(y0 + h/2*q[5], y1 + h/2*q[6], y2 + h/2*q[7], y3 + h/2*q[8], y4 + h/2*q[9], b); for (let j = 0; j < 5; j++) q[10 + j] = b[j];
    F.betaInto(y0 + h*q[10], y1 + h*q[11], y2 + h*q[12], y3 + h*q[13], y4 + h*q[14], b);
    y0 += h/6*(q[0] + 2*q[5] + 2*q[10] + b[0]); y1 += h/6*(q[1] + 2*q[6] + 2*q[11] + b[1]); y2 += h/6*(q[2] + 2*q[7] + 2*q[12] + b[2]);
    y3 += h/6*(q[3] + 2*q[8] + 2*q[13] + b[3]); y4 += h/6*(q[4] + 2*q[9] + 2*q[14] + b[4]);
    const n0 = 4 * PI / (y0 * y0), n1 = 4 * PI / (y1 * y1), n2 = 4 * PI / (y2 * y2);
    const pr = [p0, p1, p2], nw = [n0, n1, n2];
    for (const [a, c, k] of [[0, 1, "M12"], [0, 2, "M13"], [1, 2, "M23"]]) if (out[k] === null) {
      const d0 = pr[a] - pr[c], d1 = nw[a] - nw[c];
      if (d0 * d1 <= 0 && d0 !== d1) { const f = d0 / (d0 - d1); out[k] = { mu: Mt * Math.exp((i - 1 + f) * h), ainv: pr[a] + f * (nw[a] - pr[a]) }; }
    }
    p0 = n0; p1 = n1; p2 = n2;
  }
  return out;
}
module.exports = { overlaps };
