const H = require('./higgs_old.js'); const PI = Math.PI, MPL = 1.22e19; const beta = H.beta, higgsInputs = H.higgsInputs;
  function overlaps(inp = {}, steps = 6000) {
    const Mt = inp.Mt ?? 172.57, h = Math.log(MPL / Mt) / steps;
    let y = higgsInputs(inp);
    const add = (u, v, s) => u.map((x, i) => x + s * v[i]), ai = y => y.slice(0, 3).map(g => 4 * PI / (g * g));
    const pairs = [[0, 1, "M12"], [0, 2, "M13"], [1, 2, "M23"]], out = { M12: null, M13: null, M23: null };
    let prev = ai(y);
    for (let i = 1; i <= steps; i++) {
      const q1 = beta(y), q2 = beta(add(y, q1, h/2)), q3 = beta(add(y, q2, h/2)), q4 = beta(add(y, q3, h));
      y = y.map((x, j) => x + h/6 * (q1[j] + 2*q2[j] + 2*q3[j] + q4[j]));
      const now = ai(y);
      for (const [a, b, k] of pairs) if (out[k] === null) {
        const d0 = prev[a] - prev[b], d1 = now[a] - now[b];
        if (d0 * d1 <= 0 && d0 !== d1) { const f = d0 / (d0 - d1); out[k] = { mu: Mt * Math.exp((i - 1 + f) * h), ainv: prev[a] + f * (now[a] - prev[a]) }; }
      }
      prev = now;
    }
    return out;
  }

module.exports = { overlaps };
