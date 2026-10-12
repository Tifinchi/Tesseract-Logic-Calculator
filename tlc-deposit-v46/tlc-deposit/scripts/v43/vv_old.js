const PI = Math.PI, SQ2 = Math.SQRT2;
function make(mH, GF) {
    function VV(mV, GV, delta, n = 320) {
      const dq_ = (mH - 2e-4) / (n - 1), bw = s => (mV * GV / PI) / ((s - mV * mV) ** 2 + (mV * GV) ** 2); let tot = 0;
      for (let i = 0; i < n; i++) { const q1 = 1e-4 + i * dq_, x = q1 * q1 / (mH * mH), w1 = bw(q1 * q1) * 2 * q1;
        for (let j = 0; j < n; j++) { const q2 = 1e-4 + j * dq_; if (q1 + q2 > mH) break;
          const y = q2 * q2 / (mH * mH), lam = Math.max(0, (1 - x - y) ** 2 - 4 * x * y);
          tot += delta * GF * mH ** 3 / (16 * SQ2 * PI) * Math.sqrt(lam) * (lam + 12 * x * y) * w1 * bw(q2 * q2) * 2 * q2; } }
      return tot * dq_ * dq_;
    }
  return VV;
}
module.exports = { VVold: (mH, mV, GV, delta, GF) => make(mH, GF)(mV, GV, delta) };
