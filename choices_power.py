# What each open definition does to the power of the frozen set (section 8.8, Appendix M.3).
# Same scoring as blindtest.py, Part 2 (Bayes factor against today's split normal spread by the
# future error; "pass" = within 2 combined sigma), vectorised with numpy. The baseline rows must
# reproduce Table 15 before the variants are read.
import numpy as np
from math import pi, sin, radians, log, exp

rng = np.random.default_rng(20261010)
NS = 20000
L = log(100)
s15 = sin(radians(15))**2

def lnN(x, m, s):
    return -0.5*((x - m)/s)**2 - np.log(s) - 0.5*np.log(2*pi)

def std_pred_logpdf(mu, sl, sh, sf):
    lo, hi = mu - 9*max(sl, sf), mu + 9*max(sh, sf)
    xs = np.linspace(lo, hi, 6001); h = xs[1] - xs[0]
    norm = 2/((2*pi)**0.5*(sl + sh))
    w = norm*np.exp(-0.5*((xs - mu)/np.where(xs < mu, sl, sh))**2)*h
    def f(m):
        m = np.atleast_1d(m)
        out = np.empty_like(m)
        for i in range(0, len(m), 2000):
            mm = m[i:i+2000, None]
            out[i:i+2000] = np.log(np.maximum(1e-300, (w*np.exp(-0.5*((mm - xs)/sf)**2)).sum(axis=1)/((2*pi)**0.5*sf)))
        return out
    return f

def sn_draw(mu, sl, sh, n):
    z = np.abs(rng.standard_normal(n))
    up = rng.random(n) < sh/(sl + sh)
    return np.where(up, mu + z*sh, mu - z*sl)

def score(c, sc, mu, sl, sh, sf):
    s_c = np.hypot(sc, sf)
    lp = std_pred_logpdf(mu, sl, sh, sf)
    m_r = rng.normal(c, s_c, NS); lb_c = lnN(m_r, c, s_c) - lp(m_r)
    m_w = sn_draw(mu, sl, sh, NS) + rng.normal(0, sf, NS); lb_s = lnN(m_w, c, s_c) - lp(m_w)
    p_pass = np.mean(np.abs(m_w - c) <= 2*s_c)
    side = sh if c > mu else sl
    pull = (mu - c)/np.hypot(side, sc)
    return dict(pull=pull, p=p_pass, Bpass=0.9545/max(p_pass, 1e-9), lb_c=lb_c, lb_s=lb_s)

rows = {
    "sin2 theta12": ((1 - s15)/(3*(1 - s15/3)), 0.0, 0.3096, 0.0073, 0.0057, 0.0016),
    "sin2 theta13": (s15/3, 0.0, 0.0220, 0.00056, 0.00056, 0.0004),
    "alpha": (90.0, 0.0, 84.1, 3.0, 3.7, 2.0),
    "Omega_m": (pi**2/32, 0.0, 0.3045, 0.0032, 0.0032, 0.003),
    "DG 10 GeV^2": (1/6, 0.0, 0.231, 0.047, 0.047, 0.023),
    "DSigma": (1/3, 0.0, 0.352, 0.029, 0.029, 0.050),
    "top, today = direct": (172.57, 0.38, 172.60, 0.27, 0.27, 0.035),
}
variants = {
    # Delta G at mu = 2 GeV: Zhao et al. evolved (section 5.3), 0.158 +- 0.039; future error 10% of the value
    "DG 2 GeV": (1/6, 0.0, 0.158, 0.039, 0.039, 0.016),
    # top, pole scheme against the PDG 2026 average of pole determinations
    "top, today = pole average": (172.57, 0.38, 172.1, 0.6, 0.6, 0.035),
    # top, direct average carried as pole-equivalent with the MC-pole difference 0 +- 0.52
    "top, today = direct as pole (+-0.52)": (172.57, 0.38, 172.60, 0.59, 0.59, 0.035),
}
res = {k: score(*v) for k, v in {**rows, **variants}.items()}
print("%-38s %7s %8s %7s %9s %9s" % ("prediction", "pull", "P(pass|", "B if", "med lnB", "med lnB"))
print("%-38s %7s %8s %7s %9s %9s" % ("", "today", "wrong)", "pass", "if right", "if wrong"))
for k, r in res.items():
    print("%-38s %+7.2f %8.3f %7.1f %9.2f %9.2f" % (k, r["pull"], r["p"], r["Bpass"], np.median(r["lb_c"]), np.median(r["lb_s"])))

def total(keys, label):
    c = sum(res[k]["lb_c"] for k in keys); s = sum(res[k]["lb_s"] for k in keys)
    print("%-62s if right: median B = %.3g, P(B > 100) = %.2f; if wrong: median B = %.2g, P(B < 1/100) = %.2f" % (
        label, exp(np.median(c)), np.mean(c > L), exp(np.median(s)), np.mean(s < -L)))

own = ["sin2 theta12", "sin2 theta13", "alpha", "Omega_m", "DG 10 GeV^2", "DSigma", "top, today = direct"]
print()
total(own, "the seven as in Table 15 (check: 300, 0.73; 1e-13, 0.93)")
total([k if k != "DG 10 GeV^2" else "DG 2 GeV" for k in own], "Delta G at 2 GeV")
total([k if k != "top, today = direct" else "top, today = pole average" for k in own], "top against the pole average")
total([{"DG 10 GeV^2": "DG 2 GeV", "top, today = direct": "top, today = pole average"}.get(k, k) for k in own], "both changed")
for dg in ("DG 10 GeV^2", "DG 2 GeV"):
    j = res["sin2 theta12"]["p"]*res["alpha"]["p"]*res[dg]["p"]
    print("passing theta12, alpha and %s together if wrong: %.4f (1 in %.0f)" % (dg, j, 1/j))
