# What the survivors of permute.py are worth: the chance of matching this well by accident across every permutation
# tried (look-elsewhere), the test on data the choice did not use (out of sample), the digits each survivor fixes
# beyond the measurement (precision), the time each takes (speed), and the cross-family checks. Writes perm.json
# for the calculator's permutation explorer.
import math, json, time, sys
import numpy as np
sys.path.insert(0, '/tmp/claude-0/-home-claude/bbcd54b6-3962-59fe-98a5-d2cef4b7e16c/scratchpad/reality/perm')
import permute as pm

out, t_rge = pm.run()
rng = np.random.default_rng(20261011)
NT = 20000

def matrices(P, keys):
    """predictions, their sigma, and a mask of which permutations reach each key"""
    n, k = len(P), len(keys)
    V = np.full((n, k), np.nan); SC = np.zeros((n, k)); M = np.zeros((n, k), bool)
    for i, p in enumerate(P):
        for j, key in enumerate(keys):
            if key == 'mh' and 'mh_pred' in p: V[i, j] = p['mh_pred']; SC[i, j] = 0.92; M[i, j] = True
            elif key in p['pred'] and key != 'mh':
                v, sc = p['pred'][key]
                if isinstance(v, (int, float)): V[i, j] = v; SC[i, j] = sc; M[i, j] = True
    return V, SC, M
def meas_arrays(keys):
    c, lo, hi = [], [], []
    for key in keys:
        m = pm.D[key]
        if len(m) == 2: c.append(m[0]); lo.append(m[1]); hi.append(m[1])
        else: c.append(m[0]); lo.append(m[1]); hi.append(m[2])
    return np.array(c), np.array(lo), np.array(hi)
def best_chi2(V, SC, M, centres, lo, hi, other_ok):
    """for each trial row of centres: the smallest chi2 among permutations whose every reached primary key is within 2 sigma"""
    T = centres.shape[0]; res = np.full(T, np.inf); ch = max(1, int(1.5e7/(V.shape[0]*V.shape[1])))
    for t0 in range(0, T, ch):
        c = centres[t0:t0 + ch][:, None, :]                       # (t, 1, k)
        s = np.where(V[None] > c, hi[None, None, :], lo[None, None, :])
        s = np.sqrt(s**2 + SC[None]**2)
        z = np.where(M[None], (c - V[None])/s, 0.0)
        ok = np.all(np.abs(z) <= 2.0, axis=2) & other_ok[None]
        chi = np.where(ok, np.sum(z*z, axis=2), np.inf)
        res[t0:t0 + ch] = chi.min(axis=1)
    return res

LEE = {
 'charged leptons': dict(keys=['mtau'], spread='rel', width=0.25),
 'lepton mixing': dict(keys=['s13', 's12'], spread='rel', width=0.25),
 'quark mixing angle': dict(keys=['alpha_all', 'legs_lat'], spread='rel', width=0.25),
 'Higgs vacuum load': dict(keys=['mt_dir'], spread='rel', width=0.05, dep=['mh', 'mt_pole']),   # the Higgs is the same relation read backwards
 "proton's spin legs": dict(keys=['dS_hermes', 'dG_lat'], spread='rel', width=0.25),
 "proton's radius": dict(keys=['rp_codata'], spread='rel', width=0.25),
 'matter share': dict(keys=['Om_sncb', 'q0_pop'], spread='share', width=0.25),
 'vacuum energy (hum)': dict(keys=['rhoL4'], spread='log', width=5.0),
}
summary = {}
print('=== look-elsewhere: how often some permutation matches random values as well as the chosen one matches the real ones')
for fam, cfg in LEE.items():
    P = out[fam]['P']; keys = cfg['keys']
    V, SC, M = matrices(P, keys)
    # keys outside the drawn set that a permutation must still pass (e.g. the electron anchor of the load, r_mu)
    dep = cfg.get('dep', [])     # measurements that would move with the drawn ones; not held at their real values
    other_ok = np.array([all(abs(v) <= 2 for k, v in p['pulls'].items() if k not in keys and k not in dep and k not in pm.SECONDARY and k not in out[fam]['inherited']) for p in P])
    c0, lo, hi = meas_arrays(keys)
    # the reference is the calculator's current value: its permutation if it survives, otherwise the surviving reading
    # that gives the same numbers on the drawn measurements (the value kept, the failed reading dropped)
    cur = [p for p in P if p['current'] and p['survives']]
    if not cur:
        cv = [p for p in P if p['current']]
        same = [p for p in P if p['survives'] and cv and all(k in p['pred'] and k in cv[0]['pred'] and p['pred'][k][0] == cv[0]['pred'][k][0] for k in keys if k in cv[0]['pred'] and k in p['pred'])
                and any(k in p['pred'] for k in keys)]
        cur = same
    ref = cur[0] if cur else min((p for p in P if p['survives']), key=lambda p: p['chi2_live'])
    chi_obs = sum(ref['pulls'][k]**2 for k in keys if k in ref['pulls'])
    nt = NT*10 if len(keys) == 1 else NT
    if cfg['spread'] == 'rel':
        C = c0[None] * (1 + cfg['width']*(2*rng.random((nt, len(keys))) - 1))
    elif cfg['spread'] == 'share':
        om = c0[0]*(1 + cfg['width']*(2*rng.random(nt) - 1)); C = np.stack([om, 1.5*om - 1], axis=1)
    else:
        C = 10**(math.log10(c0[0]) + cfg['width']*(rng.random((nt, 1)) - 0.5))
    if cfg['spread'] == 'log':
        # relative errors scale with the drawn value
        best = np.full(nt, np.inf)
        for t0 in range(0, nt, 500):
            cc = C[t0:t0 + 500]
            s = np.sqrt((0.0046*cc)**2 + SC[None, :, 0]**2)
            z = (cc - V[None, :, 0])/s
            ok = (np.abs(z) <= 2.0) & other_ok[None]
            best[t0:t0 + 500] = np.where(ok, z*z, np.inf).min(axis=1)
    else:
        best = best_chi2(V, SC, M, C, lo, hi, other_ok)
    p_acc = float(np.mean(best <= chi_obs + 1e-12))
    # the same with the window floored at 1 sigma per measurement, the convention section 16.8 uses for the standard theories
    p_fl = float(np.mean(best <= max(chi_obs, len(keys)) + 1e-12))
    nsurv = sum(p['survives'] for p in P)
    summary[fam] = dict(n=len(P), survivors=nsurv, ref=ref['rule'], chi_obs=chi_obs, p=p_acc, bits=(-math.log2(p_acc) if p_acc > 0 else float('inf')),
                        p_floor=p_fl, bits_floor=(-math.log2(p_fl) if p_fl > 0 else float('inf')), draws=int(nt))
    print(f"  {fam:22s} {len(P):6d} permutations, {nsurv:3d} survive; reference '{ref['rule']}', chi2 {chi_obs:.3f}; P(accident) = {p_acc:.4f}, {summary[fam]['bits']:.1f} bits; with a 1-sigma floor {p_fl:.4f}, {summary[fam]['bits_floor']:.1f} bits")

# lepton rule classes separately: the balance alone, the phase alone, the phase as Q/n
P = out['charged leptons']['P']
classes = {'Q alone (balance)': [p for p in P if p['nrules'] == 1 and p['Q'] is not None and p['delta'] is None],
           'delta alone (a fixed phase)': [p for p in P if p['nrules'] == 1 and p['Q'] is None and p['delta'] is not None],
           'delta = Q/n': [p for p in P if p['rule'].startswith('delta = Q/')]}
print('\n=== charged leptons by rule class (target m_tau drawn within 25%)')
for name, ps in classes.items():
    V, SC, M = matrices(ps, ['mtau']); c0, lo, hi = meas_arrays(['mtau'])
    C = c0[None]*(1 + 0.25*(2*rng.random((10*NT, 1)) - 1))
    surv = [p for p in ps if p['survives']]
    chi_obs = min((p['pulls']['mtau']**2 for p in surv), default=float('inf'))
    best = best_chi2(V, SC, M, C, lo, hi, np.ones(len(ps), bool))
    p_acc = float(np.mean(best <= chi_obs + 1e-12)); p_fl = float(np.mean(best <= max(chi_obs, 1.0) + 1e-12))
    print(f"  {name:30s} {len(ps):6d} permutations, {len(surv)} survive (best pull {math.sqrt(chi_obs):.3f}); P(accident) = {p_acc:.5f}, {(-math.log2(p_acc) if p_acc > 0 else float('inf')):.1f} bits; 1-sigma floor {p_fl:.5f}, {(-math.log2(p_fl) if p_fl > 0 else float('inf')):.1f} bits")
    summary['leptons: ' + name] = dict(n=len(ps), survivors=len(surv), p=p_acc, p_floor=p_fl)

# Fisher's method over the calculator's own families (the inherited lepton relation kept apart)
from scipy.stats import chi2 as _chi2
own = ['lepton mixing', 'quark mixing angle', 'Higgs vacuum load', "proton's spin legs", "proton's radius", 'matter share', 'vacuum energy (hum)']
for lab, key in (('as close as observed', 'p'), ('1-sigma floor', 'p_floor')):
    ps = [max(summary[f][key], 1.0/summary[f]['draws']) for f in own]
    X = -2*sum(math.log(x) for x in ps); pf = float(_chi2.sf(X, 2*len(ps)))
    summary['_fisher_' + key] = dict(p=pf, bits=-math.log2(pf), parts={f: summary[f][key] for f in own})
    print(f"  Fisher over the calculator's own seven families ({lab}): p = {pf:.4f}, {-math.log2(pf):.1f} bits")

# ------------------------------------------------------------------ out of sample: choose on one measurement, test on the others
print('\n=== out of sample: each family chooses on its training measurement; the test measurements were not used')
SPLITS = {
 'lepton mixing': (['s13'], ['s12', 's12_juno']),
 'quark mixing angle': (['alpha_all'], ['alpha_rr', 'alpha_bg']),
 'Higgs vacuum load': (['mt_dir'], ['mt_pole']),
 "proton's spin legs": (['dS_hermes'], ['dS_lat', 'dG_lat']),
 "proton's radius": (['rp_codata'], ['rp_muH', 'rp_ep']),
 'matter share': (['Om_sncb'], ['Om_planck', 'Om_desi', 'Om_dpa', 'q0_brout', 'OL_des', 'OL_pan', 'OL_uni']),
}
oos = {}
for fam, (train, test) in SPLITS.items():
    P = out[fam]['P']
    inh = out[fam]['inherited']
    elig = [p for p in P if all(k in p['pulls'] for k in train + test)
            and all(abs(v) <= 2 for k, v in p['pulls'].items() if k not in train + test and k not in pm.SECONDARY and k not in inh and k != 'mh')]
    def tr(p): return sum(p['pulls'][k]**2 for k in train)
    chosen = min(elig, key=lambda p: (tr(p), p['cx']))
    cur = [p for p in elig if p['current']] or [p for p in P if p['current']]
    row = {}
    for lab, p in (('chosen on training', chosen), ('current', cur[0] if cur else None)):
        if p is None: continue
        tp = {k: p['pulls'].get(k) for k in test}
        row[lab] = dict(rule=p['rule'], train={k: p['pulls'][k] for k in train}, test=tp, max_test=max(abs(v) for v in tp.values() if v is not None))
        print(f"  {fam:20s} {lab:19s} {p['rule'][:60]:60s} train {', '.join(f'{v:+.2f}' for v in row[lab]['train'].values())} | test {', '.join(f'{k} {v:+.2f}' for k, v in tp.items())}")
    oos[fam] = row

# ------------------------------------------------------------------ precision: digits each survivor fixes beyond the measurement
print('\n=== precision: the survivor\'s own uncertainty against the measurement\'s (standard theories: an input, no prediction)')
prec = []
def add_prec(name, value, s_pred, meas, unit):
    s_meas = meas[1] if len(meas) == 2 else 0.5*(meas[1] + meas[2])
    gain = (s_meas/s_pred) if s_pred > 0 else float('inf')
    prec.append(dict(name=name, value=value, s_pred=s_pred, s_meas=s_meas, unit=unit, gain=gain))
    g = 'exact' if s_pred == 0 else f'{gain:.3g}x'
    print(f"  {name:34s} {value:.10g} +- {s_pred:.3g} {unit} against +- {s_meas:.3g}: {g}")
# m_tau from Q = 2/3, propagated from m_e and m_mu/m_e (CODATA)
def tau_of(r):
    me = pm.ME; mm = me*r; s = math.sqrt(me) + math.sqrt(mm); t = me + mm; Qc = 2/3
    A, B, C = Qc - 1, 2*Qc*s, Qc*s*s - t; x = (-B - math.sqrt(B*B - 4*A*C))/(2*A); return x*x
mt0 = tau_of(pm.RMU); dmt = abs(tau_of(pm.RMU + pm.D['r_mu'][1]) - mt0)
add_prec('m_tau from Q = 2/3 (MeV)', mt0, dmt, pm.D['mtau'], 'MeV')
s15 = math.sin(math.radians(15))**2
add_prec('sin2 theta13, 15 deg tilt', s15/3, 0.0, pm.D['s13'], '')
add_prec('sin2 theta12, 15 deg tilt', (1 - s15)/(3*(1 - s15/3)), 0.0, pm.D['s12'], '')
add_prec('top pole mass, load (GeV)', pm.TOP_EDGE*(1 + [p for p in out['Higgs vacuum load']['P'] if p['current']][0]['delta']), 0.38, pm.D['mt_pole'], 'GeV')
add_prec('Delta Sigma = 1/3', 1/3, 0.0, pm.D['dS_hermes'], '')
add_prec('Delta G = 1/6 at 2 GeV', 1/6, 0.0, pm.D['dG_lat'], '')
add_prec('r_p / lambdabar_p = 4', 4.0, 0.0, pm.D['rp_codata'], '')
add_prec('Omega_m = pi^2/32', math.pi**2/32, 0.0, pm.D['Om_sncb'], '')
add_prec('vacuum energy scale (meV)', 2.017, 0.20, pm.D['rhoL4'], 'meV')

# ------------------------------------------------------------------ cross-family checks
print('\n=== cross-family checks')
L = out['Higgs vacuum load']['P']; s120 = [p for p in L if p['survives'] and ' at 120 deg' in p['rule']]
print(f"  phase spacing shared: the leptons' 120 deg fixes the load's angle; load survivors at 120 deg: {len(s120)} of {sum(p['survives'] for p in L)}: {[p['rule'] for p in s120]}")
N = out['neutrino masses']['P']; nq = [p for p in N if p['survives'] and abs(p['Q'] - 2/3) < 1e-12]
print(f"  balance shared: the charged leptons' Q = 2/3 fixes the neutrinos'; survivors: {[p['rule'] for p in nq]}")
H = out['vacuum energy (hum)']['P']; Sh = [p for p in out['matter share']['P'] if p['survives']][0]
for p in [p for p in H if p['survives']]:
    h0 = pm.H0_from(p['E'], 1 - Sh['Om'])
    pl = [(k, pm.pull(pm.D[k], h0, 0.2*h0)) for k in ('H0_cmb', 'H0_ladder', 'H0_jwst')]
    print(f"  hum {p['rule']:22s} E = {p['E']:.3f} meV -> H0 = {h0:.1f} km/s/Mpc with the share; " + ', '.join(f'{k} {v:+.2f}' for k, v in pl))

# ------------------------------------------------------------------ speed
print('\n=== speed: time per evaluation (this machine, Python)')
for fam, r in out.items():
    print(f"  {fam:22s} {1e6*r['t_eval']:9.1f} us per permutation")
print(f"  the vacuum's edge by the two-loop running (the standard method the load needs once): {t_rge/2:.3f} s per edge")

# ------------------------------------------------------------------ the explorer's data
def slim(p):
    return dict(rule=p['rule'], survives=bool(p['survives']), current=bool(p['current']), chi2=(round(p['chi2_live'], 4) if math.isfinite(p['chi2_live']) else None),
                pulls={k: (round(v, 3) if math.isfinite(v) else None) for k, v in p['pulls'].items()})
J = {}
for fam, r in out.items():
    P = r['P']
    keep = [p for p in P if p['survives'] or p['current']]
    rest = sorted([p for p in P if not (p['survives'] or p['current'])], key=lambda p: p['maxpull'])[:12]
    J[fam] = dict(n=len(P), survivors=sum(p['survives'] for p in P), inherited=r['inherited'], shown=[slim(p) for p in keep + rest],
                  lee=summary.get(fam), oos=oos.get(fam))
J['_meta'] = dict(date='2026-10-10 data; run 2026-10-11', threshold=2.0, product=float(math.prod(len(r['P']) for r in out.values())),
                  product_survivors=float(math.prod(max(1, sum(p['survives'] for p in r['P'])) for r in out.values())), precision=prec,
                  lepton_classes={k[9:]: v for k, v in summary.items() if k.startswith('leptons: ')},
                  fisher={k[8:]: v for k, v in summary.items() if k.startswith('_fisher_')})
def clean(x):
    if isinstance(x, dict): return {k: clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [clean(v) for v in x]
    if isinstance(x, (float, np.floating)): return float(x) if math.isfinite(x) else None
    if isinstance(x, (np.integer,)): return int(x)
    if isinstance(x, np.bool_): return bool(x)
    return x
json.dump(clean(J), open('/tmp/claude-0/-home-claude/bbcd54b6-3962-59fe-98a5-d2cef4b7e16c/scratchpad/reality/perm/perm.json', 'w'), separators=(',', ':'), allow_nan=False)
print('\nwrote perm.json')
