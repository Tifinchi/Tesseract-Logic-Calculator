# The look-elsewhere row for Delta G at its deposited alignment (2 GeV): alphabet members near 0.158 +- 0.039.
import importlib.util, io, contextlib
spec = importlib.util.spec_from_file_location("ab", "alphabet_counts.py"); ab = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()): spec.loader.exec_module(ab)
m, s = 0.158, 0.039
w2 = [v for v in ab.merged if m - 2*s <= v <= m + 2*s]; w3 = [v for v in ab.merged if m - 3*s <= v <= m + 3*s]
near = min(ab.merged, key=lambda v: abs(v - m))
print("Delta G at 2 GeV, 0.158 +- 0.039: within 2 sigma %d, within 3 sigma %d, nearest %.6f, pull %+.2f; 1/6 pull %+.2f" % (len(w2), len(w3), near, (m - near)/s, (m - 1/6)/s))
