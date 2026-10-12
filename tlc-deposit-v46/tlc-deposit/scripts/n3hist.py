# Reproduce the N3 null distribution of blindtest.py (same seed, same order) and print its histogram.
import collections
src = open('blindtest.py').read().split('\n')
setup = '\n'.join(src[:105])            # lines 1-105: alphabets, targets, cube hits, null generators, NR
import io, contextlib
buf = io.StringIO()
g = {}
with contextlib.redirect_stdout(buf):
    exec(setup, g)
h2 = []
for _ in range(g['NR']):
    A = g['null_N3'](); h2.append(g['hits'](A, g['inr'], 2))
c = collections.Counter(h2)
print("N3 alphabets:", g['NR'], "mean", sum(h2)/len(h2))
for k in sorted(c): print(k, c[k])
print("P(>=19) =", sum(1 for h in h2 if h >= 19)/len(h2))
