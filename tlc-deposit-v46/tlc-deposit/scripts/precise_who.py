# Which precise targets do random-ingredient alphabets hit? (same recipe and seed as decide.py, static alphabets)
import collections, io, contextlib
src = open('decide.py').read()
head = src.split('print("=== Test 3')[0]
g = {}
exec(head, g)
import random
random.seed(20261010 + 5)
cnt = collections.Counter(); NR = 400
for _ in range(NR):
    A = g['build'](g['base_from'](random.uniform(2.6, 3.7), random.uniform(1.25, 1.6), random.uniform(1.55, 1.95), random.uniform(1.45, 1.8)))
    for name in g['hits'](A, g['PRECISE']): cnt[name] += 1
for k, v in cnt.most_common(): print("%-16s hit by %d of %d" % (k, v, NR))
