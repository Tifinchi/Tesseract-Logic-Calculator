# Root-family tables (legs 1-10) in the author's factored form, with element symbols for N <= 118,
# and the element tier tally (legs unbounded and on the 10 x 10 page).
from sympy import factorint
SYM = ("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr "
       "Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt "
       "Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og").split()
assert len(SYM) == 118
sup = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
def fac(n):
    if n == 0: return "0"
    f = factorint(n)
    if n == 1: return "1"
    return "·".join(("%d" % p) + (str(e).translate(sup) if e > 1 else "") for p, e in sorted(f.items()))
def cell(n):
    s = "√%d" % n if n != fac(n) or True else ""
    body = "√(%s)" % fac(n) if n > 1 else ("√1" if n == 1 else "0")
    if 1 <= n <= 118: body += " " + SYM[n-1]
    return body
print("SUM")
print("| a \\ b | " + " | ".join(str(b) for b in range(1, 11)) + " |")
print("| --- " * 11 + "|")
for a in range(1, 11):
    row = [cell(a*a + b*b) if b <= a else "" for b in range(1, 11)]
    print("| **%d** | " % a + " | ".join(row) + " |")
print("DIFF")
print("| a \\ b | " + " | ".join(str(b) for b in range(1, 11)) + " |")
print("| --- " * 11 + "|")
for a in range(1, 11):
    row = [cell(a*a - b*b) if b <= a else "" for b in range(1, 11)]
    print("| **%d** | " % a + " | ".join(row) + " |")
# tiers
import math
def is_sq(n): r = math.isqrt(n); return r*r == n
def sum2(n, cap=None):
    out = []
    for a in range(1, math.isqrt(n) + 1):
        b2 = n - a*a
        if b2 >= 1 and is_sq(b2):
            b = math.isqrt(b2)
            if b <= a and (cap is None or a <= cap): out.append((a, b))
    return out
def diff2(n, cap=None):
    out = []
    for b in range(1, n):
        a2 = n + b*b
        if is_sq(a2):
            a = math.isqrt(a2)
            if cap is None or a <= cap: out.append((a, b))
        if b > n: break
    return out
tally = {}; page = {}; three = []
for Z in range(1, 119):
    s, d = sum2(Z), diff2(Z)
    if is_sq(Z): t = "square"
    elif s and d: t = "both"
    elif s: t = "sum"
    elif d: t = "difference"
    else: t = "three"
    tally.setdefault(t, []).append(SYM[Z-1])
    sp, dp = sum2(Z, 10), diff2(Z, 10)
    if is_sq(Z) and Z <= 100: tp = "square"
    elif sp and dp: tp = "both"
    elif sp: tp = "sum"
    elif dp: tp = "difference"
    else:
        # third leg <= 10 on top of a page sum cell, or with a page difference cell
        ok = any((Z - c*c) > 0 and sum2(Z - c*c, 10) for c in range(1, 11)) or any((Z - c*c) > 0 and diff2(Z - c*c, 10) for c in range(1, 11)) or any(diff2(Z + c*c, 10) for c in range(1, 11))
        tp = "third leg" if ok else "four legs"
    page.setdefault(tp, []).append(SYM[Z-1])
for k, v in tally.items(): print(k, len(v), " ".join(v))
print("page:")
for k, v in page.items(): print(k, len(v), " ".join(v))
