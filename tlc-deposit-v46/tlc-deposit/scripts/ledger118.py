# Element ledger Z = 1..118: sum cells, difference cells, three-square cells, tier, standing on the 10x10 page.
import math, re, sys
EL="H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og".split()
assert len(EL)==118
isq=lambda n: math.isqrt(n)**2==n
def S(n): return [(a,math.isqrt(n-a*a)) for a in range(math.isqrt(n),0,-1) if n-a*a>=1 and isq(n-a*a) and a>=math.isqrt(n-a*a)]
def D(n): return sorted([(a,b) for a in range(1,n+1) for b in range(1,a) if a*a-b*b==n],key=lambda t:t[0])
def T3(n): return [(a,b,math.isqrt(n-a*a-b*b)) for a in range(math.isqrt(n),0,-1) for b in range(a,0,-1) if n-a*a-b*b>=1 and isq(n-a*a-b*b) and math.isqrt(n-a*a-b*b)<=b]
rows=[];cnt={};page={}
for Z in range(1,119):
    s=S(Z); d=D(Z); t=T3(Z)
    if isq(Z): tier=f"square ({math.isqrt(Z)}²)"; tk='square'
    elif s and d: tier='both'; tk='both'
    elif s: tier='sum'; tk='sum'
    elif d: tier='difference'; tk='difference'
    else: tier='three squares'; tk='three'
    cnt[tk]=cnt.get(tk,0)+1
    s10=[x for x in s if x[0]<=10]; d10=[x for x in d if x[0]<=10]; t10=[x for x in t if x[0]<=10]
    if isq(Z): pg='square'
    elif s10 and d10: pg='both'
    elif s10: pg='sum'
    elif d10: pg='difference'
    elif t10: pg='off the page (three legs)'
    else: pg='off the page (four legs)'
    page[pg]=page.get(pg,0)+1
    fs=", ".join(f"({a},{b})" for a,b in s) or "—"
    fd=", ".join(f"({a},{b})" for a,b in d[:3])+(" …" if len(d)>3 else "") if d else "—"
    ft=f"({t[0][0]},{t[0][1]},{t[0][2]})" if (not s and not d and not isq(Z)) else "—"
    rows.append((Z,EL[Z-1],fs,fd,ft,tier,pg))
print(cnt); print(page)
out=["| Z | Element | Sum cells a² + b² | Difference cells a² − b² | Three squares | Tier | On the 10 × 10 page |","| --- | --- | --- | --- | --- | --- | --- |"]
for r in rows: out.append("| "+" | ".join(str(x) for x in r)+" |")
open('ledger118.md','w').write("\n".join(out)+"\n")
# compare with HTUW ledger
src=open('/home/claude/obs/how-the-universe-works.md').read()
htuw={}
for line in src.splitlines():
    m=re.match(r'^\| (\d+) \| ([A-Z][a-z]?) \|(.*)\|$',line)
    if m and 1<=int(m.group(1))<=118 and m.group(2)==EL[int(m.group(1))-1]:
        htuw[int(m.group(1))]="| "+m.group(1)+" | "+m.group(2)+" |"+m.group(3)+"|"
diff=[(r[0],htuw.get(r[0])) for r,l in zip(rows,out[2:]) if htuw.get(r[0],'').replace(' ','')!=l.replace(' ','')]
print("rows differing from HTUW:",len(diff))
for z,h in diff[:12]:
    print("HTUW:",h); print("NEW :",out[1+z])
