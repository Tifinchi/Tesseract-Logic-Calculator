# Appendix A calculations: Balance Chart completed, phase-circle wedges and e values,
# Fibonacci Option 2 (fast doubling), Calc Elements tiers, Root Elements typos.
import sympy as sp, math, itertools
s=sp.symbols('s',positive=True)
sizes=[sp.Integer(1),sp.sqrt(2),sp.sqrt(3),sp.Integer(2),sp.Integer(3)]
rows=[
 ("Perimeter, square","4s",4*s),
 ("Perimeter, half-square triangle (s, s, s√2)","s(2+√2)",s*(2+sp.sqrt(2))),
 ("Circumference, circle in the square (diameter s)","πs",sp.pi*s),
 ("Circumference, circle through the corners (diameter s√2)","π√2·s",sp.pi*sp.sqrt(2)*s),
 ("½ area: half-square S","s²/2",s**2/2),
 ("⅓ area: diagonal band of the 3×3 grid","s²/3",s**2/3),
 ("⅓ area: 120° sector of the inscribed circle","πs²/12",sp.pi*s**2/12),
 ("Cube surface","6s²",6*s**2),
 ("½ volume: bi-4 prism","s³/2",s**3/2),
 ("⅓ volume: inscribed tetrahedron (edge s√2)","s³/3",s**3/3),
 ("Inscribed sphere surface","πs²",sp.pi*s**2),
 ("Inscribed ball volume","πs³/6",sp.pi*s**3/6),
]
print("== Balance chart completion")
for name,rule,expr in rows:
    vals=[sp.radsimp(sp.expand(expr.subs(s,v))) for v in sizes]
    print(name,"|",rule,"|"," | ".join(str(v) for v in vals),"|",", ".join(f"{float(v):.4f}" for v in vals))
# check tetrahedron volume = 1/3 cube via coordinates
P=[sp.Matrix(p) for p in [(0,0,0),(1,1,0),(1,0,1),(0,1,1)]]
V=abs(sp.Matrix.hstack(P[1]-P[0],P[2]-P[0],P[3]-P[0]).det())/6
print("tetra volume in unit cube",V)
# centroid of tetrahedron 1/4 up from a face; triangle 1/3
print("== Phase circle r=3: hypotenuse as chord, segment area")
r=3
for name,c in [("initial",2*sp.sqrt(3)),("container",3*sp.sqrt(2)),("change",sp.sqrt(6))]:
    th=2*sp.asin(c/(2*r))
    seg=sp.Rational(1,2)*r**2*(th-sp.sin(th))
    print(name,"chord",c,"theta deg",float(th*180/sp.pi),"cos",sp.nsimplify(sp.cos(th)),"segment",sp.simplify(seg),float(seg))
print("== N = S - M wedge (container minus initial), legs a>=b")
for name,a,b in [("initial",3,sp.sqrt(3)),("container",3,3),("change",sp.sqrt(3),sp.sqrt(3))]:
    N=sp.Rational(1,2)*a*(a-b)
    print(name,sp.simplify(N),float(N))
print("== e under three definitions")
for name,a,b,written in [("initial",3,sp.sqrt(3),"√3/2"),("container",3,3,"0"),("change",sp.sqrt(3),sp.sqrt(3),"3")]:
    e_area=sp.sqrt(b*(a-b)); h=a-b; e_half=b/2
    d2=a**2-b**2
    tanA = (d2/e_area**2-1) if e_area!=0 else None
    print(name,"written",written,"| area e=",sp.nsimplify(e_area),float(e_area),"| leftover h=",sp.nsimplify(h),float(h),"| b/2=",e_half,"| tanA from (d,e)=",sp.nsimplify(tanA) if tanA is not None else None)
print("== 120 and 60 deg segments r=3")
for deg in (60,90,120):
    th=sp.pi*deg/180; seg=sp.Rational(1,2)*9*(th-sp.sin(th)); print(deg,sp.simplify(seg),float(seg))
print("== Fibonacci option 2 (fast doubling chart)")
n={1:1}
n[2]=n[1]**2+n[1]**2
steps=[]
for k in range(2,7):
    km2 = k-2 if k-2>=1 else 1
    d=n[k]**2-n[km2]**2; sm=n[k]**2+n[k-1]**2
    n[2*k-1]=d; n[2*k]=sm
    steps.append((k,km2,k-1,2*k-1,2*k,d,sm))
for st in steps: print(f"n{st[0]}²−n{st[1]}²=n{st[3]}={st[5]}; n{st[0]}²+n{st[2]}²=n{st[4]}={st[6]}")
print([n[i] for i in range(1,13)])
F=[0,1]
for i in range(80): F.append(F[-1]+F[-2])
assert all(n[i]==F[i+1] for i in range(1,13))
# identities
for k in range(1,30):
    assert F[k+1]**2-F[k-1]**2==F[2*k]
    assert F[k+1]**2+F[k]**2==F[2*k+1]
    assert sp.Rational(1,F[k]**2)+sp.Rational(1,F[k+1]**2)==sp.Rational(F[2*k+1],(F[k]*F[k+1])**2)
phi=(1+sp.sqrt(5))/2
print("line angle atan(1/phi) deg",float(sp.atan(1/phi)*180/sp.pi),"doubled",float(2*sp.atan(1/phi)*180/sp.pi),"atan2",float(sp.atan(2)*180/sp.pi))
print("golden angle deg",float(360*(1-1/phi)))
print("== Calc Elements checks")
def is_sq(n): r=int(math.isqrt(n)); return r*r==n
def sum2(n,cap=None):
    out=[]
    for a in range(1,int(math.isqrt(n))+1):
        b2=n-a*a
        if b2>=1 and is_sq(b2):
            b=math.isqrt(b2)
            if a>=b and (cap is None or a<=cap): out.append((a,b))
    return out
def diff2(n,cap=None):
    out=[]
    for a in range(1,(n+1)//2+2):
        for b in range(1,a):
            if a*a-b*b==n and (cap is None or a<=cap): out.append((a,b))
    return out
def sum3(n,cap=None):
    out=[]
    m=int(math.isqrt(n))
    for a in range(1,m+1):
        for b in range(1,a+1):
            c2=n-a*a-b*b
            if c2>=1 and is_sq(c2) and math.isqrt(c2)<=b and (cap is None or a<=cap): out.append((a,b,math.isqrt(c2)))
    return out
for Z,el in [(23,'V'),(27,'Co'),(31,'Ga'),(39,'Y'),(43,'Tc'),(44,'Ru'),(47,'Ag'),(49,'In'),(4,'Be')]:
    print(el,Z,"sum2",sum2(Z),"diff2",diff2(Z)[:3],"diff2 cap10",diff2(Z,10),"sum3 cap10",sum3(Z,10)[:3],"square",is_sq(Z))
print("Ga written (3,3,3,1)",9+9+9+1,"Ag written (5,4,3,1)",25+16+9+1)
duals=[Z for Z in range(1,55) if sum2(Z) and diff2(Z)]
print("both sum and difference, Z<=54:",duals)
print("Root Elements (9,8):",9*9+8*8,"; (3,1):",9+1,"= Ne")
print("== why 3: n = n(n-1)/2 solutions", [n for n in range(1,20) if n==n*(n-1)//2 and n*(n-1)%2==0])
print("Manhattan gap 3D unit cube:",(1+1+1)**2-3)
