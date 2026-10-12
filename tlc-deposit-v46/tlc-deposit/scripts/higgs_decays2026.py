"""Calculate the Higgs on the board.
MASS (formulas fixed before computing):
  M1 balance (hers): the real Higgs sits just past its exact edge -> below m_H,edge(top). Side only.
  M2 the 4D corner rule (held at sqrt4, centre-to-corner = edge): lambda = 1/8 -> m_H = v/2.
  M3 semiconductive/held: the Higgs sits below the pair gap 2 m_W of its strongest boson reaction. Side only.
DECAYS as reactions of the board: coupling to each reaction = its gap / bias (y = sqrt2 m / v);
open channels (pair above its gap), held channels (one partner virtual: WW*, ZZ*), sub-reaction channels (loops: gg, gamma gamma)."""
import math, numpy as np
GF=1.1663788e-5; v=(math.sqrt(2)*GF)**-0.5; mW=80.3625; mZ=91.1880; GW=2.085; GZ=2.4955; mt=172.60; mtau=1.77686; mmu=0.1056584
mbb=4.183; mcc=1.2730; asZ=0.1180; alpha0=1/137.035999
print(f"bias of the board v = (sqrt2 G_F)^-1/2 = {v:.3f} GeV")
print("MASS")
print("  M1 balance: m_H below its exact edge 128.47 +- 1.01 GeV (top 172.56) -> measured 125.20 +- 0.11: predicted side, 45.74 deg")
m2=v/2; print(f"  M2 lambda = 1/8 -> m_H = v/2 = {m2:.2f} GeV vs 125.20 +- 0.11: {100*(m2/125.20-1):+.2f}%, pull {(125.20-m2)/0.11:+.1f} sigma -> MISS")
print(f"  M3 held below its W-pair gap: 2 m_W = {2*mW:.2f} GeV > 125.20 -> predicted side (m_H/2m_W = {125.20/(2*mW):.4f})")
print(f"  measured lambda (tree) = m_H^2/(2v^2) = {125.20**2/(2*v*v):.5f} (1/8 = 0.125)")
# alpha_s and running masses, 2-loop, a = alpha_s/pi
Z3=1.2020569031595942
def b(nf): return (11-2*nf/3)/4, (102-38*nf/3)/16, (2857/2-5033/18*nf+325/54*nf*nf)/64
def g(nf): return 1.0, (202/3-20*nf/9)/16, (1249+(-2216/27-160/3*Z3)*nf-140/81*nf*nf)/64
def run_a(a,mu0,mu1,nf,steps=4000):
    b0,b1,b2=b(nf); L=math.log(mu1**2/mu0**2); h=L/steps
    f=lambda a:-b0*a*a-b1*a**3-b2*a**4
    for _ in range(steps):
        k1=f(a);k2=f(a+h*k1/2);k3=f(a+h*k2/2);k4=f(a+h*k3);a+=h*(k1+2*k2+2*k3+k4)/6
    return a
def run_m(m,a,mu0,mu1,nf,steps=4000):
    b0,b1,b2=b(nf); g0,g1,g2=g(nf); L=math.log(mu1**2/mu0**2); h=L/steps
    y=np.array([a,math.log(m)])
    f=lambda y:np.array([-b0*y[0]**2-b1*y[0]**3-b2*y[0]**4, -(g0*y[0]+g1*y[0]**2+g2*y[0]**3)])
    for _ in range(steps):
        k1=f(y);k2=f(y+h*k1/2);k3=f(y+h*k2/2);k4=f(y+h*k3);y=y+h*(k1+2*k2+2*k3+k4)/6
    return math.exp(y[1]),y[0]
def decays(mH,verbose=False):
    aZ=asZ/math.pi; aH=run_a(aZ,mZ,mH,5); asH=aH*math.pi
    ab=run_a(aZ,mZ,mbb,5); mb,_=run_m(mbb,ab,mbb,mH,5)
    ac=run_a(ab,mbb,mcc,4); mc_b,_=run_m(mcc,ac,mcc,mbb,4); mc,_=run_m(mc_b,ab,mbb,mH,5)
    pre=GF*mH/(4*math.sqrt(2)*math.pi)
    dq=1+5.67*aH+(35.94-1.36*5)*aH**2+(164.14-25.77*5+0.259*25)*aH**3
    W={}
    W['bb']=3*pre*mb**2*dq; W['cc']=3*pre*mc**2*dq
    for k,m in (('tautau',mtau),('mumu',mmu)): W[k]=pre*m*m*(1-4*m*m/mH**2)**1.5
    # held channels: both vectors off shell, Breit-Wigner weights
    def VV(mV,GV,delta,n=500):
        q=np.linspace(1e-4,mH-1e-4,n); dq_=q[1]-q[0]
        bw=lambda s: (mV*GV/math.pi)/((s-mV**2)**2+(mV*GV)**2)
        tot=0.0
        Q1,Q2=np.meshgrid(q,q,indexing='ij'); m1=Q1<=mH-Q2
        x=Q1**2/mH**2; y=Q2**2/mH**2; lam=np.clip((1-x-y)**2-4*x*y,0,None)
        g0=delta*GF*mH**3/(16*math.sqrt(2)*math.pi)*np.sqrt(lam)*(lam+12*x*y)
        w=bw(Q1**2)*2*Q1*bw(Q2**2)*2*Q2
        return float(np.sum(np.where(m1,g0*w,0))*dq_*dq_)
    W['WW']=VV(mW,GW,2.0); W['ZZ']=VV(mZ,GZ,1.0)   # delta_Z = 1 already holds the identical-particle 1/2
    # sub-reaction (loop) channels
    def f(t): return math.asin(math.sqrt(t))**2 if t<=1 else -0.25*(math.log((1+math.sqrt(1-1/t))/(1-math.sqrt(1-1/t)))-1j*math.pi)**2
    A12=lambda t: 2*(t+(t-1)*f(t))/t**2
    A1=lambda t: -(2*t*t+3*t+3*(2*t-1)*f(t))/t**2
    tt=mH**2/(4*mt**2); tb=mH**2/(4*4.78**2); tW=mH**2/(4*mW**2)
    Ag=0.75*(A12(tt)+A12(tb)); K=1+(95/4-7*5/6)*aH+(156.81-5.7083*5)*aH**2
    W['gg']=GF*asH**2*mH**3/(36*math.sqrt(2)*math.pi**3)*abs(Ag)**2*K
    Aa=A1(tW)+3*(4/9)*A12(tt)+3*(1/9)*A12(tb)
    W['gamgam']=GF*alpha0**2*mH**3/(128*math.sqrt(2)*math.pi**3)*abs(Aa)**2
    tot=sum(W.values())
    if verbose: print(f"  alpha_s(m_H) = {asH:.4f}; m_b(m_H) = {mb:.3f} GeV, m_c(m_H) = {mc:.3f} GeV (3-loop running); gg QCD factor {K:.3f}")
    return W,tot

ref={'bb':58.2,'WW':21.4,'gg':8.18,'tautau':6.27,'cc':2.89,'ZZ':2.62,'gamgam':0.227,'mumu':0.0218}
for mH in (125.00,125.13):
    W,tot=decays(mH,True)
    print(f"m_H={mH}: total width {tot*1e3:.3f} MeV")
    for k in ref:
        print(f"  {k:8s} width {W[k]*1e6:9.2f} keV  BR {100*W[k]/tot:8.4f}%  YR4(125.0) {ref[k]:7.4f}%  ratio {100*W[k]/tot/ref[k]:.3f}")
