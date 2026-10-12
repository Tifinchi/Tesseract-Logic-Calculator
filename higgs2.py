"""Standard Model Higgs self-coupling running: 2-loop (full, top Yukawa) + 3-loop QCD for g3 and yt.
Convention: V = lambda (H^dag H - v^2/2)^2, m_h^2 = 2 lambda v^2. g1 GUT-normalized (g1^2 = 5/3 g'^2). t = ln(mu/Mt)."""
import math, numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
k=1/(16*math.pi**2); z3=1.2020569031595942
MPl=1.22e19
def beta(t,y,loops=2,three=True):
    g1,g2,g3,yt,l=y
    a1,a2,a3,T=g1*g1,g2*g2,g3*g3,yt*yt
    gp=0.6*a1  # g'^2
    b=[k*41/10*g1**3, k*(-19/6)*g2**3, k*(-7)*g3**3,
       k*yt*(4.5*T-8*a3-2.25*a2-0.85*a1),
       k*(24*l*l-6*T*T+0.375*(2*a2*a2+(a2+gp)**2)+l*(12*T-9*a2-3*gp))]
    if loops>=2:
        K=k*k
        b[0]+=K*g1**3*(199/50*a1+27/10*a2+44/5*a3-17/10*T)
        b[1]+=K*g2**3*(9/10*a1+35/6*a2+12*a3-1.5*T)
        b[2]+=K*g3**3*(11/10*a1+4.5*a2-26*a3-2*T)
        b[3]+=K*yt*(-12*T*T+T*(131/16*gp+225/16*a2+36*a3-12*l)+6*l*l
                    +1187/216*gp*gp-0.75*gp*a2+19/9*gp*a3-23/4*a2*a2+9*a2*a3-108*a3*a3)
        b[4]+=K*(-312*l**3-144*l*l*T+36*l*l*(3*a2+gp)-3*l*T*T+l*T*(80*a3+22.5*a2+85/6*gp)
                 -73/8*l*a2*a2+39/4*l*a2*gp+629/24*l*gp*gp
                 +30*T**3-32*T*T*a3-8/3*T*T*gp-9/4*T*a2*a2+10.5*T*a2*gp-19/4*T*gp*gp
                 +305/16*a2**3-289/48*a2*a2*gp-559/48*a2*gp*gp-379/48*gp**3)
    if three:
        K3=k**3
        b[2]+=K3*g3**3*(65/2*a3*a3-40*a3*T+15*T*T)
        b[3]+=K3*yt*(-4166/3+640*z3)*a3**3
    return b
def inputs(Mh=125.20,Mt=172.57,als=0.1180,MW=80.3692,dl=0.0,dy=0.0):
    a=(als-0.1184)/0.0007; w=(MW-80.384)/0.014; d=Mt-173.34
    gY=0.35830+0.00011*d-0.00020*w
    g2=0.64779+0.00004*d+0.00011*w
    g3=1.1666+0.00314*a-0.00046*d
    yt=0.93690+0.00556*d-0.00042*a+dy
    l=0.12604+0.00206*(Mh-125.15)-0.00004*d+dl
    return [math.sqrt(5/3)*gY,g2,g3,yt,l]
def run(Mt=172.57,top=MPl,loops=2,three=True,**kw):
    y0=inputs(Mt=Mt,**kw)
    return solve_ivp(lambda t,y:beta(t,y,loops,three),[0,math.log(top/Mt)],y0,method='DOP853',rtol=1e-11,atol=1e-14,dense_output=True)
def lam_profile(s,Mt,top=MPl,n=6000):
    ts=np.linspace(0,math.log(top/Mt),n); return ts,s.sol(ts)[4]
def summary(Mt=172.57,loops=2,three=True,top=MPl,**kw):
    s=run(Mt,top,loops,three,**kw); ts,l=lam_profile(s,Mt,top)
    i=int(np.argmin(l)); lmin=l[i]; mumin=Mt*math.exp(ts[i])
    neg=np.where(l<0)[0]
    cross=Mt*math.exp(brentq(lambda t:s.sol(t)[4],ts[neg[0]-1],ts[neg[0]])) if len(neg) else None
    yPl=s.sol(math.log(MPl/Mt)) if top>=MPl else None
    return dict(s=s,lmin=lmin,mumin=mumin,cross=cross,yPl=yPl,lPl=yPl[4] if yPl is not None else None)
