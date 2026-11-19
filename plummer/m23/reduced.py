"""Reduced line-piece model for Mathon's two-graph (exact; validated against the full model on q = 7).
Facts checked numerically from S: (i) S_BB preserves slope classes; (ii) S_BC maps the mean-zero slope-m2 part of block C
into a mean-zero slope-m1 function of block B for exactly one (B, m1) per (C, m2), and to 0 otherwise; (iii) means
couple through sigma (Paley(p) Seidel). Then x is an eigenvector iff (a) mean equations hold and (b) each block with a
non-constant piece of slope m has its partner block using the partner slope, with an allowed pair of types.
python reduced.py S.npy q K TLIM WORKERS MAXSOL OUT"""
import sys, time, itertools, numpy as np
from ortools.sat.python import cp_model
S=np.load(sys.argv[1]).astype(np.float64); q=int(sys.argv[2]); K=int(sys.argv[3]); TL=float(sys.argv[4]); W=int(sys.argv[5]); MAXSOL=int(sys.argv[6]); OUT=sys.argv[7]
N=len(S); b=q*q; p=(N-1)//b; rho=N-1-2*K
Qs={(x*x)%q for x in range(1,q)}; Ns=set(range(1,q))-Qs
TS={"0":{0},"Q":Qs,"N":Ns,"0Q":{0}|Qs,"0N":{0}|Ns,"*":Qs|Ns}; TN=list(TS)
slopes=list(range(q))+["inf"]
def form(m):
    return (lambda a,bb:a%q) if m=="inf" else (lambda a,bb,m=m:(bb-m*a)%q)
IND={(m,t):np.array([1. if form(m)(a,bb) in TS[t] else 0. for a in range(q) for bb in range(q)]) for m in slopes for t in TN}
def slope_of(g):
    for m in slopes:
        f=form(m); lab=np.array([f(a,bb) for a in range(q) for bb in range(q)])
        ok=all(np.ptp(g[lab==c])<1e-6 for c in range(q))
        if ok: return m
    return None
blk=lambda B: slice(B*b,(B+1)*b)
t0=time.time()
# (i) diagonal preserves slopes, (ii) partner map
partner={}
for C in range(p):
    for m2 in slopes:
        f=IND[(m2,"Q")]-IND[(m2,"Q")].mean()
        for B in range(p):
            g=S[blk(B),blk(C)]@f; g=g-g.mean()
            if np.abs(g).max()<1e-6: continue
            m1=slope_of(g)
            assert m1 is not None, ("non-line image",B,C,m2)
            if B==C: assert m1==m2, ("diag changes slope",B,m2,m1); continue
            assert (C,m2) not in partner, ("two partners",C,m2)
            partner[(C,m2)]=(B,m1)
assert len(partner)==p*len(slopes), len(partner)
# (iii) allowed type pairs
allowed={}; CHECKALL=[20]
for (C,m2),(B,m1) in partner.items():
    if (B,m1,C,m2) in allowed: continue
    ok=set()
    for t1 in TN:
        for t2 in TN:
            yB=-2*(IND[(m1,t1)]-IND[(m1,t1)].mean()); yC=-2*(IND[(m2,t2)]-IND[(m2,t2)].mean())
            zB=S[blk(B),blk(B)]@yB+S[blk(B),blk(C)]@yC-rho*yB; zC=S[blk(C),blk(C)]@yC+S[blk(C),blk(B)]@yB-rho*yC
            bad=np.abs(zB-zB.mean()).max()>1e-6 or np.abs(zC-zC.mean()).max()>1e-6
            if not bad and CHECKALL[0]>0:
                CHECKALL[0]-=1; y=np.zeros(N); y[blk(B)]=yB; y[blk(C)]=yC; z=S@y-rho*y
                for D in range(p):
                    zz=z[blk(D)]; assert np.abs(zz-zz.mean()).max()<1e-6, "leak to other block"
            if not bad: ok.add((t1,t2))
    allowed[(B,m1,C,m2)]=ok; allowed[(C,m2,B,m1)]={(u,v) for (v,u) in ok}
sig=np.array([[S[B*b,C*b+0] if False else 0 for C in range(p)] for B in range(p)])
# Paley Seidel sign of block pair: mean of S_BC entries * q / (-q) ... derive from block sums: sum(S_BC) = sigma_BC * (-q) * q
for B in range(p):
    for C in range(p):
        if B!=C: sig[B,C]=int(round(S[blk(B),blk(C)].sum()/(-q**3)))
print(f"q={q} p={p} rho={rho}: partners ok, allowed pair counts {sorted(set(len(v) for v in allowed.values()))}, sigma entries {sorted(set(sig.flatten()))} ({time.time()-t0:.0f}s)",flush=True)
# options
OPTS=[("-",None,None,0),("all",None,None,b)]+[(f"{m}:{t}",m,t,q*len(TS[t])) for m in slopes for t in TN]
nopt=len(OPTS); m_=cp_model.CpModel(); z=[[m_.NewBoolVar("") for o in range(nopt)] for B in range(p)]
for B in range(p): m_.AddExactlyOne(z[B])
w=[sum(OPTS[o][3]*z[B][o] for o in range(nopt)) for B in range(p)]
# mean eq: -q sum_{C!=B} sig_BC xbar_C + 1 = rho xbar_B, xbar = 1 - 2w/q^2 ; multiply by q^2
for B in range(p):
    m_.Add(-q*sum(int(sig[B,C])*(b-2*w[C]) for C in range(p) if C!=B) + b == rho*(b-2*w[B]))
m_.Add(sum(b-2*w[C] for C in range(p))*1 == rho*q*0 + (rho*b)//b*0 + int(round(rho*b/b)) if False else sum(b-2*w[C] for C in range(p))==rho)  # oo row: sum_v x_v = rho
for B in range(p):
    for o in range(2,nopt):
        _,m1,t1,_=OPTS[o]; C,m2=partner[(B,m1)]
        oks=[o2 for o2 in range(2,nopt) if OPTS[o2][1]==m2 and (t1,OPTS[o2][2]) in allowed[(B,m1,C,m2)]]
        m_.AddBoolOr([z[C][o2] for o2 in oks]).OnlyEnforceIf(z[B][o])
s=cp_model.CpSolver(); s.parameters.max_time_in_seconds=TL; s.parameters.num_workers=W
found=0
with open(OUT,'a') as fo:
    while found<MAXSOL:
        st=s.Solve(m_)
        if st not in (cp_model.OPTIMAL,cp_model.FEASIBLE): print("status",s.StatusName(st),"after",found,flush=True); break
        sel=[next(o for o in range(nopt) if s.Value(z[B][o])) for B in range(p)]
        x=np.ones(N)
        for B,o in enumerate(sel):
            if o==1: x[blk(B)]=-1
            elif o>1: x[blk(B)]=1-2*IND[(OPTS[o][1],OPTS[o][2])]
        ver=np.allclose(S@x,rho*x)
        found+=1; row=" | ".join(OPTS[o][0] for o in sel)
        print(f"SOL {found} verified={ver}: {row}",flush=True); fo.write(f"K={K}\tverified={ver}\t{row}\n"); fo.flush()
        m_.AddBoolOr([z[B][sel[B]].Not() for B in range(p)])
print(f"done {found} {time.time()-t0:.0f}s")
