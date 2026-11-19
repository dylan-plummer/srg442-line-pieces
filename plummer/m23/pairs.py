"""Express the exact eigen-equations of Mathon's two-graph for h-invariant x in per-slope features.
Block B, slope m: F_{B,m} = (S0, SQ, SN) = sums of x over the line through 0 and over one line with intercept in Q / in N.
For each (B,m) find partner (C,m*) and fit E_{B,m}(Q) - E_{B,m}(0), E_{B,m}(N) - E_{B,m}(0) as integer combinations of
F_{B,m}, F_{C,m*}; verify exactly on random +-1 h-invariant vectors. Also the block-mean equation.
python pairs.py S.npy q K OUT.npz"""
import sys, numpy as np, itertools
S=np.load(sys.argv[1]).astype(np.float64); q=int(sys.argv[2]); K=int(sys.argv[3]); OUT=sys.argv[4]
N=len(S); b=q*q; p=(N-1)//b; rho=N-1-2*K
Qs=sorted({(x*x)%q for x in range(1,q)}); Ns=sorted(set(range(1,q))-set(Qs)); tq,tn=Qs[0],Ns[0]
slopes=list(range(q))+['inf']
def lab(m): return np.array([(a if m=='inf' else (bb-m*a))%q for a in range(q) for bb in range(q)])
LAB={m:lab(m) for m in slopes}
def feats(xB,m):   # xB on F_q^2 (length b)
    L=LAB[m]; return np.array([xB[L==0].sum(), xB[L==tq].sum(), xB[L==tn].sum()])
blk=lambda B: slice(B*b,(B+1)*b)
# partner map via slope test functions
def slope_of(g):
    for m in slopes:
        L=LAB[m]
        if all(np.ptp(g[L==c])<1e-6 for c in range(q)): return m
partner={}
for C in range(p):
    for m2 in slopes:
        f=(LAB[m2]==tq).astype(float); f-=f.mean()
        for B in range(p):
            if B==C: continue
            g=S[blk(B),blk(C)]@f; g-=g.mean()
            if np.abs(g).max()>1e-6: partner[(B,slope_of(g))]=(C,m2)
assert len(partner)==p*(q+1)
# h-invariant random test vectors
g0=next(x for x in range(2,q) if len({pow(x,e,q) for e in range(1,q)})==q-1); hs=g0*g0%q
orbit_id=-np.ones(b,dtype=int); k=0
for v in range(b):
    if orbit_id[v]>=0: continue
    a,bb=divmod(v,q); u=(a,bb)
    for e in range(q):
        w=((pow(hs,e,q)*a)%q)*q+(pow(hs,e,q)*bb)%q; orbit_id[w]=k
    k+=1
norb=k; rng=np.random.default_rng(0)
def E_of(x):      # line-sum equations for all blocks/slopes: ((S - rho I) x)(v) + S[v,oo]
    y=S[:N-1,:N-1]@x[:N-1]-rho*x[:N-1]+S[:N-1,N-1]
    return y
T=40; samples=[]
for _ in range(T):
    x=np.ones(N); 
    for B in range(p): x[blk(B)]=rng.choice([-1.,1.],norb)[orbit_id]
    samples.append((x,E_of(x)))
coef={}
bad=0
for (B,m),(C,m2) in partner.items():
    Xf=[]; Yq=[]; Yn=[]
    for x,y in samples:
        Xf.append(np.concatenate([feats(x[blk(B)],m),feats(x[blk(C)],m2),[1]]))
        fy=feats(y[blk(B)],m); Yq.append(fy[1]-fy[0]); Yn.append(fy[2]-fy[0])
    Xf=np.array(Xf); cq,_,_,_=np.linalg.lstsq(Xf,np.array(Yq),rcond=None); cn,_,_,_=np.linalg.lstsq(Xf,np.array(Yn),rcond=None)
    if np.abs(Xf@cq-Yq).max()>1e-6 or np.abs(Xf@cn-Yn).max()>1e-6: bad+=1
    coef[(B,str(m))]=(C,str(m2),np.round(cq,6),np.round(cn,6))
print(f"q={q} K={K}: pairs {len(partner)}, non-exact fits {bad}")
seen=set()
for key in list(coef)[:6]:
    C,m2,cq,cn=coef[key]; print(key,'->',(C,m2),'\n   E_Q-E_0:',cq,'\n   E_N-E_0:',cn)
distinct=sorted({(tuple(v[2]),tuple(v[3])) for v in coef.values()})
print("distinct coefficient patterns:",len(distinct))
for d in distinct: print("  ",d)
np.save(OUT,np.array([(k[0],k[1],v[0],v[1],*v[2],*v[3]) for k,v in coef.items()],dtype=object),allow_pickle=True)
