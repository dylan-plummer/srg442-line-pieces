"""Block structure for h-invariant blocks: per slope m two half-line orbits (lambda*w_m, lambda in Q / N) + origin.
Extract sigma(m,m') (which half of slope m' the line l_m = tq meets), verify closed forms, and slope types (u/v) and signs
per (B,m) from pairs.npy. Saves struct.npz."""
import sys, numpy as np
q=int(sys.argv[1]); PR=np.load(sys.argv[2],allow_pickle=True); OUT=sys.argv[3]
Qs=sorted({(x*x)%q for x in range(1,q)}); QS=set(Qs); tq=Qs[0]; tn=sorted(set(range(1,q))-QS)[0]
slopes=[str(i) for i in range(q)]+['inf']; ns=len(slopes)
def lab(ms,a,b): return (a if ms=='inf' else (b-int(ms)*a))%q
dirv={ms:((0,1) if ms=='inf' else (1,int(ms))) for ms in slopes}   # direction vectors: slope m line through 0 is {(a, m a)}; inf: {(0,b)}
# careful: line through 0 of slope m = {l_m = 0} = {b = m a} -> direction (1,m); for 'inf' l = a = 0 -> direction (0,1)
def half(pt):   # which (slope, half) a nonzero point lies on
    a,b=pt
    for ms in slopes:
        if lab(ms,a,b)==0:
            d=dirv[ms]; lam=(a*pow(d[0],q-2,q))%q if d[0] else (b*pow(d[1],q-2,q))%q
            return ms,('Q' if lam in QS else 'N')
sig=np.zeros((ns,ns),dtype=np.int64)
for i,ms in enumerate(slopes):
    pts=[(a,b) for a in range(q) for b in range(q) if lab(ms,a,b)==tq]
    for pt in pts:
        m2,hh=half(pt); j=slopes.index(m2); assert j!=i
        sig[i,j]=1 if hh=='Q' else -1
    # verify the N-intercept line meets the opposite halves
    for pt in [(a,b) for a in range(q) for b in range(q) if lab(ms,a,b)==tn]:
        m2,hh=half(pt); j=slopes.index(m2); assert sig[i,j]==(1 if hh=='N' else -1)
print("sigma: symmetric?",np.array_equal(sig,sig.T),"antisym?",np.array_equal(sig,-sig.T)," sigma sigma^T diag/offdiag:",sorted(set(np.diag(sig@sig.T))),sorted(set((sig@sig.T)[~np.eye(ns,dtype=bool)].tolist()))[:6])
# random verification of closed forms
rng=np.random.default_rng(0)
for _ in range(200):
    e0=rng.choice([-1,1]); al=np.zeros(ns,int); de=np.zeros(ns,int)
    for i in range(ns):
        if rng.random()<0.5: al[i]=rng.choice([-1,1])
        else: de[i]=rng.choice([-1,1])
    x={}
    for a in range(q):
        for b in range(q):
            if (a,b)==(0,0): x[(a,b)]=e0; continue
            ms,hh=half((a,b)); i=slopes.index(ms); x[(a,b)]=al[i]+(de[i] if hh=='Q' else -de[i])
    A=al.sum(); P=sig@de
    for i,ms in enumerate(slopes):
        S0=sum(x[(a,b)] for a in range(q) for b in range(q) if lab(ms,a,b)==0)
        SQ=sum(x[(a,b)] for a in range(q) for b in range(q) if lab(ms,a,b)==tq)
        SN=sum(x[(a,b)] for a in range(q) for b in range(q) if lab(ms,a,b)==tn)
        assert S0==e0+(q-1)*al[i] and SQ==A-al[i]+P[i] and SN==A-al[i]-P[i], (S0,SQ,SN)
print("closed forms verified: S0 = e0 + (q-1) al_m, SQ = A - al_m + P_m, SN = A - al_m - P_m, P = sigma delta")
# slope types per block from pair equations
types={}
for r in PR:
    B,ms=int(r[0]),str(r[1]); cq=np.rint(np.array(r[4:11],float)).astype(int); cn=np.rint(np.array(r[11:18],float)).astype(int)
    C,ms2=int(r[2]),str(r[3])
    # u-type if Q-eq own part is (-k,k,0) and partner part is (0, c, -c): u_B = (c/k) w_C
    if cq[2]==0 and cq[3]==0 and cq[0]==-cq[1] and cq[4]==-cq[5]:
        types[(B,ms)]=('u',cq[4]//cq[1]*(-1)*-1, C,ms2)   # u_B = ratio * w_C with ratio = -(c/k)*(-1)...
    elif cn[1]==0 and cn[3]==0 and cn[0]==-cn[2] and cn[4]==-cn[5]:
        types[(B,ms)]=('v',None,C,ms2)
    else: raise SystemExit(f"unrecognized pattern {cq} {cn}")
# exact ratios: Q-eq: -k S0 + k SQ + c SQ_C - c SN_C = 0 -> k(SQ - S0) = -c w_C -> u_B = (c/k) w_C
U={};V={}
for r in PR:
    B,ms=int(r[0]),str(r[1]); cq=np.rint(np.array(r[4:11],float)).astype(int); cn=np.rint(np.array(r[11:18],float)).astype(int)
    if types[(B,ms)][0]=='u': U[(B,ms)]=(cq[4]/cq[1], int(r[2]), str(r[3]))
    else: V[(B,ms)]=(cn[4]/cn[2], int(r[2]), str(r[3]))     # -k S0 + k SN + c SQ_C - c SN_C = 0 -> k(SN - S0) = -c w_C -> v_B = (c/k) w_C ; here k = cn[2]
cnt=[(sum(1 for (B,ms) in U if B==b0), sum(1 for (B,ms) in V if B==b0)) for b0 in range(max(k[0] for k in types)+1)]
print("u/v slots per block:",sorted(set(cnt)),"; u ratios",sorted(set(v[0] for v in U.values())),"v ratios",sorted(set(v[0] for v in V.values())))
utype=np.array([[ (B,ms) in U for ms in slopes] for B in range(len(cnt))])
print("u-type slope sets identical across blocks?", all((utype[0]==utype[B]).all() for B in range(len(cnt))), "; u-type slopes:", [ms for i,ms in enumerate(slopes) if utype[0,i]])
np.savez(OUT,sig=sig,utype=utype,slopes=np.array(slopes),U=np.array([(k[0],k[1],*v) for k,v in U.items()],dtype=object),V=np.array([(k[0],k[1],*v) for k,v in V.items()],dtype=object),allow_pickle=True)
