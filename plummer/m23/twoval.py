# All T subset F_q with D_T(y) = sum_{x in T} chi'(y - x) (chi'(0) = -1) taking exactly two values; report gaps.
import sys, numpy as np, collections
q=int(sys.argv[1]); Qs={(x*x)%q for x in range(1,q)}
chi=lambda u: -1 if u%q==0 else (1 if u%q in Qs else -1)
M=np.array([[chi(y-x) for x in range(q)] for y in range(q)],dtype=np.float32)
res=collections.Counter(); ex={}
CH=1<<20
for start in range(0,1<<q,CH):
    idx=np.arange(start,min(start+CH,1<<q),dtype=np.int64)
    U=((idx[None,:]>>np.arange(q)[:,None])&1).astype(np.float32)
    D=M@U; mx=D.max(0); mn=D.min(0)
    nvals=np.array([len(np.unique(D[:,j])) for j in np.nonzero(mx-mn>0)[0]]) if False else None
    # two-valued: every entry is mx or mn
    two=((D==mx[None,:])|(D==mn[None,:])).all(0) & (mx>mn)
    for j in np.nonzero(two)[0]:
        g=int(mx[j]-mn[j]); T=int(idx[j]); size=bin(T).count('1')
        res[(g,size)]+=1
        ex.setdefault((g,size),T)
for (g,size),c in sorted(res.items()): print(f"gap {g} |T|={size}: {c} sets, e.g. {[i for i in range(q) if ex[(g,size)]>>i&1]}")
