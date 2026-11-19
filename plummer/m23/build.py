"""Mathon's conference graph for prime q = 4t-1, p = q+2 = s^2 (s prime), PC = Paley(p) on F_s[i]/(i^2 - n), labels s*B1 + B2.
Saves the Seidel matrix of M + {oo} (int8, oo last) and checks S^2 = (N-1) I in float32 BLAS."""
import sys, time, numpy as np
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # mathon.py, srgutil.py (copied from the srg442 repo)
from mathon import blocks
q=int(sys.argv[1]); s=int(sys.argv[2]); out=sys.argv[3]; p=s*s; assert p==q+2
nsq=next(x for x in range(2,s) if pow(x,(s-1)//2,s)==s-1)            # non-square mod s
def mul(x,y): return ((x[0]*y[1]+x[1]*y[0])%s, (x[1]*y[1]+nsq*x[0]*y[0])%s)   # (B1 i + B2)
F=[(a,b) for a in range(s) for b in range(s)]
sq={mul(x,x) for x in F if x!=(0,0)}
PC=np.array([[1 if (u!=v and ((u[0]-v[0])%s,(u[1]-v[1])%s) in sq) else 0 for v in F] for u in F],dtype=np.int64)
t,Ad,Ak=blocks(q); n=2*t; b=q*q; N=p*b+1
def L(i,j):
    d=(j-i)%p; return d if d<=n else d-p
def A(k): return Ak[k] if k>0 else Ak[-k].T
S=np.ones((N,N),dtype=np.int8); np.fill_diagonal(S,0)
for i in range(p):
    for j in range(p):
        blk=Ad if i==j else (A(L(i,j)) if PC[i,j] else 1-A(L(i,j)))
        S[i*b:(i+1)*b, j*b:(j+1)*b]=1-2*blk
np.fill_diagonal(S,0)
assert (S==S.T).all()
t0=time.time(); Sf=S.astype(np.float32); S2=Sf@Sf
ok=np.array_equal(S2, (N-1)*np.eye(N,dtype=np.float32))
print(f"q={q} p={p} N={N}: Seidel S^2 = (N-1) I: {ok}  ({time.time()-t0:.0f}s)",flush=True)
np.save(out,S)
