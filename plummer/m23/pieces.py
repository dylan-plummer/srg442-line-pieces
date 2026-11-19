"""Line-piece regular members of Mathon's two-graph (general prime q = 4t-1, p = q+2 blocks of q^2 points (a,b)).
Per block: empty, full, or {(a,b): f(a,b) in T}, f in {b - m a (m in F_q), a}, T a nonempty proper union of {0}, Q, N.
Solve (S - rho I) x = 0 (oo unswitched) with CP-SAT, enumerate with blocking clauses.
python pieces.py S.npy q K TLIM WORKERS MAXSOL OUT"""
import sys, time, numpy as np
from ortools.sat.python import cp_model
S=np.load(sys.argv[1]); q=int(sys.argv[2]); K=int(sys.argv[3]); TL=float(sys.argv[4]); W=int(sys.argv[5]); MAXSOL=int(sys.argv[6]); OUT=sys.argv[7]
N=len(S); b=q*q; p=(N-1)//b; rho=N-1-2*K
Qs={(x*x)%q for x in range(1,q)}; Ns=set(range(1,q))-Qs
TS={"0":{0},"Q":Qs,"N":Ns,"0Q":{0}|Qs,"0N":{0}|Ns,"*":Qs|Ns}
forms={f"{m}":(lambda a,bb,m=m:(bb-m*a)%q) for m in range(q)}; forms["inf"]=lambda a,bb:a%q
OPTS=[("-",np.zeros(b)),("all",np.ones(b))]
for fn,f in forms.items():
    for tn,T in TS.items():
        OPTS.append((f"{fn}:{tn}",np.array([1. if f(a,bb) in T else 0. for a in range(q) for bb in range(q)])))
V=np.stack([v for _,v in OPTS],1).astype(np.float32); nopt=V.shape[1]
t0=time.time()
Sf=S.astype(np.float32); c0=(Sf.sum(1)-rho)            # R 1, R = S - rho I
coef=np.zeros((N,p,nopt),dtype=np.float32)
for B in range(p):
    blk=Sf[:,B*b:(B+1)*b]@V
    blk[B*b:(B+1)*b]-=rho*V
    coef[:,B,:]=blk
coef=np.rint(coef).astype(np.int64); c0=np.rint(c0).astype(np.int64)
assert (c0%2==0).all()
keys={}
for i in range(N):
    keys.setdefault(coef[i].tobytes()+int(c0[i]).to_bytes(8,'little',signed=True), i)
print(f"q={q} p={p} N={N} K={K} rho={rho}: {nopt} options/block, {len(keys)} distinct equations ({time.time()-t0:.0f}s)",flush=True)
m=cp_model.CpModel(); z=[[m.NewBoolVar("") for o in range(nopt)] for B in range(p)]
for B in range(p): m.AddExactlyOne(z[B])
for i in keys.values():
    m.Add(sum(int(coef[i,B,o])*z[B][o] for B in range(p) for o in range(nopt) if coef[i,B,o])==int(c0[i]//2))
s=cp_model.CpSolver(); s.parameters.max_time_in_seconds=TL; s.parameters.num_workers=W
found=0
with open(OUT,'a') as fo:
    while found<MAXSOL:
        st=s.Solve(m)
        if st not in (cp_model.OPTIMAL,cp_model.FEASIBLE):
            print("status",s.StatusName(st),"after",found,"solutions",flush=True); break
        sel=[next(o for o in range(nopt) if s.Value(z[B][o])) for B in range(p)]
        found+=1; row=" | ".join(OPTS[o][0] for o in sel)
        print(f"SOL {found}: {row}",flush=True); fo.write(f"K={K}\t{row}\n"); fo.flush()
        m.AddBoolOr([z[B][sel[B]].Not() for B in range(p)])
print(f"done {found} solutions {time.time()-t0:.0f}s")
