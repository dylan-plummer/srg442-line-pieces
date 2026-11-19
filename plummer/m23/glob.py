"""Global exact search over local block patterns (h-invariant regular members of Mathon's two-graph).
python glob.py S.npy q K STRUCT.npz LOC.npy OUT [TGRP=0|1] [ROOTS=a:b slice of root candidates] """
import sys, os, time, numpy as np, collections
S=np.load(sys.argv[1]); q=int(sys.argv[2]); K=int(sys.argv[3]); st=np.load(sys.argv[4],allow_pickle=True); LOC=np.load(sys.argv[5]).astype(np.int64); OUT=sys.argv[6]
N=len(S); b=q*q; p=(N-1)//b; s=int(round(p**0.5)); rho=N-1-2*K
sig=st['sig'].astype(np.int64); ut=st['utype'][0]; slopes=[str(x) for x in st['slopes']]; ns=len(slopes); SI={m:i for i,m in enumerate(slopes)}
TG=os.environ.get("TGRP","0")=="1"
# pattern features
e0=LOC[:,0]; Sst=LOC[:,1:]
al=np.where(Sst==0,1,np.where(Sst==1,-1,0)); de=np.where(Sst==2,1,np.where(Sst==3,-1,0))
A=al.sum(1); P=de@sig.T
sg=np.where(ut,-1,1); F=e0[:,None]+q*al-A[:,None]+sg[None,:]*P      # u (u-slopes) or v (v-slopes)
T=e0+(q-1)*A                                                       # block total
Np=len(LOC)
# constraints: (B, i) -> (C, j, r): F[B,i] = r * 2 * P[C,j]
cons=collections.defaultdict(list)
for arr in (st['U'],st['V']):
    for B,ms,r,C,ms2 in arr: cons[int(B)].append((SI[str(ms)],int(C),SI[str(ms2)],int(round(float(r)))))
# index per slot: (F value, P value) -> pattern ids
idx=[collections.defaultdict(list) for _ in range(ns)]
for i in range(ns):
    for k,(f,pp) in enumerate(zip(F[:,i],P[:,i])): idx[i][(int(f),int(pp))].append(k)
idx=[{k:np.array(v,dtype=np.int64) for k,v in d.items()} for d in idx]
# t-orbits on blocks
orbit=list(range(p))
if TG:
    for B in range(p):
        B1,B2=divmod(B,s); orbit[B]=B2       # t: B1 -> B1+1, so orbits = fixed B2
reps=sorted(set(orbit)); members={o:[B for B in range(p) if orbit[B]==o] for o in reps}
# Paley mean condition: -q sum_{C!=B} sigP_BC T_C + q^2 = rho T_B ; sigP from S block sums
sigP=np.zeros((p,p),dtype=np.int64)
for B in range(p):
    for C in range(p):
        if B!=C: sigP[B,C]=int(round(S[B*b:(B+1)*b,C*b:(C+1)*b].astype(np.float64).sum()/(-q**3)))
def full_check(assign):
    Tv=np.array([T[assign[orbit[B]]] for B in range(p)])
    if not all(-q*int((sigP[B]*Tv).sum())+b==rho*int(Tv[B]) for B in range(p)): return False,"mean"
    return True,""
# unary filters (constraints inside one orbit)
dom={o:np.arange(Np) for o in reps}
for o in reps:
    for B in members[o]:
        for (i,C,j,r) in cons[B]:
            if orbit[C]==o:
                d=dom[o]; ok=F[d,i]==r*2*P[d,j]; dom[o]=d[ok]
print(f"q={q} K={K} TGRP={int(TG)}: {Np} local patterns, orbits {len(reps)}, domains after unary {[len(dom[o]) for o in reps]}",flush=True)
sols=[]; nodes=[0]; t0=time.time()
def dfs(assign,dom):
    nodes[0]+=1
    free=[o for o in reps if o not in assign]
    if not free:
        ok,why=full_check(assign)
        if ok: sols.append(dict(assign)); print(f"SOLUTION {len(sols)} at {time.time()-t0:.0f}s nodes {nodes[0]}",flush=True)
        return
    o=min(free,key=lambda x:len(dom[x]))
    for k in dom[o]:
        newdom={}; dead=False
        for o2 in free:
            if o2==o: continue
            d=dom[o2]
            for B in members[o]:
                for (i,C,j,r) in cons[B]:
                    if orbit[C]!=o2: continue
                    # F[B,i] = r*2*P[C,j]  and the partner's own constraint back: F[C,j] = r'*2*P[B,i]
                    need_P=F[k,i]
                    if need_P % (2*r): dead=True; break
                    needP=need_P//(2*r)
                    d=d[P[d,j]==needP]
                    for (i2,C2,j2,r2) in cons[C]:
                        if C2==B and j2==i: d=d[F[d,j]==r2*2*P[k,i]]
                    if len(d)==0: dead=True; break
                if dead: break
            if dead: break
            newdom[o2]=d
        if dead: continue
        if not assign and os.environ.get("TSYM")=="1" and not TG:
            for j in range(1,s):
                Bj=j*s
                if Bj in newdom: newdom[Bj]=newdom[Bj][newdom[Bj]>=k]
                if Bj in newdom and len(newdom[Bj])==0: dead=True
            if dead: continue
        assign[o]=int(k); nd=dict(dom); nd.update(newdom); nd[o]=np.array([k])
        dfs(assign,nd); del assign[o]
        if len(sols)>=int(os.environ.get("MAXSOL","1000000")): return
roots=os.environ.get("ROOTS")
if roots:
    a,bb=map(int,roots.split(':')); o0=min(reps,key=lambda x:len(dom[x])); dom[o0]=dom[o0][a:bb]
dfs({},dom)
print(f"done: {len(sols)} solutions, {nodes[0]} nodes, {time.time()-t0:.0f}s",flush=True)
with open(OUT,'a') as f:
    for sl in sols: f.write(f"q={q}\tK={K}\tTGRP={int(TG)}\t"+" ".join(str(sl[orbit[B]]) for B in range(p))+"\n")
