"""Rebuild x from a glob*.txt line (pattern ids per block) and check S x = rho x exactly; write graph6 if regular.
python verify_hit.py S.npy q K STRUCT.npz LOC.npy GLOB.txt"""
import sys, numpy as np
S=np.load(sys.argv[1]); q=int(sys.argv[2]); K=int(sys.argv[3]); st=np.load(sys.argv[4],allow_pickle=True); LOC=np.load(sys.argv[5]).astype(int)
N=len(S); b=q*q; p=(N-1)//b; rho=N-1-2*K; slopes=[str(x) for x in st['slopes']]
QS={(x*x)%q for x in range(1,q)}
dirv={ms:((0,1) if ms=='inf' else (1,int(ms))) for ms in slopes}
def lab(ms,a,bb): return (a if ms=='inf' else (bb-int(ms)*a))%q
def half(a,bb):
    for i,ms in enumerate(slopes):
        if lab(ms,a,bb)==0:
            d=dirv[ms]; lam=(a*pow(d[0],q-2,q))%q if d[0] else (bb*pow(d[1],q-2,q))%q
            return i,(lam in QS)
H=[[None if (a,bb)==(0,0) else half(a,bb) for bb in range(q)] for a in range(q)]
for line in open(sys.argv[6]):
    ids=list(map(int,line.split('\t')[-1].split())); x=np.ones(N)
    for B,k in enumerate(ids):
        pat=LOC[k]; e0=pat[0]; stt=pat[1:]
        for a in range(q):
            for bb in range(q):
                if (a,bb)==(0,0): val=e0
                else:
                    i,isq=H[a][bb]; s_=stt[i]
                    val={0:1,1:-1}[s_] if s_<2 else ((1 if s_==2 else -1)*(1 if isq else -1))
                x[B*b+a*q+bb]=val
    ok=np.array_equal(S.astype(np.float32)@x.astype(np.float32),rho*x.astype(np.float32))
    print("eigenvector:",ok,"|X| =",int((x<0).sum()))
