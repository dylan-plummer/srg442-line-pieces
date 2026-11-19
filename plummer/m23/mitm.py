"""Exact enumeration of local block patterns (h-invariant blocks of Mathon's two-graph, q = s^2 - 2).
Pattern = (e0, per slope state: alpha=+-1 or delta=+-1). Conditions: u-slopes: e0 + q al - A - P == 0 mod 2*(s+1)/2*... :
uses the ratios from struct.npz: u-slot u = r_u * w_partner (w even) -> u == 0 mod 2|r_u|; |P| <= floor(2q/(2|r_v|)) from
v_partner = r_v w;  v-slot v == 0 mod 2|r_v|, |P| <= floor(2q/(2|r_u|)).
python mitm.py q STRUCT.npz OUT.npy"""
import sys, time, numpy as np
q=int(sys.argv[1]); st=np.load(sys.argv[2],allow_pickle=True); OUT=sys.argv[3]
sig=st['sig'].astype(np.int64); ut=st['utype'][0]; ns=len(ut)
ru=int(abs(float(st['U'][0][2]))); rv=int(abs(float(st['V'][0][2])))
M=np.where(ut,2*ru,2*rv).astype(np.int64)          # modulus per slope
BND=np.where(ut,(2*q)//(2*rv),(2*q)//(2*ru))       # |P| bound per slope: |w_own| = |partner feature| / r <= 2q / r, w = 2P
sgn=np.where(ut,-1,1)                               # u = e0 + q al - A - P ; v = e0 + q al - A + P   (q al: S0 - (A - al) = e0+(q-1)al - A + al)
H1=list(range(ns//2)); H2=list(range(ns//2,ns))
print(f"q={q} slopes {ns} u-mod {2*ru} v-mod {2*rv} bounds u {BND[ut][0]} v {BND[~ut][0]}",flush=True)
def half(H):
    k=len(H); idx=np.arange(4**k,dtype=np.int64)
    S=((idx[:,None]//(4**np.arange(k)))%4).astype(np.int8)          # state per slope in H
    al=np.where(S==0,1,np.where(S==1,-1,0)).astype(np.int64); de=np.where(S==2,1,np.where(S==3,-1,0)).astype(np.int64)
    A=al.sum(1)
    P=de@sig[:,H].T                                                   # contributions to P_i for all i
    c=-A[:,None]+sgn[None,:]*P
    c[:,H]+=q*al
    return S,c,P
t0=time.time()
S1,c1,P1=half(H1); S2,c2,P2=half(H2)
print(f"halves built {time.time()-t0:.0f}s",flush=True)
base=np.cumprod(np.concatenate([[1],M[:-1]])).astype(np.int64)
k2=((((-c2)%M)*base).sum(1))
order=np.argsort(k2,kind='stable'); k2s=k2[order]
res=[]
for e0 in (1,-1):
    k1=((((e0+c1)%M)*base).sum(1))
    lo=np.searchsorted(k2s,k1,'left'); hi=np.searchsorted(k2s,k1,'right'); cnt=hi-lo
    tot=int(cnt.sum()); print(f" e0={e0}: residue matches {tot}",flush=True)
    i1=np.repeat(np.arange(len(k1)),cnt)
    off=np.arange(tot)-np.repeat(np.cumsum(cnt)-cnt,cnt)
    i2=order[np.repeat(lo,cnt)+off]
    CH=10**7; keep=[]
    for a in range(0,tot,CH):
        P=P1[i1[a:a+CH]]+P2[i2[a:a+CH]]
        ok=(np.abs(P)<=BND[None,:]).all(1)
        keep.append(np.stack([i1[a:a+CH][ok],i2[a:a+CH][ok]],1))
    kk=np.concatenate(keep) if keep else np.zeros((0,2),dtype=np.int64)
    print(f" e0={e0}: after |P| bounds {len(kk)}",flush=True)
    for x,y in kk: res.append(np.concatenate([[e0],S1[x],S2[y]]))
R=np.array(res,dtype=np.int8); np.save(OUT,R); print(f"total local patterns {len(R)} ({time.time()-t0:.0f}s)")
