"""q = 23, step 1: per-arc line-piece solutions of the d = 1 template, for general q (fibre GF(q)^2, q prime = 3 mod 4).
On an arc a -> b (tail direction l_a positive, head direction l_b negative, block E(u,v) = [l_a(u) - l_b(v) in T]) a
switching vector x with x = g o l_a on fibre a and x = h o l_b on fibre b (g, h : F_q -> {+1,-1}, non-constant) gives
S x = rho x on the non-constant parts of the two fibres iff, with s~ = s - mean(s),
    head:  (rho + q) h~(w) = -2q sum_{t in T} g~(w + t)
    tail:  (rho - q) g~(w) = -2q sum_{t in T} h~(w - t).
Enumerates all 2^q functions g (exhaustive), derives h from the head equation, checks both equations exactly in
integers.  Output: json {rho: {T: [[g, h], ...]}} with g, h as +-1 lists.
usage: c3_arcs.py q OUTJSON"""
import sys, json, itertools, collections
import numpy as np

q = int(sys.argv[1]); outf = sys.argv[2]
p = q + 2; r = int(round(p ** 0.5)); assert r * r == p
Q = {(x * x) % q for x in range(1, q)}
Ts = {'Q': sorted(Q), '0N': sorted({0} | (set(range(q)) - Q)), 'N': sorted(set(range(1, q)) - Q), '0Q': sorted({0} | Q)}
res = {}
for rho in (q * r, -q * r):
    res[str(rho)] = {}
    for Tn, T in Ts.items():
        Mt = np.zeros((q, q), dtype=np.float32)        # V[:, w] = sum_{t in T} g[(w + t) % q]
        for w in range(q):
            for t in T: Mt[(w + t) % q, w] += 1
        sols = []
        CH = 1 << 20
        for c0 in range(0, 1 << q, CH):
            idx = np.arange(c0, min(c0 + CH, 1 << q), dtype=np.int64)
            G = (1 - 2 * ((idx[:, None] >> np.arange(q)[None, :]) & 1)).astype(np.float32)
            V = G @ Mt
            vmin = V.min(1, keepdims=True); vmax = V.max(1, keepdims=True)
            two = ((V == vmin) | (V == vmax)).all(1) & (vmax[:, 0] > vmin[:, 0])
            for k in np.nonzero(two)[0]:
                g = G[k].astype(np.int64); v = V[k]
                # h~ = -2q u / (rho + q), u = v - |T| Sg / q: h~ is a decreasing (rho+q>0) / increasing function of v
                h = np.where(v == v.max(), 1, -1) if rho + q < 0 else np.where(v == v.min(), 1, -1)
                h = h.astype(np.int64)
                if abs(h.sum()) == q: continue
                Sg, Sh = int(g.sum()), int(h.sum())
                gq = q * g - Sg; hq = q * h - Sh                   # q * (s - mean)
                head = all((rho + q) * hq[w] == -2 * q * sum(gq[(w + t) % q] for t in T) for w in range(q))
                tail = all((rho - q) * gq[w] == -2 * q * sum(hq[(w - t) % q] for t in T) for w in range(q))
                if head and tail: sols.append([g.tolist(), h.tolist()])
        cnt = collections.Counter((sum(g), sum(h)) for g, h in sols)
        print(f'q={q} rho={rho} T={Tn}: {len(sols)} per-arc solutions; (M_tail, M_head) classes {dict(cnt)}', flush=True)
        res[str(rho)][Tn] = sols
json.dump(res, open(outf, 'w'))
