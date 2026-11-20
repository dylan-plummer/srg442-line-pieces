#!/usr/bin/env python3
"""Checks for the sporadicity section, on the Seidel matrices of Mathon's two-graphs (oo last, vertex q^2 i + q a + b).
python3 -I sporadic_check.py S7.npy [S23.npy]

1. Lemma (slope decomposition), for every fiber i and slope m, on random mean-zero g:
   S_ii (g o l_m) = -q g o l_m for m in {0} u Q, +q for m in N u {oo};
   for exactly one j != i, S_ji (g o l_m) != 0, and then S_ji (g o l_m) = q Sigma_ji (kappa * g) o l_{m*} with m* = -m
   (0* = oo), kappa = chi'(-.) if m* in N u {oo}, chi' if m* in {0} u Q; Sigma_ji = -(row sum of S_ji)/q = +-1.
2. Proposition (pair equations), on the 14 members of srg442's Table (rho = -21 for A, +21 for B), every partner pair:
   u_{j,m*} = Sigma_ij (eps s + 1) w_{i,m} / 2  and  v_{i,m} = -Sigma_ij (eps s - 1) w_{j,m*} / 2,
   (i, m) the slot with m in {0} u Q, (j, m*) its partner, rho = eps q s; u, v, w from the line sums L(0), L(Q), L(N).
   Also the closed forms of the line sums of an h-invariant fiber function (sigma = bordered Jacobsthal matrix)."""
import sys, numpy as np

def setup(q):
    Q = {(x * x) % q for x in range(1, q)}
    chi = np.array([0] + [1 if y in Q else -1 for y in range(1, q)])
    chip = chi.copy(); chip[0] = -1
    slopes = list(range(q)) + ["inf"]
    lab = {m: np.array([(a if m == "inf" else (b - m * a)) % q for a in range(q) for b in range(q)]) for m in slopes}
    return Q, chi, chip, slopes, lab

def conv(k, g, q):          # (k * g)(y) = sum_x k(y - x) g(x)
    return np.array([sum(k[(y - x) % q] * g[x] for x in range(q)) for y in range(q)])

def star(m, q): return "inf" if m == 0 else (0 if m == "inf" else (-m) % q)

def check_lemma(S, q, rng):
    Q, chi, chip, slopes, lab = setup(q)
    b = q * q; p = (len(S) - 1) // b
    minus = set(Q) | {0}
    chipr = np.array([chip[(-y) % q] for y in range(q)])
    for i in range(p):
        for j in range(p):
            if i != j:
                rs = S[i * b:(i + 1) * b, j * b:(j + 1) * b].sum(1)
                assert np.all(rs == rs[0]) and abs(rs[0]) == q
            else:
                assert np.all(S[i * b:(i + 1) * b, i * b:(i + 1) * b].sum(1) == 0)
    for i in range(p):
        Sii = S[i * b:(i + 1) * b, i * b:(i + 1) * b]
        for m in slopes:
            g = rng.standard_normal(q); g -= g.mean(); f = g[lab[m]]
            lam = -q if m in minus else q
            assert np.allclose(Sii @ f, lam * f)
            hits = []
            for j in range(p):
                if j == i: continue
                Sji = S[j * b:(j + 1) * b, i * b:(i + 1) * b]
                r = Sji @ f
                if np.abs(r).max() > 1e-8: hits.append((j, r, -Sji.sum(1)[0] / q))
            assert len(hits) == 1
            j, r, sig = hits[0]; ms = star(m, q)
            kap = chip if ms in minus else chipr
            assert np.allclose(r, q * sig * conv(kap, g, q)[lab[ms]]), (q, i, m, j)
    return True

TABLE = {   # srg442 Table 2: columns 00 01 02 10 11 12 20 21 22 (fiber i = 3 B1 + B2)
 "A1": "6:0 1:Q all 1:Q all 6:0 all 6:0 1:Q", "A2": "1:Q 6:0 all all 1:Q 6:0 6:0 all 1:Q",
 "A3": "0:0N - all 5:* 2:0N inf:0 3:0 all 4:Q", "A4": "4:0N all - 5:0 2:Q all 3:* 3:0 4:0N",
 "A5": "3:0 all 4:Q 5:* 2:0N 5:0 2:0N - all",
 "B1": "- - - inf:0Q inf:0Q 5:0Q 2:* 0:* 0:*", "B2": "5:N 2:* - 3:0Q all 4:0 5:N 2:* -",
 "B3": "5:0Q 2:0 all 3:N - 4:* 3:N - 4:*", "B4": "- 4:* 3:0Q - 4:* 3:0Q - 4:* 3:0Q",
 "B5": "2:* - 5:0Q 2:* - 5:0Q 2:* - 5:0Q", "B6": "- inf:0Q 0:* - inf:0Q 0:* - inf:0Q 0:*",
 "B7": "4:* 3:0Q - 4:* 3:0Q - 4:* 3:0Q -", "B8": "0:* - inf:0Q 0:* - inf:0Q 0:* - inf:0Q",
 "B9": "inf:0Q 0:* - inf:0Q 0:* - inf:0Q 0:* -"}

def member(row, q=7):
    Q, chi, chip, slopes, lab = setup(q); N = set(range(1, q)) - Q
    T = {"0": {0}, "Q": Q, "N": N, "0Q": {0} | Q, "0N": {0} | N, "*": Q | N}
    x = np.ones(9 * q * q + 1, dtype=np.int64)
    for i, e in enumerate(row.split()):
        blk = slice(i * q * q, (i + 1) * q * q)
        if e == "all": x[blk] = -1
        elif e != "-":
            m, t = e.split(":"); m = "inf" if m == "inf" else int(m)
            x[blk] = np.where(np.isin(lab[m], list(T[t])), -1, 1)
    return x

def check_pairs(S, q=7, s=3):
    Q, chi, chip, slopes, lab = setup(q); b = q * q; minus = set(Q) | {0}
    # sigma and closed forms
    d = {m: ((0, 1) if m == "inf" else (1, m)) for m in slopes}
    ell = lambda m, v: v[0] % q if m == "inf" else (v[1] - m * v[0]) % q
    sig = np.array([[0 if m == mm else chi[ell(m, d[mm])] for mm in slopes] for m in slopes])
    assert np.array_equal(sig @ sig.T, q * np.eye(q + 1, dtype=int))
    npairs = 0
    for name, row in TABLE.items():
        x = member(row); rho = -21 if name[0] == "A" else 21; eps = rho // (q * s)
        assert np.array_equal(S @ x, rho * x), name
        Lsum = lambda i, m: tuple(int(x[i * b:(i + 1) * b][lab[m] == c].sum()) for c in (0, 1, 3))
        for i in range(9):
            xf = x[i * b:(i + 1) * b]
            # closed forms (every line piece is h-invariant)
            e0 = int(xf[0]); al = {}; de = {}
            for m in slopes:
                hp = {int(xf[q * (t * d[m][0] % q) + t * d[m][1] % q]) for t in Q}
                hn = {int(xf[q * (t * d[m][0] % q) + t * d[m][1] % q]) for t in set(range(1, q)) - Q}
                assert len(hp) == len(hn) == 1
                fp, fn = hp.pop(), hn.pop(); al[m] = (fp + fn) // 2; de[m] = (fp - fn) // 2
            A = sum(al.values())
            for k, m in enumerate(slopes):
                P = sum(sig[k, kk] * de[mm] for kk, mm in enumerate(slopes))
                assert Lsum(i, m) == (e0 + (q - 1) * al[m], A - al[m] + P, A - al[m] - P)
            for m in minus:
                g = np.zeros(q); g[1] = 1; g -= g.mean(); f = g[lab[m]]
                js = [j for j in range(9) if j != i and np.abs(S[j * b:(j + 1) * b, i * b:(i + 1) * b] @ f).max() > 1e-8]
                assert len(js) == 1; j = js[0]; ms = star(m, q)
                sgn = -S[i * b:(i + 1) * b, j * b:(j + 1) * b].sum(1)[0] // q
                Li, Lj = Lsum(i, m), Lsum(j, ms)
                ui, vi, wi = Li[0] - Li[1], Li[0] - Li[2], Li[1] - Li[2]
                uj, vj, wj = Lj[0] - Lj[1], Lj[0] - Lj[2], Lj[1] - Lj[2]
                assert 2 * uj == sgn * (eps * s + 1) * wi, (name, i, m, j)
                assert 2 * vi == -sgn * (eps * s - 1) * wj, (name, i, m, j)
                npairs += 1
    return npairs

if __name__ == "__main__":
    rng = np.random.default_rng(1)
    S7 = np.load(sys.argv[1]).astype(np.float64)
    print("lemma q=7:", check_lemma(S7, 7, rng))
    print("pair equations on the 14 members: pairs checked", check_pairs(S7.astype(np.int64)))
    if len(sys.argv) > 2:
        S23 = np.load(sys.argv[2]).astype(np.float32)
        print("lemma q=23:", check_lemma(S23, 23, rng))
    print("ok")
