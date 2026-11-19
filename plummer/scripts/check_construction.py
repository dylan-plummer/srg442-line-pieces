#!/usr/bin/env python3
"""Standalone builder and verifier for the paper (plain Python, no packages, no code shared with the search).

Builds Mathon's conference graph M = srg(441,220,109,110) for q = 7, p = 9 from the formulas of Construction 2.2, adds an
isolated vertex oo, switches with respect to the set X of Construction 3.1, and checks by counting common neighbours
of all pairs that the result is srg(442,231,120,121) and its complement srg(442,210,99,100). With a graph6 file as
argument it also checks that the file is exactly the constructed (labelled) graph.

usage: python3 check_construction.py [../graphs/srg442_231_120_121.g6]
"""
import sys, itertools

Q = 7
QR = {(x * x) % Q for x in range(1, Q)}                  # {1, 2, 4}
def chi(x): x %= Q; return 0 if x == 0 else (1 if x in QR else -1)
G = 3                                                   # primitive root mod 7
BLOCKS = list(itertools.product(range(3), repeat=2))    # Paley(9) = K3 x K3, B = (B1, B2), index 3*B1 + B2
def pc_adj(B, C): return B != C and (B[0] == C[0] or B[1] == C[1])
def latin(i, j):                                        # circulant skew Latin square of order 9 on the indices
    d = (j - i) % 9; return d if d <= 4 else d - 9


def inner(k, a, b, a2, b2):
    """entry ((a,b),(a2,b2)) of A_k (k = 1..4; A_{-k} = A_k^T), blocks of order q^2 indexed by (a, b) in GF(7)^2"""
    if k < 0: return inner(-k, a2, b2, a, b)
    if k < 4: c = pow(G, 2 * k, Q) * (a + a2); return 1 if chi(b2 - c - b) == 1 else 0      # F P[g^{2k}(a+a')]
    return 1 if chi(b2 - a) == 1 else 0                                                        # E P[a]


def diag(a, b, a2, b2):
    """A_d: 0 if a = a', I + F if chi(a' - a) = 1, J - F otherwise"""
    if a == a2: return 0
    f = 1 if chi(b2 - b) == 1 else 0
    return (1 if b == b2 else 0) + f if chi(a2 - a) == 1 else 1 - f


VERT = [(B, a, b) for B in range(9) for a in range(Q) for b in range(Q)]   # vertex 49*B + 7*a + b; oo = 441


def mathon_adj(u, v):
    (B, a, b), (C, a2, b2) = VERT[u], VERT[v]
    if B == C: return diag(a, b, a2, b2)
    e = inner(latin(B, C), a, b, a2, b2)
    return e if pc_adj(BLOCKS[B], BLOCKS[C]) else 1 - e


def in_X(v):
    B, a, b = VERT[v]; l = (BLOCKS[B][0] + BLOCKS[B][1]) % 3
    return l == 2 or (l == 1 and chi(b - a) == 1) or (l == 0 and (a + b) % Q == 0)


def build():
    n = 442; M = [[0] * n for _ in range(n)]
    for u in range(441):
        for v in range(441):
            if u != v: M[u][v] = mathon_adj(u, v)
    X = [in_X(v) for v in range(441)] + [False]
    A = [[0] * n for _ in range(n)]
    for u in range(n):
        for v in range(n):
            if u != v: A[u][v] = M[u][v] ^ (X[u] != X[v])       # Seidel switching w.r.t. X
    return M, X, A


def srg_params(A):
    n = len(A); nb = [frozenset(j for j in range(n) if A[i][j]) for i in range(n)]
    assert all(A[i][j] == A[j][i] for i in range(n) for j in range(n)) and all(not A[i][i] for i in range(n))
    ks = {len(s) for s in nb}; lam, mu = set(), set()
    for i in range(n):
        for j in range(i + 1, n): (lam if j in nb[i] else mu).add(len(nb[i] & nb[j]))
    return n, ks, lam, mu


def graph6_decode(s):
    b = [c - 63 for c in s.strip().encode()]; n = b[0]; b = b[1:]
    if n == 63: n = (b[0] << 12) | (b[1] << 6) | b[2]; b = b[3:]
    bits = [(x >> (5 - i)) & 1 for x in b for i in range(6)]; A = [[0] * n for _ in range(n)]; t = 0
    for j in range(1, n):
        for i in range(j): A[i][j] = A[j][i] = bits[t]; t += 1
    return A


if __name__ == "__main__":
    M, X, A = build()
    print("Mathon graph:", srg_params([r[:441] for r in M[:441]]), "(expect 441, {220}, {109}, {110})")
    print("|X| =", sum(X))
    print("switched graph:", srg_params(A), "(expect 442, {231}, {120}, {121})")
    Ac = [[0 if i == j else 1 - A[i][j] for j in range(442)] for i in range(442)]
    print("complement:", srg_params(Ac), "(expect 442, {210}, {99}, {100})")
    if len(sys.argv) > 1:
        F = graph6_decode(open(sys.argv[1]).read().split()[0])
        print("identical to", sys.argv[1], ":", F == A)
