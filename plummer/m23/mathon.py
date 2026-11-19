"""Mathon's PC-graphs (conference graphs) of order p q^2, q = 4t-1 prime, p = 4t+1 the order of a PC-graph (Mathon,
Canad. J. Math. 30 (1978), Thm 4.1), as read from the paper; every variant is checked with is_srg.

GF(q) (q prime here), a_1..a_q the field elements, P[c] the translation x -> x + c, F_xy = [chi(y - x) = 1] (Paley skew
design), E_xy = [chi(y) = 1] (constant rows). Blocks of order q^2 indexed by (i, j) in GF(q)^2:
  A_d  : 0 if i = j, I + F if chi(j - i) = 1, J - F otherwise          (the Kronecker PC-graph of order q^2, (1.10))
  A_k  : F P[g^(2k) (a_i + a_j)], 1 <= k <= 2t-1;   A_2t : E P[a_i]
Big matrix, p x p blocks: A_d on the diagonal, A(l_ij) if the order-p PC-graph has the edge ij, J - A(l_ij) otherwise,
with L the circulant skew Latin square l_ij = j - i reduced to [-2t, 2t] and A(-k) = A(k)^T.

python mathon.py Q OUT.g6        (Q = 3 -> srg(45,22,10,11), Q = 7 -> srg(441,220,109,110) with Paley(9) = K3 x K3)
"""
import sys, os, itertools
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from srgutil import graph6, is_srg


def blocks(q):
    t = (q + 1) // 4; sq = {(x * x) % q for x in range(1, q)}; chi = lambda a: 0 if a % q == 0 else (1 if a % q in sq else -1)
    g = next(x for x in range(2, q) if len({pow(x, e, q) for e in range(1, q)}) == q - 1) if q > 2 else 1
    I = np.eye(q, dtype=np.int64); J = np.ones((q, q), dtype=np.int64)
    P = lambda c: np.array([[1 if (y - x) % q == c % q else 0 for y in range(q)] for x in range(q)], dtype=np.int64)
    F = np.array([[1 if chi(y - x) == 1 else 0 for y in range(q)] for x in range(q)], dtype=np.int64)
    E = np.array([[1 if chi(y) == 1 else 0 for y in range(q)] for x in range(q)], dtype=np.int64)
    def big(f): return np.block([[f(i, j) for j in range(q)] for i in range(q)])
    Ad = big(lambda i, j: 0 * J if i == j else (I + F if chi(j - i) == 1 else J - F))
    Ak = {k: big(lambda i, j, k=k: F @ P(pow(g, 2 * k, q) * (i + j))) for k in range(1, 2 * t)}
    Ak[2 * t] = big(lambda i, j: E @ P(i))
    return t, Ad, Ak


def mathon(q, PC):
    t, Ad, Ak = blocks(q); p = len(PC); n = 2 * t; JJ = np.ones_like(Ad)
    def L(i, j):
        d = (j - i) % p; return d if d <= n else d - p
    def A(k): return Ak[k] if k > 0 else Ak[-k].T
    return np.block([[Ad if i == j else (A(L(i, j)) if PC[i, j] else JJ - A(L(i, j))) for j in range(p)] for i in range(p)])


if __name__ == "__main__":
    q, out = int(sys.argv[1]), sys.argv[2]
    if q == 3:                                                  # pentagon
        PC = np.array([[1 if (j - i) % 5 in (1, 4) else 0 for j in range(5)] for i in range(5)], dtype=np.int64)
    else:                                                       # Paley(9) = K3 x K3, vertices ordered so the circulant
        V = list(itertools.product(range(3), repeat=2))         # labelling is the one used for L (any order works for Thm 4.1)
        PC = np.array([[1 if (a != b and (a[0] == b[0] or a[1] == b[1])) else 0 for b in V] for a in V], dtype=np.int64)
    A = mathon(q, PC); print("order", len(A), "symmetric", bool((A == A.T).all()), "srg", is_srg(A))
    if is_srg(A):
        with open(out, "w") as fh: fh.write(graph6(A) + "\n")
