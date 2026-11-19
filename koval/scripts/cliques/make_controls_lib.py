"""Controls: block graphs of known Steiner systems (positive), written as graph6.
Expected Delsarte-clique counts are stated with the reason; the counter must reproduce them."""
import itertools, random, sys, os


def g6(n, edges):
    adj = [[0]*n for _ in range(n)]
    for a, b in edges: adj[a][b] = adj[b][a] = 1
    bits = [adj[i][j] for j in range(1, n) for i in range(j)]
    while len(bits) % 6: bits.append(0)
    if n < 63: head = chr(n+63)
    else: head = '~' + chr(((n >> 12) & 63)+63) + chr(((n >> 6) & 63)+63) + chr((n & 63)+63)
    body = ''.join(chr(63+int(''.join(map(str, bits[i:i+6])), 2)) for i in range(0, len(bits), 6))
    return head + body

def block_graph(blocks):
    blocks = [frozenset(b) for b in blocks]
    E = [(i, j) for i, j in itertools.combinations(range(len(blocks)), 2) if blocks[i] & blocks[j]]
    return g6(len(blocks), E)

def check_design(pts, blocks, lam=1):
    c = {}
    for b in blocks:
        for p in itertools.combinations(sorted(b), 2): c[p] = c.get(p, 0)+1
    assert len(c) == len(pts)*(len(pts)-1)//2 and set(c.values()) == {lam}, "not a 2-design"

def proj_lines(q_dim, q):   # lines of PG(d,q), q prime
    vecs = [v for v in itertools.product(range(q), repeat=q_dim+1) if any(v)]
    def norm(v):
        i = next(k for k in range(len(v)) if v[k]); inv = pow(v[i], q-2, q)
        return tuple(x*inv % q for x in v)
    pts = sorted(set(map(norm, vecs))); idx = {p: i for i, p in enumerate(pts)}
    lines = set()
    for a, b in itertools.combinations(pts, 2):
        L = frozenset(idx[norm(tuple((s*x+t*y) % q for x, y in zip(a, b)))] for s in range(q) for t in range(q) if s or t)
        lines.add(L)
    return list(range(len(pts))), [sorted(L) for L in lines]

def aff_lines(dim, q):     # lines of AG(dim,q), q prime
    pts = list(itertools.product(range(q), repeat=dim)); idx = {p: i for i, p in enumerate(pts)}
    lines = set()
    for a, b in itertools.combinations(pts, 2):
        d = tuple((y-x) % q for x, y in zip(a, b))
        lines.add(frozenset(idx[tuple((x+t*dd) % q for x, dd in zip(a, d))] for t in range(q)))
    return list(range(len(pts))), [sorted(L) for L in lines]

