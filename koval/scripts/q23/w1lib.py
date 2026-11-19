"""Mathon's conference matrices of order p q^2 + 1, and the d = 1 Paley-fibre template.
Letters as in Mathon 1978: q = 4t - 1 (fibre field GF(q), prime here), p = 4t + 1 = q + 2
(base PC-graph / conference graph).  Vertices of the PC-graph: (base i, fibre index), fibre index = outer*q + inner.

build_mathon(q, base): literal transcription of Mathon (3.1)-(3.6), (4.1)-(4.3) (Canad. J. Math. 30 (1978) 325-328).
build_template(q, base): our d = 1 template (srg833 Prop. 6.6): Pal(q^2) fibres, circulant tournament i -> i + d,
    out-arc of difference d uses the d-th positive direction at the tail, the d-th negative at the head,
    T = Q_q on base edges, {0} u N_q on base non-edges.
seidel_square_check(A): exact check S^2 = (n-1) I for S = Seidel matrix of A + isolated point (float32 BLAS, every
    partial sum is an integer of absolute value <= n < 2^24, hence exact).
"""
import numpy as np

# ---------------------------------------------------------------- fields
def qr_set(q):
    return {(x * x) % q for x in range(1, q)}

def least_nonres(q):
    Q = qr_set(q)
    return min(x for x in range(2, q) if x not in Q)

def primitive_root(q):
    for g in range(2, q):
        if len({pow(g, k, q) for k in range(q - 1)}) == q - 1:
            return g

class GF:
    """GF(r) (r prime) or GF(r^2) = GF(r)[w], w^2 = least non-residue; element (x, y) <-> index r*x + y."""
    def __init__(self, order):
        self.order = order
        r = int(round(order ** 0.5))
        if r * r == order and all(r % d for d in range(2, r)):
            self.r, self.deg = r, 2
        else:
            assert all(order % d for d in range(2, order)), 'prime or prime^2 only'
            self.r, self.deg = order, 1
        r = self.r
        if self.deg == 1:
            self.els = [(x,) for x in range(r)]
            self.sq = {(x * x % r,) for x in range(1, r)}
        else:
            n0 = least_nonres(r); self.n0 = n0
            self.els = [(x, y) for x in range(r) for y in range(r)]
            mul = lambda a, b: ((a[0] * b[0] + n0 * a[1] * b[1]) % r, (a[0] * b[1] + a[1] * b[0]) % r)
            self.sq = {mul(e, e) for e in self.els if e != (0, 0)}
        assert len(self.sq) == (order - 1) // 2

    def sub(self, a, b):
        return tuple((x - y) % self.r for x, y in zip(a, b))

    def paley_adj(self):
        n = self.order
        A = np.zeros((n, n), dtype=np.uint8)
        for i, a in enumerate(self.els):
            for j, b in enumerate(self.els):
                if i != j and self.sub(a, b) in self.sq:
                    A[i, j] = 1
        return A

def check_conference_graph(Ab):
    p = len(Ab); k = (p - 1) // 2
    A = Ab.astype(np.int64); A2 = A @ A
    off = ~np.eye(p, dtype=bool)
    return bool((A.sum(1) == k).all() and (A == A.T).all() and (A2[(A == 1)] == (p - 5) // 4).all()
                and (A2[(A == 0) & off] == (p - 1) // 4).all())

# ---------------------------------------------------------------- Mathon's literal blocks
def mathon_blocks(q):
    assert q % 4 == 3 and all(q % d for d in range(2, q))
    t = (q + 1) // 4
    a = list(range(1, q)) + [0]                       # a_i = i, a_q = 0: a_{q-i} = -a_i (Theorem 1.2 numbering)
    idx = {v: i for i, v in enumerate(a)}
    Q = qr_set(q); chi1 = np.array([[1 if (x % q) in Q else 0 for x in range(q)]])[0]
    av = np.array(a)
    g = primitive_root(q)
    D = (av[None, :] - av[:, None]) % q               # a_j - a_i
    F = chi1[D].astype(np.uint8)                      # (3.3)
    E = np.tile(chi1[(av - a[0]) % q], (q, 1)).astype(np.uint8)   # (3.2): e_ij = [chi(a_j - a_1) = 1]

    def FP(c):  # (F P[c])_{r,s} = [chi(a_s - c - a_r) = 1]
        return chi1[(av[None, :] - c - av[:, None]) % q].astype(np.uint8)

    def EP(c):
        return np.tile(chi1[(av - c - a[0]) % q], (q, 1)).astype(np.uint8)
    I = np.eye(q, dtype=np.uint8); J = np.ones((q, q), dtype=np.uint8)
    n = q * q
    Ad = np.zeros((n, n), dtype=np.uint8)
    for i in range(q):
        for j in range(q):
            if i == j: continue
            Ad[i*q:(i+1)*q, j*q:(j+1)*q] = (I + F) if F[i, j] == 1 else (J - F)
    Ak = {}
    for k in range(1, 2 * t):
        M = np.zeros((n, n), dtype=np.uint8); c0 = pow(g, 2 * k, q)
        for i in range(q):
            for j in range(q):
                M[i*q:(i+1)*q, j*q:(j+1)*q] = FP(c0 * (a[i] + a[j]) % q)
        Ak[k] = M
    M = np.zeros((n, n), dtype=np.uint8)
    for i in range(q):
        for j in range(q):
            M[i*q:(i+1)*q, j*q:(j+1)*q] = EP(a[i])
    Ak[2 * t] = M
    return {'t': t, 'g': g, 'a': a, 'E': E, 'F': F, 'Ad': Ad, 'Ak': Ak}

def skew_latin_circulant(p):
    n = (p - 1) // 2
    L = np.zeros((p, p), dtype=np.int64)
    for i in range(1, p + 1):
        for j in range(1, p + 1):
            if i - j > n: L[i-1, j-1] = j - i + p
            elif j - i > n: L[i-1, j-1] = j - i - p
            else: L[i-1, j-1] = j - i
    return L

def build_mathon(q, base_adj, L=None):
    """Theorem 4.1: A_ij = A_d (i = j), A(l_ij) (base edge), J - A(l_ij) (non-edge); A(k) = A_k, A(-k) = A_k^T."""
    p = len(base_adj); assert p == q + 2
    mb = mathon_blocks(q)
    if L is None: L = skew_latin_circulant(p)
    n = q * q; N = p * n
    A = np.zeros((N, N), dtype=np.uint8)
    for i in range(p):
        for j in range(p):
            if i == j: B = mb['Ad']
            else:
                l = int(L[i, j]); B = mb['Ak'][l] if l > 0 else mb['Ak'][-l].T
                if not base_adj[i, j]: B = 1 - B
            A[i*n:(i+1)*n, j*n:(j+1)*n] = B
    return A, mb

# ---------------------------------------------------------------- d = 1 Paley-fibre template
def template_setup(q):
    Fq2 = GF(q * q)
    pts = Fq2.els                                       # (x, y), index q*x + y
    Af = Fq2.paley_adj()
    dirs = [(1, 0), (0, 1)] + [(1, b) for b in range(1, q)]
    X = np.array([u[0] for u in pts]); Y = np.array([u[1] for u in pts])
    formv = {l: (l[0] * X + l[1] * Y) % q for l in dirs}
    Sf = 1.0 - np.eye(q * q) - 2.0 * Af
    sign = {}
    z = np.exp(2j * np.pi / q)
    for l in dirs:
        v = z ** formv[l]; w = Sf @ v; r = w / v
        assert np.allclose(r, r[0]) and abs(abs(r[0]) - q) < 1e-6
        sign[l] = int(round(r[0].real / q))
    pos = [l for l in dirs if sign[l] == 1]; neg = [l for l in dirs if sign[l] == -1]
    assert len(pos) == len(neg) == (q + 1) // 2
    return {'Af': Af, 'dirs': dirs, 'formv': formv, 'pos': pos, 'neg': neg}

def build_template(q, base_adj, setup=None):
    p = len(base_adj); n = q * q; N = p * n
    st = setup or template_setup(q)
    Q = qr_set(q); Tq = np.array([1 if x in Q else 0 for x in range(q)], dtype=np.uint8)
    T0N = np.array([1 if (x == 0 or x not in Q) else 0 for x in range(q)], dtype=np.uint8)
    A = np.zeros((N, N), dtype=np.uint8)
    for i in range(p):
        A[i*n:(i+1)*n, i*n:(i+1)*n] = st['Af']
    arcs = {}
    for i in range(p):
        for d in range(1, (p - 1) // 2 + 1):
            j = (i + d) % p
            la, lb = st['pos'][d-1], st['neg'][d-1]
            adj = bool(base_adj[i, j]); T = Tq if adj else T0N
            E = T[(st['formv'][la][:, None] - st['formv'][lb][None, :]) % q]
            A[i*n:(i+1)*n, j*n:(j+1)*n] = E
            A[j*n:(j+1)*n, i*n:(i+1)*n] = E.T
            arcs[(i, j)] = (la, lb, 'Q' if adj else '0N')
    return A, st, arcs

# ---------------------------------------------------------------- verification
def seidel_f32(A):
    N = len(A); n = N + 1
    S = np.ones((n, n), dtype=np.float32)
    S[:N, :N] -= 2.0 * A.astype(np.float32)
    np.fill_diagonal(S, 0.0)
    return S

def seidel_square_check(A, B=1024, log=None):
    """Exact: all n^2 entries of S^2 compared with (n-1) I.  Returns dict."""
    S = seidel_f32(A); n = len(S)
    assert (A == A.T).all() and not A.diagonal().any()
    bad = 0; maxdev = 0.0
    for r0 in range(0, n, B):
        R = S[r0:r0+B] @ S
        R[np.arange(len(R)), np.arange(r0, r0 + len(R))] -= (n - 1)
        m = float(np.abs(R).max()); maxdev = max(maxdev, m)
        bad += int(np.count_nonzero(R))
        if log: log(f'  rows {r0}..{r0+len(R)-1}: max|S^2-(n-1)I| = {m}')
    return {'n': n, 'entries_checked': n * n, 'nonzero_deviations': bad, 'max_dev': maxdev, 'ok': bad == 0}

def srg_check(A, B=1024):
    """Exact srg check via float32 A^2 (entries <= n < 2^24)."""
    N = len(A); Af = A.astype(np.float32); k = A.sum(1)
    if not (k == k[0]).all(): return {'regular': False}
    lam = set(); mu = set()
    for r0 in range(0, N, B):
        R = Af[r0:r0+B] @ Af
        a = A[r0:r0+B].astype(bool)
        dia = np.zeros_like(a); dia[np.arange(len(a)), np.arange(r0, r0 + len(a))] = True
        if not np.array_equal(R[dia], k[r0:r0+len(a)].astype(np.float32)): return {'regular': True, 'diag': False}
        lam |= set(np.unique(R[a]).astype(int).tolist()); mu |= set(np.unique(R[~a & ~dia]).astype(int).tolist())
    return {'regular': True, 'n': N, 'k': int(k[0]), 'lambda': sorted(lam), 'mu': sorted(mu),
            'srg': len(lam) == 1 and len(mu) == 1}

def to_g6(A):
    n = len(A)
    bits = np.concatenate([A[j, :j] for j in range(1, n)]).astype(np.uint8)  # (j, i), i < j, by j then i == graph6 order
    pad = (-len(bits)) % 6
    bits = np.concatenate([bits, np.zeros(pad, dtype=np.uint8)]).reshape(-1, 6)
    vals = (bits * np.array([32, 16, 8, 4, 2, 1], dtype=np.uint8)).sum(1).astype(np.uint8) + 63
    if n < 63: head = chr(63 + n)
    else: head = '~' + chr(63 + (n >> 12)) + chr(63 + ((n >> 6) & 63)) + chr(63 + (n & 63))
    return head + vals.tobytes().decode('ascii')

def from_g6(s):
    s = s.strip()
    if s[0] == '~': n = ((ord(s[1]) - 63) << 12) | ((ord(s[2]) - 63) << 6) | (ord(s[3]) - 63); body = s[4:]
    else: n = ord(s[0]) - 63; body = s[1:]
    b = np.frombuffer(body.encode(), dtype=np.uint8) - 63
    bits = np.unpackbits(b[:, None], axis=1)[:, 2:].reshape(-1)
    A = np.zeros((n, n), dtype=np.uint8); pos = 0
    for j in range(1, n):
        A[j, :j] = bits[pos:pos + j]; pos += j
    assert not bits[pos:].any()
    A |= A.T
    return A
