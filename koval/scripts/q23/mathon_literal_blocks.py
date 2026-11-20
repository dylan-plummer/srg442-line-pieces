"""Mathon's literal blocks at q = 3, 7, 23 are d = 1 template data (Lemma 2.4 of the paper, and the reduction of
the q = 23 per-arc search on Mathon's literal build to the sets Q and {0} u N).

Uses the literal transcription of Mathon (3.1)-(3.6), (4.1)-(4.3) in scripts/q23/w1lib.py
(mathon_blocks, skew_latin_circulant; the code that reproduces Mathon's printed (4.15)/(4.16) at q = 3).
Mathon's matrix of Theorem 4.1 consists only of the blocks A_d (diagonal), A_k, A_k^T and J - A_k, J - A_k^T
(1 <= k <= 2t, t = (q+1)/4).  The script checks, exactly:
  (F) A_d is the Cayley graph on F_q^2 whose connection set is the union of (q+1)/2 lines through 0, minus 0
      (L = the set of these lines);
  (X) every A_k is a line incidence: A_k[u, v] = [f_k(u) - f'_k(v) in T_k] for non-zero linear forms f_k, f'_k, with
      ker f_k not in L (tail, positive), ker f'_k in L (head, negative) and T_k an affine image alpha*Q + beta of the
      non-zero squares Q (then J - A_k is the incidence with T = F_q minus T_k = alpha({0} u N) + beta);
  (P) the tail kernels of A_1..A_{2t} are the (q+1)/2 lines outside L, each once, and the head kernels are the lines
      of L, each once;
  (L) every row of the circulant skew-Latin square contains each of +-1, ..., +-2t exactly once, so every base point
      is the tail of exactly one A_k and the head of exactly one A_k for every k (the regular tournament of (T2)
      and the port bijections of (T3)).
Then Mathon's literal two-graph satisfies (T1)-(T3), and the per-arc equations (D2) of the paper for its sets T reduce to
those for Q or {0} u N by the affine substitution g -> g(alpha . + beta) (bijection of solutions).
Controls: q = 3 and q = 7 (where Mathon's two-graph is known to be a template instance, K833 Sec. 6.5),
and a negative control: A_d itself is NOT a line incidence (the test must fail on it).
"""
import sys
sys.path.insert(0, 'scripts/q23')
import numpy as np
import w1lib as W


def setup(q):
    pts = [(x, y) for x in range(q) for y in range(q)]          # index q*x + y (outer x, inner y), as in w1lib
    lines = [frozenset(((x * 1) % q, (c * x) % q) for x in range(q)) for c in range(q)] + \
            [frozenset((0, y) for y in range(q))]                # slope c lines and the vertical line
    forms = [(lambda c: (lambda x, y: (y - c * x) % q))(c) for c in range(q)] + [lambda x, y: x % q]
    fv = [np.array([f(x, y) for (x, y) in pts]) for f in forms]  # form j has kernel lines[j]
    for j in range(q + 1):
        assert {pts[i] for i in range(q * q) if fv[j][i] == 0} == set(lines[j])
    return pts, lines, fv


def row_kernel(B, pts, lines, q):
    """index j of the line L such that the rows of B are constant exactly on the cosets of L, or None."""
    _, cls = np.unique(B, axis=0, return_inverse=True)
    cls = np.asarray(cls).ravel()
    if len(set(cls.tolist())) != q: return None
    zero = {pts[i] for i in range(q * q) if cls[i] == cls[0]}
    for j, L in enumerate(lines):
        if zero == set(L): return j
    return None


def incidence(B, q, pts, lines, fv):
    jt = row_kernel(B, pts, lines, q); jh = row_kernel(B.T, pts, lines, q)
    if jt is None or jh is None: return None
    a, b = fv[jt], fv[jh]
    for lam in range(1, q):
        T = {int(a[u]) for u in range(q * q) if B[u, int(np.argmax(b == 0))]}   # column with f'(v) = 0
        # B[u, v] == [a(u) - lam * b(v) in T]
        D = (a[:, None] - lam * b[None, :]) % q
        if (np.isin(D, list(T)).astype(np.uint8) == B).all():
            return jt, jh, lam, sorted(T)
    return None


def affine_of_Q(T, q):
    Q = {(x * x) % q for x in range(1, q)}
    for al in range(1, q):
        for be in range(q):
            if {(al * s + be) % q for s in Q} == set(T): return al, be
    return None


def check(q):
    mb = W.mathon_blocks(q); t = mb['t']; pts, lines, fv = setup(q)
    Ad = mb['Ad']; n = q * q
    # (F) Cayley with a union of (q+1)/2 lines
    nb0 = {pts[v] for v in range(n) if Ad[0, v]}
    for u in range(n):
        x, y = pts[u]
        nbu = {((pts[v][0] - x) % q, (pts[v][1] - y) % q) for v in range(n) if Ad[u, v]}
        assert nbu == nb0, 'fiber not Cayley'
    Lset = [j for j, L in enumerate(lines) if set(L) - {(0, 0)} <= nb0]
    assert len(Lset) == (q + 1) // 2 and set().union(*[set(lines[j]) for j in Lset]) - {(0, 0)} == nb0, 'fiber lines'
    # negative control: the fiber block is not a line incidence
    assert incidence(Ad, q, pts, lines, fv) is None, 'negative control failed'
    tails, heads, out = [], [], []
    for k in range(1, 2 * t + 1):
        r = incidence(mb['Ak'][k], q, pts, lines, fv)
        assert r is not None, 'A_%d not a line incidence' % k
        jt, jh, lam, T = r
        assert jt not in Lset and jh in Lset, 'signs at A_%d' % k
        af = affine_of_Q(T, q); assert af is not None, 'T of A_%d not an affine image of Q' % k
        # J - A_k: the complementary set is an affine image of {0} u N
        Tc = sorted(set(range(q)) - set(T))
        rc = incidence(1 - mb['Ak'][k], q, pts, lines, fv); assert rc is not None and rc[:3] == (jt, jh, lam) and rc[3] == Tc
        tails.append(jt); heads.append(jh); out.append((k, jt, jh, lam, T, af))
    assert sorted(tails) == sorted(set(range(q + 1)) - set(Lset)) and sorted(heads) == sorted(Lset), 'ports'
    p = q + 2; Ls = W.skew_latin_circulant(p)
    for i in range(p):
        row = sorted(int(Ls[i, j]) for j in range(p) if j != i)
        assert row == sorted(list(range(1, 2 * t + 1)) + [-k for k in range(1, 2 * t + 1)]), 'skew-Latin row %d' % i
        assert all(Ls[i, j] == -Ls[j, i] for j in range(p) if j != i)
    slopes = sorted(('inf' if j == q else j) for j in Lset)
    print('q=%d: fiber = Cayley graph on %d lines (slopes %s); A_1..A_%d are line incidences with T an affine image of Q;'
          ' tails = the %d positive lines, heads = the %d negative lines, each once; skew-Latin rows OK' %
          (q, len(Lset), slopes, 2 * t, len(tails), len(heads)), flush=True)
    for (k, jt, jh, lam, T, af) in out:
        print('   A_%d: tail kernel %s, head kernel %s, lambda %d, T = %s = %d*Q + %d' %
              (k, 'inf' if jt == q else jt, 'inf' if jh == q else jh, lam, T, af[0], af[1]), flush=True)


if __name__ == '__main__':
    for q in (3, 7, 23):
        check(q)
    print('ALL CHECKS PASSED')
