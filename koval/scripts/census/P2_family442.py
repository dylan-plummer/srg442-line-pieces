"""P2 : the d = 1 template on 441 = 7^2 * 9 vertices, generalised.

Vertex set Z_9 x F_7^2 (fibre i, point u = (x, y), index 49 i + 7 x + y).

A CHOICE consists of
  * L[i]   : a set of 4 of the 8 lines through 0 of AG(2,7) for every fibre i; fibre i is the Cayley graph on F_7^2
             with connection set (union of L[i]) minus 0 -- an srg(49,24,11,12) (Latin-square type LS_4(7));
             Paley(49) and Mathon's Kronecker PC-graph are both of this kind;
  * arcs   : one entry (a, b, fa, fb, T) per pair of fibres, meaning E_ab(u, v) = [fa(u) - fb(v) in T]
             (fa, fb non-zero linear forms on F_7^2, T a subset of F_7; phase t and slopes are inside fa, fb, T).

Template conditions (sufficient, Prop. d1 of srg833.tex generalised; checked by `template_ok`):
  (T1) the quotient: the pairs {a,b} with |T| = 3 form a conference graph on 9 points (Paley(9) = K3 x K3);
  (T2) |T| = 3: T is a (7,3,1) difference set (c*Q + t); |T| = 4: T is the complement of one;
  (T3) at every fibre a each of the 8 lines is the kernel of exactly one form used at a, and for every arc the tail
       form has its kernel OUTSIDE L[tail] (s = +1) and the head form INSIDE L[head] (s = -1)  [=> regular tournament].
Every choice is verified exactly (A^2 = kI + lam A + mu (J-I-A)); the template conditions are only used to generate.

Two-graph invariant: the canonical certificate (pynauty) of the double cover on 2 * 442 vertices of the Seidel matrix
of G + isolated point: (x, s) ~ (y, t) iff x != y and S_xy * s * t = -1.  Isomorphic double covers <=> switching-
equivalent Seidel matrices, i.e. isomorphic two-graphs (antipodal pairs are the only pairs with no common neighbour).

Run: python3 -I scripts/census/P2_family442.py <mode> ...
"""
import sys, os, hashlib, itertools, random, time, json
import numpy as np
import pynauty

P = 7
QR = frozenset({1, 2, 4}); NR = frozenset({3, 5, 6})
INV = {x: pow(x, P - 2, P) for x in range(1, P)}
PTS = [(x, y) for x in range(P) for y in range(P)]
NV = 9 * 49
K, LAM, MU = 220, 109, 110
# lines through 0: index c in 0..6 = span (1, c); index 7 = span (0, 1)
DIRV = [(1, c) for c in range(P)] + [(0, 1)]


def line_of_vec(v):
    x, y = v[0] % P, v[1] % P
    if x: return (y * INV[x]) % P
    assert y; return 7


def kernel_line(f):
    a, b = f
    return line_of_vec((b, -a))


def form_val(f):
    """vector of f(u) over the 49 points"""
    return np.array([(f[0] * x + f[1] * y) % P for (x, y) in PTS], dtype=np.int64)


def fibre_adj(L):
    conn = set()
    for k in L:
        dx, dy = DIRV[k]
        for t in range(1, P): conn.add(((t * dx) % P, (t * dy) % P))
    M = np.zeros((49, 49), dtype=np.uint8)
    for i, (x, y) in enumerate(PTS):
        for j, (z, w) in enumerate(PTS):
            if ((x - z) % P, (y - w) % P) in conn: M[i, j] = 1
    return M


PALEY9 = np.array([[1 if a != b and (a // 3 == b // 3 or a % 3 == b % 3) else 0 for b in range(9)] for a in range(9)],
                  dtype=np.int64)       # K3 x K3 = Paley(9) with F_9 = F_3[w], w^2 = 2, index 3x + y


def build(choice):
    A = np.zeros((NV, NV), dtype=np.uint8)
    for i in range(9):
        A[49 * i:49 * i + 49, 49 * i:49 * i + 49] = fibre_adj(choice['L'][i])
    for (a, b, fa, fb, T) in choice['arcs']:
        va, vb = form_val(fa), form_val(fb)
        D = (va[:, None] - vb[None, :]) % P
        E = np.isin(D, list(T)).astype(np.uint8)
        A[49 * a:49 * a + 49, 49 * b:49 * b + 49] = E
        A[49 * b:49 * b + 49, 49 * a:49 * a + 49] = E.T
    return A


def verify(A):
    A = A.astype(np.int64)
    if not (A == A.T).all() or A.trace(): return False, 'not symmetric / loops'
    deg = A.sum(1)
    if not (deg == K).all(): return False, 'degrees %s' % sorted(set(deg.tolist()))
    A2 = A @ A
    J = np.ones_like(A); I = np.eye(NV, dtype=np.int64)
    ok = (A2 == K * I + LAM * A + MU * (J - I - A)).all()
    return bool(ok), 'srg(441,220,109,110)' if ok else 'A^2 mismatch'


def is_diffset(T):
    T = list(T)
    d = sorted((s - t) % P for s in T for t in T if s != t)
    return d == list(range(1, P))


def template_ok(choice, delta=PALEY9):
    used = {i: [] for i in range(9)}; seen = set()
    for (a, b, fa, fb, T) in choice['arcs']:
        key = frozenset((a, b)); assert key not in seen; seen.add(key)
        T = frozenset(T)
        if len(T) == 3:
            if not is_diffset(T) or not delta[a, b]: return False, 'T3 at %d,%d' % (a, b)
        elif len(T) == 4:
            if not is_diffset(frozenset(range(P)) - T) or delta[a, b]: return False, 'T4 at %d,%d' % (a, b)
        else: return False, '|T|'
        ka, kb = kernel_line(fa), kernel_line(fb)
        if ka in choice['L'][a] or kb not in choice['L'][b]: return False, 'sign at arc %d->%d' % (a, b)
        used[a].append(ka); used[b].append(kb)
    if len(seen) != 36: return False, 'pairs'
    for i in range(9):
        if sorted(used[i]) != list(range(8)): return False, 'lines at %d' % i
    return True, 'ok'


# ---------------- graph6 + certificates ----------------
def g6(A):
    n = len(A); bits = [int(A[i, j]) for j in range(n) for i in range(j)]
    while len(bits) % 6: bits.append(0)
    s = '~' + chr(63 + (n >> 12)) + chr(63 + ((n >> 6) & 63)) + chr(63 + (n & 63)) if n >= 63 else chr(63 + n)
    return s + ''.join(chr(63 + int(''.join(map(str, bits[t:t + 6])), 2)) for t in range(0, len(bits), 6))


def load_g6(line):
    s = line.strip()
    if s[0] == '~': m = ((ord(s[1]) - 63) << 12) | ((ord(s[2]) - 63) << 6) | (ord(s[3]) - 63); body = s[4:]
    else: m = ord(s[0]) - 63; body = s[1:]
    bits = []
    for ch in body:
        v = ord(ch) - 63; bits.extend((v >> (5 - i)) & 1 for i in range(6))
    A = np.zeros((m, m), dtype=np.uint8); idx = 0
    for j in range(m):
        for i in range(j):
            if bits[idx]: A[i, j] = A[j, i] = 1
            idx += 1
    return A


def seidel_plus_inf(A):
    n = len(A); S = np.ones((n + 1, n + 1), dtype=np.int64) - np.eye(n + 1, dtype=np.int64)
    S[:n, :n] -= 2 * A.astype(np.int64)
    return S      # infinity isolated (S = +1 to every vertex)


def cover_cert(S):
    """canonical certificate of the double cover of the Seidel matrix S (exact two-graph invariant)."""
    N = len(S); adj = {}
    for x in range(N):
        neg = np.nonzero(S[x] == -1)[0]; pos = np.nonzero(S[x] == 1)[0]
        # (x,+) ~ (y,+),(y,-) according to S_xy*s*t = -1
        adj[2 * x] = [2 * int(y) for y in neg] + [2 * int(y) + 1 for y in pos]
        adj[2 * x + 1] = [2 * int(y) + 1 for y in neg] + [2 * int(y) for y in pos]
    g = pynauty.Graph(2 * N, adjacency_dict=adj)
    return hashlib.sha256(pynauty.certificate(g)).hexdigest()


def cover_aut(S):
    N = len(S); adj = {}
    for x in range(N):
        neg = np.nonzero(S[x] == -1)[0]; pos = np.nonzero(S[x] == 1)[0]
        adj[2 * x] = [2 * int(y) for y in neg] + [2 * int(y) + 1 for y in pos]
        adj[2 * x + 1] = [2 * int(y) + 1 for y in neg] + [2 * int(y) for y in pos]
    gens, mm, e, orb, no = pynauty.autgrp(pynauty.Graph(2 * N, adjacency_dict=adj))
    return int(round(mm * 10 ** e)), no


def graph_aut(A):
    n = len(A); adj = {i: [int(j) for j in np.nonzero(A[i])[0]] for i in range(n)}
    gens, mm, e, orb, no = pynauty.autgrp(pynauty.Graph(n, adjacency_dict=adj))
    return int(round(mm * 10 ** e)), no


# ---------------- the two controls ----------------
def mathon_choice():
    """Mathon's conference graph of order 441 (Canad. J. Math. 30 (1978) Thm 4.1, as read by Plummer): fibre
    (i, x) ~ (j, y) iff j - i != 0 and y - x in {0} u chi(j-i) Q  -> lines (1, c), c in {0,1,2,4}; arc P -> Q iff
    (Q - P) mod 9 in {1..4} = k; k <= 3: E = [y - x - m(i+j) in Q], m = 3^(2k) mod 7; k = 4: E = [y - i in Q];
    complement for non-adjacent P, Q in Paley(9) = K3 x K3."""
    L = [frozenset({0, 1, 2, 4})] * 9; arcs = []
    for Pp in range(9):
        for k in range(1, 5):
            Qq = (Pp + k) % 9
            if k <= 3: m = pow(3, 2 * k, P); fa, fb = (m, 1), ((-m) % P, 1)
            else: fa, fb = (1, 0), (0, 1)
            T = NR if PALEY9[Pp, Qq] else frozenset({0}) | QR
            arcs.append((Pp, Qq, fa, fb, T))
    return {'L': L, 'arcs': arcs, 'name': 'mathon'}


def paley49_lines():
    # F_49 = F_7[w], w^2 = 3 (lib833: first non-residue mod 7 is 3); u = x + y w
    sq = set()
    for x in range(P):
        for y in range(P):
            if (x, y) != (0, 0): sq.add(((x * x + 3 * y * y) % P, (2 * x * y) % P))
    L = frozenset(k for k in range(8) if ((DIRV[k][0]) % P, DIRV[k][1] % P) in sq)
    assert len(L) == 4; return L


def ours_choice():
    """our d1_z9 graph (mathon_compare2.py): Paley(49) fibres, pos/neg forms in lib833.lines() order,
    i -> i+d uses pos[d-1] at i and neg[d-1] at i+d; T = QR on Paley(9)-edges, {0} u NR otherwise."""
    L = paley49_lines(); forms = [(1, 0), (0, 1)] + [(1, b) for b in range(1, P)]
    pos = [f for f in forms if kernel_line(f) not in L]; neg = [f for f in forms if kernel_line(f) in L]
    arcs = []
    for i in range(9):
        for d in (1, 2, 3, 4):
            j = (i + d) % 9
            arcs.append((i, j, pos[d - 1], neg[d - 1], QR if PALEY9[i, j] else frozenset({0}) | NR))
    return {'L': [L] * 9, 'arcs': arcs, 'name': 'ours_d1_z9'}


def sha(s): return hashlib.sha256(s.encode()).hexdigest()


def graph_cert(A):
    n = len(A); adj = {i: [int(j) for j in np.nonzero(A[i])[0]] for i in range(n)}
    return hashlib.sha256(pynauty.certificate(pynauty.Graph(n, adjacency_dict=adj))).hexdigest()


# ---------------- line-piece switching sets ----------------
from fractions import Fraction as Fr


def arc_piece_solutions(T, rho):
    """all (g, h): F_7 -> {+1,-1}, g non-constant, such that x = g(fa(u)) on the tail fibre and h(fb(v)) on the head
    fibre satisfy the non-constant part of S x = rho x for the block [fa(u) - fb(v) in T] (tail positive, head negative).
    rho = 21:  h~ = -1/2 sum_t g~(w+t)  and  g~ = -(h~ * T);   rho = -21: h~ = sum_t g~(w+t) and h~ * T = 2 g~."""
    T = sorted(T); out = []
    for bits in range(1, 2 ** P - 1):
        g = [1 if (bits >> z) & 1 else -1 for z in range(P)]
        mg = Fr(sum(g), P); gt = [gz - mg for gz in g]
        conv = [sum(gt[(w + t) % P] for t in T) for w in range(P)]
        ht = [(-c / 2) if rho == 21 else c for c in conv]
        vals = sorted(set(ht))
        if len(vals) != 2 or vals[1] - vals[0] != 2: continue
        c0 = -1 - vals[0]; h = [int(v + c0) for v in ht]; mh = Fr(sum(h), P)
        hts = [hz - mh for hz in h]
        back = [sum(hts[(z - t) % P] for t in T) for z in range(P)]
        tail = all((-back[z] == gt[z]) if rho == 21 else (back[z] == 2 * gt[z]) for z in range(P))
        if tail: out.append((tuple(g), tuple(h), mg, mh))
    return out


def piece_members(choice, rho, limit_build=None):
    """enumerate all line-piece switching vectors x (x_inf = +1) with S x = rho x.  Returns list of
    (matching, per-arc (g,h) choice, constant signs) as dicts; exact quotient condition 7 e' m + 1 = rho m."""
    eps = np.zeros((9, 9), dtype=np.int64); arcinfo = {}
    for (a, b, fa, fb, T) in choice['arcs']:
        s = 1 if len(T) == 3 else -1; eps[a, b] = eps[b, a] = s
        arcinfo[(a, b)] = (fa, fb, frozenset(T))
    sols = {ab: arc_piece_solutions(info[2], rho) for ab, info in arcinfo.items()}
    # group by integer mean pair (7 m_tail, 7 m_head)
    groups = {}
    for ab, lst in sols.items():
        gdict = {}
        for s in lst: gdict.setdefault((int(s[2] * 7), int(s[3] * 7)), []).append(s)
        groups[ab] = gdict
    arcs = list(arcinfo); r7 = rho // 7
    found = []          # entries: (matching, mean-pairs, free signs, multiplicity)

    def check(M):       # eps' M + 1 = (rho/7) M, M = 7 m
        return all(sum(int(eps[i, j]) * M[j] for j in range(9) if j != i) + 1 == r7 * M[i] for i in range(9))

    def rec(k, used, chosen):
        free = [i for i in range(9) if i not in used]
        def assign(ci, M, picks):
            if ci == len(chosen):
                for signs in itertools.product((7, -7), repeat=len(free)):
                    MM = dict(M)
                    for i, s in zip(free, signs): MM[i] = s
                    if check([MM[i] for i in range(9)]):
                        mult = 1
                        for ab, mp in zip(chosen, picks): mult *= len(groups[ab][mp])
                        found.append({'matching': list(chosen), 'means': list(picks),
                                      'free': {i: s // 7 for i, s in zip(free, signs)}, 'mult': mult})
                return
            (a, b) = chosen[ci]
            for mp in groups[(a, b)]:
                M2 = dict(M); M2[a] = mp[0]; M2[b] = mp[1]
                assign(ci + 1, M2, picks + [mp])
        assign(0, {}, [])
        for t in range(k, len(arcs)):
            a, b = arcs[t]
            if a in used or b in used or not sols[(a, b)]: continue
            rec(t + 1, used | {a, b}, chosen + [(a, b)])
    rec(0, frozenset(), [])
    return found, arcinfo, groups


def expand_member(entry, groups, rng=None, index=0):
    """pick one concrete (g,h) per matched arc (the index-th combination, or random)"""
    picks = []
    for ab, mp in zip(entry['matching'], entry['means']):
        lst = groups[ab][mp]
        picks.append(lst[rng.randrange(len(lst))] if rng else lst[index % len(lst)])
    return {'matching': entry['matching'], 'picks': picks, 'free': entry['free']}


def member_vector(sol, arcinfo):
    x = np.ones(NV + 1, dtype=np.int64)
    for (a, b), (g, h, mg, mh) in zip(sol['matching'], sol['picks']):
        fa, fb, T = arcinfo[(a, b)]
        x[49 * a:49 * a + 49] = np.array(g)[form_val(fa)]
        x[49 * b:49 * b + 49] = np.array(h)[form_val(fb)]
    for i, s in sol['free'].items(): x[49 * i:49 * i + 49] = s
    return x


def switched_graph(A, x):
    S = seidel_plus_inf(A); S2 = S * x[:, None] * x[None, :]
    B = (1 - S2) // 2; np.fill_diagonal(B, 0); return B.astype(np.uint8)


def verify_srg_general(B):
    B = B.astype(np.int64); n = len(B)
    if not (B == B.T).all() or B.trace(): return None
    deg = set(B.sum(1).tolist())
    if len(deg) != 1: return None
    k = deg.pop(); B2 = B @ B
    lam = set(B2[B == 1].tolist()); off = (B == 0) & ~np.eye(n, dtype=bool); mu = set(B2[off].tolist())
    if len(lam) == 1 and len(mu) == 1: return (n, k, lam.pop(), mu.pop())
    return None


def descendant(S, w):
    d = S[w].copy(); d[w] = 1
    Sw = S * d[:, None] * d[None, :]
    keep = [i for i in range(len(S)) if i != w]
    B = (Sw[np.ix_(keep, keep)] == -1).astype(np.uint8)
    return B


def certset(A, points=None, log=None):
    S = seidel_plus_inf(A); out = {}
    pts = range(len(S)) if points is None else points
    for w in pts:
        out[w] = graph_cert(descendant(S, w))
        if log and w % 50 == 0: print('   descendant', w, flush=True)
    return out


def form_of_line(k):
    dx, dy = DIRV[k]; return (dy % P, (-dx) % P)


def random_regular_tournament(rng, steps=400):
    out = {i: {(i + d) % 9 for d in (1, 2, 3, 4)} for i in range(9)}
    for _ in range(steps):
        a = rng.randrange(9); b = rng.choice(sorted(out[a])); cs = [c for c in out[b] if a in out[c]]
        if cs:
            c = rng.choice(cs); out[a].discard(b); out[b].discard(c); out[c].discard(a)
            out[b].add(a); out[c].add(b); out[a].add(c)
    return out


def perturb(choice, knobs, rng):
    """random template choice derived from `choice`; knobs subset of
    {'phase','scale','slope','bij','fibre','tour'} ('fibre'/'tour' rebuild the form assignment from scratch)."""
    L = list(choice['L']); arcs = [list(a) for a in choice['arcs']]
    if 'fibre' in knobs or 'tour' in knobs:
        if 'fibre' in knobs: L = [frozenset(rng.sample(range(8), 4)) for _ in range(9)]
        if 'tour' in knobs:
            out = random_regular_tournament(rng); pairs = [(a, b) for a in range(9) for b in out[a]]
        else: pairs = [(a[0], a[1]) for a in arcs]
        posl = {i: rng.sample([k for k in range(8) if k not in L[i]], 4) for i in range(9)}
        negl = {i: rng.sample(sorted(L[i]), 4) for i in range(9)}
        arcs = []
        for (a, b) in pairs:
            fa = form_of_line(posl[a].pop()); fb = form_of_line(negl[b].pop())
            arcs.append([a, b, fa, fb, QR if PALEY9[a, b] else frozenset({0}) | NR])
    if 'bij' in knobs:
        for i in range(9):
            tails = [r for r in arcs if r[0] == i]; fs = [r[2] for r in tails]; rng.shuffle(fs)
            for r, f in zip(tails, fs): r[2] = f
            heads = [r for r in arcs if r[1] == i]; fs = [r[3] for r in heads]; rng.shuffle(fs)
            for r, f in zip(heads, fs): r[3] = f
    for r in arcs:
        if 'slope' in knobs:
            s = rng.randrange(1, P); r[3] = ((s * r[3][0]) % P, (s * r[3][1]) % P)
        if 'scale' in knobs:
            c = rng.randrange(1, P); r[4] = frozenset((c * t) % P for t in r[4])
        if 'phase' in knobs:
            t0 = rng.randrange(P); r[4] = frozenset((t + t0) % P for t in r[4])
    return {'L': L, 'arcs': [tuple(r) for r in arcs], 'name': choice['name'] + '+' + '/'.join(sorted(knobs))}


def choice_json(ch):
    return {'name': ch['name'], 'L': [sorted(l) for l in ch['L']],
            'arcs': [[a, b, list(fa), list(fb), sorted(T)] for (a, b, fa, fb, T) in ch['arcs']]}


def choice_from_json(d):
    return {'name': d['name'], 'L': [frozenset(l) for l in d['L']],
            'arcs': [(a, b, tuple(fa), tuple(fb), frozenset(T)) for (a, b, fa, fb, T) in d['arcs']]}


if __name__ == '__main__':
    mode = sys.argv[1]
    root = '.'                      # run from koval/
    OUT = root + '/data/P2'
    if mode == 'base':             # base ours|mathon OUT.g6 : the descendant at infinity (srg(441,220,109,110)) as graph6
        ch = ours_choice() if sys.argv[2] == 'ours' else mathon_choice()
        A = build(ch); ok = verify(A); assert ok[0], ok
        s = g6(A); open(sys.argv[3], 'w', newline='\n').write(s + '\n')
        print(ch['name'], ok, 'template_ok', template_ok(ch), 'sha256(file)', sha(s + '\n'), flush=True)
    os.makedirs(OUT, exist_ok=True)
    if mode == 'certset':          # certset NAME  (ours | mathon | path to a choice json)
        name = sys.argv[2]
        ch = ours_choice() if name == 'ours' else mathon_choice() if name == 'mathon' else choice_from_json(json.load(open(name)))
        A = build(ch); ok = verify(A); print(ch['name'], ok, flush=True); assert ok[0]
        # Aut(G) fixes infinity inside Aut(two-graph), so its orbit representatives (+ infinity) meet every
        # Aut(two-graph)-orbit on the 442 points: the descendant classes are the certs at these points.
        n = len(A); adj = {i: [int(j) for j in np.nonzero(A[i])[0]] for i in range(n)}
        orb = pynauty.autgrp(pynauty.Graph(n, adjacency_dict=adj))[3]
        reps = sorted({int(o) for o in orb}) + [n]
        print('orbit representatives', len(reps), flush=True)
        t = time.time(); cs = certset(A, points=reps, log=True)
        tag = name if name in ('ours', 'mathon') else os.path.basename(name)[:-5]
        json.dump({'name': ch['name'], 'certs': cs}, open('%s/certset_%s.json' % (OUT, tag), 'w'))
        print('done', tag, 'distinct descendant classes', len(set(cs.values())), '%.0fs' % (time.time() - t), flush=True)
    if mode == 'members':          # members ours|mathon|choice.json RHO NPERPATTERN SEED [ref.g6 ...]
        name, rho, nper, seed = sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
        ch = ours_choice() if name == 'ours' else mathon_choice() if name == 'mathon' else choice_from_json(json.load(open(name)))
        refs = {}
        for f in sys.argv[6:]:
            for ln, line in enumerate(open(f)):
                if line.strip():
                    B = load_g6(line); refs[graph_cert(B)] = '%s:%d' % (os.path.basename(f), ln + 1)
                    Bc = (1 - B) ; np.fill_diagonal(Bc, 0); refs[graph_cert(Bc.astype(np.uint8))] = '%s:%d(compl)' % (os.path.basename(f), ln + 1)
        print('reference certs', len(refs), flush=True)
        A = build(ch); assert verify(A)[0]
        found, arcinfo, groups = piece_members(ch, rho); rng = random.Random(seed)
        print(ch['name'], 'rho', rho, 'patterns', len(found), 'labelled', sum(e['mult'] for e in found), flush=True)
        seen = {}; tag = name if name in ('ours', 'mathon') else os.path.basename(name)[:-5]
        for pi, e in enumerate(found):
            for r in range(nper):
                sol = expand_member(e, groups, None if r == 0 else rng, 0)
                B = switched_graph(A, member_vector(sol, arcinfo)); par = verify_srg_general(B)
                if par is None: print('  pattern', pi, 'NOT SRG', flush=True); continue
                c = graph_cert(B)
                if c not in seen:
                    seen[c] = (pi, r, par, graph_aut(B), refs.get(c, 'new'))
                    with open('%s/members_%s_rho%d.g6' % (OUT, tag, rho), 'a') as fh: fh.write(g6(B) + '\n')
            print('  pattern', pi, 'matching', e['matching'], 'distinct so far', len(seen), flush=True)
        for c, v in seen.items(): print('  ', c[:12], v, flush=True)
        print('SUMMARY', tag, rho, 'distinct graphs', len(seen), 'matching Plummer/ref', sum(1 for v in seen.values() if v[4] != 'new'), flush=True)
    if mode == 'members0':         # members0 ours|mathon RHO OUTTAG : the UNTRANSLATED member of every pattern
        name, rho, otag = sys.argv[2], int(sys.argv[3]), sys.argv[4]   # (g, h minus-sets unions of {0}, Q, N)
        ch = ours_choice() if name == 'ours' else mathon_choice()
        A = build(ch); found, arcinfo, groups = piece_members(ch, rho)
        cls = [frozenset({0}), QR, NR]
        unions = {frozenset().union(*[c for c, b in zip(cls, bits) if b]) for bits in itertools.product((0, 1), repeat=3)}
        def untranslated(s):
            return frozenset(z for z in range(P) if s[0][z] == -1) in unions and frozenset(z for z in range(P) if s[1][z] == -1) in unions
        with open('%s/%s.g6' % (OUT, otag), 'w') as fh:
            for pi, e in enumerate(found):
                picks = []
                for ab, mp in zip(e['matching'], e['means']):
                    u = [s for s in groups[ab][mp] if untranslated(s)]; picks.append(u[0] if u else None)
                if None in picks: print('pattern', pi, 'no untranslated pick', flush=True); continue
                B = switched_graph(A, member_vector({'matching': e['matching'], 'picks': picks, 'free': e['free']}, arcinfo))
                print('pattern', pi, verify_srg_general(B), flush=True); fh.write(g6(B) + '\n')
    if mode == 'refcheck':         # refcheck CAND.g6 REF.g6 [REF2.g6 ...]: exact srg check + iso vs references (+complements)
        refs = {}
        for f in sys.argv[3:]:
            for ln, line in enumerate(open(f)):
                if line.strip():
                    B = load_g6(line); refs[graph_cert(B)] = '%s:%d' % (os.path.basename(f), ln + 1)
                    Bc = (1 - B); np.fill_diagonal(Bc, 0); refs[graph_cert(Bc.astype(np.uint8))] = '%s:%d(compl)' % (os.path.basename(f), ln + 1)
        seen = set(); new = 0
        for ln, line in enumerate(open(sys.argv[2])):
            if not line.strip(): continue
            B = load_g6(line); par = verify_srg_general(B); c = graph_cert(B)
            dup = c in seen; seen.add(c); new += (c not in refs and not dup)
            print('cand', ln + 1, par, 'sha256(g6)', sha(line.strip())[:16], 'cert', c[:12], 'dup' if dup else '', refs.get(c, 'NOT in refs'), flush=True)
        print('SUMMARY candidates', ln + 1, 'pairwise distinct', len(seen), 'not isomorphic to any ref', new, flush=True)
    if mode == 'compare':          # compare certset_*.json  samples_*.json
        sets = {}
        for f in sorted(os.listdir(OUT)):
            if f.startswith('certset_') and f.endswith('.json'):
                d = json.load(open(OUT + '/' + f)); sets[f[8:-5]] = set(d['certs'].values())
        for a in sets: print('two-graph', a, 'descendant classes', len(sets[a]), flush=True)
        names = list(sets)
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                print('  ', names[i], 'vs', names[j], 'common classes', len(sets[names[i]] & sets[names[j]]))
        for f in sorted(os.listdir(OUT)):
            if f.startswith('samples_') and f.endswith('.json'):
                for r in json.load(open(OUT + '/' + f)):
                    hit = [a for a in sets if r['cert'] in sets[a]]
                    print(r['name'], 'srg', r['srg'], 'aut(G)', r['aut'], '-> two-graph', hit or 'NEW (not in %s)' % names)
    if mode == 'tourcount':        # tourcount SEED N : line-piece pattern counts vs the signed tournament
        seed, n = int(sys.argv[2]), int(sys.argv[3]); rng = random.Random(seed)
        for s in range(n):
            ch = perturb(ours_choice(), {'tour', 'fibre', 'bij', 'phase', 'scale', 'slope'} if s else {'fibre', 'bij', 'phase', 'scale', 'slope'}, rng)
            out = {}
            for a, b, *_ in ch['arcs']: out.setdefault(a, []).append(b)
            # tournament invariant: score sequence of 'positive' (|T|=3) out-arcs, and number of cyclic triangles by sign type
            sc = sorted(sum(1 for b in out[a] if PALEY9[a, b]) for a in range(9))
            cnt = []
            for rho in (21, -21):
                found, arcinfo, groups = piece_members(ch, rho)
                cnt.append((len(found), sum(e['mult'] for e in found)))
            print('sample', s, 'tour-out-edge-scores', sc, 'patterns/labelled rho=21:', cnt[0], 'rho=-21:', cnt[1],
                  json.dumps({str(a): sorted(out[a]) for a in range(9)}), flush=True)
    if mode == 'samples':          # samples SEED N
        seed, nsamp = int(sys.argv[2]), int(sys.argv[3]); rng = random.Random(seed)
        knobsets = [{'phase'}, {'scale'}, {'slope'}, {'bij'}, {'phase', 'scale', 'slope'}, {'fibre'}, {'tour'},
                    {'fibre', 'tour', 'bij', 'phase', 'scale', 'slope'}]
        os.makedirs(OUT + '/samples', exist_ok=True); res = []
        for s in range(nsamp):
            base = ours_choice() if s % 2 == 0 else mathon_choice()
            ch = perturb(base, knobsets[(s // 2) % len(knobsets)], rng); ch['name'] += '#%d_%d' % (seed, s)
            okt = template_ok(ch); A = build(ch); ok = verify(A)
            c = graph_cert(A); aut = graph_aut(A)
            fn = '%s/samples/s%d_%d.json' % (OUT, seed, s); json.dump(choice_json(ch), open(fn, 'w'))
            res.append({'file': fn, 'name': ch['name'], 'template_ok': okt[0], 'srg': ok[0], 'cert': c, 'aut': aut})
            print(ch['name'], okt, ok, 'aut(G)', aut, c[:12], flush=True)
        json.dump(res, open('%s/samples_%d.json' % (OUT, seed), 'w'), indent=1)
    if mode == 'timing':
        A = build(ours_choice()); print('built', verify(A), flush=True)
        t = time.time(); print('graph cert 441', graph_cert(A)[:16], '%.1fs' % (time.time() - t), flush=True)
        t = time.time(); print('graph aut 441', graph_aut(A), '%.1fs' % (time.time() - t), flush=True)
    if mode == 'pieces':
        for ch in (ours_choice(), mathon_choice()):
            A = build(ch)
            for rho in (21, -21):
                t = time.time(); found, arcinfo, groups = piece_members(ch, rho)
                tot = sum(e['mult'] for e in found)
                print(ch['name'], 'rho', rho, 'mean-patterns', len(found), 'labelled line-piece vectors (x_inf=+1)', tot,
                      'per-arc solution counts', sorted(set(sum(len(v) for v in g.values()) for g in groups.values())),
                      '%.1fs' % (time.time() - t), flush=True)
                hist = {}
                for s in found: hist[len(s['matching'])] = hist.get(len(s['matching']), 0) + s['mult']
                print('   by matching size', hist, flush=True)
                rng = random.Random(1)
                for e in found[:2] + found[-2:]:
                    x = member_vector(expand_member(e, groups, rng), arcinfo); B = switched_graph(A, x)
                    print('   check one:', len(e['matching']), verify_srg_general(B), flush=True)
    if mode == 'controls':
        # mathon441_plummer442.g6: Mathon's instance as given by D. Plummer (file data/mathon441.g6 of his repository
        # github.com/dylan-plummer/srg442); NOT the SageMath graph used by scripts/q23/c1_controls.py.
        # srg_441_220_109_110_d1_z9.g6: the Paley-fibre graph of [K833, Computation 6.7]; written by
        # `base ours data/srg_441_220_109_110_d1_z9.g6` (SHA-256 0356d1e21c3ff601157935dc52d380384314fe196d10af7585d67a1dd7fb6c80).
        mfile = root + '/data/mathon441_plummer442.g6'
        ofile = root + '/data/srg_441_220_109_110_d1_z9.g6'
        for ch, f in ((mathon_choice(), mfile), (ours_choice(), ofile)):
            t0 = time.time(); A = build(ch)
            print(ch['name'], 'template_ok', template_ok(ch), 'verify', verify(A), flush=True)
            B = load_g6(open(f).readline())
            print('  equal to file as LABELLED graphs:', bool((A == B).all()), flush=True)
            print('  |Aut(G)|, orbits', graph_aut(A), flush=True)
