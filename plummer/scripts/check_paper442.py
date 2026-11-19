"""D.P.'s checks of statements in the joint 442 paper (main.tex), 14 Nov 2026.

  A  Table tab:orders  (n = p q^2 + 1, p = q + 2 a square, member parameters)
  B  Omega_K data (base, tournament, fibres, forms, T) as printed in the paper; template conditions; srg; S^2 = p q^2 I
  C  operator facts used in the proof of Lemma decouple, block by block, on Omega_K, Omega_M and random template data:
       C1  S_ii (g o h) = -q (g o h) if ker h in L_i, +q otherwise   (g mean zero);   S_ii 1 = 0
       C2  S_ab 1 = eps_ab q 1;  S_{inf,i} = 1
       C3  S_ab (g o h_b) = 0 unless ker h_b = ker f^-_ab; then = -2q sum_{t in T} g(f^+_ab(u) - t)   (h_b = f^-_ab)
       C4  S_ba (g o h_a) = 0 unless ker h_a = ker f^+_ab; then = -2q sum_{t in T} g(f^-_ab(v) + t)   (h_a = f^+_ab)
       C5  |hat 1_T(chi)|^2 != 0 for every non-trivial character (convolution by 1_T injective on mean-zero functions)
  D  Lemma decouple (3) on line-piece members of Omega_K and Omega_M (both degrees), and S x = rho x
  E  Table tab:aut regenerated from plummer/data/census_props.tsv

Imports A.K.'s P2_family442.py read-only (not modified).  Run:
  python3 -I plummer/scripts/check_paper442.py        (numpy, pynauty)
"""
import sys, os, random, math, cmath, csv, ast
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))                     # repository root
sys.path.insert(0, os.path.join(ROOT, 'koval', 'scripts', 'census'))
import numpy as np
import P2_family442 as K2

P = 7; q = 7; NV = 441
fails = []


def ok(cond, msg):
    if not cond: fails.append(msg); print('FAIL', msg)


# ---------------- A ----------------
def is_prime_power(m):
    for r in range(2, m + 1):
        if m % r == 0:
            while m % r == 0: m //= r
            return m == 1
    return False


print('A  orders table')
rows = []
for qq in range(3, 101):
    if qq % 4 != 3 or not is_prime_power(qq): continue
    p = qq + 2; s = math.isqrt(p)
    if s * s != p: continue
    n = p * qq * qq + 1
    mem = []
    for rho in (qq * s, -qq * s):                    # rho = n - 1 - 2k
        k = (n - 1 - rho) // 2
        r, t = (-1 + qq * s) // 2, (-1 - qq * s) // 2  # restricted eigenvalues (-1 -+ rho')/2
        lam = k + r + t + r * t; mu = k + r * t
        assert k * (k - lam - 1) == (n - k - 1) * mu
        mem.append((n, k, lam, mu))
    rows.append((qq, p, n, sorted(mem, key=lambda z: z[1])))
    print('   q=%d p=%d n=%d  %s' % (qq, p, n, sorted(mem, key=lambda z: z[1])))
ok([r[0] for r in rows] == [7, 23, 47, 79], 'A: q list')
paper = {7: [(442, 210, 99, 100), (442, 231, 120, 121)],
         23: [(13226, 6555, 3248, 3249), (13226, 6670, 3363, 3364)],
         47: [(108242, 53956, 26895, 26896), (108242, 54285, 27224, 27225)],
         79: [(505522, 252405, 126024, 126025), (505522, 253116, 126735, 126736)]}
for qq, p, n, mem in rows: ok(paper[qq] == mem, 'A: row q=%d' % qq)

# ---------------- B ----------------
print('B  Omega_K data')
OKc = K2.ours_choice(); OMc = K2.mathon_choice()
Lk = OKc['L'][0]
print('   L_i (line indices; c = span(1,c), 7 = span(0,1)):', sorted(Lk))
ok(sorted(Lk) == [0, 3, 4, 7] and all(l == Lk for l in OKc['L']), 'B: Omega_K fibres = lines of slopes 0, 3, 4, inf')
ok(all(sorted(l) == [0, 1, 2, 4] for l in OMc['L']), 'B: Omega_M fibres = slopes {0} u Q')
for (a, b, fa, fb, T) in OKc['arcs']:
    d = (b - a) % 9; want = {1: ((1, 1), (1, 0)), 2: ((1, 3), (0, 1)), 3: ((1, 4), (1, 2)), 4: ((1, 6), (1, 5))}[d]
    ok((fa, fb) == want and sorted(T) == ([1, 2, 4] if K2.PALEY9[a, b] else [0, 3, 5, 6]), 'B: Omega_K arc %d->%d as printed' % (a, b))
for (a, b, fa, fb, T) in OMc['arcs']:
    d = (b - a) % 9; want = {1: ((2, 1), (5, 1)), 2: ((4, 1), (3, 1)), 3: ((1, 1), (6, 1)), 4: ((1, 0), (0, 1))}[d]
    ok((fa, fb) == want and sorted(T) == ([3, 5, 6] if K2.PALEY9[a, b] else [0, 1, 2, 4]), 'B: Omega_M arc %d->%d as printed' % (a, b))
for (a, b, fa, fb, T) in OKc['arcs'][:4]:
    print('   arc %d->%d  f+=%s (ker %d)  f-=%s (ker %d)  T=%s' % (a, b, fa, K2.kernel_line(fa), fb, K2.kernel_line(fb), sorted(T)))
tour = sorted(set((b - a) % 9 for (a, b, *_ ) in OKc['arcs']))
ok(tour == [1, 2, 3, 4], 'B: Omega_K tournament i -> i+d, d=1..4')
for ch in (OKc, OMc):
    ok(K2.template_ok(ch)[0], 'B: template_ok %s' % ch['name'])
    A = K2.build(ch); ok(K2.verify(A)[0], 'B: srg441 %s' % ch['name'])
    S = K2.seidel_plus_inf(A)
    ok((S @ S == 9 * 49 * np.eye(442, dtype=np.int64)).all(), 'B: S^2 %s' % ch['name'])
print('   Omega_K = Omega_M fibres?', sorted(OKc['L'][0]) == sorted(OMc['L'][0]))

# ---------------- C ----------------
print('C  operator facts')
PTS = K2.PTS
forms_by_line = {K2.kernel_line(f): f for f in [(1, 0), (0, 1)] + [(1, b) for b in range(1, P)]}


def fn(g, f):                     # g o f as a vector on the 49 points
    return np.array(g)[K2.form_val(f)]


def rand_mean0(rng):
    g = np.array([rng.randrange(-5, 6) for _ in range(P)], dtype=np.int64); return 7 * g - g.sum()  # integer, mean 0


def check_ops(ch, rng, tag):
    A = K2.build(ch); S = K2.seidel_plus_inf(A)
    blk = lambda i, j: S[49 * i:49 * i + 49, 49 * j:49 * j + 49]
    one = np.ones(49, dtype=np.int64); eps = {}
    for (a, b, fa, fb, T) in ch['arcs']:
        eps[(a, b)] = eps[(b, a)] = 1 if len(T) == 3 else -1
    ok((S[441, :441] == 1).all(), 'C2 inf row %s' % tag)
    for i in range(9):
        ok((blk(i, i) @ one == 0).all(), 'C1 S_ii 1 %s' % tag)
        for k in range(8):
            g = rand_mean0(rng); v = fn(g, forms_by_line[k])
            sig = -1 if k in ch['L'][i] else 1
            ok((blk(i, i) @ v == sig * q * v).all(), 'C1 %s fibre %d line %d' % (tag, i, k))
    for (a, b) in eps:
        if a != b: ok((blk(a, b) @ one == eps[(a, b)] * q * one).all(), 'C2 %s %d,%d' % (tag, a, b))
    for (a, b, fa, fb, T) in ch['arcs']:
        T = sorted(T); va = K2.form_val(fa); vb = K2.form_val(fb)
        for k in range(8):
            g = rand_mean0(rng)
            # C3: from head b to tail a
            h = forms_by_line[k]; img = blk(a, b) @ fn(g, h)
            if k == K2.kernel_line(fb):
                img2 = blk(a, b) @ fn(g, fb)
                want = -2 * q * np.array([sum(g[(va[u] - t) % P] for t in T) for u in range(49)])
                ok((img2 == want).all(), 'C3 %s arc %d->%d' % (tag, a, b))
            else: ok((img == 0).all(), 'C3 zero %s arc %d->%d line %d' % (tag, a, b, k))
            # C4: from tail a to head b
            img = blk(b, a) @ fn(g, h)
            if k == K2.kernel_line(fa):
                img2 = blk(b, a) @ fn(g, fa)
                want = -2 * q * np.array([sum(g[(vb[v] + t) % P] for t in T) for v in range(49)])
                ok((img2 == want).all(), 'C4 %s arc %d->%d' % (tag, a, b))
            else: ok((img == 0).all(), 'C4 zero %s arc %d->%d line %d' % (tag, a, b, k))
        w = cmath.exp(2j * math.pi / P)
        ok(min(abs(sum(w ** (c * t) for t in T)) ** 2 for c in range(1, P)) > 0.5, 'C5 %s' % tag)
    # cross blocks between fibres not joined by an arc do not exist (complete base)
    ok(len(ch['arcs']) == 36, 'C %s 36 arcs' % tag)


rng = random.Random(20261114)
check_ops(OKc, rng, 'Omega_K'); check_ops(OMc, rng, 'Omega_M')
for s in range(4):
    knobs = {'tour', 'fibre', 'bij', 'phase', 'scale', 'slope'} if s % 2 else {'fibre', 'bij', 'phase', 'scale', 'slope'}
    ch = K2.perturb(OKc if s < 2 else OMc, knobs, rng)
    ok(K2.template_ok(ch)[0], 'C template_ok random %d' % s)
    check_ops(ch, rng, 'random%d' % s)
print('   C done (Omega_K, Omega_M, 4 random template data)')

# ---------------- D ----------------
print('D  members: S x = rho x and the quotient equation (3)')
for ch in (OKc, OMc):
    A = K2.build(ch); S = K2.seidel_plus_inf(A)
    eps = np.zeros((9, 9), dtype=np.int64)
    for (a, b, fa, fb, T) in ch['arcs']: eps[a, b] = eps[b, a] = 1 if len(T) == 3 else -1
    for rho in (21, -21):
        found, arcinfo, groups = K2.piece_members(ch, rho)
        nchk = 0
        for e in rng.sample(found, min(12, len(found))):
            sol = K2.expand_member(e, groups, rng=rng); x = K2.member_vector(sol, arcinfo)
            ok((S @ x == rho * x).all(), 'D eigen %s %d' % (ch['name'], rho))
            M = [int(x[49 * i:49 * i + 49].sum()) // 7 for i in range(9)]           # M_i = sum_s g_i(s)
            for i in range(9):
                ok(q * (1 + sum(eps[i, j] * M[j] for j in range(9) if j != i)) == rho * M[i], 'D (3) %s' % ch['name'])
            ok(q * sum(M) == rho, 'D inf row %s' % ch['name'])
            nchk += 1
        print('   %s rho=%+d: %d member patterns, %d checked' % (ch['name'], rho, len(found), nchk))

# ---------------- E ----------------
print('E  aut table from census_props.tsv')
tsv = os.path.join(ROOT, 'plummer', 'data', 'census_props.tsv')
cnt = Counter()
with open(tsv) as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        src = row['source']; tg = 'M' if 'mathon' in src else 'K' if ('ours' in src or 'paley' in src.lower()) else src
        deg = '210' if '210' in src else '231' if '231' in src else '?'
        o = row['orbits']
        typ = o if not o.startswith('[') else ' '.join('%d^%d' % (k, v) for k, v in sorted(Counter(ast.literal_eval(o)).items()))
        cnt[(tg, deg, row['aut'], typ, bool(row['ours']))] += 1
for k in sorted(cnt): print('  ', k, cnt[k])
print('   total', sum(cnt.values()))

print('\nRESULT:', 'ALL PASS' if not fails else '%d FAILURES' % len(fails))
