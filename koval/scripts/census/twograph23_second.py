"""Second program for the pairwise non-isomorphism of the 23 template two-graphs on 442 points.

Builds each two-graph from its template data with our own builder (scripts/census/P2_family442.py; template data in
data/template_data/samples/s1_*.json and data/template_data/gauge/n_*.json, the same files as for f5_invariant.py),
checks the descendant at infinity exactly (srg(441,220,109,110) by A^2) and S^2 = 441 I on all entries, and runs descinv.c (a new C
program, independent of f5.c) on the Seidel matrix: one labeling-free hash I(D_x) per point x.  The invariant of the
two-graph is the multiset {I(D_x)}.  Different multisets => non-isomorphic two-graphs; disjoint sets => no descendant
of one is isomorphic to a descendant of the other.

Controls (same code, known answers):
  * gauge moves (gauge/g_*.json, isomorphic to Omega_K by construction) must give Omega_K's multiset;
  * re-rooted + randomly relabeled + randomly switched copies of Omega_M, Omega_K, s1_14 (trivial group),
    s1_3 (order 3) and n_phase_arc0 must give the multiset of their source.
Run from koval/ (descinv.c compiled to bin/descinv):
      python3 -I scripts/census/twograph23_second.py compute K M   (jobs with index = K mod M; 1 thread each)
      python3 -I scripts/census/twograph23_second.py merge         (reads the per-job files, writes
                                                                    data/twograph23/result.json)
Also computed: the complementary two-graphs -Omega_K, -Omega_M (Seidel matrix -S), to compare the 21 further
two-graphs with them.
"""
import sys, os, json, time, hashlib, subprocess
sys.path.insert(0, 'scripts/census')
import numpy as np
import P2_family442 as T

BIN = 'bin/descinv'
B2 = 'data/template_data'
OUT = 'data/twograph23'
os.makedirs(OUT, exist_ok=True)

TWOG = [('ours', 'ours'), ('mathon', 'mathon')] + \
       [('s1_%d' % i, '%s/samples/s1_%d.json' % (B2, i)) for i in range(16)] + \
       [(t, '%s/gauge/%s.json' % (B2, t)) for t in
        ('n_phase_arc0', 'n_scaleQN_arc0', 'n_bijswap_f0', 'n_slope_sq_arc0', 'n_slope_nsq_arc0')]
GAUGE = [(t, '%s/gauge/%s.json' % (B2, t)) for t in
         ('g_translate_f0', 'g_gl_f0', 'g_scalar3_f0', 'g_relabel_paley9', 'g_combo')]
CTRL = [('ctl_mathon', 'mathon', 12), ('ctl_ours', 'ours', 11), ('ctl_s1_14', 's1_14', 14), ('ctl_s1_3', 's1_3', 13),
        ('ctl_n_phase', 'n_phase_arc0', 7)]


def choice(src):
    if src == 'ours': return T.ours_choice()
    if src == 'mathon': return T.mathon_choice()
    return T.choice_from_json(json.load(open(src)))


def seidel(src):
    ch = choice(src)
    A = T.build(ch)
    ok = T.verify(A); assert ok[0], ok
    S = T.seidel_plus_inf(A)
    S2 = S @ S
    assert (S2 == 441 * np.eye(442, dtype=np.int64)).all(), 'S^2 != 441 I'
    return S, T.template_ok(ch)[0]


def run_inv(S, tag):
    cache = '%s/%s.txt' % (OUT, tag)
    if not os.path.exists(cache):
        n = len(S)
        inp = '%d\n' % n + '\n'.join(''.join('1' if S[i, j] == -1 else '0' for j in range(n)) for i in range(n)) + '\n'
        t = time.time()
        out = subprocess.run([BIN], input=inp, capture_output=True, text=True, check=True).stdout
        open(cache + '.tmp', 'w').write(out); os.replace(cache + '.tmp', cache)
        print('   %s: descinv %.0f s' % (tag, time.time() - t), flush=True)
    return read_inv(tag)


def read_inv(tag):
    rows = [l.split() for l in open('%s/%s.txt' % (OUT, tag)).read().split('\n') if l.strip()]
    assert len(rows) == 442 and [int(r[0]) for r in rows] == list(range(442))
    return [r[1] for r in rows], (min(int(r[2]) for r in rows), max(int(r[3]) for r in rows))


def summarise(hs):
    return {'multiset_sha': hashlib.sha256(' '.join(sorted(hs)).encode()).hexdigest()[:24],
            'distinct_descendant_hashes': len(set(hs)), 'hash_at_inf': hs[441],
            'inf_hash_multiplicity': hs.count(hs[441])}


NEG = [('neg_ours', 'ours'), ('neg_mathon', 'mathon')]


def jobs():
    J = [('plain', t, s) for t, s in TWOG + GAUGE] + [('neg', t, s) for t, s in NEG] + [('ctrl', t, (s, seed)) for t, s, seed in CTRL]
    return J


def src_of(tag):
    return dict(TWOG + GAUGE)[tag]


def compute(k, m):
    for idx, (kind, tag, src) in enumerate(jobs()):
        if idx % m != k: continue
        meta = {'tag': tag, 'kind': kind}
        if kind == 'plain':
            S, tok = seidel(src); meta.update(template_ok=tok, srg_and_S2_exact=True, src=src)
        elif kind == 'neg':
            S, tok = seidel(src_of(src)); S = -S; meta.update(of=src, srg_and_S2_exact=True)
        else:
            sname, seed = src
            S0, _ = seidel(src_of(sname)); N = len(S0); rng = np.random.default_rng(seed)
            w = int(rng.integers(0, N - 1)); D = T.descendant(S0, w)
            S2 = T.seidel_plus_inf(D); perm = rng.permutation(N); S3 = S2[np.ix_(perm, perm)]
            x = rng.choice([-1, 1], size=N); S = (S3 * x[:, None] * x[None, :]).astype(np.int64)
            meta.update(source=sname, reroot=w, seed=seed)
        run_inv(S, tag)
        json.dump(meta, open('%s/%s.json' % (OUT, tag), 'w'))
        print('done', tag, flush=True)


def merge():
    res = {}
    for kind, tag, src in jobs():
        meta = json.load(open('%s/%s.json' % (OUT, tag)))
        hs, wr = read_inv(tag)
        res[tag] = dict(summarise(hs), W_range=wr, hashes=hs, **meta)
        print(tag, {k: v for k, v in res[tag].items() if k != 'hashes'}, flush=True)
    names = [t for t, _ in TWOG]
    verdict = {'pairs_equal_multiset': [], 'pairs_sharing_a_descendant_hash': []}
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = res[names[i]], res[names[j]]
            if a['multiset_sha'] == b['multiset_sha']: verdict['pairs_equal_multiset'].append((names[i], names[j]))
            if set(a['hashes']) & set(b['hashes']): verdict['pairs_sharing_a_descendant_hash'].append((names[i], names[j]))
    verdict['gauge_equal_to_ours'] = {t: res[t]['multiset_sha'] == res['ours']['multiset_sha'] for t, _ in GAUGE}
    verdict['controls_equal_to_source'] = {t: res[t]['multiset_sha'] == res[s]['multiset_sha'] for t, s, _ in CTRL}
    verdict['all_23_pairwise_distinct'] = not verdict['pairs_equal_multiset']
    verdict['all_23_pairwise_disjoint'] = not verdict['pairs_sharing_a_descendant_hash']
    verdict['all_controls_pass'] = all(verdict['gauge_equal_to_ours'].values()) and all(verdict['controls_equal_to_source'].values())
    verdict['all_template_ok'] = all(res[t]['template_ok'] for t in names)
    verdict['distinct_descendant_hashes'] = {t: res[t]['distinct_descendant_hashes'] for t in names}
    # complements: the 21 further two-graphs against -Omega_K, -Omega_M; and +-Omega_K, +-Omega_M among themselves
    four = ['ours', 'mathon', 'neg_ours', 'neg_mathon']
    verdict['four_signed_pairwise_disjoint'] = all(not (set(res[a]['hashes']) & set(res[b]['hashes']))
                                                   for i, a in enumerate(four) for b in four[i + 1:])
    verdict['further_21_vs_neg_disjoint'] = {t: not (set(res[t]['hashes']) & (set(res['neg_ours']['hashes']) | set(res['neg_mathon']['hashes'])))
                                             for t in names if t not in ('ours', 'mathon')}
    verdict['all_21_disjoint_from_neg'] = all(verdict['further_21_vs_neg_disjoint'].values())
    print(json.dumps({k: v for k, v in verdict.items() if k != 'further_21_vs_neg_disjoint'}, indent=1), flush=True)
    json.dump({'verdict': verdict, 'per_twograph': res}, open(OUT + '/result.json', 'w'))


if __name__ == '__main__':
    if sys.argv[1] == 'compute': compute(int(sys.argv[2]), int(sys.argv[3]))
    else: merge()
