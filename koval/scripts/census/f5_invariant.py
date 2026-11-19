"""Driver for f5.c: a nauty-free invariant of a regular two-graph on 442 points (first of the two two-graph invariants).

  1. F(x,y) = #K4 through y in the descendant D_x = #5-sets {x,y,a,b,c} with all 10 triples coherent (C program f5;
     symmetric by definition -> checked).  F is an isomorphism invariant of the two-graph (edge-coloured K_442).
  2. Colour refinement on (points, F): c_0(x) = sha(sorted row F(x,.)), c_{t+1}(x) = sha(c_t(x), sorted (F(x,y), c_t(y))).
     Colours are labelling-independent strings, so they are comparable ACROSS two-graphs: an isomorphism Omega -> Omega'
     mapping x to y forces c_t(x) = c'_t(y) for every t.  Since D_x determines (Omega, x), equal descendants => equal
     colours; different colours at any fixed round => non-isomorphic descendants.
  3. If the stable colouring is discrete, the order by colour is canonical; switching-normalise (row of the first
     point all +1) -> canonical Seidel matrix -> sha = COMPLETE two-graph invariant; Aut(two-graph) = 1 and the 442
     descendants are pairwise non-isomorphic.
Run from koval/ (f5.c compiled to bin/f5, e.g. gcc -O2 -o bin/f5 scripts/census/f5.c):
  python3 -I scripts/census/f5_invariant.py run NAME TAG
        NAME = ours | mathon | a template-data json (data/template_data/...); writes data/f5/c2_TAG.json, F_TAG.npy
  python3 -I scripts/census/f5_invariant.py control NAME SEED TAG
        random re-root at a point + relabel + switching of NAME, then the same (positive control)
Then scripts/census/f5_collect.py groups the results.  ('ours' is the Paley-fibre two-graph, 'mathon' is Mathon's.)
"""
import sys, os, json, time, hashlib, subprocess
sys.path.insert(0, 'scripts/census')
import numpy as np
import P2_family442 as T

SCR = 'data/f5'; F5 = 'bin/f5'
os.makedirs(SCR, exist_ok=True)
ROUNDS_FIXED = 3


def h(s): return hashlib.sha256(s.encode()).hexdigest()[:24]


def fmat(S):
    n = len(S)
    inp = '%d\n' % n + '\n'.join(''.join('1' if S[i, j] == -1 else '0' for j in range(n)) for i in range(n)) + '\n'
    out = subprocess.run([F5], input=inp, capture_output=True, text=True, check=True).stdout
    F = np.array([[int(v) for v in l.split()] for l in out.strip().split('\n')], dtype=np.int64)
    assert F.shape == (n, n)
    return F


def refine(F):
    n = len(F); off = ~np.eye(n, dtype=bool)
    c = [h('r0|' + ','.join(map(str, sorted(F[x][off[x]].tolist())))) for x in range(n)]
    hist = [c]; t = 0
    while True:
        t += 1
        nc = [h(c[x] + '|' + ';'.join(sorted('%d:%s' % (F[x, y], c[y]) for y in range(n) if y != x))) for x in range(n)]
        hist.append(nc)
        if len(set(nc)) == len(set(c)) and t >= ROUNDS_FIXED: break
        c = nc
    return hist


def analyse(S, name):
    n = len(S); t0 = time.time()
    F = fmat(S); tf = time.time() - t0
    assert (F == F.T).all(), 'F not symmetric'
    hist = refine(F); stable = hist[-1]
    ncls = len(set(stable))
    res = {'name': name, 't_f5': round(tf, 1), 'F_min': int(F[~np.eye(n, dtype=bool)].min()), 'F_max': int(F.max()),
           'k4_desc_totals_minmax': [int(F.sum(1).min()) // 4, int(F.sum(1).max()) // 4],
           'rounds': len(hist) - 1, 'classes_per_round': [len(set(c)) for c in hist],
           'colour_multiset_round': {str(t): h(json.dumps(sorted(hist[t]))) for t in range(ROUNDS_FIXED + 1)},
           'stable_colour_multiset': h(json.dumps(sorted(stable))), 'discrete': ncls == n,
           'point_colours_round': {str(t): hist[t] for t in range(ROUNDS_FIXED + 1)}, 'stable_colours': stable}
    if ncls == n:
        order = sorted(range(n), key=lambda x: stable[x])
        Sc = S[np.ix_(order, order)]; d = Sc[0].copy(); d[0] = 1
        Sn = (Sc * d[:, None] * d[None, :]).astype(np.int8)
        res['canonical_seidel_sha'] = hashlib.sha256(Sn.tobytes()).hexdigest()
    sizes = {}
    for col in stable: sizes[col] = sizes.get(col, 0) + 1
    res['stable_class_sizes'] = sorted(sizes.values())
    res['inf_class_size'] = sizes[stable[n - 1]]
    return res, F


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'run':
        name, tag = sys.argv[2], sys.argv[3]
        ch = T.ours_choice() if name == 'ours' else T.mathon_choice() if name == 'mathon' else T.choice_from_json(json.load(open(name)))
        A = T.build(ch); ok = T.verify(A); assert ok[0], ok
        S = T.seidel_plus_inf(A)
        res, F = analyse(S, ch['name']); res.update({'src': name, 'template_ok': T.template_ok(ch)[0], 'srg_G_inf': ok[1]})
    elif mode == 'control':
        name, seed, tag = sys.argv[2], int(sys.argv[3]), sys.argv[4]
        ch = T.ours_choice() if name == 'ours' else T.mathon_choice() if name == 'mathon' else T.choice_from_json(json.load(open(name)))
        A = T.build(ch); S = T.seidel_plus_inf(A); N = len(S); rng = np.random.default_rng(seed)
        w = int(rng.integers(0, N - 1)); D = T.descendant(S, w)
        S2 = T.seidel_plus_inf(D); perm = rng.permutation(N); S3 = S2[np.ix_(perm, perm)]
        x = rng.choice([-1, 1], size=N); S4 = S3 * x[:, None] * x[None, :]
        res, F = analyse(S4, ch['name'] + '/reroot%d/perm/switch' % w); res.update({'src': name, 'seed': seed})
    np.save('%s/F_%s.npy' % (SCR, tag), F)
    json.dump(res, open('%s/c2_%s.json' % (SCR, tag), 'w'))
    print(tag, res['name'], 't_f5', res['t_f5'], 'classes/round', res['classes_per_round'], 'discrete', res['discrete'],
          'canon', res.get('canonical_seidel_sha', '-')[:16], 'stable', res['stable_colour_multiset'][:12],
          'inf class size', res['inf_class_size'], flush=True)
