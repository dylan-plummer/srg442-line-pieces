"""Check of the template data of the 23 two-graphs on 442 points (data/template_data/).  Needs numpy and pynauty
(imported by P2_family442.py; no pynauty call is made here).  About a minute on one core.

  1. the 16 files samples/s1_*.json are exactly what  P2_family442.py samples 1 16  generates (seed 1), and the 10 files
     gauge/*.json are exactly what gauge_choices.py generates (compared as parsed JSON);
  2. every one of the 26 files, and the two built-in choices 'ours' (Paley fibers, Omega_P) and 'mathon' (Omega_M),
     satisfies the template conditions checked by P2_family442.template_ok, and its graph G (the descendant at infinity)
     is srg(441,220,109,110) exactly (A^2 in int64) with S^2 = 441 I for the Seidel matrix S of G + an isolated point;
  3. the graph6 string of G for 'ours' has the SHA-256 of the file srg_441_220_109_110_d1_z9.g6 of [K833].
Run from koval/:  python3 -I scripts/checks/check_template_data.py      (prints ALL OK at the end)
"""
import sys, os, json, glob, random, hashlib
sys.path.insert(0, 'scripts/census')
import numpy as np
import P2_family442 as T
import gauge_choices as G

D = 'data/template_data'
D1Z9_SHA = '0356d1e21c3ff601157935dc52d380384314fe196d10af7585d67a1dd7fb6c80'
ok = True


def canon(d): return json.loads(json.dumps(d))


# 1. regeneration
rng = random.Random(1)
knobsets = [{'phase'}, {'scale'}, {'slope'}, {'bij'}, {'phase', 'scale', 'slope'}, {'fibre'}, {'tour'},
            {'fibre', 'tour', 'bij', 'phase', 'scale', 'slope'}]          # as in P2_family442.py, mode 'samples'
for s in range(16):
    base = T.ours_choice() if s % 2 == 0 else T.mathon_choice()
    ch = T.perturb(base, knobsets[(s // 2) % len(knobsets)], rng); ch['name'] += '#%d_%d' % (1, s)
    r = canon(T.choice_json(ch)) == json.load(open('%s/samples/s1_%d.json' % (D, s)))
    print('samples/s1_%d.json regenerated:' % s, r, ch['name'], flush=True); ok &= r
gauge = G.choices()
r = sorted(gauge) == sorted(os.path.basename(f)[:-5] for f in glob.glob(D + '/gauge/*.json'))
print('gauge file names as generated:', r); ok &= r
for name, ch in sorted(gauge.items()):
    r = canon(T.choice_json(ch)) == json.load(open('%s/gauge/%s.json' % (D, name)))
    print('gauge/%s.json regenerated:' % name, r, flush=True); ok &= r

# 2. template conditions and exact parameters
I442 = 441 * np.eye(442, dtype=np.int64)
srcs = [('ours', T.ours_choice()), ('mathon', T.mathon_choice())] + \
       [(os.path.relpath(f, D).replace(os.sep, '/'), T.choice_from_json(json.load(open(f))))
        for f in sorted(glob.glob(D + '/samples/*.json')) + sorted(glob.glob(D + '/gauge/*.json'))]
for name, ch in srcs:
    tok = T.template_ok(ch)
    A = T.build(ch); v = T.verify(A)
    S = T.seidel_plus_inf(A); np.fill_diagonal(S, 0); s2 = bool((S @ S == I442).all())
    r = tok[0] and v[0] and s2
    print('%-28s template_ok %s  %s  S^2 = 441 I: %s' % (name, tok, v[1], s2), flush=True); ok &= r
    if name == 'ours':
        h = hashlib.sha256((T.g6(A) + '\n').encode()).hexdigest()
        r = h == D1Z9_SHA
        print('   graph6 of ours = srg_441_220_109_110_d1_z9.g6 of [K833] (sha256 %s...):' % h[:16], r); ok &= r
print('checked', len(srcs), 'template choices (2 built-in + 26 files)')
print('ALL OK' if ok else 'MISMATCH')
