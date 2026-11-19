"""second program (nauty-free) for the P3 member classes.
Runs gcanon (C) on every graph (and complement), turns its canonical order into sha256 of the canonically ordered
adjacency matrix, and decides: within-file distinctness, cross-file disjointness, complements, matches with the
reference (untrusted, parsed only) files.  Positive control: random relabellings must reproduce the canonical sha.
usage, from koval/ (gcanon.c compiled to bin/gcanon; the four classify outputs in data/P3/):
       python3 -I scripts/census/members.py OUTDIR"""
import sys, os, json, hashlib, subprocess, random
sys.path.insert(0, 'scripts/census')
import numpy as np
import P2_family442 as T

P3 = 'data/P3'; GC = 'bin/gcanon'
GEN = ['ours_231_gen', 'ours_210_gen', 'mathon_231_gen', 'mathon_210_gen']
REF = ['known_mathon_210', 'known_mathon_231', 'plummer15_on_ours_231_idx0-4_14', 'plummer15_idx5-13_complements_210',
       'mathon_members14']
OUT = sys.argv[1]; os.makedirs(OUT, exist_ok=True)
# The reference files (the 29 members known before the census) are not part of this folder; any that is absent
# is skipped, and the census checks (1) and (2) run without it.
_absent = [r for r in REF if not os.path.exists('%s/%s.g6' % (P3, r))]
if _absent: print('reference files not found, skipped:', _absent, flush=True)
REF = [r for r in REF if r not in _absent]


_TRI = {}


def fast_load(s):
    if s[0] == "~": n = ((ord(s[1]) - 63) << 12) | ((ord(s[2]) - 63) << 6) | (ord(s[3]) - 63); body = s[4:]
    else: n = ord(s[0]) - 63; body = s[1:]
    b = np.frombuffer(body.encode(), dtype=np.uint8) - 63
    bits = np.unpackbits(b[:, None], axis=1)[:, 2:].reshape(-1)
    if n not in _TRI: _TRI[n] = np.tril_indices(n, -1)
    r, c = _TRI[n]; m = len(r); assert not bits[m:].any()
    A = np.zeros((n, n), dtype=np.uint8); A[c, r] = bits[:m]; A[r, c] = bits[:m]
    return A


def canon_file(path, comp=False):
    lines = [l.strip() for l in open(path) if l.strip()]
    out = subprocess.run([GC] + (['c'] if comp else []), input='\n'.join(lines) + '\n', capture_output=True, text=True,
                         check=True).stdout.strip().split('\n')
    assert len(out) == len(lines), (path, len(out), len(lines))
    res = []
    for line, o in zip(lines, out):
        head, ordtxt = o.split(' ord=')
        f = head.split(); ordv = np.array([int(v) for v in ordtxt.split()])
        A = fast_load(line)
        if comp: A = (1 - A); np.fill_diagonal(A, 0)
        assert sorted(ordv.tolist()) == list(range(len(A)))
        C = A[np.ix_(ordv, ordv)]
        sha = hashlib.sha256(np.packbits(C).tobytes()).hexdigest()
        res.append({'srg': f[3], 'r0': f[4][3:], 'ms': f[5][3:], 'cls': int(f[6][4:]), 'leaves': int(f[7][7:]),
                    'trunc': 'TRUNCATED' in head, 'sha': sha})
    return lines, res


if __name__ == "__main__":
    R = {}; rep = {}
    from concurrent.futures import ThreadPoolExecutor
    jobs = [(nm, cp) for nm in GEN for cp in (False, True)]
    with ThreadPoolExecutor(6) as ex:
        outs = list(ex.map(lambda j: canon_file("%s/%s.g6" % (P3, j[0]), j[1]), jobs))
    pre = {j: o for j, o in zip(jobs, outs)}
    for name in GEN:
        for comp in (False, True):
            tag = name + ("_COMPL" if comp else "")
            lines, res = pre[(name, comp)]; R[tag] = res
            shas = [r['sha'] for r in res]
            rep[tag] = {'graphs': len(res), 'distinct_canon': len(set(shas)), 'srg_params': sorted(set(r['srg'] for r in res)),
                        'discrete_without_individualisation': sum(r['cls'] == 442 for r in res),
                        'max_leaves': max(r['leaves'] for r in res), 'truncated': sum(r['trunc'] for r in res)}
            print(tag, json.dumps(rep[tag]), flush=True)
    for name in REF:
        for comp in (False, True):
            tag = name + ('_COMPL' if comp else '')
            lines, res = canon_file('%s/%s.g6' % (P3, name), comp); R[tag] = res
    idx = {}
    for tag in R:
        for i, r in enumerate(R[tag]): idx.setdefault(r['sha'], []).append((tag, i))
    # cross-file
    def inter(a, b): return sorted(set(r['sha'] for r in R[a]) & set(r['sha'] for r in R[b]))
    cross = {}
    pairs = [('ours_231_gen', 'mathon_231_gen'), ('ours_210_gen', 'mathon_210_gen'),
             ('ours_210_gen_COMPL', 'ours_231_gen'), ('ours_210_gen_COMPL', 'mathon_231_gen'),
             ('mathon_210_gen_COMPL', 'ours_231_gen'), ('mathon_210_gen_COMPL', 'mathon_231_gen'),
             ('ours_210_gen_COMPL', 'mathon_210_gen_COMPL'), ('ours_231_gen_COMPL', 'mathon_231_gen_COMPL')]
    for a, b in pairs:
        cross['%s & %s' % (a, b)] = len(inter(a, b)); print('common', a, b, cross['%s & %s' % (a, b)], flush=True)
    allgen = set()
    for t in ('ours_231_gen', 'ours_210_gen', 'mathon_231_gen', 'mathon_210_gen'):
        allgen |= {('231' if '231' in t else '210', r['sha']) for r in R[t]}
    union231 = {r['sha'] for t in ('ours_231_gen', 'mathon_231_gen', 'ours_210_gen_COMPL', 'mathon_210_gen_COMPL') for r in R[t]}
    print('distinct srg(442,231,120,121) classes in the union of the 231 files and complemented 210 files:', len(union231), flush=True)
    # references
    refrep = {}
    for name in REF:
        for comp in (False, True):
            tag = name + ('_COMPL' if comp else '')
            for i, r in enumerate(R[tag]):
                hits = [(t, j) for (t, j) in idx[r['sha']] if t.endswith('_gen') or t.endswith('_gen_COMPL')]
                hits = [h for h in hits if not h[0].endswith('_COMPL')]   # report matches against the gen files as written
                refrep['%s[%d]' % (tag, i)] = {'srg': r['srg'], 'gen_matches': hits}
                if hits: print('ref', tag, i, r['srg'], '->', hits, flush=True)
    matched = {(t, j) for v in refrep.values() for (t, j) in v['gen_matches']}
    for g in ('ours_231_gen', 'ours_210_gen', 'mathon_231_gen', 'mathon_210_gen'):
        m = sorted(j for (t, j) in matched if t == g)
        print('gen file', g, 'lines matching a reference (graph or complement):', m, 'not matching:', len(R[g]) - len(set(m)), flush=True)
        rep[g]['ref_matched_lines'] = m
    # positive control: random relabellings
    rng = random.Random(20261030); ctl = []
    for g in GEN:
        lines = [l.strip() for l in open('%s/%s.g6' % (P3, g)) if l.strip()]
        picks = rng.sample(range(len(lines)), 5)
        # include non-discrete ones (cls < 442) if any
        nd = [i for i, r in enumerate(R[g]) if r['cls'] < 442]
        picks += rng.sample(nd, min(3, len(nd)))
        tmp = '/tmp/ctl_%s.g6' % g
        with open(tmp, 'w') as fh:
            for i in picks:
                A = T.load_g6(lines[i]); p = list(range(len(A))); rng.shuffle(p); p = np.array(p)
                fh.write(T.g6(A[np.ix_(p, p)]) + '\n')
        _, res = canon_file(tmp)
        for i, r in zip(picks, res): ctl.append((g, i, R[g][i]['cls'], r['sha'] == R[g][i]['sha']))
    print('positive controls (file, line, cls, same canon):', ctl, 'ALL PASS' if all(c[3] for c in ctl) else 'FAIL', flush=True)
    json.dump({'summary': rep, 'cross': cross, 'union231': len(union231), 'refs': refrep, 'controls': ctl,
               'canon': {t: [r['sha'] for r in R[t]] for t in R}}, open(OUT + '/members_census.json', 'w'))
