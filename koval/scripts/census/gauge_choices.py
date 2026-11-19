"""Single-knob template choices on the Paley-fibre two-graph ('ours'): the files data/template_data/gauge/*.json.
Expected gauge (isomorphic two-graph by construction), files g_*.json: fibre translation, fibre scalar, fibre GL(2,7)
map, fibre relabelling by an automorphism of Paley(9), and a combination.  Non-gauge, files n_*.json: ONE arc changed
(phase, Q<->N scaling, head-form slope by a square / non-square, tail-form swap at one fibre).
The sixteen files data/template_data/samples/s1_*.json come from  P2_family442.py samples 1 16  (seed 1, 16 samples).
Run from koval/:  python3 -I scripts/census/gauge_choices.py [OUTDIR]   (default data/template_data/gauge; deterministic)"""
import sys, os, json, random
sys.path.insert(0, 'scripts/census')
import P2_family442 as T
P = 7


def mat_inv(g):
    a, b, c, d = g[0][0], g[0][1], g[1][0], g[1][1]; det = (a * d - b * c) % P; di = pow(det, P - 2, P)
    return [[d * di % P, -b * di % P], [-c * di % P, a * di % P]]


def apply_gl(ch, i, g):
    gi = mat_inv(g)
    def fnew(f): return ((f[0] * gi[0][0] + f[1] * gi[1][0]) % P, (f[0] * gi[0][1] + f[1] * gi[1][1]) % P)
    L = list(ch['L'])
    L[i] = frozenset(T.line_of_vec((g[0][0] * T.DIRV[k][0] + g[0][1] * T.DIRV[k][1], g[1][0] * T.DIRV[k][0] + g[1][1] * T.DIRV[k][1])) for k in ch['L'][i])
    arcs = []
    for (a, b, fa, fb, Tt) in ch['arcs']:
        arcs.append((a, b, fnew(fa) if a == i else fa, fnew(fb) if b == i else fb, Tt))
    return {'L': L, 'arcs': arcs}


def translate(ch, i, u0):
    arcs = []
    for (a, b, fa, fb, Tt) in ch['arcs']:
        if a == i: Tt = frozenset((t + fa[0] * u0[0] + fa[1] * u0[1]) % P for t in Tt)
        if b == i: Tt = frozenset((t - fb[0] * u0[0] - fb[1] * u0[1]) % P for t in Tt)
        arcs.append((a, b, fa, fb, Tt))
    return {'L': list(ch['L']), 'arcs': arcs}


def paley9_aut_perm():
    # F_9 = F_3[w], index 3x+y <-> x + y w ... any map preserving K3 x K3: swap coordinates and translate
    return [3 * ((v % 3 + 1) % 3) + (v // 3) for v in range(9)]


def relabel(ch, s):
    L = [None] * 9
    for i in range(9): L[s[i]] = ch['L'][i]
    arcs = [(s[a], s[b], fa, fb, Tt) for (a, b, fa, fb, Tt) in ch['arcs']]
    return {'L': L, 'arcs': arcs}


def one_arc(ch, fn):
    arcs = list(ch['arcs']); arcs[0] = fn(arcs[0]); return {'L': list(ch['L']), 'arcs': arcs}


def choices():
    base = T.ours_choice(); out = {}
    out['g_translate_f0'] = translate(base, 0, (1, 2))
    out['g_scalar3_f0'] = apply_gl(base, 0, [[3, 0], [0, 3]])
    out['g_gl_f0'] = apply_gl(base, 0, [[1, 2], [3, 5]])
    out['g_relabel_paley9'] = relabel(base, paley9_aut_perm())
    out['g_combo'] = translate(apply_gl(apply_gl(base, 4, [[2, 1], [1, 1]]), 7, [[0, 1], [1, 0]]), 2, (3, 3))
    out['n_phase_arc0'] = one_arc(base, lambda r: (r[0], r[1], r[2], r[3], frozenset((t + 1) % P for t in r[4])))
    out['n_scaleQN_arc0'] = one_arc(base, lambda r: (r[0], r[1], r[2], r[3], frozenset((3 * t) % P for t in r[4])))
    out['n_slope_sq_arc0'] = one_arc(base, lambda r: (r[0], r[1], r[2], (2 * r[3][0] % P, 2 * r[3][1] % P), r[4]))
    out['n_slope_nsq_arc0'] = one_arc(base, lambda r: (r[0], r[1], r[2], (3 * r[3][0] % P, 3 * r[3][1] % P), r[4]))
    # tail-form swap at fibre 0 between its first two out-arcs
    arcs = list(base['arcs']); idx = [k for k, r in enumerate(arcs) if r[0] == 0][:2]
    r0, r1 = arcs[idx[0]], arcs[idx[1]]
    arcs[idx[0]] = (r0[0], r0[1], r1[2], r0[3], r0[4]); arcs[idx[1]] = (r1[0], r1[1], r0[2], r1[3], r1[4])
    out['n_bijswap_f0'] = {'L': list(base['L']), 'arcs': arcs}
    for name, ch in out.items(): ch['name'] = 'ours+' + name
    return out


if __name__ == '__main__':
    OUTD = sys.argv[1] if len(sys.argv) > 1 else 'data/template_data/gauge'
    os.makedirs(OUTD, exist_ok=True)
    for name, ch in choices().items():
        ok = T.template_ok(ch)
        json.dump(T.choice_json(ch), open('%s/%s.json' % (OUTD, name), 'w'))
        print(name, 'template_ok', ok, flush=True)
