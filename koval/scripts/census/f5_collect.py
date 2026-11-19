"""Group the two-graphs by the f5_invariant.py invariants; exact verdicts only.
Run from koval/:  python3 -I scripts/census/f5_collect.py TAG1 TAG2 ...
(reads data/f5/c2_TAG.json; pynauty |Aut(G_inf)| for the non-discrete ones)"""
import sys, json, itertools
sys.path.insert(0, 'scripts/census')
import numpy as np, pynauty
import P2_family442 as T
SCR = 'data/f5'
tags = sys.argv[1:]; R = {t: json.load(open('%s/c2_%s.json' % (SCR, t))) for t in tags}


def key(r):     # complete invariant if discrete; else the stable colour multiset (invariant, maybe incomplete)
    return ('canon', r['canonical_seidel_sha']) if r['discrete'] else ('multiset', r['stable_colour_multiset'])


groups = {}
for t in tags: groups.setdefault(key(R[t]), []).append(t)
print('GROUPS', len(groups))
for k, v in groups.items(): print('  ', k[0], k[1][:16], v)
# cross point-colour sets at fixed rounds (shared descendant <=> equal colour at every round)
for rnd in ('0', '3'):
    sets = {t: set(R[t]['point_colours_round'][rnd]) for t in tags}
    bad = [(a, b, len(sets[a] & sets[b])) for a, b in itertools.combinations(tags, 2) if sets[a] & sets[b] and key(R[a]) != key(R[b])]
    print('round', rnd, 'pairs in different groups sharing a point colour:', bad)
for t in tags:
    r = R[t]; line = {'tag': t, 'name': r['name'], 'discrete': r['discrete'], 'classes_per_round': r['classes_per_round'],
                      'inf_class_size': r['inf_class_size'], 'k4_desc_totals_minmax': r['k4_desc_totals_minmax']}
    if not r['discrete'] and len(sys.argv) > 1 and 'ctl' not in t:
        src = r.get('src')
        ch = T.ours_choice() if src == 'ours' else T.mathon_choice() if src == 'mathon' else T.choice_from_json(json.load(open(src)))
        A = T.build(ch); n = len(A); adj = {i: [int(j) for j in np.nonzero(A[i])[0]] for i in range(n)}
        g, mm, e, orb, no = pynauty.autgrp(pynauty.Graph(n, adjacency_dict=adj))
        line['aut_G_inf'] = int(round(mm * 10 ** e)); line['orbits_G_inf_plus_inf'] = no + 1
        line['stable_classes'] = len(set(r['stable_colours']))
    print(json.dumps(line), flush=True)
