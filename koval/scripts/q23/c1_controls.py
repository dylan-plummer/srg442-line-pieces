"""Controls of the template code (w1lib.py) on known graphs.
(a) q = 3: literal Mathon blocks reproduce the printed (4.15) A_d, A_1, A_2 and the block pattern of A_I in (4.16);
    the 45-vertex graph is srg(45,22,10,11) and S^2 = 45 I.
(b) q = 7, p = 9: literal Mathon build and Paley-fibre template build; S^2 = 441 I exactly; template build compared
    byte-for-byte with our archived srg_441_220_109_110_d1_z9.g6; two-graph comparison (pynauty certificate of the
    double cover = complete two-graph invariant) with data/mathon441_sage_t2.g6 = SageMath's
    graphs.MathonPseudocyclicStronglyRegularGraph(t=2), written with graph6_string().  That graph is a third regular
    two-graph on 442 points (isomorphic to neither Mathon's literal build nor the Paley-fibre template build); it is
    NOT the file mathon441_plummer442.g6 used by scripts/census/P2_family442.py.
    Also needs data/srg_441_220_109_110_d1_z9.g6 (scripts/census/P2_family442.py base ours ...).
usage, from koval/: python3 -I scripts/q23/c1_controls.py      (W1_NAUTY=1 adds the slow pynauty part)
"""
import sys, os, json, hashlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pynauty
from w1lib import *

OUT = 'data/q23'; os.makedirs(OUT, exist_ok=True)
res = {}

def mat(rows):
    return np.array([[int(c) for c in r.replace(' ', '')] for r in rows], dtype=np.uint8)

Ad_pr = mat(['000 110 101', '000 011 110', '000 101 011', '101 000 110', '110 000 011', '011 000 101',
             '110 101 000', '011 110 000', '101 011 000'])
A1_pr = mat(['100 010 001', '010 001 100', '001 100 010', '010 001 100', '001 100 010', '100 010 001',
             '001 100 010', '100 010 001', '010 001 100'])
A2_pr = mat(['001 001 001'] * 3 + ['100 100 100'] * 3 + ['010 010 010'] * 3)
mb = mathon_blocks(3)
res['q3_printed_4.15'] = {'Ad': bool((mb['Ad'] == Ad_pr).all()), 'A1': bool((mb['Ak'][1] == A1_pr).all()),
                          'A2': bool((mb['Ak'][2] == A2_pr).all())}
print('q=3 (4.15) reproduced:', res['q3_printed_4.15'], flush=True)
# (4.16) A_I: pentagon 1-2-3-4-5, circulant skew-Latin square with first row (0,1,2,-2,-1)
C5 = np.zeros((5, 5), dtype=np.uint8)
for i in range(5): C5[i, (i+1) % 5] = C5[(i+1) % 5, i] = 1
A45, _ = build_mathon(3, C5)
J9 = np.ones((9, 9), dtype=np.uint8)
d = {'d': mb['Ad'], '1': mb['Ak'][1], '2': mb['Ak'][2], '1T': mb['Ak'][1].T, '2T': mb['Ak'][2].T}
d.update({'~' + k: J9 - v for k, v in list(d.items()) if k != 'd'})
AI = [['d', '1', '~2', '~2T', '1T'], ['1T', 'd', '1', '~2', '~2T'], ['~2T', '1T', 'd', '1', '~2'],
      ['~2', '~2T', '1T', 'd', '1'], ['1', '~2', '~2T', '1T', 'd']]
AIm = np.block([[d[x] for x in row] for row in AI])
res['q3_A_I_4.16'] = bool((AIm == A45).all())
res['q3_srg45'] = srg_check(A45)
res['q3_S2'] = seidel_square_check(A45)
print('q=3 A_I (4.16) equal:', res['q3_A_I_4.16'], res['q3_srg45'], res['q3_S2']['ok'], flush=True)

# ---------------- q = 7
F9 = GF(9); P9 = F9.paley_adj(); assert check_conference_graph(P9)
AM, _ = build_mathon(7, P9)
AT, st, arcs = build_template(7, P9)
for name, A in (('mathon_literal_q7', AM), ('template_q7', AT)):
    res[name] = {'S2': seidel_square_check(A), 'srg': srg_check(A)}
    print(name, res[name], flush=True)
ref = open('data/srg_441_220_109_110_d1_z9.g6').read().strip()
res['template_q7_equals_archived_d1_z9_g6'] = to_g6(AT) == ref
print('template q=7 g6 == archived srg_441_220_109_110_d1_z9.g6:', res['template_q7_equals_archived_d1_z9_g6'], flush=True)
if not os.path.exists('data/mathon441_sage_t2.g6'):   # optional SageMath control graph (see the docstring)
    print('data/mathon441_sage_t2.g6 not found: SageMath comparison skipped', flush=True)
    json.dump(res, open(OUT + '/c1_controls.json', 'w'), indent=1); sys.exit(0)
M441 = from_g6(open('data/mathon441_sage_t2.g6').read())
res['sage_mathon441_g6'] = {'srg': srg_check(M441), 'S2': seidel_square_check(M441)}
res['mathon_literal_q7_equals_sage_mathon441_g6_entrywise'] = bool((M441 == AM).all())
print('mathon441_sage_t2.g6:', res['sage_mathon441_g6'], ' entrywise equal to our literal Mathon build:',
      res['mathon_literal_q7_equals_sage_mathon441_g6_entrywise'], flush=True)
if os.environ.get('W1_NAUTY') != '1':      # the pynauty part (884-vertex double covers) is slow; optional
    json.dump(res, open(OUT + '/c1_controls.json', 'w'), indent=1); sys.exit(0)

def seidel_int(A):
    N = len(A); S = np.ones((N + 1, N + 1), dtype=np.int64); S[:N, :N] -= 2 * A.astype(np.int64)
    np.fill_diagonal(S, 0); return S

def cover_cert(S):
    N = len(S); adj = {}
    for x in range(N):
        neg = np.nonzero(S[x] == -1)[0]; pos = np.nonzero(S[x] == 1)[0]
        adj[2*x] = [2*int(y) for y in neg] + [2*int(y)+1 for y in pos]
        adj[2*x+1] = [2*int(y)+1 for y in neg] + [2*int(y) for y in pos]
    g = pynauty.Graph(2 * N, adjacency_dict=adj)
    return hashlib.sha256(pynauty.certificate(g)).hexdigest()

def graph_cert(A):
    n = len(A); g = pynauty.Graph(n, adjacency_dict={i: [int(j) for j in np.nonzero(A[i])[0]] for i in range(n)})
    return hashlib.sha256(pynauty.certificate(g)).hexdigest()

rng = np.random.default_rng(20261102)
certs = {}
for name, A in (('mathon_literal_q7', AM), ('template_q7', AT), ('sage_mathon441_g6', M441)):
    certs[name] = {'twograph_cover': cover_cert(seidel_int(A)), 'graph': graph_cert(A)}
# positive control for the invariant: random relabelling + random switching of the mathon literal Seidel matrix
S = seidel_int(AM); perm = rng.permutation(len(S)); x = rng.choice([-1, 1], size=len(S))
S2 = S[np.ix_(perm, perm)] * x[:, None] * x[None, :]
certs['mathon_literal_q7_relabel_switch'] = {'twograph_cover': cover_cert(S2)}
res['certs'] = certs
res['twograph_mathon_literal_eq_sage_mathon441'] = certs['mathon_literal_q7']['twograph_cover'] == certs['sage_mathon441_g6']['twograph_cover']
res['twograph_template_eq_sage_mathon441'] = certs['template_q7']['twograph_cover'] == certs['sage_mathon441_g6']['twograph_cover']
res['invariant_positive_control'] = certs['mathon_literal_q7_relabel_switch']['twograph_cover'] == certs['mathon_literal_q7']['twograph_cover']
res['graph_mathon_literal_eq_sage_mathon441'] = certs['mathon_literal_q7']['graph'] == certs['sage_mathon441_g6']['graph']
for k in ('mathon_literal_q7_equals_sage_mathon441_g6_entrywise', 'twograph_mathon_literal_eq_sage_mathon441', 'graph_mathon_literal_eq_sage_mathon441',
          'twograph_template_eq_sage_mathon441', 'invariant_positive_control'):
    print(k, res[k], flush=True)
json.dump(res, open(OUT + '/c1_controls.json', 'w'), indent=1)
