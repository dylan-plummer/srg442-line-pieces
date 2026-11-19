"""q = 23, step 2: global line-piece system on the base (CP-SAT).
Base: Paley(p) on GF(p) (p = r^2), elements indexed r*x + y; circulant tournament i -> i + d (d = 1..(p-1)/2) on the
indices, exactly as in w1lib.build_template.  Unknowns: a set of arcs forming a matching (the non-constant fibres),
for each chosen arc a class (M_tail, M_head) of per-arc solutions for its T, and a sign c_i = +-1 for every unmatched
fibre (M_i = q c_i).  Constant-part equation (exact):  1 + sum_j eps_ij M_j = (rho / q) M_i  for all i,
eps_ij = +1 on base edges, -1 on non-edges.  (The point at infinity has value +1; its equation sum_i q M_i = rho
follows by summing, since eps 1 = 0.)
usage: c4_global.py q ARCSJSON rho OUTJSON MODE [timelimit]   MODE = all (enumerate every pattern) | one"""
import sys, os, json, collections, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from ortools.sat.python import cp_model
from w1lib import GF, check_conference_graph

q = int(sys.argv[1]); arcsj = json.load(open(sys.argv[2])); rho = int(sys.argv[3]); outf = sys.argv[4]; mode = sys.argv[5]
tl = float(sys.argv[6]) if len(sys.argv) > 6 else 600.0
p = q + 2; r = int(round(p ** 0.5)); assert r * r == p and abs(rho) == q * r
Fp = GF(p); B = Fp.paley_adj(); assert check_conference_graph(B)
eps = np.where(B == 1, 1, -1); np.fill_diagonal(eps, 0)
arcs = []
for i in range(p):
    for d in range(1, (p - 1) // 2 + 1):
        j = (i + d) % p; arcs.append((i, j, 'Q' if B[i, j] else '0N'))
cls = {Tn: collections.Counter((sum(g), sum(h)) for g, h in arcsj[str(rho)][Tn]) for Tn in ('Q', '0N')}
m = cp_model.CpModel()
cp = [m.NewBoolVar(f'cp{i}') for i in range(p)]; cn = [m.NewBoolVar(f'cn{i}') for i in range(p)]
y = {}
for e, (i, j, Tn) in enumerate(arcs):
    for c in cls[Tn]: y[(e, c)] = m.NewBoolVar(f'y{e}_{c}')
inc = collections.defaultdict(list); Mterms = collections.defaultdict(list)
for (e, c), v in y.items():
    i, j, Tn = arcs[e]; inc[i].append(v); inc[j].append(v)
    Mterms[i].append((c[0], v)); Mterms[j].append((c[1], v))
M = []
for i in range(p):
    m.Add(cp[i] + cn[i] + sum(inc[i]) == 1)
    M.append(q * cp[i] - q * cn[i] + sum(a * v for a, v in Mterms[i]))
k = rho // q
for i in range(p):
    m.Add(1 + sum(int(eps[i, j]) * M[j] for j in range(p) if j != i) == k * M[i])

def pattern(sv):
    chosen = [(arcs[e][:2], list(c), cls[arcs[e][2]][c]) for (e, c), v in y.items() if sv.Value(v)]
    Mv = [int(sv.Value(x)) if not isinstance(x, int) else x for x in M]
    cons = {i: (1 if sv.Value(cp[i]) else -1) for i in range(p) if sv.Value(cp[i]) or sv.Value(cn[i])}
    mult = 1
    for _, _, mm in chosen: mult *= mm
    return {'arcs': [[list(a), c] for a, c, _ in chosen], 'const': cons, 'M': Mv, 'labelled': mult}

class CB(cp_model.CpSolverSolutionCallback):
    def __init__(s): super().__init__(); s.pats = []
    def on_solution_callback(s):
        s.pats.append(pattern(s))

sv = cp_model.CpSolver(); sv.parameters.max_time_in_seconds = tl
t0 = time.time()
if mode == 'all':
    sv.parameters.enumerate_all_solutions = True; sv.parameters.num_workers = 1
    cb = CB(); st = sv.Solve(m, cb); pats = cb.pats
else:
    sv.parameters.num_workers = 8; st = sv.Solve(m); pats = [pattern(sv)] if st in (cp_model.OPTIMAL, cp_model.FEASIBLE) else []
status = sv.StatusName(st)
out = {'q': q, 'p': p, 'rho': rho, 'mode': mode, 'status': status, 'time': round(time.time() - t0, 1),
       'per_arc_classes': {Tn: {str(c): n for c, n in cls[Tn].items()} for Tn in cls}, 'patterns': len(pats),
       'labelled_total': sum(x['labelled'] for x in pats), 'matching_sizes': dict(collections.Counter(len(x['arcs']) for x in pats)),
       'solutions': pats[:50]}
print(json.dumps({k2: v for k2, v in out.items() if k2 != 'solutions'}), flush=True)
json.dump(out, open(outf, 'w'))
