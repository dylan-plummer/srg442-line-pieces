"""q = 23: build the regular two-graph on p q^2 + 1 = 13226 points (q = 23, p = 25) in two ways and verify exactly.
  mathon : literal Mathon Theorem 4.1 (base Paley(25) on GF(25) = GF(5)[w], w^2 = 2, index 5x+y; circulant
           skew-Latin square (4.1)),
  template: our d = 1 Paley-fiber template (Pal(529) fibers, same base, circulant tournament on the indices).
Checks per build: symmetric 0/1, zero diagonal; S^2 = 13225 I on ALL 13226^2 entries (float32 BLAS, exact);
srg(13225, 6612, 3305, 3306) of the descendant at infinity via float32 A^2 (exact); writes the 13225-vertex
descendant as graph6 with SHA-256.
usage: c5_build23.py mathon|template"""
import sys, os, json, time, hashlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from w1lib import *

which = sys.argv[1]; q, p = 23, 25
OUT = 'data/q23'; os.makedirs(OUT, exist_ok=True)
t0 = time.time(); log = lambda s: print(f'[{time.time() - t0:7.1f}s] {s}', flush=True)
F25 = GF(25); B = F25.paley_adj(); assert check_conference_graph(B); log('Paley(25) base ok')
if which == 'mathon':
    A, mb = build_mathon(q, B); info = {'g': mb['g'], 't': mb['t']}
else:
    A, st, arcs = build_template(q, B); info = {'pos': st['pos'], 'neg': st['neg']}
log(f'built {which}: {A.shape}, dtype {A.dtype}, {A.nbytes / 1e6:.0f} MB')
res = {'which': which, 'q': q, 'p': p, 'info': {k: str(v) for k, v in info.items()}}
res['symmetric'] = bool((A == A.T).all()); res['zero_diag'] = not bool(A.diagonal().any())
res['is01'] = bool(A.max() <= 1)
log(f"symmetric {res['symmetric']} zero_diag {res['zero_diag']}")
res['S2'] = seidel_square_check(A, B=1024); log(f"S^2 check: {res['S2']}")
res['srg_descendant'] = srg_check(A, B=1024); log(f"descendant: {res['srg_descendant']}")
g6 = to_g6(A); fn = f'{OUT}/graphs/twograph13226_{which}_descendant_srg13225.g6'
open(fn, 'w').write(g6 + '\n')
res['g6'] = fn; res['g6_sha256'] = hashlib.sha256((g6 + '\n').encode()).hexdigest(); res['g6_bytes'] = len(g6) + 1
A2 = from_g6(g6); res['g6_roundtrip'] = bool((A2 == A).all()); del A2
log(f"g6 written {fn} sha256 {res['g6_sha256']} roundtrip {res['g6_roundtrip']}")
json.dump(res, open(f'{OUT}/c5_build23_{which}.json', 'w'), indent=1)
