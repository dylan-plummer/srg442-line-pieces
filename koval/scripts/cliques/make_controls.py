"""Controls: block graphs of known Steiner systems (positive), written as graph6.
Expected Delsarte-clique counts are stated with the reason; the counter must reproduce them."""
import itertools, random, sys, os
OUT = sys.argv[1]

def g6(n, edges):
    adj = [[0]*n for _ in range(n)]
    for a, b in edges: adj[a][b] = adj[b][a] = 1
    bits = [adj[i][j] for j in range(1, n) for i in range(j)]
    while len(bits) % 6: bits.append(0)
    if n < 63: head = chr(n+63)
    else: head = '~' + chr(((n >> 12) & 63)+63) + chr(((n >> 6) & 63)+63) + chr((n & 63)+63)
    body = ''.join(chr(63+int(''.join(map(str, bits[i:i+6])), 2)) for i in range(0, len(bits), 6))
    return head + body

def block_graph(blocks):
    blocks = [frozenset(b) for b in blocks]
    E = [(i, j) for i, j in itertools.combinations(range(len(blocks)), 2) if blocks[i] & blocks[j]]
    return g6(len(blocks), E)

def check_design(pts, blocks, lam=1):
    c = {}
    for b in blocks:
        for p in itertools.combinations(sorted(b), 2): c[p] = c.get(p, 0)+1
    assert len(c) == len(pts)*(len(pts)-1)//2 and set(c.values()) == {lam}, "not a 2-design"

def proj_lines(q_dim, q):   # lines of PG(d,q), q prime
    vecs = [v for v in itertools.product(range(q), repeat=q_dim+1) if any(v)]
    def norm(v):
        i = next(k for k in range(len(v)) if v[k]); inv = pow(v[i], q-2, q)
        return tuple(x*inv % q for x in v)
    pts = sorted(set(map(norm, vecs))); idx = {p: i for i, p in enumerate(pts)}
    lines = set()
    for a, b in itertools.combinations(pts, 2):
        L = frozenset(idx[norm(tuple((s*x+t*y) % q for x, y in zip(a, b)))] for s in range(q) for t in range(q) if s or t)
        lines.add(L)
    return list(range(len(pts))), [sorted(L) for L in lines]

def aff_lines(dim, q):     # lines of AG(dim,q), q prime
    pts = list(itertools.product(range(q), repeat=dim)); idx = {p: i for i, p in enumerate(pts)}
    lines = set()
    for a, b in itertools.combinations(pts, 2):
        d = tuple((y-x) % q for x, y in zip(a, b))
        lines.add(frozenset(idx[tuple((x+t*dd) % q for x, dd in zip(a, d))] for t in range(q)))
    return list(range(len(pts))), [sorted(L) for L in lines]

def sts_stinson(v, seed):  # Stinson hill-climbing for a random STS(v)
    rnd = random.Random(seed)
    live = {x: set(range(v)) - {x} for x in range(v)}
    other = {}; blocks = set(); nb = v*(v-1)//6
    while len(blocks) < nb:
        x = rnd.choice([x for x in range(v) if live[x]])
        y, z = rnd.sample(sorted(live[x]), 2)
        if z not in live[y]:
            w = other[(y, z)]; B = frozenset((y, z, w)); blocks.discard(B)
            for a, b in itertools.permutations(B, 2):
                live[a].add(b); other.pop((a, b), None)
            for a, b in ((y, z),):
                pass
        B = frozenset((x, y, z)); blocks.add(B)
        for a, b in itertools.permutations(B, 2):
            live[a].discard(b); other[(a, b)] = next(iter(B - {a, b}))
    return list(range(v)), [sorted(b) for b in blocks]

ctl = []
pts, L = proj_lines(3, 2); check_design(pts, L)
ctl.append(("PG32_lines_S2_3_15", L, 15+15, "srg(35,18,9,9); 7-cliques = 15 point pencils + 15 planes (Grassmann duality)", 7, 3))
pts, L = aff_lines(3, 3); check_design(pts, L)
ctl.append(("AG33_lines_S2_3_27", L, 27, "srg(117,36,15,9); 13-cliques = 27 point pencils (a plane's 12 lines contain parallels)", 13, 3))
pts, L = proj_lines(4, 2); check_design(pts, L)
ctl.append(("PG42_lines_S2_3_31", L, 31, "srg(155,42,17,9); 15-cliques = 31 point pencils (planes have 7 lines, solids 35 non-pairwise-meeting)", 15, 3))
pts, L = proj_lines(3, 3); check_design(pts, L)
ctl.append(("PG33_lines_S2_4_40", L, 40+40, "srg(130,48,20,16); 13-cliques = 40 pencils + 40 planes", 13, 4))
for s in range(3):
    pts, B = sts_stinson(19, s); check_design(pts, B)
    ctl.append((f"STS19_random_seed{s}", B, None, "srg(57,24,11,9); >= 19 pencils of size 9 (extra 9-cliques would need 9 pairwise-meeting triples outside a pencil)", 9, 3))
for s in range(2):
    pts, B = sts_stinson(15, s); check_design(pts, B)
    ctl.append((f"STS15_random_seed{s}", B, None, "srg(35,18,9,9); >= 15 pencils", 7, 3))
os.makedirs(OUT, exist_ok=True)
with open(os.path.join(OUT, "controls.tsv"), "w") as f:
    for name, B, exp, note, W, ksz in ctl:
        g = block_graph(B)
        open(os.path.join(OUT, name + ".g6"), "w").write(g + "\n")
        with open(os.path.join(OUT, name + ".blocks"), "w") as fb:
            for b in B: fb.write(" ".join(map(str, b)) + "\n")
        f.write(f"{name}\t{len(B)}\t{W}\t{exp}\t{note}\n")
        print(name, len(B), W, exp, note)
