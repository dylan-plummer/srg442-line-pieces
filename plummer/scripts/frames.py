"""Frames of the line-piece members (Kovaľ's census expanded to labeled switching sets).
python -I frames.py CENSUSDIR OUT.tsv

For each two-graph (mathon = Omega_M, ours = Omega_K) and degree, every census class is expanded under Aut(Omega) to its
labeled switching sets x (x = N(oo) of the member of that degree, oo not switched). The frame of x records per fiber:
W / E (whole / empty) or "t{clq|coc}m": x is the union of t parallel lines of the form b - m a (m = inf: the form a),
which are cliques (net direction) or cocliques of the fiber's graph. Checked and written per frame: number of labeled
sets, the decoupled pairs and the number of choices on each, the members fixed by the fiber-fixing C3 (centered
members), the orbit of the frame under Aut(Omega), and the census classes meeting it by |Aut| and orbit type."""
import sys, os, collections
import numpy as np, igraph as ig
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from census_props import g6dec

LV = {m: np.array([(b - m * a) % 7 for a in range(7) for b in range(7)]) for m in range(7)}
LV["inf"] = np.array([a for a in range(7) for b in range(7)])
I = np.eye(442, dtype=np.int64)


def fibre_type(xf, C):
    s = int(xf.sum())
    if s in (0, 49): return "W" if s else "E"
    for m, lv in LV.items():
        if all(len(set(xf[lv == c])) == 1 for c in range(7)):
            P = (lv[:, None] == lv[None, :]) & ~np.eye(49, dtype=bool)
            net = "clq" if C[P].all() else ("coc" if not C[P].any() else "mix")
            return f"{len({c for c in range(7) if xf[lv == c][0]})}{net}{m}"
    raise ValueError("not a line piece")


def group(M):
    gens = [np.array(p) for p in ig.Graph.Adjacency(M.tolist(), mode="undirected").automorphism_group()]
    G = {tuple(range(441)): np.arange(441)}
    while True:
        new = {tuple(q[p]): q[p] for q in G.values() for p in gens}
        if new.keys() <= G.keys(): return list(G.values())
        G.update(new)


if __name__ == "__main__":
    K, out = sys.argv[1], sys.argv[2]
    props = {}
    for l in open(os.path.join(os.path.dirname(out), "census_props.tsv")).read().splitlines()[1:]:
        r = l.split("\t"); props[(r[0], int(r[1]))] = (int(r[3]), r[4], r[6])
    fo = open(out, "w")
    fo.write("twograph\tdegree\tframe\torbit\tsets\tpairs\tchoices\tcentred\tcentred_x\tclasses\n")
    for tg in ("mathon", "ours"):
        M = None
        for deg, fn in ((231, f"koval442_{tg}_deg231_575.g6"), (210, f"koval442_{tg}_deg210_complemented_603.g6")):
            reps = []
            for ln, line in enumerate(open(os.path.join(K, fn)).read().split()):
                A = g6dec(line)
                if A.sum(1)[0] != deg: A = 1 - I - A                       # the member itself, of degree deg
                x = A[441, :441].astype(np.int64)
                Mp = A[:441, :441] ^ (x[:, None] != x[None, :]).astype(np.int64)
                if M is None:
                    M = Mp; G = group(M); S = 1 - np.eye(441, dtype=np.int64) - 2 * M
                    fix = [q for q in G if all(q[49 * i] // 49 == i for i in range(9))]
                    assert len(G) == 9 and len(fix) == 3
                assert np.array_equal(Mp, M)
                reps.append((x, (fn, ln)))
            ftype = lambda x: tuple(fibre_type(x[49 * i:49 * i + 49], M[49 * i:49 * i + 49, 49 * i:49 * i + 49]) for i in range(9))
            fr, cls = collections.defaultdict(dict), collections.defaultdict(set)
            for x, key in reps:
                for q in G:
                    y = np.zeros(441, np.int64); y[q] = x; f = ftype(y)
                    fr[f][y.tobytes()] = key; cls[f].add(key)
            assert sum(map(len, fr.values())) == 5145
            orbit = {}
            for f in sorted(fr):                                         # orbits of frames under Aut(Omega)
                if f in orbit: continue
                x = np.frombuffer(next(iter(fr[f])), dtype=np.int64); o = len(set(orbit.values()))
                for q in G:
                    y = np.zeros(441, np.int64); y[q] = x; orbit.setdefault(ftype(y), o)
            for f, X in sorted(fr.items(), key=lambda kv: (orbit[kv[0]], kv[0])):
                L = [np.frombuffer(b, dtype=np.int64) for b in X]; x0 = 1 - 2 * L[0]
                nc = [i for i in range(9) if f[i] not in "WE"]
                pairs = [(i, j) for i in nc for j in nc if i < j and np.ptp(
                    S[49 * i:49 * i + 49, 49 * j:49 * j + 49] @ (x0[49 * j:49 * j + 49] - x0[49 * j:49 * j + 49].mean())) > 1e-9]
                ch = [len({x[49 * i:49 * i + 49].tobytes() + x[49 * j:49 * j + 49].tobytes() for x in L}) for i, j in pairs]
                assert int(np.prod(ch)) == len(L)                            # the members are all combinations
                cen = [x for x in L if all((x[q] == x).all() for q in fix)]
                fixed = lambda k: 442 if props[k][0] == 1 else props[k][1].count("1")    # points fixed by Aut
                cc = collections.Counter(f"aut{props[k][0]}fix{fixed(k)}" + ("*" if props[k][2] else "") for k in cls[f])
                fo.write(f"{tg}\t{deg}\t{','.join(f)}\t{orbit[f]}\t{len(L)}\t{' '.join(f'{i}{j}' for i, j in pairs)}\t"
                         f"{'x'.join(map(str, ch))}\t{len(cen)}\t{np.packbits(cen[0].astype(np.uint8)).tobytes().hex() if cen else ''}\t"
                         f"{' '.join(f'{k}={v}' for k, v in sorted(cc.items()))}\n")
            print(tg, deg, len(fr), "frames in", len(set(orbit.values())), "orbits")
    fo.close()
