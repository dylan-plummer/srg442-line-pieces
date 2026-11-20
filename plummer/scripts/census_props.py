"""Properties of the 2356 line-piece members (Kovaľ's census files) and where our 29 graphs lie among them.
python -I census_props.py CENSUSDIR OUT.tsv OURS.g6 [OURS.g6 ...]
Per graph: source file, line, exact srg check, |Aut| (bliss via igraph), automorphism orbit lengths, canonical hash
(igraph canonical_permutation), and the first of our files that contains an isomorphic graph."""
import sys, os, glob, hashlib
from multiprocessing import Pool
import numpy as np, igraph as ig


def g6dec(s):
    s = s.strip(); d = [ord(c) - 63 for c in s]
    if d[0] == 63: n = (d[1] << 12) | (d[2] << 6) | d[3]; d = d[4:]
    else: n = d[0]; d = d[1:]
    bits = [(x >> (5 - i)) & 1 for x in d for i in range(6)]
    A = np.zeros((n, n), dtype=np.int64); t = 0
    for j in range(1, n):
        for i in range(j): A[i, j] = A[j, i] = bits[t]; t += 1
    return A


def work(line):
    A = g6dec(line); n = len(A); J = np.ones((n, n), dtype=np.int64); I = np.eye(n, dtype=np.int64)
    if A.sum(1)[0] != 231: A = J - I - A                         # our files hold some members as degree 210
    ok = bool((A == A.T).all() and A.trace() == 0 and np.array_equal(A @ A, 231 * I + 120 * A + 121 * (J - I - A)))
    g = ig.Graph.Adjacency(A.tolist(), mode="undirected")
    c = np.array(g.permute_vertices(g.canonical_permutation()).get_adjacency().data, dtype=np.uint8)
    gens = g.automorphism_group()
    orb = sorted(len(o) for o in ig.Graph(n=n, edges=[(i, int(p[i])) for p in gens for i in range(n)]).connected_components()) \
        if gens else [1] * n
    return ok, int(g.count_automorphisms()), orb, hashlib.sha1(np.packbits(c).tobytes()).hexdigest()


if __name__ == "__main__":
    K, out, ours = sys.argv[1], sys.argv[2], sys.argv[3:]
    files = sorted(glob.glob(os.path.join(K, "koval442_*_deg*.g6")))
    cnames = set(map(os.path.basename, files))
    jobs = [(os.path.basename(f), i, l) for f in files for i, l in enumerate(open(f).read().split())]
    jobs += [(os.path.basename(f), i, l) for f in ours for i, l in enumerate(open(f).read().split())]
    with Pool(int(os.environ.get("NPROC", 4))) as P: res = P.map(work, [l for *_, l in jobs], chunksize=4)
    mine = {}
    for (src, i, _), (ok, na, orb, h) in zip(jobs, res):
        if src in map(os.path.basename, ours): mine.setdefault(h, f"{src}:{i}")
    with open(out, "w") as fo:
        fo.write("source\tline\tsrg\taut\torbits\thash\tours\n")
        for (src, i, _), (ok, na, orb, h) in zip(jobs, res):
            if src in cnames:                                       # koval442_new15.g6 is one of ours
                fo.write(f"{src}\t{i}\t{int(ok)}\t{na}\t{orb if len(set(orb)) > 1 or na > 1 else '1^442'}\t{h}\t{mine.get(h, '')}\n")
    rows = [l.split("\t") for l in open(out)][1:]
    print(len(rows), "census graphs,", sum(r[2] == "0" for r in rows), "failing srg,", len({r[5] for r in rows}), "distinct;",
          sum(bool(r[6].strip()) for r in rows), "of ours matched,", len(mine), "of ours distinct")
