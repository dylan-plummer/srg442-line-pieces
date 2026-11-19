"""Check of the census numbers on our side against the files of this folder.  Standard library only.

What it checks (and nothing else):
  * the census archive data/koval442_linepiece_2356.zip: it holds exactly the four graph6 files and README.txt; each
    graph6 file has the SHA-256, size and number of lines listed in README.md (575, 603, 575, 603 graphs);
  * the clique table data/cliques_per_graph.tsv (made with scripts/cliques/dcliq.c): 2356 rows, one per census graph;
    0 cliques of size 22 in both passes, srg(442,231,120,121) for every graph, clique number 13 for 2192 graphs and
    14 for 164, all 164 in the file of Omega_M of degree 231;
  * the orbit-count arithmetic used in the text (5145 vectors, 575 and 603 orbits, Table 3 totals, 52 descendant
    classes, the Delsarte bound 22).
It does not check the graphs themselves (see P3_members.py, members.py, gcanon.c) or any number of the co-author's side.
Run from koval/:  python3 -I scripts/checks/check_numbers.py      (prints ALL OK at the end)
"""
import hashlib, zipfile
from collections import Counter

ZIP = 'data/koval442_linepiece_2356.zip'
FILES = [('koval442_ours_deg231_575.g6', 575, 9343175, '46fe141c53a83fffb633659c5b7a6c0c0d68d4c1e980324a15a7d483c68ebed4'),
         ('koval442_ours_deg210_complemented_603.g6', 603, 9798147, '8f91dba6aa6ac3930a3ffb8f334be4e59983fa08840888c1d71de7bda3829751'),
         ('koval442_mathon_deg231_575.g6', 575, 9343175, 'ab384ab8bf074ddf4b819db90e1fc39e9fc0fdab4751f301186c0473d232f163'),
         ('koval442_mathon_deg210_complemented_603.g6', 603, 9798147, '7e1207a4f8b8735176954fc2726a74cef7da5c438dfd6e9686d6e933d58bfde2')]
ok = True
print('archive', ZIP, 'sha256', hashlib.sha256(open(ZIP, 'rb').read()).hexdigest())
with zipfile.ZipFile(ZIP) as z:
    names = sorted(z.namelist())
    r = names == sorted([f[0] for f in FILES] + ['README.txt'])
    print('archive members:', names, 'as listed:', r); ok &= r
    for fn, cnt, size, sha in FILES:
        b = z.read(fn)
        lines = [l for l in b.split(b'\n') if l.strip()]
        r = (len(lines) == cnt, len(b) == size, hashlib.sha256(b).hexdigest() == sha)
        print(fn, 'lines/size/sha256 as listed:', r); ok &= all(r)

rows = [l.split('\t') for l in open('data/cliques_per_graph.tsv').read().split('\n')[1:] if l.strip()]
print('clique table rows', len(rows)); ok &= len(rows) == 2356
om = Counter(); om_file = Counter(); n14 = Counter(); n13 = {}
for r in rows:
    f, line, n22a, n22b, w, c13, c14, srg = r[:8]
    ok &= n22a == '0' and n22b == '0' and srg == 'srg(442,231,120,121)'
    om[w] += 1; om_file[(f, w)] += 1
    if int(c14) > 0: n14[int(c14)] += 1
    n13.setdefault(f, []).append(int(c13))
print('22-cliques: all zero in both passes, all srg:', ok)
print('clique numbers:', dict(om), 'by file:', dict(om_file))
r = (dict(om) == {'13': 2192, '14': 164}
     and {f for (f, w) in om_file if w == '14'} == {'koval442_mathon_deg231_575'})
print('clique number 13 for 2192 graphs, 14 for 164, all 164 in koval442_mathon_deg231_575:', r); ok &= r
r = Counter(f for f, *_ in rows) == Counter({f[0][:-3]: f[1] for f in FILES})
print('rows per census file 575/603/575/603:', r); ok &= r
print('graphs with 14-cliques, by number of 14-cliques:', dict(n14), 'total', sum(n14.values()))
print('13-clique counts min-max by file:', {f: (min(v), max(v)) for f, v in n13.items()})
# arithmetic used in the text
assert (5145 + 30) // 9 == 575 and (5145 + 30) % 9 == 0
assert (5145 + 30 + 6 * 42) // 9 == 603 and (5145 + 30 + 6 * 42) % 9 == 0
assert 2 * (575 + 603) == 2356 and 2356 - 29 == 2327 and 15 * 7 ** 3 == 5145
assert 570 * 9 + 5 * 3 == 5145 and 558 * 9 + 3 * 3 + 36 * 3 + 6 == 5145
assert 2 * (570 + 558) == 2256 and 2 * (5 + 3) == 16 and 2 * 36 == 72 and 2 * 6 == 12 and 2256 + 16 + 72 + 12 == 2356
assert (442 - 10) // 3 + 10 == 154 and (442 - 1) // 3 + 1 == 148 and 1 + 3 + 48 == 52
assert 1 + 231 // 11 == 22 and 22 * 210 // 420 == 11
print('arithmetic identities: OK')
print('ALL OK' if ok else 'MISMATCH')
