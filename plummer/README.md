# D.P.'s programs and data for the joint 442 paper

All programs are run with `python3 -I`. `reproduce.sh WORKDIR [q23]` re-runs everything below except the two long q = 23
searches, and compares the outputs with the files in `data/`.

## scripts/

| program | backs (paper) | input → output |
|---|---|---|
| `check_paper442.py` | Tables orders, arcs (Ω_M, Ω_P data), aut; facts (F1)–(F3) in the proof of the decoupling lemma; part (3) on 48 members | imports A.K.'s `koval/scripts/census/P2_family442.py` (unmodified) |
| `census_props.py` | Theorem (census): exact srg check, \|Aut\| (bliss via igraph), orbit lengths, canonical hash of all 2356 graphs; our 29 graphs matched | unzipped census + `data/graphs/*.g6` → `data/census_props.tsv` |
| `frames.py` | frames: 5145 = 15 · 7³ labeled members per degree and two-graph, one centered member per frame | census + `census_props.tsv` → `data/frames.tsv` |
| `make_figures.py member \| aut` | Figure 2 (member, `data/member_x.npy`, Sx = −21x checked) and Table aut | → `paper/figures/` |
| `sporadic_check.py S7.npy [S23.npy]` | Lemma (slope decomposition) at q = 7, 23; pair equations on the 14 graphs of [Plummer442] (504 pairs) | |
| `sporadic_figs.py S7.npy` | the two figures of the sporadicity section (every drawn number recomputed) | → `paper/figures/` |
| `check_construction.py` | Mathon's srg(441,220,109,110) built from his Theorem 4.1 (library for `make_figures.py`) | |

## data/

- `graphs/srg442_231_120_121_all.g6`: the 14 graphs of [Plummer442] (members of Ω_M), `graphs/srg442_210_99_100_all.g6`:
  their complements; `graphs/koval442_new15.g6`: 15 members of Ω_P found earlier with the same line-piece program
  (14 centered, 1 rigid; in the census they match 6 graphs of degree 231 and 9 of degree 210). The 14 + 15 = 29 graphs
  are the bracketed counts in Table aut.
- `census_props.tsv` (+ `.out`, the run summary), `frames.tsv`, `member_x.npy` (the member drawn in Figure 2).

## m23/ — h-invariant regular members of Mathon's two-graph (Proposition q23h, Remark "controls", Table "search", Figure "reduction")

| step | program | |
|---|---|---|
| 1 | `build.py q s S.npy` | Mathon's two-graph Ω_q (uses `mathon.py`, `srgutil.py` from the srg442 repository), checks S² = (n−1)I |
| 2 | `pairs.py S.npy q K OUT` | the exact pair equations of h-invariant vectors (fitted and verified on random vectors) |
| 3 | `structure.py q PAIRS OUT.npz` | closed forms of the h-invariant line sums (conference matrix σ), slot types |
| 4 | `mitm.py q STRUCT LOC.npy` | all local patterns per fiber (q = 7: 6256; q = 23: 206,632 for K = 6670, 189,142 for K = 6555) |
| 5 | `glob.py S q K STRUCT LOC OUT` | exact depth-first search with forward checking; `TSYM=1 ROOTS=a:b` slices, `TGRP=1` for ⟨h,t⟩ |
| 6 | `verify_hit.py` | rebuilds each hit and checks Sx = ρx exactly |
| – | `twoval.py q` | all T ⊆ F_q with χ′ ∗ 1_T two-valued (Remark two-valued; q = 7, 11, 19, 23, 31) |
| – | `reduced.py` | the reduced line-piece model (control: reproduces the 15 + 15 at q = 7) |
| – | `pieces.py S.npy q K TLIM WORKERS MAXSOL OUT` | CP-SAT on the line-piece model of [Plummer442] (fibers empty, whole, or ℓ_m⁻¹(T), T a union of {0}, Q, N); q = 23: INFEASIBLE for K = 6670 and 6555 (`results/p23_6670.log`, `results/p23_6555.log`, about 22 min each) |

`results/`: `glob7.txt` (q = 7: exactly the 15 + 15 centered members), `glob23*.txt` (q = 23: empty),
`h_slices_6670/`, `h_slices_6555/` (the 20 slices of each complete q = 23 h-invariant search: 284,271,840 and
276,655,619 nodes, 0 solutions), `twoval23.txt`, `twoval31.txt`, `red7.txt`, `red23.txt`, and `SUMMARY.txt` (the
working log of 8 November 2026, kept as a record; superseded where it says "running").

To redo a q = 23 search (about 47 CPU hours for K = 6670, 45 for K = 6555; 20 slices in parallel):
```
python3 -I m23/build.py 23 5 S23.npy
python3 -I m23/pairs.py S23.npy 23 6670 pairs23.npy && python3 -I m23/structure.py 23 pairs23.npy struct23.npz
python3 -I m23/mitm.py 23 struct23.npz loc23.npy
for i in $(seq 0 19); do a=$((i*206632/20)); b=$(((i+1)*206632/20))
  TSYM=1 ROOTS=$a:$b python3 -I m23/glob.py S23.npy 23 6670 struct23.npz loc23.npy glob23_h.txt > s$i.log & done; wait
```
(for K = 6555 replace 6670 and 206632 by 6555 and 189142).
