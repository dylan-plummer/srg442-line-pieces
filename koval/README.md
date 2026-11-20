# The census of 2356 strongly regular graphs on 442 points: data and programs (A. Kovaľ)

This folder holds the data and the programs behind the statements of the paper that rest on the second author's
side: the census of the line-piece regular members of two regular two-graphs on 442 points, their pairwise
non-isomorphism, the clique numbers, the comparison of the template two-graphs, and the computations at q = 23.
The first author's programs are in `../plummer/`.

**Names.** Ω_M is Mathon's conference two-graph for q = 7. Ω_P is the variant whose nine fibers are copies of
Paley(49) (Ω_P is called Ω_K in the draft of 14 Nov). In file names and program identifiers Ω_M is `mathon` and Ω_P is
`ours`. Some program comments still write `Omega_K` for Ω_P. "Base graph" means the descendant of a two-graph at its
point ∞. It is an srg(441,220,109,110).

## Contents

```
data/koval442_linepiece_2356.zip      the census: four graph6 files + README.txt
data/cliques_per_graph.tsv            clique table, one row per census graph (2356 rows)
data/template_data/samples/s1_*.json  template data of 16 further two-graphs on 442 points
data/template_data/gauge/n_*.json     template data of 5 further two-graphs (one arc of Ω_P changed)
data/template_data/gauge/g_*.json     5 gauge controls (isomorphic to Ω_P by construction)
scripts/census/                       template builder, census enumeration, canonical labeling, two-graph invariants
scripts/cliques/                      clique counter and its controls
scripts/q23/                          Mathon's blocks and the template at q = 3, 7, 23; the two-graphs on 13226 points
scripts/q23/results/                  saved outputs of the q = 23 per-arc and global computations
scripts/checks/                       check_numbers.py, check_template_data.py
```

## The census files

The archive `data/koval442_linepiece_2356.zip` (SHA-256
`a830a39ca8465956c6c20474389a025e48a239c667c4c123c7a0d158e117050f`, 6,848,534 bytes) contains a short README.txt
and four graph6 files, one graph per line. The hash of the archive depends on how it was packed. The hashes of the
four graph6 files do not, so cite those.

| file in the zip | graphs | two-graph, degree | SHA-256 of the graph6 file |
|---|---|---|---|
| `koval442_ours_deg231_575.g6` | 575 | Ω_P, 231 | `46fe141c53a83fffb633659c5b7a6c0c0d68d4c1e980324a15a7d483c68ebed4` |
| `koval442_ours_deg210_complemented_603.g6` | 603 | Ω_P, 210, complemented | `8f91dba6aa6ac3930a3ffb8f334be4e59983fa08840888c1d71de7bda3829751` |
| `koval442_mathon_deg231_575.g6` | 575 | Ω_M, 231 | `ab384ab8bf074ddf4b819db90e1fc39e9fc0fdab4751f301186c0473d232f163` |
| `koval442_mathon_deg210_complemented_603.g6` | 603 | Ω_M, 210, complemented | `7e1207a4f8b8735176954fc2726a74cef7da5c438dfd6e9686d6e933d58bfde2` |

All 2356 graphs are srg(442,231,120,121). A member of degree 210 is stored as its complement, which lies in the
complementary two-graph. Each file holds one graph per orbit of Aut(Ω) = C3 × C3.

## Requirements

- Python 3 with numpy. We used Python 3.12, numpy 1.26 and pynauty 2.8.
- pynauty for every program that imports `scripts/census/P2_family442.py`, and for `c1_controls.py`.
- OR-Tools (CP-SAT) for `scripts/q23/c4_global.py` and for `P3_members.py solve`.
- gcc for the C programs: `gcc -O2 -o bin/NAME scripts/.../NAME.c`.
- SageMath only for the optional SageMath control graph of `c1_controls.py`.

Run everything from this folder (`koval/`). The programs read and write relative paths under `data/` and `bin/`.

Two input graphs are not stored here:
- `data/srg_441_220_109_110_d1_z9.g6` is the base graph of Ω_P, the graph of [K833, Computation 6.7]. Write it with
  `python3 -I scripts/census/P2_family442.py base ours data/srg_441_220_109_110_d1_z9.g6`. The file has SHA-256
  `0356d1e21c3ff601157935dc52d380384314fe196d10af7585d67a1dd7fb6c80`. Use `base mathon` for the base graph of Ω_M.
- Two different graphs, used only as controls, have distinct names:
  - `data/mathon441_plummer442.g6` is Mathon's instance as given by the first author: the file `data/mathon441.g6`
    of [Plummer442], github.com/dylan-plummer/srg442. `P2_family442.py controls` uses it.
  - `data/mathon441_sage_t2.g6` is SageMath's `graphs.MathonPseudocyclicStronglyRegularGraph(t=2)`, written with
    `graph6_string()`. `scripts/q23/c1_controls.py` uses it. It gives a third regular two-graph on 442 points,
    isomorphic to neither Ω_M nor Ω_P.

## How each statement is reproduced

### 1. The template two-graphs Ω_M and Ω_P (Definition 2.1, Table 2, Lemma 2.4)

- `scripts/census/P2_family442.py` is the builder used by all census programs. `mathon_choice()` builds Ω_M and
  `ours_choice()` builds Ω_P. `template_ok` checks the template conditions, and `verify` checks
  A² = kI + λA + μ(J − I − A) exactly.
- `python3 -I scripts/checks/check_template_data.py` checks Ω_M, Ω_P and the 26 files of `data/template_data/`.
  Expected for every one: `template_ok (True, 'ok')`, `srg(441,220,109,110)`, `S^2 = 441 I: True`. The base graph of
  Ω_P must have the SHA-256 above. The last line is `ALL OK`.
- `python3 -I scripts/census/P2_family442.py controls` needs the two files above. It prints
  `equal to file as LABELED graphs: True` for Ω_M (against `mathon441_plummer442.g6`) and for Ω_P (against
  `srg_441_220_109_110_d1_z9.g6`). It also prints |Aut(G)| = 9 with 51 orbits for each.
- `python3 -I scripts/q23/mathon_literal_blocks.py` checks that Mathon's literal blocks at q = 3, 7 and 23 are
  template data. It includes a negative control: A_d itself is not a line incidence. The last line is
  `ALL CHECKS PASSED`.
- `python3 -I scripts/q23/c1_controls.py` reproduces Mathon's printed matrices (4.15) and (4.16) at q = 3. It uses the
  transcription of Mathon's blocks in `scripts/q23/w1lib.py`, the same one `mathon_literal_blocks.py` uses. Expected
  output:
  - `q=3 (4.15) reproduced: {'Ad': True, 'A1': True, 'A2': True}`
  - `q=3 A_I (4.16) equal: True`, with srg(45,22,10,11) and S² = 45 I
  - at q = 7, S² = 441 I for both builds
  - `template q=7 g6 == archived srg_441_220_109_110_d1_z9.g6: True`
  - for the SageMath graph: srg and S² = 441 I, and `entrywise equal ...: False`

  The optional part `W1_NAUTY=1` compares two-graphs by pynauty certificates of the 884-vertex double covers. It did
  not finish in reasonable time in our runs. The statement that the SageMath graph is a third two-graph was obtained
  by comparing pynauty certificates of descendants, with a program that is not in this folder.

### 2. The census: 5145 line-piece vectors, 575 and 603 classes per two-graph (Theorem 3.3, Table 3)

`scripts/census/P3_members.py` takes TG = the base graph (section above) and K = 231 or 210. For Ω_P:

```
python3 -I scripts/census/P3_members.py aut TG.g6 data/P3/aut_ours.json
python3 -I scripts/census/P3_members.py pairdfs TG.g6 231 gen data/P3/dfs_ours_231_gen.json
python3 -I scripts/census/P3_members.py classify TG.g6 data/P3/dfs_ours_231_gen.json data/P3/aut_ours.json - data/P3 ours_231_gen
```

Run the same with 210, and with the base graph of Ω_M and the tag `mathon`.

- `aut` finds the group: |Aut| = 9, fiber-preserving.
- `pairdfs` is an exhaustive enumeration of the line-piece model, with status `COMPLETE (exhaustive DFS)`. It finds
  5145 solutions in each of the four cases and takes about 15 s.
- `classify` takes one representative per Aut-orbit and checks each one exactly. It takes a few minutes and prints:
  - degree 231: `5145 solutions, 575 classes (Aut-orbit sizes {3: 5, 9: 570})`
  - degree 210: `5145 solutions, 603 classes (Aut-orbit sizes {1: 6, 3: 39, 9: 558})`

The degree-231 outputs `ours_231_gen.g6` and `mathon_231_gen.g6` are byte-identical to the archive files of
degree 231. The degree-210 outputs (SHA-256 `932c7977…5c143` for Ω_P and `82bde2fa…3882c` for Ω_M) hold the members
of degree 210. The archive stores the complement of each, line by line in the same order. `P3_members.py` also has
`--restrict h`, which is the first author's restriction to h-invariant pieces.

### 3. The 2356 graphs are pairwise non-isomorphic (Theorem 3.3)

Compile `gcanon.c` to `bin/gcanon`. Then run `python3 -I scripts/census/members.py OUTDIR`, which needs the four
`classify` outputs in `data/P3/`. `gcanon.c` is a canonical labeling that does not use nauty: pair invariant
W(u,v) = number of edges inside N(u) ∩ N(v), color refinement, and individualization–refinement over all vertices of
the target cell. `members.py` compares the canonical forms within and across the files, against complements, and
against random relabelings (positive control).

Expected output:
- 575/603/575/603 distinct canonical forms, complements likewise;
- 0 common classes for every pair of files, complements included;
- `distinct srg(442,231,120,121) classes ... : 2356`;
- `positive controls ... ALL PASS`.

The reference files of the 29 members known before the census are not in this folder. When they are absent,
`members.py` skips them and says so.

### 4. Two template two-graphs are not isomorphic; at least 23 two-graphs on 442 points (Table 5)

The 23 two-graphs are Ω_P, Ω_M, the 16 sample files `s1_*` and the 5 files `n_*`. The 16 samples change translates,
scalings, slopes, the bijections between arcs and lines, the line sets of the fibers, or the tournament. They were
generated by `P2_family442.py samples 1 16`. The 5 `n_*` change one arc of Ω_P and were generated by
`scripts/census/gauge_choices.py`. `check_template_data.py` checks that both generators reproduce the files.

Two invariants are computed, both without nauty.
- **First invariant: `f5.c` with its driver `f5_invariant.py`.**
  - F(x,y) is the number of 5-point sets containing x and y all of whose ten triples are coherent. Equivalently it
    is the number of 4-cliques through y in the descendant at x.
  - The driver refines colors on the points with F as edge colors. When the coloring is discrete, it also gives a
    complete canonical form.
  - Compile `f5.c` to `bin/f5`, then run `python3 -I scripts/census/f5_invariant.py run NAME TAG` for NAME = `ours`,
    `mathon` and every template file.
  - Controls: `f5_invariant.py control NAME SEED TAG` re-roots at a random point, relabels randomly and switches
    randomly. The original controls used `ours 11`, `mathon 12`, `s1_3 13` and `s1_14 14`.
  - Then run `python3 -I scripts/census/f5_collect.py TAG ...`.
  - Expected: 23 groups, and no pair in different groups shares a point color (rounds 0 and 3). Each control falls
    into the group of its source, and each `g_*` falls into the group of `ours`.
  - |Aut(G_∞)| is 9 for Ω_P and Ω_M, and 3 or 1 for the others. ∞ is a singleton color class, so Aut(Ω) fixes ∞.
  - One run takes about 2 minutes.
- **Second invariant: `descinv.c` with its driver `twograph23_second.py`.**
  - `descinv.c` hashes each descendant D_x labeling-free. For every vertex u it takes the sorted list of
    2 W(u,v) + A(u,v), where W(u,v) is the number of edges of D_x inside the common neighborhood of u and v. It then
    sorts these vertex hashes. The invariant of the two-graph is the multiset of the 442 hashes.
  - The driver also checks srg(441,220,109,110) and S² = 441 I for each of the 23 two-graphs, and it adds the
    complementary two-graphs of Ω_P and Ω_M.
  - Compile `descinv.c` to `bin/descinv`, then run `python3 -I scripts/census/twograph23_second.py compute 0 1` and
    `... merge`. That is 35 jobs of about 1–2 minutes each.
  - Expected in `data/twograph23/result.json`: `all_23_pairwise_distinct`, `all_23_pairwise_disjoint`,
    `all_controls_pass`, `all_template_ok`, `four_signed_pairwise_disjoint` and `all_21_disjoint_from_neg` are all
    `true`.
  - The distinct descendant hashes are 52 for Ω_P and for Ω_M, and 154 or 442 for the others.

Disjoint hash sets mean that no descendant of one two-graph is isomorphic to a descendant of the other, so the two
two-graphs are not isomorphic. A hash collision can only hide a difference, never create one. Ω_P and Ω_M have 52
descendant classes each and share none, as stated in [K833, Computation 6.7].

### 5. Clique numbers: no 22-clique; clique number 13 for 2192 graphs and 14 for 164 (Section 3)

Compile `dcliq.c` to `bin/dcliq`. For each census file F.g6, run these passes:

- `bin/dcliq 22 11 120 < F.g6` (pass A) counts 22-cliques with the Delsarte shortcut. It also gives the exact clique
  number (cap 120 s, never reached) and checks the srg parameters.
- `bin/dcliq 22 -1 < F.g6` (pass B) counts 22-cliques by plain enumeration.
- `bin/dcliq 13 -1` and `bin/dcliq 14 -1` (pass C) count the 13- and 14-cliques.

Each pass takes 0.1–1 s per graph. Controls: `python3 -I scripts/cliques/make_controls.py OUTDIR` writes block
graphs of known Steiner systems (PG(3,2), AG(3,3), PG(4,2), PG(3,3), PG(3,4)) with the expected number of Delsarte
cliques, and the counter must reproduce them.

`data/cliques_per_graph.tsv` is the result, one row per graph. Its columns are file, line, n22 for pass A, n22 for
pass B, clique number, number of 13-cliques, number of 14-cliques, and the srg check. Its SHA-256 is
`01221512e7259a9a5ef53e5c36b80f9aa672797785528f17e9e2ac0891119a4f`.

It shows 0 cliques of size 22 in both passes for all 2356 graphs. The clique number is 13 for 2192 graphs and 14 for
164, and all 164 are members of Ω_M of degree 231. 22 = 1 + 231/11 is the Delsarte bound, so no census graph is
geometric; in particular none is the block graph of a Steiner system S(2,11,221). The script that tallied the pass
outputs into the table is not included. `check_numbers.py` checks the table.

### 6. q = 23 (Proposition 4.10, the two-graphs on 13226 points)

- `python3 -I scripts/q23/c3_arcs.py q OUT.json` runs the per-arc equations (D2), exhaustively over all 2^q
  functions g. For q = 7 it gives 14 solutions per arc, for each of the sets Q, N, {0} ∪ Q, {0} ∪ N and both signs
  of ρ; see `results/arcs_q7.json`. For q = 23 it gives 0 solutions in every case; see `results/arcs_q23.json` and
  `results/arcs_q23.log`.
- `python3 -I scripts/q23/c4_global.py q ARCS.json RHO OUT.json one|all` runs the global step by CP-SAT. For q = 23
  and ρ = ±115 the result is `INFEASIBLE`; see `results/global_q23_rho*.json`. The control at q = 7, ρ = ±21 is in
  `results/global_q7_rho*.json`.
- `python3 -I scripts/q23/c5_build23.py mathon|template` builds the regular two-graph on 13226 points in two ways:
  Mathon's Theorem 4.1 literally, and the Paley(529)-fiber template. Both use the library `w1lib.py`.
  - It checks S² = 13225 I on all 13226² = 174,927,076 entries, and that the descendant at ∞ is
    srg(13225,6612,3305,3306).
  - It writes that descendant as graph6. Our runs gave SHA-256
    `bdccdee566ebe074d2b2dff2573ba7eb9dc6672c301b13db1d84bf46e58780d8` (template) and
    `4311c1384f4be2bd08cbbaf09f3a9439ed42043c2545ae7304cc0182eb9fa7e5` (Mathon).
  - It needs several GB of memory (one 13226 × 13226 float32 matrix is 0.7 GB) and a numpy with an optimized BLAS.
    It takes 1–2 minutes per build.

### 7. The numbers on our side

`python3 -I scripts/checks/check_numbers.py` uses the standard library only. It checks:
- that the archive holds exactly the four graph6 files and README.txt;
- the SHA-256, size and line count of each graph6 file against the table above;
- the clique table: 2356 rows, no 22-clique, srg for every graph, clique numbers 13 and 14 as stated;
- the orbit-count arithmetic of the text: 5145 = 15 · 7³, 575 and 603 orbits, the totals of Table 3, 52 descendant
  classes, and the Delsarte bound 22.

The last line is `ALL OK`. It does not check the graphs themselves (sections 2 and 3 do that) or any number on the
first author's side.

## Not in this folder

- The four `classify` outputs and the run logs. They are reproduced by section 2.
- The reference graphs of the 29 members known before the census.
- The script that complemented the degree-210 members for the archive.
- The tally script of the clique runs.
- The descendant comparison that identified the SageMath graph as a third two-graph.

## Checks run on this folder

On 16 November 2026, on a fresh copy of this folder, in one single-threaded process:
- `check_numbers.py`: ALL OK.
- `check_template_data.py`: ALL OK (28 template choices).
- `c3_arcs.py 7`: 14 solutions per arc for every set T and both signs of ρ, equal to `results/arcs_q7.json`.
- `mathon_literal_blocks.py`: ALL CHECKS PASSED.
- `P2_family442.py base ours`: SHA-256 as above.
- `c1_controls.py`, with that file and without the SageMath graph:
  - (4.15) and (4.16) reproduced at q = 3;
  - S² = 441 I at q = 7 for both builds;
  - `template q=7 g6 == archived srg_441_220_109_110_d1_z9.g6: True`;
  - the SageMath comparison was skipped.
- All C programs compile with `gcc -O2`, and all Python files compile.
