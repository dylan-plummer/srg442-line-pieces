# srg442-line-pieces

Paper, data and verification programs for

> D. Plummer and A. Kovaľ, *Strongly regular graphs on 442 points from line pieces in two conference two-graphs*,
> 19 November 2026.

We give a complete census of the line-piece regular members of two conference two-graphs on 442 points, Mathon's Ω_M
and the Paley(49)-fibre two-graph Ω_P: 2356 pairwise non-isomorphic strongly regular graphs with parameters
(442,231,120,121) or (442,210,99,100), organised into 15 frames per two-graph and degree. We also show that line pieces
are sporadic in Mathon's family q = s² − 2, and prove non-existence statements at q = 23 (13226 points).

## Contents

| folder | contents |
|---|---|
| `paper/` | `main.tex` (with `sections/` and `figures/`) and the compiled `main.pdf` |
| `koval/` | A. Kovaľ's census archive `data/koval442_linepiece_2356.zip` (four graph6 files; full SHA-256 in `koval/README.md`), template data, and the programs for the census, the comparison of the two-graphs, cliques and the q = 23 statements |
| `plummer/` | D. Plummer's programs and data: census properties, frames, figures, the sporadicity checks and the q = 23 h-invariant search (`m23/`); `reproduce.sh` reruns them |

The table of programs at the end of the paper lists which program supports which statement.

## Build and check

```
cd paper && pdflatex main && pdflatex main
bash plummer/reproduce.sh /path/to/scratch          # D.P.'s checks, about 10 minutes; add "q23" for the q = 23 Seidel matrix
python3 -I koval/scripts/checks/check_numbers.py   # A.K.'s number checks (see koval/README.md for the input paths)
```

Requirements: TeX Live (amsmath, booktabs, tikz, hyperref); Python 3 with numpy, python-igraph, matplotlib and
pynauty (OR-Tools for the CP-SAT programs); gcc for the C programs in `koval/scripts`.
