#!/bin/bash
# Re-run D.P.'s checks from a fresh clone.  Usage:  bash plummer/reproduce.sh WORKDIR [q23]
#   WORKDIR : scratch directory for unpacked census and Seidel matrices (created; nothing is written into the repo
#             except paper/figures/aut_table.tex, which is regenerated and must come out unchanged)
#   q23     : also build the q = 23 Seidel matrix (167 MB, a few minutes) and run sporadic_check.py on it.
# Python 3 with numpy, igraph, matplotlib, pynauty; OR-Tools only for m23/reduced.py.
# The full q = 23 h-invariant searches (47 and 45 CPU hours) are not re-run here: see plummer/README.md.
set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd); W=$(mkdir -p "$1" && cd "$1" && pwd); PY=${PYTHON:-python3}
S=$REPO/plummer/scripts; M=$REPO/plummer/m23; D=$REPO/plummer/data
echo "== 1. paper checks (tables, Omega_K data, decoupling facts, members)"
$PY -I $S/check_paper442.py | tail -1
echo "== 2. census properties (srg, |Aut|, orbits, canonical hashes, our 29) against data/census_props.tsv"
mkdir -p $W/census && (cd $W/census && unzip -oq $REPO/koval/data/koval442_linepiece_2356.zip)
NPROC=${NPROC:-8} $PY -I $S/census_props.py $W/census $W/census_props.tsv $D/graphs/srg442_231_120_121_all.g6 $D/graphs/koval442_new15.g6
cmp $W/census_props.tsv $D/census_props.tsv && echo "census_props.tsv identical"
echo "== 3. frames (5145 = 15 * 7^3 per degree and two-graph)"
cp $W/census_props.tsv $W/census/census_props.tsv
$PY -I $S/frames.py $W/census $W/census/frames.tsv | tail -3
cmp $W/census/frames.tsv $D/frames.tsv && echo "frames.tsv identical"
echo "== 4. aut table regenerated"
$PY -I $S/make_figures.py aut >/dev/null && (cd $REPO && git diff --exit-code -- paper/figures/aut_table.tex) && echo "aut_table.tex unchanged"
echo "== 5. Mathon q = 7: Seidel matrix, slope lemma + pair equations on the 14 graphs, figures' numbers"
$PY -I $M/build.py 7 3 $W/S7.npy
$PY -I $S/sporadic_check.py $W/S7.npy | tail -3
echo "== 6. exact h-invariant search at q = 7 (control: 15 members of each degree, the centred members)"
for K in 231 210; do
  $PY -I $M/pairs.py $W/S7.npy 7 $K $W/pairs7_$K.npy | tail -1
  $PY -I $M/structure.py 7 $W/pairs7_$K.npy $W/struct7_$K.npz | tail -1
  $PY -I $M/mitm.py 7 $W/struct7_$K.npz $W/loc7_$K.npy | tail -1
  rm -f $W/glob7_$K.txt; $PY -I $M/glob.py $W/S7.npy 7 $K $W/struct7_$K.npz $W/loc7_$K.npy $W/glob7_$K.txt | tail -1
  echo "q=7 K=$K: $(wc -l < $W/glob7_$K.txt) h-invariant regular members"
  $PY -I $M/verify_hit.py $W/S7.npy 7 $K $W/struct7_$K.npz $W/loc7_$K.npy $W/glob7_$K.txt | tail -1
done
if [ "${2:-}" = q23 ]; then
  echo "== 7. q = 23: Seidel matrix (S^2 = 13225 I) and the slope lemma"
  $PY -I $M/build.py 23 5 $W/S23.npy
  $PY -I $S/sporadic_check.py $W/S7.npy $W/S23.npy | tail -3
fi
echo "done"
