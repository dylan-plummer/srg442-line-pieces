#!/usr/bin/env python3
"""Figures for the joint paper on srg(442,231,120,121) (matplotlib part; the flowchart is figures/flowchart.tikz).
python3 make_figures.py [member] [aut]     -> ../figures/member.pdf (+ .png preview), aut_table.tex

member: one line-piece member of Mathon's two-graph Omega_M, drawn on the nine AG(2,7) fibers, with the fiber pairs of
        the decoupling lemma computed from the Seidel blocks (not assumed). Data: data/member_x.npy, a switching vector x
        with Sx = -21x (checked here), a non-centered member of the frame with the whole fibers on a line of the base.
aut:    table of |Aut| and orbit type of the 2356 census graphs (data/census_props.tsv from census_props.py), by two-graph
        and degree, with how many are among our 29 -> figures/aut_table.tex
Colors: gray whole fiber, white empty, blue lines in a net direction, orange lines outside the net (all figures)."""
import os, sys, json, collections
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
OUT = os.path.join(HERE, "..", "..", "paper", "figures"); os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "serif", "font.size": 8, "axes.linewidth": 0.5, "pdf.fonttype": 42,
                     "mathtext.fontset": "cm"})
BLUE, ORANGE, INK, GREY, LIGHT = "#2a78d6", "#eb6834", "#1b1f27", "#9aa1ad", "#e3e6eb"
SOL = os.path.join(HERE, "..", "data", "member_x.npy")    # a non-centered census member of Omega_M, degree 231 (see README)


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".pdf"), bbox_inches="tight", pad_inches=0.03)
    fig.savefig(os.path.join(OUT, name + ".png"), bbox_inches="tight", pad_inches=0.03, dpi=200)
    plt.close(fig)


def seidel_mathon():
    import check_construction as cc
    M = np.array(cc.build()[0], dtype=np.int64)
    n = 442; A = np.zeros((n, n), dtype=np.int64); A[:441, :441] = M[:441, :441]
    return np.ones((n, n), dtype=np.int64) - np.eye(n, dtype=np.int64) - 2 * A


def form_name(m):
    """b - m a with the sign of the smaller representative: b, b-a, b-2a, b-3a, b+3a, b+2a, b+a"""
    return "b" if m == 0 else (f"b-{m if m > 1 else ''}a" if m <= 3 else f"b+{7 - m if m < 6 else ''}a")


FORMS = [(form_name(m), lambda a, b, m=m: (b - m * a) % 7) for m in range(7)] + [("a", lambda a, b: a % 7)]


def piece(xf):
    """(form name, T) with xf = -1 exactly where form(a, b) in T; None for a constant fiber"""
    if len(set(xf)) == 1: return None
    for name, f in FORMS:
        val = {}
        if all(val.setdefault(f(a, b), xf[7 * a + b]) == xf[7 * a + b] for a in range(7) for b in range(7)):
            return name, sorted(c for c, s in val.items() if s == -1)
    raise ValueError("not a line piece")


def kind(xf, C):
    """W / E for a constant fiber, else net / non: are the lines of the piece cliques of the fiber graph C?"""
    if len(set(xf)) == 1: return "W" if xf[0] == -1 else "E"
    f = dict(FORMS)[piece(xf)[0]]; lv = np.array([f(a, b) for a in range(7) for b in range(7)])
    return "net" if C[(lv[:, None] == lv[None, :]) & ~np.eye(49, dtype=bool)].all() else "non"


KCOL = {"W": GREY, "E": "white", "net": BLUE, "non": ORANGE}       # one color code for every figure


def legend(ax, x0, y0, dy, fs=6.5, what=("W", "E", "net", "non")):
    txt = {"W": "whole fiber", "E": "empty fiber", "net": "lines in a net direction", "non": "lines outside the net"}
    for n, k in enumerate(what):
        ax.add_patch(Rectangle((x0, y0 - n * dy - 0.35 * dy), 0.7 * dy, 0.7 * dy, fc=KCOL[k], ec=INK, lw=0.4))
        ax.text(x0 + 1.1 * dy, y0 - n * dy, txt[k], va="center", fontsize=fs, color=INK)


def fig_member():
    S = seidel_mathon(); x = np.load(SOL).astype(np.int64)
    rho = -21; assert (S @ x == rho * x).all()
    if x[441] == -1: x = -x                                      # oo is not switched
    fib = [x[49 * i:49 * i + 49] for i in range(9)]
    nc = [i for i in range(9) if len(set(fib[i])) > 1]
    pairs = []
    for i in nc:
        for j in nc:
            if i < j:
                z = S[49 * i:49 * i + 49, 49 * j:49 * j + 49] @ (fib[j] - fib[j].mean())
                if np.abs(z - z.mean()).max() > 1e-9: pairs.append((i, j))
    assert len(pairs) * 2 == len(nc) and len({v for p in pairs for v in p}) == len(nc)     # a perfect matching

    const = [i for i in range(9) if i not in nc]
    C = (1 - S[:441, :441] - np.eye(441, dtype=np.int64)) // 2                     # the conference graph
    kd = [kind(fib[i], C[49 * i:49 * i + 49, 49 * i:49 * i + 49]) for i in range(9)]
    fig = plt.figure(figsize=(6.3, 2.75))
    ax0 = fig.add_axes([0.0, 0.0, 0.33, 1.0]); ax = fig.add_axes([0.35, 0.0, 0.65, 1.0])
    for a_ in (ax0, ax): a_.set_aspect("equal"); a_.axis("off")
    # (a) the base points with the tournament i -> i + {1,..,4}; pairs bold, constant fibers filled / hollow
    R = 1.0; pos = {i: np.array([R * np.sin(2 * np.pi * i / 9), R * np.cos(2 * np.pi * i / 9)]) for i in range(9)}
    def arc(a_, i, j, **kw):
        p, q = pos[i], pos[j]; u = (q - p) / np.linalg.norm(q - p)
        a_.add_patch(FancyArrowPatch(p + 0.13 * u, q - 0.13 * u, **kw))
    for i in range(9):
        for d in range(1, 5): arc(ax0, i, (i + d) % 9, arrowstyle="-|>", mutation_scale=5, lw=0.4, color="#c9ced6", zorder=1)
    for i, j in pairs:
        if (j - i) % 9 > 4: i, j = j, i
        arc(ax0, i, j, arrowstyle="-|>", mutation_scale=9, lw=1.5, color=INK, zorder=2)
    for i in range(9):
        ax0.add_patch(plt.Circle(pos[i], 0.1, fc=KCOL[kd[i]], ec=INK, lw=0.7, zorder=3))
        ax0.text(*(pos[i] * 1.24), str(i), ha="center", va="center", fontsize=8, color=INK)
    legend(ax0, -1.3, -1.36, 0.15, fs=6, what=[k for k in ("W", "E", "net", "non") if k in kd])
    ax0.set_xlim(-1.35, 1.35); ax0.set_ylim(-1.95, 1.35)
    ax0.text(-1.3, 1.3, "(a)", fontsize=8, va="top"); ax.text(-0.4, 3 * 10.4 + 0.3, "(b)", fontsize=8, va="top")
    # (b) one row per pair (in arc direction), the constant fibers in the last column
    h = 3.0; cell = 2 * h / 7
    def grid(i, x0, y0):
        for a_ in range(7):
            for b_ in range(7):
                sw = fib[i][7 * a_ + b_] == -1
                ax.add_patch(Rectangle((x0 + b_ * cell, y0 + (6 - a_) * cell), cell, cell, fc=KCOL[kd[i]] if sw else "white",
                                       ec="#c9ced6", lw=0.2, zorder=3))
        ax.add_patch(Rectangle((x0, y0), 2 * h, 2 * h, fc="none", ec=INK, lw=0.6, zorder=4))
        pc = piece(fib[i])
        lab = ("whole fiber" if fib[i][0] == -1 else "empty") if pc is None else rf"${pc[0]}\in\{{{','.join(map(str, pc[1]))}\}}$"
        ax.text(x0 + h, y0 - 0.3, f"fiber {i}\n" + lab, ha="center", va="top", fontsize=6.5, color=INK, linespacing=1.3)
    for r, (i, j) in enumerate(pairs):
        if (j - i) % 9 > 4: i, j = j, i
        y0 = (2 - r) * 10.4 + 2.4
        grid(i, 0, y0); grid(j, 10.5, y0)
        ax.add_patch(FancyArrowPatch((2 * h + 0.4, y0 + h), (10.5 - 0.4, y0 + h), arrowstyle="-|>", mutation_scale=9,
                                     lw=1.5, color=INK, zorder=2))
    for r, i in enumerate(const): grid(i, 23.0, (2 - r) * 10.4 + 2.4)
    ax.plot([21.4, 21.4], [0.4, 3 * 10.4], color="#c9ced6", lw=0.6)
    ax.set_xlim(-0.5, 29.5); ax.set_ylim(-0.6, 3 * 10.4 + 0.4)
    save(fig, "member")
    print("member:", "pairs", pairs, "constant", [i for i in range(9) if i not in nc])


def orbit_type(aut, orb):
    """automorphism orbit lengths as 1^a 3^b 9^c"""
    if orb == "1^442": return "1^{442}"
    c = collections.Counter(json.loads(orb))
    return r"\,".join(f"{k}^{{{c[k]}}}" for k in sorted(c))


def tab_aut():
    """|Aut| and orbit type of the 2356 census graphs; in brackets how many are among our 29 -> figures/aut_table.tex"""
    rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(HERE, "..", "data", "census_props.tsv"))][1:]
    groups = [("mathon_deg231", r"$\Omega_{\mathrm M}$", 231), ("mathon_deg210", r"$\Omega_{\mathrm M}$", 210),
              ("ours_deg231", r"$\Omega_{\mathrm P}$", 231), ("ours_deg210", r"$\Omega_{\mathrm P}$", 210)]
    key = lambda r: (int(r[3]), -len(json.loads(r[4])) if r[4] != "1^442" else -442)
    types = sorted({key(r) + (orbit_type(r[3], r[4]),) for r in rows})
    cnt = {g: collections.Counter(key(r) for r in rows if g in r[0]) for g, *_ in groups}
    mine = {g: collections.Counter(key(r) for r in rows if g in r[0] and r[6]) for g, *_ in groups}
    cell = lambda g, t: (f"{cnt[g][t]}" + (f" ({mine[g][t]})" if mine[g][t] else "")) if cnt[g][t] else "--"
    L = [r"% generated by scripts/make_figures.py aut from data/census_props.tsv", r"\begin{tabular}{ll" + "r" * len(types) + "r}",
         r"\toprule", r" & & " + " & ".join(rf"$|\Aut|={t[0]}$" for t in types) + r" & total \\",
         r"two-graph & degree & " + " & ".join(rf"${t[2]}$" for t in types) + r" & \\", r"\midrule"]
    for g, name, deg in groups:
        L.append(f"{name} & ${deg}$ & " + " & ".join(cell(g, t[:2]) for t in types) +
                 f" & {sum(cnt[g].values())}" + (f" ({sum(mine[g].values())})" if mine[g] else "") + r" \\")
    tot = sum((cnt[g] for g, *_ in groups), collections.Counter()); tm = sum((mine[g] for g, *_ in groups), collections.Counter())
    L += [r"\midrule", r"all & & " + " & ".join(f"{tot[t[:2]]} ({tm[t[:2]]})" if tm[t[:2]] else f"{tot[t[:2]]}" for t in types) +
          f" & {sum(tot.values())} ({sum(tm.values())})" + r" \\", r"\bottomrule", r"\end{tabular}"]
    open(os.path.join(OUT, "aut_table.tex"), "w").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    which = sys.argv[1:] or ["member", "aut"]
    if "member" in which: fig_member()
    if "aut" in which: tab_aut()
