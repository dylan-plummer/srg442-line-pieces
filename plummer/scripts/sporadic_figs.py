#!/usr/bin/env python3
"""Figures for the sporadicity section (Mathon's family q = s^2 - 2).
python3 -I sporadic_figs.py S7.npy        -> ../figures/sporadic_gaps.pdf, ../figures/sporadic_reduction.pdf (+ .png)

S7.npy: Seidel matrix of Mathon's two-graph for q = 7 (int8, 442 x 442, oo last, vertex 49 i + 7 a + b), as built by
paut_runs/m23/code/build.py; it is checked here against srg442's Construction (S^2 = 441 I and the member A_1).
Every number drawn is recomputed here: the gaps of the centered line pieces for q = 7, 23, 47, 79, 167, 359, 727, the
partner pair (fiber 1, slope 1) <-> (fiber 7, slope 6) of A_1 and its line sums, and the conference matrix sigma.
The search counts in panel (b) of the reduction figure are copied from the logs of glob.py: plummer/m23/results/glob7_231.log, glob7_210.log (q = 7) and
h_slices_6670/, h_slices_6555/ (q = 23)."""
import os, sys, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Patch

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "..", "paper", "figures"); os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "serif", "font.size": 8, "axes.linewidth": 0.5, "pdf.fonttype": 42,
                     "mathtext.fontset": "cm", "xtick.major.width": 0.5, "ytick.major.width": 0.5})
BLUE, ORANGE, INK, GREY, LIGHT, MUTED = "#2a78d6", "#eb6834", "#1b1f27", "#9aa1ad", "#e3e6eb", "#5b6270"


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".pdf"), bbox_inches="tight", pad_inches=0.03)
    fig.savefig(os.path.join(OUT, name + ".png"), bbox_inches="tight", pad_inches=0.03, dpi=200)
    plt.close(fig)


# ----------------------------------------------------------------------------------------------- arithmetic in F_q
def chi_table(q):
    Q = {(x * x) % q for x in range(1, q)}
    return np.array([0] + [1 if y in Q else -1 for y in range(1, q)]), Q


def kappa_values(q, T):
    """values of (chi' * 1_T)(y), chi' = chi - delta_0, on y = 0, y in Q, y in N; None if not constant on Q or N"""
    chi, Q = chi_table(q)
    chip = chi.copy(); chip[0] = -1
    D = np.array([sum(chip[(y - x) % q] for x in T) for y in range(q)])
    vq = {D[y] for y in range(1, q) if y in Q}; vn = {D[y] for y in range(1, q) if y not in Q}
    assert len(vq) == 1 and len(vn) == 1
    return D[0], vq.pop(), vn.pop()


def centred_gaps(q):
    """gaps of the two-valued centered line pieces T (unions of {0}, Q, N, not empty, not all), for chi' and chi'(-.)"""
    _, Q = chi_table(q); N = set(range(1, q)) - Q
    gaps = set()
    for e0 in (0, 1):
        for eq in (0, 1):
            for en in (0, 1):
                T = ({0} if e0 else set()) | (Q if eq else set()) | (N if en else set())
                if not 0 < len(T) < q: continue
                for TT in (T, {(-x) % q for x in T}):          # chi'(-.) * 1_T is chi' * 1_{-T}, reflected
                    vals = set(kappa_values(q, TT))
                    if len(vals) == 2: gaps.add(max(vals) - min(vals))
    return gaps


# ----------------------------------------------------------------------------------------------- figure 1: the gaps
def fig_gaps():
    fam = [(3, 7), (5, 23), (7, 47), (9, 79), (13, 167), (15, 223), (19, 359), (21, 439), (27, 727)]
    for s, q in fam:
        assert q == s * s - 2
        g = centred_gaps(q) if q <= 450 else {2, (q + 1) // 2}   # closed form; computed directly up to q = 439
        assert g == {2, (q + 1) // 2}, (q, g)
    fig, ax = plt.subplots(figsize=(5.2, 2.35))
    xs = np.arange(len(fam))
    for x, (s, q) in zip(xs, fam):
        need = [s - 1, s + 1]; have = [2, (q + 1) // 2]
        ax.plot([x, x], need, color=INK, lw=2.2, solid_capstyle="round", zorder=2)
        ax.scatter([x, x], need, s=22, color=INK, zorder=3)
        ax.scatter([x, x], have, s=38, facecolor="white", edgecolor=BLUE, lw=1.3, marker="D", zorder=4)
        if s == 3:
            ax.add_patch(plt.Rectangle((x - 0.32, 1.65), 0.64, 3.0, fill=False, ec=MUTED, lw=0.6, ls=(0, (2, 1.5))))
    ax.set_yscale("log", base=2)
    ax.set_yticks([2, 4, 8, 16, 32, 64, 128, 256, 512]); ax.set_yticklabels(["2", "4", "8", "16", "32", "64", "128", "256", "512"])
    ax.set_ylim(1.5, 560)
    ax.set_xticks(xs); ax.set_xticklabels([f"$q={q}$\n$s={s}$" for s, q in fam])
    ax.set_xlim(-0.6, len(fam) - 0.4)
    ax.set_ylabel("gap")
    ax.grid(axis="y", color=LIGHT, lw=0.5, zorder=0)
    for side in ("top", "right"): ax.spines[side].set_visible(False)
    ax.tick_params(axis="x", length=0)
    h = [Line2D([], [], color=INK, lw=2.2, marker="o", ms=4, label="needed by a partner pair: $s-1$ and $s+1$"),
         Line2D([], [], color=BLUE, lw=0, marker="D", ms=5, mfc="white", mew=1.3,
                label="offered by a line piece: $2$ or $(q+1)/2$")]
    ax.legend(handles=h, loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=2, frameon=False, fontsize=7.2, handlelength=1.6, columnspacing=1.2)
    save(fig, "sporadic_gaps")


# ----------------------------------------------------------------------------------------------- figure 2: reduction
A1 = {0: (6, "0"), 1: (1, "Q"), 2: "all", 3: (1, "Q"), 4: "all", 5: (6, "0"), 6: "all", 7: (6, "0"), 8: (1, "Q")}


def member_A1(q=7):
    _, Q = chi_table(q); N = set(range(1, q)) - Q
    sets = {"0": {0}, "Q": Q, "N": N}
    x = np.ones(9 * q * q + 1, dtype=np.int64)
    for i, pc in A1.items():
        for a in range(q):
            for b in range(q):
                if pc == "all" or (pc != "all" and (b - pc[0] * a) % q in sets[pc[1]]):
                    x[i * q * q + q * a + b] = -1
    return x


def line_sums(xf, m, q=7):
    """(L(0), L(1), L(3)): sums of the fiber function over the lines of slope m with intercept 0, 1 (in Q), 3 (in N)"""
    lab = np.array([(a if m == "inf" else (b - m * a)) % q for a in range(q) for b in range(q)])
    return tuple(int(xf[lab == c].sum()) for c in (0, 1, 3))


def halfline_signs(xf, q=7):
    """origin sign and, for each slope m, the signs on {t d_m : t in Q} and {t d_m : t in N}"""
    _, Q = chi_table(q)
    out = {}
    for m in list(range(q)) + ["inf"]:
        d = (0, 1) if m == "inf" else (1, m)
        sg = {}
        for t in range(1, q):
            v = xf[q * (t * d[0] % q) + (t * d[1] % q)]
            sg.setdefault(t in Q, set()).add(int(v))
        assert all(len(v) == 1 for v in sg.values())          # h-invariant
        out[m] = (sg[True].pop(), sg[False].pop())
    return int(xf[0]), out


def partner(S, i, m, q=7):
    b = q * q
    lab = lambda mm: np.array([(a if mm == "inf" else (bb - mm * a)) % q for a in range(q) for bb in range(q)])
    f = (lab(m) == 1).astype(float); f -= f.mean()
    hits = []
    for j in range(9):
        if j == i: continue
        g = S[j * b:(j + 1) * b, i * b:(i + 1) * b] @ f; g -= g.mean()
        if np.abs(g).max() > 1e-9:
            for mm in list(range(q)) + ["inf"]:
                L = lab(mm)
                if all(np.ptp(g[L == c]) < 1e-9 for c in range(q)): hits.append((j, mm))
    assert len(hits) == 1, hits
    return hits[0]


def draw_wheel(ax, cx, cy, R, e0, signs, net, title, sub):
    order = [0, 1, 2, 3, 4, 5, 6, "inf"]
    for k, m in enumerate(order):
        th = math.radians(90 - 22.5 * k)
        col = BLUE if m in net else ORANGE
        for half, sgn in ((0, signs[m][0]), (1, signs[m][1])):     # half 0: t in Q, towards the label
            dx, dy = math.cos(th + math.pi * half), math.sin(th + math.pi * half)
            if sgn == -1:
                ax.plot([cx + 0.09 * R * dx, cx + R * dx], [cy + 0.09 * R * dy, cy + R * dy], color=col, lw=2.6,
                        solid_capstyle="round", zorder=3)
            else:
                ax.plot([cx + 0.09 * R * dx, cx + R * dx], [cy + 0.09 * R * dy, cy + R * dy], color=GREY, lw=0.6,
                        zorder=2)
        lx, ly = cx + 1.2 * R * math.cos(th), cy + 1.2 * R * math.sin(th)
        ax.text(lx, ly, "$\\infty$" if m == "inf" else f"${m}$", ha="center", va="center", fontsize=7,
                color=BLUE if m in net else ORANGE)
    ax.scatter([cx], [cy], s=26, zorder=5, facecolor=INK if e0 == -1 else "white", edgecolor=INK, lw=0.8)
    ax.text(cx, cy + 1.52 * R, title, ha="center", va="bottom", fontsize=7.6, color=INK)
    ax.text(cx, cy - 1.45 * R, sub, ha="center", va="top", fontsize=7.2, color=INK, linespacing=1.35)


def fig_reduction(S):
    q, s = 7, 3
    _, Qset = chi_table(q)
    net = set(Qset) | {0}                                      # slopes in {0} u Q: net directions (multiplier -q)
    x = member_A1()
    assert np.array_equal(S @ x, -21 * x), "A_1 is not a solution"
    # sigma: conference matrix of order q+1
    slopes = list(range(q)) + ["inf"]
    chi, _ = chi_table(q)
    def ell(m, d): return d[0] if m == "inf" else (d[1] - m * d[0]) % q
    dvec = {m: ((0, 1) if m == "inf" else (1, m)) for m in slopes}
    sig = np.array([[0 if m == mm else chi[ell(m, dvec[mm])] for mm in slopes] for m in slopes])
    assert np.array_equal(sig @ sig.T, q * np.eye(q + 1, dtype=int))
    # the pair (fiber 1, slope 1) <-> (fiber 7, slope 6), its multipliers and line sums
    j, mm = partner(S.astype(float), 1, 1)
    assert (j, mm) == (7, 6)
    b = q * q; xf1, xf7 = x[b:2 * b], x[7 * b:8 * b]
    L1, L7 = line_sums(xf1, 1), line_sums(xf7, 6)
    uvw = lambda L: (L[0] - L[1], L[0] - L[2], L[1] - L[2])
    (u1, v1, w1), (u7, v7, w7) = uvw(L1), uvw(L7)
    assert u7 == ((s - 1) // 2) * w1 or u7 == -((s - 1) // 2) * w1          # rho = -qs: u_plus = +-(s-1)/2 w_minus
    assert abs(v1) == ((s + 1) // 2) * abs(w7)                               #             v_minus = +-(s+1)/2 w_plus
    e1, sg1 = halfline_signs(xf1); e7, sg7 = halfline_signs(xf7)
    # closed forms of Proposition: L(0) = e0 + (q-1) alpha_m, L(Q) = A - alpha_m + P_m, L(N) = A - alpha_m - P_m
    for e0, sg, m, L in ((e1, sg1, 1, L1), (e7, sg7, 6, L7)):
        al = {k: (v[0] + v[1]) // 2 for k, v in sg.items()}; de = {k: (v[0] - v[1]) // 2 for k, v in sg.items()}
        A = sum(al.values()); idx = slopes.index(m)
        P = sum(sig[idx, slopes.index(k)] * de[k] for k in slopes)
        assert L == (e0 + (q - 1) * al[m], A - al[m] + P, A - al[m] - P), (m, L)

    fig = plt.figure(figsize=(6.5, 2.55))
    axa = fig.add_axes([0.0, 0.0, 0.45, 1.0]); axa.set_xlim(0, 10); axa.set_ylim(0, 5.6); axa.axis("off")
    axa.set_aspect("equal")
    draw_wheel(axa, 2.2, 3.0, 1.35, e1, sg1, net, "fiber $1$: $\\ell_1^{-1}(Q)$",
               f"$(\\Lambda(0),\\Lambda(Q),\\Lambda(N))=({L1[0]},{L1[1]},{L1[2]})$\n$v_1={v1},\\ w_1={w1}$")
    draw_wheel(axa, 7.6, 3.0, 1.35, e7, sg7, net, "fiber $7$: $\\ell_6^{-1}(\\{0\\})$",
               f"$(\\Lambda(0),\\Lambda(Q),\\Lambda(N))=({L7[0]},{L7[1]},{L7[2]})$\n$u_7={u7},\\ w_7={w7}$")
    axa.add_patch(FancyArrowPatch((3.95, 3.0), (5.85, 3.0), arrowstyle="<|-|>", mutation_scale=7, color=MUTED, lw=0.7))
    axa.text(4.9, 3.25, "partners", ha="center", va="bottom", fontsize=7, color=MUTED)
    axa.text(4.9, 2.75, "$u_7=w_1$\n$v_1=-2\\,w_7$", ha="center", va="top", fontsize=7, color=INK, linespacing=1.3)
    fig.text(0.0, 0.98, "(a)", fontsize=8.5, va="top")

    # (b) the search, log10 counts (from the glob.py logs; both degrees together)
    stages = ["$h$-invariant vectors", "after the local\nconditions", "nodes of the\nexact search", "regular members"]
    c7 = [153 * math.log10(2), math.log10(2 * 6256 ** 9), math.log10(418869 + 417337), math.log10(30)]
    c23 = [1225 * math.log10(2), math.log10(206632 ** 25 + 189142 ** 25), math.log10(284271840 + 276655619), None]
    lab7 = ["$10^{%d}$" % round(c7[0]), "$10^{%d}$" % round(c7[1]), "$8\\cdot10^{5}$", "$30$"]
    lab23 = ["$10^{%d}$" % round(c23[0]), "$10^{%d}$" % round(c23[1]), "$6\\cdot10^{8}$", "$0$"]
    axb = fig.add_axes([0.6, 0.17, 0.38, 0.72])
    y = np.arange(len(stages))[::-1]; hgt = 0.34
    for k in range(len(stages)):
        axb.barh(y[k] + hgt / 2 + 0.02, c7[k], height=hgt, color=INK, zorder=3)
        axb.text(c7[k] + 6, y[k] + hgt / 2 + 0.02, lab7[k], va="center", fontsize=7, color=INK)
        if c23[k] is not None:
            axb.barh(y[k] - hgt / 2 - 0.02, c23[k], height=hgt, color=GREY, zorder=3)
            axb.text(c23[k] + 6, y[k] - hgt / 2 - 0.02, lab23[k], va="center", fontsize=7, color=INK)
        else:
            axb.text(6, y[k] - hgt / 2 - 0.02, lab23[k] + " (none exist)", va="center", fontsize=7, color=INK)
    axb.set_yticks(y); axb.set_yticklabels(stages, fontsize=7.3)
    axb.set_xlim(0, 420); axb.set_ylim(-0.62, 3.62); axb.set_xticks([0, 100, 200, 300, 400])
    axb.set_xlabel("$\\log_{10}$ (count)")
    axb.grid(axis="x", color=LIGHT, lw=0.5, zorder=0)
    for side in ("top", "right"): axb.spines[side].set_visible(False)
    axb.tick_params(axis="y", length=0)
    axb.legend(handles=[Patch(color=INK, label="$q=7$ (442 points)"), Patch(color=GREY, label="$q=23$ (13226 points)")],
               loc="lower right", frameon=False, fontsize=7)
    fig.text(0.455, 0.98, "(b)", fontsize=8.5, va="top")
    save(fig, "sporadic_reduction")
    return dict(L1=L1, L7=L7, uvw1=(u1, v1, w1), uvw7=(u7, v7, w7))


if __name__ == "__main__":
    S = np.load(sys.argv[1]).astype(np.int64)
    assert S.shape == (442, 442) and np.array_equal(S @ S, 441 * np.eye(442, dtype=np.int64))
    fig_gaps()
    print(fig_reduction(S))
    print("ok")
