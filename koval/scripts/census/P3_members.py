#!/usr/bin/env python3
"""P3 : regular members of a fibred regular two-graph on 442 points whose switching set is a
union of LINE PIECES in every fibre (the reduction used by D. Plummer, github.com/dylan-plummer/srg442, there with
h-invariant pieces only; this is our own implementation, written from the method description, none of his code run).

Setting.  Omega = two-graph of M + {oo}, M a conference graph srg(441,220,109,110) with 9 fibres of 49 vertices,
vertex 49*B + 7*a + b <-> (B, (a,b) in F_7^2), oo = 441.  Seidel S (-1 adjacent, +1 non-adjacent, oo isolated in M+oo).
A switching vector x in {+-1}^442 (x_oo = +1) gives a K-regular graph iff S x = rho x, rho = 441 - 2K.

PIECE MODEL (a RESTRICTION, not a theorem): on every fibre B, the switched set is {u : g_B(l_d(u)) = 1} for one of
the 8 parallel classes l_d of F_7^2 and an arbitrary g_B : F_7 -> {0,1} (constant g = empty/full fibre).
  --restrict h   : additionally g_B(2t) = g_B(t) (Plummer's 50 options per fibre: unions of {0}, QR, NQR)
  --restrict gen : arbitrary g_B (2 + 8*126 = 1010 options per fibre)
Encoding: one-hot y[B,d] (direction), bools g[B,t], w[B,d,t] = y[B,d] & g[B,t];  sw_u = sum_d w[B,d,l_d(u)];
equations (S - rho I) sw = (S - rho I) 1 / 2, rewritten on line sums.  Constant g is pinned to direction 0
(symmetry of the encoding).  All solutions are enumerated with blocking clauses; after every solution the whole orbit
under the given automorphism group of M (which fixes oo) is blocked, so every solution is a new Aut-orbit.
COMPLETE = CP-SAT proves the blocked model INFEASIBLE.

Subcommands (run from koval/; solve needs OR-Tools, the others numpy + pynauty):
  recover TG.g6 MEMBERS.g6 OUT.json        # switching vectors of known members (direct labelling or nauty iso)
  aut TG.g6 OUT.json                       # generators of Aut(M) (pynauty), closed into the group
  structure TG.g6                          # Fourier check: which (fibre,line) sectors S couples
  solve TG.g6 K RESTRICT AUT.json OUT.json [TLIM] [WORKERS] [MAXSOL] [CHECK.json]
  verify TG.g6 SOLS.json KNOWN.g6|- OUTDIR TAG
"""
import sys, os, json, time, hashlib, itertools, collections
import numpy as np

N0, N = 441, 442
INF = 441


def g6_decode(s):
    s = s.strip()
    if s[0] == '~':
        n = ((ord(s[1]) - 63) << 12) | ((ord(s[2]) - 63) << 6) | (ord(s[3]) - 63); data = s[4:]
    else:
        n = ord(s[0]) - 63; data = s[1:]
    bits = []
    for ch in data:
        v = ord(ch) - 63; bits.extend((v >> (5 - i)) & 1 for i in range(6))
    A = np.zeros((n, n), dtype=np.int64); k = 0
    for j in range(n):
        for i in range(j):
            if bits[k]: A[i, j] = A[j, i] = 1
            k += 1
    return A


def g6_encode(A):
    n = len(A); bits = [int(A[i, j]) for j in range(n) for i in range(j)]
    while len(bits) % 6: bits.append(0)
    s = '~' + chr(63 + (n >> 12)) + chr(63 + ((n >> 6) & 63)) + chr(63 + (n & 63)) if n > 62 else chr(63 + n)
    return s + ''.join(chr(63 + int(''.join(map(str, bits[t:t + 6])), 2)) for t in range(0, len(bits), 6))


def read_g6_file(f):
    return [g6_decode(l) for l in open(f) if l.strip()]


def seidel_plus_inf(M):
    """Seidel matrix of M + isolated oo (N x N)."""
    S = np.ones((N, N), dtype=np.int64); S[:N0, :N0] = 1 - 2 * M; np.fill_diagonal(S, 0)
    return S


def seidel(A):
    S = 1 - 2 * A; np.fill_diagonal(S, 0); return S


LINES = [(0, 1)] + [(1, m) for m in range(7)]   # l(a,b) = l0*a + l1*b ... we use the 8 directions as forms below
# The 8 linear forms on F_7^2 up to scalars: b - m a (m in F_7) and a.  Form index d: 0..6 -> b - d a, 7 -> a.
def form(d, a, b):
    return (b - d * a) % 7 if d < 7 else a % 7


FORM_NAMES = [f"b-{d}a" for d in range(7)] + ["a"]
LV = np.array([[form(d, u // 7, u % 7) for u in range(49)] for d in range(8)])   # LV[d,u] = value of form d at u


def srg_check(A):
    """exact: regular + A^2 = kI + lam A + mu (J - I - A); returns (ok, params)."""
    n = len(A); deg = A.sum(1)
    if not np.all(deg == deg[0]): return False, 'not regular'
    k = int(deg[0]); A2 = A @ A
    off = ~np.eye(n, dtype=bool)
    lam_vals = set(A2[(A == 1) & off].tolist()); mu_vals = set(A2[(A == 0) & off].tolist())
    if len(lam_vals) != 1 or len(mu_vals) > 1: return False, (k, sorted(lam_vals)[:5], sorted(mu_vals)[:5])
    return True, (n, k, lam_vals.pop(), mu_vals.pop() if mu_vals else None)


def nauty_graph(A):
    import pynauty
    n = len(A); adj = {i: np.nonzero(A[i])[0].tolist() for i in range(n)}
    return pynauty.Graph(n, adjacency_dict=adj)


def cert_hash(A):
    import pynauty
    return hashlib.sha256(pynauty.certificate(nauty_graph(A))).hexdigest()


def piece_of(xB):
    """xB in {+-1}^49 (switch indicator = (1-x)/2).  Return list of (d, g) descriptions with sw = g(form d)."""
    sw = (1 - xB) // 2
    if sw.min() == sw.max(): return [("const", int(sw[0]))]
    out = []
    for d in range(8):
        g = {}
        ok = True
        for u in range(49):
            t = LV[d, u]
            if g.setdefault(t, sw[u]) != sw[u]: ok = False; break
        if ok: out.append((FORM_NAMES[d], ''.join(str(int(g[t])) for t in range(7))))
    return out


def h_invariant(gstr):
    return all(gstr[t] == gstr[(2 * t) % 7] for t in range(7))


# ---------------------------------------------------------------- recover
def recover(tg, members, out):
    M = read_g6_file(tg)[0]; S0 = seidel_plus_inf(M); res = []
    certM = None; perms_cache = []
    for idx, G in enumerate(read_g6_file(members)):
        ok, par = srg_check(G); SG = seidel(G); how = 'direct'
        x = SG[INF] * S0[INF]; x[INF] = 1
        for pi, perm in enumerate(perms_cache):         # labellings found for earlier members
            if np.array_equal(SG, S0 * np.outer(x, x)): break
            SG = seidel(G)[np.ix_(perm, perm)]; x = SG[INF] * S0[INF]; x[INF] = 1; how = f'cached-iso#{pi}'
        if not np.array_equal(SG, S0 * np.outer(x, x)):
            SG = seidel(G); how = 'direct'
            # labelling differs: find w with descendant(G, w) ~ M, map w -> oo via canonical labellings
            import pynauty
            if certM is None:
                gM = nauty_graph(M); certM = pynauty.certificate(gM); labM = pynauty.canon_label(gM)
            found = None
            for w in range(N):
                D = SG[w].copy(); D[w] = 1
                Sw = SG * np.outer(D, D)                      # row w all +1 (w isolated)
                keep = [v for v in range(N) if v != w]
                Aw = ((1 - Sw[np.ix_(keep, keep)]) // 2); np.fill_diagonal(Aw, 0)
                gw = nauty_graph(Aw)
                if pynauty.certificate(gw) == certM:
                    labw = pynauty.canon_label(gw)
                    # canonical position c: Aw vertex labw[c] <-> M vertex labM[c]
                    perm = np.empty(N, dtype=np.int64)               # perm[new index in M+oo labelling] = G vertex
                    for c in range(N0): perm[labM[c]] = keep[labw[c]]
                    perm[INF] = w
                    SGp = SG[np.ix_(perm, perm)]
                    x = SGp[INF] * S0[INF]; x[INF] = 1
                    assert np.array_equal(SGp, S0 * np.outer(x, x)), 'iso map failed'
                    found = w; perms_cache.append(perm); break
            assert found is not None, f'member {idx}: no descendant isomorphic to M'
            how = f'iso(w={found})'
        rho = int(x @ S0 @ x) // N
        assert np.array_equal(S0 @ x, rho * x), 'recovered x is not an eigenvector'
        K = (N - 1 - rho) // 2
        pieces = [piece_of(x[49 * B:49 * B + 49]) for B in range(9)]
        in_gen = all(len(p) > 0 for p in pieces)
        in_h = in_gen and all(p[0][0] == 'const' or any(h_invariant(g) for _, g in p) for p in pieces)
        res.append(dict(idx=idx, srg=list(par) if ok else str(par), how=how, K=K, x=x.tolist(),
                        pieces=[[list(t) for t in p] for p in pieces], in_gen=in_gen, in_h=in_h,
                        sizes=[int(((1 - x[49 * B:49 * B + 49]) // 2).sum()) for B in range(9)]))
        print(idx, how, 'K', K, 'gen', in_gen, 'h', in_h, res[-1]['sizes'], flush=True)
    json.dump(res, open(out, 'w'))


# ---------------------------------------------------------------- aut
def aut(tg, out):
    import pynauty
    M = read_g6_file(tg)[0]
    gens, mantissa, exp, orbits, norb = pynauty.autgrp(nauty_graph(M))
    order = int(round(mantissa * 10 ** exp))
    gens = [tuple(g) for g in gens]; ident = tuple(range(N0)); G = {ident}; frontier = [ident]
    while frontier:
        nf = []
        for p in frontier:
            for g in gens:
                q = tuple(g[p[i]] for i in range(N0))
                if q not in G: G.add(q); nf.append(q)
        frontier = nf
        assert len(G) <= 10000
    assert len(G) == order, (len(G), order)
    for p in G:
        P = np.array(p); assert np.array_equal(M[np.ix_(P, P)], M)
    # does every element map fibres to fibres and lines to lines (piece model invariant)?
    fib = []
    for p in G:
        P = np.array(p); fib.append(all(len({P[49 * B + u] // 49 for u in range(49)}) == 1 for B in range(9)))
    print('|Aut M| =', order, 'orbits', norb, 'fibre-preserving elements', sum(fib), flush=True)
    json.dump(dict(order=order, group=[list(p) for p in G]), open(out, 'w'))


# ---------------------------------------------------------------- structure (Fourier sectors)
def structure(tg):
    M = read_g6_file(tg)[0]; S = seidel_plus_inf(M).astype(float)
    # orthonormal basis of zero-mean functions of form d on a fibre: 6 vectors per (B,d)
    Pm = {}
    for d in range(8):
        F = np.zeros((49, 7))
        for u in range(49): F[u, LV[d, u]] = 1
        F = F - F.mean(0); U, s, _ = np.linalg.svd(F, full_matrices=False); Pm[d] = U[:, :6]
    coupling = {}
    for B in range(9):
        for d in range(8):
            v = np.zeros((N, 6)); v[49 * B:49 * B + 49] = Pm[d]
            Sv = S @ v; partners = []
            for B2 in range(9):
                for d2 in range(8):
                    w = np.zeros((N, 6)); w[49 * B2:49 * B2 + 49] = Pm[d2]
                    nrm = np.linalg.norm(w.T @ Sv)
                    if nrm > 1e-6: partners.append((B2, FORM_NAMES[d2], round(float(nrm), 3)))
            coupling[(B, d)] = partners
            print(B, FORM_NAMES[d], '->', partners, flush=True)


# ---------------------------------------------------------------- solve
def solve(tg, K, restrict, autf, out, tlim=3600.0, workers=4, maxsol=100000, checkf=None):
    from ortools.sat.python import cp_model
    M = read_g6_file(tg)[0]; S0 = seidel_plus_inf(M); rho = N - 1 - 2 * K
    R = S0 - rho * np.eye(N, dtype=np.int64); c0 = R.sum(1); assert np.all(c0 % 2 == 0)
    rhs = c0 // 2
    # line-sum coefficients: C[i, B, d, t] = sum_{u in fibre B, form d(u) = t} R[i, 49B+u]
    C = np.zeros((N, 9, 8, 7), dtype=np.int64)
    for B in range(9):
        blk = R[:, 49 * B:49 * B + 49]
        for d in range(8):
            for t in range(7):
                C[:, B, d, t] = blk[:, LV[d] == t].sum(1)
    m = cp_model.CpModel()
    y = [[m.NewBoolVar(f"y{B}_{d}") for d in range(8)] for B in range(9)]
    g = [[m.NewBoolVar(f"g{B}_{t}") for t in range(7)] for B in range(9)]
    w = [[[m.NewBoolVar(f"w{B}_{d}_{t}") for t in range(7)] for d in range(8)] for B in range(9)]
    for B in range(9):
        m.AddExactlyOne(y[B])
        for d in range(8):
            for t in range(7):
                m.AddImplication(w[B][d][t], y[B][d]); m.AddImplication(w[B][d][t], g[B][t])
                m.AddBoolOr([w[B][d][t], y[B][d].Not(), g[B][t].Not()])
            m.Add(sum(w[B][d]) == 7 * 0 + sum(g[B])).OnlyEnforceIf(y[B][d])
        # constant g -> direction 0
        allsame = m.NewBoolVar(f"c{B}")
        m.AddBoolAnd([g[B][t] for t in range(7)]).OnlyEnforceIf([allsame, g[B][0]])
        m.AddBoolAnd([g[B][t].Not() for t in range(7)]).OnlyEnforceIf([allsame, g[B][0].Not()])
        # allsame <-> all g equal
        diffs = []
        for t in range(1, 7):
            e = m.NewBoolVar(""); m.Add(g[B][t] != g[B][0]).OnlyEnforceIf(e); m.Add(g[B][t] == g[B][0]).OnlyEnforceIf(e.Not())
            diffs.append(e)
        m.AddBoolOr(diffs + [allsame]); [m.AddImplication(allsame, e.Not()) for e in diffs]
        m.AddImplication(allsame, y[B][0])
        if restrict == 'h':
            for t in range(7): m.Add(g[B][t] == g[B][(2 * t) % 7])
    seen_rows = set(); neq = 0
    for i in range(N):
        key = C[i].tobytes() + int(rhs[i]).to_bytes(8, 'little', signed=True)
        if key in seen_rows: continue
        seen_rows.add(key); neq += 1
        terms = [(int(C[i, B, d, t]), w[B][d][t]) for B in range(9) for d in range(8) for t in range(7) if C[i, B, d, t]]
        m.Add(sum(c * v for c, v in terms) == int(rhs[i]))
    group = [np.array(p) for p in json.load(open(autf))['group']] if autf != '-' else [np.arange(N0)]

    def x_of(Bd):
        sw = np.zeros(N, dtype=np.int64)
        for B, (d, gs) in enumerate(Bd):
            sw[49 * B:49 * B + 49] = np.array([gs[t] for t in LV[d]])
        return 1 - 2 * sw

    def assign_of(x):
        """piece assignment (d, g) per fibre in canonical encoding, or None if x not in the model."""
        a = []
        for B in range(9):
            p = piece_of(x[49 * B:49 * B + 49])
            if not p: return None
            if p[0][0] == 'const': a.append((0, (p[0][1],) * 7)); continue
            d = FORM_NAMES.index(p[0][0]); gs = tuple(int(c) for c in p[0][1])
            if restrict == 'h' and not h_invariant(p[0][1]): return None
            a.append((d, gs))
        return a

    def block(a):
        lits = []
        for B, (d, gs) in enumerate(a):
            lits.append(y[B][d].Not()); lits += [g[B][t].Not() if gs[t] else g[B][t] for t in range(7)]
        m.AddBoolOr(lits)

    log = dict(K=K, rho=rho, restrict=restrict, neq=neq, aut_order=len(group), solutions=[], rounds=[])
    print(f"== K={K} rho={rho} restrict={restrict}: {neq} distinct equations, |Aut| {len(group)}", flush=True)
    if checkf:   # must-find membership: fix each known x and check FEASIBLE in this exact model
        for rec in json.load(open(checkf)):
            if rec['K'] != K: continue
            a = assign_of(np.array(rec['x']))
            if a is None: print(f"   control {rec['idx']}: NOT in model", flush=True); continue
            hint = []
            for B, (d, gs) in enumerate(a):
                hint += [y[B][d] == 1] + [g[B][t] == gs[t] for t in range(7)]
            mm = m.Clone(); [mm.Add(h) for h in hint]
            s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = 60; s.parameters.num_workers = 1
            print(f"   control {rec['idx']}: model {s.StatusName(s.Solve(mm))}", flush=True)
    t0 = time.time(); status = None
    while len(log['solutions']) < maxsol:
        s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = tlim; s.parameters.num_workers = workers
        t = time.time(); st = s.Solve(m); dt = time.time() - t; log['rounds'].append([s.StatusName(st), round(dt, 2)])
        if st == cp_model.INFEASIBLE: status = 'COMPLETE'; break
        if st not in (cp_model.FEASIBLE, cp_model.OPTIMAL): status = f'INCOMPLETE ({s.StatusName(st)} at {tlim:.0f}s)'; break
        a = [(next(d for d in range(8) if s.Value(y[B][d])), tuple(int(s.Value(g[B][t])) for t in range(7))) for B in range(9)]
        x = x_of(a); assert np.array_equal(S0 @ x, rho * x), 'solver output is not an eigenvector'
        orbit = set(); nblk = 0
        for P in group:
            xp = x.copy(); xp[P] = x[:N0]          # x'(P(v)) = x(v)
            key = xp.tobytes()
            if key in orbit: continue
            orbit.add(key); assert np.array_equal(S0 @ xp, rho * xp)
            ap = assign_of(xp)
            if ap is not None: block(ap); nblk += 1
        desc = ' | '.join('-' if gs.count(gs[0]) == 7 and gs[0] == 0 else ('all' if gs.count(1) == 7 else f"{FORM_NAMES[d]}:{''.join(map(str, gs))}") for d, gs in a)
        log['solutions'].append(dict(x=x.tolist(), desc=desc, orbit=len(orbit), blocked=nblk, t=round(dt, 2)))
        print(f"   {len(log['solutions'])}: {desc}  (orbit {len(orbit)}, {nblk} in model, {dt:.1f}s)", flush=True)
    log['status'] = status or f'INCOMPLETE (maxsol {maxsol})'; log['time'] = round(time.time() - t0, 1)
    print(f"   STATUS {log['status']}: {len(log['solutions'])} Aut-orbits in {log['time']}s", flush=True)
    json.dump(log, open(out, 'w'))


# ---------------------------------------------------------------- solvefull (NO piece restriction; sparse line-sum basis)
def solvefull(tg, K, autf, inv, out, tlim=600.0, workers=4):
    """Unrestricted Sx = rho x with the equations rewritten as LINE-SUMMED rows: row (B,d,t) = sum of the rows of
    (S - rho I) over the line {l_d = t} of fibre B.  By the sector lemma (structure) such a row is a function of l_d on
    fibre B, of l' on the one coupled fibre B', and constant on every other fibre, so it is written with 14 line-sum
    variables + 9 fibre totals.  inv = '-' or comma list of group-element indices whose invariance is imposed."""
    from ortools.sat.python import cp_model
    M = read_g6_file(tg)[0]; S0 = seidel_plus_inf(M); rho = N - 1 - 2 * K
    R = S0 - rho * np.eye(N, dtype=np.int64)
    m = cp_model.CpModel(); x = [m.NewBoolVar(f"s{u}") for u in range(N0)]   # s = switched (x = 1 - 2s)
    group = [np.array(p) for p in json.load(open(autf))['group']]
    if inv != '-':
        for e in map(int, inv.split(',')):
            P = group[e]
            for u in range(N0): m.Add(x[u] == x[int(P[u])])
    L = {}
    for B in range(9):
        for d in range(8):
            for t in range(7):
                L[B, d, t] = sum(x[49 * B + u] for u in range(49) if LV[d, u] == t)
    T = [sum(x[49 * B:49 * B + 49]) for B in range(9)]
    neq = 0; maxterms = 0
    rows = []
    for B in range(9):
        for d in range(8):
            for t in range(7):
                rows.append(R[[49 * B + u for u in range(49) if LV[d, u] == t]].sum(0))
    rows.append(R[INF])
    for row in rows:
        # express row . sw with sw = s (switch indicator), x = 1 - 2 sw:  row.1 - 2 row.sw = 0
        const = int(row.sum()); assert const % 2 == 0
        expr = []; nterm = 0
        for B in range(9):
            r = row[49 * B:49 * B + 49]
            if r.min() == r.max():
                if r[0]: expr.append(int(r[0]) * T[B]); nterm += 1
                continue
            done = False
            for d in range(8):
                vals = {}
                if all(vals.setdefault(LV[d, u], r[u]) == r[u] for u in range(49)):
                    expr += [int(vals[t]) * L[B, d, t] for t in range(7) if vals[t]]; nterm += 7; done = True; break
            if not done:
                expr += [int(r[u]) * x[49 * B + u] for u in range(49) if r[u]]; nterm += 49
        m.Add(sum(expr) == const // 2); neq += 1; maxterms = max(maxterms, nterm)
    print(f"== full K={K} rho={rho} inv={inv}: {neq} line-summed equations, max {maxterms} line/total terms", flush=True)
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = tlim; s.parameters.num_workers = workers
    t = time.time(); st = s.Solve(m); dt = time.time() - t
    res = dict(K=K, inv=inv, status=s.StatusName(st), time=round(dt, 1))
    if st in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        xv = np.array([1 - 2 * s.Value(x[u]) for u in range(N0)] + [1])
        assert np.array_equal(S0 @ xv, rho * xv); res['x'] = xv.tolist()
        res['pieces'] = [[list(q) for q in piece_of(xv[49 * B:49 * B + 49])] for B in range(9)]
    print(f"   {res['status']} in {dt:.1f}s", res.get('pieces', ''), flush=True)
    json.dump(res, open(out, 'w'))


# ---------------------------------------------------------------- pairdfs (independent exact enumeration of the piece model)
def pairdfs(tg, K, restrict, out):
    """Exact enumeration of the piece model WITHOUT a solver, via the sector lemma (checked numerically by `structure`
    and re-checked here on every leaf by the full S x = rho x): S preserves the trivial sector (fibre-constant vectors
    + oo) and couples every zero-mean line sector (B, l) to exactly one sector (B', l') of another fibre.  A piece fibre
    along l has its zero-mean part in (B, l) alone, so non-constant fibres come in pairs {B, B'} along the coupled lines
    (a matching of K_9); each pair must solve its own sector equation; constant fibres and oo satisfy the 10 x 10
    trivial-sector equation.  Enumerates all matchings x valid pair-patterns x constants and keeps x_oo = +1."""
    M = read_g6_file(tg)[0]; S0 = seidel_plus_inf(M); rho = N - 1 - 2 * K
    R = S0 - rho * np.eye(N, dtype=np.int64)
    # coupled line of fibre pair (B, B'): the direction d at B whose line-sum row hits fibre B' non-constantly
    coup = {}
    for B in range(9):
        for d in range(8):
            row = R[[49 * B + u for u in range(49) if LV[d, u] == 0]].sum(0)
            hit = [B2 for B2 in range(9) if B2 != B and row[49 * B2:49 * B2 + 49].min() != row[49 * B2:49 * B2 + 49].max()]
            assert len(hit) == 1, (B, d, hit); coup[B, hit[0]] = d
    assert len(coup) == 72
    gs_all = [tuple((c >> t) & 1 for t in range(7)) for c in range(1, 127)]
    if restrict == 'h': gs_all = [g for g in gs_all if all(g[t] == g[(2 * t) % 7] for t in range(7))]
    def zm(B, d, g):   # zero-mean part of x on fibre B for piece (d, g), times 7 (integer)
        x = np.array([1 - 2 * g[t] for t in LV[d]]); return 49 * x - x.sum()
    valid = {}
    for B in range(9):
        for B2 in range(B + 1, 9):
            d, d2 = coup[B, B2], coup[B2, B]; lst = []
            cols = list(range(49 * B, 49 * B + 49)) + list(range(49 * B2, 49 * B2 + 49))
            Rsub = R[:, cols]
            Z1 = {g: zm(B, d, g) for g in gs_all}; Z2 = {g: zm(B2, d2, g) for g in gs_all}
            negb = {}
            for g2 in gs_all: negb.setdefault((-(Rsub[:, 49:] @ Z2[g2])).tobytes(), []).append(g2)
            for g in gs_all:
                for g2 in negb.get((Rsub[:, :49] @ Z1[g]).tobytes(), []): lst.append((g, g2))
            valid[B, B2] = lst
    npairs = sum(len(v) for v in valid.values())
    print(f"== pairdfs K={K} rho={rho} restrict={restrict}: valid pair patterns per fibre pair:",
          sorted({len(v) for v in valid.values()}), 'total', npairs, flush=True)
    sols = []; leaves = 0
    def matchings(free):
        if not free: yield []; return
        B = free[0]; rest = free[1:]
        for m_ in matchings(rest): yield m_                     # B unmatched
        for B2 in rest:
            if valid[B, B2]:
                for m_ in matchings([b for b in rest if b != B2]): yield [(B, B2)] + m_
    # trivial sector: x-bar = (s_B / 49 on fibre B, x_oo); (S - rho I) x-bar = 0  <=>  W (s_0..s_8, 49 x_oo) = 0
    W = np.zeros((10, 10), dtype=np.int64)
    for B in range(9):
        for B2 in range(9):
            rsum = R[49 * B, 49 * B2:49 * B2 + 49].sum()
            assert np.all(R[49 * B:49 * B + 49, 49 * B2:49 * B2 + 49].sum(1) == rsum) or B == B2
            W[B, B2] = S0[49 * B, 49 * B2:49 * B2 + 49].sum()
        W[B, B] -= rho; W[B, 9] = 1                       # x_oo coefficient (times 49 below)
    W[9, :9] = 1; W[9, 9] = -rho                          # oo row: sum_u x_u - rho x_oo
    def ssum_of(B, d, g): return int(sum(1 - 2 * g[t] for t in LV[d]))
    spairs = {p: {} for p in valid}
    for (B, B2), lst in valid.items():
        for g, g2 in lst:
            spairs[B, B2].setdefault((ssum_of(B, coup[B, B2], g), ssum_of(B2, coup[B2, B], g2)), []).append((g, g2))
    t0 = time.time()
    for mt in matchings(list(range(9))):
        matched = {b for p in mt for b in p}; un = [b for b in range(9) if b not in matched]
        keys = [list(spairs[p].keys()) for p in mt]
        for skey in itertools.product(*keys, *[[(49,), (-49,)]] * len(un)):
            leaves += 1
            s = np.zeros(10, dtype=np.int64); s[9] = 49            # x_oo = +1 (times 49 in the fibre rows)
            for (B, B2), (a, b) in zip(mt, skey[:len(mt)]): s[B], s[B2] = a, b
            for B, (c,) in zip(un, skey[len(mt):]): s[B] = c
            tri = W[:9] @ s; tri_oo = s[:9].sum() - rho
            if np.any(tri) or tri_oo: continue
            for choice in itertools.product(*[spairs[p][k] for p, k in zip(mt, skey[:len(mt)])]):
                sw = np.zeros(N, dtype=np.int64)
                for (B, B2), (g, g2) in zip(mt, choice):
                    sw[49 * B:49 * B + 49] = [g[t] for t in LV[coup[B, B2]]]
                    sw[49 * B2:49 * B2 + 49] = [g2[t] for t in LV[coup[B2, B]]]
                for B, (c,) in zip(un, skey[len(mt):]): sw[49 * B:49 * B + 49] = 0 if c == 49 else 1
                x = 1 - 2 * sw
                assert np.array_equal(S0 @ x, rho * x), 'sector lemma violated'
                sols.append(x.tolist())
    print(f"   pairdfs: {leaves} leaves, {len(sols)} solutions (all x, no Aut reduction) in {time.time() - t0:.1f}s", flush=True)
    json.dump(dict(K=K, restrict=restrict, status='COMPLETE (exhaustive DFS)', solutions=[dict(x=x, desc='') for x in sols]),
              open(out, 'w'))


# ---------------------------------------------------------------- verify
def verify(tg, solf, knownf, outdir, tag, autf='-', srgcheck='all'):
    """1) group the solutions into orbits of Aut(M) (= Aut(Omega), which fixes oo; two members of one switching
    class are isomorphic iff their switching vectors are in one orbit), 2) every orbit representative: exact srg
    check via A^2 (float BLAS, entries < 2^53 so exact) and a pynauty certificate; distinct orbits must have distinct
    certificates (asserted: cross-check of the orbit argument), 3) compare with the known members."""
    M = read_g6_file(tg)[0]; S0 = seidel_plus_inf(M)
    group = [np.array(p) for p in json.load(open(autf))['group']] if autf != '-' else [np.arange(N0)]
    known = {}
    if knownf != '-':
        for i, G in enumerate(read_g6_file(knownf)): known[cert_hash(G)] = i
    sols = json.load(open(solf))['solutions']; keyidx = {}
    for i, sol in enumerate(sols): keyidx[np.array(sol['x'], dtype=np.int8).tobytes()] = i
    orbit_of = [-1] * len(sols); reps = []
    for i, sol in enumerate(sols):
        if orbit_of[i] >= 0: continue
        x = np.array(sol['x'], dtype=np.int8); o = len(reps); reps.append(i); size = 0
        for P in group:
            xp = x.copy(); xp[P] = x[:N0]; j = keyidx.get(xp.tobytes())
            assert j is not None, 'solution set not closed under Aut'
            if orbit_of[j] < 0: orbit_of[j] = o; size += 1
    print(f"{tag}: {len(sols)} solutions, {len(reps)} Aut-orbits", flush=True)
    seen = {}; lines = []; rep = []; t0 = time.time()
    for o, i in enumerate(reps):
        x = np.array(sols[i]['x']); Sx = S0 * np.outer(x, x); A = (1 - Sx) // 2; np.fill_diagonal(A, 0)
        Af = A.astype(np.float64); A2 = np.rint(Af @ Af).astype(np.int64)
        deg = A.sum(1); k = int(deg[0]); assert np.all(deg == k)
        off = ~np.eye(N, dtype=bool); lv = set(A2[(A == 1) & off].tolist()); mv = set(A2[(A == 0) & off].tolist())
        assert len(lv) == 1 and len(mv) == 1, (i, k, lv, mv); par = (N, k, lv.pop(), mv.pop())
        h = cert_hash(A); assert h not in seen, f'orbits {seen[h]} and {o} give isomorphic graphs'
        seen[h] = o; lines.append(g6_encode(A))
        sizes = [int(((1 - x[49 * B:49 * B + 49]) // 2).sum()) for B in range(9)]
        pieces = [piece_of(x[49 * B:49 * B + 49])[0] for B in range(9)]
        hinv = all(p[0] == 'const' or h_invariant(p[1]) for p in pieces)
        rep.append(dict(orbit=o, sol=i, params=list(par), cert=h, known=known.get(h), sizes=sizes, h_invariant=hinv,
                        pieces=[list(p) for p in pieces]))
        if o % 50 == 0: print(f"  {o}/{len(reps)} {time.time() - t0:.0f}s", flush=True)
    os.makedirs(outdir, exist_ok=True); f = os.path.join(outdir, f"{tag}.g6")
    open(f, 'w').write(''.join(l + '\n' for l in lines))
    sha = hashlib.sha256(open(f, 'rb').read()).hexdigest()
    kh = sorted({r['known'] for r in rep if r['known'] is not None})
    json.dump(dict(file=f, sha256=sha, solutions=len(sols), classes=len(lines), known_hit=kh, known_total=len(known),
                   h_invariant_classes=sum(r['h_invariant'] for r in rep), reports=rep), open(f + '.json', 'w'), indent=1)
    print(f"{tag}: {len(sols)} solutions -> {len(lines)} classes (h-invariant {sum(r['h_invariant'] for r in rep)}); "
          f"known hit {kh} of {len(known)}; {f} sha256 {sha}  ({time.time() - t0:.0f}s)", flush=True)


def k4_count(A):
    """number of 4-cliques (float BLAS, exact for these sizes)."""
    Af = A.astype(np.float64); tot = 0.0
    for v in range(len(A)):
        nb = np.nonzero(A[v])[0]; Bm = Af[np.ix_(nb, nb)]; tot += float((Bm * (Bm @ Bm)).sum()) / 6
    return int(round(tot / 4))


def descendant(SG, w):
    D = SG[w].copy(); D[w] = 1; Sw = SG * np.outer(D, D); keep = [v for v in range(len(SG)) if v != w]
    Aw = (1 - Sw[np.ix_(keep, keep)]) // 2; np.fill_diagonal(Aw, 0); return Aw, keep


def classify(tg, solf, autf, knownf, outdir, tag):
    """Isomorphism classes of the solutions = Aut(M)-orbits (Aut(Omega) = Aut(M) fixes oo, so two members are
    isomorphic iff their normalised switching vectors lie in one Aut(M)-orbit).  Every orbit representative is
    checked exactly (regular + A^2 identity, float BLAS on 0/1 matrices, exact below 2^53).  Known members: the
    switching vector is read off directly when the labelling agrees with M + oo; otherwise oo is located by the
    4-clique count of the descendants, the descendant is mapped to M by pynauty canonical labelling (one call)."""
    import pynauty
    M = read_g6_file(tg)[0]; S0 = seidel_plus_inf(M)
    group = [np.array(p) for p in json.load(open(autf))['group']]
    sols = json.load(open(solf))['solutions']; keyidx = {np.array(s['x'], dtype=np.int8).tobytes(): i for i, s in enumerate(sols)}
    orbit_of = [-1] * len(sols); reps = []; osize = []
    for i, sol in enumerate(sols):
        if orbit_of[i] >= 0: continue
        x = np.array(sol['x'], dtype=np.int8); o = len(reps); reps.append(i); c = 0
        for P in group:
            xp = x.copy(); xp[P] = x[:N0]; j = keyidx.get(xp.tobytes()); assert j is not None, 'not Aut-closed'
            if orbit_of[j] < 0: orbit_of[j] = o; c += 1
        osize.append(c)
    t0 = time.time(); lines = []; off = ~np.eye(N, dtype=bool); params = set()
    for i in reps:
        x = np.array(sols[i]['x']); A = (1 - S0 * np.outer(x, x)) // 2; np.fill_diagonal(A, 0)
        Af = A.astype(np.float64); A2 = np.rint(Af @ Af).astype(np.int64); deg = A.sum(1); k = int(deg[0])
        assert np.all(deg == k); lv = set(A2[(A == 1) & off].tolist()); mv = set(A2[(A == 0) & off].tolist())
        assert len(lv) == 1 and len(mv) == 1; params.add((N, k, lv.pop(), mv.pop())); lines.append(g6_encode(A))
    os.makedirs(outdir, exist_ok=True); f = os.path.join(outdir, f"{tag}.g6"); open(f, 'w').write(''.join(l + '\n' for l in lines))
    sha = hashlib.sha256(open(f, 'rb').read()).hexdigest()
    print(f"{tag}: {len(sols)} solutions, {len(reps)} classes (Aut-orbit sizes {dict(sorted(collections.Counter(osize).items()))}), "
          f"all exact srg {sorted(params)}; {f} sha256 {sha} ({time.time() - t0:.0f}s)", flush=True)
    hits = []
    if knownf != '-':
        perm = None
        for idx, G in enumerate(read_g6_file(knownf)):
            SG = seidel(G); x = SG[INF] * S0[INF]; x[INF] = 1; how = 'direct'
            if not np.array_equal(SG, S0 * np.outer(x, x)) and perm is not None:
                SG = SG[np.ix_(perm, perm)]; x = SG[INF] * S0[INF]; x[INF] = 1; how = 'relabelled'
            if not np.array_equal(SG, S0 * np.outer(x, x)):
                SG = seidel(G); gM = nauty_graph(M); certM = pynauty.certificate(gM); labM = pynauty.canon_label(gM)
                # oo is the unique fixed point of Aut(Omega) (order 9, 51 + 1 orbits), and Aut(G) <= Aut(Omega) fixes oo:
                # try the fixed points of Aut(G) first (all points if Aut(G) = 1)
                gensG = pynauty.autgrp(nauty_graph(G))[0]
                cand = [w for w in range(N) if all(gg[w] == w for gg in gensG)]
                print(f"  known {idx}: |gens Aut(G)| = {len(gensG)}, {len(cand)} fixed-point candidates for oo", flush=True)
                for w in cand:
                    Aw, keep = descendant(SG, w)
                    gw = nauty_graph(Aw)
                    if pynauty.certificate(gw) != certM: continue
                    labw = pynauty.canon_label(gw); perm = np.empty(N, dtype=np.int64)
                    for c in range(N0): perm[labM[c]] = keep[labw[c]]
                    perm[INF] = w; break
                assert perm is not None
                SG = SG[np.ix_(perm, perm)]; x = SG[INF] * S0[INF]; x[INF] = 1; how = f'relabelled (oo = vertex {w})'
                assert np.array_equal(SG, S0 * np.outer(x, x))
            j = keyidx.get(x.astype(np.int8).tobytes())
            hits.append(dict(idx=idx, how=how, in_model=j is not None, cls=orbit_of[j] if j is not None else None))
            print(f"  known {idx}: {how}, in piece model: {j is not None}, class {hits[-1]['cls']}", flush=True)
    json.dump(dict(file=f, sha256=sha, solutions=len(sols), classes=len(reps), orbit_sizes=osize, params=sorted(params),
                   known=hits, reps=[sols[i]['x'] for i in reps]), open(f + '.json', 'w'))


if __name__ == '__main__':
    cmd, a = sys.argv[1], sys.argv[2:]
    if cmd == 'recover': recover(*a)
    elif cmd == 'aut': aut(*a)
    elif cmd == 'structure': structure(*a)
    elif cmd == 'solve':
        solve(a[0], int(a[1]), a[2], a[3], a[4], float(a[5]) if len(a) > 5 else 3600, int(a[6]) if len(a) > 6 else 4,
              int(a[7]) if len(a) > 7 else 100000, a[8] if len(a) > 8 else None)
    elif cmd == 'verify': verify(*a)
    elif cmd == 'pairdfs': pairdfs(a[0], int(a[1]), a[2], a[3])
    elif cmd == 'classify': classify(*a)
    elif cmd == 'solvefull':
        solvefull(a[0], int(a[1]), a[2], a[3], a[4], float(a[5]) if len(a) > 5 else 600, int(a[6]) if len(a) > 6 else 4)
