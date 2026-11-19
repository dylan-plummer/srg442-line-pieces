"""Small helpers shared by the search scripts: graph6 encoding/decoding and an exact strongly-regular check."""
import os
import numpy as np

GAP = os.environ.get("GAP", "gap")       # GAP 4 with the GRAPE package
GP = os.environ.get("GP", "gp")          # PARI/GP


def graph6(A):
    """graph6 string of a 0/1 adjacency matrix"""
    v = A.shape[0]
    bits = []
    for j in range(1, v):
        for i in range(j):
            bits.append(int(A[i, j]))
    bits += [0] * (-len(bits) % 6)
    if v < 63: head = [63 + v]
    elif v < 258048: head = [126, 63 + (v >> 12 & 63), 63 + (v >> 6 & 63), 63 + (v & 63)]
    out = head + [63 + int("".join(map(str, bits[t:t + 6])), 2) for t in range(0, len(bits), 6)]
    return bytes(out).decode()


def g6(s):
    """adjacency matrix (numpy int64) of a graph6 string"""
    s = s.strip(); b = [c - 63 for c in s.encode()]
    if b[0] == 63: n = (b[1] << 12) | (b[2] << 6) | b[3]; b = b[4:]
    else: n = b[0]; b = b[1:]
    bits = [(x >> (5 - i)) & 1 for x in b for i in range(6)]
    A = np.zeros((n, n), dtype=np.int64); t = 0
    for j in range(1, n):
        for i in range(j): A[i, j] = A[j, i] = bits[t]; t += 1
    return A


def is_srg(A):
    """(v, k, lambda, mu) if A is strongly regular, else None"""
    v = len(A); k = int(A[0].sum()); A2 = A @ A; I = np.eye(v, dtype=bool)
    if not (A.sum(1) == k).all(): return None
    la = set(A2[(A == 1)]); mu = set(A2[(A == 0) & ~I])
    return (v, k, la.pop(), mu.pop()) if len(la) == 1 and len(mu) == 1 else None
