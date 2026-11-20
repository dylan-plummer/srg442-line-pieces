/* nauty-free canonical form of (strongly regular) graphs, n <= 512.
   Pair invariant W(u,v) = #edges inside N(u) cap N(v); edge color code(u,v) = 2 W(u,v) + A(u,v).
   Color refinement with 64-bit multiset hashing (deterministic => an isomorphism invariant; a hash collision can
   only MERGE colors, never separate isomorphic vertices).  Canonical form by individualization-refinement:
   target cell = smallest non-singleton cell (ties: smallest color), branch over ALL its vertices, keep the
   lexicographically smallest color-ordered adjacency matrix  =>  canonical labeling (exact, collision-proof:
   equal outputs are equal matrices; isomorphic inputs give equal outputs).
   usage: gcanon [c] < file.g6     ('c' = use the complement of every graph)
   output per graph: "G idx n srg=(k,l,m)|NOTSRG r0=<hash> ms=<hash> cls=<#stable classes w/o individualization>
                      leaves=<#discrete leaves> ord= v0 v1 ... "                                                   */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#define MAXN 512
#define WD 8
typedef uint64_t u64;
static int n; static u64 ADJ[MAXN][WD]; static int A[MAXN][MAXN]; static int CODE[MAXN][MAXN];
static long long leaves, maxleaves = 200000;
static u64 mix(u64 z) { z += 0x9e3779b97f4a7c15ULL; z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
                        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; return z ^ (z >> 31); }
static int cmpu(const void *a, const void *b) { u64 x = *(const u64 *)a, y = *(const u64 *)b; return x < y ? -1 : x > y; }
static int ncls(const u64 *c) { u64 s[MAXN]; memcpy(s, c, n * 8); qsort(s, n, 8, cmpu);
                                int k = 1; for (int i = 1; i < n; i++) k += s[i] != s[i - 1]; return k; }
static void refine(u64 *c) {
    int k = ncls(c); u64 nc[MAXN];
    for (int round = 0; round < 1000; round++) {
        for (int x = 0; x < n; x++) { u64 s = 0;
            for (int y = 0; y < n; y++) if (y != x) s += mix(c[y] * 0x100000001b3ULL + (u64)CODE[x][y] * 0x51ULL + 7);
            nc[x] = mix(c[x] ^ mix(s)); }
        int k2 = ncls(nc); memcpy(c, nc, n * 8);
        if (k2 == k) return; k = k2;
    }
}
static int best_have; static int best_ord[MAXN]; static unsigned char best_mat[MAXN][MAXN / 8 + 1];
static int ord_of(const u64 *c, int *ord) {   /* discrete coloring -> order by color */
    for (int i = 0; i < n; i++) ord[i] = i;
    for (int i = 1; i < n; i++) { int v = ord[i], j = i; while (j > 0 && c[ord[j - 1]] > c[v]) { ord[j] = ord[j - 1]; j--; } ord[j] = v; }
    return 0;
}
static void leaf(const u64 *c) {
    int ord[MAXN]; ord_of(c, ord); leaves++;
    static unsigned char mat[MAXN][MAXN / 8 + 1]; memset(mat, 0, sizeof(mat));
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (A[ord[i]][ord[j]]) mat[i][j >> 3] |= 1 << (j & 7);
    int better = !best_have;
    if (!better) { int r = memcmp(mat, best_mat, sizeof(mat)); better = r < 0; }
    if (better) { best_have = 1; memcpy(best_mat, mat, sizeof(mat)); memcpy(best_ord, ord, sizeof(ord)); }
}
static void search(u64 *c) {
    refine(c);
    if (ncls(c) == n) { leaf(c); return; }
    if (leaves > maxleaves) return;
    /* target cell: smallest size > 1, ties smallest color */
    u64 s[MAXN]; memcpy(s, c, n * 8); qsort(s, n, 8, cmpu);
    int bsz = n + 1; u64 bcol = 0;
    for (int i = 0; i < n;) { int j = i; while (j < n && s[j] == s[i]) j++;
        if (j - i > 1 && j - i < bsz) { bsz = j - i; bcol = s[i]; } i = j; }
    for (int v = 0; v < n; v++) if (c[v] == bcol) {
        u64 c2[MAXN]; memcpy(c2, c, n * 8); c2[v] = mix(c[v] ^ 0xabcdef12345ULL); search(c2);
    }
}
static int readg6(char *s) {
    int p = 0; if (s[0] == '~') { n = ((s[1] - 63) << 12) | ((s[2] - 63) << 6) | (s[3] - 63); p = 4; } else { n = s[0] - 63; p = 1; }
    if (n > MAXN) return -1;
    int bitpos = 0, len = strlen(s);
    memset(A, 0, sizeof(A));
    for (int j = 1; j < n; j++) for (int i = 0; i < j; i++) {
        int byte = p + bitpos / 6; if (byte >= len) return -2;
        int b = ((s[byte] - 63) >> (5 - bitpos % 6)) & 1; bitpos++;
        A[i][j] = A[j][i] = b; }
    return 0;
}
static char line[1 << 20];
int main(int argc, char **argv) {
    int comp = argc > 1 && argv[1][0] == 'c'; int idx = 0;
    while (fgets(line, sizeof(line), stdin)) {
        int L = strlen(line); while (L && (line[L - 1] == '\n' || line[L - 1] == '\r' || line[L - 1] == ' ')) line[--L] = 0;
        if (!L) continue;
        if (readg6(line)) { printf("G %d PARSEERROR\n", idx++); continue; }
        if (comp) for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (i != j) A[i][j] ^= 1;
        memset(ADJ, 0, sizeof(ADJ));
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (A[i][j]) ADJ[i][j >> 6] |= 1ULL << (j & 63);
        /* srg check + W */
        int k = -1, lam = -1, mu = -1, ok = 1;
        for (int i = 0; i < n && ok; i++) { int d = 0; for (int w = 0; w < WD; w++) d += __builtin_popcountll(ADJ[i][w]);
            if (k < 0) k = d; else if (d != k) ok = 0; }
        for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) {
            u64 C[WD]; int cn = 0; for (int w = 0; w < WD; w++) { C[w] = ADJ[u][w] & ADJ[v][w]; cn += __builtin_popcountll(C[w]); }
            if (A[u][v]) { if (lam < 0) lam = cn; else if (cn != lam) ok = 0; } else { if (mu < 0) mu = cn; else if (cn != mu) ok = 0; }
            long long e2 = 0;
            for (int w = 0; w < WD; w++) { u64 m = C[w]; while (m) { int c = w * 64 + __builtin_ctzll(m); m &= m - 1;
                for (int x = 0; x < WD; x++) e2 += __builtin_popcountll(C[x] & ADJ[c][x]); } }
            CODE[u][v] = CODE[v][u] = (int)(e2 / 2) * 2 + A[u][v];
        }
        u64 c[MAXN]; u64 s0 = 0;
        for (int x = 0; x < n; x++) { u64 s = 0; for (int y = 0; y < n; y++) if (y != x) s += mix((u64)CODE[x][y] * 0x51ULL + 3); c[x] = mix(s); }
        { u64 t[MAXN]; memcpy(t, c, n * 8); qsort(t, n, 8, cmpu); for (int i = 0; i < n; i++) s0 = mix(s0 ^ t[i]); }
        u64 cs[MAXN]; memcpy(cs, c, n * 8); refine(cs); int cls = ncls(cs);
        u64 ms = 0; { u64 t[MAXN]; memcpy(t, cs, n * 8); qsort(t, n, 8, cmpu); for (int i = 0; i < n; i++) ms = mix(ms ^ t[i]); }
        best_have = 0; leaves = 0; search(cs);
        if (ok) printf("G %d %d srg=%d,%d,%d", idx, n, k, lam, mu); else printf("G %d %d NOTSRG", idx, n);
        printf(" r0=%016llx ms=%016llx cls=%d leaves=%lld%s ord=", (unsigned long long)s0, (unsigned long long)ms, cls, leaves,
               leaves > maxleaves ? " TRUNCATED" : "");
        for (int i = 0; i < n; i++) printf(i ? " %d" : "%d", best_ord[i]);
        printf("\n"); fflush(stdout); idx++;
    }
    return 0;
}
