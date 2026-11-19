/* descinv.c -- second, nauty-free invariant program for regular two-graphs on n points.
 *
 * Independent of f5.c (which counts K4 per vertex of each descendant).
 * Input (stdin): n, then n lines of n characters '0'/'1' ('1' where the Seidel entry is -1; diagonal ignored).
 * For every point x the descendant D_x is formed (switch so that x is isolated, delete x).  For every pair u != v of
 * D_x let W(u,v) = number of edges of D_x inside N(u) cap N(v), and code(u,v) = 2 W(u,v) + A(u,v).
 * vertex hash  h(u)  = H(sorted list of code(u,v), v != u)
 * descendant hash I(D_x) = H(sorted list of h(u), u in D_x)
 * Output: one line per point x: "x I(D_x) Wmin Wmax".
 * I(D_x) depends only on the isomorphism class of D_x (all steps are labelling-free), so an isomorphism of two-graphs
 * Omega -> Omega' mapping x to y forces I(D_x) = I'(D_y).  Hash collisions can only make different descendants look
 * equal, never make isomorphic descendants look different; so "different multisets of I" is a proof of
 * non-isomorphism, while equality proves nothing.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define MAXN 512
#define WORDS 8
typedef uint64_t u64;

static int n;
static signed char S[MAXN][MAXN];
static u64 adj[MAXN][WORDS];
static unsigned short code[MAXN][MAXN];

static u64 mix(u64 h, u64 v) { /* splitmix64-style combine */
    h ^= v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    h ^= h >> 30; h *= 0xbf58476d1ce4e5b9ULL; h ^= h >> 27; h *= 0x94d049bb133111ebULL; h ^= h >> 31;
    return h;
}
static int cmp_us(const void *a, const void *b) { return (int)*(const unsigned short *)a - (int)*(const unsigned short *)b; }
static int cmp_u64(const void *a, const void *b) { u64 x = *(const u64 *)a, y = *(const u64 *)b; return x < y ? -1 : x > y; }

int main(void) {
    if (scanf("%d", &n) != 1 || n < 3 || n > MAXN) { fprintf(stderr, "bad n\n"); return 1; }
    static char line[MAXN + 8];
    for (int i = 0; i < n; i++) {
        if (scanf("%s", line) != 1 || (int)strlen(line) != n) { fprintf(stderr, "bad row %d\n", i); return 1; }
        for (int j = 0; j < n; j++) S[i][j] = (i == j) ? 0 : (line[j] == '1' ? -1 : 1);
    }
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
        if (S[i][j] != S[j][i]) { fprintf(stderr, "not symmetric\n"); return 1; }
    int m = n - 1;
    static int idx[MAXN];
    static unsigned short row[MAXN];
    static u64 vh[MAXN];
    for (int x = 0; x < n; x++) {
        int k = 0;
        for (int v = 0; v < n; v++) if (v != x) idx[k++] = v;
        memset(adj, 0, sizeof(adj));
        for (int a = 0; a < m; a++) for (int b = 0; b < m; b++) {
            if (a == b) continue;
            int u = idx[a], v = idx[b];
            int s = S[u][v] * S[x][u] * S[x][v];          /* switched so that row x is all +1 */
            if (s == -1) adj[a][b >> 6] |= 1ULL << (b & 63);
        }
        int wmin = 1 << 30, wmax = -1;
        for (int a = 0; a < m; a++) {
            code[a][a] = 0;
            for (int b = a + 1; b < m; b++) {
                u64 C[WORDS]; int w2 = 0;
                for (int t = 0; t < WORDS; t++) C[t] = adj[a][t] & adj[b][t];
                for (int t = 0; t < WORDS; t++) {
                    u64 word = C[t];
                    while (word) {
                        int bit = __builtin_ctzll(word); word &= word - 1;
                        int w = (t << 6) | bit;
                        for (int s = 0; s < WORDS; s++) w2 += __builtin_popcountll(adj[w][s] & C[s]);
                    }
                }
                int W = w2 / 2;
                if (W < wmin) wmin = W;
                if (W > wmax) wmax = W;
                int A = (adj[a][b >> 6] >> (b & 63)) & 1;
                code[a][b] = code[b][a] = (unsigned short)(2 * W + A);
            }
        }
        for (int a = 0; a < m; a++) {
            int r = 0;
            for (int b = 0; b < m; b++) if (b != a) row[r++] = code[a][b];
            qsort(row, r, sizeof(row[0]), cmp_us);
            u64 h = 0x12345678abcdefULL;
            for (int t = 0; t < r; t++) h = mix(h, row[t]);
            vh[a] = h;
        }
        qsort(vh, m, sizeof(vh[0]), cmp_u64);
        u64 H = 0xfeedfacecafebeefULL;
        for (int a = 0; a < m; a++) H = mix(H, vh[a]);
        printf("%d %016llx %d %d\n", x, (unsigned long long)H, wmin, wmax);
        fflush(stdout);
    }
    return 0;
}
