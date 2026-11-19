/* F(x,y) = number of K4 through y in the descendant D_x of a two-graph
   = number of 5-sets {x,y,a,b,c} all of whose 10 triples are coherent (so F is symmetric: checked by the caller).
   input (stdin): n, then n lines of n chars, '1' iff S_xy = -1 (diagonal '0').
   output (stdout): n lines of n integers F(x,y) (F(x,x) = 0).  Bitsets, exact integer arithmetic. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#define MAXN 1024
#define W ((MAXN + 63) / 64)
static uint64_t NEG[MAXN][W], NX[MAXN][W], ALL[W];
static int n, nw;
static inline int bit(const uint64_t *s, int i) { return (s[i >> 6] >> (i & 63)) & 1; }
int main(void) {
    if (scanf("%d", &n) != 1 || n > MAXN) return 1;
    nw = (n + 63) / 64;
    char *line = malloc(n + 8);
    for (int i = 0; i < n; i++) {
        if (scanf("%s", line) != 1 || (int)strlen(line) != n) return 2;
        for (int j = 0; j < n; j++) if (line[j] == '1') NEG[i][j >> 6] |= 1ULL << (j & 63);
    }
    for (int j = 0; j < n; j++) ALL[j >> 6] |= 1ULL << (j & 63);
    long long *row = malloc(sizeof(long long) * n);
    for (int x = 0; x < n; x++) {
        const uint64_t *nx = NEG[x];
        for (int a = 0; a < n; a++) {
            int fl = bit(nx, a);
            for (int k = 0; k < nw; k++) NX[a][k] = NEG[a][k] ^ nx[k] ^ (fl ? ALL[k] : 0);
            NX[a][x >> 6] &= ~(1ULL << (x & 63));      /* safety: x is not a vertex of D_x */
            NX[a][a >> 6] &= ~(1ULL << (a & 63));
        }
        for (int v = 0; v < n; v++) {
            row[v] = 0; if (v == x) continue;
            const uint64_t *Nv = NX[v]; long long cnt = 0;
            uint64_t C[W], D[W];
            for (int a = 0; a < n; a++) {
                if (!bit(Nv, a)) continue;
                for (int k = 0; k < nw; k++) C[k] = Nv[k] & NX[a][k];
                /* keep b > a */
                for (int k = 0; k < nw; k++) {
                    if (k < (a >> 6)) C[k] = 0;
                    else if (k == (a >> 6)) C[k] &= ((a & 63) == 63) ? 0 : (~0ULL << ((a & 63) + 1));
                }
                for (int kb = 0; kb < nw; kb++) {
                    uint64_t wb = C[kb];
                    while (wb) {
                        int b = kb * 64 + __builtin_ctzll(wb); wb &= wb - 1;
                        /* c > b, c in C & N(b) */
                        for (int k = 0; k < nw; k++) D[k] = C[k] & NX[b][k];
                        for (int k = 0; k < nw; k++) {
                            if (k < (b >> 6)) continue;
                            uint64_t m = D[k];
                            if (k == (b >> 6)) m &= ((b & 63) == 63) ? 0 : (~0ULL << ((b & 63) + 1));
                            cnt += __builtin_popcountll(m);
                        }
                    }
                }
            }
            row[v] = cnt;
        }
        for (int v = 0; v < n; v++) printf(v ? " %lld" : "%lld", row[v]);
        printf("\n"); fflush(stdout);
    }
    return 0;
}
