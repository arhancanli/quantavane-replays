/* |U+U| and |U-U| for a set U of distinct non-negative integers (our replay of the two AlphaEvolve rows of 3a).
   Input on stdin: n, then n integers. Sums u_i + u_j (i <= j) and positive differences u_i - u_j (u_i > u_j) are
   marked in bitsets of sizes 2 max + 1 and max + 1; |U-U| = 2 #(positive differences) + 1. Output: "sums diffs max". */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
static int cmp(const void *a, const void *b) { int64_t x = *(const int64_t *)a, y = *(const int64_t *)b; return (x > y) - (x < y); }
int main(void) {
    long n; if (scanf("%ld", &n) != 1 || n < 1) return 2;
    int64_t *u = malloc(n * sizeof(int64_t));
    for (long i = 0; i < n; i++) if (scanf("%lld", (long long *)&u[i]) != 1 || u[i] < 0) return 2;
    qsort(u, n, sizeof(int64_t), cmp);
    for (long i = 1; i < n; i++) if (u[i] == u[i - 1]) { fprintf(stderr, "repeated element\n"); return 2; }
    int64_t mx = u[n - 1];
    size_t sbits = (size_t)(2 * mx + 1), dbits = (size_t)(mx + 1);
    uint64_t *S = calloc(sbits / 64 + 1, 8), *D = calloc(dbits / 64 + 1, 8);
    if (!S || !D) { fprintf(stderr, "out of memory\n"); return 3; }
    for (long i = 0; i < n; i++) {
        int64_t ui = u[i];
        for (long j = 0; j <= i; j++) {
            uint64_t s = (uint64_t)(ui + u[j]); S[s >> 6] |= 1ULL << (s & 63);
            uint64_t d = (uint64_t)(ui - u[j]); D[d >> 6] |= 1ULL << (d & 63);
        }
    }
    long long ns = 0, nd = 0;
    for (size_t k = 0; k <= sbits / 64; k++) ns += __builtin_popcountll(S[k]);
    D[0] &= ~1ULL;                                  /* d = 0 is counted once, below */
    for (size_t k = 0; k <= dbits / 64; k++) nd += __builtin_popcountll(D[k]);
    printf("%lld %lld %lld\n", ns, 2 * nd + 1, (long long)mx);
    return 0;
}
