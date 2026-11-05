/* c3a_hist3.c: the min-cost residue histogram of the block-length-3 controlled-carry construction (3a page, [Syl2026]
 * k = 3), over all Q = q^3 residues, without storing them.
 *
 * M = <G> cap [0, B]; for d in D = M - M, kappa(d) = min{a + b : a, b in M, a - b = d}. For r mod Q,
 * c(r) = min{kappa(d0) + kappa(d1) + kappa(d2) : d_i in D, d0 + d1 q + d2 q^2 = r (mod Q)}. Since D lies in [-B, B] and
 * B < q, writing r = r0 + r1 q + r2 q^2 (digits in [0, q)) every admissible (d0, d1, d2) has the form
 *   d0 = r0 + q g1,  d1 = r1 - g1 + q g2,  d2 = r2 - g2 + q g3   with g1, g2, g3 in {0, -1}
 * (the carries of the subtraction), and every such vector with all d_i in D is admissible. So with
 *   X_g(r0, r1) = min over g1 of kappa(r0 + q g1) + kappa(r1 - g1 + q g),   T_g(r2) = min over g3 of kappa(r2 - g + q g3),
 * c(r) = min(X_0 + T_0(r2), X_-1 + T_-1(r2)). For each row r0 we fold the q pairs (X_0, X_-1) over r1 into distinct
 * pairs with multiplicities, then sweep r2 once per distinct pair.
 *
 * usage: c3a_hist3 B q PART NPARTS g1 g2 ...           rows r0 = PART, PART + NPARTS, ... (summed by c3a_check.py)
 *        c3a_hist3 B q 0 1 g1 g2 ... --plain          the same histogram by the triple loop over (r0, r1, r2), no folding
 *        c3a_hist3 B q 0 1 g1 g2 ... --scatter        by scattering every triple of D^3 into an array of all Q costs,
 *                                                     straight from the definition (small q only)
 * prints "M |M|", "MM |M+M|", "c COST COUNT" lines, "h COST" (the cost of Q/2 for even q, part 0 only),
 * "pairs DISTINCT" (the folded pairs swept) and "unreached COUNT".
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define INF 0x3fffffff

static long B, q;
static int *kap; /* kap[d] for -B <= d <= B */

static inline int K(long d) { return (d < -B || d > B) ? INF : kap[d]; }
static inline int min2(int a, int b) { return a < b ? a : b; }
static inline int add(int a, int b) { return (a >= INF || b >= INF) ? INF : a + b; }

static int cmp64(const void *a, const void *b) {
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return x < y ? -1 : x > y;
}

int main(int argc, char **argv) {
    int plain = argc > 1 && !strcmp(argv[argc - 1], "--plain");
    int scatter = argc > 1 && !strcmp(argv[argc - 1], "--scatter");
    int ng = argc - 5 - plain - scatter;
    if (ng < 1) { fprintf(stderr, "usage: c3a_hist3 B q PART NPARTS g1 ... [--plain|--scatter]\n"); return 2; }
    B = atol(argv[1]); q = atol(argv[2]);
    long part = atol(argv[3]), nparts = atol(argv[4]);
    if (B < 1 || q <= B || nparts < 1 || part < 0 || part >= nparts) { fprintf(stderr, "need 0 < B < q\n"); return 2; }
    long *g = malloc(sizeof(long) * ng);
    for (int i = 0; i < ng; i++) g[i] = atol(argv[5 + i]);
    char *in = calloc(B + 1, 1);
    in[0] = 1;
    for (long n = 1; n <= B; n++)
        for (int i = 0; i < ng; i++)
            if (g[i] <= n && in[n - g[i]]) { in[n] = 1; break; }
    long nm = 0, nmm = 0;
    for (long n = 0; n <= B; n++) nm += in[n];
    char *mm = calloc(2 * B + 1, 1);
    for (long a = 0; a <= B; a++)
        if (in[a])
            for (long b = a; b <= B; b++)
                if (in[b]) mm[a + b] = 1;
    for (long s = 0; s <= 2 * B; s++) nmm += mm[s];
    int *kap0 = malloc(sizeof(int) * (2 * B + 1));
    kap = kap0 + B;
    for (long d = 0; d <= B; d++) {
        kap[d] = INF;
        for (long b = 0; b + d <= B; b++)
            if (in[b] && in[b + d]) { kap[d] = (int)(2 * b + d); break; }
        kap[-d] = kap[d];
    }
    long Q = q * q * q, cmax = 6 * B + 3;
    uint64_t *hist = calloc(cmax + 1, sizeof(uint64_t)), unreached = 0, distinct = 0;
    /* T_g(r2) for g = 0, -1 */
    int *T0 = malloc(sizeof(int) * q), *T1 = malloc(sizeof(int) * q);
    for (long r2 = 0; r2 < q; r2++) {
        T0[r2] = min2(K(r2), K(r2 - q));
        T1[r2] = min2(K(r2 + 1), K(r2 + 1 - q));
    }
    if (scatter) {
        uint32_t *cost = malloc(sizeof(uint32_t) * Q);
        if (!cost) { fprintf(stderr, "no memory for the scatter\n"); return 1; }
        for (long r = 0; r < Q; r++) cost[r] = INF;
        for (long d0 = -B; d0 <= B; d0++) {
            if (kap[d0] == INF) continue;
            for (long d1 = -B; d1 <= B; d1++) {
                if (kap[d1] == INF) continue;
                for (long d2 = -B; d2 <= B; d2++) {
                    if (kap[d2] == INF) continue;
                    long r = ((d0 + d1 * q + d2 * q * q) % Q + Q) % Q;
                    uint32_t c = (uint32_t)(kap[d0] + kap[d1] + kap[d2]);
                    if (c < cost[r]) cost[r] = c;
                }
            }
        }
        for (long r = 0; r < Q; r++) {
            if (cost[r] == INF) unreached++;
            else hist[cost[r]]++;
        }
        if (q % 2 == 0) printf("h %u\n", cost[Q / 2]);
    } else {
        uint64_t *keys = malloc(sizeof(uint64_t) * q);
        for (long r0 = part; r0 < q; r0 += nparts) {
            for (long r1 = 0; r1 < q; r1++) {
                /* g1 = 0: d0 = r0, d1 = r1 + q g;  g1 = -1: d0 = r0 - q, d1 = r1 + 1 + q g */
                int a0 = K(r0), a1 = K(r0 - q);
                int x0 = min2(add(a0, K(r1)), add(a1, K(r1 + 1)));
                int x1 = min2(add(a0, K(r1 - q)), add(a1, K(r1 + 1 - q)));
                if (plain) {
                    for (long r2 = 0; r2 < q; r2++) {
                        int c = min2(add(x0, T0[r2]), add(x1, T1[r2]));
                        if (c >= INF) unreached++;
                        else hist[c]++;
                    }
                } else {
                    keys[r1] = ((uint64_t)(uint32_t)x0 << 32) | (uint32_t)x1;
                }
            }
            if (plain) continue;
            qsort(keys, q, sizeof(uint64_t), cmp64);
            for (long i = 0; i < q;) {
                long j = i;
                while (j < q && keys[j] == keys[i]) j++;
                uint64_t w = (uint64_t)(j - i);
                int x0 = (int)(keys[i] >> 32), x1 = (int)(keys[i] & 0xffffffffu);
                distinct++;
                for (long r2 = 0; r2 < q; r2++) {
                    int c = min2(add(x0, T0[r2]), add(x1, T1[r2]));
                    if (c >= INF) unreached += w;
                    else hist[c] += w;
                }
                i = j;
            }
        }
        if (part == 0 && q % 2 == 0) { /* h = Q/2 has digits (0, 0, q/2) */
            int a0 = K(0), a1 = K(-q);
            int x0 = min2(add(a0, K(0)), add(a1, K(1))), x1 = min2(add(a0, K(-q)), add(a1, K(1 - q)));
            printf("h %d\n", min2(add(x0, T0[q / 2]), add(x1, T1[q / 2])));
        }
    }
    printf("M %ld\nMM %ld\n", nm, nmm);
    for (long c = 0; c <= cmax; c++)
        if (hist[c]) printf("c %ld %llu\n", c, (unsigned long long)hist[c]);
    printf("pairs %llu\nunreached %llu\n", (unsigned long long)distinct, (unsigned long long)unreached);
    return 0;
}
