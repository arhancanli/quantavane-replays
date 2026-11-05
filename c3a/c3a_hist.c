/* c3a_hist.c: the min-cost residue histogram of a block-length-2 controlled-carry construction (3a page, [Syl2026]
 * k = 2), counted residue by residue without storing them.
 *
 * M = <G> cap [0, B]; for d in D = M - M, kappa(d) = min{a + b : a, b in M, a - b = d}. For a residue r mod Q = q^2,
 * c(r) = min{kappa(d0) + kappa(d1) : d0, d1 in D, d0 + d1 q = r (mod Q)}. Every d lies in [-B, B] and B < q, so
 * d0 = r (mod q) leaves at most two choices, rho = r mod q and rho - q, and each fixes d1 modulo q, which again leaves
 * at most two: c(r) is the least of at most four sums.
 *
 * usage: c3a_hist B q PART NPARTS g1 g2 ...            residues r = PART, PART + NPARTS, ... (summed by c3a_check.py)
 *        c3a_hist B q 0 1 g1 g2 ... --scatter         the same histogram by scattering every pair of D x D into an
 *                                                      array of all Q costs (a check for small q only)
 * prints "M |M|", "MM |M+M|", "c COST COUNT" lines, "h COST" (the cost of Q/2, part 0 only), "unreached COUNT".
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define INF 0x7fffffff

int main(int argc, char **argv) {
    int scatter = argc > 1 && !strcmp(argv[argc - 1], "--scatter");
    int ng = argc - 5 - scatter;
    if (ng < 1) { fprintf(stderr, "usage: c3a_hist B q PART NPARTS g1 ... [--scatter]\n"); return 2; }
    long B = atol(argv[1]), q = atol(argv[2]);
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
    int *kap0 = malloc(sizeof(int) * (2 * B + 1)), *kap = kap0 + B; /* kap[d] for -B <= d <= B */
    for (long d = 0; d <= B; d++) {
        kap[d] = INF;
        for (long b = 0; b + d <= B; b++)
            if (in[b] && in[b + d]) { kap[d] = (int)(2 * b + d); break; }
        kap[-d] = kap[d];
    }
    long Q = q * q, cmax = 4 * B + 2;
    uint64_t *hist = calloc(cmax + 1, sizeof(uint64_t)), unreached = 0;
    if (scatter) {
        uint32_t *cost = malloc(sizeof(uint32_t) * Q);
        for (long r = 0; r < Q; r++) cost[r] = INF;
        for (long d0 = -B; d0 <= B; d0++) {
            if (kap[d0] == INF) continue;
            for (long d1 = -B; d1 <= B; d1++) {
                if (kap[d1] == INF) continue;
                long r = ((d0 + d1 * q) % Q + Q) % Q;
                uint32_t c = (uint32_t)(kap[d0] + kap[d1]);
                if (c < cost[r]) cost[r] = c;
            }
        }
        for (long r = 0; r < Q; r++) {
            if (cost[r] == INF) unreached++;
            else hist[cost[r]]++;
        }
        printf("h %u\n", cost[Q / 2]);
    } else {
        for (long r = part; r < Q; r += nparts) {
            long rho = r % q, s0 = r / q;
            int best = INF;
            for (int a = 0; a < 2; a++) {
                long d0 = a ? rho - q : rho, s = a ? (s0 + 1) % q : s0;
                if (d0 < -B || d0 > B || kap[d0] == INF) continue;
                for (int b = 0; b < 2; b++) {
                    long d1 = b ? s - q : s;
                    if (d1 < -B || d1 > B || kap[d1] == INF) continue;
                    int c = kap[d0] + kap[d1];
                    if (c < best) best = c;
                }
            }
            if (best == INF) unreached++;
            else hist[best]++;
            if (r == Q / 2) printf("h %d\n", best);
        }
    }
    printf("M %ld\nMM %ld\n", nm, nmm);
    for (long c = 0; c <= cmax; c++)
        if (hist[c]) printf("c %ld %llu\n", c, (unsigned long long)hist[c]);
    printf("unreached %llu\n", (unsigned long long)unreached);
    return 0;
}
