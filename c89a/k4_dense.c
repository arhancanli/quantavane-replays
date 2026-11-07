/* k4_dense.c -- the second program: tr((X D_a)^3) for single roots a, by dense matrix products modulo 2^61 - 1.
 *
 * Quantavane replays, our own code.  For a root a, u = X_a and Y = X D_a restricted to the support S of u (the
 * other rows and columns contribute nothing), so tr((X D_a)^3) = sum_{b,c,d in S} Y_bc Y_cd Y_db.  It uses none of
 * the coarse structure that k4_block.c relies on: no constant blocks, no triangles, no decomposition.  Each row of
 * Y*Y is accumulated in 128-bit integers from reduced entries (products below 2^122; at most 32 are added between
 * reductions).
 *
 * usage: k4_dense MATRIX.bin n den colour root [root ...]     (colour 0 = red X = M, 1 = blue X = den - M)
 * prints "DENSE colour c a <root> <residue mod 2^61-1>" per root.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

typedef unsigned __int128 u128;
static const uint64_t P = (1ULL << 61) - 1;
static inline uint64_t redp(u128 x) {
    while (x >> 61) x = (x >> 61) + (x & P);
    uint64_t r = (uint64_t)x;
    return r >= P ? r - P : r;
}

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: k4_dense MATRIX.bin n den colour root...\n"); return 2; }
    const long n = atol(argv[2]); const uint64_t den = strtoull(argv[3], NULL, 10); const int colour = atoi(argv[4]);
    uint32_t *M = malloc((size_t)n * n * 4);
    FILE *fp = fopen(argv[1], "rb");
    if (!fp || fread(M, 4, (size_t)n * n, fp) != (size_t)n * n) { fprintf(stderr, "cannot read\n"); return 2; }
    fclose(fp);
    long *S = malloc(sizeof(long) * n);
    for (int ri = 5; ri < argc; ri++) {
        const long a = atol(argv[ri]);
        long k = 0;
        for (long b = 0; b < n; b++) {
            uint64_t v = colour ? den - M[a * n + b] : M[a * n + b];
            if (v) S[k++] = b;
        }
        uint64_t *Y = malloc(sizeof(uint64_t) * k * k);
        u128 *acc = malloc(sizeof(u128) * k);
        for (long i = 0; i < k; i++)
            for (long j = 0; j < k; j++) {
                uint64_t x = colour ? den - M[S[i] * n + S[j]] : M[S[i] * n + S[j]];
                uint64_t uj = colour ? den - M[a * n + S[j]] : M[a * n + S[j]];
                Y[i * k + j] = redp((u128)x * uj);
            }
        uint64_t tot = 0;
        for (long i = 0; i < k; i++) {
            for (long j = 0; j < k; j++) acc[j] = 0;
            int cnt = 0;
            for (long l = 0; l < k; l++) {
                uint64_t y = Y[i * k + l];
                if (y) {
                    const uint64_t *row = Y + l * k;
                    for (long j = 0; j < k; j++) acc[j] += (u128)y * row[j];
                    if (++cnt == 32) { for (long j = 0; j < k; j++) acc[j] = redp(acc[j]); cnt = 0; }
                }
            }
            u128 t = 0; int c2 = 0;
            for (long j = 0; j < k; j++) {
                t += (u128)redp(acc[j]) * Y[j * k + i];
                if (++c2 == 32) { t = redp(t); c2 = 0; }
            }
            tot = redp((u128)tot + redp(t));
        }
        printf("DENSE colour %d a %ld %llu\n", colour, a, (unsigned long long)tot);
        fflush(stdout);
        free(Y); free(acc);
    }
    return 0;
}
