/* The three lattice sums of c16_trunc_check.py in C, for large cutoffs (our own code; no certificate program was read or
   run). For c_k = q/(n d), n = |k|^2, d = q + p n (mu = p/q), over the disc 0 < x^2 + y^2 <= R^2:
       S0 = sum q/(n d),  S1 = sum |k|^2 c_k^2 = sum q^2/(n d^2),  S2 = sum |k|^4 c_k^2 = sum q^2/d^2.
   Each term is floored and ceiled in integers at the scale 2^E (unsigned 128-bit division), and the disc is covered by
   its 8-fold symmetry (0 <= y <= x, multiplicity 4 on the axes and the diagonal, 8 elsewhere), exactly as the Python
   function sums() does. Output: three lines "lo hi" (decimal integers, the sums times 2^E rounded down and up).
   Refuses inputs where a numerator or denominator could exceed 127 bits.
   usage: c16_trunc_fast R p q E */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

typedef unsigned __int128 u128;

static void put(u128 x) {
    char b[64]; int i = 63; b[i] = 0;
    if (!x) { putchar('0'); return; }
    while (x) { b[--i] = (char)('0' + (int)(x % 10)); x /= 10; }
    fputs(b + i, stdout);
}

static int bits(u128 x) { int k = 0; while (x) { k++; x >>= 1; } return k; }

int main(int argc, char **argv) {
    if (argc != 5) { fprintf(stderr, "usage: c16_trunc_fast R p q E\n"); return 2; }
    long long R = atoll(argv[1]), p = atoll(argv[2]), q = atoll(argv[3]); int E = atoi(argv[4]);
    if (R < 1 || R > 100000 || p < 1 || q < 1 || p > 1000000 || q > 1000000 || E < 1 || E > 110) { fprintf(stderr, "bad input\n"); return 2; }
    u128 R2 = (u128)R * R, dmax = (u128)q + (u128)p * R2;
    u128 T0 = ((u128)1 << E) * (u128)q, T1 = ((u128)1 << E) * (u128)q * (u128)q;
    if (bits(R2) + 2 * bits(dmax) > 126 || bits(T1) > 126) { fprintf(stderr, "would overflow 128 bits\n"); return 2; }
    u128 lo[3] = {0, 0, 0}, hi[3] = {0, 0, 0};
    long long terms = 0;
    for (long long y = 0; y <= R; y++) {
        for (long long x = (y > 1 ? y : 1); x <= R; x++) {
            u128 n = (u128)x * x + (u128)y * y;
            if (n > R2) break;
            u128 m = (y == 0 || y == x) ? 4 : 8;
            u128 d = (u128)q + (u128)p * n;
            u128 den[3] = { n * d, n * d * d, d * d }, num[3] = { T0, T1, T1 };
            for (int i = 0; i < 3; i++) {
                u128 f = num[i] / den[i], r = num[i] - f * den[i];
                lo[i] += m * f;
                hi[i] += m * (f + (r ? 1 : 0));
            }
            terms++;
        }
    }
    for (int i = 0; i < 3; i++) { put(lo[i]); putchar(' '); put(hi[i]); putchar('\n'); }
    fprintf(stderr, "terms (one eighth of the disc): %lld\n", terms);
    return 0;
}
