/* diff2d.c (Quantavane replays, 3a): |U-U| for a capped digit construction, modulo each prime given.
   U = { sum a_i B^i : a in A^d, sum a_i <= T } with 2 max A < B. Every digit of u - v is a_i - b_i (balanced, no
   carries), and a digit vector t is realised iff the coordinatewise least choice x_i = min{x in A : x - t_i in A},
   y_i = x_i - t_i, has sum x <= T and sum y <= T. So |U-U| is the box sum over [0,T]^2 of the coefficients of
   F(u,v)^d, F = sum over the digit differences t of u^x(t) v^y(t). The terms come on stdin as lines "x y m"; the
   two-dimensional dynamic program runs in place (from high to low indices, so each step reads old values) and the
   program prints "p residue" for each prime p; the caller reconstructs the integer by the Chinese remainder theorem.
   usage: diff2d d T p1 p2 ... < terms */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: diff2d d T p1 p2 ... < terms\n"); return 2; }
    int d = atoi(argv[1]), T = atoi(argv[2]);
    int nt = 0, tx[256], ty[256];
    uint64_t tm[256];
    int x, y; unsigned long long m;
    while (scanf("%d %d %llu", &x, &y, &m) == 3) {
        if (nt >= 256 || x < 0 || y < 0 || x > T || y > T) { fprintf(stderr, "bad term\n"); return 2; }
        tx[nt] = x; ty[nt] = y; tm[nt] = m; nt++;
    }
    if (d < 0 || T < 0) return 2;
    size_t W = (size_t)T + 1;
    uint64_t *S = malloc(W * W * sizeof(uint64_t));
    if (!S) { fprintf(stderr, "no memory\n"); return 2; }
    for (int a = 3; a < argc; a++) {
        uint64_t p = strtoull(argv[a], NULL, 10);
        if (p < 3 || p >= ((uint64_t)1 << 62)) { fprintf(stderr, "bad prime\n"); return 2; }
        for (size_t k = 0; k < W * W; k++) S[k] = 0;
        S[0] = 1;
        for (int step = 0; step < d; step++) {
            for (int i = T; i >= 0; i--) {
                for (int j = T; j >= 0; j--) {
                    unsigned __int128 acc = 0;
                    for (int t = 0; t < nt; t++) {
                        int ii = i - tx[t], jj = j - ty[t];
                        if (ii >= 0 && jj >= 0) acc += (unsigned __int128)S[(size_t)ii * W + jj] * tm[t];
                    }
                    S[(size_t)i * W + j] = (uint64_t)(acc % p);
                }
            }
        }
        unsigned __int128 tot = 0;
        for (size_t k = 0; k < W * W; k++) tot += S[k];
        printf("%llu %llu\n", (unsigned long long)p, (unsigned long long)(tot % p));
        fflush(stdout);
    }
    free(S);
    return 0;
}
