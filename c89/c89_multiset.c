/* c89_multiset.c: exact counts behind the monochromatic K_t density of a weight-class Cayley step graphon on Z_2^n,
 * by multisets of column types (method A; c89_refine.c is method B and shares no code with this file).
 *
 * The graphon has 2^n equal blocks indexed by x in Z_2^n; the red density between blocks x and y depends only on the
 * Hamming weight w = |x xor y| (w = 0 inside a block). Each weight has a class: R (red density 1), B (red density 0)
 * or a fractional index k in 0..3 (red density f_k strictly between 0 and 1). By translation invariance the first of
 * the t vertices is fixed at 0; the other t-1 vertices are any (2^n)^(t-1) tuples. A coordinate's column is the t-1
 * bits of those vertices there (T = 2^(t-1) column types), the distance of a pair is the number of coordinates whose
 * column separates it, and a multiset of n columns stands for n!/prod(m!) tuples.
 *
 * For every tuple, red is possible only if no pair has a B weight and blue only if no pair has an R weight; the red
 * product is then prod_k f_k^(e_k) and the blue one prod_k (1 - f_k)^(e_k), where e_k counts the pairs of class k.
 * This program prints the exact number of tuples behind each (colour, e_0..e_3) and the total, which must be
 * 2^(n(t-1)). The rational density is formed from these counts by c89_check.py.
 *
 * usage: c89_multiset t n CLASSES PART NPARTS   (CLASSES: n+1 comma-separated codes R, B or 0..3, for w = 0..n)
 * The search is split by the counts of the first three column types; part PART of NPARTS takes every NPARTS-th split.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXT 6
#define MAXE 15
#define MAXN 16
#define CELLS 65536 /* 16^4: four fractional indices, each exponent 0..15 */

static int t, n, T, E, part, nparts;
static int cls[MAXN + 1];          /* -1 = R, -2 = B, k >= 0 fractional index */
static int emask[1 << (MAXT - 1)]; /* bit e set when column type separates pair e */
static uint64_t binom[MAXN + 1][MAXN + 1];
static uint64_t red[CELLS], blue[CELLS], total;
static int dist[MAXE];
static long split_no;

static void leaf(uint64_t w) {
    int red_ok = 1, blue_ok = 1, idx = 0;
    for (int e = 0; e < E; e++) {
        int c = cls[dist[e]];
        if (c == -2) red_ok = 0;
        else if (c == -1) blue_ok = 0;
        else idx += 1 << (4 * c);
    }
    if (red_ok) red[idx] += w;
    if (blue_ok) blue[idx] += w;
    total += w;
}

static void dfs(int tau, int rem, uint64_t w) {
    if (tau == 3 && (split_no++ % nparts) != part) return;
    if (tau == T - 1) { /* the last type takes every remaining coordinate */
        for (int e = 0; e < E; e++) if (emask[tau] >> e & 1) dist[e] += rem;
        leaf(w);
        for (int e = 0; e < E; e++) if (emask[tau] >> e & 1) dist[e] -= rem;
        return;
    }
    for (int m = 0; m <= rem; m++) {
        if (m) for (int e = 0; e < E; e++) if (emask[tau] >> e & 1) dist[e] += 1;
        dfs(tau + 1, rem - m, w * binom[rem][m]);
    }
    for (int e = 0; e < E; e++) if (emask[tau] >> e & 1) dist[e] -= rem;
}

int main(int argc, char **argv) {
    if (argc != 6) { fprintf(stderr, "usage: c89_multiset t n CLASSES PART NPARTS\n"); return 2; }
    t = atoi(argv[1]); n = atoi(argv[2]); part = atoi(argv[4]); nparts = atoi(argv[5]);
    if (t < 4 || t > MAXT || n < 1 || n > MAXN || nparts < 1 || part < 0 || part >= nparts) {
        fprintf(stderr, "bad t, n or part\n"); return 2;
    }
    T = 1 << (t - 1); E = t * (t - 1) / 2;
    char *s = strdup(argv[3]), *tok = strtok(s, ",");
    for (int w = 0; w <= n; w++) {
        if (!tok) { fprintf(stderr, "need %d classes\n", n + 1); return 2; }
        if (!strcmp(tok, "R")) cls[w] = -1;
        else if (!strcmp(tok, "B")) cls[w] = -2;
        else if (tok[0] >= '0' && tok[0] <= '3' && !tok[1]) cls[w] = tok[0] - '0';
        else { fprintf(stderr, "bad class %s\n", tok); return 2; }
        tok = strtok(NULL, ",");
    }
    if (tok) { fprintf(stderr, "too many classes\n"); return 2; }
    for (int a = 0; a <= MAXN; a++) {
        binom[a][0] = 1;
        for (int b = 1; b <= a; b++) binom[a][b] = binom[a - 1][b - 1] + (b <= a - 1 ? binom[a - 1][b] : 0);
    }
    /* pairs (a, b), a < b, of the t vertices; vertex 0 is the fixed origin, vertex i >= 1 is bit i-1 of a type */
    int e = 0;
    for (int a = 0; a < t; a++)
        for (int b = a + 1; b < t; b++, e++)
            for (int tau = 0; tau < T; tau++) {
                int ba = a ? (tau >> (a - 1)) & 1 : 0, bb = (tau >> (b - 1)) & 1;
                if (ba != bb) emask[tau] |= 1 << e;
            }
    dfs(0, n, 1);
    for (int i = 0; i < CELLS; i++) {
        if (red[i]) printf("R %d %d %d %d %llu\n", i & 15, i >> 4 & 15, i >> 8 & 15, i >> 12 & 15,
                           (unsigned long long)red[i]);
        if (blue[i]) printf("B %d %d %d %d %llu\n", i & 15, i >> 4 & 15, i >> 8 & 15, i >> 12 & 15,
                            (unsigned long long)blue[i]);
    }
    printf("TOTAL %llu\n", (unsigned long long)total);
    return 0;
}
