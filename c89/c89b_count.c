/* c89b_count.c: exact counts behind the monochromatic K_t density of a step graphon on m equal blocks whose red
 * densities are 0, one fraction p, or 1 (PR #307's 256-block binary Cayley graphon for K_5).
 *
 * Input (stdin): m, then the m x m matrix of codes, 0 (red density 0), 1 (red density p) or 2 (red density 1),
 * diagonal included (a pair of vertices in the same block). Repeated blocks are counted, as for a graphon.
 *
 *   c89b_count t root   counts the tuples (0, x_2, ..., x_t): the first vertex in block 0 (valid when the matrix is
 *                       invariant under a group acting transitively on the blocks, which c89b_check.py checks first)
 *   c89b_count t all PART NPARTS   counts every ordered tuple (x_1, ..., x_t), with no symmetry assumed, for the
 *                       first vertices x_1 = PART, PART + NPARTS, ... (the parts are summed by c89b_check.py)
 *
 * A tuple is red-possible when no pair has code 0 and blue-possible when no pair has code 2; its red product is
 * p^k and its blue one (1 - p)^k, k the number of pairs with code 1. Prints "R k count" and "B k count" for every
 * k with a nonzero count; the search extends a tuple one vertex at a time and drops a colour as soon as a pair rules
 * it out, so the counts are exact integers (at most m^t).
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXM 1024
#define MAXT 6

static int m, t, part = 0, nparts = 1;
static unsigned char code[MAXM][MAXM];
static uint64_t red[MAXT * MAXT], blue[MAXT * MAXT];
static int x[MAXT];

static void ext(int v, int red_ok, int blue_ok, int k) {
    if (v == t) {
        if (red_ok) red[k]++;
        if (blue_ok) blue[k]++;
        return;
    }
    for (int y = 0; y < m; y++) {
        int r = red_ok, b = blue_ok, kk = k;
        for (int u = 0; u < v && (r || b); u++) {
            int c = code[x[u]][y];
            if (c == 0) r = 0;
            else if (c == 2) b = 0;
            else kk++;
        }
        if (!r && !b) continue;
        x[v] = y;
        ext(v + 1, r, b, kk);
    }
}

int main(int argc, char **argv) {
    if (!(argc == 3 && !strcmp(argv[2], "root")) && !(argc == 5 && !strcmp(argv[2], "all"))) {
        fprintf(stderr, "usage: c89b_count t root < matrix | c89b_count t all PART NPARTS < matrix\n");
        return 2;
    }
    if (argc == 5) {
        part = atoi(argv[3]);
        nparts = atoi(argv[4]);
        if (nparts < 1 || part < 0 || part >= nparts) { fprintf(stderr, "bad part\n"); return 2; }
    }
    t = atoi(argv[1]);
    if (t < 3 || t > MAXT || scanf("%d", &m) != 1 || m < 1 || m > MAXM) { fprintf(stderr, "bad t or m\n"); return 2; }
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++) {
            int c;
            if (scanf("%d", &c) != 1 || c < 0 || c > 2) { fprintf(stderr, "bad matrix entry\n"); return 2; }
            code[i][j] = (unsigned char)c;
        }
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
            if (code[i][j] != code[j][i]) { fprintf(stderr, "matrix not symmetric\n"); return 2; }
    if (!strcmp(argv[2], "root")) {
        x[0] = 0;
        ext(1, 1, 1, 0);
    } else {
        for (int y = part; y < m; y += nparts) {
            x[0] = y;
            ext(1, 1, 1, 0);
        }
    }
    for (int k = 0; k <= t * (t - 1) / 2; k++) {
        if (red[k]) printf("R %d %llu\n", k, (unsigned long long)red[k]);
        if (blue[k]) printf("B %d %llu\n", k, (unsigned long long)blue[k]);
    }
    return 0;
}
