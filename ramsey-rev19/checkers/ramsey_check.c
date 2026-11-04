/* An independent check of a two-colouring of K_n: does it hold a red K_s or a blue K_t? It is written apart from the
   searches (ramsey_local.c, ramsey_circ.c, ramsey_block.c), shares no code with them and uses another algorithm:
   Bron-Kerbosch with Tomita's pivot. At each node it branches only on the candidates that are not neighbours of the
   pivot, the candidate with the most neighbours among the candidates, and stops when the clique so far plus the
   number of colours a greedy colouring of the candidates needs cannot reach the size sought. No exclusion set is
   kept: a clique of the size sought lies in a maximal clique, and each maximal clique within the candidates holds the
   pivot or a vertex branched on, so the search finds one whenever one exists.

   ramsey_check S T FILE
   FILE holds the n x n adjacency matrix, one line of n characters 0/1 per vertex, 1 red and 0 blue; it must be
   symmetric with a zero diagonal. Prints "n N red_K<S> none|found blue_K<T> none|found nodes X seconds Y", then the
   vertices of each clique found. With OMEGA=1 it also prints "omega_red A omega_blue B", the largest clique of each
   colour. Exit status 0 means neither clique exists, 1 that one does, 2 a bad input. */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAXV 2048
#define WORDS (MAXV / 64)

static int nv, nw;
static uint64_t g[2][MAXV][WORDS];   /* g[0]: red neighbours, g[1]: blue neighbours */
static uint64_t nodes;
static int clique[MAXV];

static int count(const uint64_t *s) { int k = 0; for (int i = 0; i < nw; i++) k += __builtin_popcountll(s[i]); return k; }

/* colours a greedy colouring of S needs, counting no further than enough */
static int greedy_colours(int c, const uint64_t *S, int enough) {
    uint64_t left[WORDS], avail[WORDS];
    int k = 0;
    memcpy(left, S, sizeof(uint64_t) * nw);
    while (k < enough) {
        int empty = 1;
        for (int i = 0; i < nw; i++) if (left[i]) { empty = 0; break; }
        if (empty) break;
        k++;
        memcpy(avail, left, sizeof(uint64_t) * nw);
        for (int i = 0; i < nw; i++)
            while (avail[i]) {
                int x = i * 64 + __builtin_ctzll(avail[i]);
                avail[i] &= avail[i] - 1;
                left[i] &= ~(1ULL << (x & 63));
                for (int j = 0; j < nw; j++) avail[j] &= ~g[c][x][j];
            }
    }
    return k;
}

/* is there a clique of `want` vertices in colour c holding clique[0..size-1] and otherwise inside P? */
static int bk(int c, uint64_t *P, int size, int want) {
    nodes++;
    if (size >= want) return 1;
    int left = count(P);
    if (size + left < want) return 0;
    if (size + greedy_colours(c, P, want - size) < want) return 0;
    int pivot = -1, most = -1;
    for (int i = 0; i < nw; i++)
        for (uint64_t m = P[i]; m; m &= m - 1) {
            int x = i * 64 + __builtin_ctzll(m), k = 0;
            for (int j = 0; j < nw; j++) k += __builtin_popcountll(P[j] & g[c][x][j]);
            if (k > most) { most = k; pivot = x; }
        }
    uint64_t branch[WORDS], next[WORDS];
    for (int i = 0; i < nw; i++) branch[i] = P[i] & ~g[c][pivot][i];
    for (int i = 0; i < nw; i++)
        while (branch[i]) {
            int x = i * 64 + __builtin_ctzll(branch[i]);
            branch[i] &= branch[i] - 1;
            for (int j = 0; j < nw; j++) next[j] = P[j] & g[c][x][j];
            clique[size] = x;
            if (bk(c, next, size + 1, want)) return 1;
            P[i] &= ~(1ULL << (x & 63));
            if (size + count(P) < want) return 0;
        }
    return 0;
}

static int has(int c, int want) {
    uint64_t all[WORDS] = {0};
    for (int x = 0; x < nv; x++) all[x >> 6] |= 1ULL << (x & 63);
    return bk(c, all, 0, want);
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: ramsey_check S T FILE\n"); return 2; }
    int want[2] = {atoi(argv[1]), atoi(argv[2])};
    if (want[0] < 2 || want[1] < 2) { fprintf(stderr, "bad sizes\n"); return 2; }
    FILE *f = fopen(argv[3], "r");
    if (!f) { fprintf(stderr, "cannot read %s\n", argv[3]); return 2; }
    static char row[MAXV + 8];
    static uint8_t a[MAXV][MAXV];
    int width = -1;
    while (fgets(row, sizeof row, f)) {
        int m = (int)strspn(row, "01");
        if (m == 0) continue;
        if (width < 0) width = m;
        if (m != width || nv >= MAXV) { fprintf(stderr, "not a square 0/1 matrix\n"); return 2; }
        for (int y = 0; y < m; y++) a[nv][y] = (uint8_t)(row[y] - '0');
        nv++;
    }
    fclose(f);
    if (nv != width || nv < 1) { fprintf(stderr, "%d rows of %d\n", nv, width); return 2; }
    nw = (nv + 63) / 64;
    for (int x = 0; x < nv; x++) {
        if (a[x][x]) { fprintf(stderr, "1 on the diagonal at %d\n", x); return 2; }
        for (int y = 0; y < nv; y++) {
            if (a[x][y] != a[y][x]) { fprintf(stderr, "not symmetric at %d,%d\n", x, y); return 2; }
            if (x != y) g[a[x][y] ? 0 : 1][x][y >> 6] |= 1ULL << (y & 63);
        }
    }
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    int found[2], verts[2][MAXV];
    for (int c = 0; c < 2; c++) {
        found[c] = has(c, want[c]);
        if (found[c]) memcpy(verts[c], clique, sizeof(int) * want[c]);
    }
    clock_gettime(CLOCK_MONOTONIC, &t1);
    printf("n %d red_K%d %s blue_K%d %s nodes %llu seconds %.2f\n", nv, want[0], found[0] ? "found" : "none", want[1],
           found[1] ? "found" : "none", (unsigned long long)nodes,
           (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) * 1e-9);
    for (int c = 0; c < 2; c++)
        if (found[c]) {
            printf("%s:", c ? "blue" : "red");
            for (int i = 0; i < want[c]; i++) printf(" %d", verts[c][i]);
            printf("\n");
        }
    if (getenv("OMEGA") && atoi(getenv("OMEGA"))) {
        int omega[2];
        for (int c = 0; c < 2; c++) {
            /* without a clique of want[c] the largest is below it: try sizes downwards, else upwards */
            int k = want[c];
            if (found[c]) { while (has(c, k + 1)) k++; }
            else { k--; while (k > 1 && !has(c, k)) k--; }
            omega[c] = k;
        }
        printf("omega_red %d omega_blue %d\n", omega[0], omega[1]);
    }
    return found[0] || found[1];
}
