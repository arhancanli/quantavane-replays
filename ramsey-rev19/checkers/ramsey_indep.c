/* Independent clique check of a two-colouring of K_n given as n lines of n characters '0'/'1' (written apart from the
   engine's checkers). It confirms the matrix is symmetric with a zero diagonal and circulant; translations are then
   automorphisms of both colour graphs, so each colour's clique number is 1 + the clique number of vertex 0's
   neighbourhood in that colour, found exactly by branch and bound with greedy-colouring bounds (MCQ).

   cc -O2 -o ramsey_indep ramsey_indep.c && ./ramsey_indep matrix.txt
   Circulant colourings are searched from vertex 0 only; any other symmetric colouring is searched over all vertices. */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 1024
#define W 16

typedef struct { uint64_t w[W]; } bits;

static int n, best;
static long nodes;
static bits adj[MAXN];
static char row[MAXN][MAXN + 1];

static int empty(const bits *s) { for (int i = 0; i < W; i++) if (s->w[i]) return 0; return 1; }
static int low(const bits *s) { for (int i = 0; i < W; i++) if (s->w[i]) return 64 * i + __builtin_ctzll(s->w[i]); return -1; }
static void clr(bits *s, int v) { s->w[v >> 6] &= ~(1ULL << (v & 63)); }

static void expand(bits P, int size) {
    int order[MAXN], col[MAXN], m = 0, c = 0;
    bits U = P;
    nodes++;
    while (!empty(&U)) {  /* greedy colouring of P: each class is an independent set, so a clique uses one vertex per class */
        c++;
        bits Q = U;
        while (!empty(&Q)) {
            int v = low(&Q);
            clr(&Q, v); clr(&U, v);
            for (int k = 0; k < W; k++) Q.w[k] &= ~adj[v].w[k];
            order[m] = v; col[m++] = c;
        }
    }
    for (int i = m - 1; i >= 0; i--) {  /* P holds order[0..i], which needs at most col[i] colours */
        if (size + col[i] <= best) return;
        int v = order[i];
        bits NP;
        for (int k = 0; k < W; k++) NP.w[k] = P.w[k] & adj[v].w[k];
        if (empty(&NP)) { if (size + 1 > best) best = size + 1; }
        else expand(NP, size + 1);
        clr(&P, v);
    }
}

int main(int argc, char **argv) {
    char line[4096];
    FILE *f = argc > 1 ? fopen(argv[1], "r") : NULL;
    if (!f) { fprintf(stderr, "usage: ramsey_indep matrix.txt\n"); return 2; }
    while (fgets(line, sizeof line, f)) {
        size_t len = strcspn(line, "\r\n");
        if (!len) continue;
        if (n >= MAXN || len > MAXN) { fprintf(stderr, "matrix larger than %d\n", MAXN); return 2; }
        memcpy(row[n], line, len); row[n++][len] = 0;
    }
    fclose(f);
    int symmetric = 1, circulant = 1;
    for (int i = 0; i < n; i++) {
        if ((int)strlen(row[i]) != n) { fprintf(stderr, "row %d is not %d long\n", i, n); return 2; }
        if (row[i][i] != '0') symmetric = 0;
    }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            if (row[i][j] != '0' && row[i][j] != '1') { fprintf(stderr, "entry %d,%d is not 0/1\n", i, j); return 2; }
            if (row[i][j] != row[j][i]) symmetric = 0;
            if (row[i][j] != row[0][(j - i + n) % n]) circulant = 0;
        }
    printf("n %d symmetric %d circulant %d\n", n, symmetric, circulant);
    if (!symmetric) return 1;  /* a colouring must be symmetric with a zero diagonal */
    for (int c = 1; c >= 0; c--) {
        int deg = 0;
        memset(adj, 0, sizeof adj);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                if (i != j && row[i][j] == '0' + c) adj[i].w[j >> 6] |= 1ULL << (j & 63);
        for (int j = 0; j < n; j++) deg += (int)(adj[0].w[j >> 6] >> (j & 63) & 1);
        best = 1; nodes = 0;
        if (circulant) expand(adj[0], 1);  /* by symmetry a largest clique may contain vertex 0 */
        else {  /* no symmetry: search every vertex set */
            bits all;
            memset(&all, 0, sizeof all);
            for (int v = 0; v < n; v++) all.w[v >> 6] |= 1ULL << (v & 63);
            expand(all, 0);
        }
        printf("entries '%d': degree%s %d, largest clique %d (%ld branch nodes)\n", c, circulant ? "" : " of vertex 0", deg, best, nodes);
        fflush(stdout);
    }
    return 0;
}
