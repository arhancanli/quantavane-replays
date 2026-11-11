/* Exact independence test for a circulant graph Cay(Z_n, S), written for the 17b replay (no code of the certificate's
   checkers or search was read or run). Input (stdin), produced by c17b_check.py:
       n K
       s_1 ... s_m              (the connection set S, closed under negation, 0 not in S)
       r                        (the number of orbits of the stabilizer of 0 on the non-neighbours of 0)
       then for each of the n vertices v: orb(v)   (-1 for 0 and the neighbours of 0, else 0..r-1)
       then for each orbit i: rep(i)               (a member of orbit i)
   Question: is there an independent set of size K? By vertex-transitivity it may contain 0, so the rest is a clique
   of size K-1 in R(0), the non-neighbours of 0, in the red graph (the complement). Orbital branching: let i be the
   first orbit (in index order) that the clique meets; a stabilizer element maps the clique to one containing rep(i)
   and still avoiding orbits 0..i-1 (they are unions of orbits). So the answer is no exactly when, for every i, the
   red graph on P_i = R(0) ∩ R(rep i) ∩ {v : orb(v) >= i} has no clique of size K-2.
   Each P_i is searched exactly by branch and bound with greedy-colouring bounds (Tomita's MCQ), as a decision: the
   bound starts at K-3, and the search stops at the first clique of size K-2, which is then printed and checked.
   Output: one line per orbit (its size, nodes, seconds, and NONE or the clique found), then RESULT.
   Optional arguments FIRST STEP run only the orbits FIRST, FIRST + STEP, ... (each subproblem still excludes every
   orbit with a smaller index, so the union over FIRST = 0 .. STEP-1 is the whole search). */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAXN 512
#define W (MAXN / 64)

typedef struct { uint64_t w[W]; } bits;

static int n, K, target, found, foundsz, clique[MAXN], best_clique[MAXN];
static long long nodes;
static bits red[MAXN];
static int inS[MAXN];

static int empty(const bits *s) { for (int i = 0; i < W; i++) if (s->w[i]) return 0; return 1; }
static int low(const bits *s) { for (int i = 0; i < W; i++) if (s->w[i]) return 64 * i + __builtin_ctzll(s->w[i]); return -1; }
static void clr(bits *s, int v) { s->w[v >> 6] &= ~(1ULL << (v & 63)); }
static void set(bits *s, int v) { s->w[v >> 6] |= 1ULL << (v & 63); }
static int has(const bits *s, int v) { return (s->w[v >> 6] >> (v & 63)) & 1; }

/* clique of size `size` in clique[0..size-1]; P = common red neighbours not yet decided */
static void expand(bits P, int size) {
    int order[MAXN], col[MAXN], m = 0, c = 0;
    bits U = P;
    if (found) return;
    nodes++;
    while (!empty(&U)) {  /* greedy colouring: each class is independent in red, so a red clique uses one per class */
        c++;
        bits Q = U;
        while (!empty(&Q)) {
            int v = low(&Q);
            clr(&Q, v); clr(&U, v);
            for (int k = 0; k < W; k++) Q.w[k] &= ~red[v].w[k];
            order[m] = v; col[m++] = c;
        }
    }
    for (int i = m - 1; i >= 0 && !found; i--) {
        if (size + col[i] < target) return;  /* cannot reach the target size */
        int v = order[i];
        clique[size] = v;
        if (size + 1 >= target) { found = 1; foundsz = size + 1; memcpy(best_clique, clique, sizeof(int) * (size + 1)); return; }
        bits NP;
        for (int k = 0; k < W; k++) NP.w[k] = P.w[k] & red[v].w[k];
        if (!empty(&NP)) expand(NP, size + 1);
        clr(&P, v);
    }
}

int main(int argc, char **argv) {
    int m, r;
    int first = argc > 1 ? atoi(argv[1]) : 0, step = argc > 2 ? atoi(argv[2]) : 1;   /* run orbits first, first+step, ... */
    if (first < 0 || step < 1) { fprintf(stderr, "bad orbit range\n"); return 2; }
    if (scanf("%d %d", &n, &K) != 2 || n < 3 || n > MAXN || K < 3) { fprintf(stderr, "bad header\n"); return 2; }
    if (scanf("%d", &m) != 1) return 2;
    for (int i = 0; i < m; i++) { int s; if (scanf("%d", &s) != 1 || s <= 0 || s >= n) return 2; inS[s] = 1; }
    for (int s = 1; s < n; s++) if (inS[s] != inS[n - s]) { fprintf(stderr, "S not symmetric\n"); return 2; }
    for (int u = 0; u < n; u++) {
        memset(&red[u], 0, sizeof(bits));
        for (int v = 0; v < n; v++) if (v != u && !inS[((v - u) % n + n) % n]) set(&red[u], v);
    }
    if (scanf("%d", &r) != 1 || r < 1 || r > MAXN) return 2;
    int orb[MAXN], rep[MAXN];
    for (int v = 0; v < n; v++) if (scanf("%d", &orb[v]) != 1) return 2;
    for (int i = 0; i < r; i++) if (scanf("%d", &rep[i]) != 1) return 2;
    /* consistency: the orbits partition R(0) exactly, and each representative lies in its orbit */
    for (int v = 0; v < n; v++) {
        int inR = (v != 0) && has(&red[0], v);
        if (inR != (orb[v] >= 0) || orb[v] >= r) { fprintf(stderr, "orbit labels do not partition R(0)\n"); return 2; }
    }
    for (int i = 0; i < r; i++) if (orb[rep[i]] != i) { fprintf(stderr, "representative not in its orbit\n"); return 2; }
    target = K - 2;
    long long total_nodes = 0;
    int any = 0;
    clock_t t_all = clock();
    for (int i = first; i < r && !any; i += step) {
        bits P;
        memset(&P, 0, sizeof P);
        for (int v = 0; v < n; v++)
            if (orb[v] >= i && v != rep[i] && has(&red[0], v) && has(&red[rep[i]], v)) set(&P, v);
        int sz = 0;
        for (int v = 0; v < n; v++) sz += has(&P, v);
        found = 0; nodes = 0;
        clock_t t0 = clock();
        if (target <= 0) found = 1;
        else if (sz >= target) expand(P, 0);
        double sec = (double)(clock() - t0) / CLOCKS_PER_SEC;
        total_nodes += nodes;
        if (found) {
            any = 1;
            printf("orbit %2d rep %3d |P| = %3d nodes %lld %.1f s: FOUND independent set {0, %d", i, rep[i], sz, nodes, sec, rep[i]);
            for (int j = 0; j < foundsz; j++) printf(", %d", best_clique[j]);
            printf("}\n");
        } else {
            printf("orbit %2d rep %3d |P| = %3d nodes %lld %.1f s: NONE\n", i, rep[i], sz, nodes, sec);
        }
        fflush(stdout);
    }
    double sec = (double)(clock() - t_all) / CLOCKS_PER_SEC;
    printf("RESULT %s: an independent set of size %d %s (nodes %lld, %.1f s)\n", any ? "FOUND" : "NONE", K,
           any ? "exists" : "does not exist", total_nodes, sec);
    return any ? 1 : 0;
}
