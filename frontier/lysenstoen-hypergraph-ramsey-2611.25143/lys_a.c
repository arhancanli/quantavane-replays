/* Program A (C): exact check of the six colourings of Lysenstoen, arXiv 2611.25143.
 * Written from the paper's definitions and the ancillary README's description of the data formats;
 * the ancillary verifiers were neither read nor run.
 *
 *   lys_a tri  <n> <s0> <s1> <colouring.txt>     3-subsets of n points, colour string in lexicographic order.
 *        Exhaustive clique search with NO symmetry: every s0-set with all triples colour 0 and every s1-set with all
 *        triples colour 1 is searched for.
 *   lys_a five <short_form.txt> <s>               5-subsets of Z_83 given by the colours of the 5-sets through 0, 1.
 *        Builds the full colouring with the normalisation (min, second min) -> (0, 1), checks that all 20 ordered
 *        normalisations of every 5-set give the same colour (so the colouring is well defined and affine-invariant),
 *        then searches EVERY s-set (no symmetry) for one with all 5-subsets of one colour.
 *   lys_a proj <q> <table.txt> <s0> <s1> <pgl|psl> 4-subsets of the projective line over F_q (points 0..q-1, q = inf).
 *        pgl: colour = table[x], x the cross-ratio of the sorted 4-set; psl: table[x] if delta is a square, else
 *        table[1/x]. Checks: (pgl) the table is constant on the six-element cross-ratio classes; (psl, and pgl when
 *        q < 300) invariance of the colour of EVERY 4-subset under x -> x+1 and x -> -1/x (and x -> r x, r a
 *        primitive root, for pgl); the orbit of {0, 1, inf} under these maps is ALL 3-subsets (breadth-first search),
 *        so every set can be moved to one through 0, 1, inf. Then every s0-set and s1-set through {0, 1, inf}.
 * Output lines "RESULT ... found=0" mean no forbidden set. Exit status 0 iff nothing forbidden was found and all
 * structural checks passed.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static char *readstr(const char *path, long *len) {
    FILE *f = fopen(path, "rb");
    if (!f) { perror(path); exit(2); }
    fseek(f, 0, SEEK_END); long n = ftell(f); fseek(f, 0, SEEK_SET);
    char *s = malloc(n + 1); long k = 0; int c;
    while ((c = fgetc(f)) != EOF) if (c == '0' || c == '1') s[k++] = (char)(c - '0');
    fclose(f); s[k] = 0; *len = k; return s;
}

/* ------------------------------------------------------------------ triples */
#define MAXN 128
static unsigned char C3[MAXN][MAXN][MAXN];
typedef struct { uint64_t w[2]; } bs;
static bs NB[2][MAXN][MAXN];
static long long found3, nodes3;
static int tn;

static void rec3(int col, int *S, int k, int need, bs cand) {
    nodes3++;
    if (k == need) { found3++; if (found3 <= 3) { printf("  FOUND colour %d:", col); for (int i = 0; i < k; i++) printf(" %d", S[i]); printf("\n"); } return; }
    for (int w = 0; w < 2; w++) {
        uint64_t m = cand.w[w];
        while (m) {
            int b = __builtin_ctzll(m); m &= m - 1; int p = w * 64 + b;
            bs nc; nc.w[0] = cand.w[0]; nc.w[1] = cand.w[1];
            /* keep only points > p */
            if (p < 64) { nc.w[0] &= (p == 63) ? 0 : (~0ULL << (p + 1)); }
            else { nc.w[0] = 0; nc.w[1] &= (p - 64 == 63) ? 0 : (~0ULL << (p - 64 + 1)); }
            for (int i = 0; i < k; i++) for (int j = i + 1; j < k; j++) { /* nothing: pairs already enforced */ }
            /* adding p: every triple {S_i, S_j, p} must be colour col (enforced by cand); for later points l,
               need {S_i, p, l} colour col for all i: intersect with NB[col][S_i][p] */
            for (int i = 0; i < k; i++) { nc.w[0] &= NB[col][S[i]][p].w[0]; nc.w[1] &= NB[col][S[i]][p].w[1]; }
            S[k] = p;
            if (k + 1 == need) { nodes3++; found3++; if (found3 <= 3) { printf("  FOUND colour %d:", col); for (int i = 0; i <= k; i++) printf(" %d", S[i]); printf("\n"); } continue; }
            if (k + 1 >= 2 && !(nc.w[0] | nc.w[1])) continue;
            rec3(col, S, k + 1, need, nc);
        }
    }
}

static int do_tri(int n, int s0, int s1, const char *path) {
    long len; char *s = readstr(path, &len);
    long expect = (long)n * (n - 1) * (n - 2) / 6;
    if (len != expect) { printf("RESULT tri n=%d BAD LENGTH %ld != %ld\n", n, len, expect); return 1; }
    long idx = 0, c0 = 0;
    for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) for (int c = b + 1; c < n; c++) {
        unsigned char v = (unsigned char)s[idx++]; c0 += (v == 0);
        C3[a][b][c] = C3[a][c][b] = C3[b][a][c] = C3[b][c][a] = C3[c][a][b] = C3[c][b][a] = v;
    }
    memset(NB, 0, sizeof NB);
    for (int col = 0; col < 2; col++) for (int a = 0; a < n; a++) for (int b = 0; b < n; b++) if (a != b)
        for (int c = 0; c < n; c++) if (c != a && c != b && C3[a][b][c] == col) NB[col][a][b].w[c >> 6] |= 1ULL << (c & 63);
    tn = n; int bad = 0;
    printf("tri n=%d triples=%ld colour0=%ld colour1=%ld\n", n, len, c0, len - c0);
    int sz[2] = {s0, s1};
    for (int col = 0; col < 2; col++) {
        found3 = 0; nodes3 = 0; int S[16];
        /* start from every pair (a<b): candidates = points > b with {a,b,c} colour col */
        for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) {
            bs cand = NB[col][a][b];
            if (b < 63) cand.w[0] &= ~0ULL << (b + 1); else if (b == 63) cand.w[0] = 0; else { cand.w[0] = 0; cand.w[1] &= (b - 64 == 63) ? 0 : (~0ULL << (b - 64 + 1)); }
            S[0] = a; S[1] = b;
            if (sz[col] == 2) { found3++; continue; }
            rec3(col, S, 2, sz[col], cand);
        }
        printf("RESULT tri n=%d colour=%d size=%d found=%lld (search nodes %lld)\n", n, col, sz[col], found3, nodes3);
        if (found3) bad = 1;
    }
    return bad;
}

/* ------------------------------------------------------------------ 5-subsets of Z_83 */
static long long binom[100][8];
static int MODQ = 83;
static int inv83[83];
static long PRE[5][90]; /* PRE[i][x] = sum_{v<x} C(82-v, 4-i) */
static long rank5(int *x) { /* lexicographic rank of a sorted 5-set among combinations(range(83),5) */
    long r = 0; int prev = -1;
    for (int i = 0; i < 5; i++) { r += PRE[i][x[i]] - PRE[i][prev + 1]; prev = x[i]; }
    return r;
}
static long rank3from2(int x, int y, int z) { /* rank of {0,1,x,y,z} among 5-sets through 0,1 = rank of (x,y,z) in combinations(range(2,83),3) */
    int a[5] = {0, 1, x, y, z}; return rank5(a); /* the 5-sets through 0 and 1 come first in lex order */
}
static void sort3(int *a) { int t; if (a[0] > a[1]) { t = a[0]; a[0] = a[1]; a[1] = t; } if (a[1] > a[2]) { t = a[1]; a[1] = a[2]; a[2] = t; } if (a[0] > a[1]) { t = a[0]; a[0] = a[1]; a[1] = t; } }
static unsigned char *FULL5;
static int norm_colour(const char *sh, int *S, int u, int v) { /* map S[u] -> 0, S[v] -> 1 */
    int d = ((S[v] - S[u]) % 83 + 83) % 83; int iv = inv83[d]; int img[3], k = 0;
    for (int i = 0; i < 5; i++) if (i != u && i != v) img[k++] = (int)((((long)(S[i] - S[u]) % 83 + 83) % 83) * iv % 83);
    sort3(img);
    return sh[rank3from2(img[0], img[1], img[2])];
}
static int do_five(const char *path, int s) {
    long len; char *sh = readstr(path, &len);
    for (int i = 0; i < 100; i++) for (int j = 0; j < 8; j++) binom[i][j] = (j == 0) ? 1 : (i == 0 ? 0 : binom[i - 1][j - 1] + binom[i - 1][j]);
    for (int i = 0; i < 5; i++) { PRE[i][0] = 0; for (int x = 0; x < 89; x++) PRE[i][x + 1] = PRE[i][x] + (x <= 82 ? binom[82 - x][4 - i] : 0); }
    if (len != binom[81][3]) { printf("RESULT five BAD SHORT LENGTH %ld\n", len); return 1; }
    for (int a = 1; a < 83; a++) for (int b = 1; b < 83; b++) if (a * b % 83 == 1) inv83[a] = b;
    long total = binom[83][5];
    FULL5 = malloc(total);
    long idx = 0, incons = 0, c0 = 0;
    int S[5];
    for (S[0] = 0; S[0] < 83; S[0]++) for (S[1] = S[0] + 1; S[1] < 83; S[1]++) for (S[2] = S[1] + 1; S[2] < 83; S[2]++)
    for (S[3] = S[2] + 1; S[3] < 83; S[3]++) for (S[4] = S[3] + 1; S[4] < 83; S[4]++) {
        int c = norm_colour(sh, S, 0, 1);
        for (int u = 0; u < 5; u++) for (int v = 0; v < 5; v++) if (u != v && norm_colour(sh, S, u, v) != c) incons++;
        FULL5[idx++] = (unsigned char)c; c0 += (c == 0);
    }
    printf("five n=83 5-subsets=%ld colour0=%ld inconsistent normalisations=%ld\n", idx, c0, incons);
    /* the first 85320 entries (sets through 0,1) must equal the short form itself */
    long mism = 0; for (long i = 0; i < len; i++) if (FULL5[i] != (unsigned char)sh[i]) mism++;
    printf("five: full colouring restricted to sets through 0,1 differs from the short form in %ld places\n", mism);
    long long found[2] = {0, 0}; long long checked = 0;
    int T[6];
    for (T[0] = 0; T[0] < 83; T[0]++) for (T[1] = T[0] + 1; T[1] < 83; T[1]++) for (T[2] = T[1] + 1; T[2] < 83; T[2]++)
    for (T[3] = T[2] + 1; T[3] < 83; T[3]++) for (T[4] = T[3] + 1; T[4] < 83; T[4]++) {
        int x[5] = {T[0], T[1], T[2], T[3], T[4]}; int c = FULL5[rank5(x)];
        for (T[5] = T[4] + 1; T[5] < 83; T[5]++) {
            checked++; int ok = 1;
            for (int omit = 0; omit < 5 && ok; omit++) { int y[5], k = 0; for (int i = 0; i < 6; i++) if (i != omit) y[k++] = T[i]; if (FULL5[rank5(y)] != c) ok = 0; }
            if (ok) { found[c]++; if (found[c] <= 3) printf("  FOUND colour %d: %d %d %d %d %d %d\n", c, T[0], T[1], T[2], T[3], T[4], T[5]); }
        }
    }
    printf("RESULT five n=83 size=%d 6-sets examined=%lld found colour0=%lld colour1=%lld\n", s, checked, found[0], found[1]);
    return (incons || mism || found[0] || found[1]) ? 1 : 0;
}

/* ------------------------------------------------------------------ 4-subsets of the projective line */
static int Q, INF; static int *INV, *ISSQ; static char *TAB; static int PSL;
static long md(long a) { a %= Q; return a < 0 ? a + Q : a; }
static int colour4(int p0, int p1, int p2, int p3) {
    int a[4] = {p0, p1, p2, p3}, t;
    for (int i = 0; i < 4; i++) for (int j = i + 1; j < 4; j++) if (a[j] < a[i]) { t = a[i]; a[i] = a[j]; a[j] = t; }
    long x;
    if (a[3] == INF) x = md((long)(a[0] - a[2])) * INV[md(a[1] - a[2])] % Q;
    else x = md(md(a[0] - a[2]) * md(a[1] - a[3])) * INV[md(md(a[0] - a[3]) * md(a[1] - a[2]))] % Q;
    if (!PSL) return TAB[x];
    long delta = md(md(a[1] - a[0]) * md(a[2] - a[0]) % Q * md(a[2] - a[1]));
    return ISSQ[delta] ? TAB[x] : TAB[INV[x]];
}
static int g_shift(int x) { return x == INF ? INF : (x + 1) % Q; }
static int g_inv(int x) { if (x == INF) return 0; if (x == 0) return INF; return (int)md(-(long)INV[x]); }
static int RPRIM;
static int g_mul(int x) { return x == INF ? INF : (int)((long)x * RPRIM % Q); }
static long long fnd; static int SZ, COL; static int CH[16];
static void recp(int k, int start) {
    if (k == SZ) { fnd++; if (fnd <= 3) { printf("  FOUND colour %d:", COL); for (int i = 0; i < k; i++) printf(" %d", CH[i]); printf("\n"); } return; }
    for (int p = start; p < Q; p++) {
        int ok = 1;
        for (int i = 0; i < k && ok; i++) for (int j = i + 1; j < k && ok; j++) for (int l = j + 1; l < k && ok; l++)
            if (colour4(CH[i], CH[j], CH[l], p) != COL) ok = 0;
        if (!ok) continue;
        CH[k] = p; recp(k + 1, p + 1);
    }
}
static long tri_index(int a, int b, int c, int n) { /* a<b<c, simple 3D index */ return ((long)a * n + b) * n + c; }
static int do_proj(int q, const char *path, int s0, int s1, const char *kind) {
    Q = q; INF = q; PSL = (strcmp(kind, "psl") == 0); int n = q + 1, bad = 0;
    long len; TAB = readstr(path, &len);
    if (len != q) { printf("RESULT proj BAD TABLE LENGTH %ld\n", len); return 1; }
    for (int d = 2; d < q; d++) if (q % d == 0) { printf("q not prime\n"); return 1; }
    INV = calloc(q, sizeof(int)); ISSQ = calloc(q, sizeof(int));
    for (int a = 1; a < q; a++) for (int b = 1; b < q; b++) if ((long)a * b % q == 1) { INV[a] = b; break; }
    for (int a = 1; a < q; a++) ISSQ[(long)a * a % q] = 1;
    for (RPRIM = 2; RPRIM < q; RPRIM++) { long x = 1; int ord = 0; do { x = x * RPRIM % q; ord++; } while (x != 1); if (ord == q - 1) break; }
    if (!PSL) {
        long viol = 0;
        for (int x = 2; x < q; x++) {
            long cls[6] = {x, INV[x], md(1 - x), INV[md(1 - x)], md((long)x * INV[md(x - 1)]), md((long)md(x - 1) * INV[x])};
            for (int i = 0; i < 6; i++) if (TAB[cls[i]] != TAB[x]) viol++;
        }
        printf("proj q=%d: table constant on cross-ratio classes: %s (%ld violations)\n", q, viol ? "NO" : "yes", viol);
        if (viol) bad = 1;
    }
    if (PSL || q < 300) {
        long viol = 0, cnt = 0;
        for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) for (int c = b + 1; c < n; c++) for (int d = c + 1; d < n; d++) {
            int col = colour4(a, b, c, d); cnt++;
            if (colour4(g_shift(a), g_shift(b), g_shift(c), g_shift(d)) != col) viol++;
            if (colour4(g_inv(a), g_inv(b), g_inv(c), g_inv(d)) != col) viol++;
            if (!PSL && colour4(g_mul(a), g_mul(b), g_mul(c), g_mul(d)) != col) viol++;
        }
        printf("proj q=%d: invariance of all %ld 4-subsets under x+1, -1/x%s: %ld violations\n", q, cnt, PSL ? "" : ", r*x", viol);
        if (viol) bad = 1;
    }
    /* orbit of {0,1,inf} on 3-subsets */
    {
        long N3 = (long)n * n * n; unsigned char *seen = calloc(N3 / 8 + 1, 1);
        int *qa = malloc(sizeof(int) * ((long)n * (n - 1) * (n - 2) / 6 + 10)); long head = 0, tail = 0, cnt = 0;
#define SEEN(i) (seen[(i) >> 3] & (1 << ((i) & 7)))
#define MARK(i) (seen[(i) >> 3] |= (1 << ((i) & 7)))
        int st[3] = {0, 1, INF}; long id = tri_index(0, 1, INF, n); MARK(id); qa[0] = (int)id; tail = 1; cnt = 1;
        while (head < tail) {
            long pk = qa[head++]; int z = (int)(pk % n), y = (int)((pk / n) % n), x = (int)(pk / ((long)n * n));
            for (int g = 0; g < (PSL ? 2 : 3); g++) {
                int u[3];
                for (int i = 0; i < 3; i++) { int v = (i == 0 ? x : i == 1 ? y : z); u[i] = g == 0 ? g_shift(v) : g == 1 ? g_inv(v) : g_mul(v); }
                int t; if (u[0] > u[1]) { t = u[0]; u[0] = u[1]; u[1] = t; } if (u[1] > u[2]) { t = u[1]; u[1] = u[2]; u[2] = t; } if (u[0] > u[1]) { t = u[0]; u[0] = u[1]; u[1] = t; }
                long k = tri_index(u[0], u[1], u[2], n);
                if (!SEEN(k)) { MARK(k); qa[tail++] = (int)k; cnt++; }
            }
        }
        long all3 = (long)n * (n - 1) * (n - 2) / 6;
        printf("proj q=%d: orbit of {0,1,inf} has %ld of %ld 3-subsets\n", q, cnt, all3);
        if (cnt != all3) bad = 1;
        (void)st;
    }
    int sz[2] = {s0, s1};
    for (COL = 0; COL < 2; COL++) {
        SZ = sz[COL]; fnd = 0;
        CH[0] = 0; CH[1] = 1; CH[2] = INF;
        /* sets through 0, 1, inf: the other points from 2..q-1 */
        recp(3, 2);
        printf("RESULT proj q=%d colour=%d size=%d sets through {0,1,inf} found=%lld\n", q, COL, SZ, fnd);
        if (fnd) bad = 1;
    }
    return bad;
}

int main(int argc, char **argv) {
    if (argc >= 6 && !strcmp(argv[1], "tri")) return do_tri(atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), argv[5]);
    if (argc >= 4 && !strcmp(argv[1], "five")) return do_five(argv[2], atoi(argv[3]));
    if (argc >= 7 && !strcmp(argv[1], "proj")) return do_proj(atoi(argv[2]), argv[3], atoi(argv[4]), atoi(argv[5]), argv[6]);
    fprintf(stderr, "usage: see header\n"); return 2;
}
