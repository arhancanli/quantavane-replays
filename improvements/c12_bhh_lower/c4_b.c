/* c4_b.c -- certified lower bound for Dong and Cao's C_4 (arXiv 2610.02809, Theorem 3.11), program B.
 * Written separately from c4_a.c; same mathematics, same refinement policy (so the two must agree to rounding):
 *   C_4 = (5 pi^{3/2}/32) * int r2 g(G_1(y)) A(y)^{-7/2}, y1 = 0, y2 = (r2, 0), y3 = (xi3, eta3 >= 0), y4 = (xi4, eta4),
 *   cubes of side 1/H in r2, eta3 in [0,2], xi3, xi4, eta4 in [-2,2] (DC (20), (25));
 *   a cube is split into 32 halves while depth < ADEPTH, or while depth < DEPTH and the graphs compatible with its
 *   distance bounds do not all have the same g; a final cube contributes g_min * Amax^{-7/2} * int r2 over the cube.
 * Differences from program A: box distances by clamping (gap / far), g by union-find over all 2^4 deletion sets, the
 * compatible graphs by scanning all 64 edge masks, the rounded-up table index by integer bisection, an explicit work
 * stack instead of recursion, the table stored as a full symmetric cube, pairwise (blocked) summation.
 * Build: cc -O2 -o c4_b c4_b.c -lm     Usage: C4_DEPTH=d C4_ADEPTH=e c4_b TABLE H [claim]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int NT, MT;
static double *TAB3, *TABC, *TABL = NULL;
static double G[64];
static const int PU[6] = {0, 0, 0, 1, 1, 2}, PV[6] = {1, 2, 3, 2, 3, 3};

static int findp(int *par, int x) { while (par[x] != x) x = par[x] = par[par[x]]; return x; }

static int components(int E, int keep) {      /* components of the subgraph induced on 'keep' */
    int par[4] = {0, 1, 2, 3}, n = 0;
    for (int e = 0; e < 6; e++)
        if (((E >> e) & 1) && ((keep >> PU[e]) & 1) && ((keep >> PV[e]) & 1)) {
            int r1 = findp(par, PU[e]), r2 = findp(par, PV[e]);
            if (r1 != r2) par[r1] = r2;
        }
    for (int v = 0; v < 4; v++) if (((keep >> v) & 1) && findp(par, v) == v) n++;
    return n;
}

static double gfun(int E) {
    if (components(E, 15) != 1) return 0.0;
    int deg[4] = {0};
    for (int e = 0; e < 6; e++) if ((E >> e) & 1) { deg[PU[e]]++; deg[PV[e]]++; }
    int twiceD = 0, leaves = 0;
    for (int v = 0; v < 4; v++) { if (deg[v] < 2) twiceD += 2 - deg[v]; if (deg[v] == 1) leaves++; }
    int q = 1, half = leaves / 2 + (leaves % 2);
    if (half > q) q = half;
    for (int Sd = 1; Sd < 16; Sd++) {
        int keep = 15 ^ Sd, removed = 0;
        for (int v = 0; v < 4; v++) removed += (Sd >> v) & 1;
        int val = components(E, keep) - removed;
        if (val > q) q = val;
    }
    return q - twiceD / 2.0;
}

static inline double tab3(int i, int j, int k) {
    if (i < 1) i = 1; if (j < 1) j = 1; if (k < 1) k = 1;
    if (i > MT) i = MT; if (j > MT) j = MT; if (k > MT) k = MT;
    return TAB3[((size_t)i * (MT + 1) + j) * (MT + 1) + k];
}

static long long FU;   /* finest unit: 1/FU */

static inline int kup(long long d2) {        /* least k >= 0 with k^2 FU^2 >= d2 NT^2, by bisection, capped */
    long long lo = 0, hi = MT;               /* if even MT is too small, return MT (the table saturates at 2) */
    long long X = d2 * (long long)NT * NT, F2 = FU * FU;
    if (hi * hi * F2 < X) return MT;
    while (lo < hi) { long long mid = (lo + hi) / 2; if (mid * mid * F2 >= X) hi = mid; else lo = mid + 1; }
    return (int)lo;
}

static inline int kdown(long long d2) {      /* largest k with k^2 FU^2 <= d2 NT^2, by bisection, capped */
    long long X = d2 * (long long)NT * NT, F2 = FU * FU, lo = 0, hi = MT;
    if (hi * hi * F2 <= X) return MT;
    while (lo < hi) { long long mid = (lo + hi + 1) / 2; if (mid * mid * F2 <= X) lo = mid; else hi = mid - 1; }
    return (int)lo;
}

static inline void gapfar(long long alo, long long ahi, long long blo, long long bhi, long long *g, long long *f) {
    long long t1 = blo - ahi, t2 = alo - bhi;
    long long gap = t1 > t2 ? t1 : t2; if (gap < 0) gap = 0;
    long long f1 = bhi - alo, f2 = ahi - blo;
    *g = gap; *f = f1 > f2 ? f1 : f2;
}

typedef struct { long long lo[5]; long long s; int depth; } Cube;

int main(int argc, char **argv) {
    if (argc < 3) return 1;
    FILE *fp = fopen(argv[1], "r"); if (!fp) return 1;
    int H = atoi(argv[2]);
    int DEPTH = getenv("C4_DEPTH") ? atoi(getenv("C4_DEPTH")) : 0;
    int ADEPTH = getenv("C4_ADEPTH") ? atoi(getenv("C4_ADEPTH")) : 0;
    double TAU = getenv("C4_TAU") ? atof(getenv("C4_TAU")) : -1.0;   /* adaptive splitting (decisions only) */
    if (ADEPTH > DEPTH) DEPTH = ADEPTH;
    char t[8], val[64];
    if (fscanf(fp, "%7s %d", t, &NT) != 2) return 2;
    MT = 2 * NT;
    TAB3 = malloc(sizeof(double) * (size_t)(MT + 1) * (MT + 1) * (MT + 1));
    TABC = malloc(sizeof(double) * (MT + 1));
    for (size_t z = 0; z < (size_t)(MT + 1) * (MT + 1) * (MT + 1); z++) TAB3[z] = NAN;
    for (int z = 0; z <= MT; z++) TABC[z] = NAN;
    while (fscanf(fp, "%7s", t) == 1) {
        if (t[0] == 'C') { int k; if (fscanf(fp, "%d %63s", &k, val) != 2) return 2; TABC[k] = nextafter(atof(val), HUGE_VAL); }
        else {
            int i, j, k; if (fscanf(fp, "%d %d %d %63s", &i, &j, &k, val) != 4) return 2;
            double v = nextafter(atof(val), HUGE_VAL);
            int p[3] = {i, j, k};
            static const int perm[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
            for (int q = 0; q < 6; q++) TAB3[((size_t)p[perm[q][0]] * (MT + 1) + p[perm[q][1]]) * (MT + 1) + p[perm[q][2]]] = v;
        }
    }
    fclose(fp);
    for (int i = 1; i <= MT; i++) for (int j = 1; j <= MT; j++) for (int k = 1; k <= MT; k++)
        if (isnan(TAB3[((size_t)i * (MT + 1) + j) * (MT + 1) + k])) { fprintf(stderr, "table hole %d %d %d\n", i, j, k); return 3; }
    for (int z = 0; z <= MT; z++) if (isnan(TABC[z])) { fprintf(stderr, "pair table hole %d\n", z); return 3; }
    double eqv = tab3(NT, NT, NT), ref = 3 * M_PI - 3 * (2 * M_PI / 3 - sqrt(3.0) / 2) + (M_PI - sqrt(3.0)) / 2;
    if (eqv < ref - 1e-12 || eqv > ref + 1e-9) { fprintf(stderr, "table check failed\n"); return 4; }
    if (getenv("C4_LTAB")) {
        FILE *q = fopen(getenv("C4_LTAB"), "r"); int n2;
        if (!q || fscanf(q, "%7s %d", t, &n2) != 2 || n2 != NT) { fprintf(stderr, "bad L table\n"); return 6; }
        TABL = malloc(sizeof(double) * (MT + 1));
        for (int z = 0; z <= MT; z++) TABL[z] = NAN;
        while (fscanf(q, "%7s", t) == 1) { int k; if (fscanf(q, "%d %63s", &k, val) != 2) return 6; TABL[k] = nextafter(atof(val), -HUGE_VAL); }
        fclose(q);
        for (int z = 0; z <= MT; z++) if (isnan(TABL[z])) { fprintf(stderr, "L hole\n"); return 6; }
    }
    for (int E = 0; E < 64; E++) G[E] = gfun(E);
    printf("g: path %.1f star %.1f cycle %.1f K4 %.1f\n", G[(1<<0)|(1<<3)|(1<<5)], G[7], G[(1<<0)|(1<<3)|(1<<5)|(1<<2)], G[63]);

    FU = (long long)H << DEPTH;
    const long long F2 = FU * FU, S0 = 1LL << DEPTH;
    double blocks_total = 0, blocks_comp = 0;
    long long visited = 0, positive = 0, split = 0;
    Cube *stk = malloc(sizeof(Cube) * 4096);
    const int R = 2 * H;
    for (int a = 0; a < R; a++) {
        double row = 0;          /* each row summed in blocks of at most 2^20 terms, blocks summed with compensation */
        double blk = 0; long blkn = 0;
        for (int b = -R; b < R; b++) for (int c = 0; c < R; c++) for (int d = -R; d < R; d++) for (int e = -R; e < R; e++) {
            int sp = 0;
            Cube c0 = {{a * S0, b * S0, c * S0, d * S0, e * S0}, S0, 0};
            stk[sp++] = c0;
            while (sp) {
                Cube Q = stk[--sp];
                visited++;
                long long X[4][2], Y[4][2];         /* boxes of y1..y4 */
                X[0][0] = X[0][1] = 0; Y[0][0] = Y[0][1] = 0;
                X[1][0] = Q.lo[0]; X[1][1] = Q.lo[0] + Q.s; Y[1][0] = Y[1][1] = 0;
                X[2][0] = Q.lo[1]; X[2][1] = Q.lo[1] + Q.s; Y[2][0] = Q.lo[2]; Y[2][1] = Q.lo[2] + Q.s;
                X[3][0] = Q.lo[3]; X[3][1] = Q.lo[3] + Q.s; Y[3][0] = Q.lo[4]; Y[3][1] = Q.lo[4] + Q.s;
                long long dmin2[6], dmax2[6];
                for (int k = 0; k < 6; k++) {
                    long long gx, fx, gy, fy;
                    gapfar(X[PU[k]][0], X[PU[k]][1], X[PV[k]][0], X[PV[k]][1], &gx, &fx);
                    gapfar(Y[PU[k]][0], Y[PU[k]][1], Y[PV[k]][0], Y[PV[k]][1], &gy, &fy);
                    dmin2[k] = gx * gx + gy * gy; dmax2[k] = fx * fx + fy * fy;
                }
                if (dmin2[1] > 4 * F2 || dmin2[2] > 4 * F2) continue;     /* |y3| or |y4| > 2: zero */
                int must = 0, may = 0;
                for (int k = 0; k < 6; k++) { if (dmax2[k] <= F2) must |= 1 << k; if (dmin2[k] <= F2) may |= 1 << k; }
                double glo = 1e9, ghi = -1e9;
                for (int E = 0; E < 64; E++) if ((E & must) == must && (E & ~may) == 0) {
                    if (G[E] < glo) glo = G[E];
                    if (G[E] > ghi) ghi = G[E];
                }
                if (ghi <= 0) continue;
                if (TAU < 0 && (Q.depth < ADEPTH || (glo < ghi && Q.depth < DEPTH))) {
                    split++;
                    long long h = Q.s / 2;
                    for (int m = 0; m < 32; m++) {
                        Cube C; C.s = h; C.depth = Q.depth + 1;
                        for (int z = 0; z < 5; z++) C.lo[z] = Q.lo[z] + (((m >> z) & 1) ? h : 0);
                        stk[sp++] = C;
                    }
                    continue;
                }
                if (TAU < 0 && glo <= 0) continue;
                int K[6];
                for (int k = 0; k < 6; k++) K[k] = kup(dmax2[k]);
                /* pair index helper */
                #define PI_(u, v) ((u) < (v) ? ((u) == 0 ? (v) - 1 : ((u) == 1 ? (v) + 1 : 5)) : ((v) == 0 ? (u) - 1 : ((v) == 1 ? (u) + 1 : 5)))
                double best = HUGE_VAL;
                for (int out = 0; out < 4; out++) {              /* the point left out of the triple */
                    int tri[3], nt = 0;
                    for (int v = 0; v < 4; v++) if (v != out) tri[nt++] = v;
                    double A3v = tab3(K[PI_(tri[0], tri[1])], K[PI_(tri[0], tri[2])], K[PI_(tri[1], tri[2])]);
                    for (int z = 0; z < 3; z++) {
                        double v = nextafter(A3v + TABC[K[PI_(out, tri[z])]], HUGE_VAL);
                        if (v < best) best = v;
                    }
                    if (TABL) for (int z = 0; z < 3; z++) {     /* the left-out disc outside two discs of the triple */
                        int u1 = tri[z], u2 = tri[(z + 1) % 3];
                        double A3o = tab3(K[PI_(out, u1)], K[PI_(out, u2)], K[PI_(u1, u2)]);
                        double v = nextafter(nextafter(A3v + A3o, HUGE_VAL) - TABL[kdown(dmin2[PI_(u1, u2)])], HUGE_VAL);
                        if (v < best) best = v;
                    }
                }
                long long r0 = X[1][0], r1 = X[1][1];
                double wr = (double)(r1 * r1 - r0 * r0) * 0.5;
                double s4 = (double)Q.s * (double)Q.s * (double)Q.s * (double)Q.s;
                double wq = wr * s4 / ((double)F2 * (double)F2 * (double)F2);
                if (TAU >= 0 && Q.depth < DEPTH) {
                    double Alow = M_PI;                           /* heuristic estimate, decisions only */
                    for (int out = 0; out < 4; out++) {
                        int tri[3], nt = 0;
                        for (int v = 0; v < 4; v++) if (v != out) tri[nt++] = v;
                        double v = tab3(kdown(dmin2[PI_(tri[0], tri[1])]), kdown(dmin2[PI_(tri[0], tri[2])]), kdown(dmin2[PI_(tri[1], tri[2])]));
                        if (v > Alow) Alow = v;
                    }
                    double est = wq * (ghi * pow(Alow, -3.5) - (glo > 0 ? glo : 0) * pow(best, -3.5));
                    if (est > TAU) {
                        split++;
                        long long h = Q.s / 2;
                        for (int m = 0; m < 32; m++) {
                            Cube C; C.s = h; C.depth = Q.depth + 1;
                            for (int z = 0; z < 5; z++) C.lo[z] = Q.lo[z] + (((m >> z) & 1) ? h : 0);
                            stk[sp++] = C;
                        }
                        continue;
                    }
                }
                if (glo <= 0) continue;
                double pw = pow(best, -3.5) * (1.0 - 1e-13);
                double term = glo * pw * wq * (1.0 - 1e-15);
                blk += term; blkn++;
                if (blkn == (1L << 20)) { row += blk; blk = 0; blkn = 0; }
                positive++;
            }
        }
        row += blk;
        double y = row - blocks_comp, tt = blocks_total + y;
        blocks_comp = (tt - blocks_total) - y; blocks_total = tt;
        if ((a + 1) % 8 == 0) { fprintf(stderr, "rows %d/%d: %.12e\n", a + 1, R, blocks_total); fflush(stderr); }
    }
    /* summation error: <= (2^20 + rows + ...) * eps relative on positive terms; lower the total by 1e-9 */
    double S = blocks_total * (1.0 - 1e-9);
    double pil = 3.14159265358979;
    double C4 = 5.0 * pil * sqrt(pil) / 32.0 * (1.0 - 1e-15) * S * (1.0 - 1e-15);
    printf("H %d depth %d adepth %d table N %d: %lld cubes visited, %lld positive, %lld split\n", H, DEPTH, ADEPTH, NT, visited, positive, split);
    printf("C_4 >= %.12f   [Dong-Cao: 0.005040573276]\n", C4);
    if (argc > 3) printf("claim C_4 >= %s: %s\n", argv[3], C4 >= atof(argv[3]) ? "HOLDS" : "NOT IMPLIED");
    return 0;
}
