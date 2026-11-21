/* c4_a.c -- certified lower bound for Dong and Cao's C_4 (arXiv 2610.02809, Theorem 3.11), program A.
 *
 * C_4 = (2 pi Gamma(7/2) / 4!) * int r2 g(G_1(y)) A(y)^{-7/2} dr2 dxi3 deta3 dxi4 deta4   (DC (20)/(25), eta3 >= 0)
 *     = (5 pi^{3/2} / 32) * int ...,
 * y1 = 0, y2 = (r2, 0), y3 = (xi3, eta3), y4 = (xi4, eta4); the integrand vanishes unless every |y_j| <= 2 (DC's
 * claim, k - 2 = 2), so the domain is r2, eta3 in [0, 2], xi3, xi4, eta4 in [-2, 2], cut into cubes of side 1/H.
 * On each cube Q (DC's Lemma 3.15, with a sharper area bound):
 *   - pair distances: exact integer bounds d2min, d2max of |y_i - y_j|^2 in units of 1/H^2;
 *   - edges: present throughout if d2max <= H^2, absent throughout if d2min > H^2, else undetermined;
 *     g_min = min of g over every graph between the two (g = q - D for connected graphs, 0 otherwise; DC Def. 3.8);
 *   - area: A(y) <= A3(i, j, k) + c(d_lm) for every triple {i, j, k} and every attachment of the fourth point l to a
 *     point m of the triple (the fourth disc adds at most its part outside disc m), where A3 is an UPPER bound of the
 *     union of three unit discs at the rounded-up distances (planar Kneser-Poulsen theorem, Bezdek-Connelly 2002:
 *     the union area does not decrease when no distance decreases) and c(d) = pi - lens(d) is increasing; both come
 *     from a rigorous table (c4_tables.py, rounded up), distances rounded UP to the table grid 1/N;
 *   - contribution >= g_min * Amax^{-7/2} * int_Q r2 = g_min * Amax^{-7/2} * (2a+1) / (2 H^6).
 * Floating point: table values are read with strtod and moved one ulp up; each A bound is a sum of two such values
 * moved one ulp up; pow() is lowered by a relative 1e-13; every product is lowered by a relative 1e-15; the sum uses
 * Neumaier compensation and the final total is lowered by a relative 1e-12 (far above every rounding error made).
 *
 * Build: cc -O2 -o c4_a c4_a.c -lm      Usage: c4_a TABLE H [claim]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

static int N, M;                      /* table grid 1/N, sides up to M = 2N */
static double *A3;                     /* A3[(i*(M+1)+j)*(M+1)+k], i <= j <= k */
static double *C1;
static double *L2 = NULL;                 /* optional: LOWER bounds of the union of two discs at k/N */
static double gtab[64];                /* g of each edge mask on 4 vertices */
static double gmin_tab[64][64];        /* [present][undetermined] */
static double gmax_tab[64][64];
static int DIAG = 0;                   /* diagnostic: use g_max (NOT a bound) */

static const int EU[6] = {0, 0, 0, 1, 1, 2}, EV[6] = {1, 2, 3, 2, 3, 3};

static int edge_index(int u, int v) {
    if (u > v) { int t = u; u = v; v = t; }
    for (int e = 0; e < 6; e++) if (EU[e] == u && EV[e] == v) return e;
    return -1;
}

/* number of connected components of the graph induced on vertex set vs (bit mask) with edge mask E */
static int cc_count(int E, int vs) {
    int seen = 0, comps = 0;
    for (int s = 0; s < 4; s++) {
        if (!((vs >> s) & 1) || ((seen >> s) & 1)) continue;
        comps++;
        int stack = 1 << s; seen |= 1 << s;
        while (stack) {
            int u = __builtin_ctz(stack); stack &= stack - 1;
            for (int v = 0; v < 4; v++) {
                if (!((vs >> v) & 1) || ((seen >> v) & 1) || v == u) continue;
                if ((E >> edge_index(u, v)) & 1) { seen |= 1 << v; stack |= 1 << v; }
            }
        }
    }
    return comps;
}

static double g_of(int E) {
    if (cc_count(E, 15) != 1) return 0.0;
    int deg[4] = {0, 0, 0, 0};
    for (int e = 0; e < 6; e++) if ((E >> e) & 1) { deg[EU[e]]++; deg[EV[e]]++; }
    double D = 0; int n1 = 0;
    for (int v = 0; v < 4; v++) { if (deg[v] < 2) D += 0.5 * (2 - deg[v]); if (deg[v] == 1) n1++; }
    int q = 1;
    if ((n1 + 1) / 2 > q) q = (n1 + 1) / 2;
    for (int S = 0; S < 16; S++) {
        int rest = 15 & ~S, sz = __builtin_popcount(S);
        int val = (rest ? cc_count(E, rest) : 0) - sz;
        if (val > q) q = val;
    }
    return (double)q - D;
}

static inline double A3get(int i, int j, int k) {     /* any order, indices clamped to [1, M] */
    if (i < 1) i = 1; if (j < 1) j = 1; if (k < 1) k = 1;
    if (i > M) i = M; if (j > M) j = M; if (k > M) k = M;
    int t;
    if (i > j) { t = i; i = j; j = t; }
    if (j > k) { t = j; j = k; k = t; }
    if (i > j) { t = i; i = j; j = t; }
    return A3[((size_t)i * (M + 1) + j) * (M + 1) + k];
}

/* smallest integer k with (k/N)^2 >= d2max / H^2, i.e. k^2 * H^2 >= d2max * N^2 (exact integers) */
static inline int up_index(long long d2max, int H) {
    long long X = d2max * (long long)N * N;
    long long HH = (long long)H * H;
    long long k = (long long)floor(sqrt((double)X / (double)HH));
    while (k > 0 && (k - 1) * (k - 1) * HH >= X) k--;
    while (k * k * HH < X) k++;
    if (k > M) k = M;
    return (int)k;
}

/* largest k with (k/N)^2 <= d2min / HF^2 (finest units), exact integers */
static long long HFg;
static inline int down_index_f(long long d2min) {
    long long X = d2min * (long long)N * N;
    long long HH = HFg * HFg;
    long long k = (long long)floor(sqrt((double)X / (double)HH));
    if (k < 0) k = 0;
    while (k * k * HH > X) k--;
    while ((k + 1) * (k + 1) * HH <= X) k++;
    if (k > M) k = M;
    return (int)k;
}

/* smallest k with (k/N)^2 >= d2max / HF^2 (finest units), exact integers */
static inline int up_index_f(long long d2max) {
    long long X = d2max * (long long)N * N;
    long long HH = HFg * HFg;
    long long k = (long long)floor(sqrt((double)X / (double)HH));
    if (k < 0) k = 0;
    while (k > 0 && (k - 1) * (k - 1) * HH >= X) k--;
    while (k * k * HH < X) k++;
    if (k > M) k = M;
    return (int)k;
}

static inline void range_abs(long long lo, long long hi, long long *mn, long long *mx) {
    long long a = llabs(lo), b = llabs(hi);
    *mx = a > b ? a : b;
    *mn = (lo <= 0 && hi >= 0) ? 0 : (a < b ? a : b);
}


static long long HF, HF2;
static double sum = 0, comp = 0;
static long long ncells = 0, npos = 0, nsplit = 0;
static int ADEPTH = 0;                 /* split every cube with g_max > 0 down to this depth (sharper area bounds) */
static double TAU = -1;                /* adaptive mode (C4_TAU): split while the estimated gain exceeds TAU */

/* one cube: lo[0] = r2 low, lo[1] = xi3, lo[2] = eta3, lo[3] = xi4, lo[4] = eta4 (finest units), side s */
static void cube(const long long *lo, long long s, int depth, int maxdepth) {
    ncells++;
    long long mn, mx, mny, mxy, d2min[6], d2max[6];
    long long x1l = lo[0], x1h = lo[0] + s;                 /* y2 = (x1, 0) */
    long long x2l = lo[1], x2h = lo[1] + s, y2l = lo[2], y2h = lo[2] + s;   /* y3 */
    long long x3l = lo[3], x3h = lo[3] + s, y3l = lo[4], y3h = lo[4] + s;   /* y4 */
    range_abs(-x1h, -x1l, &mn, &mx); d2min[0] = mn * mn; d2max[0] = mx * mx;
    range_abs(-x2h, -x2l, &mn, &mx); range_abs(-y2h, -y2l, &mny, &mxy); d2min[1] = mn * mn + mny * mny; d2max[1] = mx * mx + mxy * mxy;
    range_abs(-x3h, -x3l, &mn, &mx); range_abs(-y3h, -y3l, &mny, &mxy); d2min[2] = mn * mn + mny * mny; d2max[2] = mx * mx + mxy * mxy;
    if (d2min[1] > 4 * HF2 || d2min[2] > 4 * HF2) return;  /* some |y_j| > 2 throughout: integrand zero */
    range_abs(x1l - x2h, x1h - x2l, &mn, &mx); range_abs(-y2h, -y2l, &mny, &mxy); d2min[3] = mn * mn + mny * mny; d2max[3] = mx * mx + mxy * mxy;
    range_abs(x1l - x3h, x1h - x3l, &mn, &mx); range_abs(-y3h, -y3l, &mny, &mxy); d2min[4] = mn * mn + mny * mny; d2max[4] = mx * mx + mxy * mxy;
    range_abs(x2l - x3h, x2h - x3l, &mn, &mx); range_abs(y2l - y3h, y2h - y3l, &mny, &mxy); d2min[5] = mn * mn + mny * mny; d2max[5] = mx * mx + mxy * mxy;
    int P = 0, U = 0;
    for (int k = 0; k < 6; k++) {
        if (d2max[k] <= HF2) P |= 1 << k;
        else if (d2min[k] <= HF2) U |= 1 << k;
    }
    double gmx = gmax_tab[P][U];
    if (gmx <= 0) return;
    double gm = DIAG ? gmx : gmin_tab[P][U];
    if (TAU < 0 && ((gm < gmx && depth < maxdepth) || depth < ADEPTH)) {
        nsplit++;
        long long h = s / 2, sub[5];
        for (int m = 0; m < 32; m++) {
            for (int t = 0; t < 5; t++) sub[t] = lo[t] + (((m >> t) & 1) ? h : 0);
            cube(sub, h, depth + 1, maxdepth);
        }
        return;
    }
    if (TAU < 0 && gm <= 0) return;
    int ki[6];
    for (int k = 0; k < 6; k++) ki[k] = up_index_f(d2max[k]);
    double Amax = INFINITY;
    static const int TR[4][3] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};
    static const int LO4[4] = {3, 2, 1, 0};
    for (int t = 0; t < 4; t++) {
        int i0 = TR[t][0], i1 = TR[t][1], i2 = TR[t][2], l = LO4[t];
        double a3 = A3get(ki[edge_index(i0, i1)], ki[edge_index(i0, i2)], ki[edge_index(i1, i2)]);
        for (int q = 0; q < 3; q++) {
            int m = TR[t][q];
            double val = nextafter(a3 + C1[ki[edge_index(l, m)]], INFINITY);
            if (val < Amax) Amax = val;
        }
        if (L2) {   /* the fourth disc outside the union of two discs m1, m2 of the triple: A3(l,m1,m2) - A2(m1,m2) */
            for (int q = 0; q < 3; q++) {
                int m1 = TR[t][q], m2 = TR[t][(q + 1) % 3];
                double a3l = A3get(ki[edge_index(l, m1)], ki[edge_index(l, m2)], ki[edge_index(m1, m2)]);
                int kd = down_index_f(d2min[edge_index(m1, m2)]);
                double val = nextafter(nextafter(a3 + a3l, INFINITY) - L2[kd], INFINITY);
                if (val < Amax) Amax = val;
            }
        }
    }
    double w = ((double)(x1h * x1h - x1l * x1l) / 2.0) * (double)s * (double)s * (double)s * (double)s;
    w = w / ((double)HF2 * (double)HF2 * (double)HF2);
    if (TAU >= 0 && depth < maxdepth) {
        /* heuristic only (the decision does not affect validity): a lower ESTIMATE of the area from the triples at
         * the smallest distances, and the estimated gain of splitting */
        double Alo = M_PI;
        for (int t = 0; t < 4; t++) {
            int kd0 = down_index_f(d2min[edge_index(TR[t][0], TR[t][1])]);
            int kd1 = down_index_f(d2min[edge_index(TR[t][0], TR[t][2])]);
            int kd2 = down_index_f(d2min[edge_index(TR[t][1], TR[t][2])]);
            double v = A3get(kd0, kd1, kd2);
            if (v > Alo) Alo = v;
        }
        double est = w * (gmx * pow(Alo, -3.5) - (gm > 0 ? gm : 0) * pow(Amax, -3.5));
        if (est > TAU) {
            nsplit++;
            long long h = s / 2, sub[5];
            for (int m = 0; m < 32; m++) {
                for (int t = 0; t < 5; t++) sub[t] = lo[t] + (((m >> t) & 1) ? h : 0);
                cube(sub, h, depth + 1, maxdepth);
            }
            return;
        }
    }
    if (gm <= 0) return;
    double p = pow(Amax, -3.5) * (1.0 - 1e-13);
    /* w = int over the cube of r2 = ((x1h^2 - x1l^2)/2) * s^4 / HF^6 (exact in double: dyadic, < 2^53) */
    double term = gm * p * w * (1.0 - 1e-15);
    double t2 = sum + term;
    if (fabs(sum) >= fabs(term)) comp += (sum - t2) + term; else comp += (term - t2) + sum;
    sum = t2;
    npos++;
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: c4_a TABLE H [claim]\n"); return 1; }
    FILE *f = fopen(argv[1], "r");
    if (!f) { perror("table"); return 1; }
    int H = atoi(argv[2]);
    if (getenv("C4_DIAG")) { DIAG = 1; fprintf(stderr, "DIAGNOSTIC MODE: g_max used, the output is NOT a lower bound\n"); }
    char tag[8];
    if (fscanf(f, "%7s %d", tag, &N) != 2 || strcmp(tag, "N")) { fprintf(stderr, "bad table header\n"); return 1; }
    M = 2 * N;
    A3 = calloc((size_t)(M + 1) * (M + 1) * (M + 1), sizeof(double));
    C1 = calloc(M + 1, sizeof(double));
    char vs[64];
    long nA = 0, nC = 0;
    while (fscanf(f, "%7s", tag) == 1) {
        if (!strcmp(tag, "C")) {
            int k; if (fscanf(f, "%d %63s", &k, vs) != 2) return 2;
            C1[k] = nextafter(strtod(vs, NULL), INFINITY); nC++;
        } else if (!strcmp(tag, "A")) {
            int i, j, k; if (fscanf(f, "%d %d %d %63s", &i, &j, &k, vs) != 4) return 2;
            A3[((size_t)i * (M + 1) + j) * (M + 1) + k] = nextafter(strtod(vs, NULL), INFINITY); nA++;
        } else { fprintf(stderr, "bad tag %s\n", tag); return 2; }
    }
    fclose(f);
    long expectA = 0;
    for (int i = 1; i <= M; i++) for (int j = i; j <= M; j++) expectA += M - j + 1;
    if (nA != expectA || nC != M + 1) { fprintf(stderr, "table incomplete: %ld/%ld triples, %ld/%d pairs\n", nA, expectA, nC, M + 1); return 3; }
    /* self-checks of the tables: equilateral side 1 (if on the grid) and monotonicity spot checks */
    if (N % 1 == 0) {
        double eq = A3get(N, N, N);
        double L1 = 2 * M_PI / 3 - sqrt(3.0) / 2, T = (M_PI - sqrt(3.0)) / 2;
        double ref = 3 * M_PI - 3 * L1 + T;
        if (!(eq >= ref - 1e-12 && eq <= ref + 1e-9)) { fprintf(stderr, "table self-check failed: %.15f vs %.15f\n", eq, ref); return 4; }
    }
    for (int k = 1; k <= M; k++) if (C1[k] < C1[k - 1]) { fprintf(stderr, "C1 not monotone at %d\n", k); return 4; }
    if (getenv("C4_LTAB")) {
        FILE *g = fopen(getenv("C4_LTAB"), "r"); int n2; char tg[8], vv[64];
        if (!g || fscanf(g, "%7s %d", tg, &n2) != 2 || n2 != N) { fprintf(stderr, "bad L table\n"); return 6; }
        L2 = calloc(M + 1, sizeof(double)); int cnt = 0;
        while (fscanf(g, "%7s", tg) == 1) { int k; if (fscanf(g, "%d %63s", &k, vv) != 2) return 6; L2[k] = nextafter(strtod(vv, NULL), -INFINITY); cnt++; }
        fclose(g);
        if (cnt != M + 1) { fprintf(stderr, "L table incomplete\n"); return 6; }
        for (int k = 1; k <= M; k++) if (L2[k] < L2[k - 1]) { fprintf(stderr, "L table not monotone\n"); return 6; }
        fprintf(stderr, "using the two-disc refinement of the area bound (L table)\n");
    }
    for (int E = 0; E < 64; E++) gtab[E] = g_of(E);
    for (int P = 0; P < 64; P++) for (int U = 0; U < 64; U++) {
        if (P & U) { gmin_tab[P][U] = -1; continue; }
        double m = 1e9, mx = -1;
        for (int S = U;; S = (S - 1) & U) { double v = gtab[P | S]; if (v < m) m = v; if (v > mx) mx = v; if (S == 0) break; }
        gmin_tab[P][U] = m; gmax_tab[P][U] = mx;
    }
    /* report g of named graphs */
    int Epath = (1 << edge_index(0, 1)) | (1 << edge_index(1, 2)) | (1 << edge_index(2, 3));
    int Estar = (1 << edge_index(0, 1)) | (1 << edge_index(0, 2)) | (1 << edge_index(0, 3));
    int Ecyc = Epath | (1 << edge_index(0, 3));
    printf("g: path %.1f  star %.1f  4-cycle %.1f  K4 %.1f  paw %.1f\n", gtab[Epath], gtab[Estar], gtab[Ecyc], gtab[63],
           gtab[Estar | (1 << edge_index(1, 2))]);
    if (gtab[Epath] != 0 || gtab[Estar] != 0.5 || gtab[Ecyc] != 1 || gtab[63] != 1) { fprintf(stderr, "g self-check failed\n"); return 5; }

    /* Refinement: base cubes of side 1/H; a cube whose compatible graphs do not all have the same g (some
     * undetermined edge matters) is split into 32 sub-cubes, down to DEPTH levels. All coordinates are integers in
     * the finest unit 1/HF, HF = H * 2^DEPTH. */
    int DEPTH = getenv("C4_DEPTH") ? atoi(getenv("C4_DEPTH")) : 0;
    ADEPTH = getenv("C4_ADEPTH") ? atoi(getenv("C4_ADEPTH")) : 0;
    if (getenv("C4_TAU")) { TAU = atof(getenv("C4_TAU")); fprintf(stderr, "adaptive mode: tau %g, max depth from C4_DEPTH\n", TAU); }
    if (ADEPTH > DEPTH) DEPTH = ADEPTH;
    HF = (long long)H << DEPTH; HFg = HF;
    HF2 = HF * HF;
    const long long S0 = 1LL << DEPTH;
    const int R = 2 * H;
    for (int a = 0; a < R; a++) {
        for (int b = -R; b < R; b++) {
            for (int c = 0; c < R; c++) {
                long long lo[5] = {a * S0, b * S0, c * S0, 0, 0};
                /* quick skip: |y3| > 2 throughout */
                {
                    long long mn, mx, mny, mxy;
                    range_abs(-lo[1] - S0, -lo[1], &mn, &mx); range_abs(-lo[2] - S0, -lo[2], &mny, &mxy);
                    if (mn * mn + mny * mny > 4 * HF2) continue;
                }
                for (int d = -R; d < R; d++) {
                    for (int e = -R; e < R; e++) {
                        lo[3] = d * S0; lo[4] = e * S0;
                        cube(lo, S0, 0, DEPTH);
                    }
                }
            }
        }
        if ((a + 1) % 8 == 0) { fprintf(stderr, "r2 rows %d/%d: partial integral %.12e (%lld positive cubes, %lld splits)\n", a + 1, R, sum + comp, npos, nsplit); fflush(stderr); }
    }
    double S = (sum + comp) * (1.0 - 1e-12);
    /* constant 5 pi^{3/2} / 32, lowered: pi >= 3.14159265358979 */
    double pi_lo = 3.14159265358979;
    double K = 5.0 * pi_lo * sqrt(pi_lo) / 32.0 * (1.0 - 1e-15);
    double C4 = K * S * (1.0 - 1e-15);
    printf("%s", L2 ? "area bound with the two-disc refinement; " : "");
    printf("H %d (base cube side 1/%d), refinement depth %d (every cube to depth %d), table N %d: %lld cubes visited, %lld positive, %lld split\n", H, H, DEPTH, ADEPTH, N, ncells, npos, nsplit);
    printf("C_4 >= %.12f   [Dong-Cao: 0.005040573276]\n", C4);
    if (argc > 3) {
        double claim = atof(argv[3]);
        printf("claim C_4 >= %s: %s\n", argv[3], C4 >= claim ? "HOLDS" : "NOT IMPLIED");
    }
    return 0;
}
