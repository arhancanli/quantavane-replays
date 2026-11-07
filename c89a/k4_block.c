/* k4_block.c -- the monochromatic K4 count of an equal-weight step graphon, by a coarse block decomposition.
 *
 * Quantavane replays, our own code, written from the definition on the 89a page of teorth/optimizationproblems:
 * for a symmetric n x n matrix W with entries in [0,1] and classes of equal weight 1/n, the density of
 * monochromatic K4 is n^-4 * sum over all n^4 ordered (a,b,c,d), repeats included, of
 * prod_{6 pairs} W + prod_{6 pairs} (1 - W).  With W = M/den (integers M in [0, den]) this is
 * (N_red + N_blue) / (n^4 den^6), where for a colour X (red: X = M, blue: X = den - M)
 *
 *     N_X = sum_{a,b,c,d} X_ab X_ac X_ad X_bc X_bd X_cd = sum_a tr((X D_a)^3),   D_a = diag(X_a1, ..., X_an).
 *
 * The method.  The n classes are cut into m = n/s coarse classes of s consecutive classes.  Every coarse block is
 * checked: it is "constant" if all its s*s entries are equal, "fractional" otherwise.  We require (and check) that
 * every constant block holds 0 or den and that every diagonal coarse block is constant; otherwise we stop.  For a
 * colour X write X = P K P^T + F, where K = den*E, E is the 0/1 matrix of the coarse blocks constant at den in
 * this colour, P the class-to-coarse incidence, and F equals X on the fractional blocks and 0 elsewhere.  Fix a
 * root a, u = X_a (its row), D = diag(u), q = P^T u (q_I = sum of u over coarse class I).  With A = P K P^T D and
 * B = F D, tr((A+B)^3) = tr A^3 + 3 tr A^2 B + 3 tr A B^2 + tr B^3 (cyclicity), and
 *
 *   tr A^3   = den^3 sum_{I,J,L} E_IJ E_JL E_LI q_I q_J q_L = den^3 sum_{E_IJ=1} q_I q_J m_IJ,   m = E diag(q) E,
 *   tr A^2 B = den^2 sum_{(I,J) fractional} m_IJ G_JI,    G_JI = sum_{b in J, c in I} u_b F_bc u_c,
 *   tr A B^2 = den   sum_c u_c sum_{I,J} E_IJ x_c[I] x_c[J],   x_c[I] = sum_{b in I} u_b F_bc,
 *   tr B^3   = sum_{b,c,d} u_b u_c u_d F_bc F_cd F_db.
 *
 * In tr B^3 the three coarse classes of b, c, d form a triangle of the fractional graph (diagonal blocks are
 * constant, so the three are distinct), and the inner sum is the same for the six orders of a triangle (F is
 * symmetric), so tr B^3 = 6 * sum over unordered triangles {I,K,J} of S_IKJ(u),
 * S_IKJ(u) = sum_{b in I, c in K, d in J} u_b u_c u_d F_bc F_cd F_db.  On a coarse class V that is not a
 * fractional partner of the root's coarse class A0, u is the constant den*E_{A0,V}; so S needs only
 *   (0 fractional vertices) den^3 E E E tr(F_IK F_KJ F_JI),
 *   (1, at V)               den^2 E E sum_{b in V} u_b (F F F)_bb,
 *   (2, constant vertex W)  den E sum_{b in V1, c in V2} u_b u_c F_bc (F_{V2 W} F_{W V1})_cb,
 *   (3)                     the full sum.
 *
 * Arithmetic: every first-level sum of products of entries (q, m, x, G, the triangle traces and products) is
 * formed exactly in unsigned 64- or 128-bit integers; the bounds, for den <= 2^30, s <= 40, n <= 7680, m <= 192:
 * q <= 2^35.4, m_IJ <= 2^43, sum_J q_J m_IJ <= 2^86.1, x and (F F)_cb <= 2^65.4, G <= 2^100.7,
 * (F F F)_bb <= 2^101.7, tr(F F F) <= 2^106.3, sum_d u_d F_cd F_db <= 2^95.4.  All further products and sums are
 * taken modulo each of six moduli: 2^64 and the primes 2^61-1, 2^62-57, 2^60-93, 2^59-55, 2^63-25.  The caller
 * reconstructs N_X by the Chinese remainder theorem from the first four (product about 2^247, above
 * n^4 den^6 <= 2^229) and checks it against the last two.
 *
 * usage: k4_block MATRIX.bin n s den [root ...]
 *   MATRIX.bin: n*n little-endian uint32, row-major.  Prints, for each colour, the residues of N_X and of its four
 *   parts, and for each listed root a the residues of tr((X D_a)^3) (for the second program, k4_dense.c).
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef unsigned __int128 u128;
#define NMOD 6
static const int MK[NMOD] = {64, 61, 62, 60, 59, 63};
static const uint64_t MC[NMOD] = {0, 1, 57, 93, 55, 25};
static uint64_t MP[NMOD];

static inline uint64_t red(u128 x, int j) {
    if (j == 0) return (uint64_t)x;
    const int k = MK[j];
    const u128 mask = (((u128)1) << k) - 1;
    while (x >> k) x = (x >> k) * MC[j] + (x & mask);
    uint64_t r = (uint64_t)x;
    if (r >= MP[j]) r -= MP[j];
    return r;
}
static inline uint64_t mulm(uint64_t a, uint64_t b, int j) { return red((u128)a * b, j); }
static inline uint64_t addm(uint64_t a, uint64_t b, int j) { return red((u128)a + b, j); }

static void die(const char *msg) { fprintf(stderr, "k4_block: %s\n", msg); exit(2); }

static void print_u128(const char *tag, u128 v) {
    char buf[64]; int i = 63; buf[i] = 0;
    if (v == 0) buf[--i] = '0';
    while (v) { buf[--i] = (char)('0' + (int)(v % 10)); v /= 10; }
    printf("%s %s\n", tag, buf + i);
}

int main(int argc, char **argv) {
    if (argc < 5) die("usage: k4_block MATRIX.bin n s den [root ...]");
    for (int j = 1; j < NMOD; j++) MP[j] = (((uint64_t)1) << MK[j]) - MC[j];
    const long n = atol(argv[2]), s = atol(argv[3]);
    const uint64_t den = strtoull(argv[4], NULL, 10);
    if (n <= 0 || s <= 0 || n % s) die("s must divide n");
    if (den == 0 || den > (1ULL << 30)) die("den must be in [1, 2^30]");
    if (s > 40 || n > 7680) die("bounds of the comment assume s <= 40, n <= 7680");
    const long m = n / s;
    if (m > 192) die("coarse classes > 192");
    const int nroots = argc - 5;
    long *roots = calloc(nroots + 1, sizeof(long));
    for (int i = 0; i < nroots; i++) { roots[i] = atol(argv[5 + i]); if (roots[i] < 0 || roots[i] >= n) die("bad root"); }

    uint32_t *M = malloc((size_t)n * n * sizeof(uint32_t));
    FILE *fp = fopen(argv[1], "rb");
    if (!fp || fread(M, sizeof(uint32_t), (size_t)n * n, fp) != (size_t)n * n) die("cannot read matrix");
    { uint32_t extra; if (fread(&extra, 1, 1, fp) != 0) die("matrix file longer than n*n entries"); }
    fclose(fp);
    for (long i = 0; i < n; i++)
        for (long k = 0; k < n; k++) {
            if (M[i * n + k] > den) die("entry above den");
            if (M[i * n + k] != M[k * n + i]) die("matrix not symmetric");
        }

    /* coarse structure (the same for both colours) */
    uint8_t *frac = calloc(m * m, 1);
    uint32_t *cval = calloc(m * m, sizeof(uint32_t));
    for (long I = 0; I < m; I++)
        for (long J = 0; J < m; J++) {
            uint32_t v = M[(I * s) * n + J * s]; int fr = 0;
            for (long i = I * s; i < I * s + s && !fr; i++)
                for (long k = J * s; k < J * s + s; k++) if (M[i * n + k] != v) { fr = 1; break; }
            frac[I * m + J] = (uint8_t)fr; cval[I * m + J] = v;
            if (!fr && v != 0 && v != den) die("a constant coarse block holds a value other than 0 or den");
            if (I == J && fr) die("a diagonal coarse block is not constant");
        }
    long *deg = calloc(m, sizeof(long)), *fr = calloc(m * m, sizeof(long)), *frpos = malloc(m * m * sizeof(long));
    for (long I = 0; I < m; I++)
        for (long J = 0; J < m; J++) { frpos[I * m + J] = -1; if (frac[I * m + J]) { frpos[I * m + J] = deg[I]; fr[I * m + deg[I]++] = J; } }
    long ntri = 0;
    long (*tri)[3] = malloc(sizeof(long[3]) * (size_t)(m * m * m / 6 + 1));
    for (long I = 0; I < m; I++)
        for (long a = 0; a < deg[I]; a++) {
            long K = fr[I * m + a]; if (K <= I) continue;
            for (long b = 0; b < deg[K]; b++) {
                long J = fr[K * m + b]; if (J <= K || !frac[J * m + I]) continue;
                tri[ntri][0] = I; tri[ntri][1] = K; tri[ntri][2] = J; ntri++;
            }
        }
    long nfrac = 0; for (long I = 0; I < m; I++) nfrac += deg[I];
    printf("n %ld s %ld m %ld den %llu fractional_blocks %ld triangles %ld\n", n, s, m, (unsigned long long)den, nfrac, ntri);

    uint32_t *X = malloc((size_t)n * n * sizeof(uint32_t));
    const long W = (m + 63) / 64;
    uint64_t *Ebits = calloc(m * W, sizeof(uint64_t));
    uint8_t *E = calloc(m * m, 1);
    u128 *T0 = malloc(sizeof(u128) * (ntri + 1));
    uint64_t *q = malloc(sizeof(uint64_t) * m);
    uint64_t *mIJ = malloc(sizeof(uint64_t) * m * m);
    u128 *x = malloc(sizeof(u128) * n * 192);               /* x[c*maxdeg + idx], idx over fr(cc(c)) */
    long maxdeg = 0; for (long I = 0; I < m; I++) if (deg[I] > maxdeg) maxdeg = deg[I];
    if (maxdeg > 192) die("degree");
    u128 *wtmp = malloc(sizeof(u128) * s * s);
    /* per-A0 caches for triangles with 1 or 2 fractional vertices, reduced modulo each modulus */
    int *tf = malloc(sizeof(int) * (ntri + 1));
    uint64_t *Pv = malloc(sizeof(uint64_t) * (ntri + 1) * s * NMOD);
    uint64_t *Qw = malloc(sizeof(uint64_t) * (ntri + 1) * s * s * NMOD);
    if (!X || !x || !Pv || !Qw || !T0) die("out of memory");

    for (int colour = 0; colour < 2; colour++) {
        for (long i = 0; i < n * n; i++) X[i] = colour == 0 ? M[i] : (uint32_t)(den - M[i]);
        memset(Ebits, 0, m * W * sizeof(uint64_t));
        for (long I = 0; I < m; I++)
            for (long J = 0; J < m; J++) {
                uint32_t v = colour == 0 ? cval[I * m + J] : (uint32_t)(den - cval[I * m + J]);
                E[I * m + J] = (uint8_t)(!frac[I * m + J] && v == den);
                if (E[I * m + J]) Ebits[I * W + J / 64] |= 1ULL << (J % 64);
            }
#define XX(i, k) ((uint64_t)X[(size_t)(i) * n + (k)])
        /* triangle traces tr(F_IK F_KJ F_JI), exact */
        for (long t = 0; t < ntri; t++) {
            long I = tri[t][0], K = tri[t][1], J = tri[t][2]; u128 tr = 0;
            for (long b = I * s; b < I * s + s; b++)
                for (long d = J * s; d < J * s + s; d++) {
                    u128 w = 0;
                    for (long c = K * s; c < K * s + s; c++) w += (u128)(XX(b, c) * XX(c, d));
                    tr += w * XX(d, b);
                }
            T0[t] = tr;
        }
        uint64_t tot[4][NMOD]; memset(tot, 0, sizeof tot);    /* sum over roots of the four parts */
        uint64_t den1[NMOD], den2[NMOD], den3[NMOD];
        for (int j = 0; j < NMOD; j++) { den1[j] = red(den, j); den2[j] = mulm(den1[j], den1[j], j); den3[j] = mulm(den2[j], den1[j], j); }

        for (long A0 = 0; A0 < m; A0++) {
            /* classify triangles relative to A0 and build the caches */
            uint64_t C0[NMOD]; memset(C0, 0, sizeof C0);
            for (long t = 0; t < ntri; t++) {
                long V[3] = {tri[t][0], tri[t][1], tri[t][2]}; int f = 0, ok = 1;
                for (int i = 0; i < 3; i++) { if (frac[A0 * m + V[i]]) f++; else if (!E[A0 * m + V[i]]) ok = 0; }
                tf[t] = ok ? f : -1;              /* -1: a constant vertex with u = 0, so S = 0 */
                if (!ok) continue;
                if (f == 0) { for (int j = 0; j < NMOD; j++) C0[j] = addm(C0[j], red(T0[t], j), j); }
                else if (f == 1) {
                    /* anchor V (fractional), the other two constant: (F_VY F_YZ F_ZV)_bb for b in V */
                    int ia = frac[A0 * m + V[0]] ? 0 : frac[A0 * m + V[1]] ? 1 : 2;
                    long Va = V[ia], Y = V[(ia + 1) % 3], Z = V[(ia + 2) % 3];
                    for (long bi = 0; bi < s; bi++) {
                        long b = Va * s + bi; u128 acc = 0;
                        for (long d = Z * s; d < Z * s + s; d++) {
                            u128 w = 0;
                            for (long c = Y * s; c < Y * s + s; c++) w += (u128)(XX(b, c) * XX(c, d));
                            acc += w * XX(d, b);
                        }
                        for (int j = 0; j < NMOD; j++) Pv[((size_t)t * s + bi) * NMOD + j] = red(acc, j);
                    }
                } else if (f == 2) {
                    /* constant vertex Wc; fractional V1 < V2 (in triangle order): (F_{V2 Wc} F_{Wc V1})_{c b} */
                    int iw = !frac[A0 * m + V[0]] ? 0 : !frac[A0 * m + V[1]] ? 1 : 2;
                    long Wc = V[iw], V1 = V[(iw + 1) % 3], V2 = V[(iw + 2) % 3];
                    for (long bi = 0; bi < s; bi++)
                        for (long ci = 0; ci < s; ci++) {
                            long b = V1 * s + bi, c = V2 * s + ci; u128 w = 0;
                            for (long d = Wc * s; d < Wc * s + s; d++) w += (u128)(XX(c, d) * XX(d, b));
                            for (int j = 0; j < NMOD; j++) Qw[(((size_t)t * s + bi) * s + ci) * NMOD + j] = red(w, j);
                        }
                }
            }
            for (long r = 0; r < s; r++) {
                const long a = A0 * s + r;
                const uint32_t *u = X + (size_t)a * n;
                for (long L = 0; L < m; L++) { uint64_t acc = 0; for (long b = L * s; b < L * s + s; b++) acc += u[b]; q[L] = acc; }
                uint64_t part[4][NMOD]; memset(part, 0, sizeof part);
                /* m_IJ for E_IJ = 1 and for fractional (I,J); tr A^3 */
                for (long I = 0; I < m; I++) {
                    u128 rowacc = 0;
                    for (long J = 0; J < m; J++) {
                        if (!E[I * m + J] && !frac[I * m + J]) continue;
                        uint64_t acc = 0;
                        for (long w = 0; w < W; w++) {
                            uint64_t bits = Ebits[I * W + w] & Ebits[J * W + w];
                            while (bits) { long L = w * 64 + __builtin_ctzll(bits); acc += q[L]; bits &= bits - 1; }
                        }
                        mIJ[I * m + J] = acc;
                        if (E[I * m + J]) rowacc += (u128)q[J] * acc;
                    }
                    if (rowacc) {
                        u128 v = rowacc * q[I];    /* <= 2^121.5 */
                        for (int j = 0; j < NMOD; j++) part[0][j] = addm(part[0][j], red(v, j), j);
                    }
                }
                /* x_c[I] = sum_{b in I} u_b F_bc for I in fr(cc(c)) */
                for (long c = 0; c < n; c++) {
                    long C = c / s;
                    for (long idx = 0; idx < deg[C]; idx++) {
                        long I = fr[C * m + idx]; u128 acc = 0;
                        for (long b = I * s; b < I * s + s; b++) acc += (u128)((uint64_t)u[b] * XX(b, c));
                        x[c * maxdeg + idx] = acc;
                    }
                }
                /* tr A^2 B: G_JI = sum_{b in J} u_b x_b[I], I in fr(J) */
                for (long J = 0; J < m; J++)
                    for (long idx = 0; idx < deg[J]; idx++) {
                        long I = fr[J * m + idx]; u128 G = 0;
                        for (long b = J * s; b < J * s + s; b++) G += x[b * maxdeg + idx] * (u128)u[b];
                        for (int j = 0; j < NMOD; j++) part[1][j] = addm(part[1][j], mulm(red(G, j), red(mIJ[I * m + J], j), j), j);
                    }
                /* tr A B^2 */
                for (long c = 0; c < n; c++) {
                    if (!u[c]) continue;
                    long C = c / s, dC = deg[C];
                    for (int j = 0; j < NMOD; j++) {
                        uint64_t xr[192]; uint64_t t = 0;
                        for (long i = 0; i < dC; i++) xr[i] = red(x[c * maxdeg + i], j);
                        for (long i = 0; i < dC; i++) {
                            long I = fr[C * m + i]; u128 w = 0;
                            for (long k = 0; k < dC; k++) if (E[I * m + fr[C * m + k]]) w += xr[k];
                            t = addm(t, mulm(xr[i], red(w, j), j), j);
                        }
                        part[2][j] = addm(part[2][j], mulm(t, red(u[c], j), j), j);
                    }
                }
                /* tr B^3 = 6 sum_t S_t(u) */
                uint64_t sb[NMOD]; for (int j = 0; j < NMOD; j++) sb[j] = mulm(C0[j], den3[j], j);
                for (long t = 0; t < ntri; t++) {
                    if (tf[t] <= 0) continue;
                    long V[3] = {tri[t][0], tri[t][1], tri[t][2]};
                    if (tf[t] == 1) {
                        int ia = frac[A0 * m + V[0]] ? 0 : frac[A0 * m + V[1]] ? 1 : 2; long Va = V[ia];
                        for (int j = 0; j < NMOD; j++) {
                            uint64_t acc = 0;
                            for (long bi = 0; bi < s; bi++) acc = addm(acc, mulm(red(u[Va * s + bi], j), Pv[((size_t)t * s + bi) * NMOD + j], j), j);
                            sb[j] = addm(sb[j], mulm(acc, den2[j], j), j);
                        }
                    } else if (tf[t] == 2) {
                        int iw = !frac[A0 * m + V[0]] ? 0 : !frac[A0 * m + V[1]] ? 1 : 2;
                        long V1 = V[(iw + 1) % 3], V2 = V[(iw + 2) % 3];
                        for (long bi = 0; bi < s; bi++)
                            for (long ci = 0; ci < s; ci++) {
                                long b = V1 * s + bi, c = V2 * s + ci;
                                u128 coef = (u128)((uint64_t)u[b] * u[c]) * XX(b, c);   /* <= 2^90 */
                                wtmp[bi * s + ci] = coef;
                            }
                        for (int j = 0; j < NMOD; j++) {
                            uint64_t acc = 0;
                            for (long bi = 0; bi < s; bi++)
                                for (long ci = 0; ci < s; ci++)
                                    acc = addm(acc, mulm(red(wtmp[bi * s + ci], j), Qw[(((size_t)t * s + bi) * s + ci) * NMOD + j], j), j);
                            sb[j] = addm(sb[j], mulm(acc, den1[j], j), j);
                        }
                    } else {  /* three fractional vertices */
                        long I = V[0], K = V[1], J = V[2];
                        for (int j = 0; j < NMOD; j++) {
                            uint64_t acc = 0;
                            for (long b = I * s; b < I * s + s; b++)
                                for (long c = K * s; c < K * s + s; c++) {
                                    u128 z = 0;
                                    for (long d = J * s; d < J * s + s; d++) z += (u128)((uint64_t)u[d] * XX(c, d)) * XX(d, b);
                                    u128 coef = (u128)((uint64_t)u[b] * u[c]) * XX(b, c);
                                    acc = addm(acc, mulm(red(coef, j), red(z, j), j), j);
                                }
                            sb[j] = addm(sb[j], acc, j);
                        }
                    }
                }
                for (int j = 0; j < NMOD; j++) part[3][j] = mulm(sb[j], 6, j);
                uint64_t Na[NMOD];
                for (int j = 0; j < NMOD; j++) {
                    uint64_t v = mulm(part[0][j], den3[j], j);
                    v = addm(v, mulm(mulm(part[1][j], den2[j], j), 3, j), j);
                    v = addm(v, mulm(mulm(part[2][j], den1[j], j), 3, j), j);
                    v = addm(v, part[3][j], j);
                    Na[j] = v;
                    tot[0][j] = addm(tot[0][j], mulm(part[0][j], den3[j], j), j);
                    tot[1][j] = addm(tot[1][j], mulm(mulm(part[1][j], den2[j], j), 3, j), j);
                    tot[2][j] = addm(tot[2][j], mulm(mulm(part[2][j], den1[j], j), 3, j), j);
                    tot[3][j] = addm(tot[3][j], part[3][j], j);
                }
                for (int i = 0; i < nroots; i++)
                    if (roots[i] == a) {
                        printf("ROOT colour %d a %ld", colour, a);
                        for (int j = 0; j < NMOD; j++) printf(" %llu", (unsigned long long)Na[j]);
                        printf("\n");
                    }
            }
            if (A0 % 32 == 31) { fprintf(stderr, "colour %d: coarse class %ld of %ld done\n", colour, A0 + 1, m); }
        }
        const char *nm[4] = {"AAA", "AAB", "ABB", "BBB"};
        for (int p = 0; p < 4; p++) {
            printf("PART colour %d %s", colour, nm[p]);
            for (int j = 0; j < NMOD; j++) printf(" %llu", (unsigned long long)tot[p][j]);
            printf("\n");
        }
        printf("TOTAL colour %d", colour);
        for (int j = 0; j < NMOD; j++) {
            uint64_t v = 0; for (int p = 0; p < 4; p++) v = addm(v, tot[p][j], j);
            printf(" %llu", (unsigned long long)v);
        }
        printf("\n");
        fflush(stdout);
    }
    (void)print_u128;
    return 0;
}
