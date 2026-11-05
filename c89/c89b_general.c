/* c89b_general.c: the monochromatic K_5 count numerators of a step graphon on m equal blocks with an arbitrary
 * symmetric matrix of red densities R/D (PR #277's 192-block certificate), modulo primes below 2^28.
 *
 * With x_1, ..., x_5 the blocks of five vertices (repeats allowed), N(X) = sum over all m^5 tuples of
 * prod_{i<j} X[x_i][x_j], and the density is (N(R) + N(D - R)) / (m^5 D^10). Two organisations of the same sum,
 * computed separately:
 *   pairs:   N = sum_{a,b} X_ab tr(S^3), S = diag(X_a o X_b) X          (two vertices, then a weighted triangle)
 *   triples: N = sum_{c,d,e} X_cd X_de X_ec v^T X v, v_x = X_xc X_xd X_xe (a triangle, then a weighted edge)
 * Each uses its own symmetry: a <= b with multiplicity 2 off the diagonal, resp. c <= d <= e with multiplicity
 * 6, 3 or 1. Products of two residues are below 2^56, so up to 256 of them are summed in 64 bits before reducing.
 *
 * usage: c89b_general pairs|triples PART NPARTS p1 [p2 ...] < matrix   (matrix: m, D, then the m x m entries R)
 * prints "p R residue" and "p B residue" (B = D - R) for each prime, for this part of the outer loop.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXM 256

static int m;
static long D;
static long R[MAXM][MAXM];
static uint32_t X[MAXM][MAXM];

static uint64_t pairs(uint64_t p, int part, int nparts) {
    static uint32_t S[MAXM][MAXM];
    static uint64_t acc[MAXM];
    uint32_t u[MAXM];
    uint64_t N = 0;
    long idx = 0;
    for (int a = 0; a < m; a++)
        for (int b = a; b < m; b++, idx++) {
            if (idx % nparts != part) continue;
            for (int d = 0; d < m; d++) u[d] = (uint32_t)((uint64_t)X[a][d] * X[b][d] % p);
            for (int c = 0; c < m; c++)
                for (int d = 0; d < m; d++) S[c][d] = (uint32_t)((uint64_t)u[c] * X[c][d] % p);
            uint64_t tr = 0;
            for (int c = 0; c < m; c++) {
                for (int e = 0; e < m; e++) acc[e] = 0;
                for (int d = 0; d < m; d++) {
                    uint64_t s = S[c][d];
                    const uint32_t *row = S[d];
                    for (int e = 0; e < m; e++) acc[e] += s * row[e];
                }
                uint64_t t = 0;
                for (int e = 0; e < m; e++) t += (acc[e] % p) * S[e][c] % p;
                tr = (tr + t % p) % p;
            }
            uint64_t f = (uint64_t)X[a][b] * tr % p;
            N = (N + (a == b ? f : 2 * f)) % p;
        }
    return N;
}

static uint64_t triples(uint64_t p, int part, int nparts) {
    uint32_t v[MAXM];
    uint64_t N = 0;
    long idx = 0;
    for (int c = 0; c < m; c++)
        for (int d = c; d < m; d++)
            for (int e = d; e < m; e++, idx++) {
                if (idx % nparts != part) continue;
                uint64_t T = (uint64_t)X[c][d] * X[d][e] % p * X[e][c] % p;
                for (int x = 0; x < m; x++) v[x] = (uint32_t)((uint64_t)X[x][c] * X[x][d] % p * X[x][e] % p);
                uint64_t q = 0;
                for (int a = 0; a < m; a++) {
                    uint64_t w = 0;
                    for (int b = 0; b < m; b++) w += (uint64_t)X[a][b] * v[b];
                    q = (q + (w % p) * v[a]) % p;
                }
                int mult = (c == d && d == e) ? 1 : (c == d || d == e) ? 3 : 6;
                N = (N + T * q % p * mult) % p;
            }
    return N;
}

int main(int argc, char **argv) {
    if (argc < 5 || (strcmp(argv[1], "pairs") && strcmp(argv[1], "triples"))) {
        fprintf(stderr, "usage: c89b_general pairs|triples PART NPARTS p1 [p2 ...] < matrix\n");
        return 2;
    }
    int part = atoi(argv[2]), nparts = atoi(argv[3]);
    if (nparts < 1 || part < 0 || part >= nparts) { fprintf(stderr, "bad part\n"); return 2; }
    if (scanf("%d %ld", &m, &D) != 2 || m < 1 || m > MAXM || D < 1 || D >= (1L << 28)) {
        fprintf(stderr, "bad m or D\n");
        return 2;
    }
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
            if (scanf("%ld", &R[i][j]) != 1 || R[i][j] < 0 || R[i][j] > D) { fprintf(stderr, "bad entry\n"); return 2; }
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
            if (R[i][j] != R[j][i]) { fprintf(stderr, "matrix not symmetric\n"); return 2; }
    for (int k = 4; k < argc; k++) {
        uint64_t p = strtoull(argv[k], NULL, 10);
        if (p <= (uint64_t)D || p >= (1ULL << 28)) { fprintf(stderr, "prime must lie in (D, 2^28)\n"); return 2; }
        for (int colour = 0; colour < 2; colour++) {
            for (int i = 0; i < m; i++)
                for (int j = 0; j < m; j++) X[i][j] = (uint32_t)((colour ? D - R[i][j] : R[i][j]) % (long)p);
            uint64_t N = strcmp(argv[1], "pairs") ? triples(p, part, nparts) : pairs(p, part, nparts);
            printf("%llu %c %llu\n", (unsigned long long)p, colour ? 'B' : 'R', (unsigned long long)N);
            fflush(stdout);
        }
    }
    return 0;
}
