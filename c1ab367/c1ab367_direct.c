/* Program B (core) for the step-function certificates of PR #367 (1a) and PR #368 (1b): plain double loops in exact
 * integer arithmetic, no product of big integers.
 *
 * Input on stdin: a first line "1a n" or "1b n Q", then n non-negative decimal integers, one per line.
 *  1a: c_k = sum_{i+j=k} a_i a_j for k = 0..2n-2, each as 2 * sum_{i<j} a_i a_j + [k even] a_{k/2}^2, accumulated in
 *      unsigned 128-bit integers; larger heights are split into the fewest equal limbs (two to four) for which the
 *      size bounds hold, the partial sums of each degree are kept apart and recombined in 384 bits. Every bound on the sizes is checked before the loop, so no
 *      accumulator can overflow.
 *  1b: v_m = sum_{0 <= i, i+m < n} k_i (Q - k_{i+m}) for m = -(n-1)..n-1, straight from the definition.
 * Output: n, S = sum of the heights, the largest value (hex, 384 bits), every index attaining it, and the fingerprint
 * sum_t v_t 3^t mod 2^61 - 1 of the whole sequence (t = k, or t = m + n - 1). Exit status 0, or 2 on bad input.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef unsigned __int128 u128;
#define NW 6
typedef struct { uint64_t w[NW]; } u256;         /* 384 bits, little-endian words (name kept short) */

static const uint64_t PR = (((uint64_t)1) << 61) - 1;

static uint64_t mod_u128(u128 x) { return (uint64_t)(x % PR); }

static uint64_t mod_u256(u256 x) {               /* Horner in base 2^64 */
    u128 r = 0;
    for (int i = NW - 1; i >= 0; i--) r = ((r << 64) + x.w[i]) % PR;
    return (uint64_t)r;
}

static u256 from128(u128 x) { u256 r; memset(&r, 0, sizeof r); r.w[0] = (uint64_t)x; r.w[1] = (uint64_t)(x >> 64); return r; }

static u256 add_shift(u256 a, u128 x, int sh) {  /* a + x * 2^sh, 0 <= sh <= 256; no carry may leave the top */
    uint64_t t[NW + 3];
    memset(t, 0, sizeof t);
    int q = sh / 64, r = sh % 64;
    uint64_t lo = (uint64_t)x, hi = (uint64_t)(x >> 64);
    if (r == 0) { t[q] = lo; t[q + 1] = hi; }
    else {
        t[q] = lo << r;
        t[q + 1] = (lo >> (64 - r)) | (hi << r);
        t[q + 2] = hi >> (64 - r);
    }
    for (int i = NW; i < NW + 3; i++) if (t[i]) { fprintf(stderr, "overflow in add_shift\n"); exit(2); }
    u128 carry = 0;
    for (int i = 0; i < NW; i++) {
        u128 s2 = (u128)a.w[i] + t[i] + carry;
        a.w[i] = (uint64_t)s2;
        carry = s2 >> 64;
    }
    if (carry) { fprintf(stderr, "overflow in add_shift\n"); exit(2); }
    return a;
}

static int cmp256(u256 a, u256 b) {
    for (int i = NW - 1; i >= 0; i--) {
        if (a.w[i] != b.w[i]) return a.w[i] < b.w[i] ? -1 : 1;
    }
    return 0;
}

static int bits128(u128 x) { int b = 0; while (x) { b++; x >>= 1; } return b; }

static int read_u128(u128 *out) {
    char buf[128];
    if (scanf("%127s", buf) != 1) return 0;
    u128 v = 0;
    for (char *p = buf; *p; p++) {
        if (*p < '0' || *p > '9') return 0;
        u128 nv = v * 10 + (u128)(*p - '0');
        if (nv / 10 != v) return 0;              /* overflow */
        v = nv;
    }
    *out = v;
    return 1;
}

static void print_u128(u128 x) {
    char buf[64]; int i = 63; buf[i] = 0;
    if (x == 0) buf[--i] = '0';
    while (x) { buf[--i] = (char)('0' + (int)(x % 10)); x /= 10; }
    printf("%s", buf + i);
}

static void print256(u256 x) {
    printf("0x");
    for (int i = NW - 1; i >= 0; i--) printf("%016llx", (unsigned long long)x.w[i]);
}

#define MAXARG 64

int main(void) {
    char mode[8];
    long n;
    if (scanf("%7s %ld", mode, &n) != 2 || n < 1 || n > (1L << 22)) { fprintf(stderr, "bad header\n"); return 2; }
    int is1a = strcmp(mode, "1a") == 0, is1b = strcmp(mode, "1b") == 0;
    if (!is1a && !is1b) { fprintf(stderr, "bad mode\n"); return 2; }
    u128 Q = 0;
    if (is1b && !read_u128(&Q)) { fprintf(stderr, "bad Q\n"); return 2; }
    u128 *a = malloc(sizeof(u128) * (size_t)n);
    u128 mx = 0, S = 0;
    for (long i = 0; i < n; i++) {
        if (!read_u128(&a[i])) { fprintf(stderr, "bad height %ld\n", i); return 2; }
        if (a[i] > mx) mx = a[i];
        S += a[i];
        if (S < a[i]) { fprintf(stderr, "sum overflow\n"); return 2; }
    }
    int nb = bits128((u128)n), mb = bits128(mx);
    u256 best;
    memset(&best, 0, sizeof best);
    long arg[MAXARG]; int narg = 0, first = 1;
    u128 fp = 0, pw = 1;                         /* fingerprint, 3^t mod p */
    long T = 2 * n - 1;
    if (is1a) {
        int split = (2 * mb + nb + 1 > 127);
        /* NL limbs of LB bits: the fewest limbs (2 to 4) with NL 2^(2 LB) n 2 < 2^127, so that no accumulator below
           can overflow (each limb product < 2^(2 LB), at most NL per degree and pair, at most n/2 pairs, doubled,
           plus the diagonal) */
        int NL = 0, LB = 0;
        if (split) {
            for (int t = 2; t <= 4 && !NL; t++) {
                int lb = (mb + t - 1) / t;
                if (2 * lb + 3 + nb + 1 <= 127 && lb <= 63) { NL = t; LB = lb; }
            }
            if (!NL) { fprintf(stderr, "heights too large for the limb bounds\n"); return 2; }
        }
        uint64_t *lim = malloc(8 * (size_t)n * 4);
        for (long i = 0; i < n; i++)
            for (int p = 0; p < 4; p++)
                lim[4 * i + p] = (split && p < NL) ? (uint64_t)((a[i] >> (LB * p)) & ((((u128)1) << LB) - 1)) : 0;
        printf("B: 1a, n = %ld, largest height %d bits, %s", n, mb, split ? "limbs (count, bits): " : "one 128-bit accumulator");
        if (split) printf("%d, %d", NL, LB);
        printf("\n");
        for (long k = 0; k < T; k++) {
            long i0 = k - (n - 1) > 0 ? k - (n - 1) : 0;
            u256 c;
            if (!split) {
                u128 acc = 0;
                for (long i = i0; 2 * i < k; i++) acc += a[i] * a[k - i];
                acc *= 2;
                if (k % 2 == 0) acc += a[k / 2] * a[k / 2];
                c = from128(acc);
            } else {
                u128 acc[7] = {0, 0, 0, 0, 0, 0, 0};
                for (long i = i0; 2 * i < k; i++) {
                    const uint64_t *x = lim + 4 * i, *y = lim + 4 * (k - i);
                    for (int p = 0; p < NL; p++)
                        for (int q = 0; q < NL; q++) acc[p + q] += (u128)x[p] * y[q];
                }
                for (int d = 0; d < 2 * NL - 1; d++) acc[d] *= 2;
                if (k % 2 == 0) {
                    const uint64_t *x = lim + 4 * (k / 2);
                    for (int p = 0; p < NL; p++)
                        for (int q = 0; q < NL; q++) acc[p + q] += (u128)x[p] * x[q];
                }
                c = from128(acc[0]);
                for (int d = 1; d < 2 * NL - 1; d++) c = add_shift(c, acc[d], LB * d);
            }
            int s = first ? 1 : cmp256(c, best);
            if (s > 0) { best = c; narg = 0; first = 0; }
            if (s >= 0 && narg < MAXARG) arg[narg++] = k;
            fp = (fp + (u128)mod_u256(c) * pw) % PR;
            pw = (pw * 3) % PR;
        }
    } else {
        if (mx > Q) { printf("B: a height exceeds Q\nB: VERDICT REFUSED\n"); return 0; }
        if (2 * bits128(Q) + nb + 1 > 127) { fprintf(stderr, "Q too large for 128 bits\n"); return 2; }
        printf("B: 1b, n = %ld, Q = ", n); print_u128(Q); printf("\n");
        for (long t = 0; t < T; t++) {
            long m = t - (n - 1);
            long i0 = m < 0 ? -m : 0, i1 = m > 0 ? n - m : n;
            u128 acc = 0;
            for (long i = i0; i < i1; i++) acc += a[i] * (Q - a[i + m]);
            u256 c = from128(acc);
            int s = first ? 1 : cmp256(c, best);
            if (s > 0) { best = c; narg = 0; first = 0; }
            if (s >= 0 && narg < MAXARG) arg[narg++] = m;
            fp = (fp + (u128)mod_u128(acc) * pw) % PR;
            pw = (pw * 3) % PR;
        }
    }
    printf("B: n %ld\nB: S ", n); print_u128(S); printf("\nB: max "); print256(best); printf("\nB: argmax");
    for (int i = 0; i < narg; i++) printf(" %ld", arg[i]);
    printf("\nB: fingerprint %llu\n", (unsigned long long)fp);
    return 0;
}
