/* Program B (core) for the 1d step function of PR #372: its own reader, then plain loops over every knot.
 *
 *   c1d372_b FILE [BUMP_INDEX BUMP_DELTA]
 *
 * Reads the certificate's JSON text {"values":[v_0,...,v_{n-1}]} with its own scanner (no JSON library: the exact
 * prefix {"values":[ , signed decimal integers separated by commas, the exact suffix ]} , whitespace allowed between
 * tokens; anything else is an error), then for every
 * knot m = 0 .. 2n-2 forms c_m = sum_{i+j=m} v_i v_j directly in signed 128-bit integers (each pair i < j once,
 * doubled, plus the square at i = j), keeping max |c_m|, every m where it is attained, and the fingerprint
 * sum_m c_m 3^m mod 2^61 - 1. No product of big integers and no transform. Before the loops it checks that
 * sum_i v_i^2 < 2^125: by Cauchy-Schwarz every |c_m| and every partial sum of the loop (at most half of sum_i v_i^2
 * before the doubling) is then below 2^125, so nothing can overflow (n max|v|^2, the guard used for #371, is above
 * 2^127 here). Output: plain "key value" lines read by c1d372_b.py.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

typedef __int128 i128;
typedef unsigned __int128 u128;
static const unsigned long long P61 = (1ULL << 61) - 1;

static void put128(const char *key, i128 x) {
    char buf[64]; int k = 63; int neg = x < 0; u128 u = neg ? (u128)(-x) : (u128)x;
    buf[k] = 0;
    do { buf[--k] = (char)('0' + (int)(u % 10)); u /= 10; } while (u);
    if (neg) buf[--k] = '-';
    printf("%s %s\n", key, buf + k);
}

static unsigned long long mod61(i128 x) {
    i128 r = x % (i128)P61;
    if (r < 0) r += P61;
    return (unsigned long long)r;
}

int main(int argc, char **argv) {
    if (argc != 2 && argc != 4) { fprintf(stderr, "usage: c1d372_b FILE [I D]\n"); return 2; }
    FILE *fh = fopen(argv[1], "rb");
    if (!fh) { perror("open"); return 2; }
    size_t cap = 1 << 16, n = 0, tlen = 0, tcap = 1 << 20;
    long long *v = malloc(cap * sizeof *v);
    char *t = malloc(tcap);
    size_t got;
    while ((got = fread(t + tlen, 1, tcap - tlen, fh)) > 0) {
        tlen += got;
        if (tlen == tcap) { tcap *= 2; t = realloc(t, tcap); }
    }
    fclose(fh);
    size_t p = 0;
#define SKIPWS while (p < tlen && (t[p] == ' ' || t[p] == '\n' || t[p] == '\r' || t[p] == '\t')) p++
    const char *pre = "{\"values\":[";
    for (const char *q = pre; *q; q++) {
        SKIPWS;
        if (p >= tlen || t[p] != *q) { fprintf(stderr, "unexpected text before the heights\n"); return 2; }
        p++;
    }
    for (;;) {
        SKIPWS;
        int neg = 0;
        if (p < tlen && t[p] == '-') { neg = 1; p++; }
        size_t d0 = p;
        unsigned long long u = 0;
        while (p < tlen && t[p] >= '0' && t[p] <= '9') {
            if (u > (1ULL << 62) / 10) { fprintf(stderr, "height out of range\n"); return 2; }
            u = u * 10 + (unsigned long long)(t[p] - '0');
            p++;
        }
        if (p == d0) { fprintf(stderr, "not an integer at byte %zu\n", p); return 2; }
        if (u > (1ULL << 62)) { fprintf(stderr, "height out of range\n"); return 2; }
        if (n == cap) { cap *= 2; v = realloc(v, cap * sizeof *v); }
        v[n++] = neg ? -(long long)u : (long long)u;
        SKIPWS;
        if (p < tlen && t[p] == ',') { p++; continue; }
        if (p < tlen && t[p] == ']') { p++; break; }
        fprintf(stderr, "unexpected character at byte %zu\n", p); return 2;
    }
    SKIPWS;
    if (p >= tlen || t[p] != '}') { fprintf(stderr, "unexpected text after the heights\n"); return 2; }
    p++;
    SKIPWS;
    if (p != tlen) { fprintf(stderr, "trailing text\n"); return 2; }
    free(t);
    if (argc == 4) {
        long long i = atoll(argv[2]), d = atoll(argv[3]);
        if (i < 0 || (size_t)i >= n) { fprintf(stderr, "bad index\n"); return 2; }
        v[i] += d;
        printf("control height %lld changed by %lld\n", i, d);
    }
    i128 S = 0; unsigned long long vmax = 0; size_t neg = 0;
    for (size_t i = 0; i < n; i++) {
        S += v[i];
        unsigned long long a = v[i] < 0 ? (unsigned long long)(-v[i]) : (unsigned long long)v[i];
        if (a > vmax) vmax = a;
        if (v[i] < 0) neg++;
    }
    /* overflow guard: sum_i v_i^2 < 2^125 (each square < 2^124 since |v_i| <= 2^62; the running sum is checked
       before every addition, so it cannot wrap) */
    u128 ss = 0, lim = ((u128)1) << 125;
    for (size_t i = 0; i < n; i++) {
        unsigned long long a = v[i] < 0 ? (unsigned long long)(-v[i]) : (unsigned long long)v[i];
        u128 sq = (u128)a * a;
        if (sq >= lim - ss) { fprintf(stderr, "sizes too large for 128 bits\n"); return 2; }
        ss += sq;
    }
    printf("n %zu\nnegative %zu\nvmax %llu\n", n, neg, vmax);
    put128("S", S);
    i128 best = -1, bestval = 0; size_t nbest = 0; size_t *where = malloc(64 * sizeof *where);
    unsigned long long fp = 0, p3 = 1;
    for (size_t m = 0; m + 1 < 2 * n; m++) {
        size_t lo = m >= n ? m - n + 1 : 0;
        i128 acc = 0;
        size_t i = lo, j = m - lo;
        for (; i < j; i++, j--) acc += (i128)v[i] * v[j];
        acc *= 2;
        if (i == j) acc += (i128)v[i] * v[i];
        fp = (unsigned long long)(((u128)fp + (u128)mod61(acc) * p3) % P61);
        p3 = (unsigned long long)((u128)p3 * 3 % P61);
        i128 a = acc < 0 ? -acc : acc;
        if (a > best) { best = a; bestval = acc; nbest = 0; }
        if (a == best && nbest < 64) where[nbest++] = m;
    }
    put128("maxabs", best);
    put128("valueatfirstmax", bestval);
    printf("argmax");
    for (size_t k = 0; k < nbest; k++) printf(" %zu", where[k]);
    printf("\nfingerprint %llu\n", fp);
    return 0;
}
