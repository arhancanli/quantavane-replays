/* c1c_ntt.c: the exact autocorrelation minimum of a C_1c step-function certificate with millions of cells.
   Written for our replays of PR #354's K = 2^22 and K = 2^24 certificates (no code of the certificate's repository
   was read or run).

   Input: a text file "K\n w_0 w_1 ... w_K\n" (non-negative decimal integers below 2^63, separated by blanks; the
   format of the certificate files). For f = sum_i w_i 1_[i/K, (i+1)/K) the bound is K min_{0<=k<=K} A_k / S^2 with
   A_k = sum_i w_i w_{i+k} and S = sum_i w_i (the 1c page's exact rule).

   Method: the cyclic autocorrelation of length L (the least power of two >= 2K) modulo three primes p below 2^62
   with 2^27 | p - 1, by a number-theoretic transform: C_j = W_j W_{L-j} (W = NTT(w)), then the inverse transform.
   With L >= 2K no lag 0..K wraps around, except when L = 2K, where lag K picks up the one extra term w_K w_0, which
   is subtracted. Every A_k is at most A_0 = sum w_i^2 (Cauchy-Schwarz), which we sum exactly and require to be below
   p1 p2, so A_k is the CRT value from p1 and p2; the third prime is an identity check on every lag (the CRT value
   reduced mod p3 must equal the third transform's value). The minimum and its lags are found exactly in 128 bits,
   and every minimizing lag, lag 0, lag K and 64 random lags are summed again directly.

   Usage: c1c_ntt FILE [--bump J] [--scale C] [--fault]   (--bump adds 1 to w_J and --scale multiplies every weight
          by C, both after reading, for red controls; --fault corrupts one residue of p1 to show that the identity
          catches it)
          c1c_ntt --selftest   (random small cases against the direct double sum, L = 2K and L > 2K both)
   Output: lines "key value"; exit 0 on success, 1 on any failed check. */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef uint64_t u64;
typedef unsigned __int128 u128;

static const u64 P[3] = {4611686009971671041ULL, 4611686007555751937ULL, 4611686004066091009ULL};

static u64 mulmod(u64 a, u64 b, u64 p) { return (u64)(((u128)a * b) % p); }
static u64 powmod(u64 a, u64 e, u64 p) {
    u64 r = 1 % p;
    a %= p;
    while (e) { if (e & 1) r = mulmod(r, a, p); a = mulmod(a, a, p); e >>= 1; }
    return r;
}
static int is_prime(u64 n) {  /* deterministic Miller-Rabin for 64-bit n */
    static const u64 B[12] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    if (n < 2) return 0;
    for (int i = 0; i < 12; i++) if (n % B[i] == 0) return n == B[i];
    u64 d = n - 1; int s = 0;
    while (!(d & 1)) { d >>= 1; s++; }
    for (int i = 0; i < 12; i++) {
        u64 x = powmod(B[i], d, n);
        if (x == 1 || x == n - 1) continue;
        int ok = 0;
        for (int r = 1; r < s; r++) { x = mulmod(x, x, n); if (x == n - 1) { ok = 1; break; } }
        if (!ok) return 0;
    }
    return 1;
}
static void u128_str(u128 x, char *out) {
    char buf[64]; int n = 0;
    if (x == 0) { strcpy(out, "0"); return; }
    while (x) { buf[n++] = (char)('0' + (int)(x % 10)); x /= 10; }
    for (int i = 0; i < n; i++) out[i] = buf[n - 1 - i];
    out[n] = 0;
}

/* forward transform of a (length n, a power of two) with tw[j] = omega^j, j < n/2, omega of order n */
static void ntt(u64 *a, size_t n, const u64 *tw, u64 p) {
    for (size_t i = 1, j = 0; i < n; i++) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) { u64 t = a[i]; a[i] = a[j]; a[j] = t; }
    }
    for (size_t len = 2; len <= n; len <<= 1) {
        size_t half = len >> 1, step = n / len;
        for (size_t i = 0; i < n; i += len) {
            for (size_t j = 0; j < half; j++) {
                u64 u = a[i + j], v = mulmod(a[i + j + half], tw[j * step], p);
                u64 s = u + v; if (s >= p) s -= p;
                a[i + j] = s;
                a[i + j + half] = (u >= v) ? u - v : u + p - v;
            }
        }
    }
}

/* cyclic autocorrelation c_k = sum_i w_i w_{(i+k) mod L}, k = 0..K, modulo p, into out[0..K] (with the L = 2K fix) */
static int autocorr_mod(const u64 *w, size_t K, size_t L, u64 p, u64 *out, u64 *a, u64 *tw) {
    u64 g = 2;
    while (powmod(g, (p - 1) / 2, p) != p - 1) g++;          /* a quadratic non-residue */
    u64 om = powmod(g, (p - 1) / L, p);                        /* of order L, since 2^27 | p - 1 and L | 2^27 */
    if (powmod(om, L / 2, p) != p - 1) return 0;
    tw[0] = 1;
    for (size_t j = 1; j < L / 2; j++) tw[j] = mulmod(tw[j - 1], om, p);
    for (size_t i = 0; i < L; i++) a[i] = (i <= K) ? w[i] % p : 0;
    ntt(a, L, tw, p);
    for (size_t j = 1; j < L / 2; j++) { u64 t = mulmod(a[j], a[L - j], p); a[j] = t; a[L - j] = t; }
    a[0] = mulmod(a[0], a[0], p);
    a[L / 2] = mulmod(a[L / 2], a[L / 2], p);
    ntt(a, L, tw, p);                                          /* inverse = forward, then reverse 1..L-1, / L */
    u64 linv = powmod(L % p, p - 2, p);
    out[0] = mulmod(a[0], linv, p);
    for (size_t k = 1; k <= K; k++) out[k] = mulmod(a[L - k], linv, p);
    if (L == 2 * K) {                                           /* lag K: remove the wrapped term w_K w_0 */
        u64 t = mulmod(w[K] % p, w[0] % p, p);
        out[K] = (out[K] >= t) ? out[K] - t : out[K] + p - t;
    }
    return 1;
}

static u128 direct(const u64 *w, size_t K, size_t k, int *ovf) {
    u128 s = 0;
    for (size_t i = 0; i + k <= K; i++) {
        u128 t = (u128)w[i] * w[i + k];
        if (s + t < s) *ovf = 1;
        s += t;
    }
    return s;
}

static u64 rng_state = 88172645463325252ULL;
static u64 rng(void) { rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7; rng_state ^= rng_state << 17; return rng_state; }

/* the whole computation: returns 0 on success and prints the results */
static int run(u64 *w, size_t K, int fault, int verbose) {
    size_t L = 1;
    while (L < 2 * K) L <<= 1;
    if (L > ((size_t)1 << 27)) { printf("error L too large\n"); return 1; }
    for (int t = 0; t < 3; t++) if (!is_prime(P[t]) || ((P[t] - 1) & ((1ULL << 27) - 1))) { printf("error prime\n"); return 1; }
    /* exact S and A_0 */
    u128 S = 0; int ovf = 0;
    for (size_t i = 0; i <= K; i++) { if (S + w[i] < S) ovf = 1; S += w[i]; }
    u128 A0 = direct(w, K, 0, &ovf);
    u128 P12 = (u128)P[0] * P[1];
    if (ovf || A0 >= P12) { printf("error A_0 does not fit below p1 p2 (overflow %d)\n", ovf); return 1; }
    u64 *a = malloc(L * sizeof(u64)), *tw = malloc((L / 2 + 1) * sizeof(u64));
    u64 *r[3];
    for (int t = 0; t < 3; t++) r[t] = malloc((K + 1) * sizeof(u64));
    if (!a || !tw || !r[0] || !r[1] || !r[2]) { printf("error memory\n"); return 1; }
    for (int t = 0; t < 3; t++) {
        if (!autocorr_mod(w, K, L, P[t], r[t], a, tw)) { printf("error root of unity\n"); return 1; }
        if (verbose) { fprintf(stderr, "prime %d done\n", t + 1); }
    }
    free(a); free(tw);
    size_t faultk = (size_t)-1;
    if (fault) { faultk = (size_t)(rng() % (K + 1)); r[0][faultk] = (r[0][faultk] + 1) % P[0]; }
    u64 inv = powmod(P[0] % P[1], P[1] - 2, P[1]);
    u128 best = 0; size_t nbest = 0, firstbest = 0, bad = 0, firstbad = 0;
    size_t *argmins = malloc(64 * sizeof(size_t));
    for (size_t k = 0; k <= K; k++) {
        u64 r1 = r[0][k], r2 = r[1][k] % P[1];
        u64 d = (r2 >= r1 % P[1]) ? r2 - r1 % P[1] : r2 + P[1] - r1 % P[1];
        u64 tt = mulmod(d, inv, P[1]);
        u128 x = (u128)r1 + (u128)P[0] * tt;                   /* the unique value below p1 p2 */
        if ((u64)(x % P[2]) != r[2][k]) { if (!bad) firstbad = k; bad++; }
        if (k == 0 || x < best) { best = x; nbest = 0; firstbest = k; }
        if (x == best) { if (nbest < 64) argmins[nbest] = k; nbest++; }
    }
    char sbuf[64];
    u128_str(S, sbuf); printf("S %s\n", sbuf);
    u128_str(A0, sbuf); printf("A0 %s\n", sbuf);
    printf("K %zu\nL %zu\n", K, L);
    printf("identity_failures %zu\n", bad);
    if (bad) printf("first_identity_failure %zu\n", firstbad);
    if (fault) printf("fault_lag %zu\n", faultk);
    u128_str(best, sbuf); printf("Amin %s\n", sbuf);
    printf("argmin_count %zu\nargmin", nbest);
    for (size_t i = 0; i < nbest && i < 64; i++) printf(" %zu", argmins[i]);
    printf("\n");
    /* direct sums: every minimizing lag (up to 64), lags 0 and K, and 64 random lags */
    int dok = 1; ovf = 0;
    for (size_t i = 0; i < nbest && i < 64; i++) if (direct(w, K, argmins[i], &ovf) != best) dok = 0;
    size_t checked = (nbest < 64 ? nbest : 64);
    size_t lags[66]; int nl = 0;
    lags[nl++] = 0; lags[nl++] = K;
    for (int i = 0; i < 64; i++) lags[nl++] = (size_t)(rng() % (K + 1));
    for (int i = 0; i < nl; i++) {
        size_t k = lags[i];
        u128 x = direct(w, K, k, &ovf);
        u64 r1 = r[0][k], r2 = r[1][k];
        u64 d = (r2 >= r1 % P[1]) ? r2 - r1 % P[1] : r2 + P[1] - r1 % P[1];
        u128 y = (u128)r1 + (u128)P[0] * mulmod(d, inv, P[1]);
        if (x != y) dok = 0;
        checked++;
    }
    printf("direct_sums_checked %zu\ndirect_sums_agree %d\n", checked, dok && !ovf);
    (void)firstbest;
    for (int t = 0; t < 3; t++) free(r[t]);
    free(argmins);
    return (bad || !dok || ovf) ? 1 : 0;
}

static int selftest(void) {
    int ok = 1, cases = 0;
    for (int c = 0; c < 200; c++) {
        size_t K = (c % 4 == 0) ? ((size_t)1 << (1 + rng() % 11)) : (size_t)(1 + rng() % 3000);
        u64 *w = malloc((K + 1) * sizeof(u64));
        for (size_t i = 0; i <= K; i++) {
            u64 m = rng() % 4;
            w[i] = (m == 0) ? 0 : (m == 1 ? rng() % 1000 : (m == 2 ? rng() >> 12 : (rng() >> 30)));
        }
        if (c % 7 == 0) w[0] = (u64)1 << 52;
        size_t L = 1; while (L < 2 * K) L <<= 1;
        u64 *a = malloc(L * sizeof(u64)), *tw = malloc((L / 2 + 1) * sizeof(u64)), *r[3];
        for (int t = 0; t < 3; t++) { r[t] = malloc((K + 1) * sizeof(u64)); autocorr_mod(w, K, L, P[t], r[t], a, tw); }
        u64 inv = powmod(P[0] % P[1], P[1] - 2, P[1]);
        int ovf = 0;
        for (size_t k = 0; k <= K; k++) {
            u128 x = direct(w, K, k, &ovf);
            u64 r1 = r[0][k], r2 = r[1][k];
            u64 d = (r2 >= r1 % P[1]) ? r2 - r1 % P[1] : r2 + P[1] - r1 % P[1];
            u128 y = (u128)r1 + (u128)P[0] * mulmod(d, inv, P[1]);
            if (x != y || (u64)(x % P[2]) != r[2][k]) ok = 0;
        }
        cases++;
        free(a); free(tw); free(w); for (int t = 0; t < 3; t++) free(r[t]);
    }
    printf("selftest %d cases (K up to 3000 and powers of two up to 2048, so L = 2K and L > 2K): %s\n",
           cases, ok ? "all equal to the direct double sum" : "MISMATCH");
    return ok ? 0 : 1;
}

int main(int argc, char **argv) {
    if (argc >= 2 && !strcmp(argv[1], "--selftest")) return selftest();
    if (argc < 2) { fprintf(stderr, "usage: c1c_ntt FILE [--bump J] [--scale C] [--fault] | --selftest\n"); return 2; }
    long long bump = -1; u64 scale = 1; int fault = 0;
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--bump") && i + 1 < argc) bump = atoll(argv[++i]);
        else if (!strcmp(argv[i], "--scale") && i + 1 < argc) scale = strtoull(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "--fault")) fault = 1;
        else { fprintf(stderr, "unknown argument %s\n", argv[i]); return 2; }
    }
    FILE *f = fopen(argv[1], "rb");
    if (!f) { printf("error cannot open\n"); return 1; }
    /* strict parser: tokens of decimal digits separated by spaces, tabs, CR or LF; first token K, then K+1 weights */
    size_t K = 0, cap = 0, n = 0; u64 *w = NULL; int first = 1, intok = 0, toolong = 0; u64 cur = 0;
    int ch;
    for (;;) {
        ch = getc_unlocked(f);
        if (ch >= '0' && ch <= '9') {
            if (cur > (((u64)1 << 63) - 1 - (u64)(ch - '0')) / 10) toolong = 1;
            cur = cur * 10 + (u64)(ch - '0'); intok = 1;
            continue;
        }
        if (ch == ' ' || ch == '\t' || ch == '\r' || ch == '\n' || ch == EOF) {
            if (intok) {
                if (toolong) { printf("error a number is not below 2^63\n"); return 1; }
                if (first) {
                    K = (size_t)cur; first = 0; cap = K + 1;
                    if (K < 1 || K > ((size_t)1 << 26)) { printf("error K out of range\n"); return 1; }
                    w = malloc(cap * sizeof(u64));
                    if (!w) { printf("error memory\n"); return 1; }
                } else {
                    if (n >= cap) { printf("error more than K+1 weights\n"); return 1; }
                    w[n++] = cur;
                }
                cur = 0; intok = 0;
            }
            if (ch == EOF) break;
            continue;
        }
        printf("error unexpected byte %d\n", ch);
        return 1;
    }
    fclose(f);
    if (first || n != K + 1) { printf("error expected K+1 = %zu weights, read %zu\n", K + 1, n); return 1; }
    if (bump >= 0) { if ((size_t)bump > K) { printf("error bump index\n"); return 1; } w[bump] += 1; }
    if (scale != 1) for (size_t i = 0; i <= K; i++) w[i] *= scale;   /* only for the controls; checked below */
    u64 wmax = 0; size_t nz = 0;
    for (size_t i = 0; i <= K; i++) { if (w[i] > wmax) wmax = w[i]; if (w[i]) nz++; }
    printf("weights %zu\nnonzero %zu\nwmax %llu\n", n, nz, (unsigned long long)wmax);
    if (wmax >= ((u64)1 << 62)) { printf("error a weight is not below 2^62 after scaling\n"); return 1; }
    return run(w, K, fault, 1);
}
