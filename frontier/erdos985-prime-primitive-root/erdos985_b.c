/* Program B (C, exact 64-bit integers) for the finite part of Erdos problem 985 as claimed by
 * Howe, Huffman, Kohlmeyer, Kruer (github.com/H0W4RD-dev/erdos-985-prime-primitive-root @ 99fd0a7).
 * Written separately from program A (no shared code), without reading the claimants' verifier.
 *
 * (a) every odd prime p <= 2791 (by trial division): the multiplicative order of each prime q
 *     is found by walking the powers q, q^2, ... until 1 (no factorisation used); the least prime
 *     q of order p - 1 must exist and be < p; the table's q must be prime, < p and of order p - 1.
 * (b) companion certificate: one segmented sieve over n in [2792, 20443^2) gives omega(n - 1)
 *     (by dividing out every prime <= 20443) and the primality of n (Eratosthenes on the same
 *     segment); n is a rank-k candidate when n is odd, omega(n-1) = k (3 <= k <= 9) and n < u_k^2;
 *     the candidates of each rank are matched in increasing order with the words of strip k-3;
 *     for prime n the word must be type 2 with b prime, b < n and b^((n-1)/r) != 1 for every prime
 *     r | n - 1; the least prime primitive root of n is also computed (same exponent test);
 *     composite n must not carry a type-2 word.
 * (c) extension: every odd prime p <= BOUND has a prime primitive root < p (bit sieve, trial
 *     division of p - 1, exponent test).
 *
 * Usage: erdos985_b SMALL_TSV STRIP_DIR BOUND [control]   control: witness | strip
 * Prints 'B: VERDICT PASS' (exit 0) or 'B: VERDICT REFUSED' (exit 1).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef uint64_t u64;
static int fails = 0;
#define FAIL(...) do { if (fails < 20) { printf("B: FAIL "); printf(__VA_ARGS__); printf("\n"); } fails++; } while (0)

static u64 mulmod(u64 a, u64 b, u64 m) { return (u64)((unsigned __int128)a * b % m); }
static u64 powmod(u64 a, u64 e, u64 m) { u64 r = 1 % m; a %= m; while (e) { if (e & 1) r = mulmod(r, a, m); a = mulmod(a, a, m); e >>= 1; } return r; }

static int *pr; static int npr;          /* primes up to 1<<16 */
static int isprime_small(u64 n) {        /* trial division, n < 2^32 */
  if (n < 2) return 0;
  for (int i = 0; i < npr; i++) { u64 r = pr[i]; if (r * r > n) return 1; if (n % r == 0) return n == r; }
  fprintf(stderr, "prime list too short\n"); exit(2);
}
static int distinct(u64 m, u64 *rs) {     /* distinct prime factors by trial division */
  int c = 0;
  for (int i = 0; i < npr && (u64)pr[i] * pr[i] <= m; i++) if (m % pr[i] == 0) { rs[c++] = pr[i]; while (m % pr[i] == 0) m /= pr[i]; }
  if (m > 1) rs[c++] = m;
  return c;
}
static int primroot(u64 q, u64 p, const u64 *rs, int c) {
  if (q % p == 0) return 0;
  for (int i = 0; i < c; i++) if (powmod(q, (p - 1) / rs[i], p) == 1) return 0;
  return 1;
}
static u64 least_root(u64 p, const u64 *rs, int c) {
  for (int i = 0; i < npr; i++) { u64 q = pr[i]; if (q >= p) return 0; if (primroot(q, p, rs, c)) return q; }
  return 0;
}
static u64 order_walk(u64 q, u64 p) {     /* multiplicative order by walking powers */
  if (q % p == 0) return 0;
  u64 x = q % p, k = 1;
  while (x != 1) { x = x * q % p; k++; if (k > p) return 0; }
  return k;
}

static const u64 UP[10] = {0, 0, 0, 143ULL * 143, 363ULL * 363, 871ULL * 871, 2019ULL * 2019, 4321ULL * 4321, 9693ULL * 9693, 20443ULL * 20443};

int main(int argc, char **argv) {
  if (argc < 4) { fprintf(stderr, "usage: %s SMALL_TSV STRIP_DIR BOUND [witness|strip]\n", argv[0]); return 2; }
  const char *control = argc > 4 ? argv[4] : "";
  u64 bound = strtoull(argv[3], 0, 10);
  /* small primes */
  int L = 1 << 16; char *c = calloc(L + 1, 1); pr = malloc(sizeof(int) * L); npr = 0;
  for (int i = 2; i <= L; i++) if (!c[i]) { pr[npr++] = i; for (long j = (long)i * i; j <= L; j += i) c[j] = 1; }

  /* (a) */
  FILE *f = fopen(argv[1], "r"); if (!f) { perror("tsv"); return 2; }
  char line[4096]; int nrows = 0; static int tq[3000]; memset(tq, 0, sizeof tq);
  if (!fgets(line, sizeof line, f) || strncmp(line, "p\tq\t", 4)) { printf("B: bad header\n"); return 2; }
  while (fgets(line, sizeof line, f)) {
    long p, q; if (sscanf(line, "%ld\t%ld", &p, &q) != 2) continue;
    if (p < 3 || p > 2791) { FAIL("(a) row with p=%ld outside 3..2791", p); continue; }
    if (tq[p]) FAIL("(a) duplicate row p=%ld", p);
    tq[p] = (int)q; nrows++;
  }
  fclose(f);
  if (!strcmp(control, "witness")) { tq[7] = 2; printf("B: CONTROL witness: table row p=7 given q=2 (order 3 mod 7)\n"); }
  int nsmall = 0, nleast = 0, noroot = 0; u64 gmax = 0, gmaxp = 0;
  for (u64 p = 3; p <= 2791; p += 2) {
    if (!isprime_small(p)) { if (tq[p]) FAIL("(a) row for composite p=%llu", (unsigned long long)p); continue; }
    nsmall++;
    u64 g = 0;
    for (int i = 0; i < npr && (u64)pr[i] < p; i++) if (order_walk(pr[i], p) == p - 1) { g = pr[i]; break; }
    if (!g) { noroot++; FAIL("(a) no prime primitive root below p=%llu", (unsigned long long)p); continue; }
    if (g > gmax) { gmax = g; gmaxp = p; }
    u64 q = tq[p];
    if (!q) { FAIL("(a) no table row for p=%llu", (unsigned long long)p); continue; }
    if (!(q < p && isprime_small(q) && order_walk(q, p) == p - 1)) FAIL("(a) table row p=%llu: q=%llu is not a prime of order p-1 below p", (unsigned long long)p, (unsigned long long)q);
    if (q == g) nleast++;
  }
  printf("B (a): odd primes p <= 2791: %d; table rows %d; least prime primitive root < p for all: %s; table q is the least one in %d rows; largest least prime primitive root %llu (p = %llu)\n",
         nsmall, nrows, noroot ? "NO" : "yes", nleast, (unsigned long long)gmax, (unsigned long long)gmaxp);

  /* (b) */
  uint16_t *W[10]; long nw[10], pos[10];
  for (int k = 3; k <= 9; k++) {
    char path[4096]; snprintf(path, sizeof path, "%s/strip-%02d.u16", argv[2], k - 3);
    FILE *g = fopen(path, "rb"); if (!g) { perror(path); return 2; }
    fseek(g, 0, SEEK_END); long sz = ftell(g); fseek(g, 0, SEEK_SET);
    unsigned char *buf = malloc(sz); if (fread(buf, 1, sz, g) != (size_t)sz) return 2; fclose(g);
    nw[k] = sz / 2; W[k] = malloc(sizeof(uint16_t) * (nw[k] + 1)); pos[k] = 0;
    for (long i = 0; i < nw[k]; i++) W[k][i] = (uint16_t)(buf[2 * i] | (buf[2 * i + 1] << 8));
    free(buf);
  }
  int strip_control_done = 0;
  long np[10] = {0}, ty[10][4]; memset(ty, 0, sizeof ty);
  u64 bmax = 0, gbm = 0, gbn = 0;
  const u64 LO = 2792, HI = UP[9];
  const u64 S = 1 << 20;
  uint32_t *rem = malloc(sizeof(uint32_t) * S); unsigned char *om = malloc(S); unsigned char *comp = malloc(S);
  int lim = 0; while (lim < npr && (u64)pr[lim] * pr[lim] < HI) lim++;   /* primes r with r^2 < HI */
  for (u64 base = LO; base < HI; base += S) {
    u64 top = base + S < HI ? base + S : HI, len = top - base;
    /* N = n - 1 for n in [base, top): N in [base-1, top-1) */
    for (u64 i = 0; i < len; i++) { rem[i] = (uint32_t)(base - 1 + i); om[i] = 0; comp[i] = 0; }
    for (int t = 0; t < lim; t++) {
      u64 r = pr[t];
      u64 st = ((base - 1 + r - 1) / r) * r;               /* first multiple of r that is >= base-1 */
      for (u64 m = st; m < top - 1; m += r) { u64 i = m - (base - 1); om[i]++; uint32_t x = rem[i] / r; while (x % r == 0) x /= r; rem[i] = x; }
      u64 s2 = ((base + r - 1) / r) * r; if (s2 < r * r) s2 = r * r;   /* composite n */
      for (u64 m = s2; m < top; m += r) comp[m - base] = 1;
    }
    for (u64 i = 0; i < len; i++) {
      u64 n = base + i;
      if (!(n & 1)) continue;
      int k = om[i] + (rem[i] > 1);
      if (k < 3 || k > 9 || n >= UP[k]) continue;
      if (pos[k] >= nw[k]) { FAIL("(b) rank %d: more candidates than words", k); continue; }
      long idx = pos[k]++;
      uint16_t wd = W[k][idx];
      int a = wd >> 14; u64 b = wd & 0x3FFF;
      int isp = !comp[i];
      if (isp && !strcmp(control, "strip") && k == 3 && !strip_control_done) {
        u64 rs0[16]; int c0 = distinct(n - 1, rs0);
        for (int j = 0; j < npr; j++) if (!primroot(pr[j], n, rs0, c0)) { printf("B: CONTROL strip: rank 3, n=%llu: witness %llu replaced by the prime %d (not a primitive root)\n", (unsigned long long)n, (unsigned long long)b, pr[j]); b = pr[j]; a = 2; break; }
        strip_control_done = 1;
      }
      ty[k][a]++;
      if (isp) {
        np[k]++;
        u64 rs[16]; int cc = distinct(n - 1, rs);
        if (cc != k) FAIL("(b) n=%llu: omega mismatch", (unsigned long long)n);
        u64 g = least_root(n, rs, cc);
        if (!g) FAIL("(b) no prime primitive root below n=%llu", (unsigned long long)n);
        else if (g > gbm) { gbm = g; gbn = n; }
        if (!(a == 2 && b > 1 && b < n && isprime_small(b) && primroot(b, n, rs, cc)))
          FAIL("(b) prime n=%llu: word (%d, %llu) is not a prime primitive root", (unsigned long long)n, a, (unsigned long long)b);
        if (b > bmax) bmax = b;
      } else if (a == 2) FAIL("(b) composite n=%llu carries a type-2 word", (unsigned long long)n);
    }
  }
  long tc = 0, tp = 0, t0 = 0, t1 = 0, t2 = 0;
  for (int k = 3; k <= 9; k++) {
    if (pos[k] != nw[k]) FAIL("(b) rank %d: %ld candidates but %ld words", k, pos[k], nw[k]);
    printf("B (b): rank %d, 2792 <= n < %llu: candidates %ld, words %ld, prime moduli %ld, word types 0/1/2 = %ld/%ld/%ld\n",
           k, (unsigned long long)UP[k], pos[k], nw[k], np[k], ty[k][0], ty[k][1], ty[k][2]);
    tc += pos[k]; tp += np[k]; t0 += ty[k][0]; t1 += ty[k][1]; t2 += ty[k][2];
  }
  printf("B (b): total candidates %ld, prime moduli %ld, word types 0/1/2 = %ld/%ld/%ld, largest certificate witness %llu, largest least prime primitive root %llu (n = %llu)\n",
         tc, tp, t0, t1, t2, (unsigned long long)bmax, (unsigned long long)gbm, (unsigned long long)gbn);

  /* (c) */
  u64 nb = bound / 2 + 1; unsigned char *bits = calloc(nb / 8 + 1, 1);   /* bit i: 2i+1 composite */
  for (u64 i = 1; (2 * i + 1) * (2 * i + 1) <= bound; i++) if (!(bits[i >> 3] >> (i & 7) & 1)) { u64 p = 2 * i + 1; for (u64 j = p * p; j <= bound; j += 2 * p) bits[(j / 2) >> 3] |= 1 << ((j / 2) & 7); }
  long cnt = 0, bad = 0; u64 gs = 0, gm = 0, gmp = 0;
  for (u64 p = 3; p <= bound; p += 2) {
    u64 i = p / 2; if (bits[i >> 3] >> (i & 7) & 1) continue;
    cnt++;
    u64 rs[16]; int cc = distinct(p - 1, rs);
    u64 g = least_root(p, rs, cc);
    if (!g) { bad++; FAIL("(c) no prime primitive root below p=%llu", (unsigned long long)p); continue; }
    gs += g; if (g > gm) { gm = g; gmp = p; }
  }
  printf("B (c): odd primes p <= %llu: %ld; each has a prime primitive root < p: %s; sum of least prime primitive roots %llu; largest %llu (first at p = %llu)\n",
         (unsigned long long)bound, cnt, bad ? "NO" : "yes", (unsigned long long)gs, (unsigned long long)gm, (unsigned long long)gmp);
  if (fails) { printf("B: VERDICT REFUSED (%d failures)\n", fails); return 1; }
  printf("B: VERDICT PASS\n");
  return 0;
}
