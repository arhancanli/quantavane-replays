/* 14a, program A: the determinations behind the rows C_14 >= 4 [LS1965] (BB(3)) and C_14 >= 5 [B1983] (BB(4)).
 *
 * Every n-state, 2-symbol Turing machine (n = 2, 3, 4), started in state A on a blank tape, is enumerated in tree
 * normal form: a transition is defined only when the simulation first reaches it, with the next state at most one
 * more than the largest state used so far, and the first move fixed to the right (mirror symmetry). When the
 * simulation reaches an undefined transition the machine with that transition as its halting transition is
 * recorded (it halts after that step, the halting transition writing 1); the tree then continues with every
 * definition of it, as long as one transition stays undefined. A machine that runs past the step limit without
 * reaching an undefined transition must be proved never to halt: here by the n-gram closed position set test
 * (radius 1 to 6): the set of local views (left n cells, state, head symbol, right n cells) reachable from the blank
 * tape, with every n consecutive cells on either side drawn from a closed set of n-grams, is computed to a fixed
 * point; if it contains no undefined transition, no halting configuration is reachable. A machine that no radius
 * decides is printed as a holdout. Output: S(n) (most steps, the halting step included) and Sigma(n) (most 1s).
 * Usage: c14a_a N [LIMIT]   (red controls: c14a_a N LIMIT --maxrad R)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 5
static int N, LIMIT, MAXRAD = 12;
static int W[MAXN][2], D[MAXN][2], NS[MAXN][2];   /* NS = -1: undefined */
static long long n_halt = 0, n_nonhalt = 0, n_holdout = 0, by_rad[16], by_depth[16];
static int HDEPTH = 9, SELFTEST = 0;
static long long n_selftest_fail = 0;
static int HISTRAD = 4;
static long long bestS = 0, bestSigma = 0;
static unsigned char *tape;
static int TAPELEN;

/* ---------------- n-gram closed position set ----------------
 * Alphabet of K tape symbols; the machine reads bit = sym / T and writes w * T + tag, where T = 1 (plain test,
 * K = 2) or T = N + 1 and tag = 1 + the state that wrote the cell (the history variant, K = 2(N + 1); the blank is
 * 0). n-grams are base-K integers, digit 0 next to the head. Local views (L, s, c, R) and the two n-gram sets are
 * grown to a fixed point. */
typedef unsigned long long u64;
typedef unsigned __int128 u128;
typedef struct { u128 *key; unsigned *gen; unsigned g; size_t cap, n; } hset;   /* a slot is empty unless gen == g */
static void hs_init(hset *h, size_t cap) { h->cap = cap; h->n = 0; h->g = 1; h->key = malloc(cap * sizeof(u128)); h->gen = calloc(cap, sizeof(unsigned)); }
static void hs_clear(hset *h) { h->g++; h->n = 0; }
static size_t hs_slot(hset *h, u128 k) {
  u64 z = (u64)k * 0x9E3779B97F4A7C15ULL ^ (u64)(k >> 64) * 0xC2B2AE3D27D4EB4FULL;
  size_t i = (size_t)(z ^ (z >> 29)) & (h->cap - 1);
  while (h->gen[i] == h->g && h->key[i] != k) i = (i + 1) & (h->cap - 1);
  return i;
}
static int hs_has(hset *h, u128 k) { return h->gen[hs_slot(h, k)] == h->g; }
static int hs_add(hset *h, u128 k) {   /* 1 if new; 0 if present; -1 if full */
  size_t i = hs_slot(h, k);
  if (h->gen[i] == h->g) return 0;
  if (h->n * 2 >= h->cap) return -1;
  h->key[i] = k; h->gen[i] = h->g; h->n++; return 1;
}
static hset LS, RS, LOC;
static u128 *locv; static size_t nlocv, caplocv = 1 << 22;
static int ngram_cps(int r, int hist) {
  u64 T = 1, K, P = 1;
  for (int i = 0; i < hist; i++) T *= (u64)(N + 1);    /* the last `hist` writer states of each cell */
  K = 2 * T;
  for (int i = 0; i < r; i++) { if (P > (1ULL << 56) / K) return 0; P *= K; }   /* P = K^r <= 2^56, checked before each product */
  if ((u128)P * P > (((u128)1) << 126) / ((u128)N * K)) return 0;   /* every view key ((L N + s) K + c) P + R < 2^126 */
  hs_clear(&LS); hs_clear(&RS); hs_clear(&LOC); nlocv = 0;
  hs_add(&LS, 0); hs_add(&RS, 0);
  u128 k0 = 0; hs_add(&LOC, k0); locv[nlocv++] = k0;   /* key = ((L * N + s) * K + c) * P + R */
  int changed = 1;
  while (changed) {
    changed = 0;
    for (size_t i = 0; i < nlocv; i++) {
      u128 key = locv[i];
      u64 R = (u64)(key % P); u128 t = key / P; u64 c = (u64)(t % K); t /= K; int s = (int)(t % N); u64 L = (u64)(t / N);
      int bit = (int)(c / T);
      if (NS[s][bit] < 0) return 0;
      u64 wsym = (u64)W[s][bit] * T + (hist == 0 ? 0 : ((c % T) * (N + 1) + (u64)(s + 1)) % T);
      int s2 = NS[s][bit];
      u64 news[2]; int nn = 0;
      if (D[s][bit]) {
        u64 L2 = (u64)(((u128)L * K + wsym) % P);
        int a = hs_add(&LS, L2); if (a < 0) return 0; if (a) changed = 1;
        u64 c2 = R % K;
        for (u64 x = 0; x < K; x++) {
          u64 R2 = R / K + x * (P / K);
          if (!hs_has(&RS, R2)) continue;
          u128 k = (((u128)L2 * N + s2) * K + c2) * P + R2;
          int b = hs_add(&LOC, k); if (b < 0) return 0;
          if (b) { if (nlocv == caplocv) return 0; locv[nlocv++] = k; changed = 1; }
        }
      } else {
        u64 R2 = (u64)(((u128)R * K + wsym) % P);
        int a = hs_add(&RS, R2); if (a < 0) return 0; if (a) changed = 1;
        u64 c2 = L % K;
        for (u64 x = 0; x < K; x++) {
          u64 L2 = L / K + x * (P / K);
          if (!hs_has(&LS, L2)) continue;
          u128 k = (((u128)L2 * N + s2) * K + c2) * P + R2;
          int b = hs_add(&LOC, k); if (b < 0) return 0;
          if (b) { if (nlocv == caplocv) return 0; locv[nlocv++] = k; changed = 1; }
        }
      }
      (void)news; (void)nn;
    }
  }
  return 1;
}


static void print_machine(FILE *f) {
  for (int s = 0; s < N; s++) {
    if (s) fputc('_', f);
    for (int c = 0; c < 2; c++) {
      if (NS[s][c] < 0) fputs("---", f);
      else fprintf(f, "%d%c%c", W[s][c], D[s][c] ? 'R' : 'L', 'A' + NS[s][c]);
    }
  }
}

/* simulate from blank; return -1 if past the limit, else the step at which an undefined transition is read */
static int simulate(int *hs, int *hc, long long *ones, int *maxstate) {
  memset(tape, 0, TAPELEN);
  int h = TAPELEN / 2, s = 0, steps = 0, ms = 0;
  while (steps < LIMIT) {
    int c = tape[h];
    if (NS[s][c] < 0) {
      long long o = 0;
      for (int i = 0; i < TAPELEN; i++) o += tape[i];
      if (c == 0) o += 1;                         /* the halting transition writes 1 */
      *hs = s; *hc = c; *ones = o; *maxstate = ms;
      return steps;
    }
    tape[h] = (unsigned char)W[s][c];
    h += D[s][c] ? 1 : -1;
    s = NS[s][c];
    if (s > ms) ms = s;
    steps++;
    if (h <= 0 || h >= TAPELEN - 1) { fprintf(stderr, "tape overflow\n"); exit(2); }
  }
  *maxstate = ms;
  return -1;
}

static void explore(int defined) {
  int hs, hc, ms; long long ones;
  int t = simulate(&hs, &hc, &ones, &ms);
  if (t < 0) {                                    /* must be proved non-halting */
    for (int r = 1; r <= MAXRAD; r++)
      if (ngram_cps(r, 0)) { n_nonhalt++; by_rad[r]++; return; }
    for (int h = 1; h <= HDEPTH; h++)
      for (int r = 1; r <= HISTRAD; r++)
        if (ngram_cps(r, h)) { n_nonhalt++; by_depth[h]++; return; }
    n_holdout++;
    printf("HOLDOUT "); print_machine(stdout); printf("\n");
    return;
  }
  /* the machine with (hs, hc) as halting transition halts after step t + 1 */
  n_halt++;
  if (SELFTEST) {                                 /* no test may claim that this halting machine never halts */
    for (int r = 1; r <= MAXRAD; r++) if (ngram_cps(r, 0)) n_selftest_fail++;
    for (int h = 1; h <= HDEPTH; h++) for (int r = 1; r <= HISTRAD; r++) if (ngram_cps(r, h)) n_selftest_fail++;
  }
  if (t + 1 > bestS) bestS = t + 1;
  if (ones > bestSigma) bestSigma = ones;
  if (defined == 2 * N - 1) return;               /* the last undefined transition must stay a halt */
  int top = ms + 1 < N ? ms + 1 : N - 1;
  for (int ns = 0; ns <= top; ns++)
    for (int w = 0; w < 2; w++)
      for (int d = 0; d < 2; d++) {
        if (defined == 0 && d == 0) continue;     /* the first move goes right (mirror symmetry) */
        W[hs][hc] = w; D[hs][hc] = d; NS[hs][hc] = ns;
        explore(defined + 1);
        NS[hs][hc] = -1;
      }
}

int main(int argc, char **argv) {
  if (argc < 2) { fprintf(stderr, "usage: c14a_a N [LIMIT] [--maxrad R]\n"); return 1; }
  N = atoi(argv[1]);
  LIMIT = argc > 2 ? atoi(argv[2]) : 5000;
  for (int i = 3; i < argc; i++) { if (!strcmp(argv[i], "--maxrad") && i + 1 < argc) MAXRAD = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--histrad") && i + 1 < argc) HISTRAD = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--hdepth") && i + 1 < argc) HDEPTH = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--selftest")) SELFTEST = 1;
    else if (!strcmp(argv[i], "--machine")) i++; }
  if (N < 1 || N > MAXN || MAXRAD < 1 || MAXRAD > 12 || HISTRAD < 0 || HISTRAD > 5 || HDEPTH < 0 || HDEPTH > 9) return 1;
  TAPELEN = 2 * LIMIT + 16;
  tape = malloc(TAPELEN);
  hs_init(&LS, 1 << 20); hs_init(&RS, 1 << 20); hs_init(&LOC, 1 << 23); locv = malloc(caplocv * sizeof(u128));
  for (int s = 0; s < MAXN; s++) for (int c = 0; c < 2; c++) NS[s][c] = -1;
  for (int i = 3; i < argc; i++) if (!strcmp(argv[i], "--machine") && i + 1 < argc) {
    /* decide one machine given as e.g. 1RB1LB_1LA0LC_---1LD_1RD0RA (--- = undefined, i.e. halting) */
    const char *m = argv[i + 1];
    for (int s = 0; s < N; s++) for (int c = 0; c < 2; c++) {
      const char *t = m + s * 7 + c * 3;
      if (t[0] == '-') NS[s][c] = -1; else { W[s][c] = t[0] - '0'; D[s][c] = t[1] == 'R'; NS[s][c] = t[2] - 'A'; }
    }
    int hs, hc, ms; long long ones;
    int t = simulate(&hs, &hc, &ones, &ms);
    printf("MACHINE %s: ", m);
    if (t >= 0) { printf("halts after %d steps with %lld ones\n", t + 1, ones); return 0; }
    for (int r = 1; r <= MAXRAD; r++) if (ngram_cps(r, 0)) { printf("proved non-halting (radius %d)\n", r); return 0; }
    for (int h = 1; h <= HDEPTH; h++) for (int r = 1; r <= HISTRAD; r++)
      if (ngram_cps(r, h)) { printf("proved non-halting (history %d, radius %d)\n", h, r); return 0; }
    printf("not decided within %d steps\n", LIMIT);
    return 0;
  }
  explore(0);
  printf("n = %d, step limit %d, radius <= %d: halting machines %lld, proved non-halting %lld (by radius:", N, LIMIT,
         MAXRAD, n_halt, n_nonhalt);
  for (int r = 1; r <= MAXRAD; r++) printf(" %lld", by_rad[r]);
  printf("; with the last h writers, h = 1..%d:", HDEPTH);
  for (int h = 1; h <= HDEPTH; h++) printf(" %lld", by_depth[h]);
  printf("), holdouts %lld\n", n_holdout);
  if (SELFTEST) printf("SELFTEST: every test applied to each of the %lld halting machines: %lld false proofs\n", n_halt, n_selftest_fail);
  if (n_holdout == 0) printf("RESULT n = %d: S = %lld, Sigma = %lld (every machine decided)\n", N, bestS, bestSigma);
  else printf("RESULT n = %d: undecided (%lld holdouts); the halting machines give S >= %lld, Sigma >= %lld\n", N,
              n_holdout, bestS, bestSigma);
  return 0;
}
