/* 14a, program B (written separately from program A): BB(n) for n = 2, 3, 4 by enumeration with proofs.
 *
 * Tree normal form, built with an explicit stack of saved runs: each node keeps its tape, head, state and step
 * count and resumes the run (nothing is re-simulated). The first transition moves right; a new transition may go
 * to any state already used or to the next unused one. A node whose run reads an undefined transition gives one
 * halting machine (that transition halts, writing 1) and, while two or more transitions are undefined, children
 * for every definition of it. A run that passes the step limit must be shown never to halt, by the first of:
 *   1. cycler: the full configuration recurs (Brent's method: compared with the one saved at each power of 2);
 *   2. translated cycler: two records of the head on the same side, in the same state, with the cells from the
 *      record back to the deepest point the head reached in between equal at both times (and blank beyond);
 *   3. closed sets of local views: windows of r cells on each side of the head, each cell carrying its last h
 *      writer states (h = 0..3), with the sets of windows seen on each side closed under every step; no view may
 *      reach an undefined transition.
 * Output: S(n), Sigma(n), and how many machines each method settled.
 * Usage: c14a_b N [LIMIT] [--no-translated] [--hmax H] [--rmax R]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MS 5
static int n, LIM = 5000, HMAX = 9, RMAX = 12, TRANSLATED = 1, SELFTEST = 0;
static long long selftest_bad = 0;
static int tr[MS][2];                     /* -1 undefined; else write | move << 1 | next << 2 (move 1 = right) */
static long long cnt_halt, cnt_cyc, cnt_tc, cnt_cps[10][16], cnt_left;
static long long S_best, Sigma_best;
static int TL;

#define WR(t) ((t) & 1)
#define MV(t) (((t) >> 1) & 1)
#define NX(t) ((t) >> 2)

typedef struct { unsigned char *tape; int head, st, steps, maxst, ndef; } node;

/* ---------------------------------------------------------------- 1. cycler */
static unsigned char *buf1, *buf2;
static int is_cycler(void) {
  memset(buf1, 0, TL); int h = TL / 2, s = 0, lo = h, hi = h;
  int sh = h, ss = s, slo = lo, shi = hi; long long mark = 1;
  memcpy(buf2, buf1, TL);
  for (long long t = 1; t <= LIM; t++) {
    int x = tr[s][buf1[h]]; if (x < 0) return 0;
    buf1[h] = WR(x); h += MV(x) ? 1 : -1; s = NX(x);
    if (h < lo) lo = h; if (h > hi) hi = h;
    if (s == ss && h == sh) {
      int a = lo < slo ? lo : slo, b = hi > shi ? hi : shi;
      if (!memcmp(buf1 + a, buf2 + a, b - a + 1)) return 1;
    }
    if (t == mark) { memcpy(buf2, buf1, TL); sh = h; ss = s; slo = lo; shi = hi; mark *= 2; }
  }
  return 0;
}

/* ---------------------------------------------------------------- 2. translated cycler */
#define WIN 200
#define MAXREC 4000
typedef struct { int t, pos, st; unsigned char w[WIN + 1]; } rec;
static rec *recs[2]; static int nrec[2]; static int *poshist;
static int is_translated(void) {
  memset(buf1, 0, TL); int h = TL / 2, s = 0, lo = h, hi = h;
  nrec[0] = nrec[1] = 0; poshist[0] = h;
  for (int t = 1; t <= LIM; t++) {
    int x = tr[s][buf1[h]]; if (x < 0) return 0;
    buf1[h] = WR(x); h += MV(x) ? 1 : -1; s = NX(x); poshist[t] = h;
    int side = -1;
    if (h > hi) { hi = h; side = 1; }
    if (h < lo) { lo = h; side = 0; }
    if (side < 0) continue;
    /* compare with earlier records on this side in the same state */
    for (int i = nrec[side] - 1; i >= 0; i--) {
      rec *r = &recs[side][i];
      if (r->st != s) continue;
      int D = 0;                                  /* deepest excursion back from r->pos during [r->t, t] */
      for (int u = r->t; u <= t; u++) {
        int dd = side ? r->pos - poshist[u] : poshist[u] - r->pos;
        if (dd > D) D = dd;
      }
      if (D > WIN) continue;
      int same = 1;
      for (int k = 0; k <= D && same; k++) {
        unsigned char now = side ? buf1[h - k] : buf1[h + k];
        if (now != r->w[k]) same = 0;
      }
      if (same) return 1;
    }
    if (nrec[side] < MAXREC) {
      rec *r = &recs[side][nrec[side]++];
      r->t = t; r->pos = h; r->st = s;
      for (int k = 0; k <= WIN; k++) {
        int p = side ? h - k : h + k;
        r->w[k] = (p >= 0 && p < TL) ? buf1[p] : 0;
      }
    }
  }
  return 0;
}

/* ---------------------------------------------------------------- 3. closed sets of local views */
typedef unsigned __int128 key_t128;
typedef struct { key_t128 *k; uint32_t *g; uint32_t gen; size_t cap, used; } set_t;
static set_t SL, SR, SV;
static void set_new(set_t *s, size_t cap) { s->k = malloc(cap * sizeof(key_t128)); s->g = calloc(cap, 4); s->gen = 1; s->cap = cap; s->used = 0; }
static void set_reset(set_t *s) { s->gen++; s->used = 0; }
static size_t slot(set_t *s, key_t128 key) {
  uint64_t z = (uint64_t)key + 0x632BE59BD9B4E019ULL + (uint64_t)(key >> 64) * 0xD6E8FEB86659FD93ULL;
  z ^= z >> 31; z *= 0x9E3779B97F4A7C15ULL; z ^= z >> 29;
  size_t i = z & (s->cap - 1);
  while (s->g[i] == s->gen && s->k[i] != key) i = (i + 1) & (s->cap - 1);
  return i;
}
static int set_in(set_t *s, key_t128 key) { return s->g[slot(s, key)] == s->gen; }
static int set_put(set_t *s, key_t128 key) {          /* 1 new, 0 present, -1 full */
  size_t i = slot(s, key);
  if (s->g[i] == s->gen) return 0;
  if (4 * (s->used + 1) > 3 * s->cap) return -1;
  s->g[i] = s->gen; s->k[i] = key; s->used++; return 1;
}
static key_t128 *views; static size_t nviews, capviews = 1u << 22;

static int closed_views(int r, int hdepth) {
  uint64_t base = 1;                                   /* tags: last hdepth writers, each 0 (none) or state+1 */
  for (int i = 0; i < hdepth; i++) base *= (uint64_t)(n + 1);
  uint64_t A = 2 * base;                               /* alphabet: bit * base + tag */
  uint64_t Ar = 1;
  for (int i = 0; i < r; i++) { if (Ar > (1ULL << 56) / A) return 0; Ar *= A; }   /* A^r <= 2^56, checked before multiplying */
  if ((key_t128)Ar * Ar > (((key_t128)1) << 126) / ((key_t128)n * A)) return 0;   /* view keys stay below 2^126 */
  set_reset(&SL); set_reset(&SR); set_reset(&SV); nviews = 0;
  set_put(&SL, 0); set_put(&SR, 0);
  /* view key = ((left * n + state) * A + head) * Ar + right; window digit 0 is the cell next to the head */
  set_put(&SV, 0); views[nviews++] = 0;
  for (int again = 1; again; ) {
    again = 0;
    for (size_t i = 0; i < nviews; i++) {
      key_t128 v = views[i];
      uint64_t right = (uint64_t)(v % Ar); v /= Ar; uint64_t head = (uint64_t)(v % A); v /= A; int st = (int)(v % n);
      uint64_t left = (uint64_t)(v / n);
      int bit = (int)(head / base);
      int x = tr[st][bit];
      if (x < 0) return 0;
      uint64_t tag = hdepth ? ((head % base) * (uint64_t)(n + 1) + (uint64_t)(st + 1)) % base : 0;
      uint64_t written = (uint64_t)WR(x) * base + tag;
      int nst = NX(x);
      uint64_t *grow_w = MV(x) ? &left : &right, *shrink_w = MV(x) ? &right : &left;
      set_t *grow_s = MV(x) ? &SL : &SR, *shrink_s = MV(x) ? &SR : &SL;
      uint64_t g2 = (uint64_t)(((key_t128)*grow_w * A + written) % Ar);
      int p = set_put(grow_s, g2); if (p < 0) return 0; if (p) again = 1;
      uint64_t nhead = *shrink_w % A;
      for (uint64_t far = 0; far < A; far++) {
        uint64_t s2 = *shrink_w / A + far * (Ar / A);
        if (!set_in(shrink_s, s2)) continue;
        uint64_t L = MV(x) ? g2 : s2, R = MV(x) ? s2 : g2;
        key_t128 key = (((key_t128)L * n + nst) * A + nhead) * Ar + R;
        int q = set_put(&SV, key); if (q < 0) return 0;
        if (q) { if (nviews == capviews) return 0; views[nviews++] = key; again = 1; }
      }
    }
  }
  return 1;
}

/* ---------------------------------------------------------------- enumeration */
static void show(void) {
  for (int s = 0; s < n; s++) {
    if (s) putchar('_');
    for (int c = 0; c < 2; c++)
      if (tr[s][c] < 0) printf("---"); else printf("%d%c%c", WR(tr[s][c]), MV(tr[s][c]) ? 'R' : 'L', 'A' + NX(tr[s][c]));
  }
}

static void decide(void) {
  if (is_cycler()) { cnt_cyc++; return; }
  if (TRANSLATED && is_translated()) { cnt_tc++; return; }
  for (int hd = 0; hd <= HMAX; hd++)
    for (int r = 1; r <= (hd ? 4 : RMAX); r++)
      if (closed_views(r, hd)) { cnt_cps[hd][r]++; return; }
  cnt_left++; printf("UNDECIDED "); show(); putchar('\n');
}

static void grow(node *nd) {
  unsigned char *tp = nd->tape;
  while (nd->steps < LIM) {
    int c = tp[nd->head], x = tr[nd->st][c];
    if (x < 0) {
      long long ones = 0;
      for (int i = 0; i < TL; i++) ones += tp[i];
      ones += (c == 0);
      cnt_halt++;
      if (SELFTEST) {                               /* a halting machine must defeat every test */
        int sv = tr[nd->st][c];
        if (is_cycler()) selftest_bad++;
        if (TRANSLATED && is_translated()) selftest_bad++;
        for (int hd = 0; hd <= HMAX; hd++) for (int r = 1; r <= (hd ? 4 : RMAX); r++) if (closed_views(r, hd)) selftest_bad++;
        (void)sv;
      }
      if (nd->steps + 1 > S_best) S_best = nd->steps + 1;
      if (ones > Sigma_best) Sigma_best = ones;
      if (nd->ndef + 1 >= 2 * n) return;               /* keep one transition undefined */
      int hs = nd->st, top = nd->maxst + 1 < n - 1 ? nd->maxst + 1 : n - 1;
      for (int nx = 0; nx <= top; nx++)
        for (int mv = 1; mv >= 0; mv--) {
          if (nd->ndef == 0 && mv == 0) continue;
          for (int w = 0; w <= 1; w++) {
            node ch = *nd;
            ch.tape = malloc(TL); memcpy(ch.tape, tp, TL);
            tr[hs][c] = w | mv << 1 | nx << 2;
            ch.ndef = nd->ndef + 1;
            grow(&ch);
            tr[hs][c] = -1;
            free(ch.tape);
          }
        }
      return;
    }
    tp[nd->head] = WR(x); nd->head += MV(x) ? 1 : -1; nd->st = NX(x); nd->steps++;
    if (nd->st > nd->maxst) nd->maxst = nd->st;
    if (nd->head < 1 || nd->head > TL - 2) { fprintf(stderr, "tape too short\n"); exit(2); }
  }
  decide();
}

int main(int argc, char **argv) {
  if (argc < 2) return 1;
  n = atoi(argv[1]);
  for (int i = 2; i < argc; i++) {
    if (!strcmp(argv[i], "--no-translated")) TRANSLATED = 0;
    else if (!strcmp(argv[i], "--selftest")) SELFTEST = 1;
    else if (!strcmp(argv[i], "--hmax")) HMAX = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--rmax")) RMAX = atoi(argv[++i]);
    else LIM = atoi(argv[i]);
  }
  if (n < 1 || n > MS || HMAX < 0 || HMAX > 9 || RMAX < 1 || RMAX > 15) return 1;
  TL = 2 * LIM + 8;
  buf1 = malloc(TL); buf2 = malloc(TL); poshist = malloc((LIM + 1) * sizeof(int));
  recs[0] = malloc(MAXREC * sizeof(rec)); recs[1] = malloc(MAXREC * sizeof(rec));
  set_new(&SL, 1u << 20); set_new(&SR, 1u << 20); set_new(&SV, 1u << 23); views = malloc(capviews * sizeof(key_t128));
  for (int s = 0; s < MS; s++) tr[s][0] = tr[s][1] = -1;
  node root; root.tape = calloc(TL, 1); root.head = TL / 2; root.st = 0; root.steps = 0; root.maxst = 0; root.ndef = 0;
  grow(&root);
  long long ncps = 0, byh[10] = {0};
  for (int a = 0; a < 10; a++) for (int c = 0; c < 16; c++) { ncps += cnt_cps[a][c]; byh[a] += cnt_cps[a][c]; }
  printf("n = %d, limit %d: halting %lld; never halting: cyclers %lld, translated cyclers %lld, closed view sets %lld"
         " (by writer-history depth 0..%d:", n, LIM, cnt_halt, cnt_cyc, cnt_tc, ncps, HMAX);
  for (int a = 0; a <= HMAX; a++) printf(" %lld", byh[a]);
  printf("); undecided %lld\n", cnt_left);
  if (SELFTEST) printf("SELFTEST: every test applied to each of the %lld halting machines: %lld false proofs\n", cnt_halt, selftest_bad);
  if (cnt_left == 0) printf("RESULT n = %d: S = %lld, Sigma = %lld (every machine decided)\n", n, S_best, Sigma_best);
  else printf("RESULT n = %d: undecided (%lld machines); halting machines give S >= %lld, Sigma >= %lld\n", n, cnt_left, S_best, Sigma_best);
  return 0;
}
