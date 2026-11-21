/* Program B (C, exact 128-bit integer orientation) for the upper-bound side of Rohan Arun's Erdos problem
 * 838 minima (github.com/rohanarun/erdos-838-convex-polygon-minima @ c301008).  Written separately from
 * program A (no shared code), without reading or running the claimant's programs.
 *
 * Method B: for every subset S of size >= 3 the convex hull is built by Andrew's monotone chain with exact
 * integer orientation (points sorted by x then y; a point is popped when the turn is not strictly left);
 * S is in convex position iff its hull has |S| vertices.  No three collinear: every triple's orientation
 * is nonzero.  Counts by size and the totals after deleting each single point are reported.
 *
 * Usage: erdos838_b [--control collinear|move] SPEC ...
 *   SPEC = PATH:FORMAT:INDEX=TOTAL[,b3,b4,...][@d1/d2/...] with FORMAT xy | jsonl | lean
 *   (jsonl: the INDEX-th line of a JSON Lines file whose record ends in a list of [x, y] pairs;
 *    lean: lines '| i => (x, y)').
 * Prints 'B: VERDICT PASS' (exit 0) or 'B: VERDICT REFUSED' (exit 1).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef __int128 i128;
typedef struct { long long x, y; } P;
static int fails = 0;
#define FAIL(...) do { if (fails < 20) { printf("B: FAIL "); printf(__VA_ARGS__); printf("\n"); } fails++; } while (0)

static i128 cross(P o, P a, P b) { return (i128)(a.x - o.x) * (b.y - o.y) - (i128)(a.y - o.y) * (b.x - o.x); }
static int cmpP(const void *u, const void *v) {
  const P *a = u, *b = v;
  if (a->x != b->x) return a->x < b->x ? -1 : 1;
  if (a->y != b->y) return a->y < b->y ? -1 : 1;
  return 0;
}
static int hull_size(P *q, int k) {   /* q sorted; returns number of strict hull vertices */
  P h[64]; int m = 0;
  for (int i = 0; i < k; i++) { while (m >= 2 && cross(h[m - 2], h[m - 1], q[i]) <= 0) m--; h[m++] = q[i]; }
  for (int i = k - 2, t = m + 1; i >= 0; i--) { while (m >= t && cross(h[m - 2], h[m - 1], q[i]) <= 0) m--; h[m++] = q[i]; }
  return m - 1;
}

static int read_ints_line(const char *s, long long *v, int max) {   /* all integers on the line after any quoted string */
  int c = 0; const char *p = s;
  const char *q = strrchr(s, '"'); if (q) p = q + 1;
  while (*p && c < max) {
    if ((*p == '-' && p[1] >= '0' && p[1] <= '9') || (*p >= '0' && *p <= '9')) { char *e; v[c++] = strtoll(p, &e, 10); p = e; }
    else p++;
  }
  return c;
}

static int load(const char *path, const char *fmt, int index, P *pts) {
  FILE *f = fopen(path, "r"); if (!f) { perror(path); exit(2); }
  static char line[1 << 16]; int n = 0, li = 0;
  while (fgets(line, sizeof line, f)) {
    long long v[64];
    if (!strcmp(fmt, "xy")) { int c = read_ints_line(line, v, 64); if (c >= 2) { pts[n].x = v[0]; pts[n].y = v[1]; n++; } }
    else if (!strcmp(fmt, "jsonl")) {
      if (strlen(line) < 3) continue;
      if (li++ != index) continue;
      int c = read_ints_line(line, v, 64);
      for (int i = 0; i + 1 < c; i += 2) { pts[n].x = v[i]; pts[n].y = v[i + 1]; n++; }
    } else if (!strcmp(fmt, "lean")) {
      char *a = strstr(line, "=>"); char *bar = strchr(line, '|');
      if (a && bar && strchr(a, '(')) { int c = read_ints_line(bar, v, 64); if (c == 3) { if (v[0] != n) { printf("B: lean index out of order\n"); exit(2); } pts[n].x = v[1]; pts[n].y = v[2]; n++; } }
    }
  }
  fclose(f);
  return n;
}

int main(int argc, char **argv) {
  const char *control = "";
  int ai = 1;
  if (argc > 2 && !strcmp(argv[1], "--control")) { control = argv[2]; ai = 3; }
  for (int w = ai; w < argc; w++) {
    char spec[4096]; strncpy(spec, argv[w], sizeof spec - 1); spec[sizeof spec - 1] = 0;
    char *eq = strrchr(spec, '='); *eq = 0; char *rhs = eq + 1;
    char *path = spec, *fmt = "xy"; int index = 0;
    char *c1 = strchr(spec, ':'); if (c1) { *c1 = 0; fmt = c1 + 1; char *c2 = strchr(fmt, ':'); if (c2) { *c2 = 0; index = atoi(c2 + 1); } }
    long long claim[32]; int nclaim = 0; long long dels[32]; int ndel = 0;
    char *at = strchr(rhs, '@'); if (at) { *at = 0; for (char *t = strtok(at + 1, "/"); t; t = strtok(0, "/")) dels[ndel++] = atoll(t); }
    for (char *t = strtok(rhs, ","); t; t = strtok(0, ",")) claim[nclaim++] = atoll(t);
    P pts[64]; int n = load(path, fmt, index, pts);
    const char *base = strrchr(path, '/'); base = base ? base + 1 : path;
    char tag[512]; if (!strcmp(fmt, "jsonl")) snprintf(tag, sizeof tag, "%s[%d]", base, index); else snprintf(tag, sizeof tag, "%s", base);
    if (w == ai && !strcmp(control, "collinear")) { pts[2].x = 2 * pts[1].x - pts[0].x; pts[2].y = 2 * pts[1].y - pts[0].y; printf("B: CONTROL collinear: %s point 2 replaced by 2*P1 - P0 = (%lld, %lld)\n", tag, pts[2].x, pts[2].y); }
    if (w == ai && !strcmp(control, "move")) { pts[0].x -= 1000; pts[0].y += 37; printf("B: CONTROL move: %s point 0 moved to (%lld, %lld)\n", tag, pts[0].x, pts[0].y); }
    int ncol = 0, nrep = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
      if (pts[i].x == pts[j].x && pts[i].y == pts[j].y) nrep++;
      for (int k = j + 1; k < n; k++) if (cross(pts[i], pts[j], pts[k]) == 0) ncol++;
    }
    if (nrep) FAIL("%s: %d repeated points", tag, nrep);
    if (ncol) FAIL("%s: %d collinear triples", tag, ncol);
    long long bysize[65] = {0}, dele[64] = {0};
    for (uint32_t S = 0; S < (1u << n); S++) {
      int k = __builtin_popcount(S); if (k < 3) continue;
      P q[64]; int m = 0;
      for (int i = 0; i < n; i++) if (S >> i & 1) q[m++] = pts[i];
      qsort(q, m, sizeof(P), cmpP);
      if (hull_size(q, m) == k) { bysize[k]++; for (int i = 0; i < n; i++) if (!(S >> i & 1)) dele[i]++; }
    }
    long long tot = 0; for (int k = 3; k <= n; k++) tot += bysize[k];
    /* sort deletion totals */
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) if (dele[j] < dele[i]) { long long t = dele[i]; dele[i] = dele[j]; dele[j] = t; }
    printf("B: %s: n=%d, collinear triples %d, convex subsets by size 3..%d =", tag, n, ncol, n);
    for (int k = 3; k <= n; k++) printf("%s%lld", k == 3 ? " " : ",", bysize[k]);
    printf(", total %lld (claimed %lld)", tot, claim[0]);
    if (ndel) { printf(", deletion totals "); for (int i = 0; i < n; i++) printf("%s%lld", i ? "/" : "", dele[i]); }
    printf("\n");
    if (tot != claim[0]) FAIL("%s: total %lld != claimed %lld", tag, tot, claim[0]);
    if (nclaim > 1) {
      for (int k = 3; k <= n; k++) { long long want = (k - 3 + 1 < nclaim) ? claim[k - 3 + 1] : 0; if (bysize[k] != want) FAIL("%s: size %d count %lld != claimed %lld", tag, k, bysize[k], want); }
    }
    if (ndel) { if (ndel != n) FAIL("%s: %d deletion totals claimed for %d points", tag, ndel, n); else for (int i = 0; i < n; i++) if (dele[i] != dels[i]) FAIL("%s: deletion totals differ at rank %d", tag, i); }
  }
  if (fails) { printf("B: VERDICT REFUSED (%d failures)\n", fails); return 1; }
  printf("B: VERDICT PASS\n");
  return 0;
}
