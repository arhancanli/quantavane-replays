/* Homomorphism counts of K_4 into a graph and into its looped complement, for the 89a row [PPSS2024].
 * Input: one graph in graph6 format (n <= 2048). Output (one line each):
 *   n, edges
 *   routeA red  : sum over roots a of tr((X D_a)^3), X = adjacency (no loops)
 *   routeA blue : the same for X = looped complement (X_vv = 1)
 *   routeB red  : number of 4-cliques c4 (u < v < w < x) of the graph
 *   routeB blue : c1, c2, c3, c4 of the complement (as a simple graph)
 *   roots       : the distinct per-root values of tr((X D_a)^3) for each colour (vertex-transitivity check)
 * Route A counts ordered 4-tuples (a, b, c, d), repeats allowed, with all six pairs coloured X; route B counts
 * cliques, and hom(K_4 looped, X looped) = c1 + 14 c2 + 36 c3 + 24 c4 (surjections of 4 points onto k points),
 * hom(K_4, X) = 24 c4 without loops. The two routes share no counting code.
 * Usage: k4hom FILE [flip u v]   (flip: toggle the edge uv first, for the red control) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXN 2048
#define W (MAXN / 64)
static uint64_t R[MAXN][W], B[MAXN][W], Bl[MAXN][W];
static int n, nw;

static int pc(uint64_t x) { return __builtin_popcountll(x); }
static void set(uint64_t *m, int i) { m[i >> 6] |= 1ULL << (i & 63); }
static int get(const uint64_t *m, int i) { return (m[i >> 6] >> (i & 63)) & 1; }

static long long trace3(uint64_t (*X)[W], int a) {
  /* sum over b, c, d in S = X[a] of X_bc X_cd X_db */
  const uint64_t *S = X[a];
  long long t = 0;
  uint64_t Sb[W], Sbc[W];
  for (int b = 0; b < n; b++) {
    if (!get(S, b)) continue;
    for (int k = 0; k < nw; k++) Sb[k] = S[k] & X[b][k];
    for (int c = 0; c < n; c++) {
      if (!get(Sb, c)) continue;
      int s = 0;
      for (int k = 0; k < nw; k++) { Sbc[k] = Sb[k] & X[c][k]; s += pc(Sbc[k]); }
      t += s;
    }
  }
  return t;
}

static void cliques(uint64_t (*X)[W], long long *c) {
  /* c[1..4]: numbers of cliques of sizes 1..4 of the simple graph X */
  uint64_t U[W], UV[W], UVW[W];
  c[1] = n; c[2] = c[3] = c[4] = 0;
  for (int u = 0; u < n; u++) {
    memset(U, 0, sizeof U);
    for (int k = 0; k < nw; k++) U[k] = X[u][k];
    for (int i = 0; i <= u; i++) U[i >> 6] &= ~(1ULL << (i & 63));         /* keep v > u */
    for (int v = u + 1; v < n; v++) {
      if (!get(U, v)) continue;
      c[2]++;
      for (int k = 0; k < nw; k++) UV[k] = U[k] & X[v][k];
      for (int i = 0; i <= v; i++) UV[i >> 6] &= ~(1ULL << (i & 63));
      for (int w = v + 1; w < n; w++) {
        if (!get(UV, w)) continue;
        c[3]++;
        long long s = 0;
        for (int k = 0; k < nw; k++) UVW[k] = UV[k] & X[w][k];
        for (int i = 0; i <= w; i++) UVW[i >> 6] &= ~(1ULL << (i & 63));
        for (int k = 0; k < nw; k++) s += pc(UVW[k]);
        c[4] += s;
      }
    }
  }
}

int main(int argc, char **argv) {
  if (argc < 2) { fprintf(stderr, "usage: k4hom FILE [flip u v]\n"); return 2; }
  FILE *f = fopen(argv[1], "rb");
  if (!f) { perror("open"); return 2; }
  static unsigned char buf[1 << 22];
  size_t len = fread(buf, 1, sizeof buf, f);
  fclose(f);
  while (len && (buf[len - 1] == '\n' || buf[len - 1] == '\r')) len--;
  size_t pos;
  if (buf[0] < 126) { n = buf[0] - 63; pos = 1; }
  else if (buf[1] < 126) { n = ((buf[1] - 63) << 12) | ((buf[2] - 63) << 6) | (buf[3] - 63); pos = 4; }
  else { fprintf(stderr, "n too large\n"); return 2; }
  if (n > MAXN || n < 1) { fprintf(stderr, "bad n\n"); return 2; }
  nw = (n + 63) / 64;
  long long need_bits = (long long)n * (n - 1) / 2;
  if ((long long)(len - pos) != (need_bits + 5) / 6) { fprintf(stderr, "length mismatch\n"); return 2; }
  long long bit = 0, edges = 0;
  for (int j = 1; j < n; j++)
    for (int i = 0; i < j; i++, bit++) {
      int byte = buf[pos + bit / 6];
      if (byte < 63 || byte > 126) { fprintf(stderr, "bad byte\n"); return 2; }
      if (((byte - 63) >> (5 - bit % 6)) & 1) { set(R[i], j); set(R[j], i); edges++; }
    }
  /* trailing padding bits must be zero */
  for (long long b2 = bit; b2 % 6; b2++) if (((buf[pos + b2 / 6] - 63) >> (5 - b2 % 6)) & 1) { fprintf(stderr, "padding\n"); return 2; }
  if (argc == 5 && !strcmp(argv[2], "flip")) {
    int u = atoi(argv[3]), v = atoi(argv[4]);
    R[u][v >> 6] ^= 1ULL << (v & 63); R[v][u >> 6] ^= 1ULL << (u & 63);
    edges += get(R[u], v) ? 1 : -1;
  }
  for (int v = 0; v < n; v++)
    for (int w = 0; w < n; w++)
      if (w != v && !get(R[v], w)) { set(B[v], w); set(Bl[v], w); }
  for (int v = 0; v < n; v++) set(Bl[v], v);          /* loops in the complement */
  long long Nr = 0, Nb = 0, rmin = -1, rmax = -1, bmin = -1, bmax = -1;
  for (int a = 0; a < n; a++) {
    long long tr = trace3(R, a), tb = trace3(Bl, a);
    Nr += tr; Nb += tb;
    if (rmin < 0 || tr < rmin) rmin = tr;
    if (tr > rmax) rmax = tr;
    if (bmin < 0 || tb < bmin) bmin = tb;
    if (tb > bmax) bmax = tb;
  }
  long long cr[5], cb[5];
  cliques(R, cr);
  cliques(B, cb);
  printf("n %d edges %lld\n", n, edges);
  printf("routeA red %lld blue %lld\n", Nr, Nb);
  printf("routeB red %lld %lld %lld %lld blue %lld %lld %lld %lld\n", cr[1], cr[2], cr[3], cr[4], cb[1], cb[2], cb[3], cb[4]);
  printf("roots red %lld %lld blue %lld %lld\n", rmin, rmax, bmin, bmax);
  return 0;
}
