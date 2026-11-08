/* flag_count.c -- our own counting code for the 89c lower-bound replay (PR #308).
 *
 * Mode "sums": for every 8-vertex graph H in GRAPHS (one 28-bit edge mask per line; edge (x,y), x < y, is bit
 * y(y-1)/2 + x), print "idx mono S6 S4 S2 S0", where mono is the number of 6-sets that are cliques or independent
 * sets of H, and S_t is the sum, over the blocks with t roots (both colours; colour 1 reads the complement of H), over
 * all ordered root tuples R of size t whose labelled induced graph equals the block's root graph, and over all
 * ordered pairs (A, B) of disjoint m-sets outside R, of G_b[i][j] = sum_r L_ri L_rj, where i and j are the listed
 * flags equal to the induced graphs on R+A and R+B up to permutations of the extension vertices (pairs with an
 * unlisted flag add 0). Exact: 128-bit integers, overflow checked.
 * Mode "aut": for every graph print "idx canon aut": the least mask over the 8! relabellings, and the number of
 * relabellings that fix the mask.
 * The certificate is read in the format stated in PROOF.md of the packet (header n k D blocks bnum bden mode; per
 * block t m root nflags nrows rootmode colour, the flag codes, then the integer rows).
 * Usage: flag_count sums CERT GRAPHS | flag_count aut GRAPHS
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef __int128 i128;
#define PIDX(x, y) ((y) * ((y) - 1) / 2 + (x))

static void die(const char *m) { fprintf(stderr, "flag_count: %s\n", m); exit(2); }

static void print128(i128 v) {
  char buf[64]; int p = 63; int neg = v < 0; unsigned __int128 u = neg ? -(unsigned __int128)v : (unsigned __int128)v;
  buf[p] = 0;
  if (u == 0) buf[--p] = '0';
  while (u) { buf[--p] = '0' + (int)(u % 10); u /= 10; }
  if (neg) buf[--p] = '-';
  fputs(buf + p, stdout);
}

typedef struct { int t, m, colour, root, nf, nr, extbits; long long *codes; long long *rows; i128 *G; int *table; } Block;

static int nblocks;
static Block *B;
static int lookup[2][7][1 << 15];         /* [colour][t][rootcode] -> block or -1 */

static int pair_count(int v) { return v * (v - 1) / 2; }

static long long permute_ext(long long code, int t, int m, const int *sg) {
  /* relabel extension vertices t+k -> t+sg[k], roots fixed */
  int v = t + m; long long out = 0;
  for (int y = 1; y < v; y++) for (int x = 0; x < y; x++) if (code >> PIDX(x, y) & 1) {
    int x2 = x < t ? x : t + sg[x - t], y2 = y < t ? y : t + sg[y - t];
    if (x2 > y2) { int tmp = x2; x2 = y2; y2 = tmp; }
    out |= 1LL << PIDX(x2, y2);
  }
  return out;
}

static int next_perm(int *a, int n) {
  int i = n - 2; while (i >= 0 && a[i] >= a[i + 1]) i--; if (i < 0) return 0;
  int j = n - 1; while (a[j] <= a[i]) j--; int tmp = a[i]; a[i] = a[j]; a[j] = tmp;
  for (int l = i + 1, r = n - 1; l < r; l++, r--) { tmp = a[l]; a[l] = a[r]; a[r] = tmp; }
  return 1;
}

static void read_cert(const char *path) {
  FILE *f = fopen(path, "r"); if (!f) die("cannot open certificate");
  long long n, k, D, nb, bn, bd, mode;
  if (fscanf(f, "%lld %lld %lld %lld %lld %lld %lld", &n, &k, &D, &nb, &bn, &bd, &mode) != 7) die("header");
  if (n != 8 || k != 6 || mode != -2) die("unexpected header");
  nblocks = (int)nb; B = calloc(nblocks, sizeof(Block));
  memset(lookup, -1, sizeof lookup);
  for (int b = 0; b < nblocks; b++) {
    long long t, m, root, nf, nr, rmode, col;
    if (fscanf(f, "%lld %lld %lld %lld %lld %lld %lld", &t, &m, &root, &nf, &nr, &rmode, &col) != 7) die("block header");
    if (t + 2 * m != 8 || rmode != 1 || (col != 0 && col != 1)) die("unexpected block type");
    Block *X = &B[b]; X->t = (int)t; X->m = (int)m; X->colour = (int)col; X->root = (int)root;
    X->nf = (int)nf; X->nr = (int)nr; X->extbits = pair_count(t + m) - pair_count(t);
    X->codes = malloc(nf * sizeof(long long)); X->rows = malloc((nr ? nr : 1) * nf * sizeof(long long));
    for (int i = 0; i < nf; i++) if (fscanf(f, "%lld", &X->codes[i]) != 1) die("flags");
    for (long long r = 0; r < nr * nf; r++) if (fscanf(f, "%lld", &X->rows[r]) != 1) die("rows");
    if (lookup[col][t][root] != -1) die("two blocks with the same colour and root graph");
    lookup[col][t][root] = b;
    /* Gram matrix */
    X->G = calloc((size_t)nf * nf, sizeof(i128));
    for (int r = 0; r < nr; r++) {
      long long *L = X->rows + (size_t)r * nf;
      for (int i = 0; i < nf; i++) { if (!L[i]) continue;
        for (int j = 0; j < nf; j++) {
          i128 p = (i128)L[i] * (i128)L[j];
          X->G[(size_t)i * nf + j] += p;
        } }
    }
    /* flag table over extension codes, all extension permutations */
    X->table = malloc(sizeof(int) << X->extbits); for (int e = 0; e < (1 << X->extbits); e++) X->table[e] = -1;
    int rmask = (1 << pair_count(t)) - 1;
    for (int i = 0; i < nf; i++) {
      if ((X->codes[i] & rmask) != root) die("a flag whose root part is not the block's root graph");
      if (X->codes[i] >> pair_count(t + m)) die("a flag code out of range");
      int sg[4] = {0, 1, 2, 3};
      do {
        long long pc = permute_ext(X->codes[i], (int)t, (int)m, sg);
        int e = (int)(pc >> pair_count(t));
        if (X->table[e] != -1 && X->table[e] != i) die("two listed flags in one orbit");
        X->table[e] = i;
      } while (next_perm(sg, (int)m));
    }
  }
  int extra; if (fscanf(f, "%d", &extra) == 1) die("trailing data in certificate");
  fclose(f);
}

static int ngraphs; static long long *gm;
static void read_graphs(const char *path) {
  FILE *f = fopen(path, "r"); if (!f) die("cannot open graphs");
  int cap = 16384; gm = malloc(cap * sizeof(long long)); long long v;
  while (fscanf(f, "%lld", &v) == 1) { if (ngraphs == cap) { cap *= 2; gm = realloc(gm, cap * sizeof(long long)); } gm[ngraphs++] = v; }
  fclose(f);
}

/* ordered tuples and extension pairs, independent of H */
typedef struct { int r[6]; int rootpair[15]; int npairs; int A[70][4], Bv[70][4]; } Tuple;

static int adjm(long long mask, int x, int y) { if (x > y) { int t = x; x = y; y = t; } return (int)(mask >> PIDX(x, y) & 1); }

static int mode_sums(const char *cert, const char *graphs) {
  read_cert(cert); read_graphs(graphs);
  static const int TS[4] = {6, 4, 2, 0}, MS[4] = {1, 2, 3, 4};
  Tuple *tup[4]; int ntup[4];
  for (int q = 0; q < 4; q++) {
    int t = TS[q], m = MS[q], cnt = 0; tup[q] = malloc(20160 * sizeof(Tuple));
    int idx[8];
    /* all ordered t-tuples: enumerate 8^t with distinctness (once) */
    long long total = 1; for (int i = 0; i < t; i++) total *= 8;
    for (long long code = 0; code < total; code++) {
      long long c = code; int used = 0, ok = 1;
      for (int i = 0; i < t; i++) { idx[i] = (int)(c % 8); c /= 8; if (used >> idx[i] & 1) { ok = 0; break; } used |= 1 << idx[i]; }
      if (!ok) continue;
      Tuple *T = &tup[q][cnt++];
      for (int i = 0; i < t; i++) T->r[i] = idx[i];
      for (int y = 1; y < t; y++) for (int x = 0; x < y; x++) {
        int a = idx[x], b = idx[y]; if (a > b) { int tmp = a; a = b; b = tmp; }
        T->rootpair[PIDX(x, y)] = PIDX(a, b);
      }
      int rem[8], nr = 0; for (int v = 0; v < 8; v++) if (!(used >> v & 1)) rem[nr++] = v;
      /* ordered pairs (A, B) of disjoint m-subsets of rem (|rem| = 2m here, so B is the complement) */
      T->npairs = 0;
      for (int s = 0; s < (1 << nr); s++) { if (__builtin_popcount(s) != m) continue;
        int a = 0, bb = 0;
        for (int i = 0; i < nr; i++) { if (s >> i & 1) T->A[T->npairs][a++] = rem[i]; else T->Bv[T->npairs][bb++] = rem[i]; }
        T->npairs++; }
    }
    ntup[q] = cnt;
  }
  /* the 28 six-sets as edge masks */
  long long six[28]; int ns = 0;
  for (int s = 0; s < 256; s++) if (__builtin_popcount(s) == 6) { long long em = 0;
    for (int y = 1; y < 8; y++) for (int x = 0; x < y; x++) if ((s >> x & 1) && (s >> y & 1)) em |= 1LL << PIDX(x, y);
    six[ns++] = em; }
  const long long FULL = (1LL << 28) - 1;
  for (int g = 0; g < ngraphs; g++) {
    long long H = gm[g]; if (H < 0 || H > FULL) die("graph mask out of range");
    int mono = 0; for (int i = 0; i < 28; i++) { long long x = H & six[i]; if (x == six[i] || x == 0) mono++; }
    i128 S[4] = {0, 0, 0, 0};
    for (int col = 0; col < 2; col++) {
      long long Gm = col ? (~H & FULL) : H;
      for (int q = 0; q < 4; q++) {
        int t = TS[q], m = MS[q], cp = pair_count(t);
        for (int u = 0; u < ntup[q]; u++) {
          Tuple *T = &tup[q][u];
          int rc = 0; for (int k = 0; k < cp; k++) rc |= (int)(Gm >> T->rootpair[k] & 1) << k;
          int b = lookup[col][t][rc]; if (b < 0) continue;
          Block *X = &B[b];
          for (int p = 0; p < T->npairs; p++) {
            int fi[2];
            for (int side = 0; side < 2; side++) {
              const int *E = side ? T->Bv[p] : T->A[p];
              int v[8]; for (int i = 0; i < t; i++) v[i] = T->r[i]; for (int i = 0; i < m; i++) v[t + i] = E[i];
              long long code = 0;
              for (int y = t; y < t + m; y++) for (int x = 0; x < y; x++) if (adjm(Gm, v[x], v[y])) code |= 1LL << PIDX(x, y);
              fi[side] = X->table[code >> cp];
            }
            if (fi[0] >= 0 && fi[1] >= 0) {
              i128 add = X->G[(size_t)fi[0] * X->nf + fi[1]];
              i128 before = S[q]; S[q] += add;
              if ((add > 0 && S[q] < before) || (add < 0 && S[q] > before)) die("128-bit overflow");
            }
          }
        }
      }
    }
    printf("%d %d ", g, mono); print128(S[0]); putchar(' '); print128(S[1]); putchar(' '); print128(S[2]); putchar(' ');
    print128(S[3]); putchar('\n');
  }
  return 0;
}

static int mode_aut(const char *graphs) {
  read_graphs(graphs);
  /* all 8! relabellings as maps of the 28 edge positions */
  static unsigned char emap[40320][28]; int np = 0; int a[8] = {0, 1, 2, 3, 4, 5, 6, 7};
  do { for (int y = 1; y < 8; y++) for (int x = 0; x < y; x++) {
         int u = a[x], w = a[y]; if (u > w) { int tmp = u; u = w; w = tmp; } emap[np][PIDX(x, y)] = (unsigned char)PIDX(u, w); }
       np++; } while (next_perm(a, 8));
  if (np != 40320) die("permutation count");
  for (int g = 0; g < ngraphs; g++) {
    long long H = gm[g], best = -1; int aut = 0; int bitsv[28], nbits = 0;
    for (int e = 0; e < 28; e++) if (H >> e & 1) bitsv[nbits++] = e;
    for (int p = 0; p < np; p++) {
      long long h2 = 0; for (int i = 0; i < nbits; i++) h2 |= 1LL << emap[p][bitsv[i]];
      if (h2 == H) aut++;
      if (best < 0 || h2 < best) best = h2;
    }
    printf("%d %lld %d\n", g, best, aut);
  }
  return 0;
}

int main(int argc, char **argv) {
  if (argc == 4 && !strcmp(argv[1], "sums")) return mode_sums(argv[2], argv[3]);
  if (argc == 3 && !strcmp(argv[1], "aut")) return mode_aut(argv[2]);
  fprintf(stderr, "usage: flag_count sums CERT GRAPHS | flag_count aut GRAPHS\n"); return 2;
}
