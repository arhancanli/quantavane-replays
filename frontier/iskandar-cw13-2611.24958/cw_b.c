/* Program B (C, exact 64-bit bitsets): the thirteen binary constant-weight codes of P. R. A. Iskandar,
 * arXiv 2611.24958v1.  Written separately from program A (no shared code); the paper's verifiers were
 * not read or run.  Reads ORBIT_MANIFEST.csv for (n, d, w, M, blocks, rep_file) and the code file
 * codes/A{n}_{d}_{w}_{M}.txt.
 *
 * For each code: every line is a 0/1 string of length n (n <= 64) and weight w (popcount); the words
 * are sorted and must be distinct; there are exactly M; ALL PAIRS x, y satisfy |x & y| <= w - d/2
 * (distance 2(w - |x & y|) >= d); the code is rebuilt from the representatives under g (each block of
 * length > 1 shifted by one place) and, sorted, must equal the sorted word list, with orbit sizes
 * summing to M.
 * Usage: cw_b ANC_DIR [distance|repeat|weight]   (the control edits the first code after reading)
 * Prints 'B: VERDICT PASS' (exit 0) or 'B: VERDICT REFUSED' (exit 1).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef uint64_t u64;
static int fails = 0;
#define FAIL(...) do { if (fails < 20) { printf("B: FAIL "); printf(__VA_ARGS__); printf("\n"); } fails++; } while (0)
static int cmpu(const void *a, const void *b) { u64 x = *(const u64 *)a, y = *(const u64 *)b; return x < y ? -1 : x > y; }

static u64 *read_words(const char *path, int n, long *cnt, long *badlen) {
  FILE *f = fopen(path, "r"); if (!f) { perror(path); exit(2); }
  long cap = 1 << 17, k = 0; u64 *v = malloc(sizeof(u64) * cap); char line[512]; *badlen = 0;
  while (fgets(line, sizeof line, f)) {
    int L = (int)strcspn(line, "\r\n"); if (L == 0) continue;
    if (L != n) { (*badlen)++; continue; }
    u64 x = 0; int ok = 1;
    for (int i = 0; i < n; i++) { if (line[i] == '1') x |= 1ULL << i; else if (line[i] != '0') ok = 0; }
    if (!ok) { (*badlen)++; continue; }
    if (k == cap) { cap *= 2; v = realloc(v, sizeof(u64) * cap); }
    v[k++] = x;
  }
  fclose(f); *cnt = k; return v;
}

int main(int argc, char **argv) {
  if (argc < 2) { fprintf(stderr, "usage: %s ANC_DIR [control]\n", argv[0]); return 2; }
  const char *control = argc > 2 ? argv[2] : "";
  char path[4096]; snprintf(path, sizeof path, "%s/ORBIT_MANIFEST.csv", argv[1]);
  FILE *mf = fopen(path, "r"); if (!mf) { perror(path); return 2; }
  char line[4096]; if (!fgets(line, sizeof line, mf)) return 2;   /* header */
  int ci = 0; long totw = 0, totr = 0;
  while (fgets(line, sizeof line, mf)) {
    int n, d, w, M, old, ord; char blocks[512], rest[2048];
    if (sscanf(line, "%d,%d,%d,%d,%d,%d,%511[^,],%2047[^\n]", &n, &d, &w, &M, &old, &ord, blocks, rest) != 8) continue;
    /* rest = rep_count,"orbit sizes" or plain,rep_file,sha,sha : take the field ending in .txt */
    char *rf = strstr(rest, "reps/"); if (!rf) { FAIL("manifest row without rep file"); continue; }
    char repfile[512]; sscanf(rf, "%511[^,]", repfile);
    snprintf(path, sizeof path, "%s/codes/A%d_%d_%d_%d.txt", argv[1], n, d, w, M);
    long k, badlen; u64 *v = read_words(path, n, &k, &badlen);
    if (ci == 0 && *control) {
      u64 x = v[0];
      if (!strcmp(control, "distance")) { int i1 = __builtin_ctzll(x); int i0 = __builtin_ctzll(~x); v[1] = (x & ~(1ULL << i1)) | (1ULL << i0); }
      else if (!strcmp(control, "repeat")) v[1] = x;
      else if (!strcmp(control, "weight")) v[1] = x | (1ULL << __builtin_ctzll(~x));
      printf("B: CONTROL %s: A(%d,%d,%d) word 2 replaced (bitmask %llx)\n", control, n, d, w, (unsigned long long)v[1]);
    }
    if (badlen) FAIL("A(%d,%d,%d): %ld lines not 0/1 strings of length %d", n, d, w, badlen, n);
    long badw = 0; for (long i = 0; i < k; i++) if (__builtin_popcountll(v[i]) != w) badw++;
    if (badw) FAIL("A(%d,%d,%d): %ld words of weight != %d", n, d, w, badw, w);
    if (k + badlen != M) FAIL("A(%d,%d,%d): %ld words, claimed %d", n, d, w, k + badlen, M);
    u64 *s = malloc(sizeof(u64) * k); memcpy(s, v, sizeof(u64) * k); qsort(s, k, sizeof(u64), cmpu);
    long dup = 0; for (long i = 1; i < k; i++) if (s[i] == s[i - 1]) dup++;
    if (dup) FAIL("A(%d,%d,%d): %ld repeated words", n, d, w, dup);
    /* pairwise intersections.  A pair with |x & y| > lim >= 1 shares at least two coordinates, so it is met in
       the bucket of words containing any two of its common coordinates: every word is put in the bucket of each
       pair (i, j) of its coordinates and all pairs inside each bucket are intersected.  For M <= 30000 every
       pair of words is also intersected directly (all pairs) and the two maxima must agree. */
    int lim = w - d / 2; long viol = 0; int maxint = 0;
    if (lim < 1) { FAIL("bucket method needs w - d/2 >= 1"); }
    long *cnt = calloc(64 * 64 + 1, sizeof(long));
    for (long i = 0; i < k; i++) { u64 x = v[i]; for (int a = 0; a < n; a++) if (x >> a & 1) for (int b = a + 1; b < n; b++) if (x >> b & 1) cnt[a * 64 + b + 1]++; }
    for (int q = 1; q <= 64 * 64; q++) cnt[q] += cnt[q - 1];
    long *fill = calloc(64 * 64, sizeof(long)); u64 *bk = malloc(sizeof(u64) * (cnt[64 * 64] + 1));
    for (long i = 0; i < k; i++) { u64 x = v[i]; for (int a = 0; a < n; a++) if (x >> a & 1) for (int b = a + 1; b < n; b++) if (x >> b & 1) { int q = a * 64 + b; bk[cnt[q] + fill[q]++] = x; } }
    long bpairs = 0;
    for (int q = 0; q < 64 * 64; q++) {
      long lo = cnt[q], hi = cnt[q + 1];
      for (long a = lo; a < hi; a++) for (long b = a + 1; b < hi; b++) {
        int c = __builtin_popcountll(bk[a] & bk[b]); bpairs++;
        if (c > maxint) maxint = c;
        if (c > lim) viol++;
      }
    }
    if (maxint < 1 && k > 1) maxint = 1;   /* pairs outside every bucket meet in at most one point */
    free(cnt); free(fill); free(bk);
    int allmax = -1;
    if (k <= 30000) {
      allmax = 0;
      for (long i = 0; i < k; i++) { u64 x = v[i]; for (long j = i + 1; j < k; j++) { int c = __builtin_popcountll(x & v[j]); if (c > allmax) allmax = c; } }
      if (allmax > maxint && allmax >= 2) FAIL("A(%d,%d,%d): all-pairs maximum %d differs from bucket maximum %d", n, d, w, allmax, maxint);
      if (allmax > lim) viol++;
    }
    if (viol) FAIL("A(%d,%d,%d): %ld word pairs found (counted once per shared bucket) meeting in more than %d points (distance < %d)", n, d, w, viol, lim, d);
    /* orbit rebuild */
    int perm[64], st = 0; char bcopy[512]; strcpy(bcopy, blocks);
    for (char *t = strtok(bcopy, " "); t; t = strtok(0, " ")) { int b = atoi(t); for (int j = 0; j < b; j++) perm[st + j] = st + (j + 1) % b; st += b; }
    if (st != n) FAIL("A(%d,%d,%d): blocks sum to %d", n, d, w, st);
    snprintf(path, sizeof path, "%s/%s", argv[1], repfile);
    long nr, badr; u64 *r = read_words(path, n, &nr, &badr);
    u64 *o = malloc(sizeof(u64) * (M + 64 * nr + 1)); long no = 0;
    for (long i = 0; i < nr; i++) {
      u64 x = r[i];
      do { o[no++] = x; u64 y = 0; for (int b = 0; b < n; b++) if (x >> b & 1) y |= 1ULL << perm[b]; x = y; } while (x != r[i] && no < M + 64 * nr);
    }
    qsort(o, no, sizeof(u64), cmpu);
    long ou = 0; for (long i = 0; i < no; i++) if (i == 0 || o[i] != o[i - 1]) o[ou++] = o[i];
    int same = (ou == k - dup) && (no == ou) && (no == M);
    if (same) { long u = 0; for (long i = 0; i < k; i++) if (i == 0 || s[i] != s[i - 1]) { if (s[i] != o[u]) { same = 0; break; } u++; } }
    if (!same) FAIL("A(%d,%d,%d): orbit rebuild (%ld words, %ld distinct) differs from the list", n, d, w, no, ou);
    printf("B: A(%d,%d,%d) >= %d: words %ld, repeats %ld, weight ok %s, largest intersection of two words %d (bound %d; minimum distance %d; %ld bucket pairs%s), orbit rebuild from %ld reps under blocks [%s]: %ld words, equals list: %s\n",
           n, d, w, M, k, dup, badw ? "NO" : "yes", maxint, lim, 2 * (w - maxint), bpairs, allmax >= 0 ? ", all pairs agree" : "", nr, blocks, no, same ? "yes" : "NO");
    fflush(stdout);
    totw += k; totr += nr; ci++;
    free(v); free(s); free(r); free(o);
  }
  fclose(mf);
  printf("B: codes %d, total words %ld, total representatives %ld\n", ci, totw, totr);
  if (fails) { printf("B: VERDICT REFUSED (%d failures)\n", fails); return 1; }
  printf("B: VERDICT PASS\n");
  return 0;
}
