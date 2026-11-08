/* Combination step of our C45 replay (Romanoff's constant), written from the obstruction inequality on the 45a page.

   The base primes are split into five groups by the modulus their orders divide: G1 | 12, G2 | 36, G3 | 60, G4 | 48 and
   GA | 28, so that every order divides 5040 = lcm(12, 36, 60, 48, 28). For r in Z/5040 put t = r mod 12. Given t, the
   residues r = t (mod 12) correspond one to one (CRT) to tuples (r mod 12, r mod 36, r mod 60, r mod 48, r mod 28) with
   the given t, of sizes 1, 3, 5, 4, 7 (1*3*5*4*7 = 420 = 5040/12), and r is uncovered iff each component is uncovered by
   its group. So nu0 = sum_t prod_g u_g(t), where u_g(t) counts the uncovered residues of group g above t (for GA it
   depends on t mod 4 only). python's c45_check.py computes the exact distribution of each group's vector (u_g(t))_t;
   this program combines them into the histogram of nu0:
     p12  = G1 x G2 elementwise products, p123 = p12 x G3, both reduced to a canonical form under the 48 affine maps
            x -> u x + s of Z/12 (every remaining distribution is invariant under them: they lift to affine maps of
            Z/5040, which permute each group's configurations with their weights);
     V(a) = sum_{b mod 3} p123(t) v4(t), t = CRT(a mod 4, b mod 3), for every pair (p123, G4 vector), reduced under
            the 8 affine maps of Z/4 (GA's distribution is invariant under them);
     nu0  = sum_{a mod 4} V(a) vA(a), for every pair (V, GA vector), accumulated into hist[nu0].
   Weights are exact unsigned 128-bit integers (the caller checks that every total fits). Tables have at most 16M
   slots of 24 bytes; the V table is folded into hist (with GA) and cleared whenever it is 70% full.
   Input (stdin): five blocks "G1 n" / "G2 n" / "G3 n" / "G4 n" (12 values and a decimal weight per line) and "GA n"
   (4 values and a weight). Output: lines "nu count", and statistics on stderr. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef unsigned __int128 u128;

static u128 parse_u128(const char *s) { u128 x = 0; while (*s >= '0' && *s <= '9') x = x * 10 + (u128)(*s++ - '0'); return x; }
static void print_u128(FILE *f, u128 x) {
    char buf[64]; int i = 63; buf[i] = 0;
    if (x == 0) { fputs("0", f); return; }
    while (x) { buf[--i] = (char)('0' + (int)(x % 10)); x /= 10; }
    fputs(buf + i, f);
}

/* open addressing with keys and weights in separate arrays (24 bytes a slot); key 0 marks an empty slot, so keys are
   stored plus one */
typedef struct { uint64_t *k; u128 *w; size_t cap, n; } map;
static void map_init(map *m, size_t cap) {
    m->cap = cap; m->n = 0; m->k = calloc(cap, sizeof(uint64_t)); m->w = calloc(cap, sizeof(u128));
    if (!m->k || !m->w) { fprintf(stderr, "out of memory\n"); exit(2); }
}
static void map_free(map *m) { free(m->k); free(m->w); m->k = NULL; m->w = NULL; m->n = 0; }
static uint64_t mix(uint64_t x) { x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33; return x; }
static size_t CAP_LIMIT = (size_t)1 << 24;   /* 16M slots of 24 bytes: 400 MB per table (the memory cap) */
static void map_put(map *m, uint64_t key, u128 w) {   /* no growth */
    uint64_t k1 = key + 1;
    size_t i = mix(key) & (m->cap - 1);
    while (m->k[i] && m->k[i] != k1) i = (i + 1) & (m->cap - 1);
    if (!m->k[i]) { m->k[i] = k1; m->w[i] = 0; m->n++; }
    m->w[i] += w;
}
static void map_add(map *m, uint64_t key, u128 w) {   /* grows at load 0.7, up to CAP_LIMIT */
    if (10 * (m->n + 1) > 7 * m->cap) {
        if (m->cap * 2 > CAP_LIMIT) { fprintf(stderr, "a table would exceed %zu slots; stopping\n", CAP_LIMIT); exit(3); }
        map o = *m; map_init(m, o.cap * 2);
        for (size_t i = 0; i < o.cap; i++) if (o.k[i]) map_put(m, o.k[i] - 1, o.w[i]);
        map_free(&o);
    }
    map_put(m, key, w);
}

typedef struct { unsigned char v[12]; u128 w; } vec;
static vec *readblock(const char *name, int len, int *count) {
    char tag[8]; int n;
    if (scanf("%7s %d", tag, &n) != 2 || strcmp(tag, name)) { fprintf(stderr, "expected block %s\n", name); exit(2); }
    vec *a = calloc((size_t)n + 1, sizeof(vec));
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < len; k++) { int x; if (scanf("%d", &x) != 1 || x < 0 || x > 255) { fprintf(stderr, "bad value\n"); exit(2); } a[i].v[k] = (unsigned char)x; }
        char s[64]; if (scanf("%63s", s) != 1) { fprintf(stderr, "bad weight\n"); exit(2); }
        a[i].w = parse_u128(s);
    }
    *count = n; return a;
}

static int MAPS[48][12];
static uint64_t pack12(const unsigned char *v) { uint64_t k = 0; for (int i = 0; i < 12; i++) k = (k << 4) | v[i]; return k; }
static void unpack12(uint64_t k, unsigned char *v) { for (int i = 11; i >= 0; i--) { v[i] = k & 15; k >>= 4; } }
static uint64_t canon12(const unsigned char *v) {
    uint64_t best = UINT64_MAX; unsigned char w[12];
    for (int m = 0; m < 48; m++) { for (int x = 0; x < 12; x++) w[MAPS[m][x]] = v[x]; uint64_t k = pack12(w); if (k < best) best = k; }
    return best;
}

static int MAPS4[8][4];
static uint64_t canon4(const int *V) {   /* the 8 affine maps x -> u x + s of Z/4 (GA's distribution is invariant) */
    uint64_t best = UINT64_MAX; int w[4];
    for (int m = 0; m < 8; m++) {
        for (int x = 0; x < 4; x++) w[MAPS4[m][x]] = V[x];
        uint64_t k = 0; for (int x = 0; x < 4; x++) k = (k << 8) | (uint64_t)w[x];
        if (k < best) best = k;
    }
    return best;
}

static map combine(map *A, vec *B, int nb) {   /* canonical products of every entry of A with every vector of B */
    map out; map_init(&out, 1 << 16);
    unsigned char a[12], c[12];
    for (size_t i = 0; i < A->cap; i++) {
        if (!A->k[i]) continue;
        unpack12(A->k[i] - 1, a);
        for (int j = 0; j < nb; j++) {
            for (int x = 0; x < 12; x++) { int y = a[x] * B[j].v[x]; if (y > 15) { fprintf(stderr, "value above 15\n"); exit(2); } c[x] = (unsigned char)y; }
            map_add(&out, canon12(c), A->w[i] * B[j].w);
        }
    }
    return out;
}

static void flush(map *Vm, vec *GA, int nA, u128 *hist) {   /* fold the V table into hist with GA, then clear it */
    for (size_t k = 0; k < Vm->cap; k++) {
        if (!Vm->k[k]) continue;
        uint64_t key = Vm->k[k] - 1; int V[4];
        for (int a = 3; a >= 0; a--) { V[a] = (int)(key & 255); key >>= 8; }
        for (int j = 0; j < nA; j++) {
            int nu = 0;
            for (int a = 0; a < 4; a++) nu += V[a] * GA[j].v[a];
            if (nu > 5040) { fprintf(stderr, "nu above 5040\n"); exit(2); }
            hist[nu] += Vm->w[k] * GA[j].w;
        }
        Vm->k[k] = 0;
    }
    Vm->n = 0;
}

int main(void) {
    int n1, n2, n3, n4, nA;
    for (int u = 0, m = 0; u < 12; u++) {
        if (u != 1 && u != 5 && u != 7 && u != 11) continue;
        for (int s = 0; s < 12; s++, m++) for (int x = 0; x < 12; x++) MAPS[m][x] = (u * x + s) % 12;
    }
    for (int u = 1, m = 0; u < 4; u += 2) for (int s = 0; s < 4; s++, m++) for (int x = 0; x < 4; x++) MAPS4[m][x] = (u * x + s) % 4;
    vec *G1 = readblock("G1", 12, &n1), *G2 = readblock("G2", 12, &n2), *G3 = readblock("G3", 12, &n3);
    vec *G4 = readblock("G4", 12, &n4), *GA = readblock("GA", 4, &nA);
    map p1; map_init(&p1, 1 << 10);
    for (int i = 0; i < n1; i++) map_add(&p1, canon12(G1[i].v), G1[i].w);
    map p12 = combine(&p1, G2, n2);
    fprintf(stderr, "canonical states: p1 %zu, p12 %zu\n", p1.n, p12.n);
    map p123 = combine(&p12, G3, n3);
    fprintf(stderr, "canonical states: p123 %zu\n", p123.n);
    map_free(&p12);
    /* G4 rows: for each a mod 4 the three values at t = CRT(a, b), b = 0, 1, 2 */
    int T[4][3];
    for (int t = 0; t < 12; t++) T[t % 4][t % 3] = t;
    u128 *hist = calloc(5041, sizeof(u128));
    map Vm; map_init(&Vm, CAP_LIMIT);
    size_t flush_at = 7 * (CAP_LIMIT / 10), total_v = 0;
    int flushes = 0;
    unsigned char p[12];
    for (size_t i = 0; i < p123.cap; i++) {
        if (!p123.k[i]) continue;
        unpack12(p123.k[i] - 1, p);
        u128 wp = p123.w[i];
        for (int j = 0; j < n4; j++) {
            int V[4];
            for (int a = 0; a < 4; a++) {
                int s = 0;
                for (int b = 0; b < 3; b++) s += p[T[a][b]] * G4[j].v[T[a][b]];
                V[a] = s;
            }
            map_put(&Vm, canon4(V), wp * G4[j].w);
        }
        if (Vm.n > flush_at) { total_v += Vm.n; flush(&Vm, GA, nA, hist); flushes++; }
    }
    total_v += Vm.n; flush(&Vm, GA, nA, hist); flushes++;
    fprintf(stderr, "V vectors folded: %zu in %d flushes\n", total_v, flushes);
    map_free(&p123); map_free(&Vm);
    for (int nu = 0; nu <= 5040; nu++) if (hist[nu]) { printf("%d ", nu); print_u128(stdout, hist[nu]); printf("\n"); }
    return 0;
}
