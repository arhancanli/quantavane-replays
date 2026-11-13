/* Base histogram of nu0 on Z/T0, T0 = 12252240 = 2^4 3^2 5 7 11 13 17, for our C45 replays of the 29- to 34-prime
   bases (rows [S2026] of PR #290 and [P2026b], [P2026c] of PR #281), written from the obstruction inequality on the
   45a page. Our own code; no program of any certificate was read or run.

   The base primes are grouped by the part of their order outside 12. Groups G1 (orders | 12), G9 (| 36), G5 (| 60)
   and G16 (| 48) see t = r mod 12; groups whose orders are prime to 3 (G7 | 28, G13 | 52, G11 | 44, G17 | 68) see
   only s = r mod 4. Given t, the residues r mod T0 correspond (CRT) to independent private coordinates of the
   groups, and r is uncovered iff each group leaves its coordinate uncovered, so
       nu0 = sum_t prod_g u_g(t) = sum_{s mod 4} V(s) B(s),   V(s) = sum_{t = s (mod 4)} (u1 u9 u5 u16)(t),
   B(s) = the product of the s-groups' u_g(s). c45_frame_check.py computes each group's exact law; this program
     p123 = canonical products of G1, G9, G5 under the 48 affine maps of Z/12 (as in c45_comb.c);
     V    = for every (p123, G16 vector) the 4-vector V, reduced under the 8 affine maps x -> ux + c of Z/4;
     B    = the full law of the s-groups' product (every vector; B's law is invariant under the same 8 maps, which the
            program checks, so reducing V alone is exact);
     hist[nu] += w_V w_B for every pair (V, B) with nu = V.B <= C, with 256-bit counts.
   Pairs with nu > C are not enumerated: a k-d tree on B stores each node's componentwise minimum and total weight,
   and a node is skipped when V.(minimum) > C; its weight is still added to the V's skipped weight, and for every V
   the program checks accepted + skipped = the total weight of B. It also returns the exact first and second moments
   of nu0 over all pairs (from the moments of B), which the Python driver compares with closed forms.
   V vectors with entry sum < SIG (default 200), whose pairs are most of the work, are kept in their own table and
   processed once, after the stream; the others are processed whenever their table is 70% full. To bound memory the
   low vectors are split into PASSES classes by a hash of the vector; pass k replays the stream and keeps class k only
   (the high vectors are processed in pass 0). Work is split over THREADS threads with private histograms. State is
   checkpointed in STATEDIR at the end of each pass's stream and every five minutes while its low table is processed,
   so a stopped run resumes there (a pass's stream is replayed if it was stopped during it).

   Usage: c45_frame C THREADS STATEDIR [SIG [PASSES]] < input      input: blocks "G1 n" "G9 n" "G5 n" "G16 n" (12 values and a
   decimal weight per line), then "NB k" and k blocks "name n" (4 values and a weight). Output (stdout): "nu count"
   for 0 <= nu <= C, then lines "# above count", "# moment1 ...", "# moment2 ...", "# total ..."; stderr: progress. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <pthread.h>
#include <unistd.h>

typedef unsigned __int128 u128;
typedef struct { uint64_t w[4]; } u256;

static u128 parse_u128(const char *s) { u128 x = 0; while (*s >= '0' && *s <= '9') x = x * 10 + (u128)(*s++ - '0'); return x; }
static void print_u256(FILE *f, const u256 *x) {
    uint64_t t[4]; memcpy(t, x->w, sizeof t); char buf[100]; int i = 99; buf[i] = 0;
    if (!(t[0] | t[1] | t[2] | t[3])) { fputs("0", f); return; }
    while (t[0] | t[1] | t[2] | t[3]) {           /* divide by 10 */
        u128 r = 0;
        for (int k = 3; k >= 0; k--) { u128 cur = (r << 64) | t[k]; t[k] = (uint64_t)(cur / 10); r = cur % 10; }
        buf[--i] = (char)('0' + (int)r);
    }
    fputs(buf + i, f);
}
static inline void add_u256(u256 *a, const u256 *b) {
    u128 c = (u128)a->w[0] + b->w[0]; a->w[0] = (uint64_t)c;
    c = (c >> 64) + a->w[1] + b->w[1]; a->w[1] = (uint64_t)c;
    c = (c >> 64) + a->w[2] + b->w[2]; a->w[2] = (uint64_t)c;
    c = (c >> 64) + a->w[3] + b->w[3]; a->w[3] = (uint64_t)c;
    if (c >> 64) { fprintf(stderr, "256-bit overflow\n"); exit(4); }
}
/* a += x * y, x and y below 2^128 */
static inline void addmul(u256 *a, u128 x, u128 y) {
    uint64_t x0 = (uint64_t)x, x1 = (uint64_t)(x >> 64), y0 = (uint64_t)y, y1 = (uint64_t)(y >> 64);
    u128 p00 = (u128)x0 * y0, p01 = (u128)x0 * y1, p10 = (u128)x1 * y0, p11 = (u128)x1 * y1;
    u128 mid = p01 + p10; uint64_t midc = (uint64_t)(mid < p01);
    u128 c = (u128)a->w[0] + (uint64_t)p00; a->w[0] = (uint64_t)c;
    c = (c >> 64) + a->w[1] + (uint64_t)(p00 >> 64) + (uint64_t)mid; a->w[1] = (uint64_t)c;
    c = (c >> 64) + a->w[2] + (uint64_t)(mid >> 64) + (uint64_t)p11; a->w[2] = (uint64_t)c;
    c = (c >> 64) + a->w[3] + (uint64_t)(p11 >> 64) + midc; a->w[3] = (uint64_t)c;
    if (c >> 64) { fprintf(stderr, "256-bit overflow\n"); exit(4); }
}
/* a += x * y, x a 256-bit number, y below 2^128 (the product must fit 256 bits) */
static void addmul_wide(u256 *a, const u256 *x, u128 y) {
    uint64_t yy[2] = {(uint64_t)y, (uint64_t)(y >> 64)}, r[6] = {0, 0, 0, 0, 0, 0};
    for (int i = 0; i < 4; i++) {
        u128 c = 0;
        for (int j = 0; j < 2; j++) { c += (u128)x->w[i] * yy[j] + r[i + j]; r[i + j] = (uint64_t)c; c >>= 64; }
        for (int k = i + 2; c && k < 6; k++) { c += r[k]; r[k] = (uint64_t)c; c >>= 64; }
    }
    if (r[4] | r[5]) { fprintf(stderr, "256-bit overflow\n"); exit(4); }
    u256 t = {{r[0], r[1], r[2], r[3]}}; add_u256(a, &t);
}
static void mul_u256_small(u256 *a, uint64_t m) {   /* a *= m */
    u128 c = 0;
    for (int k = 0; k < 4; k++) { c += (u128)a->w[k] * m; a->w[k] = (uint64_t)c; c >>= 64; }
    if (c) { fprintf(stderr, "256-bit overflow\n"); exit(4); }
}

/* open addressing, keys stored plus one */
typedef struct { uint64_t *k; u128 *w; size_t cap, n; } map;
static void map_init(map *m, size_t cap) {
    m->cap = cap; m->n = 0; m->k = calloc(cap, sizeof(uint64_t)); m->w = calloc(cap, sizeof(u128));
    if (!m->k || !m->w) { fprintf(stderr, "out of memory\n"); exit(2); }
}
static void map_free(map *m) { free(m->k); free(m->w); m->k = NULL; m->w = NULL; m->n = 0; }
static uint64_t mix(uint64_t x) { x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33; return x; }
static void map_put(map *m, uint64_t key, u128 w) {
    uint64_t k1 = key + 1; size_t i = mix(key) & (m->cap - 1);
    while (m->k[i] && m->k[i] != k1) i = (i + 1) & (m->cap - 1);
    if (!m->k[i]) { m->k[i] = k1; m->w[i] = 0; m->n++; }
    m->w[i] += w;
}
static u128 map_get(const map *m, uint64_t key) {
    uint64_t k1 = key + 1; size_t i = mix(key) & (m->cap - 1);
    while (m->k[i]) { if (m->k[i] == k1) return m->w[i]; i = (i + 1) & (m->cap - 1); }
    return 0;
}
static size_t GROW_LIMIT = (size_t)1 << 25;
static void map_add(map *m, uint64_t key, u128 w) {
    if (10 * (m->n + 1) > 7 * m->cap) {
        if (m->cap * 2 > GROW_LIMIT) { fprintf(stderr, "a table would exceed %zu slots; stopping\n", GROW_LIMIT); exit(3); }
        map o = *m; map_init(m, o.cap * 2);
        for (size_t i = 0; i < o.cap; i++) if (o.k[i]) map_put(m, o.k[i] - 1, o.w[i]);
        map_free(&o);
    }
    map_put(m, key, w);
}

typedef struct { unsigned char v[12]; u128 w; } vec;
static vec *readblock(const char *name, int len, int *count) {
    char tag[16]; int n;
    if (scanf("%15s %d", tag, &n) != 2 || (name && strcmp(tag, name))) { fprintf(stderr, "expected block %s\n", name ? name : "(any)"); exit(2); }
    vec *a = calloc((size_t)n + 1, sizeof(vec));
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < len; k++) { int x; if (scanf("%d", &x) != 1 || x < 0 || x > 255) { fprintf(stderr, "bad value\n"); exit(2); } a[i].v[k] = (unsigned char)x; }
        char s[64]; if (scanf("%63s", s) != 1) { fprintf(stderr, "bad weight\n"); exit(2); }
        a[i].w = parse_u128(s);
    }
    *count = n; return a;
}

static int MAPS[48][12], MAPS4[8][4];
static uint64_t pack12(const unsigned char *v) { uint64_t k = 0; for (int i = 0; i < 12; i++) k = (k << 4) | v[i]; return k; }
static void unpack12(uint64_t k, unsigned char *v) { for (int i = 11; i >= 0; i--) { v[i] = k & 15; k >>= 4; } }
static uint64_t canon12(const unsigned char *v) {
    uint64_t best = UINT64_MAX; unsigned char w[12];
    for (int m = 0; m < 48; m++) { for (int x = 0; x < 12; x++) w[MAPS[m][x]] = v[x]; uint64_t k = pack12(w); if (k < best) best = k; }
    return best;
}
static uint64_t canon4(const int *V) {
    uint64_t best = UINT64_MAX; int w[4];
    for (int m = 0; m < 8; m++) {
        for (int x = 0; x < 4; x++) w[MAPS4[m][x]] = V[x];
        uint64_t k = 0; for (int x = 0; x < 4; x++) k = (k << 8) | (uint64_t)w[x];
        if (k < best) best = k;
    }
    return best;
}
static map combine(map *A, vec *B, int nb) {
    map out; map_init(&out, 1 << 16); unsigned char a[12], c[12];
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

/* ---- B side: full law and k-d tree ---- */
static int NB; static uint32_t *BX[4]; static u128 *BW;
typedef struct { uint32_t mn[4]; int lo, hi, left, right; u128 sw; } node;
static node *NODES; static int NNODES;
static u128 WB_TOTAL;
static int cmp_dim;
static int *ORD;
static int cmp_idx(const void *a, const void *b) {
    uint32_t x = BX[cmp_dim][*(const int *)a], y = BX[cmp_dim][*(const int *)b];
    return x < y ? -1 : x > y;
}
static int build(int lo, int hi) {   /* over ORD[lo..hi) */
    int id = NNODES++;
    node *nd = &NODES[id];
    nd->lo = lo; nd->hi = hi; nd->left = nd->right = -1; nd->sw = 0;
    uint32_t mn[4] = {UINT32_MAX, UINT32_MAX, UINT32_MAX, UINT32_MAX}, mx[4] = {0, 0, 0, 0};
    for (int i = lo; i < hi; i++) for (int d = 0; d < 4; d++) { uint32_t x = BX[d][ORD[i]]; if (x < mn[d]) mn[d] = x; if (x > mx[d]) mx[d] = x; }
    for (int d = 0; d < 4; d++) nd->mn[d] = mn[d];
    for (int i = lo; i < hi; i++) nd->sw += BW[ORD[i]];
    if (hi - lo > 16) {
        int bd = 0; for (int d = 1; d < 4; d++) if (mx[d] - mn[d] > mx[bd] - mn[bd]) bd = d;
        cmp_dim = bd; qsort(ORD + lo, (size_t)(hi - lo), sizeof(int), cmp_idx);
        int mid = (lo + hi) / 2;
        int L = build(lo, mid), R = build(mid, hi);
        NODES[id].left = L; NODES[id].right = R;
    }
    return id;
}

/* ---- workers ---- */
static long long CUT; static int THREADS = 1; static int SIG = 200; static int PASSES = 1; static int STOP_AFTER = -1;
typedef struct { uint64_t key; u128 w; } ventry;
typedef struct {
    const ventry *v; size_t n;
    u256 *hist;                      /* private histogram, C + 1 bins */
    u256 below, above;              /* weight of pairs with nu <= C and nu > C */
    u256 m1, m2;                     /* sum w nu and sum w nu^2 over all pairs, from B's moments */
    long long bad;                   /* V's whose accepted + skipped weight differs from W_B */
    long long pairs;
} job;
static u256 MB1;                     /* sum_B w_B B(0) (the same for every s) */
static u256 MB2[4][4];               /* sum_B w_B B(s) B(t) */
static uint32_t ROOTMIN;

static void *worker(void *arg) {
    job *J = arg;
    int *stack = malloc(sizeof(int) * 4096);
    for (size_t e = 0; e < J->n; e++) {
        uint64_t key = J->v[e].key; u128 wv = J->v[e].w; int V[4];
        for (int a = 3; a >= 0; a--) { V[a] = (int)(key & 255); key >>= 8; }
        long long sig = (long long)V[0] + V[1] + V[2] + V[3];
        /* moments of nu over the pairs of this V, from B's moments (exact, B's law invariant) */
        {
            u256 t1 = MB1; mul_u256_small(&t1, (uint64_t)sig); addmul_wide(&J->m1, &t1, wv);
            u256 acc = {{0, 0, 0, 0}};
            for (int s = 0; s < 4; s++) for (int t = 0; t < 4; t++) { u256 x = MB2[s][t]; mul_u256_small(&x, (uint64_t)(V[s] * V[t])); add_u256(&acc, &x); }
            addmul_wide(&J->m2, &acc, wv);
        }
        u128 acc_lo = 0, acc_hi = 0;
        if ((long long)sig * ROOTMIN > CUT) { acc_hi = WB_TOTAL; }
        else {
            int sp = 0; stack[sp++] = 0;
            while (sp) {
                const node *nd = &NODES[stack[--sp]];
                long long lb = (long long)V[0] * nd->mn[0] + (long long)V[1] * nd->mn[1] + (long long)V[2] * nd->mn[2] + (long long)V[3] * nd->mn[3];
                if (lb > CUT) { acc_hi += nd->sw; continue; }
                if (nd->left < 0) {
                    for (int i = nd->lo; i < nd->hi; i++) {
                        long long nu = (long long)V[0] * BX[0][i] + (long long)V[1] * BX[1][i] + (long long)V[2] * BX[2][i] + (long long)V[3] * BX[3][i];
                        if (nu <= CUT) { addmul(&J->hist[nu], wv, BW[i]); acc_lo += BW[i]; J->pairs++; }
                        else acc_hi += BW[i];
                    }
                } else { stack[sp++] = nd->left; stack[sp++] = nd->right; if (sp > 4090) { fprintf(stderr, "stack\n"); exit(5); } }
            }
        }
        if (acc_lo + acc_hi != WB_TOTAL) J->bad++;
        addmul(&J->below, wv, acc_lo);
        addmul(&J->above, wv, acc_hi);
    }
    free(stack);
    return NULL;
}

static u256 *HIST; static u256 BELOW, ABOVE, M1, M2; static long long BAD, PAIRS;
static u256 *THIST[64];
static void run_batch(const ventry *v, size_t n) {
    pthread_t th[64]; job J[64];
    size_t per = (n + THREADS - 1) / THREADS;
    for (int t = 0; t < THREADS; t++) {
        size_t a = t * per, b = a + per > n ? n : a + per; if (a > n) a = n;
        memset(&J[t], 0, sizeof(job)); J[t].v = v + a; J[t].n = b - a; J[t].hist = THIST[t];
        pthread_create(&th[t], NULL, worker, &J[t]);
    }
    for (int t = 0; t < THREADS; t++) {
        pthread_join(th[t], NULL);
        add_u256(&BELOW, &J[t].below); add_u256(&ABOVE, &J[t].above); add_u256(&M1, &J[t].m1); add_u256(&M2, &J[t].m2);
        BAD += J[t].bad; PAIRS += J[t].pairs;
    }
}
static void merge_hists(void) {
    for (int t = 0; t < THREADS; t++) {
        for (long long nu = 0; nu <= CUT; nu++) if (THIST[t][nu].w[0] | THIST[t][nu].w[1] | THIST[t][nu].w[2] | THIST[t][nu].w[3]) { add_u256(&HIST[nu], &THIST[t][nu]); memset(&THIST[t][nu], 0, sizeof(u256)); }
    }
}

/* ---- state ---- */
static char STATEDIR[1024];
static uint64_t INPUT_HASH;
static void save_state(const char *phase, long long chunk) {
    char tmp[1200], fin[1200];
    snprintf(tmp, sizeof tmp, "%s/state.tmp", STATEDIR); snprintf(fin, sizeof fin, "%s/state.bin", STATEDIR);
    FILE *f = fopen(tmp, "wb"); if (!f) { fprintf(stderr, "cannot write state\n"); exit(6); }
    char ph[8] = {0}; strncpy(ph, phase, 7);
    fwrite(&INPUT_HASH, 8, 1, f); fwrite(&CUT, 8, 1, f); fwrite(ph, 8, 1, f); fwrite(&chunk, 8, 1, f);
    fwrite(&BELOW, sizeof(u256), 1, f); fwrite(&ABOVE, sizeof(u256), 1, f); fwrite(&M1, sizeof(u256), 1, f); fwrite(&M2, sizeof(u256), 1, f);
    fwrite(&BAD, 8, 1, f); fwrite(&PAIRS, 8, 1, f);
    fwrite(HIST, sizeof(u256), (size_t)CUT + 1, f);
    if (fclose(f)) { fprintf(stderr, "cannot write state\n"); exit(6); }
    rename(tmp, fin);
}
static int load_state(char *phase, long long *chunk) {
    char fin[1200]; snprintf(fin, sizeof fin, "%s/state.bin", STATEDIR);
    FILE *f = fopen(fin, "rb"); if (!f) return 0;
    uint64_t h; long long c; char ph[8];
    if (fread(&h, 8, 1, f) != 1 || fread(&c, 8, 1, f) != 1 || h != INPUT_HASH || c != CUT) { fclose(f); fprintf(stderr, "state is for another input; ignored\n"); return 0; }
    if (fread(ph, 8, 1, f) != 1 || fread(chunk, 8, 1, f) != 1) { fclose(f); return 0; }
    int ok = fread(&BELOW, sizeof(u256), 1, f) == 1 && fread(&ABOVE, sizeof(u256), 1, f) == 1 && fread(&M1, sizeof(u256), 1, f) == 1 && fread(&M2, sizeof(u256), 1, f) == 1
             && fread(&BAD, 8, 1, f) == 1 && fread(&PAIRS, 8, 1, f) == 1 && fread(HIST, sizeof(u256), (size_t)CUT + 1, f) == (size_t)CUT + 1;
    fclose(f);
    if (!ok) return 0;
    memcpy(phase, ph, 8);
    return 1;
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: c45_frame C THREADS STATEDIR [SIG] < input\n"); return 2; }
    CUT = atoll(argv[1]); THREADS = atoi(argv[2]); snprintf(STATEDIR, sizeof STATEDIR, "%s", argv[3]);
    if (argc > 4) SIG = atoi(argv[4]);
    if (argc > 5) PASSES = atoi(argv[5]);
    if (argc > 6) STOP_AFTER = atoi(argv[6]);      /* testing: stop after this pass, as if killed */
    if (PASSES < 1 || PASSES > 64 || (PASSES & (PASSES - 1))) { fprintf(stderr, "PASSES must be a power of 2 up to 64\n"); return 2; }
    if (THREADS < 1 || THREADS > 64 || CUT < 1) { fprintf(stderr, "bad arguments\n"); return 2; }
    time_t T0 = time(NULL);
    for (int u = 0, m = 0; u < 12; u++) {
        if (u != 1 && u != 5 && u != 7 && u != 11) continue;
        for (int s = 0; s < 12; s++, m++) for (int x = 0; x < 12; x++) MAPS[m][x] = (u * x + s) % 12;
    }
    for (int u = 1, m = 0; u < 4; u += 2) for (int s = 0; s < 4; s++, m++) for (int x = 0; x < 4; x++) MAPS4[m][x] = (u * x + s) % 4;
    int n1, n2, n3, n4;
    vec *G1 = readblock("G1", 12, &n1), *G2 = readblock("G9", 12, &n2), *G3 = readblock("G5", 12, &n3), *G4 = readblock("G16", 12, &n4);
    int ng; char tag[16];
    if (scanf("%15s %d", tag, &ng) != 2 || strcmp(tag, "NB")) { fprintf(stderr, "expected NB\n"); return 2; }
    /* input hash (for the state file): sizes and every weight */
    INPUT_HASH = 1469598103934665603ULL;
    vec *GS[4] = {G1, G2, G3, G4}; int NS[4] = {n1, n2, n3, n4};
    for (int g = 0; g < 4; g++) for (int i = 0; i < NS[g]; i++) { INPUT_HASH = mix(INPUT_HASH ^ (uint64_t)GS[g][i].w ^ pack12(GS[g][i].v)); }
    /* B law */
    map Bm; map_init(&Bm, 1 << 10); map_put(&Bm, ((uint64_t)1 << 48) | ((uint64_t)1 << 32) | ((uint64_t)1 << 16) | 1, 1);
    for (int g = 0; g < ng; g++) {
        int n; vec *G = readblock(NULL, 4, &n); map nb; map_init(&nb, 1 << 10);
        for (int i = 0; i < n; i++) INPUT_HASH = mix(INPUT_HASH ^ (uint64_t)G[i].w ^ ((uint64_t)G[i].v[0] << 24 | G[i].v[1] << 16 | G[i].v[2] << 8 | G[i].v[3]) ^ ((uint64_t)g << 56));
        for (size_t i = 0; i < Bm.cap; i++) {
            if (!Bm.k[i]) continue;
            uint64_t k = Bm.k[i] - 1, b[4]; for (int s = 3; s >= 0; s--) { b[s] = k & 0xffff; k >>= 16; }
            for (int j = 0; j < n; j++) {
                uint64_t c = 0;
                for (int s = 0; s < 4; s++) { uint64_t y = b[s] * G[j].v[s]; if (y > 0xffff) { fprintf(stderr, "B value above 65535\n"); return 2; } c = (c << 16) | y; }
                map_add(&nb, c, Bm.w[i] * G[j].w);
            }
        }
        map_free(&Bm); Bm = nb; free(G);
    }
    /* invariance of B's law under the 8 affine maps of Z/4 */
    long long noninv = 0;
    NB = (int)Bm.n;
    for (int d = 0; d < 4; d++) BX[d] = malloc(sizeof(uint32_t) * NB);
    BW = malloc(sizeof(u128) * NB); ORD = malloc(sizeof(int) * NB);
    uint32_t *tx[4]; for (int d = 0; d < 4; d++) tx[d] = malloc(sizeof(uint32_t) * NB);
    u128 *tw = malloc(sizeof(u128) * NB);
    {
        int j = 0;
        for (size_t i = 0; i < Bm.cap; i++) {
            if (!Bm.k[i]) continue;
            uint64_t k = Bm.k[i] - 1; uint32_t b[4];
            for (int s = 3; s >= 0; s--) { b[s] = (uint32_t)(k & 0xffff); k >>= 16; }
            for (int m = 0; m < 8; m++) {
                uint32_t w[4]; for (int x = 0; x < 4; x++) w[MAPS4[m][x]] = b[x];
                uint64_t c = ((uint64_t)w[0] << 48) | ((uint64_t)w[1] << 32) | ((uint64_t)w[2] << 16) | w[3];
                if (map_get(&Bm, c) != Bm.w[i]) noninv++;
            }
            for (int d = 0; d < 4; d++) tx[d][j] = b[d];
            tw[j] = Bm.w[i]; ORD[j] = j; j++;
        }
    }
    for (int d = 0; d < 4; d++) BX[d] = tx[d];
    BW = tw;
    map_free(&Bm);
    WB_TOTAL = 0; memset(&MB1, 0, sizeof MB1); memset(MB2, 0, sizeof MB2);
    for (int i = 0; i < NB; i++) {
        if (WB_TOTAL + BW[i] < WB_TOTAL) { fprintf(stderr, "B weight overflow\n"); return 4; }
        WB_TOTAL += BW[i]; addmul(&MB1, BW[i], BX[0][i]);
        for (int s = 0; s < 4; s++) for (int t = 0; t < 4; t++) addmul(&MB2[s][t], BW[i], (u128)BX[s][i] * BX[t][i]);
    }
    NODES = malloc(sizeof(node) * (size_t)(NB / 4 + 64)); NNODES = 0;
    build(0, NB);
    /* store B in tree order */
    {
        uint32_t *nx[4]; for (int d = 0; d < 4; d++) nx[d] = malloc(sizeof(uint32_t) * NB);
        u128 *nw = malloc(sizeof(u128) * NB);
        for (int i = 0; i < NB; i++) { for (int d = 0; d < 4; d++) nx[d][i] = BX[d][ORD[i]]; nw[i] = BW[ORD[i]]; }
        for (int d = 0; d < 4; d++) { free(BX[d]); BX[d] = nx[d]; }
        free(BW); BW = nw;
    }
    ROOTMIN = NODES[0].mn[0]; for (int d = 1; d < 4; d++) if (NODES[0].mn[d] < ROOTMIN) ROOTMIN = NODES[0].mn[d];
    fprintf(stderr, "B law: %d vectors, %d tree nodes, least entry %u; invariant under the 8 affine maps of Z/4: %s\n", NB, NNODES, ROOTMIN, noninv ? "NO" : "yes");
    if (noninv) return 7;
    HIST = calloc((size_t)CUT + 1, sizeof(u256));
    for (int t = 0; t < THREADS; t++) { THIST[t] = calloc((size_t)CUT + 1, sizeof(u256)); if (!THIST[t]) { fprintf(stderr, "out of memory\n"); return 2; } }
    if (!HIST) { fprintf(stderr, "out of memory\n"); return 2; }
    char lowpath[1200]; snprintf(lowpath, sizeof lowpath, "%s/low.bin", STATEDIR);
    INPUT_HASH = mix(INPUT_HASH ^ ((uint64_t)SIG << 32) ^ (uint64_t)PASSES);
    char phase[8] = {0}; long long st = 0;
    int resumed = load_state(phase, &st);   /* st = (pass << 32) | chunk; phase "low": that pass's low table is saved */
    int pass0 = 0; long long chunk0 = 0, have_low = 0;
    if (resumed && !strcmp(phase, "low")) { pass0 = (int)(st >> 32); chunk0 = st & 0xffffffffLL; have_low = 1; }
    else if (resumed && !strcmp(phase, "pass")) { pass0 = (int)(st >> 32); }
    else { memset(&BELOW, 0, sizeof BELOW); memset(&ABOVE, 0, sizeof ABOVE); memset(&M1, 0, sizeof M1); memset(&M2, 0, sizeof M2); BAD = 0; PAIRS = 0; }
    if (resumed) fprintf(stderr, "resumed at pass %d, phase %s, chunk %lld\n", pass0, phase, chunk0);
    int T[4][3]; for (int t = 0; t < 12; t++) T[t % 4][t % 3] = t;
    for (int pass = pass0; pass < PASSES; pass++) {
        ventry *low = NULL; size_t nlow = 0;
        if (have_low && pass == pass0) {
            FILE *f = fopen(lowpath, "rb");
            if (!f) { fprintf(stderr, "state without its low table; start again from scratch\n"); return 6; }
            if (fread(&nlow, sizeof nlow, 1, f) != 1) return 6;
            low = malloc(sizeof(ventry) * (nlow + 1));
            if (fread(low, sizeof(ventry), nlow, f) != nlow) { fprintf(stderr, "short low table\n"); return 6; }
            fclose(f);
        } else {
            chunk0 = 0;
            map p1; map_init(&p1, 1 << 10);
            for (int i = 0; i < n1; i++) map_add(&p1, canon12(G1[i].v), G1[i].w);
            map p12 = combine(&p1, G2, n2); map_free(&p1);
            GROW_LIMIT = (size_t)1 << 25;
            map p123 = combine(&p12, G3, n3);
            map_free(&p12);
            if (pass == 0) fprintf(stderr, "canonical states: p123 %zu (%ld s)\n", p123.n, (long)(time(NULL) - T0));
            size_t lcap = ((size_t)1 << 24) / (size_t)PASSES; if (lcap < (1 << 16)) lcap = 1 << 16;
            map lowm; map_init(&lowm, lcap); GROW_LIMIT = 2 * lcap;
            size_t HCAP = pass == 0 ? (size_t)1 << 20 : 16; map highm; map_init(&highm, HCAP);
            size_t flush_at = 7 * (HCAP / 10), total_high = 0; int flushes = 0;
            ventry *buf = malloc(sizeof(ventry) * HCAP);          /* a flush can hold up to the table capacity */
            unsigned char p[12]; size_t done = 0;
            for (size_t i = 0; i < p123.cap; i++) {
                if (!p123.k[i]) continue;
                unpack12(p123.k[i] - 1, p); u128 wp = p123.w[i];
                for (int j = 0; j < n4; j++) {
                    int V[4];
                    for (int a = 0; a < 4; a++) { int s = 0; for (int b = 0; b < 3; b++) s += p[T[a][b]] * G4[j].v[T[a][b]]; V[a] = s; }
                    int sig = V[0] + V[1] + V[2] + V[3];
                    uint64_t key = canon4(V);
                    if (sig < SIG) { if ((int)(mix(key ^ 0x9e3779b97f4a7c15ULL) & (uint64_t)(PASSES - 1)) == pass) map_add(&lowm, key, wp * G4[j].w); }
                    else if (pass == 0) map_put(&highm, key, wp * G4[j].w);
                }
                if (pass == 0 && highm.n > flush_at) {
                    size_t n = 0; for (size_t k = 0; k < highm.cap; k++) if (highm.k[k]) { buf[n].key = highm.k[k] - 1; buf[n].w = highm.w[k]; n++; highm.k[k] = 0; }
                    highm.n = 0; total_high += n; flushes++;
                    run_batch(buf, n);
                }
                if (++done % 250000 == 0) fprintf(stderr, "pass %d stream: %zu of %zu p123 states, low table %zu, %d flushes (%ld s)\n", pass, done, p123.n, lowm.n, flushes, (long)(time(NULL) - T0));
            }
            if (pass == 0) {
                size_t n = 0; for (size_t k = 0; k < highm.cap; k++) if (highm.k[k]) { buf[n].key = highm.k[k] - 1; buf[n].w = highm.w[k]; n++; }
                total_high += n; flushes++; run_batch(buf, n);
            }
            free(buf); map_free(&highm); map_free(&p123);
            nlow = lowm.n; low = malloc(sizeof(ventry) * (nlow + 1));
            { size_t n = 0; for (size_t k = 0; k < lowm.cap; k++) if (lowm.k[k]) { low[n].key = lowm.k[k] - 1; low[n].w = lowm.w[k]; n++; } }
            map_free(&lowm);
            fprintf(stderr, "pass %d stream done: %zu high V entries in %d flushes, %zu low V vectors (%ld s)\n", pass, total_high, flushes, nlow, (long)(time(NULL) - T0));
            merge_hists();
            FILE *f = fopen(lowpath, "wb");
            if (!f || fwrite(&nlow, sizeof nlow, 1, f) != 1 || fwrite(low, sizeof(ventry), nlow, f) != nlow || fclose(f)) { fprintf(stderr, "cannot write low table\n"); return 6; }
            save_state("low", ((long long)pass << 32));
        }
        size_t CH = 32768;
        size_t nchunks = (nlow + CH - 1) / CH;
        time_t last_save = time(NULL);
        for (size_t c = (size_t)chunk0; c < nchunks; c++) {
            size_t a = c * CH, b = a + CH > nlow ? nlow : a + CH;
            run_batch(low + a, b - a);
            if (time(NULL) - last_save >= 300 && c + 1 < nchunks) { merge_hists(); save_state("low", ((long long)pass << 32) | (long long)(c + 1)); last_save = time(NULL); }
            if ((c + 1) % 16 == 0 || c + 1 == nchunks) fprintf(stderr, "pass %d low table: chunk %zu of %zu, %lld pairs so far (%ld s)\n", pass, c + 1, nchunks, PAIRS, (long)(time(NULL) - T0));
        }
        merge_hists();
        save_state("pass", ((long long)(pass + 1) << 32));
        free(low); have_low = 0; chunk0 = 0;
        if (pass == STOP_AFTER) { fprintf(stderr, "stopped after pass %d (test)\n", pass); return 9; }
    }
    for (long long nu = 0; nu <= CUT; nu++) {
        if (HIST[nu].w[0] | HIST[nu].w[1] | HIST[nu].w[2] | HIST[nu].w[3]) { printf("%lld ", nu); print_u256(stdout, &HIST[nu]); printf("\n"); }
    }
    printf("# below "); print_u256(stdout, &BELOW); printf("\n");
    printf("# above "); print_u256(stdout, &ABOVE); printf("\n");
    printf("# moment1 "); print_u256(stdout, &M1); printf("\n");
    printf("# moment2 "); print_u256(stdout, &M2); printf("\n");
    printf("# pairs %lld\n", PAIRS);
    printf("# partition_failures %lld\n", BAD);
    fprintf(stderr, "done: %lld pairs with nu <= C, %lld partition failures (%ld s)\n", PAIRS, BAD, (long)(time(NULL) - T0));
    return BAD ? 1 : 0;
}
