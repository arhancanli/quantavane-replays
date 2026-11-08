/* 9a: independence of a vertex set in the strong power C_q^{boxtimes k} (Quantavane replays; our own code).
 *
 * Two vertices x != y of C_q^{boxtimes k} are adjacent when every coordinate differs by 0 or +-1 mod q.  A set is
 * independent when no two of its (distinct) members are adjacent.  Two routes, sharing only the input parser:
 *   route 1: a bitmap of all q^k vertices; for every member, each of the 3^k - 1 neighbours is looked up;
 *   route 2: members grouped by their first h = k/2 coordinates (counting sort); for every member, the members of
 *            the 3^h groups whose prefix is adjacent or equal are compared coordinate by coordinate on the rest.
 * Both count adjacent ordered pairs; an independent set gives 0 in both.  Duplicates are refused.
 *
 * usage: strong_check Q K FILE [--inject-neighbour | --inject-duplicate | --wrong-graph]
 *   --inject-neighbour   appends member 0 with its first coordinate moved by +1 (a red control: must find edges)
 *   --inject-duplicate   appends a copy of member 0 (a red control: must be refused as a duplicate)
 *   --wrong-graph        counts pairs that differ by at most 2 in every coordinate (a red control: must find edges);
 *                        with it, only the first 500 members are used (5^k - 1 neighbours each in route 1)
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int Q, K, H, SPREAD = 1;
static uint64_t NV;            /* q^k */
static uint8_t *pts;           /* n x K */
static uint32_t *idx;          /* base-q index of each member */
static uint8_t *bits;          /* bitmap of members */
static size_t n;

static int near(int a, int b) { /* |a - b| mod q <= SPREAD (cyclic) */
    int d = a - b; if (d < 0) d = -d; if (Q - d < d) d = Q - d; return d <= SPREAD;
}

static uint64_t pw[24];

static uint64_t hits1;
static void walk(const uint8_t *x, int i, uint64_t acc, int moved) {
    if (i == K) {
        if (moved && (bits[acc >> 3] >> (acc & 7) & 1)) hits1++;
        return;
    }
    for (int d = -SPREAD; d <= SPREAD; d++) {
        int y = ((int)x[i] + d + Q) % Q;
        walk(x, i + 1, acc + (uint64_t)y * pw[i], moved || d != 0);
    }
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: strong_check Q K FILE [control]\n"); return 2; }
    Q = atoi(argv[1]); K = atoi(argv[2]); H = K / 2;
    const char *ctl = argc > 4 ? argv[4] : "";
    if (!strcmp(ctl, "--wrong-graph")) SPREAD = 2;
    pw[0] = 1;
    for (int i = 1; i <= K; i++) pw[i] = pw[i - 1] * (uint64_t)Q;
    NV = pw[K];
    if (NV > (1ULL << 32)) { fprintf(stderr, "q^k too large\n"); return 2; }
    FILE *f = fopen(argv[3], "r");
    if (!f) { perror("open"); return 2; }
    size_t cap = 1 << 16; n = 0;
    pts = malloc(cap * K + 2 * K);
    int col = 0;
    int c;
    long long line = 1;
    int innum = 0, cur = 0;
    /* strict parser: lines of exactly K integers in [0, Q), separated by single spaces */
    while ((c = fgetc(f)) != EOF) {
        if (c >= '0' && c <= '9') { cur = innum ? cur * 10 + (c - '0') : (c - '0'); innum = 1; if (cur >= Q) { fprintf(stderr, "line %lld: value out of range\n", line); return 1; } continue; }
        if (c == ' ' || c == '\n') {
            if (innum) { if (col >= K) { fprintf(stderr, "line %lld: too many values\n", line); return 1; } pts[n * K + col++] = (uint8_t)cur; innum = 0; }
            if (c == '\n') {
                if (col == 0) { fprintf(stderr, "line %lld: empty line\n", line); return 1; }
                if (col != K) { fprintf(stderr, "line %lld: %d values\n", line, col); return 1; }
                col = 0; n++; line++;
                if (n + 2 >= cap) { cap *= 2; pts = realloc(pts, cap * K + 2 * K); }
            }
            continue;
        }
        fprintf(stderr, "line %lld: unexpected byte %d\n", line, c); return 1;
    }
    if (innum || col) { fprintf(stderr, "last line not terminated\n"); return 1; }
    fclose(f);
    size_t n_file = n;
    if (SPREAD == 2 && n > 500) n = 500;
    if (!strcmp(ctl, "--inject-neighbour") || !strcmp(ctl, "--inject-duplicate")) {
        memcpy(pts + n * K, pts, K);
        if (!strcmp(ctl, "--inject-neighbour")) pts[n * K] = (uint8_t)((pts[n * K] + 1) % Q);
        n++;
    }
    idx = malloc(n * sizeof(uint32_t));
    bits = calloc(NV / 8 + 1, 1);
    size_t dup = 0;
    for (size_t a = 0; a < n; a++) {
        uint64_t s = 0;
        for (int i = 0; i < K; i++) s += pts[a * K + i] * pw[i];
        idx[a] = (uint32_t)s;
        if (bits[s >> 3] >> (s & 7) & 1) dup++;
        bits[s >> 3] |= (uint8_t)(1 << (s & 7));
    }
    printf("members read: %zu (file %zu), C_%d^%d, duplicates: %zu\n", n, n_file, Q, K, dup);
    if (dup) { printf("RESULT: REFUSED (duplicate members)\n"); return 1; }
    /* route 1 */
    hits1 = 0;
    for (size_t a = 0; a < n; a++) walk(pts + a * K, 0, 0, 0);
    unsigned long long nb = 1; for (int i = 0; i < K; i++) nb *= (unsigned long long)(2 * SPREAD + 1);
    printf("route 1 (bitmap, all %llu neighbours of each member): adjacent ordered pairs %llu\n",
           nb - 1, (unsigned long long)hits1);
    /* route 2: counting sort by prefix (first H coordinates) */
    uint64_t G = pw[H];
    uint32_t *start = calloc(G + 1, sizeof(uint32_t));
    uint32_t *order = malloc(n * sizeof(uint32_t));
    for (size_t a = 0; a < n; a++) { uint64_t p = idx[a] % G; start[p + 1]++; }
    for (uint64_t g = 0; g < G; g++) start[g + 1] += start[g];
    uint32_t *fill = malloc((G + 1) * sizeof(uint32_t));
    memcpy(fill, start, (G + 1) * sizeof(uint32_t));
    for (size_t a = 0; a < n; a++) order[fill[idx[a] % G]++] = (uint32_t)a;
    uint64_t hits2 = 0;
    int span = 2 * SPREAD + 1;
    int combos = 1; for (int i = 0; i < H; i++) combos *= span;
    for (size_t a = 0; a < n; a++) {
        const uint8_t *x = pts + a * K;
        for (int cmb = 0; cmb < combos; cmb++) {
            int t = cmb; uint64_t key = 0;
            for (int i = 0; i < H; i++) { int d = t % span - SPREAD; t /= span; key += (uint64_t)(((int)x[i] + d + Q) % Q) * pw[i]; }
            for (uint32_t r = start[key]; r < start[key + 1]; r++) {
                uint32_t b = order[r];
                if (b == a) continue;
                const uint8_t *y = pts + (size_t)b * K;
                int ok = 1;
                for (int i = H; i < K && ok; i++) ok = near(x[i], y[i]);
                if (ok) hits2++;
            }
        }
    }
    printf("route 2 (grouped by the first %d coordinates, %d prefix groups per member): adjacent ordered pairs %llu\n",
           H, combos, (unsigned long long)hits2);
    int agree = hits1 == hits2;
    printf("routes agree: %s\n", agree ? "yes" : "NO");
    if (agree && hits1 == 0) printf("RESULT: INDEPENDENT (%zu members, no adjacent pair)\n", n);
    else printf("RESULT: NOT INDEPENDENT (%llu adjacent ordered pairs)\n", (unsigned long long)hits1);
    return (agree && hits1 == 0) ? 0 : 1;
}
