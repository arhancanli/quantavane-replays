/* c12a_check.c: an exact replay of the computation behind beta <= 0.8618 (12a, the Euclidean travelling salesman
 * constant; teorth/optimizationproblems PR #222, D. O'Keefe), written from RECIPE.md of
 * devinokeefe/travelling-sales-person-constant-bounds @ 69b864b9 alone: the repository's checkers were neither read
 * nor run.
 *
 *   c12a_check SPEC [-j N] [--orders] [--hk G1,G2,...]
 *
 * SPEC is a spec file (RECIPE.md, Section 10) and the values file is the one it names. Every integer of Sections 2
 * to 8 is recomputed from the parameters and compared with the values file: the masses, the tail term, each table
 * line (key and entries), each cell sum and the total. Then the test of Section 8, A_0 h_2 theta_2 <= theta_1
 * n^(k+1) F k h_1. Exit status 0 if every number agrees and the test holds, 1 if a number differs or the test fails,
 * 2 if a file is malformed or a value would leave the widths used.
 *
 * Widths. Sections 2 to 6 and 8 use 128-bit integers; every multiplication, addition and subtraction there is
 * checked, and any overflow stops the program with status 2. The argument of ceil_sqrt is required to be below
 * 2^126 and the weights, Dist, G and Tab below 2^64. Section 7 uses 64-bit costs after checking that
 * k max(Tab) < 2^60 (a path has at most k edges; 2^61 marks an impossible path) and n^3 k max(Tab) < 2^64 (the
 * largest sum formed in 64 bits), and adds into 128-bit cell sums.
 *
 * Cell sums (Section 7), for k >= 3, with L = k - 1 interior points. cost(g, a) = min over interior u != v of
 * E_0u + H(u, v) + E_vk, where H(u, v) is the least cost of a path from u to v through all interior points (the
 * shared-paths identity); H depends on the middle (g_2..g_L; a_1..a_L) only. The heights a_1, ..., a_L are fixed in
 * nested loops. In the loop of a_j the least costs f(u, S, w) of paths from u to w through S, for the sets S whose
 * largest point is j, are formed from those of smaller sets (Held-Karp), and all inner loops reuse them. In the loop
 * of a_(L-1) the two-path covers of {1..L-1} are formed: K(u, p, q, v) = min over S of f(u, S, p) + f(q, S', v), S'
 * the complement. Inside, H(u, v) = min over p, q of K(u, p, q, v) + E_pL + E_Lq for u, v < L (point L between p
 * and q), and H(u, L) = min over p of f(u, {1..L-1}, p) + E_pL. Then z_v = min over u of E_0u + H(u, v) for each
 * (g_1, a_0), and cost = min over v of z_v + E_vk for each (g_k, a_k).
 * Reversal: S_g' = S_g for the reversed tuple g', so a middle gap tuple and its reverse are done once; for a
 * palindromic middle both S_(g1, mid, gk) and S_(gk, mid, g1) are computed and must agree.
 * Reflection, on a_1 only: a_1 = y and a_1 = n - 1 - y contribute equally, so a_1 < (n - 1)/2 is counted twice.
 *
 * --orders: every cell sum is instead the literal sum over all n^(k+1) height tuples of the least cost over all
 *           (k - 1)! orders, with none of the three identities.
 * --hk:     only the listed cell sums (indices g_1 + g_2 m + ...) are recomputed, each cell by Held-Karp from point
 *           0 (no identities), and compared; no total. For spot checks where all orders are too many.
 */
#include <math.h>
#include <pthread.h>
#include <stdarg.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef unsigned __int128 u128;
typedef uint64_t u64;

#define MAXM 8
#define MAXK 9
#define MAXN 16
#define QQ ((u128)1 << 62)
#define DD ((u128)1 << 48)
#define FF ((u128)1 << 32)
#define TWO64 ((u128)1 << 64)

static void die(int status, const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    fprintf(stderr, "c12a_check: ");
    vfprintf(stderr, fmt, ap);
    fprintf(stderr, "\n");
    va_end(ap);
    fflush(stdout);
    exit(status);
}

static u128 add(u128 a, u128 b) { u128 r; if (__builtin_add_overflow(a, b, &r)) die(2, "overflow (addition)"); return r; }
static u128 mul(u128 a, u128 b) { u128 r; if (__builtin_mul_overflow(a, b, &r)) die(2, "overflow (multiplication)"); return r; }
static u128 sub(u128 a, u128 b) { if (a < b) die(2, "negative difference"); return a - b; }
static u128 fdiv(u128 a, u128 b) { if (b == 0) die(2, "division by zero"); return a / b; }
static u128 cdiv(u128 a, u128 b) { if (b == 0) die(2, "division by zero"); return a / b + (a % b != 0); }
static u64 to64(u128 a, const char *what) { if (a >= TWO64) die(2, "%s is 2^64 or more", what); return (u64)a; }
static int iabs(int x) { return x < 0 ? -x : x; }

static const char *u128s(u128 v, char *buf)
{
    char t[48];
    int n = 0;
    do { t[n++] = (char)('0' + (int)(v % 10)); v /= 10; } while (v);
    for (int i = 0; i < n; i++) buf[i] = t[n - 1 - i];
    buf[n] = 0;
    return buf;
}

static double now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

/* ---------- threads ---------- */
static int nthreads = 1;
typedef struct { void (*fn)(long, int, void *); void *ctx; long count; atomic_long next; } Job;
typedef struct { Job *job; int tid; } Arg;

static void *worker(void *p)
{
    Arg *a = p;
    for (;;) {
        long i = atomic_fetch_add(&a->job->next, 1);
        if (i >= a->job->count) break;
        a->job->fn(i, a->tid, a->job->ctx);
    }
    return NULL;
}

static void parallel_for(long count, void (*fn)(long, int, void *), void *ctx)
{
    Job job;
    job.fn = fn; job.ctx = ctx; job.count = count;
    atomic_init(&job.next, 0);
    int nt = nthreads < count ? nthreads : (int)(count > 0 ? count : 1);
    pthread_t th[64];
    Arg args[64];
    for (int t = 0; t < nt; t++) {
        args[t].job = &job; args[t].tid = t;
        if (t > 0 && pthread_create(&th[t], NULL, worker, &args[t])) die(2, "cannot start a thread");
    }
    worker(&args[0]);
    for (int t = 1; t < nt; t++) pthread_join(th[t], NULL);
}

/* ---------- files (Section 10) ---------- */
typedef struct { int nf; char **f; } Line;
typedef struct { int n; Line *l; const char *path; } Text;

static Text read_text(const char *path)
{
    FILE *fp = fopen(path, "rb");
    if (!fp) die(2, "cannot open %s", path);
    size_t cap = 1 << 16, len = 0;
    char *buf = malloc(cap + 1);
    for (;;) {
        len += fread(buf + len, 1, cap - len, fp);
        if (len < cap) { if (ferror(fp)) die(2, "cannot read %s", path); break; }
        cap *= 2;
        buf = realloc(buf, cap + 1);
    }
    fclose(fp);
    buf[len] = 0;
    for (size_t i = 0; i < len; i++) {
        unsigned char c = (unsigned char)buf[i];
        if (c == 10 || c == 9 || (c >= 32 && c <= 126)) continue;
        if (c == 13 && i + 1 < len && buf[i + 1] == 10) continue;
        die(2, "%s: byte %u at offset %zu is not allowed", path, c, i);
    }
    Text t = {0, NULL, path};
    int capl = 64;
    t.l = malloc(capl * sizeof(Line));
    size_t s = 0;
    while (s < len) {
        size_t e = s;
        while (e < len && buf[e] != 10) e++;
        size_t end = e;
        if (end > s && buf[end - 1] == 13) end--;
        buf[end] = 0;
        Line ln = {0, NULL};
        int capf = 4;
        ln.f = malloc(capf * sizeof(char *));
        size_t p = s;
        while (p < end) {
            while (p < end && (buf[p] == ' ' || buf[p] == '\t')) buf[p++] = 0;
            if (p >= end) break;
            if (ln.nf == capf) { capf *= 2; ln.f = realloc(ln.f, capf * sizeof(char *)); }
            ln.f[ln.nf++] = buf + p;
            while (p < end && buf[p] != ' ' && buf[p] != '\t') p++;
        }
        if (ln.nf == 0) die(2, "%s: line %d has no fields", path, t.n + 1);
        if (t.n == capl) { capl *= 2; t.l = realloc(t.l, capl * sizeof(Line)); }
        t.l[t.n++] = ln;
        s = e + 1;
    }
    if (t.n == 0) die(2, "%s: no lines", path);
    return t;
}

static u128 tonum(const char *s, const char *path, int line)
{
    if (!*s) die(2, "%s: line %d: empty number", path, line);
    u128 v = 0;
    for (const char *c = s; *c; c++) {
        if (*c < '0' || *c > '9') die(2, "%s: line %d: '%s' is not a number", path, line, s);
        u128 t;
        if (__builtin_mul_overflow(v, (u128)10, &t) || __builtin_add_overflow(t, (u128)(*c - '0'), &v))
            die(2, "%s: line %d: a number of 2^128 or more", path, line);
    }
    return v;
}

/* the seven parameter lines shared by spec and values files */
static const char *pname[7] = {"k", "h", "mesh", "gaps", "heights", "height_mesh", "tail"};
static const int pcnt[7] = {1, 2, 1, -1, 1, 1, 1};

/* ---------- parameters (Section 1) ---------- */
static int k, m, n, nu, T;
static u64 h1, h2, M, th1, th2, NT;
static u64 b[MAXM + 1], ncell[MAXM];
static int ngaps;
static u128 specnum[8][8];
static int speccnt[8];

static u64 gcd64(u64 a, u64 c) { while (c) { u64 t = a % c; a = c; c = t; } return a; }

static char *read_spec(const char *path)
{
    Text t = read_text(path);
    if (t.n != 9) die(2, "%s: %d lines, not 9", path, t.n);
    static const char *nm[8] = {"k", "h", "mesh", "gaps", "heights", "height_mesh", "tail", "threshold"};
    static const int cnt[8] = {1, 2, 1, -1, 1, 1, 1, 2};
    for (int i = 0; i < 8; i++) {
        Line *ln = &t.l[i];
        if (strcmp(ln->f[0], nm[i])) die(2, "%s: line %d is not '%s'", path, i + 1, nm[i]);
        if (cnt[i] >= 0 && ln->nf - 1 != cnt[i]) die(2, "%s: line %d needs %d numbers", path, i + 1, cnt[i]);
        if (ln->nf - 1 > 7) die(2, "%s: more than 7 gap boundaries", path);
        speccnt[i] = ln->nf - 1;
        for (int j = 1; j < ln->nf; j++) specnum[i][j - 1] = tonum(ln->f[j], path, i + 1);
    }
    Line *vl = &t.l[8];
    if (strcmp(vl->f[0], "values") || vl->nf != 2) die(2, "%s: line 9 is not 'values PATH'", path);
    const char *vp = vl->f[1];
    if (vp[0] == '/' || strchr(vp, '\\') || strchr(vp, ':')) die(2, "%s: values path not allowed", path);
    /* ranges */
    if (specnum[0][0] < 2 || specnum[0][0] > 9) die(2, "k out of range");
    k = (int)specnum[0][0];
    if (specnum[1][0] == 0 || specnum[1][1] == 0 || specnum[1][0] >= TWO64 || specnum[1][1] >= TWO64)
        die(2, "h out of range");
    h1 = (u64)specnum[1][0]; h2 = (u64)specnum[1][1];
    if (gcd64(h1, h2) != 1) die(2, "h is not in lowest terms");
    if (specnum[2][0] < 2 || specnum[2][0] > 1024) die(2, "mesh out of range");
    M = (u64)specnum[2][0];
    ngaps = speccnt[3];
    m = ngaps + 1;
    b[0] = 0;
    for (int i = 0; i < ngaps; i++) {
        if (specnum[3][i] >= TWO64) die(2, "gap boundary out of range");
        b[i + 1] = (u64)specnum[3][i];
        if (b[i + 1] <= b[i]) die(2, "gap boundaries not increasing from 1");
    }
    if (specnum[4][0] < 1 || specnum[4][0] > 16) die(2, "heights out of range");
    n = (int)specnum[4][0];
    if (specnum[5][0] < 1 || specnum[5][0] > 256) die(2, "height_mesh out of range");
    nu = (int)specnum[5][0];
    if (specnum[6][0] < 1 || specnum[6][0] > 40) die(2, "tail out of range");
    T = (int)specnum[6][0];
    if (specnum[7][0] == 0 || specnum[7][1] == 0 || specnum[7][0] >= TWO64 || specnum[7][1] >= TWO64)
        die(2, "threshold out of range");
    th1 = (u64)specnum[7][0]; th2 = (u64)specnum[7][1];
    NT = (u64)T * M;
    for (int c = 0; c + 1 < m; c++) ncell[c] = b[c + 1] - b[c];
    ncell[m - 1] = NT;
    /* the values path is relative to the folder of the spec path as given */
    const char *slash = strrchr(path, '/');
    size_t dl = slash ? (size_t)(slash - path + 1) : 0;
    char *full = malloc(dl + strlen(vp) + 1);
    memcpy(full, path, dl);
    strcpy(full + dl, vp);
    return full;
}

/* ---------- Section 2 ---------- */
static u128 lo1, up1, *Rup, *Rlo;
static u64 Imax;

static void section2(void)
{
    u128 lo[14], up[14];
    lo[13] = 0;
    up[13] = QQ;
    for (int j = 12; j >= 1; j--) {
        u128 dv = mul(M, (u128)j);
        up[j] = sub(QQ, fdiv(lo[j + 1], dv));
        lo[j] = sub(QQ, cdiv(up[j + 1], dv));
    }
    lo1 = lo[1];
    up1 = up[1];
    Imax = b[m - 1] > NT ? b[m - 1] : NT;
    if (Imax > ((u64)1 << 22)) die(2, "a boundary or T*M above 2^22 needs more memory than this program provides");
    Rup = malloc((Imax + 1) * sizeof(u128));
    Rlo = malloc((Imax + 1) * sizeof(u128));
    Rup[0] = Rlo[0] = QQ;
    for (u64 i = 1; i <= Imax; i++) {
        Rup[i] = cdiv(mul(Rup[i - 1], up1), QQ);
        Rlo[i] = fdiv(mul(Rlo[i - 1], lo1), QQ);
    }
}

/* ---------- Section 3 ---------- */
static u128 Clam, Cmu, P[MAXM], ET;
static u64 *W[MAXM];

static void section3(void)
{
    Clam = sub(QQ, mul(M, sub(QQ, up1)));
    Cmu = sub(mul(M, sub(QQ, lo1)), lo1);
    if (Clam == 0 || Cmu == 0) die(2, "C_lambda or C_mu is not positive");
    for (int c = 0; c < m; c++) {
        u64 nc = ncell[c];
        W[c] = malloc((nc + 1) * sizeof(u64));
        for (u64 i = 0; i <= nc; i++) {
            u128 num = 0, w;
            if (i < nc) num = add(num, mul(Rup[i], Clam));
            if (i >= 1) num = add(num, mul(Rup[i - 1], Cmu));
            if (c < m - 1) {
                w = cdiv(cdiv(num, (u128)1 << 14), sub(QQ, Rup[nc]));
            } else {
                w = cdiv(num, (u128)1 << 76);
                if (i == NT) w = add(w, cdiv(Rup[NT], (u128)1 << 14));
            }
            W[c][i] = to64(w, "a weight W");
        }
    }
    for (int c = 0; c + 1 < m; c++) P[c] = cdiv(sub(Rup[b[c]], Rlo[b[c + 1]]), (u128)1 << 14);
    P[m - 1] = cdiv(Rup[b[m - 1]], (u128)1 << 14);
    ET = cdiv(Rup[NT], (u128)1 << 30);
}

/* ---------- Section 4 ---------- */
static int NK;
static unsigned char (*key)[MAXM];
static int *kpar, *kcell, *ksz;
static u64 *Voff, *Vlen, *Vmax, *Vbuf, *bkey;

static u128 binom(int a, int c)
{
    u128 r = 1;
    for (int i = 1; i <= c; i++) r = r * (u128)(a - c + i) / (u128)i;
    return r;
}

/* index of a key in key order: the vectors with sum <= k before it, less the zero vector */
static int key_rank(const int *kv)
{
    u128 r = 0;
    int rem = k;
    for (int c = 0; c < m; c++) {
        for (int x = 0; x < kv[c]; x++) r += binom(m - c - 1 + rem - x, rem - x);
        rem -= kv[c];
    }
    if (rem < 0) die(2, "internal: key with more than k gaps");
    return (int)r - 1;
}

static void gen_keys(int c, int rem, int *cur, int *cnt)
{
    if (c == m) {
        int s = 0;
        for (int i = 0; i < m; i++) s += cur[i];
        if (s == 0) return;
        for (int i = 0; i < m; i++) key[*cnt][i] = (unsigned char)cur[i];
        (*cnt)++;
        return;
    }
    for (int x = 0; x <= rem; x++) { cur[c] = x; gen_keys(c + 1, rem - x, cur, cnt); }
    cur[c] = 0;
}

static int *vlist, nvlist;

static void vjob(long idx, int tid, void *ctx)
{
    (void)tid; (void)ctx;
    int t = vlist[idx], c = kcell[t], tp = kpar[t];
    u64 nc = ncell[c], len = Vlen[t], lenp = Vlen[tp], mw = 0;
    const u64 *Vp = Vbuf + Voff[tp], *Wc = W[c];
    u64 *V = Vbuf + Voff[t];
    for (u64 j = 0; j <= nc; j++) if (Wc[j] > mw) mw = Wc[j];
    /* every sum below is at most (n_c + 1) max(V') max(W): check that it stays below 2^128 */
    mul(mul(Vmax[tp], mw), nc + 1);
    u64 mx = 0;
    for (u64 s = 0; s < len; s++) {
        u64 jlo = s + 1 > lenp ? s + 1 - lenp : 0, jhi = s < nc ? s : nc;
        u128 acc = 0;
        for (u64 j = jlo; j <= jhi; j++) acc += (u128)Vp[s - j] * Wc[j];
        u64 v = to64(cdiv(acc, DD), "a weight V");
        V[s] = v;
        if (v > mx) mx = v;
    }
    Vmax[t] = mx;
}

static void section4(void)
{
    NK = (int)(binom(m + k, k) - 1);
    key = malloc((size_t)NK * sizeof *key);
    int cur[MAXM] = {0}, cnt = 0;
    gen_keys(0, k, cur, &cnt);
    if (cnt != NK) die(2, "internal: %d keys, expected %d", cnt, NK);
    kpar = malloc(NK * sizeof(int)); kcell = malloc(NK * sizeof(int)); ksz = malloc(NK * sizeof(int));
    Voff = malloc(NK * sizeof(u64)); Vlen = malloc(NK * sizeof(u64)); Vmax = malloc(NK * sizeof(u64));
    bkey = malloc(NK * sizeof(u64));
    u64 tot = 0;
    for (int t = 0; t < NK; t++) {
        int kv[MAXM], s = 0, c = -1;
        u128 len = 1, bk = 0;
        for (int i = 0; i < m; i++) {
            kv[i] = key[t][i];
            s += kv[i];
            if (kv[i] > 0) c = i;
            len = add(len, mul(kv[i], ncell[i]));
            bk = add(bk, mul(kv[i], b[i]));
        }
        if (key_rank(kv) != t) die(2, "internal: key rank");
        ksz[t] = s;
        kcell[t] = c;
        kv[c]--;
        kpar[t] = s == 1 ? -1 : key_rank(kv);
        Vlen[t] = (u64)len;
        bkey[t] = to64(bk, "b_kappa");
        Voff[t] = tot;
        tot += (u64)len;
        if (tot > ((u64)1 << 26)) die(2, "the weight vectors need more memory than this program provides");
    }
    Vbuf = malloc(tot * sizeof(u64));
    for (int t = 0; t < NK; t++)
        if (ksz[t] == 1) {
            int c = kcell[t];
            u64 mx = 0;
            for (u64 i = 0; i <= ncell[c]; i++) { Vbuf[Voff[t] + i] = W[c][i]; if (W[c][i] > mx) mx = W[c][i]; }
            Vmax[t] = mx;
        }
    vlist = malloc(NK * sizeof(int));
    for (int s = 2; s <= k; s++) {
        nvlist = 0;
        for (int t = NK - 1; t >= 0; t--) if (ksz[t] == s) vlist[nvlist++] = t;  /* longer vectors first */
        parallel_for(nvlist, vjob, NULL);
    }
}

/* ---------- Section 5 ---------- */
static u64 *G, J;
static u128 alpha, gam, aM;
static u64 chi[2 * 256 + 1];

static u64 ceil_sqrt(u128 A)    /* A < 2^126 */
{
    if (A == 0) return 0;
    u64 s = (u64)sqrt((double)A);
    if (s == 0) s = 1;
    for (int it = 0; it < 2; it++) s = (u64)(((u128)s + A / s) / 2);
    while ((u128)s * s > A) s--;
    while ((u128)(s + 1) * (s + 1) <= A) s++;
    return (u128)s * s == A ? s : s + 1;
}

static void gjob(long blk, int tid, void *ctx)
{
    (void)tid; (void)ctx;
    int ntau = n * nu + 1;
    u64 *dist = malloc((size_t)ntau * sizeof(u64));
    u64 j0 = (u64)blk * 64, j1 = j0 + 64 > J + 1 ? J + 1 : j0 + 64;
    for (u64 j = j0; j < j1; j++) {
        u128 aj2 = mul(mul(alpha, j), mul(alpha, j));
        for (int tau = 0; tau < ntau; tau++) {
            u128 X = add(aj2, mul(mul(gam, (u128)tau), mul(gam, (u128)tau)));
            if (X >= ((u128)1 << 62)) die(2, "the argument of ceil_sqrt is 2^126 or more");
            dist[tau] = to64(cdiv(ceil_sqrt(X << 64), aM), "Dist");
        }
        for (int d = 0; d < n; d++) {
            u128 s = 0;
            for (int t = -nu; t <= nu; t++) s = add(s, mul(chi[t + nu], dist[iabs(d * nu + t)]));
            G[j * (u64)n + (u64)d] = to64(cdiv(s, (u128)4 * nu * nu), "G");
        }
    }
    free(dist);
}

static void section5(void)
{
    for (int t = -nu; t <= nu; t++) {
        u64 s = 0;
        for (int i = 0; i <= nu; i++) {
            int j = i - t;
            if (j < 0 || j > nu) continue;
            u64 ei = (i == 0 || i == nu) ? 1 : 2, ej = (j == 0 || j == nu) ? 1 : 2;
            s += ei * ej;
        }
        chi[t + nu] = s;
    }
    alpha = mul(mul(mul(h2, h2), n), nu);
    gam = mul(M, mul(h1, h1));
    aM = mul(alpha, M);
    u128 Jw = mul(k, add(b[m - 1], NT));
    if (mul(add(Jw, 1), n) > ((u128)1 << 25)) die(2, "the table of G needs more memory than this program provides");
    J = (u64)Jw;
    G = malloc((J + 1) * (u64)n * sizeof(u64));
    parallel_for((long)(J / 64 + 1), gjob, NULL);
}

/* ---------- Section 6 ---------- */
static u64 *Tab, maxTab;

static void tjob(long t, int tid, void *ctx)
{
    (void)tid; (void)ctx;
    u128 acc[MAXN] = {0};
    const u64 *V = Vbuf + Voff[t];
    for (u64 s = 0; s < Vlen[t]; s++) {
        const u64 *g = G + (bkey[t] + s) * (u64)n;
        for (int d = 0; d < n; d++) acc[d] = add(acc[d], mul(V[s], g[d]));
    }
    for (int d = 0; d < n; d++)
        Tab[t * n + d] = to64(add(cdiv(acc[d], DD), mul(key[t][m - 1], ET)), "Tab");
}

static void section6(void)
{
    Tab = malloc((size_t)NK * n * sizeof(u64));
    parallel_for(NK, tjob, NULL);
    maxTab = 0;
    for (long i = 0; i < (long)NK * n; i++) if (Tab[i] > maxTab) maxTab = Tab[i];
}

/* ---------- Section 7 ---------- */
#define CINF ((u64)1 << 61)
#define NV 8                 /* lanes for the interior points 1..L (v - 1), padded with CINF */
static long mpw[MAXK + 2];
static u128 *Ssum;
static unsigned char *Sset;
static int L, L1;
static int npair, pid[MAXK][MAXK], nmem[128], mem[128][8];
static atomic_long done7;
static long todo7;
static double t7start;
static pthread_mutex_t plock = PTHREAD_MUTEX_INITIALIZER;

typedef struct {
    const u64 *er[MAXK][MAXK], *x0[MAXM][MAXK], *yk[MAXM][MAXK];
    int a[MAXK + 1];
    u64 weight;
    u64 E[MAXK][MAXK], eL[MAXK];
    u64 f[128][MAXK][MAXK];
    u64 Kc[28][28];
    u64 KK[21][NV][NV];
    u64 Fall[MAXK][MAXK], H[MAXK][NV];
    u64 z[MAXM][MAXN][NV], Y[MAXM][NV][MAXN];
    u64 acc[MAXM][MAXM];
    u128 S16[MAXM][MAXM];
} Work;
static Work *works[64];

static const u64 *row_of(const int *kv) { return Tab + (size_t)key_rank(kv) * n; }

/* f(u, S, x) for the sets S with largest point j, by removing the end x; f(x, S, u) = f(u, S, x) */
static void hk_level(Work *w, int j)
{
    int top = 1 << (j - 1);
    w->f[top][j][j] = 0;
    for (int S = top + 1; S < 2 * top; S++)
        for (int ix = 0; ix < nmem[S]; ix++) {
            int x = mem[S][ix], Pm = S & ~(1 << (x - 1));
            for (int iu = 0; iu < ix; iu++) {
                int u = mem[S][iu];
                u64 best = CINF;
                for (int iy = 0; iy < nmem[Pm]; iy++) {
                    int y = mem[Pm][iy];
                    u64 c = w->f[Pm][u][y] + w->E[y][x];
                    best = c < best ? c : best;
                }
                w->f[S][u][x] = w->f[S][x][u] = best;
            }
        }
}

static void covers(Work *w)
{
    int full = (1 << L1) - 1;
    for (int i = 0; i < npair; i++) for (int i2 = 0; i2 < npair; i2++) w->Kc[i][i2] = CINF;
    for (int S = 1; S < full; S += 2) {
        int C = full ^ S;
        for (int iu = 0; iu < nmem[S]; iu++)
            for (int ip = iu; ip < nmem[S]; ip++) {
                if (ip == iu && nmem[S] > 1) continue;
                int u = mem[S][iu], p = mem[S][ip], i1 = pid[u][p];
                u64 fa = w->f[S][u][p];
                for (int iq = 0; iq < nmem[C]; iq++)
                    for (int iv = iq; iv < nmem[C]; iv++) {
                        if (iv == iq && nmem[C] > 1) continue;
                        int q = mem[C][iq], v = mem[C][iv], i2 = pid[q][v];
                        u64 c = fa + w->f[C][q][v];
                        if (c < w->Kc[i1][i2]) w->Kc[i1][i2] = w->Kc[i2][i1] = c;
                    }
            }
    }
    for (int u = 1; u <= L1; u++) for (int p = 1; p <= L1; p++) w->Fall[u][p] = w->f[full][u][p];
    /* KK[u < v][p - 1][q - 1] = K(u, p, q, v), CINF where p = v, q = u or q = p */
    for (int u = 1, uv = 0; u <= L1; u++)
        for (int v = u + 1; v <= L1; v++, uv++)
            for (int p = 1; p <= L1; p++)
                for (int q = 1; q <= L1; q++)
                    w->KK[uv][p - 1][q - 1] = (p == v || q == u || q == p) ? CINF : w->Kc[pid[u][p]][pid[q][v]];
}

/* H(u, v) for the current a_L, then the endpoint loops; NN is the number of heights */
static inline __attribute__((always_inline)) void inner_n(Work *w, const int NN)
{
    u64 el[NV], wpq[NV][NV];
    for (int p = 0; p < NV; p++) el[p] = p < L1 ? w->E[p + 1][L] : 0;
    for (int p = 0; p < NV; p++) for (int q = 0; q < NV; q++) wpq[p][q] = el[p] + el[q];
    for (int u = 1, uv = 0; u <= L1; u++)
        for (int v = u + 1; v <= L1; v++, uv++) {
            u64 bv[NV];
            for (int q = 0; q < NV; q++) bv[q] = CINF;
            for (int p = 0; p < L1; p++)
                for (int q = 0; q < NV; q++) { u64 c = w->KK[uv][p][q] + wpq[p][q]; bv[q] = c < bv[q] ? c : bv[q]; }
            u64 best = CINF;
            for (int q = 0; q < NV; q++) best = bv[q] < best ? bv[q] : best;
            w->H[u][v - 1] = w->H[v][u - 1] = best;
        }
    for (int p = 1; p <= L1; p++) w->eL[p] = w->E[p][L];
    for (int u = 1; u <= L1; u++) {
        u64 best = CINF;
        for (int p = 1; p <= L1; p++) { u64 c = w->Fall[u][p] + w->eL[p]; best = c < best ? c : best; }
        w->H[u][L - 1] = w->H[L][u - 1] = best;
    }
    for (int g1 = 0; g1 < m; g1++)
        for (int a0 = 0; a0 < NN; a0++) {
            u64 zr[NV];
            for (int v = 0; v < NV; v++) zr[v] = CINF;
            for (int u = 1; u <= L; u++) {
                u64 xu = w->x0[g1][u][iabs(a0 - w->a[u])];
                const u64 *hr = w->H[u];
                for (int v = 0; v < NV; v++) { u64 c = xu + hr[v]; zr[v] = c < zr[v] ? c : zr[v]; }
            }
            for (int v = 0; v < NV; v++) w->z[g1][a0][v] = zr[v];
        }
    for (int gk = 0; gk < m; gk++)
        for (int v = 1; v <= L; v++)
            for (int ak = 0; ak < NN; ak++) w->Y[gk][v - 1][ak] = w->yk[gk][v][iabs(w->a[v] - ak)];
    for (int g1 = 0; g1 < m; g1++)
        for (int a0 = 0; a0 < NN; a0++)
            for (int gk = 0; gk < m; gk++) {
                u64 best[MAXN];
                for (int ak = 0; ak < NN; ak++) best[ak] = CINF;
                for (int v = 0; v < L; v++) {
                    u64 zv = w->z[g1][a0][v];
                    const u64 *yr = w->Y[gk][v];
                    for (int ak = 0; ak < NN; ak++) { u64 c = zv + yr[ak]; best[ak] = c < best[ak] ? c : best[ak]; }
                }
                u64 sum = 0;
                for (int ak = 0; ak < NN; ak++) sum += best[ak];
                w->acc[g1][gk] += sum;
            }
}

static void inner6(Work *w) { inner_n(w, 6); }
static void inner_any(Work *w) { inner_n(w, n); }

static void lev(Work *w, int j)
{
    int hi = j == 1 ? (n + 1) / 2 : n;
    if (j == L) memset(w->acc, 0, sizeof w->acc);
    for (int y = 0; y < hi; y++) {
        if (j == 1) w->weight = 2 * y < n - 1 ? 2 : 1;
        w->a[j] = y;
        for (int i = 1; i < j; i++) w->E[i][j] = w->E[j][i] = w->er[i][j][iabs(w->a[i] - y)];
        if (j < L) {
            hk_level(w, j);
            if (j == L1) covers(w);
            lev(w, j + 1);
        } else {
            if (n == 6) inner6(w); else inner_any(w);
        }
    }
    if (j == L)
        for (int g1 = 0; g1 < m; g1++)
            for (int gk = 0; gk < m; gk++) w->S16[g1][gk] += (u128)w->weight * w->acc[g1][gk];
}

static void progress(void)
{
    long d = atomic_fetch_add(&done7, 1) + 1;
    long step = todo7 / 50 > 0 ? todo7 / 50 : 1;
    if (todo7 >= 100 && (d % step == 0 || d == todo7)) {
        double el = now() - t7start;
        pthread_mutex_lock(&plock);
        fprintf(stderr, "section 7: %ld of %ld done, %.0f s, about %.0f s left\n", d, todo7, el,
                el / (double)d * (double)(todo7 - d));
        pthread_mutex_unlock(&plock);
    }
}

static void task7(long ti, int tid, void *ctx)
{
    (void)ctx;
    Work *w = works[tid];
    int g[MAXK + 2];
    long r = ti;
    for (int l = 2; l <= L; l++) { g[l] = (int)(r % m); r /= m; }
    int cmp = 0;
    for (int l = 2; l <= L && !cmp; l++) if (g[l] != g[L + 2 - l]) cmp = g[l] < g[L + 2 - l] ? -1 : 1;
    if (cmp > 0) return;
    int kv[MAXM];
    for (int i = 1; i <= L; i++)
        for (int j = i + 1; j <= L; j++) {
            memset(kv, 0, sizeof kv);
            for (int l = i + 1; l <= j; l++) kv[g[l]]++;
            w->er[i][j] = row_of(kv);
        }
    for (int c = 0; c < m; c++)
        for (int u = 1; u <= L; u++) {
            memset(kv, 0, sizeof kv);
            kv[c]++;
            for (int l = 2; l <= u; l++) kv[g[l]]++;
            w->x0[c][u] = row_of(kv);
            memset(kv, 0, sizeof kv);
            kv[c]++;
            for (int l = u + 1; l <= L; l++) kv[g[l]]++;
            w->yk[c][u] = row_of(kv);
        }
    memset(w->S16, 0, sizeof w->S16);
    lev(w, 1);
    long mid = 0, rmid = 0;
    for (int l = 2; l <= L; l++) { mid += g[l] * mpw[l - 1]; rmid += g[L + 2 - l] * mpw[l - 1]; }
    for (int g1 = 0; g1 < m; g1++)
        for (int gk = 0; gk < m; gk++) {
            long i1 = g1 + mid + gk * mpw[k - 1], i2 = gk + rmid + g1 * mpw[k - 1];
            Ssum[i1] = w->S16[g1][gk];
            Sset[i1]++;
            if (cmp < 0) { Ssum[i2] = w->S16[g1][gk]; Sset[i2]++; }
            else if (w->S16[g1][gk] != w->S16[gk][g1]) die(2, "internal: reversal fails for a palindromic middle");
        }
    progress();
}

static void setup7(void)
{
    L = k - 1;
    L1 = L - 1;
    npair = 0;
    for (int u = 1; u <= L1; u++) for (int p = u; p <= L1; p++) pid[u][p] = pid[p][u] = npair++;
    for (int S = 0; S < (1 << L1); S++) {
        nmem[S] = 0;
        for (int p = 1; p <= L1; p++) if (S >> (p - 1) & 1) mem[S][nmem[S]++] = p;
    }
    for (int t = 0; t < nthreads; t++) {
        works[t] = malloc(sizeof(Work));
        Work *w = works[t];
        for (int S = 0; S < 128; S++) for (int u = 0; u < MAXK; u++) for (int x = 0; x < MAXK; x++) w->f[S][u][x] = CINF;
        for (int u = 0; u < MAXK; u++) for (int x = 0; x < NV; x++) w->H[u][x] = CINF;
        for (int i = 0; i < 21; i++) for (int p = 0; p < NV; p++) for (int q = 0; q < NV; q++) w->KK[i][p][q] = CINF;
    }
}

/* k = 2: one order */
static void sums_k2(void)
{
    for (int g1 = 0; g1 < m; g1++)
        for (int g2 = 0; g2 < m; g2++) {
            int kv[MAXM] = {0};
            kv[g1] = 1;
            const u64 *r1 = row_of(kv);
            kv[g1] = 0; kv[g2] = 1;
            const u64 *r2 = row_of(kv);
            u128 s = 0;
            for (int a0 = 0; a0 < n; a0++)
                for (int a1 = 0; a1 < n; a1++)
                    for (int a2 = 0; a2 < n; a2++) s += (u128)r1[iabs(a0 - a1)] + r2[iabs(a1 - a2)];
            Ssum[g1 + g2 * m] = s;
            Sset[g1 + g2 * m]++;
        }
}

/* --orders: all (k-1)! orders in every cell */
static int nperm, *perm;

static void make_perms(void)
{
    int q = k - 1, a[MAXK];
    nperm = 1;
    for (int i = 2; i <= q; i++) nperm *= i;
    perm = malloc((size_t)nperm * q * sizeof(int));
    for (int i = 0; i < q; i++) a[i] = i + 1;
    for (int t = 0; t < nperm; t++) {
        memcpy(perm + (size_t)t * q, a, q * sizeof(int));
        int i = q - 2;                       /* next permutation in lexicographic order */
        while (i >= 0 && a[i] > a[i + 1]) i--;
        if (i < 0) break;
        int j = q - 1;
        while (a[j] < a[i]) j--;
        int tmp = a[i]; a[i] = a[j]; a[j] = tmp;
        for (int x = i + 1, y = q - 1; x < y; x++, y--) { tmp = a[x]; a[x] = a[y]; a[y] = tmp; }
    }
}

static void gap_rows(long gi, const u64 *row[MAXK + 1][MAXK + 1])
{
    int g[MAXK + 1], kv[MAXM];
    long r = gi;
    for (int l = 1; l <= k; l++) { g[l] = (int)(r % m); r /= m; }
    for (int i = 0; i <= k; i++)
        for (int j = i + 1; j <= k; j++) {
            memset(kv, 0, sizeof kv);
            for (int l = i + 1; l <= j; l++) kv[g[l]]++;
            row[i][j] = row_of(kv);
        }
}

static void orders_job(long gi, int tid, void *ctx)
{
    (void)tid; (void)ctx;
    const u64 *row[MAXK + 1][MAXK + 1];
    gap_rows(gi, row);
    int a[MAXK + 1] = {0}, q = k - 1;
    u64 E[MAXK + 1][MAXK + 1];
    u128 S = 0;
    for (;;) {
        for (int i = 0; i <= k; i++)
            for (int j = i + 1; j <= k; j++) E[i][j] = E[j][i] = row[i][j][iabs(a[i] - a[j])];
        u64 best = UINT64_MAX;
        for (int t = 0; t < nperm; t++) {
            const int *pi = perm + (size_t)t * q;
            u64 c = E[0][pi[0]] + E[pi[q - 1]][k];
            for (int i = 0; i + 1 < q; i++) c += E[pi[i]][pi[i + 1]];
            if (c < best) best = c;
        }
        S += best;
        int i = 0;
        while (i <= k && ++a[i] == n) a[i++] = 0;
        if (i > k) break;
    }
    Ssum[gi] = S;
    Sset[gi]++;
    progress();
}

/* --hk: Held-Karp from point 0 in every cell of the listed gap tuples */
static long *hkg;
static int nhkg;
static u128 *hkpart;

static void hk_job(long ci, int tid, void *ctx)
{
    (void)tid; (void)ctx;
    long gi = hkg[ci / (n * n)];
    int chunk = (int)(ci % (n * n));
    const u64 *row[MAXK + 1][MAXK + 1];
    gap_rows(gi, row);
    int a[MAXK + 1] = {0}, q = k - 1, full = (1 << q) - 1;
    u64 E[MAXK + 1][MAXK + 1];
    static _Thread_local u64 dp[256][MAXK];
    u128 S = 0;
    a[0] = chunk % n;
    a[1] = chunk / n;
    for (;;) {
        for (int i = 0; i <= k; i++)
            for (int j = i + 1; j <= k; j++) E[i][j] = E[j][i] = row[i][j][iabs(a[i] - a[j])];
        for (int S1 = 1; S1 <= full; S1++)
            for (int v = 1; v <= q; v++) {
                if (!(S1 >> (v - 1) & 1)) continue;
                int R = S1 & ~(1 << (v - 1));
                if (!R) { dp[S1][v] = E[0][v]; continue; }
                u64 best = UINT64_MAX;
                for (int x = 1; x <= q; x++)
                    if (R >> (x - 1) & 1) { u64 c = dp[R][x] + E[x][v]; if (c < best) best = c; }
                dp[S1][v] = best;
            }
        u64 best = UINT64_MAX;
        for (int v = 1; v <= q; v++) { u64 c = dp[full][v] + E[v][k]; if (c < best) best = c; }
        S += best;
        int i = 2;
        while (i <= k && ++a[i] == n) a[i++] = 0;
        if (i > k) break;
    }
    hkpart[ci] = S;
}

/* ---------- main ---------- */
int main(int argc, char **argv)
{
    double t0 = now();
    if (argc < 2) die(2, "usage: c12a_check SPEC [-j N] [--orders] [--hk G1,G2,...]");
    const char *spec = argv[1];
    int orders = 0;
    const char *hklist = NULL;
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "-j") && i + 1 < argc) {
            nthreads = atoi(argv[++i]);
            if (nthreads < 1 || nthreads > 64) die(2, "-j N needs 1 <= N <= 64");
        } else if (!strcmp(argv[i], "--orders")) {
            orders = 1;
        } else if (!strcmp(argv[i], "--hk") && i + 1 < argc) {
            hklist = argv[++i];
        } else {
            die(2, "unknown argument %s", argv[i]);
        }
    }

    /* steps 1 to 4 of Section 10 */
    char *vpath = read_spec(spec);
    Text v = read_text(vpath);
    const char *vp = vpath;
    if (v.l[0].nf != 2 || strcmp(v.l[0].f[0], "tsp-block-values") || strcmp(v.l[0].f[1], "1"))
        die(2, "%s: line 1 is not 'tsp-block-values 1'", vp);
    if (v.n < 12) die(2, "%s: too few lines", vp);
    for (int i = 0; i < 7; i++) {
        Line *ln = &v.l[i + 1];
        if (strcmp(ln->f[0], pname[i])) die(2, "%s: line %d is not '%s'", vp, i + 2, pname[i]);
        if (pcnt[i] >= 0 && ln->nf - 1 != pcnt[i]) die(2, "%s: line %d needs %d numbers", vp, i + 2, pcnt[i]);
        for (int j = 1; j < ln->nf; j++) tonum(ln->f[j], vp, i + 2);
    }
    if (strcmp(v.l[8].f[0], "masses")) die(2, "%s: line 9 is not 'masses'", vp);
    for (int j = 1; j < v.l[8].nf; j++) tonum(v.l[8].f[j], vp, 9);
    if (strcmp(v.l[9].f[0], "tail_term") || v.l[9].nf != 2) die(2, "%s: line 10 is not 'tail_term N'", vp);
    u128 sET = tonum(v.l[9].f[1], vp, 10);
    if (strcmp(v.l[10].f[0], "table") || v.l[10].nf != 2) die(2, "%s: line 11 is not 'table N'", vp);
    u128 ntab = tonum(v.l[10].f[1], vp, 11);
    if (ntab > (u128)(v.n - 11)) die(2, "%s: fewer table lines than stated", vp);
    int tab0 = 11, sl = 11 + (int)ntab;
    for (int i = tab0; i < sl; i++) for (int j = 0; j < v.l[i].nf; j++) tonum(v.l[i].f[j], vp, i + 1);
    if (sl >= v.n || strcmp(v.l[sl].f[0], "sums") || v.l[sl].nf != 2) die(2, "%s: no 'sums N' line", vp);
    u128 nsum = tonum(v.l[sl].f[1], vp, sl + 1);
    if (nsum > (u128)(v.n - sl - 1)) die(2, "%s: fewer sums than stated", vp);
    int sum0 = sl + 1, tl = sum0 + (int)nsum;
    for (int i = sum0; i < tl; i++) {
        if (v.l[i].nf != 1) die(2, "%s: line %d is not one number", vp, i + 1);
        tonum(v.l[i].f[0], vp, i + 1);
    }
    if (tl >= v.n || strcmp(v.l[tl].f[0], "total") || v.l[tl].nf != 2) die(2, "%s: no 'total N' line", vp);
    u128 stot = tonum(v.l[tl].f[1], vp, tl + 1);
    if (tl + 1 != v.n) die(2, "%s: lines after 'total'", vp);

    printf("spec %s\nvalues %s\n", spec, vpath);
    printf("parameters: k %d, h %llu/%llu, mesh %llu, gaps", k, (unsigned long long)h1, (unsigned long long)h2,
           (unsigned long long)M);
    for (int i = 1; i < m; i++) printf(" %llu", (unsigned long long)b[i]);
    printf(" (m = %d), heights %d, height_mesh %d, tail %d, threshold %llu/%llu\n", m, n, nu, T,
           (unsigned long long)th1, (unsigned long long)th2);
    for (int i = 0; i < 7; i++) {
        Line *ln = &v.l[i + 1];
        int same = ln->nf - 1 == speccnt[i];
        for (int j = 1; same && j < ln->nf; j++) same = tonum(ln->f[j], vp, i + 2) == specnum[i][j - 1];
        if (!same) {
            printf("REJECTED: parameter line '%s' of the values file differs from the spec file\n", pname[i]);
            return 1;
        }
    }
    mpw[0] = 1;
    for (int i = 1; i <= k + 1; i++) mpw[i] = mpw[i - 1] * m;
    NK = (int)(binom(m + k, k) - 1);
    if (v.l[8].nf - 1 != m) die(2, "%s: %d masses, not m = %d", vp, v.l[8].nf - 1, m);
    if (ntab != (u128)NK) die(2, "%s: table count is not C(m + k, k) - 1 = %d", vp, NK);
    for (int i = tab0; i < sl; i++) if (v.l[i].nf != m + n) die(2, "%s: line %d has not m + n numbers", vp, i + 1);
    if (nsum != (u128)mpw[k]) die(2, "%s: sums count is not m^k = %ld", vp, mpw[k]);

    /* Sections 2 to 6 */
    section2();
    section3();
    double t1 = now();
    section4();
    double t2 = now();
    section5();
    double t3 = now();
    section6();
    double t4 = now();
    fprintf(stderr, "sections 2-3 %.2f s, 4 %.2f s, 5 %.2f s, 6 %.2f s\n", t1 - t0, t2 - t1, t3 - t2, t4 - t3);

    long bad = 0, compared = 0;
    char s1[48], s2[48];
    long badm = 0;
    for (int c = 0; c < m; c++) {
        u128 sv = tonum(v.l[8].f[c + 1], vp, 9);
        compared++;
        if (sv != P[c]) {
            if (badm++ < 5) printf("differs: mass P_%d stored %s, recomputed %s\n", c, u128s(sv, s1), u128s(P[c], s2));
        }
    }
    printf("masses: %d compared, %d equal\n", m, m - (int)badm);
    compared++;
    int badet = sET != ET;
    if (badet) printf("differs: tail_term stored %s, recomputed %s\n", u128s(sET, s1), u128s(ET, s2));
    printf("tail_term: 1 compared, %s (%s)\n", badet ? "differs" : "equal", u128s(ET, s1));
    long badt = 0, badtl = 0;
    for (int t = 0; t < NK; t++) {
        Line *ln = &v.l[tab0 + t];
        int lbad = 0;
        for (int c = 0; c < m; c++) { compared++; if (tonum(ln->f[c], vp, tab0 + t + 1) != key[t][c]) { lbad = 1; badt++; } }
        for (int d = 0; d < n; d++) {
            compared++;
            u128 sv = tonum(ln->f[m + d], vp, tab0 + t + 1);
            if (sv != Tab[(size_t)t * n + d]) {
                badt++;
                if (badtl < 5)
                    printf("differs: table line %d entry %d stored %s, recomputed %s\n", t + 1, d, u128s(sv, s1),
                           u128s(Tab[(size_t)t * n + d], s2));
                lbad = 1;
            }
        }
        badtl += lbad;
    }
    printf("table: %d lines of m + n = %d numbers (key and Tab(key, 0..n-1)), %ld numbers compared, %ld lines differ\n",
           NK, m + n, (long)NK * (m + n), badtl);
    bad += badm + badet + badt;

    /* Section 7 */
    Ssum = calloc((size_t)mpw[k], sizeof(u128));
    Sset = calloc((size_t)mpw[k], 1);
    if (k * (u128)maxTab >= ((u128)1 << 60) || (u128)n * n * n * k * maxTab >= TWO64)
        die(2, "Tab too large for the 64-bit costs of Section 7");
    double t5 = now();
    atomic_init(&done7, 0);
    t7start = t5;
    if (hklist) {
        nhkg = 0;
        hkg = malloc(strlen(hklist) * sizeof(long));
        for (const char *p = hklist; *p;) {
            char *e;
            long gi = strtol(p, &e, 10);
            if (e == p || gi < 0 || gi >= mpw[k]) die(2, "bad --hk list");
            hkg[nhkg++] = gi;
            p = *e == ',' ? e + 1 : e;
            if (*e && *e != ',') die(2, "bad --hk list");
        }
        hkpart = calloc((size_t)nhkg * n * n, sizeof(u128));
        parallel_for((long)nhkg * n * n, hk_job, NULL);
        long hbad = 0;
        for (int i = 0; i < nhkg; i++) {
            u128 s = 0;
            for (int c = 0; c < n * n; c++) s += hkpart[(size_t)i * n * n + c];
            u128 sv = tonum(v.l[sum0 + hkg[i]].f[0], vp, sum0 + (int)hkg[i] + 1);
            printf("cell sum S_%ld (Held-Karp in every cell): %s, stored %s, %s\n", hkg[i], u128s(s, s1),
                   u128s(sv, s2), s == sv ? "equal" : "DIFFERS");
            hbad += s != sv;
        }
        printf("section 7 run time %.1f s, %d threads\n", now() - t5, nthreads);
        printf("RESULT: %s: %d cell sums recomputed by Held-Karp in every cell, %ld differ (masses, tail term and "
               "table: %ld numbers differ)\n", hbad || bad ? "REJECTED" : "AGREE", nhkg, hbad, bad);
        return hbad || bad ? 1 : 0;
    }
    const char *method;
    if (orders) {
        make_perms();
        todo7 = mpw[k];
        parallel_for(mpw[k], orders_job, NULL);
        method = "all (k-1)! orders in every cell, no identities";
    } else if (k == 2) {
        sums_k2();
        method = "one order (k = 2)";
    } else {
        setup7();
        long nt = mpw[k - 2], canon = 0;
        for (long ti = 0; ti < nt; ti++) {
            long r = ti;
            int g[MAXK + 2], cmp = 0;
            for (int l = 2; l <= L; l++) { g[l] = (int)(r % m); r /= m; }
            for (int l = 2; l <= L && !cmp; l++) if (g[l] != g[L + 2 - l]) cmp = g[l] < g[L + 2 - l] ? -1 : 1;
            canon += cmp <= 0;
        }
        todo7 = canon;
        parallel_for(nt, task7, NULL);
        method = "shared paths, reversal, reflection on a_1";
    }
    for (long i = 0; i < mpw[k]; i++) if (Sset[i] != 1) die(2, "internal: cell sum %ld set %d times", i, Sset[i]);
    double t6 = now();
    fprintf(stderr, "section 7 %.1f s\n", t6 - t5);
    long bads = 0;
    for (long i = 0; i < mpw[k]; i++) {
        u128 sv = tonum(v.l[sum0 + i].f[0], vp, sum0 + (int)i + 1);
        compared++;
        if (sv != Ssum[i]) {
            if (bads < 5) printf("differs: cell sum S_%ld stored %s, recomputed %s\n", i, u128s(sv, s1), u128s(Ssum[i], s2));
            bads++;
        }
    }
    printf("sums: %ld compared, %ld equal (%s)\n", mpw[k], mpw[k] - bads, method);
    bad += bads;

    /* Section 8 */
    u128 *A = malloc((size_t)mpw[k] * sizeof(u128));
    memcpy(A, Ssum, (size_t)mpw[k] * sizeof(u128));
    for (int j = k; j >= 1; j--) {
        long cnt = mpw[j - 1];
        for (long i = 0; i < cnt; i++) {
            u128 s = 0;
            for (int c = 0; c < m; c++) s = add(s, mul(P[c], A[i + c * cnt]));
            A[i] = cdiv(s, DD);
        }
    }
    u128 A0 = A[0];
    compared++;
    int badtot = stot != A0;
    if (badtot) printf("differs: total stored %s, recomputed %s\n", u128s(stot, s1), u128s(A0, s2));
    printf("total: 1 compared, %s (%s)\n", badtot ? "differs" : "equal", u128s(A0, s1));
    bad += badtot;
    u128 den = mul(mul(FF, (u128)k), (u128)h1);
    for (int i = 0; i <= k; i++) den = mul(den, (u128)n);
    u128 bnd = cdiv(mul(mul(A0, h2), (u128)1000000000000ULL), den);
    u128 lhs = mul(mul(A0, h2), th2), rhs = mul(th1, den);
    printf("bound %s.%012llu (rounded up)\n", u128s(bnd / 1000000000000ULL, s1),
           (unsigned long long)(bnd % 1000000000000ULL));
    printf("test: A_0 h_2 theta_2 = %s %s theta_1 n^(k+1) F k h_1 = ", u128s(lhs, s1), lhs <= rhs ? "<=" : ">");
    printf("%s\n", u128s(rhs, s2));
    printf("numbers compared: %ld, differ: %ld\n", compared, bad);
    printf("run time %.1f s (section 7 %.1f s), %d threads\n", now() - t0, t6 - t5, nthreads);
    if (bad) { printf("RESULT: REJECTED: %ld stored numbers differ from the recomputed ones\n", bad); return 1; }
    if (lhs > rhs) {
        printf("RESULT: REJECTED: the stored values agree, but B > %llu/%llu\n", (unsigned long long)th1,
               (unsigned long long)th2);
        return 1;
    }
    printf("RESULT: ACCEPTED: all %ld stored numbers equal the recomputed ones and B <= %llu/%llu\n", compared,
           (unsigned long long)th1, (unsigned long long)th2);
    return 0;
}
