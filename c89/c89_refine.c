/* c89_refine.c: the same exact counts as c89_multiset.c, computed another way (method B; no shared code).
 *
 * The vertices are placed one at a time. Before vertex v is placed, the n coordinates fall into cells by the bits of
 * the vertices already placed (vertex 0 is the origin, so it never splits a cell); placing v means choosing, in each
 * cell of size s, how many j of its coordinates get bit 1, in binom(s, j) ways. The distances from v to the earlier
 * vertices are then final, so a branch stops as soon as both colours are impossible (a pair of class B rules out
 * red, a pair of class R rules out blue); the weight it stands for is added to the total with the factor 2^(n(t-1-v))
 * for the vertices not yet placed, and the total must again be 2^(n(t-1)).
 *
 * usage: c89_refine t n CLASSES   (CLASSES as for c89_multiset); prints the same lines as c89_multiset.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXT 6
#define MAXN 16

static int t, n;
static int cls[MAXN + 1]; /* -1 = R, -2 = B, k >= 0 fractional index */
static uint64_t C[MAXN + 1][MAXN + 1];
static uint64_t red_count[16][16][16][16], blue_count[16][16][16][16], total;

/* place vertex v: cells s[0..L-1] (L = 2^(v-1)); c = the cell being split; jv = choices made so far */
static void place(int v, const int *s, int c, int *jv, uint64_t w, int red_ok, int blue_ok, const int *ex) {
    int L = 1 << (v - 1);
    if (c == L) {
        int e2[4] = {ex[0], ex[1], ex[2], ex[3]}, r = red_ok, b = blue_ok;
        for (int u = 0; u < v; u++) {
            int d = 0;
            for (int k = 0; k < L; k++) {
                int bit_u = u ? (k >> (u - 1)) & 1 : 0;
                d += bit_u ? s[k] - jv[k] : jv[k];
            }
            int cl = cls[d];
            if (cl == -2) r = 0;
            else if (cl == -1) b = 0;
            else e2[cl]++;
        }
        if (!r && !b) {
            total += w << (n * (t - 1 - v));
            return;
        }
        if (v == t - 1) {
            if (r) red_count[e2[0]][e2[1]][e2[2]][e2[3]] += w;
            if (b) blue_count[e2[0]][e2[1]][e2[2]][e2[3]] += w;
            total += w;
            return;
        }
        int s2[1 << (MAXT - 1)], j2[1 << (MAXT - 1)];
        for (int k = 0; k < L; k++) {
            s2[k] = s[k] - jv[k];
            s2[k + L] = jv[k];
        }
        place(v + 1, s2, 0, j2, w, r, b, e2);
        return;
    }
    for (int j = 0; j <= s[c]; j++) {
        jv[c] = j;
        place(v, s, c + 1, jv, w * C[s[c]][j], red_ok, blue_ok, ex);
    }
}

int main(int argc, char **argv) {
    if (argc != 4) { fprintf(stderr, "usage: c89_refine t n CLASSES\n"); return 2; }
    t = atoi(argv[1]); n = atoi(argv[2]);
    if (t < 4 || t > MAXT || n < 1 || n > MAXN) { fprintf(stderr, "bad t or n\n"); return 2; }
    char *s = strdup(argv[3]), *tok = strtok(s, ",");
    for (int w = 0; w <= n; w++) {
        if (!tok) { fprintf(stderr, "need %d classes\n", n + 1); return 2; }
        if (!strcmp(tok, "R")) cls[w] = -1;
        else if (!strcmp(tok, "B")) cls[w] = -2;
        else if (tok[0] >= '0' && tok[0] <= '3' && !tok[1]) cls[w] = tok[0] - '0';
        else { fprintf(stderr, "bad class %s\n", tok); return 2; }
        tok = strtok(NULL, ",");
    }
    if (tok) { fprintf(stderr, "too many classes\n"); return 2; }
    for (int a = 0; a <= MAXN; a++)
        for (int b = 0; b <= MAXN; b++)
            C[a][b] = (b == 0) ? 1 : (a == 0) ? 0 : C[a - 1][b - 1] + C[a - 1][b];
    int s0[1] = {n}, j0[1], ex0[4] = {0, 0, 0, 0};
    place(1, s0, 0, j0, 1, 1, 1, ex0);
    for (int a = 0; a < 16; a++)
        for (int b = 0; b < 16; b++)
            for (int c = 0; c < 16; c++)
                for (int d = 0; d < 16; d++) {
                    if (red_count[a][b][c][d]) printf("R %d %d %d %d %llu\n", a, b, c, d,
                                                      (unsigned long long)red_count[a][b][c][d]);
                    if (blue_count[a][b][c][d]) printf("B %d %d %d %d %llu\n", a, b, c, d,
                                                       (unsigned long long)blue_count[a][b][c][d]);
                }
    printf("TOTAL %llu\n", (unsigned long long)total);
    return 0;
}
