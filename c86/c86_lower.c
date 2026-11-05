/* c86_lower.c: the dual inequality behind a lower bound for the Schur-Siegel-Smyth trace constant (86a page, the
 * 1.80220* row of PR #315, arXiv:2609.29298): for every x >= 0,
 *     g(x) = x - sum_Q lambda_Q log|Q(x)| - lambda_0 (2 U(x) - I) >= lambda,
 * with eta a probability measure with constant density on each cell [a_j, a_{j+1}] of [0, 8] (mass m_j), U its
 * logarithmic potential and I its energy.
 *
 *   c86_lower energy < data      a rigorous lower bound for I(eta) (as given, masses unnormalised)
 *   c86_lower check X0 X1 ILO LAM < data
 *                                certify g >= LAM on [X0, X1] by bisection, with ILO a lower bound for I
 *
 * Arithmetic is binary64. Basic operations are correctly rounded, and log() is taken to be within 1 ulp (macOS libm
 * documents this); every quantity carries an explicit error bound that covers 8 ulp per log and the standard
 * bounds for sums, and the certified margin (2.1e-5) is about eleven orders of magnitude above these bounds.
 *
 * Potentials. Cell j has width w and centre c; u_j(x) = (1/w) int log|x - y| dy over the cell. If the cell is far,
 * |x - c| >= 64 w, u_j(x) = log|x - c| - w^2/(24 (x - c)^2) + r with |r| <= w^4/(36 (x - c)^4) (the series of
 * log(1 + z) averaged over a symmetric z with |z| <= 1/128); otherwise u_j(x) = (Psi(x - a) - Psi(x - b)) / w with
 * Psi(t) = t log|t| - t. u_j is unimodal with its minimum at c, so on [s, t] it is at most max(u_j(s), u_j(t)).
 * Energy. I = sum_ij (m_i/w_i)(m_j/w_j) E_ij, E_ii = w^2 (log w - 3/2), E_ij for i != j by the double antiderivative
 * Phi(t) = t^2/2 log|t| - 3t^2/4 when the cells are near and by the same series (z^2 terms) when they are far.
 * Polynomials. max |Q| on [s, t] from the centred Taylor form of degree 3 with a crude fourth-order remainder.
 */
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define U53 (DBL_EPSILON / 2)   /* unit roundoff 2^-53 */

static int n, nq;
static double *edge, *mass, *wid, *cen;
static int *deg;
static double **coef, *lam, lam0;
static int *active, nact;     /* cells with positive mass */

static double up(double x) { return nextafter(x, INFINITY); }
static double dn(double x) { return nextafter(x, -INFINITY); }

/* u_j(x) with an error bound err (absolute) */
static double ucell(int j, double x, double *err) {
    double w = wid[j], c = cen[j], d = x - c, ad = fabs(d);
    if (ad >= 64 * w) {
        double L = log(ad), q = w / d, q2 = q * q;
        double v = L - q2 / 24;
        /* remainder <= q^4/36; rounding: log within 1 ulp, a few ops */
        *err = q2 * q2 / 36 + 8 * U53 * (fabs(L) + 1) + 1e-300;
        return v;
    }
    double a = edge[j] - x, b = edge[j + 1] - x;   /* int_{a}^{b} log|t| dt over t = y - x */
    double Pb = (b == 0) ? 0 : b * log(fabs(b)) - b;
    double Pa = (a == 0) ? 0 : a * log(fabs(a)) - a;
    double v = (Pb - Pa) / w;
    double mag = (fabs(b) * (fabs(log(fabs(b) + 1e-300)) + 1) + fabs(a) * (fabs(log(fabs(a) + 1e-300)) + 1)) / w;
    *err = 16 * U53 * mag + 4 * U53 * fabs(v) + 1e-300;
    return v;
}

static double Phi(double t) { return t == 0 ? 0 : t * t * (0.5 * log(fabs(t)) - 0.75); }

/* a lower bound for I(eta), as sum and error */
static void energy(double *Iest, double *Ierr) {
    double s = 0, e = 0, sabs = 0;
    for (int ii = 0; ii < nact; ii++) {
        int i = active[ii];
        double mi = mass[i], wi = wid[i];
        double tii = mi * mi * (log(wi) - 1.5);
        s += tii;
        sabs += fabs(tii);
        e += 8 * U53 * mi * mi * (fabs(log(wi)) + 1.5);
        for (int jj = ii + 1; jj < nact; jj++) {
            int j = active[jj];
            double mj = mass[j], wj = wid[j], D = cen[j] - cen[i], aD = fabs(D);
            double t, te;
            if (aD >= 64 * (wi + wj)) {
                double L = log(aD), z2 = (wi * wi + wj * wj) / (D * D);
                double avg = L - z2 / 24;
                double zmax = (wi + wj) / (2 * aD);           /* <= 1/128 */
                double rem = (8.0 / 3.0) * zmax * z2 / 12;     /* bound on the averaged cubic-and-higher terms */
                t = 2 * mi * mj * avg;
                te = 2 * mi * mj * (rem + 8 * U53 * (fabs(L) + 1));
            } else {
                double a = edge[i], b = edge[i + 1], c = edge[j], d = edge[j + 1];
                double E = -Phi(b - d) + Phi(b - c) + Phi(a - d) - Phi(a - c);
                double mag = fabs(Phi(b - d)) + fabs(Phi(b - c)) + fabs(Phi(a - d)) + fabs(Phi(a - c));
                double f = mi * mj / (wi * wj);
                t = 2 * f * E;
                te = 2 * f * (16 * U53 * (mag + 4 * fmax(fmax(fabs(b - d), fabs(a - c)), 1e-300) *
                                          fmax(fmax(fabs(b - d), fabs(a - c)), 1e-300))) + 8 * U53 * fabs(t);
            }
            s += t;
            sabs += fabs(t);
            e += te;
        }
    }
    double nt = (double)nact * nact;
    *Iest = s;
    *Ierr = e + nt * U53 / (1 - nt * U53) * sabs;     /* gamma_n bound for the recursive sum */
}

/* upper bound of max |Q| on [s, t] by the centred Taylor form: with m the midpoint and r the half-width,
 * |Q(x)| <= sum_{j<4} (|c_j| + e_j) r^j + r^4 sum_i |a_i| C(i, 4) M^(i-4), c_j = Q^(j)(m)/j! by repeated synthetic
 * division, e_j their rounding bounds (from the same division on |a_i| and |m|), M = max(|s|, |t|). Interval Horner
 * would overestimate by about sum |a_i| i M^(i-1) (t - s), hopeless for degree 10 near x = 5. */
static double qmax(int k, double s, double t) {
    int d = deg[k];
    double m = 0.5 * (s + t), r = up(0.5 * (t - s)) + fabs(m) * U53;   /* r covers the rounding of m */
    double b[32], A[32];
    for (int i = 0; i <= d; i++) { b[i] = coef[k][i]; A[i] = fabs(coef[k][i]); }
    double am = fabs(m), tot = 0, rp = 1;
    int nn = d;
    for (int j = 0; j < 4 && nn >= 0; j++) {
        /* synthetic division of b[0..nn] by (x - m): quotient in b[0..nn-1], remainder c_j */
        double acc = b[0], aac = A[0];
        for (int i = 1; i <= nn; i++) {
            double nb = acc * m + b[i], nA = aac * am + A[i];
            b[i - 1] = acc; A[i - 1] = aac;
            acc = nb; aac = nA;
        }
        double ej = 4 * (d + 1) * U53 * aac + 1e-300;
        tot += (fabs(acc) + ej) * rp;
        rp *= r;
        nn--;
    }
    if (d >= 4) {
        double M = fmax(fabs(s), fabs(t)), R = 0;
        for (int i = 4; i <= d; i++) {          /* coefficient of x^i is coef[k][d - i] */
            double cb = (double)i * (i - 1) * (i - 2) * (i - 3) / 24.0;
            R += fabs(coef[k][d - i]) * cb * pow(M, i - 4);
        }
        tot += R * rp;                          /* rp = r^4 here */
    }
    return tot * (1 + 16 * U53);
}

/* a rigorous lower bound for g on [s, t] (0 <= s <= t). A cell whose centre is outside (s, t) is monotone on [s, t],
 * so its maximum there is at the endpoint farther from the centre; otherwise both endpoints are evaluated. */
static double glow(double s, double t, double ILO, double *Uup) {
    double Us = 0, err = 0;
    for (int jj = 0; jj < nact; jj++) {
        int j = active[jj];
        double c = cen[j], e1, e2, v;
        if (c <= s) v = ucell(j, t, &e1);
        else if (c >= t) v = ucell(j, s, &e1);
        else {
            double v1 = ucell(j, s, &e1), v2 = ucell(j, t, &e2);
            v = v1 > v2 ? v1 : v2;
            e1 = fmax(e1, e2);
        }
        Us += mass[j] * v;
        err += mass[j] * e1 + 4 * U53 * mass[j] * fabs(v);
    }
    err += nact * U53 * 2 * fabs(Us) + 1e-15;
    double Uhi = Us + err;
    *Uup = Uhi;
    double P = 0, perr = 0;
    for (int k = 0; k < nq; k++) {
        if (lam[k] == 0) continue;
        double M = qmax(k, s, t);
        if (M <= 0) M = 1e-300;
        double L = log(M);
        P += lam[k] * L;
        perr += lam[k] * (8 * U53 * fabs(L) + 1e-300) + 4 * U53 * fabs(lam[k] * L);
    }
    double g = s - (P + perr) - lam0 * (2 * Uhi) + lam0 * ILO;
    return g - 64 * U53 * (fabs(s) + fabs(P) + 2 * lam0 * fabs(Uhi) + lam0 * fabs(ILO) + 1);
}

static int readall(void) {
    if (scanf("%d", &n) != 1 || n < 1) return 0;
    edge = malloc(sizeof(double) * (n + 1));
    mass = malloc(sizeof(double) * n);
    wid = malloc(sizeof(double) * n);
    cen = malloc(sizeof(double) * n);
    for (int i = 0; i <= n; i++) if (scanf("%la", &edge[i]) != 1) return 0;
    for (int i = 0; i < n; i++) if (scanf("%la", &mass[i]) != 1) return 0;
    if (scanf("%d", &nq) != 1 || nq < 1) return 0;
    deg = malloc(sizeof(int) * nq);
    coef = malloc(sizeof(double *) * nq);
    lam = malloc(sizeof(double) * nq);
    for (int k = 0; k < nq; k++) {
        if (scanf("%d", &deg[k]) != 1 || deg[k] < 1) return 0;
        coef[k] = malloc(sizeof(double) * (deg[k] + 1));
        for (int i = 0; i <= deg[k]; i++) if (scanf("%la", &coef[k][i]) != 1) return 0;
        if (scanf("%la", &lam[k]) != 1 || lam[k] < 0) return 0;
    }
    if (scanf("%la", &lam0) != 1 || lam0 < 0) return 0;
    active = malloc(sizeof(int) * n);
    nact = 0;
    for (int i = 0; i < n; i++) {
        wid[i] = edge[i + 1] - edge[i];
        cen[i] = edge[i] + wid[i] / 2;
        if (!(wid[i] > 0) || mass[i] < 0) return 0;
        if (mass[i] > 0) active[nact++] = i;
    }
    return 1;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: c86_lower energy | check X0 X1 ILO LAM < data\n"); return 2; }
    if (!readall()) { fprintf(stderr, "bad data\n"); return 2; }
    if (!strcmp(argv[1], "energy")) {
        double I, e;
        energy(&I, &e);
        printf("I %.17g\nerr %.3e\nlo %a\n", I, e, dn(dn(I - e)));
        return 0;
    }
    if (argc != 6 || strcmp(argv[1], "check")) { fprintf(stderr, "usage\n"); return 2; }
    double X0 = atof(argv[2]), X1 = atof(argv[3]), ILO = strtod(argv[4], NULL), LAM = strtod(argv[5], NULL);
    /* work stack of intervals over [X0, X1]; the cell edges and the polynomial roots need no special handling: the
     * bounds hold on any interval, they are only looser there, and bisection tightens them */
    long cap = 1 << 22, top = 0, done = 0;
    double *S = malloc(sizeof(double) * cap), *T = malloc(sizeof(double) * cap);
    double worst = INFINITY, worst_s = 0, best_mid = INFINITY, best_mid_x = 0;
    int nstart = 1 << 10;
    for (int i = nstart - 1; i >= 0; i--) {
        S[top] = X0 + (X1 - X0) * i / nstart;
        T[top] = X0 + (X1 - X0) * (i + 1) / nstart;
        top++;
    }
    T[0] = X1;
    S[top - 1] = X0;
    while (top > 0) {
        top--;
        double s = S[top], t = T[top], Uup;
        double g = glow(s, t, ILO, &Uup);
        if (g >= LAM) {
            done++;
            if (g < worst) { worst = g; worst_s = s; }
            if (t - s < 1e-5) {      /* a point value for the report: the bound at the midpoint alone */
                double m = 0.5 * (s + t), u2;
                double gm = glow(m, m, ILO, &u2);
                if (gm < best_mid) { best_mid = gm; best_mid_x = m; }
            }
            if ((done & ((1 << 20) - 1)) == 0) fprintf(stderr, "  %ld intervals, at x = %.6f\n", done, s);
            continue;
        }
        if (t - s < 1e-13) {
            printf("FAIL at [%.17g, %.17g]: lower bound %.12g < %.12g\n", s, t, g, LAM);
            return 1;
        }
        if (top + 2 > cap) { fprintf(stderr, "stack overflow\n"); return 3; }
        double m = 0.5 * (s + t);
        S[top] = m; T[top] = t; top++;
        S[top] = s; T[top] = m; top++;
    }
    printf("PASS intervals %ld\nsmallest_interval_bound %.15f at %.12f\nsmallest_point_value %.15f at %.12f\n", done,
           worst, worst_s, best_mid, best_mid_x);
    return 0;
}
