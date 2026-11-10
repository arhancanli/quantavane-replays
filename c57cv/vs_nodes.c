/* vs_nodes.c: interval enclosures of f_i(x_k) at the negative zeros x_k of T_M, for the 57c [vS2026] replay
   (vs_conformal_check.py writes the input and adds the quadrature remainders). Our own code, written from the paper's
   formula (4.4):
     f_i(x) = -a_i x + sum_n c_{i,n} (pi/n) T_n(x) - sum_j (pi/N_ij) sum_l log|G~_ij(x, y_l)| psi_j(y_l)
   with y_l = cos((2l-1) pi/(2N)), G~_ij = G_ij = 2 sinh((w_i(x) - w_j(y))/2) for j != i, and for j = i
   |K_i|^2 = alpha^2 shc(tau alpha)^2 + beta^2 sinc(tau beta)^2, tau = (x - y)/2, alpha = a_i,
   beta = B1 (x + y) + B2 (x^3 + x^2 y + x y^2 + y^3) (B1 = b1 a^2, B2 = b2 a^4), which is |Delta shc(tau Delta)|^2.
   For j != i (x, y real): |G|^2 = 4 sinh^2(U/2) + 4 sin^2(V/2), U = a_i x - a_j y,
   V = Theta_i - Theta_j + B1_i x^2 + B2_i x^4 - B1_j y^2 - B2_j y^4.
   Interval arithmetic in binary64: every + - * / result is widened by one ulp each way (IEEE operations are correctly
   rounded); libm's sin, cos, sinh, log are taken to be within one ulp and widened by four. Angles of the cosines are
   reduced exactly in integers (n (2k+1) mod 4M).
   Input on stdin (hex floats, exact), output one line per node: k  lo  hi  (hex).  Usage: vs_nodes < in > out */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct { double lo, hi; } I;
static double dn(double x) { return nextafter(x, -INFINITY); }
static double up(double x) { return nextafter(x, INFINITY); }
static double dn4(double x) { return dn(dn(dn(dn(x)))); }
static double up4(double x) { return up(up(up(up(x)))); }
static I mk(double a, double b) { I r = {a, b}; return r; }
static I pt(double a) { return mk(a, a); }
static I add(I a, I b) { return mk(dn(a.lo + b.lo), up(a.hi + b.hi)); }
static I sub(I a, I b) { return mk(dn(a.lo - b.hi), up(a.hi - b.lo)); }
static I neg(I a) { return mk(-a.hi, -a.lo); }
static I mul(I a, I b) {
  double p1 = a.lo * b.lo, p2 = a.lo * b.hi, p3 = a.hi * b.lo, p4 = a.hi * b.hi;
  double lo = fmin(fmin(p1, p2), fmin(p3, p4)), hi = fmax(fmax(p1, p2), fmax(p3, p4));
  return mk(dn(lo), up(hi));
}
static I sqr(I a) {
  double l = fabs(a.lo), h = fabs(a.hi);
  if (a.lo <= 0 && a.hi >= 0) return mk(0.0, up(fmax(l, h) * fmax(l, h)));
  double m = fmin(l, h), M = fmax(l, h);
  return mk(dn(m * m), up(M * M));
}
static I dvd(I a, I b) {  /* b must not contain 0 */
  if (!(b.lo > 0 || b.hi < 0)) { fprintf(stderr, "division by an interval containing 0\n"); exit(2); }
  double p1 = a.lo / b.lo, p2 = a.lo / b.hi, p3 = a.hi / b.lo, p4 = a.hi / b.hi;
  return mk(dn(fmin(fmin(p1, p2), fmin(p3, p4))), up(fmax(fmax(p1, p2), fmax(p3, p4))));
}
static I ilog(I a) {
  if (!(a.lo > 0)) { fprintf(stderr, "log of a non-positive interval\n"); exit(2); }
  return mk(dn4(log(a.lo)), up4(log(a.hi)));
}
static I isinh(I a) { return mk(dn4(sinh(a.lo)), up4(sinh(a.hi))); }   /* increasing */
/* sin^2 on [a.lo, a.hi], |a| <= 64: exact for monotone pieces; conservative ([0, ...] or [..., 1]) near extrema */
static I PI_;
static I isin2(I a) {
  /* sin^2 on [a.lo, a.hi] with a.hi - a.lo < 1: its extrema are at multiples k pi/2; only the multiple nearest the
     midpoint can lie inside (the others are at least pi/4 - 1/2 away from the interval) */
  if (!(a.hi - a.lo < 1.0)) return mk(0.0, 1.0);
  double mid = 0.5 * (a.lo + a.hi);
  long k = (long)floor(mid / (0.5 * PI_.lo) + 0.5);
  double e1 = (double)k * (0.5 * PI_.lo), e2 = (double)k * (0.5 * PI_.hi);
  double elo = dn(dn(fmin(e1, e2))), ehi = up(up(fmax(e1, e2)));
  double s1 = sin(a.lo), s2 = sin(a.hi);
  double v1 = s1 * s1, v2 = s2 * s2;
  double lo = fmin(v1, v2), hi = fmax(v1, v2);
  if (ehi >= a.lo && elo <= a.hi) {           /* an extremum may lie inside */
    if (k % 2 == 0) lo = 0.0; else hi = 1.0;
  }
  /* libm sin within 4 ulp: sin^2 within 2 * 4 ulp relative; widen by 4e-15 relative */
  lo = dn(lo * (1.0 - 4e-15)); hi = up(hi * (1.0 + 4e-15));
  if (lo < 0) lo = 0;
  if (hi > 1) hi = 1;
  return mk(lo, hi);
}
/* shc(s) = sinh(s)/s and sinc(s) = sin(s)/s for |s| <= 2 on an interval: even series in s^2, with remainder */
static I ishc(I s, int sinc) {
  I s2 = sqr(s);                       /* [0, S] */
  /* sum_{k=0}^{K} (+-1)^k s^{2k}/(2k+1)!, monotone pieces handled by evaluating at both ends of s2 (each term is
     monotone in s2 >= 0; the alternating sum is enclosed by summing per-term enclosures) */
  I sum = pt(1.0), term = pt(1.0);
  int K = 18;
  for (int k = 1; k <= K; k++) {
    term = dvd(mul(term, s2), pt((double)((2 * k) * (2 * k + 1))));
    sum = (sinc && (k % 2 == 1)) ? sub(sum, term) : add(sum, term);
  }
  /* remainder: |sum_{k>K}| <= S^{K+1}/(2K+3)! * 2 for S <= 4 */
  double S = s2.hi, r = 1.0;
  for (int k = 0; k < K + 1; k++) r *= S;
  for (int k = 1; k <= 2 * K + 3; k++) r /= k;
  r = up(2 * r) ;
  return mk(dn(sum.lo - r), up(sum.hi + r));
}
static I icos_m(long m, long fourN) {   /* cos(m pi / (2N)) with fourN = 4N, m reduced mod 4N */
  m %= fourN; if (m < 0) m += fourN;
  double twoN = (double)(fourN / 2);
  I ang = dvd(mul(pt((double)m), PI_), pt(twoN));    /* in [0, 2 pi) */
  /* cos on [lo, hi] of tiny width: Lipschitz 1 */
  double mid = 0.5 * (ang.lo + ang.hi);
  double c = cos(mid), w = up(fmax(mid - ang.lo, ang.hi - mid));
  return mk(dn(dn4(c) - w), up(up4(c) + w));
}

int main(void) {
  long M, J;
  if (scanf("%ld %ld", &M, &J) != 2) return 1;
  double a, b;
  scanf("%la %la", &a, &b); PI_ = mk(a, b);
  long ii; scanf("%ld", &ii);           /* index of the representative slit i */
  double cf[5][64];
  for (int c = 0; c < 5; c++) for (int n = 0; n < 64; n++) scanf("%la", &cf[c][n]);
  int cls[64]; long N[64]; I A[64], B1[64], B2[64], TH[64];
  for (long j = 0; j < J; j++) {
    double l1, h1, l2, h2, l3, h3, l4, h4;
    scanf("%d %la %la %la %la %la %la %la %la %ld", &cls[j], &l1, &h1, &l2, &h2, &l3, &h3, &l4, &h4, &N[j]);
    A[j] = mk(l1, h1); B1[j] = mk(l2, h2); B2[j] = mk(l3, h3); TH[j] = mk(l4, h4);
  }
  int ci = cls[ii];
  /* per slit j: the nodes y_l, psi_j(y_l), and the y-parts of U and V */
  I *Y[64], *PSI[64], *UY[64], *VY[64];
  for (long j = 0; j < J; j++) {
    long n_ = N[j];
    Y[j] = malloc(n_ * sizeof(I)); PSI[j] = malloc(n_ * sizeof(I)); UY[j] = malloc(n_ * sizeof(I)); VY[j] = malloc(n_ * sizeof(I));
    for (long l = 1; l <= n_; l++) {
      I y = icos_m(2 * l - 1, 4 * n_);
      I ps = pt(0.0);
      for (int n = 1; n < 64; n += 2) ps = add(ps, mul(pt(cf[cls[j]][n]), icos_m((long)n * (2 * l - 1), 4 * n_)));
      I y2 = sqr(y);
      Y[j][l - 1] = y; PSI[j][l - 1] = ps;
      UY[j][l - 1] = mul(A[j], y);
      VY[j][l - 1] = add(add(TH[j], mul(B1[j], y2)), mul(B2[j], sqr(y2)));
    }
  }
  for (long k = M / 2; k < M; k++) {          /* the zeros x_k = cos((2k+1) pi/(2M)) < 0 */
    I x = icos_m(2 * k + 1, 4 * M);
    I x2 = sqr(x);
    I ux = mul(A[ii], x);
    I vx = add(add(TH[ii], mul(B1[ii], x2)), mul(B2[ii], sqr(x2)));
    /* -a_i x + sum_n c_{i,n} (pi/n) T_n(x) */
    I f = neg(ux);
    for (int n = 1; n < 64; n += 2)
      f = add(f, dvd(mul(mul(pt(cf[ci][n]), PI_), icos_m((long)n * (2 * k + 1), 4 * M)), pt((double)n)));
    for (long j = 0; j < J; j++) {
      I acc = pt(0.0);
      for (long l = 0; l < N[j]; l++) {
        I g2;
        if (j == ii) {
          I y = Y[j][l];
          I tau = dvd(sub(x, y), pt(2.0));
          I xy = add(x, y);
          I cub = add(add(mul(x2, x), mul(x2, y)), add(mul(x, sqr(y)), mul(sqr(y), y)));
          I beta = add(mul(B1[j], xy), mul(B2[j], cub));
          I alpha = A[j];
          I s1 = ishc(mul(tau, alpha), 0), s2 = ishc(mul(tau, beta), 1);
          g2 = add(mul(sqr(alpha), sqr(s1)), mul(sqr(beta), sqr(s2)));
        } else {
          I U = sub(ux, UY[j][l]);
          I V = sub(vx, VY[j][l]);
          I sh = isinh(dvd(U, pt(2.0)));
          g2 = mul(pt(4.0), add(sqr(sh), isin2(dvd(V, pt(2.0)))));
        }
        I lg = dvd(ilog(g2), pt(2.0));                 /* log|G~| */
        acc = add(acc, mul(lg, PSI[j][l]));
      }
      /* minus (pi/N) acc */
      f = sub(f, dvd(mul(PI_, acc), pt((double)N[j])));
    }
    printf("%ld %a %a\n", k, f.lo, f.hi);
  }
  return 0;
}
