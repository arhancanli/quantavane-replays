#!/usr/bin/env python3
"""10a [P2026] K_G <= 1.7813*: program A. A replay of the computational lemma of Peng, arXiv 2609.20074v1
(Lemma 5.1: b_1 > 0.881850816, sum_{3<=n<=301} |b_n| < 3e-8, B_3 < 107) and of the closing inequality, from the
paper's text and its exact parameters (Appendix A) only. No program of the source was read or run.

Notation (paper, Sections 2-5). h_j = He_j/sqrt(j!) (orthonormal for the Gaussian measure gamma);
P(x) = sum alpha_j h_j(x) (j = 1, 3, ..., 11); rho(t) = sum c_j t^j (j = 1, ..., 51 odd); f(w, x) = sgn(w + P(x));
A_ab = E[f(W, X) h_a(W) h_b(X)] = E[h_b(X) q_a(P(X))], q_0(z) = erf(z/sqrt 2), q_a(z) = 2 phi(z) h_(a-1)(-z)/sqrt a;
H(t) = (pi/2) sum A_ab^2 rho(t)^a (-t)^b = sum b_n t^n; gamma = 5000 pi/17813.
The criterion (Lemma 2.2 and Proposition 2.3): gamma + sum_{n>=3} |b_n| < b_1 implies K_G <= pi/(2 gamma) = 1.7813.
The tail (Proposition 4.2): sum_{n>301} |b_n| <= B_3 / sqrt(10 * 301^5), B_3^2 <= C_a * I with
C_a = 1/2 + pi/(2 sqrt a) coth(pi/sqrt a), a = 1/1000, and I the integral of the integrand of Section 5.2 over x, y in R, s in [0, 1].

Part "head" (this program): A_ab for all a + b <= 301 with a + b odd (22,952 pairs) at once. A_ab = 2 int_0^oo of
an even integrand. On [0, 11/2]: composite Gauss-Legendre (rigorous nodes and weights from Arb), 24 nodes per panel;
the quadrature error on each panel bounded by Trefethen's theorem (ATAP, Thm 19.3): h (64/15) M rho^(-2n)/(rho^2-1),
with M bounded on a complex box containing the Bernstein ellipse rho = 4 (separately for |h_b phi| and |q_a(P)|).
On [11/2, oo): P is increasing with P >= 20, and by Cramer's inequality |He_n(x)| <= 1.0865 sqrt(n!) e^(x^2/4)
(Abramowitz-Stegun 22.14.17), |h_b phi q_a(P)| <= (2 (1.0865)^2/(2 pi sqrt a)) e^(-x^2/4) e^(-P^2/4) for a >= 1;
for a = 0, int_c^oo h_b phi = h_(b-1)(c) phi(c)/sqrt b exactly and erfc(P/sqrt 2) <= e^(-P^2/2).
Then b_n (n <= 301) by ball polynomial arithmetic.

Part "b3": the integral I of the integrand displayed in Section 5.2 by an adaptive partition of [0,8] x [-8,8] x [0,1] (times 2 by the reflection
(x, y) -> (-x, -y)) into boxes, each bounded above in ball arithmetic exactly as specified in Section 5.2 (the F_ij
recurrence, the Bell polynomials, the coefficients kappa_pa of Table 1, the quadratic form split
(A-B)(u+v)^2/2 + (A+B)(u-v)^2/2 with both coefficients at least 1/(1 + R) > 1/2); the exterior |x| > 8 or |y| > 8
by an explicit bound. The partition is written out for program B.

Part "final": gamma + head + B_3/sqrt(10 301^5) < b_1.
Requires python-flint (Arb); see run_all.sh for the pinned version. Usage: p26_a.py head|b3|final [options]
"""
import sys, json, os, time, heapq
from fractions import Fraction as Q
from flint import arb, acb, arb_poly, ctx

ALPHA = {1: "0.0033816437099488204", 3: "0.1158676765085016", 5: "0.0095104711257198767",
         7: "-0.0080278072994477648", 9: "0.0058210183960585591", 11: "0.012821610768478974"}
CJ = {1: "0.89128717292541904", 3: "-0.10638551869595701", 5: "0.0022804376833522166",
      7: "-3.7814263579902146e-8", 9: "-1.0300867559781962e-5", 11: "5.3495808526173962e-10",
      13: "8.3881520203661931e-6", 15: "-1.2307543439537263e-6", 17: "5.3773250452303132e-7",
      19: "1.6847785460476277e-6", 21: "3.7811663671652489e-7", 23: "3.2246923243582631e-8",
      25: "1.8981462348684652e-7", 27: "6.4099421934521302e-8", 29: "-5.1544356066976797e-8",
      31: "6.9279841157102393e-9", 33: "5.553232634844159e-8", 35: "1.3619423272928659e-8",
      37: "-3.1796367719425194e-8", 39: "-2.2737533175642567e-8", 41: "8.4604728136655195e-9",
      43: "1.7763396879945166e-8", 45: "4.1161558338437501e-9", 47: "-8.7507870739030589e-9",
      49: "-7.9288491194999501e-9", 51: "4.9151067348141944e-10"}
N = 301

def qd(s):
    """an exact decimal string as an exact rational"""
    return Q(s)

def arbq(q):
    q = Q(q)
    return arb(q.numerator) / arb(q.denominator)

def hermite_vals(x, n):
    """[h_0(x), ..., h_n(x)] by h_(j+1) = (x h_j - sqrt j h_(j-1))/sqrt(j+1) (x an arb or acb)"""
    h = [x * 0 + 1, x]
    for j in range(1, n):
        h.append((x * h[j] - SQ[j] * h[j - 1]) / SQ[j + 1])
    return h[:n + 1]

def setup(prec):
    global SQ, ALP, PI, SQ2PI, IS2
    ctx.prec = prec
    SQ = [arb(j).sqrt() for j in range(N + 3)]
    ALP = {j: arbq(qd(s)) for j, s in ALPHA.items()}
    PI = arb.pi(); SQ2PI = (2 * PI).sqrt(); IS2 = 1 / arb(2).sqrt()

def P_of(x):
    h = hermite_vals(x, 11)
    return sum((ALP[j] * h[j] for j in ALP), x * 0)

def phi(z):
    return (-(z * z) / 2).exp() / SQ2PI

def q_vals(p, amax):
    """[q_0(p), ..., q_amax(p)]"""
    out = [(p * IS2).erf()]
    hm = hermite_vals(-p, amax)
    ph = phi(p)
    for a in range(1, amax + 1):
        out.append(2 * ph * hm[a - 1] / SQ[a])
    return out

def exact(x):
    """an exact dyadic number (radius 0 arb) as the string 'm*2^e'"""
    m, e = x.man_exp()
    return "%d*2^%d" % (int(m), int(e))

def from_exact(s):
    m, e = s.split("*2^")
    return Q(int(m)) * Q(2) ** int(e)

def mag(z):
    """an upper bound for |z| (arb or acb) as an arb"""
    return arb(abs(z).upper())

# ------------------------------------------------------------------ the head
PANELS = [(Q(i, 4), Q(i + 1, 4)) for i in range(18)] + [(Q(45 + i, 10), Q(46 + i, 10)) for i in range(5)] \
    + [(Q(200 + i, 40), Q(201 + i, 40)) for i in range(20)]           # [0, 4.5], [4.5, 5], [5, 5.5]
XC = Q(11, 2)
NGL = 24
RHO = 4

def head(out_path):
    setup(160)
    t0 = time.time()
    # P increasing with P >= 20 on [11/2, oo): P'(x) > 0 on [11/2, 16] (interval check on 1/8-cells) and beyond 16 by
    # the leading term; P(11/2) >= 20.
    xs = arb_poly([0, 1])
    Hp = [arb_poly([1]), xs]
    for j in range(1, 11):
        Hp.append((xs * Hp[j] - arb_poly([SQ[j]]) * Hp[j - 1]) * arb_poly([1 / SQ[j + 1]]))
    Ppoly = sum((Hp[j] * arb_poly([ALP[j]]) for j in ALP), arb_poly([0]))
    dP = Ppoly.derivative()
    ok_mono = all(dP(arb(44 + k) / 8 + arb(0, 0.0625)) > 0 for k in range(0, 85))
    pc = dP.coeffs()
    lead = pc[10]; rest = sum((abs(pc[j]) * arb(16) ** (j - 10) for j in range(10)), arb(0))
    ok_far = bool(lead - rest > 0)            # P'(x) >= x^10 (p'_10 - sum |p'_j| 16^(j-10)) > 0 for x >= 16
    Pxc = Ppoly(arbq(XC))
    ok_val = bool(Pxc > 20)
    print("P increasing on [11/2, oo): %s; P(11/2) = %s > 20: %s" % (ok_mono and ok_far, Pxc.str(8), ok_val))
    if not (ok_mono and ok_far and ok_val):
        raise SystemExit("head: the tail argument fails")
    pairs = [(a, b) for a in range(N + 1) for b in range(N + 1 - a) if (a + b) % 2 == 1]
    A = [[arb(0)] * (N + 1 - a) for a in range(N + 1)]
    Err = [[arb(0)] * (N + 1 - a) for a in range(N + 1)]
    gl = [arb.legendre_p_root(NGL, k, weight=True) for k in range(NGL)]
    rr = arb(RHO)
    efac = (arb(64) / 15) * rr ** (-2 * NGL) / (rr * rr - 1)
    for (lo, hi) in PANELS:
        c = arbq((lo + hi) / 2); hw = arbq((hi - lo) / 2)
        for (xk, wk) in gl:
            x = c + hw * xk
            H = hermite_vals(x, N)
            wphi = wk * hw * phi(x)
            Hw = [hb * wphi for hb in H]
            Qv = q_vals(P_of(x), N)
            for a in range(N + 1):
                qa = Qv[a]; row = A[a]
                for b in range(1 - a % 2, N + 1 - a, 2):
                    row[b] += qa * Hw[b]
        # bounds on a box containing the Bernstein ellipse: real half-width hw (rho + 1/rho)/2, imaginary hw (rho - 1/rho)/2.
        # |h_n(x+iy)| <= K e^(x^2/4) e^(|y| sqrt n) (Cramer, K = 1.0865, with He_n(x+iy) = sum C(n,k) He_k(x) (iy)^(n-k));
        # |phi(x+iy)| = e^(-(x^2-y^2)/2)/sqrt(2 pi); P on the box lies in the disk P(c) +- sum_k |P^(k)(c)/k!| r^k;
        # |q_a(s+i tau)| <= (2K/sqrt(2 pi a)) e^(-s^2/4 + tau^2/2 + |tau| sqrt(a-1)) for a >= 1; erf on the disk for a = 0.
        rre = hw * (rr + 1 / rr) / 2; rim = hw * (rr - 1 / rr) / 2
        K = arb("1.0865")
        xmax = c + rre; xmin = c - rre
        xmin2 = xmin * xmin if xmin > 0 else arb(0)
        ymax = rim
        hbase = K * (xmax * xmax / 4).exp() * (-(xmin2 - ymax * ymax) / 2).exp() / SQ2PI
        HB = [arb((hbase * (ymax * SQ[b]).exp()).upper()) for b in range(N + 1)]
        shifted = Ppoly(arb_poly([c, 1]))          # P(c + w) as a polynomial in w
        sc = shifted.coeffs()
        r = (rre * rre + rim * rim).sqrt()
        rp = sum((abs(sc[k]) * r ** k for k in range(1, len(sc))), arb(0))
        pc = sc[0]
        smin = pc - rp; smax = pc + rp
        smin2 = arb(0) if (smin < 0 < smax or not (smin > 0 or smax < 0)) else (smin * smin if smin > 0 else smax * smax)
        tau = rp
        disk = acb(arb(pc, rp.upper()), arb(0, rp.upper()))
        QA = [arb(abs((disk * IS2).erf()).upper())]
        qbase = 2 * K / SQ2PI * (-smin2 / 4 + tau * tau / 2).exp()
        QA += [arb((qbase / SQ[a] * (tau * SQ[a - 1]).exp()).upper()) for a in range(1, N + 1)]
        cf = hw * efac
        for a in range(N + 1):
            qa = QA[a] * cf; row = Err[a]
            for b in range(1 - a % 2, N + 1 - a, 2):
                row[b] += qa * HB[b]
    print("Gauss-Legendre on %d panels done (%.0f s)" % (len(PANELS), time.time() - t0))
    # tails on [11/2, oo)
    xc = arbq(XC); phic = phi(xc); Hc = hermite_vals(xc, N)
    tail_a = 2 * arb("1.0865") ** 2 / (2 * PI) * PI.sqrt() * (-(arb(20) ** 2) / 4).exp()   # / sqrt a, a >= 1
    tail_0 = arb("1.0865") / SQ2PI * PI.sqrt() * (-(arb(20) ** 2) / 2).exp()
    maxerr = 0.0; worst = None
    Afin = {}
    for a in range(N + 1):
        for b in range(1 - a % 2, N + 1 - a, 2):
            v = A[a][b]
            if a == 0:
                v = v + Hc[b - 1] * phic / SQ[b]
                e = Err[a][b] + tail_0
            else:
                e = Err[a][b] + tail_a / SQ[a]
            v = 2 * (v + arb(0, e.upper()))
            Afin[(a, b)] = v
            if float(e.upper()) > maxerr:
                maxerr = float(e.upper()); worst = (a, b)
    print("pairs: %d; largest quadrature + tail remainder bound %.3e at (a, b) = %s" % (len(Afin), maxerr, worst))
    print("A_10 = %s, A_01 = %s" % (Afin[(1, 0)].str(25), Afin[(0, 1)].str(25)))
    # Parseval check (information): sum of A_ab^2 over the computed pairs is below 1
    s2 = sum((v * v for v in Afin.values()), arb(0))
    print("information: sum of A_ab^2 over a + b <= 301 = %s (Parseval: at most 1)" % s2.str(12))
    # b_n for n <= 301
    rho = arb_poly([0] + [arbq(qd(CJ[j])) if j in CJ else 0 for j in range(1, 52)])
    Rsum = sum((abs(arbq(qd(s))) for s in CJ.values()), arb(0))
    tot = arb_poly([0])
    ra = arb_poly([1])
    for a in range(N + 1):
        S = [arb(0)] * (N + 1 - a)
        for b in range(1 - a % 2, N + 1 - a, 2):
            S[b] = Afin[(a, b)] ** 2 * (-1) ** b
        tot = tot + ra * arb_poly(S)
        ra = arb_poly((ra * rho).coeffs()[:N + 1])
    bc = [(c * PI / 2) for c in tot.coeffs()[:N + 1]]
    b1 = bc[1]
    headsum = sum((abs(bc[n]) for n in range(3, N + 1) if n < len(bc)), arb(0))
    even = max((abs(bc[n]).upper() for n in range(0, N + 1, 2) if n < len(bc)), default=0)
    print("R = sum |c_j| = %s (< 1)" % Rsum.str(30))
    print("b_1 in %s" % b1.str(25))
    print("sum_{3<=n<=301} |b_n| <= %s" % arb(headsum.upper()).str(10))
    print("information: largest |b_n| over even n <= %s (zero by parity)" % arb(even).str(3))
    print("head time %.0f s" % (time.time() - t0))
    json.dump({"bn_mid": [float(c.mid()) for c in bc], "b1_lo": exact(b1.lower()), "b1": b1.str(30),
               "headsum_hi": exact(headsum.upper()), "headsum": arb(headsum.upper()).str(10), "R": Rsum.str(30)},
              open(out_path, "w"), indent=1)

# ------------------------------------------------------------------ the B_3 integrand (Section 5.2)
AW = Q(1, 1000)                                   # the trace weight a

def setup_b3(prec):
    global PP, DPP, BELL, CJA, RSUM, DELTA, DAMP, AWA
    setup(prec)
    xs = arb_poly([0, 1])
    Hp = [arb_poly([1]), xs]
    for j in range(1, 11):
        Hp.append((xs * Hp[j] - arb_poly([SQ[j]]) * Hp[j - 1]) * arb_poly([1 / SQ[j + 1]]))
    PP = sum((Hp[j] * arb_poly([ALP[j]]) for j in ALP), arb_poly([0]))
    DPP = PP.derivative()
    # Bell polynomials: d^k/dx^k g(P(x)) = sum_i B_(k,i)(x) g^(i)(P(x)); B_(0,0) = 1, B_(k+1,i) = B_(k,i)' + P' B_(k,i-1)
    BELL = [[arb_poly([1])]]
    for k in range(4):
        row = []
        for i in range(k + 2):
            t = arb_poly([0])
            if i <= k:
                t = t + BELL[k][i].derivative()
            if i >= 1:
                t = t + DPP * BELL[k][i - 1]
            row.append(t)
        BELL.append(row)
    CJA = {j: arbq(qd(v)) for j, v in CJ.items()}
    RSUM = sum((abs(v) for v in CJA.values()), arb(0))
    DELTA = 1 - RSUM * RSUM                        # |d| >= 1 - |r|^2 >= DELTA on the circle
    DAMP = 1 / (1 + RSUM)                          # Re 1/(1 +- r) >= 1/(1 + |r|) >= DAMP
    AWA = arbq(AW)

def rho_q(t):
    """rho(t), q1 = D rho, q2 = D^2 rho, q3 = D^3 rho (D = t d/dt), by direct sums over the odd powers"""
    t2 = t * t
    tp = t
    r = q1 = q2 = q3 = t * 0
    for j in range(1, 52, 2):
        cj = CJA[j] * tp
        r += cj; q1 += j * cj; q2 += (j * j) * cj; q3 += (j * j * j) * cj
        tp = tp * t2
    return r, q1, q2, q3

KSIGN = 1                                         # red control: the sign of kappa_(1,1) flipped

def kappa(t, q1, q2, q3):
    """Table 1: (p, a) -> kappa_pa"""
    return [(1, 0, q3), (2, 0, 3 * q1 * q2), (3, 0, q1 ** 3), (0, 1, -t), (0, 2, 3 * t * t), (0, 3, -(t ** 3)),
            (1, 1, -3 * KSIGN * t * (q2 + q1)), (2, 1, -3 * t * q1 * q1), (1, 2, 3 * t * t * q1)]

def F_table(u, v, r, dinv):
    lu = (-u + r * v) * dinv; lv = (-v + r * u) * dinv
    F = [[None] * 4 for _ in range(4)]
    F[0][0] = lu * 0 + 1
    for j in range(1, 4):
        F[0][j] = lv * F[0][j - 1] - ((j - 1) * dinv * F[0][j - 2] if j >= 2 else 0)
    for i in range(1, 4):
        for j in range(0, 4):
            val = lu * F[i - 1][j]
            if i >= 2:
                val -= (i - 1) * dinv * F[i - 2][j]
            if j >= 1:
                val += j * r * dinv * F[i - 1][j - 1]
            F[i][j] = val
    return F

def R_pref(Bx, By, F, kap, e):
    """R_e of Section 5.2"""
    tot = 0
    for (p, a, k) in kap:
        m = a + e
        acc = 0
        for i in range(m + 1):
            bi = Bx[m][i]
            if p == 0 and i == 0:
                continue
            for j in range(m + 1):
                if p == 0 and j == 0:
                    continue
                bj = By[m][j]
                if (m > 0 and (i == 0 or j == 0)):
                    continue                     # B_(m,0) = 0 for m >= 1
                g = F[i + p - 1][j + p - 1] if p >= 1 else F[i - 1][j - 1]
                acc = acc + bi * bj * g
        tot = tot + k * acc
    return tot

XMONO = arb(11) / 2                    # P is increasing on [11/2, oo) (checked in "head") and odd

def P_box(X):
    """an enclosure of P on the interval X: exact by monotonicity on |x| >= 11/2, else the tighter of Horner and
    the mean-value form P(m) + P'(X)(X - m)"""
    lo = arb(X.lower()); hi = arb(X.upper())
    if lo >= XMONO or hi <= -XMONO:
        a = PP(lo); b = PP(hi)
        return arb((a.lower() + b.upper()) / 2, ((b.upper() - a.lower()) / 2).upper()) if True else None
    return derivs_box(X)[0]

def derivs_box(X):
    """[P(X), and its derivatives of order 1 to 4 on X]: the exact Taylor shift s_k at c = mid(X), plus the remainder
    sum_{k>m} C(k,m) |s_k| r^(k-m) on |w| <= r"""
    c = arb(X.mid()); r = arb(X.rad())
    sc = PP(arb_poly([c, 1])).coeffs()
    out = []
    for m in range(5):
        cen = sc[m] if m < len(sc) else arb(0)
        tail = sum((comb_int(k, m) * abs(sc[k]) * r ** (k - m) for k in range(m + 1, len(sc))), arb(0))
        out.append(fact_int(m) * (cen + arb(0, tail.upper())))
    return out

def comb_int(n, k):
    from math import comb
    return comb(n, k)

def fact_int(n):
    from math import factorial
    return factorial(n)

def bell_at(x):
    """Faa di Bruno with the derivative enclosures: B_(k,i) for k <= 4"""
    p0, p1, p2, p3, p4 = derivs_box(x)
    one = arb(1); z = arb(0)
    return [[one], [z, p1], [z, p2, p1 * p1], [z, p3, 3 * p1 * p2, p1 ** 3],
            [z, p4, 4 * p1 * p3 + 3 * p2 * p2, 6 * p1 * p1 * p2, p1 ** 4]]

def sq_lower(I):
    """an exact lower bound for min of z^2 over the interval I (no dependency loss)"""
    lo = arb(I.lower().lower()); hi = arb(I.upper().upper())
    if lo > 0:
        return arb((lo * lo).lower())
    if hi < 0:
        return arb((hi * hi).lower())
    return arb(0)

def integrand_upper(X, Y, S):
    """an upper bound of the integrand (Section 5.2) on the box X x Y x S (arb intervals), s = theta/pi"""
    t = acb(0, PI * S).exp()
    r, q1, q2, q3 = rho_q(t)
    u = P_box(X); v = P_box(Y)
    d = 1 - r * r
    dinv = 1 / d
    F = F_table(u, v, r, dinv)
    kap = kappa(t, q1, q2, q3)
    Bx = bell_at(X); By = bell_at(Y)
    R0 = R_pref(Bx, By, F, kap, 0)
    R1 = R_pref(Bx, By, F, kap, 1)
    num = arb((abs(R0) ** 2 + AWA * abs(R1) ** 2).upper())
    dl = abs(d).lower()
    den = 2 * PI * (dl if dl > DELTA else DELTA)
    amb = (1 / (1 + r)).real.lower(); apb = (1 / (1 - r)).real.lower()
    amb = amb if amb > DAMP else DAMP; apb = apb if apb > DAMP else DAMP
    sp = u + v; sm = u - v
    sp2 = sq_lower(sp); sm2 = sq_lower(sm)
    x2 = sq_lower(X); y2 = sq_lower(Y)
    # two valid lower bounds of A(u^2+v^2) - 2Buv: the split form, and min(A-B, A+B)(u^2 + v^2) (as 2|uv| <= u^2+v^2)
    e1 = amb * sp2 / 2 + apb * sm2 / 2
    e2 = (amb if amb < apb else apb) * (sq_lower(u) + sq_lower(v))
    ql = arb(e1.lower()); q2 = arb(e2.lower())
    expo = (ql if ql > q2 else q2) + (x2 + y2) / 2
    return arb((num / den * (-expo).exp()).upper())

def integrand_point(x, y, s):
    """the integrand at a point (for the self-test and estimates): a ball"""
    return integrand_upper(arb(x), arb(y), arb(s))


# ------------------------------------------------------------------ the partition
XU, YU, SU = 16, 16, 24              # coordinates in units of 2^-16 (x, y) and 2^-24 (s)

def box_arb(b):
    x0, x1, y0, y1, s0, s1 = b
    X = arb(x0 + x1) / 2 ** (XU + 1); X = arb(X, (x1 - x0) / 2 ** (XU + 1))
    Y = arb(y0 + y1) / 2 ** (YU + 1); Y = arb(Y, (y1 - y0) / 2 ** (YU + 1))
    S = arb(s0 + s1) / 2 ** (SU + 1); S = arb(S, (s1 - s0) / 2 ** (SU + 1))
    return X, Y, S

def vol(b):
    x0, x1, y0, y1, s0, s1 = b
    return Q((x1 - x0) * (y1 - y0) * (s1 - s0), 2 ** (XU + YU + SU))

def contrib(b):
    """the box's upper bound times its volume times 2 (the reflection (x, y) -> (-x, -y)), as a float rounded up"""
    from math import nextafter, inf
    X, Y, S = box_arb(b)
    U = integrand_upper(X, Y, S)
    c = U * arbq(vol(b)) * 2
    if not c.is_finite():
        return inf, inf
    f = nextafter(float(c.upper()), inf)
    return f, f

def split(b, dim):
    x0, x1, y0, y1, s0, s1 = b
    if dim == 0:
        m = (x0 + x1) // 2; return (x0, m, y0, y1, s0, s1), (m, x1, y0, y1, s0, s1)
    if dim == 1:
        m = (y0 + y1) // 2; return (x0, x1, y0, m, s0, s1), (x0, x1, m, y1, s0, s1)
    m = (s0 + s1) // 2; return (x0, x1, y0, y1, s0, m), (x0, x1, y0, y1, m, s1)

def exterior_bound():
    """|x| > 8 or |y| > 8, every s: an explicit bound (see the module docstring and RESULT.md)"""
    pc = PP.coeffs()
    c8 = pc[11] - sum((abs(pc[j]) * arb(8) ** (j - 11) for j in range(11)), arb(0))
    Phat1 = sum((abs(c) for c in pc), arb(0))                # |P(x)| <= Phat1 * max(1,|x|)^11
    D = 1 / DELTA
    L = 2 * Phat1 * D                                        # |l_u|, |l_v| <= L * (max(1,|x|) max(1,|y|))^11
    Fh = [[arb(0)] * 4 for _ in range(4)]
    Fh[0][0] = arb(1)
    for j in range(1, 4):
        Fh[0][j] = L * Fh[0][j - 1] + ((j - 1) * D * Fh[0][j - 2] if j >= 2 else 0)
    for i in range(1, 4):
        for j in range(4):
            v = L * Fh[i - 1][j]
            if i >= 2: v += (i - 1) * D * Fh[i - 2][j]
            if j >= 1: v += j * D * Fh[i - 1][j - 1]
            Fh[i][j] = v
    qh = [sum((arb(n) ** k * abs(c) for n, c in CJA.items()), arb(0)) for k in range(4)]
    kh = [(1, 0, qh[3]), (2, 0, 3 * qh[1] * qh[2]), (3, 0, qh[1] ** 3), (0, 1, arb(1)), (0, 2, arb(3)), (0, 3, arb(1)),
          (1, 1, 3 * (qh[2] + qh[1])), (2, 1, 3 * qh[1] ** 2), (1, 2, 3 * qh[1])]
    Bh = [[sum((abs(c) for c in BELL[k][i].coeffs()), arb(0)) for i in range(k + 1)] for k in range(5)]
    def Rh(e):
        tot = arb(0)
        for (p, a, k) in kh:
            m = a + e; acc = arb(0)
            for i in range(m + 1):
                for j in range(m + 1):
                    if (p == 0 and (i == 0 or j == 0)) or (m > 0 and (i == 0 or j == 0)):
                        continue
                    g = Fh[i + p - 1][j + p - 1] if p >= 1 else Fh[i - 1][j - 1]
                    acc += Bh[m][i] * Bh[m][j] * g
            tot += k * acc
        return tot
    C0, C1 = Rh(0), Rh(1)
    deg = 40 + 66                                            # Bell degree <= 40, F degree <= 6 * 11 in each variable
    K = 2 * deg                                              # the square
    mom = arb(2) ** ((arb(K) - 1) / 2) * ((arb(K) + 1) / 2).gamma()   # int_0^oo x^K e^(-x^2/2) dx
    Mx = mom; My = 2 + 2 * mom
    u8 = c8 * arb(8) ** 11
    ext = 4 * (C0 * C0 + AWA * C1 * C1) / (2 * PI * DELTA) * (-(u8 * u8) / 4).exp() * Mx * My
    return c8, u8, C0, C1, ext

def b3(out_json, out_boxes, target, max_splits, max_seconds):
    import gzip
    setup_b3(96)
    t0 = time.time()
    c8, u8, C0, C1, ext = exterior_bound()
    print("exterior: c_8 = %s > 0, c_8 8^11 = %s; C_0 <= %s, C_1 <= %s; exterior integral <= %s"
          % (c8.str(10), u8.str(10), C0.str(5), C1.str(5), arb(ext.upper()).str(3)))
    if not (c8 > 0):
        raise SystemExit("b3: c_8 not positive")
    heap = []; total = 0.0; k = 0
    for i in range(8):
        for j in range(16):
            for m in range(64):
                b = (i << XU, (i + 1) << XU, (j - 8) << YU, (j - 7) << YU, m << (SU - 6), (m + 1) << (SU - 6))
                f, c = contrib(b)
                heapq.heappush(heap, (-f, k, b, c)); k += 1; total += f
    from math import fsum
    total = fsum(-h[0] for h in heap)
    print("initial 8 x 16 x 64 boxes: upper sum %.6g" % total)
    splits = 0
    from math import fsum
    while splits < max_splits and time.time() - t0 < max_seconds:
        if total <= target or splits % 1000 == 0:
            total = fsum(-h[0] for h in heap)        # exact re-summation (the running float total cancels badly)
            if total <= target:
                break
        nf, _, b, c = heapq.heappop(heap)
        best = None
        for dim in range(3):
            if (dim == 0 and b[1] - b[0] < 2) or (dim == 1 and b[3] - b[2] < 2) or (dim == 2 and b[5] - b[4] < 2):
                continue
            ch = split(b, dim)
            cs = [contrib(x) for x in ch]
            sm = cs[0][0] + cs[1][0]
            if best is None or sm < best[0]:
                best = (sm, ch, cs)
        if best is None:
            raise SystemExit("b3: a box cannot be split further")
        total += best[0] + nf
        for x, (f, cc) in zip(best[1], best[2]):
            heapq.heappush(heap, (-f, k, x, cc)); k += 1
        splits += 1
        if splits % 20000 == 0:
            print("  %d splits, %d boxes, upper sum %.6g (%.0f s)" % (splits, len(heap), total, time.time() - t0))
    top = sorted(heap)[:8]
    for (nf, _, b, c) in top:
        x0, x1, y0, y1, s0, s1 = b
        print("  top box x [%.4f, %.4f] y [%.4f, %.4f] s [%.6f, %.6f]: %.3e" % (x0 / 2 ** XU, x1 / 2 ** XU, y0 / 2 ** YU,
              y1 / 2 ** YU, s0 / 2 ** SU, s1 / 2 ** SU, -nf))
    tot_q = sum((Q(c) for (_, _, b, c) in heap), Q(0))          # exact sum of the rounded-up floats
    exact_total = arbq(tot_q)
    I_up = arb((exact_total + ext).upper())
    print("partition: %d boxes after %d splits (%.0f s); finite part <= %s; I <= %s"
          % (len(heap), splits, time.time() - t0, arb(exact_total.upper()).str(10), I_up.str(10)))
    with gzip.open(out_boxes, "wt") as fh:
        for (_, _, b, c) in heap:
            fh.write("%d %d %d %d %d %d\n" % b)
    json.dump({"I_hi": exact(I_up), "I": I_up.str(12), "boxes": len(heap), "ext": arb(ext.upper()).str(5)},
              open(out_json, "w"), indent=1)

def phi3_point(tr, x, y, e):
    """d_x^e d_y^e Phi_3(t; x, y) at a real t (for the self-test): d^(-1/2) R_e exp(-(u^2 - 2ruv + v^2)/(2d))"""
    t = acb(tr)
    r, q1, q2, q3 = rho_q(t)
    X = arb(x); Y = arb(y)
    u = PP(X); v = PP(Y)
    d = 1 - r * r; dinv = 1 / d
    F = F_table(u, v, r, dinv)
    Re = R_pref(bell_at(X), bell_at(Y), F, kappa(t, q1, q2, q3), e)
    return (Re * (-(u * u - 2 * r * u * v + v * v) * dinv / 2).exp() / d.sqrt()).real

def gauss_expect(tr, sc, e, hstep=0.1, L=8):
    """E[d^e Phi_3(t; X, Y)] with corr(X, Y) = sc, by the trapezoidal rule in (x, z), Y = sc X + sqrt(1 - sc^2) Z"""
    from math import sqrt as fs, exp as fe, pi as fp
    w = fs(1 - sc * sc)
    n = int(L / hstep)
    tot = 0.0
    for i in range(-n, n + 1):
        x = i * hstep
        for k in range(-n, n + 1):
            z = k * hstep
            g = float(phi3_point(tr, x, sc * x + w * z, e).mid())
            tot += g * fe(-(x * x + z * z) / 2)
    return tot * hstep * hstep / (2 * fp)

def selftest(head_json):
    setup_b3(96)
    bn = json.load(open(head_json))["bn_mid"]
    ok = True
    for tr in (0.5, -0.3):
        lhs = sum(n ** 3 * bn[n] * tr ** n for n in range(len(bn)))
        rhs = gauss_expect(tr, -tr, 0)
        good = abs(lhs - rhs) < 1e-6 * max(1.0, abs(lhs))
        ok &= good
        print("self-test D^3 H(%s) = sum n^3 b_n t^n = %.10f; T_(-t) Phi_3 (Gaussian expectation) = %.10f: %s" % (tr, lhs, rhs, good))
    tr = 0.5; sc = -0.5; ds = 1e-3
    deriv = (gauss_expect(tr, sc + ds, 0) - gauss_expect(tr, sc - ds, 0)) / (2 * ds)
    e1 = gauss_expect(tr, sc, 1)
    good = abs(deriv - e1) < 1e-5 * max(1.0, abs(e1))
    ok &= good
    print("self-test E_s[d_xy Phi_3] = %.8f against d/ds E_s[Phi_3] = %.8f (t = 0.5, s = -0.5): %s" % (e1, deriv, good))
    return ok

def final(head_json, b3_json, claim):
    """gamma + sum_{3<=n<=301} |b_n| + B_3/sqrt(10 * 301^5) < b_1, with B_3^2 <= C_a I"""
    setup(160)
    hj = json.load(open(head_json)); bj = json.load(open(b3_json))
    b1 = arbq(from_exact(hj["b1_lo"])); hs = arbq(from_exact(hj["headsum_hi"])); I = arbq(from_exact(bj["I_hi"]))
    K = arbq(Q(claim))
    gam = PI / (2 * K)                                       # gamma = pi/(2 K); K = 17813/10000 gives 5000 pi/17813
    aw = arbq(AW)
    Ca = arb(1) / 2 + PI / (2 * aw.sqrt()) * (1 / (PI / aw.sqrt()).tanh())
    B3 = (Ca * I).sqrt()
    tail = B3 / (10 * arb(N) ** 5).sqrt()
    lhs = gam + hs + tail
    margin = b1 - lhs
    print("claim K_G <= %s: gamma = pi/(2K) = %s" % (claim, gam.str(15)))
    print("b_1 > %s; sum_{3<=n<=301} |b_n| < %s" % (b1.str(15), hs.str(10)))
    print("C_a = %s (a = 1/1000); I <= %s; B_3 <= sqrt(C_a I) = %s" % (Ca.str(10), I.str(10), arb(B3.upper()).str(8)))
    print("tail sum_{n>301} |b_n| <= B_3/sqrt(10 * 301^5) = %s" % arb(tail.upper()).str(8))
    print("lemma 5.1 as stated: b_1 > 0.881850816: %s; head < 3e-8: %s; B_3 < 107: %s; C_a < 50.173: %s"
          % (b1 > arb("0.881850816"), hs < arb("3e-8"), B3 < 107, Ca < arb("50.173")))
    ok = margin > 0
    print("gamma + head + tail = %s %s b_1: margin %s" % (lhs.str(15), "<" if ok else "is NOT below", margin.str(8)))
    if ok:
        print("RESULT: K_G <= %s (with the analytic lemmas of the paper, not re-proved)" % claim)
        print("VERDICT PASS")
    else:
        print("VERDICT REFUSED: the criterion gamma + sum |b_n| < b_1 is not established")

if __name__ == "__main__":
    part = sys.argv[1]
    if "--c1-shift" in sys.argv:
        k = sys.argv.index("--c1-shift"); sh = Q(sys.argv[k + 1]); del sys.argv[k:k + 2]
        CJ[1] = str(Q(CJ[1]) + sh)
        print("CONTROL: c_1 shifted by %s" % sh)
    if "--kappa-sign" in sys.argv:
        sys.argv.remove("--kappa-sign"); KSIGN = -1
        print("CONTROL: kappa_(1,1) with its sign flipped")
    if part == "head":
        head(sys.argv[2] if len(sys.argv) > 2 else "head_a.json")
    elif part == "b3":
        b3(sys.argv[2], sys.argv[3], float(sys.argv[4]), int(sys.argv[5]), float(sys.argv[6]))
    elif part == "selftest":
        print("VERDICT SELFTEST PASS" if selftest(sys.argv[2]) else "VERDICT REFUSED (self-test fails)")
    elif part == "final":
        final(sys.argv[2], sys.argv[3], sys.argv[4] if len(sys.argv) > 4 else "17813/10000")
