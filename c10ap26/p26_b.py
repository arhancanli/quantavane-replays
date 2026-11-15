#!/usr/bin/env python3
"""10a [P2026] K_G <= 1.7813*: program B, written separately from program A (no code shared).

Same statements as program A (see its docstring); different routes:
  head:  P in the monomial basis from the explicit coefficients of He_j; h_n(x) from Arb's Hermite function
         (He_n(x) = 2^(-n/2) H_n(x/sqrt 2)) at 320 bits, not from the three-term recurrence; composite Gauss-Legendre with
         30 nodes on a different panel layout and the cut at 28/5; the quadrature error by the derivative form of the
         Gauss-Legendre remainder, (b-a)^(2n+1) (n!)^4 / ((2n+1) ((2n)!)^3) |f^(2n)|, with |f^(2n)| <= (2n)! M / R^(2n)
         (Cauchy) and M bounded on a rectangle containing the R-neighbourhood of the panel; b_n by Horner in rho:
         T_b(t) = sum_a A_ab^2 rho(t)^a, H = (pi/2) sum_b (-t)^b T_b.
  b3:    reads program A's partition, checks that it tiles [0,8] x [-8,8] x [0,1] exactly (recursive reconstruction of the
         bisections inside every initial box, and the dyadic volumes), and bounds the integrand on every box with its
         own code: phi_r^(-1) d_u^i d_v^j phi_r by the closed sum
         sum_k C(i,k) C(j,k) k! (r/d)^k G_(i-k)(l_u, -1/d) G_(j-k)(l_v, -1/d),  G_m(L, c) = sum_l m!/(l! (m-2l)!) (c/2)^l L^(m-2l),
         and the Bell coefficients B_(k,i)(x) = (k!/i!) [w^k] (P(x+w) - P(x))^i by series composition; its own exterior bound.
  self-test: the identity D^3 H(t) = T_(-t) Phi_3(t) at real t with this program's own b_n and integrand.
Requires python-flint (Arb). Usage: p26_b.py head OUT | selftest HEADJSON | b3 BOXES.gz OUT | final ...
"""
import sys, json, time, gzip
from fractions import Fraction as Q
from math import comb, factorial
from flint import arb, acb, arb_poly, acb_poly, ctx

PAR = {
    "alpha": {1: "0.0033816437099488204", 3: "0.1158676765085016", 5: "0.0095104711257198767",
              7: "-0.0080278072994477648", 9: "0.0058210183960585591", 11: "0.012821610768478974"},
    "c": {1: "0.89128717292541904", 3: "-0.10638551869595701", 5: "0.0022804376833522166",
          7: "-3.7814263579902146e-8", 9: "-1.0300867559781962e-5", 11: "5.3495808526173962e-10",
          13: "8.3881520203661931e-6", 15: "-1.2307543439537263e-6", 17: "5.3773250452303132e-7",
          19: "1.6847785460476277e-6", 21: "3.7811663671652489e-7", 23: "3.2246923243582631e-8",
          25: "1.8981462348684652e-7", 27: "6.4099421934521302e-8", 29: "-5.1544356066976797e-8",
          31: "6.9279841157102393e-9", 33: "5.553232634844159e-8", 35: "1.3619423272928659e-8",
          37: "-3.1796367719425194e-8", 39: "-2.2737533175642567e-8", 41: "8.4604728136655195e-9",
          43: "1.7763396879945166e-8", 45: "4.1161558338437501e-9", 47: "-8.7507870739030589e-9",
          49: "-7.9288491194999501e-9", 51: "4.9151067348141944e-10"}}
NMAX = 301

def rat(s):
    return Q(s)

def ball(q):
    q = Q(q)
    return arb(q.numerator) / arb(q.denominator)

def init(prec):
    global PREC, PIb, S2PI, PMON, DPMON, CB, KCR
    PREC = prec; ctx.prec = prec
    PIb = arb.pi(); S2PI = (2 * PIb).sqrt()
    # P in the monomial basis: He_j(x) = sum_m (-1)^m j!/(m! (j-2m)! 2^m) x^(j-2m); h_j = He_j/sqrt(j!)
    coef = [arb(0)] * 12
    for j, a in PAR["alpha"].items():
        aj = ball(rat(a)) / arb(factorial(j)).sqrt()
        for m in range(j // 2 + 1):
            c = Q((-1) ** m * factorial(j), factorial(m) * factorial(j - 2 * m) * 2 ** m)
            coef[j - 2 * m] = coef[j - 2 * m] + aj * ball(c)
    PMON = coef
    DPMON = [coef[k] * k for k in range(1, 12)]
    CB = {j: ball(rat(v)) for j, v in PAR["c"].items()}
    KCR = arb("1.0865")

def poly_eval(cf, x):
    acc = x * 0
    for c in reversed(cf):
        acc = acc * x + c
    return acc

def hfun(n, x):
    """h_n(x) = He_n(x)/sqrt(n!) through Arb's Hermite function H_n (physicists'): He_n(x) = 2^(-n/2) H_n(x/sqrt 2)"""
    s2 = arb(2).sqrt()
    return (x / s2).hermite_h(n) * arb(2) ** (-arb(n) / 2) / arb(factorial(n)).sqrt()

def gpdf(z):
    return (-(z * z) / 2).exp() / S2PI

def qfun(a, p):
    if a == 0:
        return (p / arb(2).sqrt()).erf()
    return 2 * gpdf(p) * hfun(a - 1, -p) / arb(a).sqrt()

def upper_float(x):
    return float(arb(abs(x).upper()).upper()) * (1 + 1e-12)

# ------------------------------------------------------------------ head
XCB = Q(28, 5)
PANB = [(Q(k, 8), Q(k + 1, 8)) for k in range(36)] + [(Q(450 + 2 * k, 100), Q(452 + 2 * k, 100)) for k in range(55)]
NGB = 30

def head(out):
    init(320)
    t0 = time.time()
    # P increasing on [28/5, oo) with P(28/5) >= 25: P' > 0 on [28/5, 20] by interval cells, beyond 20 by the leading term
    ok = all(poly_eval(DPMON, arb((Q(28, 5) + Q(k, 10) + Q(1, 20)).numerator) / (Q(28, 5) + Q(k, 10) + Q(1, 20)).denominator
                       + arb(0, 0.0501)) > 0 for k in range(144))
    lead = DPMON[10]; rest = sum((abs(DPMON[j]) * arb(20) ** (j - 10) for j in range(10)), arb(0))
    pxc = poly_eval(PMON, ball(XCB))
    ok = ok and lead - rest > 0 and pxc > 25
    print("P increasing on [28/5, oo) and P(28/5) = %s > 25: %s" % (pxc.str(8), ok))
    if not ok:
        raise SystemExit("head: tail argument fails")
    acc = {}
    err = {}
    nodes = [arb.legendre_p_root(NGB, k, weight=True) for k in range(NGB)]
    n = NGB
    nfac = arb(factorial(n)); n2fac = arb(factorial(2 * n))
    for (lo, hi) in PANB:
        a0 = ball(lo); h = ball(hi - lo)
        for (xk, wk) in nodes:
            x = a0 + h * (xk + 1) / 2
            w = wk * h / 2 * gpdf(x)
            hb = [hfun(b, x) * w for b in range(NMAX + 1)]
            p = poly_eval(PMON, x)
            qa = [qfun(a, p) for a in range(NMAX + 1)]
            for a in range(NMAX + 1):
                for b in range(1 - a % 2, NMAX + 1 - a, 2):
                    acc[(a, b)] = acc.get((a, b), 0) + qa[a] * hb[b]
        # remainder: (b-a)^(2n+1) (n!)^4/((2n+1)((2n)!)^3) * (2n)! M / R^(2n), M = sup |f| on [lo-R, hi+R] x [-R, R]
        R = h * 0.75
        xl = a0 - R; xr = a0 + h + R
        x2max = (xr * xr); x2min = arb(0) if xl < 0 else xl * xl
        # |h_n(x+iy)| <= K e^(x^2/4) e^(|y| sqrt n); |phi(x+iy)| <= e^(-(x^2-y^2)/2)/sqrt(2 pi)
        base_h = KCR * (x2max / 4).exp() * (-(x2min - R * R) / 2).exp() / S2PI
        # P on the rectangle: centre value plus sum |P^(k)(c)/k!| rho^k (exact Taylor shift), rho the corner distance
        c = a0 + h / 2
        sh = arb_poly(PMON)(arb_poly([c, 1])).coeffs()
        rho = ((h / 2 + R) ** 2 + R * R).sqrt()
        rp = sum((abs(sh[k]) * rho ** k for k in range(1, len(sh))), arb(0))
        pc = sh[0]
        smin2 = arb(0) if (pc - rp < 0 and pc + rp > 0) or not (pc - rp > 0 or pc + rp < 0) else \
            (((pc - rp) ** 2) if pc - rp > 0 else ((pc + rp) ** 2))
        qbase = 2 * KCR / S2PI * (-smin2 / 4 + rp * rp / 2).exp()
        disk = acb(arb(pc, rp.upper()), arb(0, rp.upper()))
        Qb = [upper_float((disk / arb(2).sqrt()).erf())] + \
             [upper_float(qbase / arb(a).sqrt() * (rp * arb(a - 1).sqrt()).exp()) for a in range(1, NMAX + 1)]
        Hb = [upper_float(base_h * (R * arb(b).sqrt()).exp()) for b in range(NMAX + 1)]
        cf = upper_float(h ** (2 * n + 1) * nfac ** 4 / ((2 * n + 1) * n2fac ** 3) * n2fac / R ** (2 * n))
        for a in range(NMAX + 1):
            ca = cf * Qb[a]
            for b in range(1 - a % 2, NMAX + 1 - a, 2):
                err[(a, b)] = err.get((a, b), 0.0) + ca * Hb[b]
    print("Gauss-Legendre, %d panels of %d nodes (%.0f s)" % (len(PANB), NGB, time.time() - t0))
    xc = ball(XCB); phic = gpdf(xc)
    # tails on [28/5, oo): a >= 1: |h_b phi q_a(P)| <= 2 K^2/(2 pi sqrt a) e^(-x^2/4) e^(-P^2/4), P >= 25;
    # a = 0: int_c^oo h_b phi = h_(b-1)(c) phi(c)/sqrt b, and |erfc(P/sqrt 2)| <= e^(-P^2/2)
    ta = upper_float(2 * KCR * KCR / (2 * PIb) * PIb.sqrt() * (-(arb(25) ** 2) / 4).exp())
    t0b = upper_float(KCR / S2PI * PIb.sqrt() * (-(arb(25) ** 2) / 2).exp())
    A = {}
    worst = 0.0
    for (a, b), v in acc.items():
        if a == 0:
            v = v + hfun(b - 1, xc) * phic / arb(b).sqrt()
            e = err[(a, b)] + t0b
        else:
            e = err[(a, b)] + ta / a ** 0.5
        e = e * (1 + 1e-9)
        worst = max(worst, e)
        A[(a, b)] = 2 * (v + arb(0, e))
    print("pairs %d; largest remainder bound %.3e" % (len(A), worst))
    print("A_10 = %s, A_01 = %s" % (A[(1, 0)].str(25), A[(0, 1)].str(25)))
    rho_p = arb_poly([0] + [CB.get(j, 0) for j in range(1, 52)])
    tot = arb_poly([0])
    for b in range(NMAX + 1):
        T = arb_poly([0])
        for a in range(NMAX - b, -1, -1):           # Horner in rho
            T = arb_poly((T * rho_p).coeffs()[:NMAX + 1 - b])
            if (a + b) % 2 == 1:
                T = T + arb_poly([A[(a, b)] ** 2])
        sh = arb_poly([0] * b + [(-1) ** b])
        tot = tot + arb_poly((T * sh).coeffs()[:NMAX + 1])
    bn = [c * PIb / 2 for c in tot.coeffs()[:NMAX + 1]]
    b1 = bn[1]
    hs = sum((abs(bn[k]) for k in range(3, len(bn))), arb(0))
    print("b_1 in %s" % b1.str(25))
    print("sum_{3<=n<=301} |b_n| <= %s" % arb(hs.upper()).str(10))
    print("head time %.0f s" % (time.time() - t0))
    m, e = b1.lower().man_exp(); m2, e2 = hs.upper().man_exp()
    json.dump({"b1_lo": "%d*2^%d" % (int(m), int(e)), "head_hi": "%d*2^%d" % (int(m2), int(e2)),
               "b1": b1.str(30), "head": arb(hs.upper()).str(10), "bn_mid": [float(c.mid()) for c in bn]},
              open(out, "w"), indent=1)

# ------------------------------------------------------------------ the integrand (own derivation of eq. (Section 5.2)-(Section 5.2))
def setup_int(prec):
    init(prec)
    global RB, DELTAB, DAMPB, AWB
    RB = sum((abs(v) for v in CB.values()), arb(0))
    DELTAB = 1 - RB * RB
    DAMPB = 1 / (1 + RB)
    AWB = ball(Q(1, 1000))

def Fclosed(i, j, Lu, Lv, r, dinv):
    """phi_r^(-1) d_u^i d_v^j phi_r for phi_r ~ exp(l), l quadratic: l_uu = l_vv = -1/d, l_uv = r/d"""
    c = -dinv
    def G(m, L):
        return sum((L ** (m - 2 * l) * (c / 2) ** l * ball(Q(factorial(m), factorial(l) * factorial(m - 2 * l)))
                    for l in range(m // 2 + 1)), L * 0)
    luv = r * dinv
    return sum((comb(i, k) * comb(j, k) * factorial(k) * luv ** k * G(i - k, Lu) * G(j - k, Lv)
                for k in range(min(i, j) + 1)), Lu * 0)

def taylor_box(X, upto=4):
    """e_m enclosing P^(m)(x)/m! for every x in X, m <= upto: s_m at c = mid(X) plus the bound
    sum_{k>m} C(k,m) |s_k| r^(k-m) (Horner in r)"""
    c = arb(X.mid()); r = arb(X.rad())
    sh = arb_poly(PMON)(arb_poly([c, 1])).coeffs()
    sh = list(sh) + [arb(0)] * (12 - len(sh))
    out = []
    for m in range(upto + 1):
        acc = arb(0)
        for k in range(11, m, -1):
            acc = acc * r + abs(sh[k]) * comb(k, m)
        acc = acc * r
        out.append(sh[m] + arb(0, acc.upper()))
    return out

def bell_series(X, kmax=4):
    """B_(k,i)(x) = (k!/i!) [w^k] (P(x+w) - P(x))^i for all x in X, from the enclosures of the Taylor coefficients"""
    e = taylor_box(X, kmax)
    inc = arb_poly([arb(0)] + e[1:kmax + 1])
    out = [[arb(0)] * (k + 1) for k in range(kmax + 1)]
    out[0][0] = arb(1)
    pw = arb_poly([1])
    for i in range(1, kmax + 1):
        pw = arb_poly((pw * inc).coeffs()[:kmax + 1])
        cf = pw.coeffs()
        for k in range(i, kmax + 1):
            if k < len(cf):
                out[k][i] = cf[k] * ball(Q(factorial(k), factorial(i)))
    return out

def Pint(X):
    """P on an interval: exact endpoints where P is monotone (|x| >= 28/5), else the Taylor-shift (mean value) form"""
    lo = arb(X.lower()); hi = arb(X.upper())
    if lo >= arb(28) / 5 or hi <= -arb(28) / 5:
        return poly_eval(PMON, lo).union(poly_eval(PMON, hi))
    return taylor_box(X, 0)[0]

KSG = 1

def kap_b(t, q):
    q1, q2, q3 = q
    return {(1, 0): q3, (2, 0): 3 * q1 * q2, (3, 0): q1 * q1 * q1, (0, 1): -t, (0, 2): 3 * t * t,
            (0, 3): -t * t * t, (1, 1): -3 * KSG * t * (q2 + q1), (2, 1): -3 * t * q1 * q1, (1, 2): 3 * t * t * q1}

def Rpref_b(e, BX, BY, u, v, r, dinv, kap):
    Lu = (-u + r * v) * dinv; Lv = (-v + r * u) * dinv
    cache = {}
    def F(i, j):
        if (i, j) not in cache:
            cache[(i, j)] = Fclosed(i, j, Lu, Lv, r, dinv)
        return cache[(i, j)]
    tot = 0
    for (p, a), kv in kap.items():
        m = a + e
        s = 0
        rng = range(0, m + 1) if m == 0 else range(1, m + 1)
        for i in rng:
            for j in rng:
                if p == 0:
                    g = F(i - 1, j - 1)
                else:
                    g = F(i + p - 1, j + p - 1)
                s = s + BX[m][i] * BY[m][j] * g
        tot = tot + kv * s
    return tot

def rq(t):
    r = q1 = q2 = q3 = t * 0
    for j, cj in CB.items():
        tj = t ** j * cj
        r += tj; q1 += tj * j; q2 += tj * (j * j); q3 += tj * (j ** 3)
    return r, (q1, q2, q3)

def sqlo(I):
    lo = arb(I.lower().lower()); hi = arb(I.upper().upper())
    if lo > 0: return arb((lo * lo).lower())
    if hi < 0: return arb((hi * hi).lower())
    return arb(0)

def upper_box(X, Y, S):
    t = acb(0, PIb * S).exp()
    r, q = rq(t)
    u = Pint(X); v = Pint(Y)
    d = 1 - r * r; dinv = 1 / d
    kap = kap_b(t, q)
    BX = bell_series(X); BY = bell_series(Y)
    R0 = Rpref_b(0, BX, BY, u, v, r, dinv, kap)
    R1 = Rpref_b(1, BX, BY, u, v, r, dinv, kap)
    num = abs(R0) ** 2 + AWB * abs(R1) ** 2
    dabs = abs(d).lower()
    den = 2 * PIb * (dabs if dabs > DELTAB else DELTAB)
    cm = (1 / (1 + r)).real.lower(); cp = (1 / (1 - r)).real.lower()
    cm = cm if cm > DAMPB else DAMPB; cp = cp if cp > DAMPB else DAMPB
    e1 = cm * sqlo(u + v) / 2 + cp * sqlo(u - v) / 2
    e2 = (cm if cm < cp else cp) * (sqlo(u) + sqlo(v))
    e = arb(e1.lower()) if e1.lower() > e2.lower() else arb(e2.lower())
    return arb((num / den * (-(e + (sqlo(X) + sqlo(Y)) / 2)).exp()).upper())

def phi3_real(tr, x, y, e):
    t = acb(tr); r, q = rq(t)
    X = arb(x); Y = arb(y)
    u = poly_eval(PMON, X); v = poly_eval(PMON, Y)
    d = 1 - r * r; dinv = 1 / d
    Re = Rpref_b(e, bell_series(X), bell_series(Y), u, v, r, dinv, kap_b(t, q))
    return (Re * (-(u * u - 2 * r * u * v + v * v) * dinv / 2).exp() / d.sqrt()).real

def selftest(headjson):
    from math import sqrt, exp, pi
    setup_int(96)
    bn = json.load(open(headjson))["bn_mid"]
    ok = True
    for tr in (0.4, -0.6):
        lhs = sum(n ** 3 * bn[n] * tr ** n for n in range(len(bn)))
        sc = -tr; w = sqrt(1 - sc * sc); h = 0.08; L = 8; nn = int(L / h); tot = 0.0
        for i in range(-nn, nn + 1):
            x = i * h
            for k in range(-nn, nn + 1):
                z = k * h
                tot += float(phi3_real(tr, x, sc * x + w * z, 0).mid()) * exp(-(x * x + z * z) / 2)
        rhs = tot * h * h / (2 * pi)
        good = abs(lhs - rhs) < 1e-6 * max(1.0, abs(lhs)); ok &= good
        print("self-test D^3 H(%s): sum n^3 b_n t^n = %.10f, Gaussian expectation of Phi_3 = %.10f: %s" % (tr, lhs, rhs, good))
    # the e = 1 part: E_s[d_xy Phi_3] = d/ds E_s[Phi_3] (Mehler), at t = 0.4, s = -0.4
    def ex(tr, sc, e, h=0.08, L=8):
        w = sqrt(1 - sc * sc); nn = int(L / h); tot = 0.0
        for i in range(-nn, nn + 1):
            x = i * h
            for k in range(-nn, nn + 1):
                z = k * h
                tot += float(phi3_real(tr, x, sc * x + w * z, e).mid()) * exp(-(x * x + z * z) / 2)
        return tot * h * h / (2 * pi)
    d1 = ex(0.4, -0.4, 1); dd = (ex(0.4, -0.399, 0) - ex(0.4, -0.401, 0)) / 0.002
    good = abs(d1 - dd) < 1e-5 * max(1.0, abs(d1)); ok &= good
    print("self-test E_s[d_xy Phi_3] = %.8f against d/ds E_s[Phi_3] = %.8f: %s" % (d1, dd, good))
    return ok

# ------------------------------------------------------------------ b3: coverage and bounds on program A's partition
XU, YU, SU = 16, 16, 24

def coverage(boxes):
    """the boxes tile [0,8] x [-8,8] x [0,1] exactly: every initial cell (8 x 16 x 64) is reconstructed by bisections"""
    cells = {}
    for b in boxes:
        x0, x1, y0, y1, s0, s1 = b
        if not (0 <= x0 < x1 <= 8 << XU and -(8 << YU) <= y0 < y1 <= 8 << YU and 0 <= s0 < s1 <= 1 << SU):
            return False, "a box leaves the domain"
        key = (x0 >> XU, (y0 + (8 << YU)) >> YU, s0 >> (SU - 6))
        cells.setdefault(key, []).append(b)
    if len(cells) != 8 * 16 * 64:
        return False, "an initial cell has no box"
    sys.setrecursionlimit(10000)
    def rec(cell, bs):
        if len(bs) == 1:
            return bs[0] == cell
        for dim in range(3):
            lo, hi = cell[2 * dim], cell[2 * dim + 1]
            if hi - lo < 2:
                continue
            m = (lo + hi) // 2
            left = [b for b in bs if b[2 * dim + 1] <= m]
            right = [b for b in bs if b[2 * dim] >= m]
            if len(left) + len(right) == len(bs) and left and right:
                c1 = list(cell); c1[2 * dim + 1] = m
                c2 = list(cell); c2[2 * dim] = m
                return rec(tuple(c1), left) and rec(tuple(c2), right)
        return False
    for (i, j, k), bs in cells.items():
        cell = (i << XU, (i + 1) << XU, (j - 8) << YU, (j - 7) << YU, k << (SU - 6), (k + 1) << (SU - 6))
        if not rec(cell, bs):
            return False, "cell %s is not tiled" % (cell,)
    vol = sum((b[1] - b[0]) * (b[3] - b[2]) * (b[5] - b[4]) for b in boxes)
    if vol != (8 << XU) * (16 << YU) * (1 << SU):
        return False, "volumes do not add up"
    return True, "%d boxes tile the domain" % len(boxes)

def exterior_b():
    """own exterior bound, |x| > 8 or |y| > 8: integrand <= (|R_0|^2 + a|R_1|^2)/(2 pi delta) e^(-(u^2+v^2)/2 - (x^2+y^2)/2);
    |R_e| <= C_e max(1,|x|)^106 max(1,|y|)^106 with C_e from absolute values; |P(x)| >= c8 |x|^11 for |x| >= 8"""
    c8 = PMON[11] - sum((abs(PMON[j]) * arb(8) ** (j - 11) for j in range(11)), arb(0))
    Ph = sum((abs(c) for c in PMON), arb(0))
    Dd = 1 / DELTAB
    L = 2 * Ph * Dd
    def Fh(i, j):
        G = lambda m: sum((L ** (m - 2 * l) * (Dd / 2) ** l * ball(Q(factorial(m), factorial(l) * factorial(m - 2 * l)))
                           for l in range(m // 2 + 1)), arb(0))
        return sum((comb(i, k) * comb(j, k) * factorial(k) * Dd ** k * G(i - k) * G(j - k) for k in range(min(i, j) + 1)), arb(0))
    qh = [sum((arb(j) ** k * abs(c) for j, c in CB.items()), arb(0)) for k in range(4)]
    kh = {(1, 0): qh[3], (2, 0): 3 * qh[1] * qh[2], (3, 0): qh[1] ** 3, (0, 1): arb(1), (0, 2): arb(3), (0, 3): arb(1),
          (1, 1): 3 * (qh[2] + qh[1]), (2, 1): 3 * qh[1] ** 2, (1, 2): 3 * qh[1]}
    # Bell coefficients bounded by the Bell polynomials of |P|-coefficient series at 1
    absP = [abs(c) for c in PMON]
    sh = [sum((abs(absP[n]) * comb(n, k) for n in range(k, 12)), arb(0)) for k in range(12)]   # |coeffs| of P(1+w)
    inc = arb_poly([0] + sh[1:5])
    Bh = [[arb(0)] * (k + 1) for k in range(5)]; Bh[0][0] = arb(1)
    pw = arb_poly([1])
    for i in range(1, 5):
        pw = arb_poly((pw * inc).coeffs()[:5]); cf = pw.coeffs()
        for k in range(i, 5):
            if k < len(cf):
                Bh[k][i] = cf[k] * ball(Q(factorial(k), factorial(i)))
    def Ch(e):
        tot = arb(0)
        for (p, a), kv in kh.items():
            m = a + e
            rng = range(0, m + 1) if m == 0 else range(1, m + 1)
            s = arb(0)
            for i in rng:
                for j in rng:
                    s += Bh[m][i] * Bh[m][j] * (Fh(i - 1, j - 1) if p == 0 else Fh(i + p - 1, j + p - 1))
            tot += kv * s
        return tot
    C0, C1 = Ch(0), Ch(1)
    K = 212
    mom = arb(2) ** ((arb(K) - 1) / 2) * ((arb(K) + 1) / 2).gamma()
    u8 = c8 * arb(8) ** 11
    ext = 4 * (C0 * C0 + AWB * C1 * C1) / (2 * PIb * DELTAB) * (-(u8 * u8) / 4).exp() * mom * (2 + 2 * mom)
    return c8, u8, C0, C1, ext

def b3check(boxes_path, out, part=0, parts=1):
    setup_int(96)
    t0 = time.time()
    boxes = [tuple(int(v) for v in line.split()) for line in gzip.open(boxes_path, "rt")]
    if part == 0:
        ok, msg = coverage(boxes)
        print("coverage: %s (%s) (%.0f s)" % (ok, msg, time.time() - t0))
        if not ok:
            raise SystemExit("b3: the partition does not tile the domain")
    boxes = boxes[part::parts]
    c8, u8, C0, C1, ext = exterior_b()
    print("exterior (own constants): c_8 = %s, c_8 8^11 = %s, C_0 <= %s, C_1 <= %s, bound %s"
          % (c8.str(10), u8.str(10), C0.str(4), C1.str(4), arb(ext.upper()).str(3)))
    if part != 0:
        ext = arb(0)
        print("part %d of %d: the exterior is counted in part 0" % (part, parts))
    tot = arb(0)
    for n, b in enumerate(boxes):
        x0, x1, y0, y1, s0, s1 = b
        X = (arb(x0) / 2 ** XU).union(arb(x1) / 2 ** XU)
        Y = (arb(y0) / 2 ** YU).union(arb(y1) / 2 ** YU)
        S = (arb(s0) / 2 ** SU).union(arb(s1) / 2 ** SU)
        U = upper_box(X, Y, S)
        vol = arb((x1 - x0) * (y1 - y0) * (s1 - s0)) / 2 ** (XU + YU + SU)
        tot += 2 * U * vol
        if (n + 1) % 50000 == 0:
            print("  %d boxes, running sum %s (%.0f s)" % (n + 1, arb(tot.upper()).str(6), time.time() - t0))
    I = arb((tot + ext).upper())
    print("program B on program A's partition: I <= %s (%.0f s)" % (I.str(10), time.time() - t0))
    m, e = I.man_exp()
    json.dump({"I_hi": "%d*2^%d" % (int(m), int(e)), "I": I.str(12), "boxes": len(boxes)}, open(out, "w"), indent=1)

def final(headjson, b3json, claim):
    setup_int(160)
    def rd(st):
        m, e = st.split("*2^"); return ball(Q(int(m)) * Q(2) ** int(e))
    hj = json.load(open(headjson)); bj = json.load(open(b3json))
    b1 = rd(hj["b1_lo"]); hs = rd(hj["head_hi"]); I = rd(bj["I_hi"])
    gam = PIb / (2 * ball(Q(claim)))
    aw = ball(Q(1, 1000)); sa = aw.sqrt()
    Ca = arb(1) / 2 + PIb / (2 * sa) * ((PIb / sa).exp() + (-(PIb / sa)).exp()) / ((PIb / sa).exp() - (-(PIb / sa)).exp())
    B3 = (Ca * I).sqrt()
    tail = B3 / (arb(10) * arb(301) ** 5).sqrt()
    margin = b1 - gam - hs - tail
    print("gamma = %s; b_1 > %s; head < %s" % (gam.str(15), b1.str(15), hs.str(10)))
    print("C_a = %s; I <= %s; B_3 <= %s; tail <= %s" % (Ca.str(10), I.str(10), arb(B3.upper()).str(8), arb(tail.upper()).str(8)))
    print("lemma 5.1: b_1 > 0.881850816: %s; head < 3e-8: %s; B_3 < 107: %s" % (b1 > arb("0.881850816"), hs < arb("3e-8"), B3 < 107))
    if margin > 0:
        print("margin b_1 - gamma - head - tail = %s > 0" % margin.str(8))
        print("RESULT: K_G <= %s (with the analytic lemmas of the paper, not re-proved)" % claim)
        print("VERDICT PASS")
    else:
        print("margin %s: VERDICT REFUSED: the criterion is not established" % margin.str(8))

if __name__ == "__main__":
    part = sys.argv[1]
    if "--c1-shift" in sys.argv:
        k = sys.argv.index("--c1-shift"); sh = Q(sys.argv[k + 1]); del sys.argv[k:k + 2]
        PAR["c"][1] = str(Q(PAR["c"][1]) + sh)
        print("CONTROL: c_1 shifted by %s" % sh)
    if "--kappa-sign" in sys.argv:
        sys.argv.remove("--kappa-sign"); KSG = -1
        print("CONTROL: kappa_(1,1) with its sign flipped")
    if part == "head":
        head(sys.argv[2])
    elif part == "selftest":
        print("VERDICT SELFTEST PASS" if selftest(sys.argv[2]) else "VERDICT REFUSED (self-test fails)")
    elif part == "b3":
        if len(sys.argv) > 4:
            i, n = sys.argv[4].split("/")
            b3check(sys.argv[2], sys.argv[3], int(i), int(n))
        else:
            b3check(sys.argv[2], sys.argv[3])
    elif part == "b3sum":
        ctx.prec = 128
        tot = arb(0); nb = 0
        for f in sys.argv[3:]:
            j = json.load(open(f)); m, e = j["I_hi"].split("*2^")
            tot += ball(Q(int(m)) * Q(2) ** int(e)); nb += j["boxes"]
        m, e = tot.upper().man_exp()
        json.dump({"I_hi": "%d*2^%d" % (int(m), int(e)), "I": arb(tot.upper()).str(12), "boxes": nb}, open(sys.argv[2], "w"), indent=1)
        print("program B, all parts: %d boxes, I <= %s" % (nb, arb(tot.upper()).str(10)))
    elif part == "final":
        final(sys.argv[2], sys.argv[3], sys.argv[4] if len(sys.argv) > 4 else "17813/10000")
