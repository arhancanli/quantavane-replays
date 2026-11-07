"""Independent check of the C_42 upper bounds 0.6906538 [Gri26], 0.6888863712 [Pin26] and 0.688983 [R2026a] (Turan's
pure power-sum constant, constants/42a.md; the last is the C42 row of PR #184), from the criterion as Griego's note
(section 4), Pinho's paper (Proposition 2) and Roehrig's notes/turan.md state it, and the exact inputs of their
certificates alone: verify_42a_certificate.py, verify_primary.py, independent/check_bounds.py, certify_turan.py and
audit.py were neither read nor run.

  ./fetch.sh && python3 -I c42_check.py [CERTS_ROOT] [--controls] [--finite] [--rows=gri,pin,roe]

The criterion, for a first-block cut-off tau >= 1/3 (then at most two corrections reach b_n): with
tau = t_0 < t_1 < ... < t_M = 1 - tau, alpha = a + ib (a > 0), s = 1 - alpha, middle values eta_0, ..., eta_{M-1},
delta_0 = eta_0 - s and delta_j = eta_j - eta_{j-1},
  I_A(x) = int_0^x v^(A-1) / (1 - v) dv,   K = I_alpha(tau),   D = I_a(tau),
  J(t, r) = int_0^(1-t-r) v^(alpha-1) / (1 - v) log((1-t-v)(1-r-v) / (t r)) dv   (0 if t + r >= 1),
  Y = 1 + s K - sum_j delta_j (I_alpha(1 - t_j) - K) + 1/2 sum_(i,j) delta_i delta_j J(t_i, t_j),
and if |s| < C, every |eta_j| < C and |Y| < C D, then C_42 <= C. Griego's note is the case M = 1 (w = delta_0,
A_1 = I_alpha(1 - tau) - K, A_2 = J(tau, tau)). We check the radius conditions in exact rationals, and evaluate I and J
by their power series, I_A(x) = sum_r x^(A+r) / (A+r) and J = sum_m c_m L^(alpha+m) / (alpha+m) with L = 1 - t - r and
c_m = log((1-t)(1-r)/(t r)) - sum_(1<=q<=m) ((1-t)^-q + (1-r)^-q) / q, with explicit tail bounds, in decimal interval
arithmetic rounded outward (x^alpha = e^(a ln x) (cos(b ln x) + i sin(b ln x)), sine and cosine by Taylor series with
the remainder bounded); then |Y|^2 < C^2 D^2 from the corners of Y's rectangle. --finite adds a floating-point look at
the construction itself, |Y_n| / sum_F |P_m| for finite n, which the criterion says tends to |Y|/D: a sanity check of
the limit, not part of the proof. Python 3.9+, standard library only."""
import cmath
import decimal
import json
import math
import os
import sys
import time
from decimal import Decimal
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
PREC = 50
EPS_DIGITS = PREC - 8                      # every series is cut where its tail is below 10^-EPS_DIGITS
LO = decimal.Context(prec=PREC, rounding=decimal.ROUND_FLOOR)
HI = decimal.Context(prec=PREC, rounding=decimal.ROUND_CEILING)
NEAR = decimal.Context(prec=PREC, rounding=decimal.ROUND_HALF_EVEN)
HARCOS = Fraction(69368, 100000)           # the page's best unstarred upper bound [Bir00]


# ---- interval arithmetic on decimals (as in our c3a and c3e replays) ------------------------------------------------
class I:
    """A closed interval [lo, hi] of decimals; every operation rounds outward."""

    def __init__(self, lo, hi=None):
        self.lo, self.hi = Decimal(lo), Decimal(lo if hi is None else hi)
        assert self.lo <= self.hi

    def __add__(self, o):
        o = o if isinstance(o, I) else I(o)
        return I(LO.add(self.lo, o.lo), HI.add(self.hi, o.hi))

    def __sub__(self, o):
        o = o if isinstance(o, I) else I(o)
        return I(LO.subtract(self.lo, o.hi), HI.subtract(self.hi, o.lo))

    def __neg__(self):
        # copy_negate is exact; unary minus would round in decimal's default context (28 digits, half-even)
        return I(self.hi.copy_negate(), self.lo.copy_negate())

    def __mul__(self, o):
        o = o if isinstance(o, I) else I(o)
        ps = [(a, b) for a in (self.lo, self.hi) for b in (o.lo, o.hi)]
        return I(min(LO.multiply(a, b) for a, b in ps), max(HI.multiply(a, b) for a, b in ps))

    def __truediv__(self, o):
        o = o if isinstance(o, I) else I(o)
        assert o.lo > 0 or o.hi < 0, "division by an interval containing 0"
        ps = [(a, b) for a in (self.lo, self.hi) for b in (o.lo, o.hi)]
        return I(min(LO.divide(a, b) for a, b in ps), max(HI.divide(a, b) for a, b in ps))

    def mag(self):
        """An upper bound for |x| over the interval (copy_abs is exact; abs() would round to 28 digits)."""
        return max(self.lo.copy_abs(), self.hi.copy_abs())

    def widen(self, r):
        """[lo - r, hi + r] for a decimal r >= 0."""
        return I(LO.subtract(self.lo, r), HI.add(self.hi, r))


def q(x):
    """The rational x as an interval: the exact quotient rounded down and up."""
    x = Fraction(x)
    return I(LO.divide(Decimal(x.numerator), Decimal(x.denominator)),
             HI.divide(Decimal(x.numerator), Decimal(x.denominator)))


def iln(x):
    """ln of a positive interval: decimal's ln is correctly rounded, so one unit either way encloses the true value."""
    assert x.lo > 0
    return I(NEAR.next_minus(NEAR.ln(x.lo)), NEAR.next_plus(NEAR.ln(x.hi)))


def iexp(x):
    """exp of an interval (decimal's exp is correctly rounded)."""
    return I(NEAR.next_minus(NEAR.exp(x.lo)), NEAR.next_plus(NEAR.exp(x.hi)))


def isincos(x):
    """(sin x, cos x) for an interval x with |x| < 8: Taylor series at the point x.lo with the remainder bounded, then
    widened by the width of x (both functions are 1-Lipschitz)."""
    p = I(x.lo)
    assert p.mag() < 8
    tol = Decimal(10) ** -(EPS_DIGITS + 4)
    s, c = I(0), I(1)
    term = I(1)                                            # p^k / k!
    k = 0
    while True:
        k += 1
        term = term * p / I(k)
        if k % 2:
            s = s + term if k % 4 == 1 else s - term
        else:
            c = c + term if k % 4 == 0 else c - term
        if k > 20 and term.mag() < tol:
            break
    # the remaining terms are at most |p|^(k+1)/(k+1)! (1 + 1/2 + 1/4 + ...) once k + 2 > 2|p|
    rem = HI.multiply(HI.multiply(term.mag(), HI.divide(p.mag(), Decimal(k + 1))), Decimal(2))
    w = HI.subtract(x.hi, x.lo)
    return s.widen(rem).widen(w), c.widen(rem).widen(w)


class Z:
    """A complex number as a rectangle re + i im of two intervals."""

    def __init__(self, re, im=None):
        self.re, self.im = re, (im if im is not None else I(0))

    def __add__(self, o):
        return Z(self.re + o.re, self.im + o.im)

    def __sub__(self, o):
        return Z(self.re - o.re, self.im - o.im)

    def __mul__(self, o):
        if isinstance(o, I):
            return Z(self.re * o, self.im * o)
        return Z(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)

    def widen(self, r):
        return Z(self.re.widen(r), self.im.widen(r))

    def abs2_hi(self):
        """An upper bound for |z|^2 over the rectangle (x^2 + y^2 is largest at a corner)."""
        x, y = self.re.mag(), self.im.mag()
        return HI.add(HI.multiply(x, x), HI.multiply(y, y))


def zq(re, im):
    """The complex rational re + i im as a rectangle."""
    return Z(q(re), q(im))


def zinv(re, im):
    """1/(re + i im) for rationals, exactly, as a rectangle."""
    n = Fraction(re) ** 2 + Fraction(im) ** 2
    return zq(Fraction(re) / n, -Fraction(im) / n)


def cpow(x, ar, ai):
    """x^(ar + i ai) for a rational 0 < x < 1 (real logarithm): e^(ar ln x) (cos(ai ln x) + i sin(ai ln x))."""
    lx = iln(q(x))
    mod = iexp(q(ar) * lx)
    s, c = isincos(q(ai) * lx)
    return Z(mod * c, mod * s)


def terms_for(ratio):
    """The number of terms after which ratio^R < 10^-EPS_DIGITS (0 < ratio < 1)."""
    return int(math.ceil(EPS_DIGITS * math.log(10) / -math.log(float(ratio)))) + 2


# ---- the integrals --------------------------------------------------------------------------------------------------
def I_series(x, ar, ai):
    """I_A(x) = sum_r x^(A+r) / (A+r) for A = ar + i ai (ai = 0: real), with the tail bound
    |sum_{r>=R}| <= x^(ar+R) / ((ar+R)(1-x))."""
    x = Fraction(x)
    R = terms_for(x)
    xa = cpow(x, ar, ai) if ai else Z(iexp(q(ar) * iln(q(x))))
    xr = I(1)
    total = Z(I(0))
    for r in range(R):
        total = total + xa * xr * zinv(Fraction(ar) + r, ai)
        xr = xr * q(x)
    tail = (iexp(q(ar) * iln(q(x))) * xr / (q(Fraction(ar) + R) * q(1 - x))).hi
    return total.widen(tail)


class Pairs:
    """J(t, r) for many pairs, sharing the parts that depend on one break or on L = 1 - t - r only."""

    def __init__(self, ar, ai):
        self.ar, self.ai = Fraction(ar), Fraction(ai)
        self.W = {}       # L -> [L^(alpha+m) / (alpha+m)] as rectangles, and L^a as an interval
        self.V = {}       # t -> [log((1-t)/t) - sum_{q<=m} (1-t)^-q / q]
        self.inv = []     # 1/(alpha+m)

    def _inv(self, m):
        while len(self.inv) <= m:
            self.inv.append(zinv(self.ar + len(self.inv), self.ai))
        return self.inv[m]

    def _W(self, L, R):
        w = self.W.get(L)
        if w is None or len(w[0]) < R:
            La = cpow(L, self.ar, self.ai)
            Lm, out = I(1), []
            for m in range(R):
                out.append(La * Lm * self._inv(m))
                Lm = Lm * q(L)
            w = self.W[L] = (out, iexp(q(self.ar) * iln(q(L))))
        return w

    def _V(self, t, R):
        v = self.V.get(t)
        if v is None or len(v) < R:
            rho = q(1 / (1 - t))
            cur, pw, out = iln(q((1 - t) / t)), I(1), []
            for m in range(R):
                if m:
                    pw = pw * rho
                    cur = cur - pw / I(m)
                out.append(cur)
            v = self.V[t] = out
        return v

    def J(self, t, r):
        t, r = Fraction(t), Fraction(r)
        L = 1 - t - r
        if L <= 0:
            return Z(I(0))
        rho = max(1 / (1 - t), 1 / (1 - r))
        R = terms_for(L * rho)
        W, La = self._W(L, R)
        Vt, Vr = self._V(t, R), self._V(r, R)
        total = Z(I(0))
        for m in range(R):
            total = total + W[m] * (Vt[m] + Vr[m])
        # tail: |c_m| <= |lambda| + sum_rho rho^(m+1)/(rho-1), so sum_{m>=R} |c_m| L^(a+m)/(a+m) is at most
        # L^a/(a+R) (|lambda| L^R/(1-L) + sum_rho rho/(rho-1) (L rho)^R/(1 - L rho))
        lam = iln(q((1 - t) * (1 - r) / (t * r))).mag()
        tail = HI.divide(HI.multiply(lam, ipow(q(L), R).hi), q(1 - L).lo)
        for p in (1 / (1 - t), 1 / (1 - r)):
            g = HI.divide(HI.multiply(q(p / (p - 1)).hi, ipow(q(L * p), R).hi), q(1 - L * p).lo)
            tail = HI.add(tail, g)
        tail = HI.multiply(tail, HI.divide(La.hi, q(self.ar + R).lo))
        return total.widen(tail)


def ipow(x, n):
    out = I(1)
    for _ in range(n):
        out = out * x
    return out


# ---- the criterion ----------------------------------------------------------------------------------------------------
def load(path):
    """(name, C, alpha, breaks, etas) as exact rationals from either certificate's JSON."""
    d = json.load(open(path))
    if "tau" in d:                                         # Griego: one middle block
        tau = Fraction(d["tau"])
        return ("Gri26", Fraction(d["bound"]), (Fraction(d["alpha"]["re"]), Fraction(d["alpha"]["im"])),
                [tau, 1 - tau], [(Fraction(d["eta"]["re"]), Fraction(d["eta"]["im"]))], d)
    name = "R2026a" if "construction" in d else "Pin26"   # Roehrig's turan.json has Pinho's layout plus a note
    return (name, Fraction(d["C"]), (Fraction(d["alpha"][0]), Fraction(d["alpha"][1])),
            [Fraction(b) for b in d["breaks"]], [(Fraction(e[0]), Fraction(e[1])) for e in d["eta"]], d)


def evaluate(C, alpha, breaks, etas, out=print, failed=None):
    """Check the criterion; returns (ok, details). The names of failed checks are appended to `failed`."""
    ar, ai = alpha
    tau, M = breaks[0], len(etas)
    ok = True

    def check(what, cond, detail=""):
        nonlocal ok
        ok &= bool(cond)
        if not cond and failed is not None:
            failed.append(what)
        out(f"  check {what}: {'PASS' if cond else 'FAIL'}{' (' + detail + ')' if detail else ''}")

    s = (1 - ar, -ai)
    check("the cut-off tau >= 1/3, so at most two corrections reach b_n", tau >= Fraction(1, 3), f"tau = {tau}")
    check("breaks increase from tau to 1 - tau, one more than the middle values",
          len(breaks) == M + 1 and breaks[-1] == 1 - tau and all(x < y for x, y in zip(breaks, breaks[1:])))
    check("a = Re alpha > 0", ar > 0)
    m_s = C * C - (s[0] ** 2 + s[1] ** 2)
    m_eta = min(C * C - (e[0] ** 2 + e[1] ** 2) for e in etas)
    check("|s| < C and every |eta_j| < C (exact)", m_s > 0 and m_eta > 0,
          f"least C^2 - |.|^2 = {float(min(m_s, m_eta)):.6e}")
    check(f"C < {float(HARCOS)} (Harcos, the best unstarred upper bound)", C < HARCOS)
    if not ok:
        return False, None
    deltas = [(etas[0][0] - s[0], etas[0][1] - s[1])] + \
             [(etas[j][0] - etas[j - 1][0], etas[j][1] - etas[j - 1][1]) for j in range(1, M)]
    K = I_series(tau, ar, ai)
    D = I_series(tau, ar, 0).re
    Y = Z(I(1)) + zq(*s) * K
    single = Z(I(0))
    for j in range(M):
        single = single + zq(*deltas[j]) * (I_series(1 - breaks[j], ar, ai) - K)
    Y = Y - single
    P = Pairs(ar, ai)
    double = Z(I(0))
    npairs = 0
    for i in range(M):
        for j in range(i, M):
            if breaks[i] + breaks[j] >= 1:
                break
            npairs += 1
            dij = (deltas[i][0] * deltas[j][0] - deltas[i][1] * deltas[j][1],
                   deltas[i][0] * deltas[j][1] + deltas[i][1] * deltas[j][0])
            w = Fraction(1, 2) if i == j else 1        # 1/2 sum over ordered pairs = sum over i < j + 1/2 diagonal
            double = double + zq(dij[0] * w, dij[1] * w) * P.J(breaks[i], breaks[j])
    Y = Y + double
    y2 = Y.abs2_hi()
    c2d2 = (q(C * C) * D * D).lo
    # displayed ends rounded outward to 30 decimals (in a context wide enough for the quantize)
    dn = lambda x: decimal.Context(prec=PREC + 10, rounding=decimal.ROUND_FLOOR).quantize(x, Decimal(10) ** -30)
    up = lambda x: decimal.Context(prec=PREC + 10, rounding=decimal.ROUND_CEILING).quantize(x, Decimal(10) ** -30)
    out(f"  K = I_alpha(tau) in Re [{dn(K.re.lo)}, {up(K.re.hi)}], Im [{dn(K.im.lo)}, {up(K.im.hi)}]")
    out(f"  D = I_a(tau) in [{dn(D.lo)}, {up(D.hi)}]")
    out(f"  Y in Re [{dn(Y.re.lo)}, {up(Y.re.hi)}]")
    out(f"       Im [{dn(Y.im.lo)}, {up(Y.im.hi)}]   ({M} middle blocks, {npairs} pairs i <= j with t_i + t_j < 1)")
    ratio = NEAR.divide(NEAR.sqrt(y2), D.lo)
    out(f"  |Y|/D <= {ratio:.25f}   (C = {float(C):.10f}; C^2 D^2 - |Y|^2 >= {HI.subtract(c2d2, y2):.6e})")
    check("|Y|^2 < C^2 D^2, so C_42 <= C", y2 < c2d2)
    return ok, {"K": K, "D": D, "Y": Y, "ratio": ratio}


def compare_griego(d, res, out=print):
    """Our enclosures against the ones Griego's certificate prints (each pair must overlap)."""
    enc = d["enclosures"]

    def overlap(ours, theirs):
        lo, hi = Fraction(theirs[0]), Fraction(theirs[1])
        return Fraction(ours.lo) <= hi and lo <= Fraction(ours.hi)
    good = overlap(res["K"].re, enc["K"]["re"]) and overlap(res["K"].im, enc["K"]["im"]) and \
        overlap(res["Y"].re, enc["Y"]["re"]) and overlap(res["Y"].im, enc["Y"]["im"]) and \
        Fraction(res["D"].hi) > Fraction(enc["D_lower"])
    out(f"  our K, D and Y agree with the certificate's printed enclosures: {'yes' if good else 'NO'}")
    return good


# the enclosures Roehrig's notes/turan.md prints (30 decimals; the 42a row of PR #184 prints 17 of them)
ROEHRIG_PRINTED = {"D": ("0.888190650452511954178223193231", "0.888190650452511954178223193232"),
                   "ratio": ("0.688982098504060748505120997583", "0.688982098504060748505120997584"),
                   "gap": ("0.000000979973624636182168032566", "0.000000979973624636182168032567")}


def compare_roehrig(C, res, out=print):
    """Our D, |Y|/D and C^2 D^2 - |Y|^2 against the enclosures Roehrig's note prints (each pair must overlap)."""
    P = {k: (Fraction(a), Fraction(b)) for k, (a, b) in ROEHRIG_PRINTED.items()}
    D, Y = res["D"], res["Y"]
    y2_hi = Fraction(Y.abs2_hi())
    xr, xi = Y.re, Y.im
    # the least |Y|^2 over the rectangle (copy_abs is exact; abs() would round to 28 digits)
    y2_lo = min(Fraction(v.copy_abs()) for v in (xr.lo, xr.hi)) ** 2 * (0 if xr.lo <= 0 <= xr.hi else 1) + \
        min(Fraction(v.copy_abs()) for v in (xi.lo, xi.hi)) ** 2 * (0 if xi.lo <= 0 <= xi.hi else 1)
    d_lo, d_hi = Fraction(D.lo), Fraction(D.hi)
    gap_lo, gap_hi = C * C * d_lo * d_lo - y2_hi, C * C * d_hi * d_hi - y2_lo
    # |Y|/D lies in [sqrt(y2_lo)/d_hi, sqrt(y2_hi)/d_lo]; compare squares to stay exact
    r_lo2, r_hi2 = y2_lo / (d_hi * d_hi), y2_hi / (d_lo * d_lo)
    good = d_lo <= P["D"][1] and P["D"][0] <= d_hi and r_lo2 <= P["ratio"][1] ** 2 and P["ratio"][0] ** 2 <= r_hi2 \
        and gap_lo <= P["gap"][1] and P["gap"][0] <= gap_hi
    out(f"  C^2 D^2 - |Y|^2 in [{float(gap_lo):.15e}, {float(gap_hi):.15e}]")
    out(f"  our D, |Y|/D and C^2 D^2 - |Y|^2 overlap the note's printed 30-decimal enclosures: {'yes' if good else 'NO'}")
    return good


def finite(C, alpha, breaks, etas, ns, out=print):
    """Floating point, not a proof: |Y_n| / sum_F |P_m| for the construction at finite n (A = floor(tau n))."""
    a = complex(float(alpha[0]), float(alpha[1]))
    s = 1 - a
    tau = breaks[0]
    deltas = [complex(float(etas[0][0]), float(etas[0][1])) - s] + \
             [complex(float(etas[j][0] - etas[j - 1][0]), float(etas[j][1] - etas[j - 1][1])) for j in range(1, len(etas))]
    for n in ns:
        beta = [1 + 0j]
        for l in range(1, n + 1):
            beta.append(beta[-1] * (a + l - 1) / l)
        A = (tau * n).__floor__()
        mid = range(A + 1, n - A)
        w = [0j] * (n + 1)
        for m in mid:
            w[m] = sum(dj for dj, tj in zip(deltas, breaks) if m > (tj * n).__floor__())
        u = [w[m] / m if A < m < n - A else 0j for m in range(n + 1)]
        y = n * beta[n] - n * sum(w[m] * beta[n - m] / m for m in mid)
        dbl = 0j
        for k in range(2 * (A + 1), n + 1):
            conv = sum(u[x] * u[k - x] for x in range(max(A + 1, k - (n - A - 1)), min(n - A - 1, k - A - 1) + 1))
            dbl += conv * beta[n - k]
        y += n / 2 * dbl
        P = [n * beta[n - m] / m for m in range(n - A, n + 1)]
        y += s * sum(P)
        out(f"  n = {n}: |Y_n| / sum_F |P_m| = {abs(y) / sum(abs(p) for p in P):.10f}  (floating point; C = {float(C):.10f})")


def run(path, controls, fin):
    t0 = time.time()
    name, C, alpha, breaks, etas, d = load(path)
    print(f"{name}: {os.path.relpath(path, HERE)}; C = {C} = {float(C):.10f}; alpha = {float(alpha[0])} + "
          f"{float(alpha[1])}i; tau = {breaks[0]}; {len(etas)} middle block(s)")
    ok, res = evaluate(C, alpha, breaks, etas)
    if ok and name == "Gri26":
        ok &= compare_griego(d, res)
    if ok and name == "R2026a":
        ok &= compare_roehrig(C, res)
    if controls:
        print("Red controls (the same check, on altered data):")
        quiet = lambda *a, **k: None
        tests = []
        tests.append((f"C lowered to {float(Fraction(res['ratio']) - Fraction(1, 10 ** 10)):.10f}, just below our |Y|/D",
                      (Fraction(res["ratio"]) - Fraction(1, 10 ** 10), alpha, breaks, etas)))
        tests.append(("eta_0 replaced by s (the first step removed)",
                      (C, alpha, breaks, [(1 - alpha[0], -alpha[1])] + etas[1:])))
        tests.append(("the cut-off moved to tau = 3/10 (three corrections would then reach b_n)",
                      (C, alpha, [Fraction(3, 10)] + breaks[1:-1] + [Fraction(7, 10)], etas)))
        for what, args in tests:
            failed = []
            refused = not evaluate(*args, out=quiet, failed=failed)[0]
            print(f"  {what}: {'refused, as expected, by: ' + '; '.join(failed) if refused else 'ACCEPTED (control failed)'}")
            ok &= refused
    if fin:
        print("Finite n (floating point, not part of the proof):")
        finite(C, alpha, breaks, etas, (500, 1000, 2000, 4000))
    print(f"RESULT {name}: {'PASS' if ok else 'FAIL'}  ({time.time() - t0:.1f} s)")
    return ok


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    root = args[0] if args else os.path.join(HERE, "certs")
    rows = {"gri": os.path.join(root, "griego", "turan42_certificate.json"),
            "pin": os.path.join(root, "pinho", "c42_256.json"),
            "roe": os.path.join(root, "roehrig", "turan.json")}
    pick = [a.split("=", 1)[1].split(",") for a in sys.argv[1:] if a.startswith("--rows=")]
    names = pick[0] if pick else list(rows)
    ok = all([run(rows[n], "--controls" in sys.argv, "--finite" in sys.argv) for n in names])
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
