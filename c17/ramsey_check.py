"""Replay of the 17a row [WMR2026], C_17 <= 3.6960839126* (PR #313), inside the theorem it relies on.

The row is exp(F(1)) for F(l) = h(l) + p(l) e^(-l), h(l) = (1+l)log(1+l) - l log l, with the quintic
p(l) = -0.25 l + 0.033 l^2 + 0.08 l^3 + 0 l^4 - 0.0778 l^5 of Wang et al. (arXiv 2603.15617v2, Appendix A.3),
which they report as passing the conditions of Gupta, Ndiaye, Norin and Wei (arXiv 2407.19026; Theorem 13 of v1,
Theorem 14 of v2): for every l in (0,1] there are M(l), Y(l) in (0,1) with
    X(l) = (1 - e^(-F'(l)))^(1/(1-M(l))) (1 - M(l)),   (X(l), Y(l)) in R,
    F(l) > -(1/2) (log X(l) + l log M(l) + l log Y(l)).
R is defined (both versions) by R(k,l) <= x^(-k) y^(-l) for ALL k, l with k + l large, so a bound
R(k,l) <= exp(f(l/k) k + o(k)) for l <= k puts (x, y) in R only if, for every s in (0,1],
    (a) f(s) + log x + s log y <= 0      (the pairs with l <= k)   and
    (b) f(s) + s log x + log y <= 0      (the pairs with k < l, through R(k,l) = R(l,k)).
The benchmark accepts (x, y) if (a) holds for (x, y) or for (y, x) with f = U (its R_0); that is one of (a), (b),
not both. Here we check l = 1 exactly in interval arithmetic of our own (Python decimal, 60 digits, outward rounding;
exp and ln are correctly rounded in decimal, and every result is widened by one unit in the last place).

Our own code, standard library only, from the two papers' TeX sources (read as data, compared with our
transcription before anything is computed). Usage: python3 -I ramsey_check.py CERTS_ROOT [--controls]
"""
import os
import sys
from decimal import Decimal as D, Context, ROUND_FLOOR, ROUND_CEILING, ROUND_HALF_EVEN
import hashlib

PREC = 60
CF = Context(prec=PREC, rounding=ROUND_FLOOR)
CC = Context(prec=PREC, rounding=ROUND_CEILING)
CN = Context(prec=PREC + 5, rounding=ROUND_HALF_EVEN)


class Iv:
    """A closed interval [lo, hi] of Decimals; every operation rounds outward."""
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        lo = D(lo)
        hi = lo if hi is None else D(hi)
        if lo > hi:
            raise ValueError("empty interval")
        self.lo, self.hi = lo, hi

    def __add__(a, b):
        b = b if isinstance(b, Iv) else Iv(b)
        return Iv(CF.add(a.lo, b.lo), CC.add(a.hi, b.hi))

    __radd__ = __add__

    def __neg__(a):
        return Iv(-a.hi, -a.lo)

    def __sub__(a, b):
        b = b if isinstance(b, Iv) else Iv(b)
        return Iv(CF.subtract(a.lo, b.hi), CC.subtract(a.hi, b.lo))

    def __rsub__(a, b):
        return Iv(b) - a

    def __mul__(a, b):
        b = b if isinstance(b, Iv) else Iv(b)
        lo = min(CF.multiply(x, y) for x in (a.lo, a.hi) for y in (b.lo, b.hi))
        hi = max(CC.multiply(x, y) for x in (a.lo, a.hi) for y in (b.lo, b.hi))
        return Iv(lo, hi)

    __rmul__ = __mul__

    def __truediv__(a, b):
        b = b if isinstance(b, Iv) else Iv(b)
        if b.lo <= 0 <= b.hi:
            raise ZeroDivisionError("interval contains 0")
        lo = min(CF.divide(x, y) for x in (a.lo, a.hi) for y in (b.lo, b.hi))
        hi = max(CC.divide(x, y) for x in (a.lo, a.hi) for y in (b.lo, b.hi))
        return Iv(lo, hi)

    def __rtruediv__(a, b):
        return Iv(b) / a

    def __repr__(self):
        return "[%s, %s]" % (fmt(self.lo, 25), fmt(self.hi, 25))


def _down(r):
    return CF.plus(CN.next_minus(r))


def _up(r):
    return CC.plus(CN.next_plus(r))


def iexp(a):
    return Iv(_down(CN.exp(a.lo)), _up(CN.exp(a.hi)))


def iln(a):
    if a.lo <= 0:
        raise ValueError("log of a non-positive interval")
    return Iv(_down(CN.ln(a.lo)), _up(CN.ln(a.hi)))


def fmt(x, digits=20):
    return format(Context(prec=digits).plus(D(x)), "f") if abs(D(x)) < 10 ** 6 else str(x)


LN2 = iln(Iv(2))
EINV = iexp(Iv(-1))


def poly(coeffs, s):
    """p(s) = sum_j coeffs[j] s^(j+1) (no constant term), Horner, interval s."""
    acc = Iv(0)
    for c in reversed(coeffs):
        acc = acc * s + Iv(c)
    return acc * s


def dpoly(coeffs, s):
    """p'(s)."""
    acc = Iv(0)
    n = len(coeffs)
    for j in range(n - 1, -1, -1):
        acc = acc * s + Iv(D(j + 1) * coeffs[j])
    return acc


def h_point(s):
    """h(s) = (1+s) ln(1+s) - s ln s at an exact point s in (0, 1]."""
    s = Iv(s)
    return (Iv(1) + s) * iln(Iv(1) + s) - s * iln(s)


def f_point(coeffs, s):
    """f(s) = h(s) + p(s) e^(-s) at an exact point."""
    si = Iv(s)
    return h_point(s) + poly(coeffs, si) * iexp(-si)


def f_cell(coeffs, a, b):
    """An enclosure of f on the cell [a, b] (0 < a < b <= 1): h is increasing, p e^(-s) by interval extension."""
    s = Iv(a, b)
    hpart = Iv(h_point(a).lo, h_point(b).hi)
    return hpart + poly(coeffs, s) * iexp(-s)


def F1(coeffs):
    """F(1) = 2 ln 2 + p(1)/e."""
    return Iv(2) * LN2 + poly(coeffs, Iv(1)) * EINV


def Fp1(coeffs):
    """F'(1) = ln 2 + (p'(1) - p(1))/e."""
    one = Iv(1)
    return LN2 + (dpoly(coeffs, one) - poly(coeffs, one)) * EINV


def logX(lnP, mu):
    """log X = log(1-mu) + log P/(1-mu) at an exact mu, with log P an interval (log P < 0)."""
    om = Iv(1) - Iv(mu)
    return iln(om) + lnP / om


# ----------------------------------------------------------------------------------------------- the sources
QUINTIC = [D("-0.25"), D("0.033"), D("0.08"), D("0.0"), D("-0.0778")]   # Appendix A.3 of [WMR2026]
U_COEFFS = [D("-0.25"), D("0.033"), D("0.08")]                          # the benchmark's R_0 rate function
F003 = [D("-0.25"), D("0.03"), D("0.08")]                               # GNNW v2, Theorem 1 (F_{0.03})
F008 = [D("-0.25"), D("0.08"), D("0.08")]                               # GNNW v2, F_{0.08} (eq. base-bound)

TRANSCRIPTION = {
    "wmr2026/sections/appendix.tex": [
        r"with $(c_1,\; c_2,\; c_3,\; c_4,\; c_5) = (-0.25,\,0.033,\,0.08,\,0.0,\,-0.0778)$.",
        r"Upper bound base $c$ & $3.7992027396\ldots$ & $\mathbf{3.6960839126\ldots}$ \\",
        r"F(\lambda) = (1+\lambda)\log(1+\lambda) - \lambda\log\lambda + p(\lambda)\,e^{-\lambda},",
        r"U(\mu) = G(\mu) + (1+\mu)\log(1+\mu) - \mu\log\mu, \quad G(\mu) = (-0.25\mu + 0.033\mu^2 + 0.08\mu^3)\,e^{-\mu}.",
        r"the pair $(x,y)$ is accepted if either $(x,y) \in \mathcal{R}_0$ or $(y,x) \in \mathcal{R}_0$.",
        r"A pair $(x,y) \in \mathcal{R}_0$ iff $-\log x - \mu\log y \ge U(\mu)$ for all $\mu \in (0,1]$.",
        r"X(\lambda) = \bigl(1 - e^{-F'(\lambda)}\bigr)^{1/(1-M(\lambda))}\,(1 - M(\lambda)).",
        r"\item $F(\lambda) > -\tfrac{1}{2}\bigl(\log X(\lambda) + \lambda\log M(\lambda) + \lambda\log Y(\lambda)\bigr)$",
    ],
    "gnnw_v2/main.tex": [
        r"R(k,\ell) \leq x^{-k}y^{-\ell}  \qquad \mathrm{for \; all} \; k, \ell \in \bb{N} \; \mathrm{such\; that \;} k+\ell \geq N .",
        r"& and  \qquad F(\lambda) > -\frac{1}{2} \s{\log X(\lambda)+ \lambda \log M(\lambda) + \lambda\log Y(\lambda)} \end{align*}",
        r"&F'(\lambda) > 0, \qquad  X(\lambda) = (1 - e^{-F'(\lambda)})^{\frac{1}{1-M(\lambda)}}(1-M(\lambda)), \qquad (X(\lambda),Y(\lambda)) \in \mc{R},",
        r"\bigl(-\tfrac14\lambda+b\lambda^2",
        r"\exp\!\left(-\frac23\lambda-\frac3{10}\lambda^2\right),",
        r"we expect that lowering the base of the exponent below $3.7$",
    ],
    "gnnw_v1/RamseyArXiv.tex": [
        r"such that $R(k,\ell) \leq x^{-k}y^{-\ell}$ for all $k+\ell \geq N_0$.",
        r"\begin{thm}\label{t:general}",
    ],
}


def sha256(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def check_sources(root, out):
    ok = True
    for rel, needles in TRANSCRIPTION.items():
        path = os.path.join(root, rel)
        with open(path, encoding="utf-8") as fh:
            text = " ".join(fh.read().split())
        out("%s sha256 %s" % (rel, sha256(path)))
        for nd in needles:
            hit = " ".join(nd.split()) in text
            ok &= hit
            out("  our transcription appears verbatim: %s   %s" % (hit, nd[:110]))
    return ok


# ----------------------------------------------------------------------------------------------- the checks
S_LIST = [D(k) / 100 for k in range(5, 101, 5)]   # the s values tried in (b); s = 1 is (a) and (b) at once


def refute_lambda1(F, region, out, tag, mu1=D("0.10"), mu2=D("0.70"), ncells=600, maxdepth=12):
    """Try to show: for every mu in (0,1) and every y certified by the bound `region` (both (a) and (b)),
    F(1) + (1/2)(log X(mu) + log mu + log y) < 0. Returns (refuted, worst upper bound, s values used)."""
    f1 = F1(F)
    lnP = iln(Iv(1) - iexp(-Fp1(F)))
    fs = {s: f_point(region, s) for s in S_LIST}
    used = set()
    worst = None

    def ub_cell(a, b, first=False):
        lxa = lnP if first else logX(lnP, a)        # sup of log X on the cell (log X decreases in mu)
        lnb = iln(Iv(b))
        best = f1 + (lxa + lnb) * Iv(D("0.5"))       # y <= 1
        best_s = None
        for s, fv in fs.items():                     # (b): log y <= -f(s) - s log x
            cand = f1 + ((Iv(1) - Iv(s)) * lxa + lnb - fv) * Iv(D("0.5"))
            if cand.hi < best.hi:
                best, best_s = cand, s
        return best.hi, best_s

    # (0, mu1]: log mu <= log mu1, log X <= log P, y <= 1
    u0, s0 = ub_cell(D(0), mu1, first=True)
    refuted = u0 < 0
    worst = u0
    out("[%s] cell (0, %s]: upper bound %s (s = %s)" % (tag, mu1, fmt(u0, 8), s0))
    # [mu2, 1): log mu <= 0, log X <= log X(mu2), y <= 1
    lx2 = logX(lnP, mu2)
    u1 = (f1 + lx2 * Iv(D("0.5"))).hi
    refuted &= u1 < 0
    worst = max(worst, u1)
    out("[%s] cell [%s, 1): upper bound %s (y <= 1)" % (tag, mu2, fmt(u1, 8)))
    stack = []
    w = (mu2 - mu1) / ncells
    for i in range(ncells):
        stack.append((mu1 + i * w, mu1 + (i + 1) * w, 0))
    ncell = 0
    open_cells = []
    while stack:
        a, b, d = stack.pop()
        u, s = ub_cell(a, b)
        if u < 0:
            ncell += 1
            worst = max(worst, u)
            if s is not None:
                used.add(s)
        elif d < maxdepth:
            m = (a + b) / 2
            stack.append((a, m, d + 1))
            stack.append((m, b, d + 1))
        else:
            open_cells.append((a, b, u))
    refuted &= not open_cells
    out("[%s] cells in [%s, %s] closed: %d; cells left open: %d%s" % (
        tag, mu1, mu2, ncell, len(open_cells),
        "" if not open_cells else " (first: [%s, %s], bound %s)" % (
            fmt(open_cells[0][0], 8), fmt(open_cells[0][1], 8), fmt(open_cells[0][2], 6))))
    return refuted, worst, sorted(used)


def region_member(region, x_hi, y, out, tag, ncells=2000, maxdepth=14):
    """Check (a) and (b) for every s in (0,1] at the point (x, y) for every x <= x_hi (monotone in x):
    max over s of f(s) + log x + s log y and of f(s) + s log x + log y must be <= 0."""
    lx = iln(Iv(x_hi)).hi
    ly = iln(Iv(y)).hi
    worst_a = worst_b = None
    stack = [(D(i) / ncells, D(i + 1) / ncells, 0) for i in range(ncells)]
    bad = 0
    while stack:
        a, b, d = stack.pop()
        if a == 0:
            # (0, b]: h increases from h(0+) = 0, so f(s) <= h(b) + max of p(s) e^(-s) on [0, b]
            fhi = CC.add(h_point(b).hi, (poly(region, Iv(0, b)) * iexp(-Iv(0, b))).hi)
            va = CC.add(fhi, lx)                          # s log y <= 0
            vb = CC.add(fhi, ly)                          # s log x <= 0
        else:
            fhi = f_cell(region, a, b).hi
            va = CC.add(CC.add(fhi, lx), CC.multiply(a, ly))   # s log y is largest at s = a (log y < 0)
            vb = CC.add(CC.add(fhi, CC.multiply(a, lx)), ly)   # s log x is largest at s = a (log x < 0)
        if va <= 0 and vb <= 0:
            worst_a = va if worst_a is None else max(worst_a, va)
            worst_b = vb if worst_b is None else max(worst_b, vb)
        elif d < maxdepth and a > 0:
            m = (a + b) / 2
            stack.append((a, m, d + 1))
            stack.append((m, b, d + 1))
        else:
            bad += 1
    out("[%s] (x, y) with x <= %s, y = %s: (a) max <= %s, (b) max <= %s, cells failing: %d" % (
        tag, fmt(x_hi, 12), fmt(y, 12), fmt(worst_a, 6) if worst_a is not None else "-",
        fmt(worst_b, 6) if worst_b is not None else "-", bad))
    return bad == 0


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    root = sys.argv[1]
    controls = "--controls" in sys.argv
    lines = []

    def out(s):
        print(s)
        lines.append(s)

    status = True
    out("17a row [WMR2026] (PR #313): C_17 <= 3.6960839126*")
    src_ok = check_sources(root, out)
    out("sources contain our transcription: %s" % src_ok)
    status &= src_ok

    # 1. the printed value
    f1 = F1(QUINTIC)
    c = iexp(f1)
    out("p(1) = %s exactly; F(1) = 2 log 2 + p(1)/e in %s" % (sum(QUINTIC), f1))
    out("exp(F(1)) = 4 exp(-0.2148/e) in %s" % c)
    trunc = D("3.6960839126") <= c.lo and c.hi < D("3.6960839127")
    out("printed 3.6960839126... is exp(F(1)) truncated at 10 decimals: %s" % trunc)
    status &= trunc

    # 2. the benchmark's rule at lambda = 1
    one = Iv(1)
    lnP = iln(Iv(1) - iexp(-Fp1(QUINTIC)))
    out("F'(1) in %s; log(1 - e^(-F'(1))) in %s" % (Fp1(QUINTIC), lnP))
    # U is increasing on (0,1]: U'(s) = log((1+s)/s) + (g'(s) - g(s)) e^(-s) >= log 2 - 1/4 > 0, because
    # g' - g = -0.25 + 0.316 s + 0.207 s^2 - 0.08 s^3 >= -0.25 on [0,1] (0.316 s - 0.08 s^3 >= 0.236 s >= 0)
    a1, a2, a3 = D("0.25") + D("2") * D("0.033"), D("3") * D("0.08") - D("0.033"), -D("0.08")
    mono = a1 == D("0.316") and a2 == D("0.207") and a3 == D("-0.08") and a1 + a3 >= 0 and a2 >= 0 \
        and (LN2 - Iv(D("0.25"))).lo > 0
    u1 = f_point(U_COEFFS, D(1))
    out("U is increasing on (0,1] (U' >= log 2 - 1/4 > 0): %s; U(1) in %s" % (mono, u1))
    mu0 = D("0.3225")
    lx0 = logX(lnP, mu0)
    any_y = mono and (lx0 + u1).hi < 0
    out("at mu = %s: log X in %s; log X + U(1) < 0, so (X, y) is in the benchmark's R_0 for EVERY y in (0,1): %s"
        % (mu0, lx0, any_y))
    y_near = D(1) - D("1e-9")
    slack_or = f1 + (lx0 + iln(Iv(mu0)) + iln(Iv(y_near))) * Iv(D("0.5"))
    out("so the benchmark's rule accepts lambda = 1: with M = %s, Y = 1 - 1e-9 the slack is in %s (> 0: %s)"
        % (mu0, slack_or, slack_or.lo > 0))
    status &= any_y and slack_or.lo > 0
    cap = iexp(-f_point(U_COEFFS, D("0.3")) - Iv(D("0.3")) * lx0)
    out("but (b) at s = 0.3 caps y at exp(-U(0.3) - 0.3 log X) in %s for that X" % cap)

    # 3. the theorem's conditions at lambda = 1, with the region a bound certifies (both (a) and (b))
    found = []
    for name, region in (("U, the benchmark's R_0 bound", U_COEFFS), ("F_0.03, GNNW Theorem 1", F003),
                         ("the quintic F itself", QUINTIC)):
        r, worst, used = refute_lambda1(QUINTIC, region, out, "quintic / " + name)
        out("[quintic / %s] for every M(1) in (0,1) and every Y(1) certified by this bound, condition (3) fails at "
            "lambda = 1: %s; left side minus right side <= %s; s values used: %s"
            % (name, r, fmt(worst, 6), ", ".join(str(s) for s in used)))
        found.append(r)
    finding = all(found)

    ctl_ok = True
    if controls:
        out("")
        out("controls (each must come out as stated):")
        # (i) GNNW's own diagonal step passes with a region certified by F_0.08 (their Lemma 16 (2), lambda = 1)
        M1 = iexp(Iv(D(-2) / 3) - Iv(D("0.3")))
        lnP3 = iln(Iv(1) - iexp(-Fp1(F003)))
        # log X at the interval M1: log(1-M) + log P/(1-M) is decreasing in M, so use the ends
        om = Iv(1) - M1
        lX1 = Iv((iln(Iv(om.lo)) + lnP3 / Iv(om.lo)).lo, (iln(Iv(om.hi)) + lnP3 / Iv(om.hi)).hi)
        X1hi = iexp(Iv(lX1.hi)).hi
        y0 = D("0.8113")   # just below Y_{F_0.08}(X_1) = 0.81158... (GNNW Lemma 15 frontier, our float estimate)
        mem = region_member(F008, X1hi, y0, out, "ctl GNNW F_0.08 region")
        sl = F1(F003) + (lX1 + iln(M1) + iln(Iv(y0))) * Iv(D("0.5"))
        ok1 = mem and sl.lo > 0
        out("[ctl GNNW F_0.03 at lambda = 1, M = exp(-2/3 - 3/10), Y = %s certified by F_0.08] slack in %s: passes: %s"
            % (y0, sl, ok1))
        ctl_ok &= ok1
        # (ii) the same refutation run on that passing instance must leave cells open
        r2, _, _ = refute_lambda1(F003, F008, out, "ctl F_0.03 / F_0.08", maxdepth=6)
        out("[ctl] the refutation does not close on GNNW's passing instance: %s" % (not r2))
        ctl_ok &= not r2
        # (iii) the quintic's point above violates (b): at s = 0.3, U(0.3) + 0.3 log X + log 0.95 > 0
        v = f_point(U_COEFFS, D("0.3")) + Iv(D("0.3")) * lx0 + iln(Iv(D("0.95")))
        out("[ctl] (X(0.3225), 0.95) violates (b) for U at s = 0.3: U(0.3) + 0.3 log X + log 0.95 in %s (> 0: %s)"
            % (v, v.lo > 0))
        ctl_ok &= v.lo > 0
        # (iv) printed value moved by one unit in the 10th decimal is refused as a truncation
        moved = D("3.6960839127")
        r4 = not (moved <= c.lo and c.hi < moved + D("1e-10"))
        out("[ctl] 3.6960839127 refused as the truncation: %s" % r4)
        ctl_ok &= r4
        out("controls all as expected: %s" % ctl_ok)
        status &= ctl_ok

    if status and finding:
        out("RESULT: FINDING. The printed value is exp(F(1)) for the quintic, truncated, but the quintic does not satisfy "
            "the conditions of GNNW (Theorem 13 of v1, 14 of v2) at lambda = 1 for any M(1) with Y(1) from the region "
            "U, F_0.03 or F itself certifies; the benchmark's PASS rests on accepting either orientation.")
        return 0
    out("RESULT: %s" % ("FAIL (a check did not come out as stated)" if not status else "no finding"))
    return 1


if __name__ == "__main__":
    sys.exit(main())
