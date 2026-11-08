"""Replay of the 18a row [GT2009], C_18 >= 1.46601 (Green and Tao, arXiv math/0703668, Section 5: a Hamming ball
A = H_{alpha n}(n) has (1/n) log2 sigma[A] = h(2 alpha) - h(alpha) + o(1) and needs
(1/n) log2 tau[A] = (1 - h(alpha)) h(alpha/(1 - h(alpha))) + o(1) translates; the paper maximizes the ratio
N(alpha)/D(alpha) over alpha in (0, 1/10] with Mathematica and reports 1.46601... at alpha = 0.0939288...).

We enclose the maximum of that ratio rigorously, in interval arithmetic of our own (Python decimal, 50 digits,
outward rounding; ln is correctly rounded in decimal and every result is widened by one unit in the last place):
a lower bound at alpha0 = 0.0939288 (which is what the row needs), and an upper bound on [0.05, 0.1] from cells, using
that D is increasing there (D' > 0 for alpha < (2 - sqrt 2)/4), that 1 - h is decreasing and that h(alpha/(1 - h(alpha)))
is increasing. The size formulas need 2 alpha <= 1/2 and alpha/(1 - h(alpha)) <= 1/2, which we check.
Not checked: the paper's lemma that an SMD covered by t translates of an m-dimensional subspace is covered by t
translates of F_2^m, and the passage from this example to C_18.

Usage: python3 -I hamming_check.py CERTS_ROOT [--controls]
"""
import hashlib
import os
import sys
from decimal import Decimal as D, Context, ROUND_FLOOR, ROUND_CEILING, ROUND_HALF_EVEN

P = 50
CF = Context(prec=P, rounding=ROUND_FLOOR)
CC = Context(prec=P, rounding=ROUND_CEILING)
CN = Context(prec=P + 5, rounding=ROUND_HALF_EVEN)
LN2 = (CF.plus(CN.next_minus(CN.ln(D(2)))), CC.plus(CN.next_plus(CN.ln(D(2)))))


def ln_lo(x):
    return CF.plus(CN.next_minus(CN.ln(x)))


def ln_hi(x):
    return CC.plus(CN.next_plus(CN.ln(x)))


def h_iv(lo, hi):
    """enclosure of the binary entropy h on [lo, hi], 0 < lo <= hi <= 1/2 (h is increasing there)"""
    def h_lo(x):  # lower bound of h(x) = -(x ln x + (1-x) ln(1-x))/ln 2
        t = CC.add(CC.multiply(x, ln_hi(x)), CC.multiply(CC.subtract(1, x), ln_hi(CC.subtract(1, x))))
        return CF.divide(CF.minus(t), LN2[1])

    def h_hi(x):
        t = CF.add(CF.multiply(x, ln_lo(x)), CF.multiply(CF.subtract(1, x), ln_lo(CF.subtract(1, x))))
        return CC.divide(CC.minus(t), LN2[0])
    return h_lo(lo), h_hi(hi)


def ratio_point(a):
    """an enclosure of N(a)/D(a) at an exact point a"""
    ha = h_iv(a, a)
    m_lo, m_hi = CF.subtract(1, ha[1]), CC.subtract(1, ha[0])           # 1 - h(a)
    u_lo, u_hi = CF.divide(a, m_hi), CC.divide(a, m_lo)                 # a/(1 - h(a))
    hu = h_iv(u_lo, u_hi)
    N = (CF.multiply(m_lo, hu[0]), CC.multiply(m_hi, hu[1]))
    h2 = h_iv(2 * a, 2 * a)
    Dv = (CF.subtract(h2[0], ha[1]), CC.subtract(h2[1], ha[0]))
    return CF.divide(N[0], Dv[1]), CC.divide(N[1], Dv[0]), u_hi


def ratio_cell_upper(a, b):
    """an upper bound of N/D on [a, b] in (0, 1/10]: N <= (1 - h(a)) h(b/(1 - h(b))), D >= D(a)"""
    ha = h_iv(a, a)
    hb = h_iv(b, b)
    m_a_hi = CC.subtract(1, ha[0])
    u_b_hi = CC.divide(b, CF.subtract(1, hb[1]))
    hu = h_iv(u_b_hi, u_b_hi)
    N_hi = CC.multiply(m_a_hi, hu[1])
    h2a = h_iv(2 * a, 2 * a)
    D_lo = CF.subtract(h2a[0], ha[1])
    return CC.divide(N_hi, D_lo), u_b_hi


def sup_upper(lo, hi, cells, target, maxdepth=30):
    w = (hi - lo) / cells
    stack = [(lo + i * w, hi if i == cells - 1 else lo + (i + 1) * w, 0) for i in range(cells)]
    worst = D(0)
    maxu = D(0)
    n = 0
    while stack:
        a, b, d = stack.pop()
        ub, u = ratio_cell_upper(a, b)
        maxu = max(maxu, u)
        if ub <= target or d >= maxdepth:
            worst = max(worst, ub)
            n += 1
        else:
            m = (a + b) / 2
            stack += [(a, m, d + 1), (m, b, d + 1)]
    return worst, n, maxu


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    root = sys.argv[1]
    controls = "--controls" in sys.argv

    def out(x):
        print(x, flush=True)

    out("18a row [GT2009]: C_18 >= 1.46601 (Green and Tao, arXiv math/0703668, Section 5)")
    path = os.path.join(root, "gt2009/main.tex")
    with open(path, "rb") as fh:
        raw = fh.read()
    out("gt2009/main.tex sha256 %s" % hashlib.sha256(raw).hexdigest())
    text = " ".join(raw.decode("latin-1").split())
    needles = [r"\frac{1}{n}\log_2 \sigma[A] = h(2\alpha) - h(\alpha) + o(1).",
               r"\frac{1}{n}\log_2 \tau[A] = (1 - h(\alpha))h\big( \frac{\alpha}{1-h(\alpha)}\big) + o(1).",
               r"Let $\alpha \in (0,1/10]$ be a constant to be optimised later.",
               r"is largest when $\alpha = 0.0939288\ldots$, in which case it is $1.46601\ldots + o(1)$."]
    status = True
    for nd in needles:
        hit = " ".join(nd.split()) in text
        status &= hit
        out("  our transcription appears verbatim: %s   %s" % (hit, nd))
    a0 = D("0.0939288")
    lo0, hi0, u0 = ratio_point(a0)
    out("ratio at alpha = 0.0939288 in [%s, %s]; alpha/(1 - h(alpha)) = %s... <= 1/2 and 2 alpha <= 1/2: %s"
        % (str(lo0)[:22], str(hi0)[:22], str(u0)[:10], u0 <= D("0.5")))
    bound = lo0 >= D("1.46601")
    out("so the example gives the exponent %s... >= 1.46601: %s" % (str(lo0)[:14], bound))
    status &= bound and u0 <= D("0.5")
    worst, ncell, maxu = sup_upper(D("0.05"), D("0.1"), 500, D("1.46602"), maxdepth=24)
    local = worst <= D("1.46602") and maxu <= D("0.5")
    out("on [0.05, 0.1] the ratio is at most %s (%d cells; D is increasing there), so the paper's 'largest at "
        "0.0939288..., 1.46601...' is the truncation of a value in [%s, 1.46602): %s"
        % (str(worst)[:12], ncell, str(lo0)[:14], local))
    status &= local
    if controls:
        out("")
        out("controls (each must come out as stated):")
        ctl = True
        r1 = not (sup_upper(D("0.09"), D("0.1"), 50, D("1.466012"), maxdepth=14)[0] <= D("1.466012"))
        out("[ctl] an upper bound 1.466012 (below the value at alpha0) cannot be certified on [0.09, 0.1]: %s" % r1)
        ctl &= r1
        lo2, hi2, _ = ratio_point(D("0.12"))
        r2 = hi2 < D("1.46601")
        out("[ctl] alpha = 0.12 gives only %s... < 1.46601: %s" % (str(hi2)[:10], r2))
        ctl &= r2
        r3 = not (lo0 >= D("1.466013"))
        out("[ctl] the value 1.466013 is refused at alpha0: %s" % r3)
        ctl &= r3
        out("controls all as expected: %s" % ctl)
        status &= ctl
    out("RESULT: %s" % ("PASS (the Hamming-ball example gives 1.4660123711... >= 1.46601 at alpha = 0.0939288; the "
                        "largest value on [0.05, 0.1] is below 1.46602)" if status else "FAIL"))
    return 0 if status else 1


if __name__ == "__main__":
    sys.exit(main())
