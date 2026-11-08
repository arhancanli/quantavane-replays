"""Replay of the 30a row [N2026], C_30 = gr(Av(1324)) >= 10.412263 (PR #313; Norton, arXiv 2608.20292v3, Theorem 1.1),
from the paper's LaTeX source alone: the evaluation behind Table 1 that the theorem rests on.

1. The strip functions. H(z,q) = 2/D with s = sqrt(1-4z), D = 2 - q + q s (eq. (H)), and H_j = Omega_j[H] with
   Omega_j[F] = (1/j!) d^j/dz^j (z^j F) (eq. (Omega)). We differentiate exactly, as rational functions in s and q with
   rational coefficients (z = (1 - s^2)/4, d/dz = (-2/s) d/ds), and compare with the closed forms (H2), (H3) printed in
   Section 5 by cross-multiplication.
2. The profile f = (beta, 3 - 4 alpha - 3 beta, 3 alpha + 2 beta - 2) (eq. (fprofile)) at the table's alpha, beta, and
   the hypotheses of Proposition 4.5 (transfer) on alpha, beta.
3. At z0 = 9604060/10^8, kappa = 516/1000: the bracket B(q) = -kappa log q + sum_j f_j log H_j(z0, q); its slope in
   w = log q at q = 1 and 9 (a sign change, inside the radius 1/Q(z0)); B and dB/dw at the table's q0; the tangent
   bound Lambda = min B >= B(q0) - |dB/dw(q0)| max|log q - log q0| over [1, 9] (Lemma 4.1: B is convex in log q);
   Phi(z0, 1, kappa) = 2 log(27 z0/4) + 2 Lambda - kappa log Q(z0) > 0 (Section 5.1); and 1/z0 > 10.412263 exactly.
Interval arithmetic of our own (Python decimal, 60 digits, outward rounding; exp, ln and sqrt are correctly rounded in
decimal and every result is widened by one unit in the last place).

Not checked: the mathematics that turns Phi > 0 into the bound (the second relaxation, the Harris separation, the
transfer proposition, BBEP's results it cites).

Our own code, standard library only. Usage: python3 -I norton_check.py CERTS_ROOT [--controls]
"""
import hashlib
import os
import sys
from decimal import Decimal as D, Context, ROUND_FLOOR, ROUND_CEILING, ROUND_HALF_EVEN
from fractions import Fraction as Fr
from math import factorial

PREC = 60
CF = Context(prec=PREC, rounding=ROUND_FLOOR)
CC = Context(prec=PREC, rounding=ROUND_CEILING)
CN = Context(prec=PREC + 5, rounding=ROUND_HALF_EVEN)


class Iv:
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        lo = D(lo) if not isinstance(lo, Fr) else frac_lo(lo)
        hi = lo if hi is None else (D(hi) if not isinstance(hi, Fr) else frac_hi(hi))
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
        return "[%s, %s]" % (fmt(self.lo), fmt(self.hi))


def frac_lo(q):
    return CF.divide(D(q.numerator), D(q.denominator))


def frac_hi(q):
    return CC.divide(D(q.numerator), D(q.denominator))


def ivq(q):
    q = Fr(q)
    return Iv(frac_lo(q), frac_hi(q))


def _down(r):
    return CF.plus(CN.next_minus(r))


def _up(r):
    return CC.plus(CN.next_plus(r))


def iln(a):
    if a.lo <= 0:
        raise ValueError("log of a non-positive interval")
    return Iv(_down(CN.ln(a.lo)), _up(CN.ln(a.hi)))


def isqrt(a):
    if a.lo < 0:
        raise ValueError("sqrt of a negative interval")
    return Iv(_down(CN.sqrt(a.lo)), _up(CN.sqrt(a.hi)))


def fmt(x, digits=22):
    return format(Context(prec=digits).plus(D(x)), "f")


# ------------------------------------------------------------- polynomials in s and q with rational coefficients
def padd(*ps):
    out = {}
    for p in ps:
        for k, v in p.items():
            out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if v != 0}


def pmul(p, r):
    out = {}
    for (i, j), v in p.items():
        for (k, l), w in r.items():
            key = (i + k, j + l)
            out[key] = out.get(key, 0) + v * w
    return {k: v for k, v in out.items() if v != 0}


def pscale(p, c):
    return {k: v * c for k, v in p.items() if v * c != 0}


def pds(p):
    return {(i - 1, j): v * i for (i, j), v in p.items() if i > 0}


def pdq(p):
    return {(i, j - 1): v * j for (i, j), v in p.items() if j > 0}


def ppow(p, n):
    out = {(0, 0): Fr(1)}
    for _ in range(n):
        out = pmul(out, p)
    return out


S1 = {(1, 0): Fr(1)}
Q1 = {(0, 1): Fr(1)}
DD = {(0, 0): Fr(2), (0, 1): Fr(-1), (1, 1): Fr(1)}        # D = 2 - q + q s
ZZ = {(0, 0): Fr(1, 4), (2, 0): Fr(-1, 4)}                 # z = (1 - s^2)/4


def rd_z(F):
    """d/dz of P/(s^a D^b) = (-2/s) d/ds: -2 (P_s s D - a P D - b q s P) / (s^(a+2) D^(b+1))."""
    P, a, b = F
    num = padd(pmul(pmul(pds(P), S1), DD), pscale(pmul(P, DD), -a), pscale(pmul(pmul(Q1, S1), P), -b))
    return (pscale(num, -2), a + 2, b + 1)


def rd_q(F):
    """d/dq of P/(s^a D^b) = (P_q D - b P (s - 1)) / (s^a D^(b+1))."""
    P, a, b = F
    num = padd(pmul(pdq(P), DD), pscale(pmul(P, {(1, 0): Fr(1), (0, 0): Fr(-1)}), -b))
    return (num, a, b + 1)


def omega(j, F):
    P, a, b = F
    G = (pmul(ppow(ZZ, j), P), a, b)
    for _ in range(j):
        G = rd_z(G)
    return (pscale(G[0], Fr(1, factorial(j))), G[1], G[2])


def same(F, G):
    """P/(s^a D^b) == R/(s^c D^d) as rational functions: P s^c D^d == R s^a D^b."""
    (P, a, b), (R, c, d) = F, G
    return pmul(pmul(P, ppow(S1, c)), ppow(DD, d)) == pmul(pmul(R, ppow(S1, a)), ppow(DD, b))


def poly_from(terms_s_by_q):
    """{k: [c0, c1, ...]} meaning sum_k q^k sum_i c_i s^i."""
    out = {}
    for k, cs in terms_s_by_q.items():
        for i, c in enumerate(cs):
            if c:
                out[(i, k)] = Fr(c)
    return out


# the closed forms printed in Section 5.1, (H2) and (H3), coefficients of s^0, s^1, ... for each power of q
H2_PAPER = (padd(poly_from({2: [-1, 3, -6, 10, -9, 3]}), poly_from({1: [2, 0, 12, -32, 18]}),
                 poly_from({0: [0, 0, 0, 32]})), 3, 3)
H2_DEN = Fr(1, 4)
H3_NUM = padd(poly_from({3: [1, -4, 8, -12, 18, -28, 32, -20, 5]}),
              pscale(poly_from({2: [-1, 2, -3, 6, -15, 30, -29, 10]}), 4),
              pscale(poly_from({1: [1, 0, 3, 0, 15, -48, 29]}), 4),
              poly_from({0: [0, 0, 0, 0, 0, 128]}))
H3_DEN = Fr(1, 8)

TEX = [
    r"H_2 \;=\; \frac{q^{2}\bigl(3s^{5} - 9s^{4} + 10s^{3} - 6s^{2} + 3s - 1\bigr) \;+\; 2q\bigl(9s^{4} - 16s^{3} + 6s^{2} + 1\bigr) \;+\; 32s^{3}} {4\,s^{3}D^{3}}",
    r"& q^{3}\bigl(5s^{8} - 20s^{7} + 32s^{6} - 28s^{5} + 18s^{4} - 12s^{3} + 8s^{2} - 4s + 1\bigr) \\ &+\; 4q^{2}\bigl(10s^{7} - 29s^{6} + 30s^{5} - 15s^{4} + 6s^{3} - 3s^{2} + 2s - 1\bigr) \\ &+\; 4q\bigl(29s^{6} - 48s^{5} + 15s^{4} + 3s^{2} + 1\bigr) \;+\; 128s^{5}\;\Bigr] .",
    r"H_3 \;=\; \frac{1}{8\,s^{5}D^{4}}\Bigl[\;",
    r"\Omega_j[z^n] \;=\; \binom{n+j}{j} z^n, \qquad\text{equivalently}\qquad \Omega_j[F(z)] \;=\; \frac{1}{j!} \frac{\partial^j}{\partial z^j}\bigl(z^j F(z)\bigr).",
    r"H(z,q) \;=\; \frac{1}{1 - q\Qc(z)} \;=\; \frac{2}{2 - q\bigl(1 - \sqrt{1-4z}\bigr)},",
    r"f_0 = \beta, \qquad f_2 = 3 - 4\alpha - 3\beta, \qquad f_3 = 3\alpha + 2\beta - 2 ,",
    r"\Phi(z,1,\kappa) \;=\; 2\log\frac{27z}{4} \;+\; 2\,\Lam(z,f,\kappa) \;-\; \kappa\log\Qc(z),",
    r"B(q) \;=\; -\kappa\log q \;+\!\! \sum_{j \in \{0,2,3\}}\!\! f_j\log H_j(z,q).",
    r"evaluation point $z_0$ & $9604060/10^{8}$ \\",
    r"component scale $\kappa$ & $516/1000$ \\",
    r"leaf density $\alpha$ & $55555555/10^{8}$ \\",
    r"empty-strip density $\beta$ & $18518518/10^{8}$ \\",
    r"tangent point $q_0$ & $3.046850759320$ \\",
    r"$\gr(\Av(1324)) \ge 10.412263$.",
]


def check_tex(root, out):
    path = os.path.join(root, "n2026/paper.tex")
    with open(path, "rb") as fh:
        raw = fh.read()
    out("n2026/paper.tex sha256 %s" % hashlib.sha256(raw).hexdigest())
    text = " ".join(raw.decode("utf-8").split())
    ok = True
    for t in TEX:
        hit = " ".join(t.split()) in text
        ok &= hit
        out("  our transcription appears verbatim: %s   %s" % (hit, t[:100]))
    return ok


def evalF(F, s, q):
    """Evaluate P/(s^a D^b) at interval s, q."""
    P, a, b = F
    acc = Iv(0)
    for (i, j), c in P.items():
        term = ivq(c)
        for _ in range(i):
            term = term * s
        for _ in range(j):
            term = term * q
        acc = acc + term
    Dv = Iv(2) - q + q * s
    den = Iv(1)
    for _ in range(a):
        den = den * s
    for _ in range(b):
        den = den * Dv
    return acc / den


def setup(z0, kappa, alpha, beta, Hs):
    s = isqrt(Iv(1) - Iv(4) * ivq(z0))
    Qz = (Iv(1) - s) / Iv(2)
    f = {0: beta, 2: 3 - 4 * alpha - 3 * beta, 3: 3 * alpha + 2 * beta - 2}
    return s, Qz, f


def bracket(Hs, dHs, f, kappa, s, q):
    """B(q) and dB/dw = q B'(q) at an interval q."""
    B = -ivq(kappa) * iln(q)
    Bw = -ivq(kappa)
    vals = {}
    for j in (0, 2, 3):
        h = evalF(Hs[j], s, q)
        dh = evalF(dHs[j], s, q)
        vals[j] = h
        B = B + ivq(f[j]) * iln(h)
        Bw = Bw + ivq(f[j]) * q * dh / h
    return B, Bw, vals


def phi_lower(z0, kappa, alpha, beta, q0, Hs, dHs, out=None):
    s, Qz, f = setup(z0, kappa, alpha, beta, Hs)
    q0i = ivq(q0)
    B0, Bw0, vals = bracket(Hs, dHs, f, kappa, s, q0i)
    _, Bw1, _ = bracket(Hs, dHs, f, kappa, s, Iv(1))
    _, Bw9, _ = bracket(Hs, dHs, f, kappa, s, Iv(9))
    radius = Iv(1) / Qz
    sign_change = Bw1.hi < 0 < Bw9.lo and Iv(9).hi < radius.lo
    lq0 = iln(q0i)
    span = Iv(max(lq0.hi, (iln(Iv(9)) - lq0).hi))
    absBw = max(abs(Bw0.lo), abs(Bw0.hi))
    lam_lo = (B0 - Iv(absBw) * span).lo
    phi_lo = (Iv(2) * iln(Iv(27) * ivq(z0) / Iv(4)) + Iv(2) * Iv(lam_lo) - ivq(kappa) * iln(Qz)).lo
    if out:
        out("s = sqrt(1 - 4 z0) in %s" % s)
        out("Q(z0) in %s; radius 1/Q(z0) in %s" % (Qz, radius))
        out("profile f = (f0, f2, f3) = (%s, %s, %s)" % (f[0], f[2], f[3]))
        for j in (0, 2, 3):
            out("H_%d(z0, q0) in %s" % (j, vals[j]))
        out("B'(1) = dB/dw at q = 1 in %s; B'(9) = (dB/dw)/9 in %s (B' changes sign on [1, 9], inside the radius: "
            "%s)" % (Bw1, Bw9 / Iv(9), sign_change))
        out("B(q0) in %s; dB/dw(q0) in %s" % (B0, Bw0))
        out("max |log q - log q0| over [1, 9] <= %s" % fmt(span.hi, 12))
        out("Lambda = min B >= %s" % fmt(lam_lo, 16))
        out("Phi(z0, 1, kappa) >= %s" % fmt(phi_lo, 12))
    return phi_lo, sign_change, B0, Bw0, vals, s, Qz, f


def phi_upper(z0, kappa, alpha, beta, q, Hs, dHs):
    """Lambda <= B(q) at any q in the range, so Phi <= 2 log(27 z/4) + 2 B(q) - kappa log Q."""
    s, Qz, f = setup(z0, kappa, alpha, beta, Hs)
    B, _, _ = bracket(Hs, dHs, f, kappa, s, ivq(q))
    return (Iv(2) * iln(Iv(27) * ivq(z0) / Iv(4)) + Iv(2) * B - ivq(kappa) * iln(Qz)).hi


def float_free_entropy(g, k):
    """E(gamma, kappa) = (g+k) log(g+k) - g log g - k log k (eq. (ent)), as an interval."""
    gi, ki = ivq(g), ivq(k)
    return (gi + ki) * iln(gi + ki) - gi * iln(gi) - ki * iln(ki)


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    root = sys.argv[1]
    controls = "--controls" in sys.argv

    def out(x):
        print(x)

    out("30a row [N2026] (PR #313): gr(Av(1324)) >= 10.412263 (Norton, arXiv 2608.20292v3, Theorem 1.1)")
    status = check_tex(root, out)
    out("sources contain our transcription: %s" % status)

    # 1. strip functions, exactly
    H = ({(0, 0): Fr(2)}, 0, 1)
    Hs = {0: H, 2: omega(2, H), 3: omega(3, H)}
    ok0 = same(omega(0, H), H)
    ok2 = same(Hs[2], (pscale(H2_PAPER[0], H2_DEN), 3, 3))
    ok3 = same(Hs[3], (pscale(H3_NUM, H3_DEN), 5, 4))
    out("Omega_0[H] = H: %s; Omega_2[H] equals the printed (H2): %s; Omega_3[H] equals the printed (H3): %s"
        % (ok0, ok2, ok3))
    status &= ok0 and ok2 and ok3
    dHs = {j: rd_q(Hs[j]) for j in Hs}

    # 2. the profile and the hypotheses of the transfer proposition
    z0, kappa = Fr(9604060, 10 ** 8), Fr(516, 1000)
    alpha, beta = Fr(55555555, 10 ** 8), Fr(18518518, 10 ** 8)
    q0 = Fr("3.046850759320")
    f = (beta, 3 - 4 * alpha - 3 * beta, 3 * alpha + 2 * beta - 2)
    prof = f == (Fr("0.18518518"), Fr("0.22222226"), Fr("0.03703701"))
    hyp = Fr(11, 20) <= alpha < Fr(5, 9) and Fr(7, 40) <= beta < Fr(5, 27) and 3 * alpha + 2 * beta > 2 \
        and 0 < z0 < Fr(1, 4) and kappa > 0
    out("profile from (fprofile): f = (%s, %s, %s), as printed: %s" % (f[0], f[1], f[2], prof))
    out("alpha in [11/20, 5/9), beta in [7/40, 5/27), 3 alpha + 2 beta = %s > 2, z0 in (0, 1/4): %s"
        % (3 * alpha + 2 * beta, hyp))
    status &= prof and hyp

    # 3. the evaluation
    phi_lo, sc, B0, Bw0, vals, s, Qz, fv = phi_lower(z0, kappa, alpha, beta, q0, Hs, dHs, out)
    def trunc_of(v, printed):
        p = D(printed)
        ulp = D(1).scaleb(p.as_tuple().exponent)
        return p <= v.lo and v.hi < p + ulp
    checks = [("Q(z0) = 0.107623395192832", trunc_of(Qz, "0.107623395192832")),
              ("H_0(z0, q0) = 1.487901331259", trunc_of(vals[0], "1.487901331259")),
              ("H_2(z0, q0) = 3.725658932481", trunc_of(vals[2], "3.725658932481")),
              ("H_3(z0, q0) = 6.161763006829", trunc_of(vals[3], "6.161763006829")),
              ("B(q0) > -0.1416704254", B0.lo > D("-0.1416704254")),
              ("|dB/dw(q0)| < 1e-13", max(abs(Bw0.lo), abs(Bw0.hi)) < D("1e-13"))]
    for name, ok in checks:
        out("Table 1, %s (decimals read as truncations, as the caption says): %s" % (name, ok))
    h3_rounded = (not checks[3][1]) and trunc_of(vals[3], "6.161763006828")
    out("  note: H_3(z0, q0) = %s..., so the printed 6.161763006829 is rounded, not truncated: %s"
        % (fmt(vals[3].lo, 16), h3_rounded))
    table = all(ok for name, ok in checks if not name.startswith("H_3")) and h3_rounded
    pos = phi_lo > 0
    bound = 10412263 * 9604060 < 10 ** 14 and Fr(1) / z0 > Fr("10.412263")
    out("Phi(z0, 1, kappa) > 0: %s (the table's 1.2762e-7 bound: %s); 1/z0 = %s... > 10.412263: %s"
        % (pos, phi_lo >= D("1.2762e-7"), fmt(D(1) / D("0.09604060"), 16), bound))
    status &= sc and table and pos and bound

    if controls:
        out("")
        out("controls (each must come out as stated):")
        ctl = True
        # (i) a perturbed closed form is refused
        bad = padd(H2_PAPER[0], {(3, 0): Fr(-1)})
        r1 = not same(Hs[2], (pscale(bad, H2_DEN), 3, 3))
        out("[ctl] (H2) with 31 s^3 in place of 32 s^3 is refused: %s" % r1)
        ctl &= r1
        # (ii) a stronger claim fails: at z = 0.0960405 (1/z = 10.41227...) Phi has a negative upper bound
        zb = Fr(960405, 10 ** 7)
        ub = phi_upper(zb, kappa, alpha, beta, q0, Hs, dHs)
        out("[ctl] z = 0.0960405 (1/z = %s...): Phi <= %s < 0, so the claim 10.41227 is refused: %s"
            % (fmt(D(1) / D("0.0960405"), 10), fmt(ub, 8), ub < 0))
        ctl &= ub < 0
        # (iii) the corner (alpha, beta) = (11/20, 7/40) has f_3 = 0 and fails the hypothesis 3 alpha + 2 beta > 2
        a2, b2 = Fr(11, 20), Fr(7, 40)
        r3 = not (3 * a2 + 2 * b2 > 2)
        out("[ctl] (alpha, beta) = (11/20, 7/40): 3 alpha + 2 beta = %s, hypothesis refused: %s" % (3 * a2 + 2 * b2, r3))
        ctl &= r3
        # (iv) BBEP's Theorem 9 exactly: at z = 8/81, Q = 1/9, and 2 log(27z/4) + 2 E(1, 1/2) + (1/2) log Q = 0,
        #      i.e. (2/3)^2 ((3/2)^(3/2) 2^(1/2))^2 (1/9)^(1/2) = 1
        s81 = isqrt(Iv(1) - Iv(4) * ivq(Fr(8, 81)))
        Q81 = (Iv(1) - s81) / Iv(2)
        r4 = frac_lo(Fr(1, 9)) <= Q81.hi and Q81.lo <= frac_hi(Fr(1, 9))
        prod = Fr(2, 3) ** 2 * (Fr(27, 8) * 2) * Fr(1, 3)
        out("[ctl] BBEP Theorem 9 reading: Q(8/81) in %s contains 1/9: %s; (2/3)^2 (27/8)(2)(1/3) = %s"
            % (Q81, r4, prod))
        ctl &= r4 and prod == 1
        # (v) the tangent point moved to q = 2.9: the tangent bound alone no longer certifies Phi > 0
        plo2 = phi_lower(z0, kappa, alpha, beta, Fr(29, 10), Hs, dHs)[0]
        out("[ctl] with the tangent taken at q = 2.9 the bound is Phi >= %s (no longer > 0: %s)" % (
            fmt(plo2, 8), plo2 <= 0))
        ctl &= plo2 <= 0
        # (vi) BBEP's own exponent (eq. (phi)) at their published (gamma, kappa) = (0.951509, 0.496339) with the
        #      equitable profile at (5/9, 5/27), saddle q solving sum f_j q dH_j/H_j = kappa: Norton's Section 5.2
        #      reports 1/z* = 10.2710129282 and q0 = 2.9170545 (BBEP: 10.27101292824530 and 2.917054)
        g, k2 = Fr("0.951509"), Fr("0.496339")
        fb = {0: Fr(5, 27), 2: Fr(2, 9), 3: Fr(1, 27)}
        ent = float_free_entropy(g, k2)

        def mid(v):
            return (v.lo + v.hi) / 2

        def saddle(zz):
            sz = isqrt(Iv(1) - Iv(4) * ivq(zz))
            lo, hi = D("1.0"), mid(Iv(1) / ((Iv(1) - sz) / Iv(2))) * D("0.999999")
            for _ in range(90):
                m = (lo + hi) / 2
                qm = Iv(m)
                gsum = -ivq(k2)
                for j in (0, 2, 3):
                    gsum = gsum + ivq(fb[j]) * qm * evalF(dHs[j], sz, qm) / evalF(Hs[j], sz, qm)
                if mid(gsum) < 0:
                    lo = m
                else:
                    hi = m
            return sz, (lo + hi) / 2

        def phi_bbep(zz):
            sz, qs = saddle(zz)
            qi = Iv(qs)
            val = (Iv(1) + ivq(g)) * iln(Iv(27) * ivq(zz) / Iv(4)) - ivq(k2) * iln(qi) + ent
            for j in (0, 2, 3):
                val = val + ivq(fb[j]) * iln(evalF(Hs[j], sz, qi))
            return mid(val), qs

        zlo, zhi = Fr("0.0973"), Fr("0.0975")
        for _ in range(45):
            zm = (zlo + zhi) / 2
            if phi_bbep(zm)[0] < 0:
                zlo = zm
            else:
                zhi = zm
        zs = (zlo + zhi) / 2
        qs = phi_bbep(zs)[1]
        inv = D(zs.denominator) / D(zs.numerator)
        r6 = abs(inv - D("10.2710129282")) < D("1e-9") and abs(qs - D("2.9170545")) < D("1e-6")
        out("[ctl] BBEP's exponent at their (gamma, kappa): 1/z* = %s, saddle q = %s (Norton: 10.2710129282, "
            "2.9170545): %s" % (fmt(inv, 14), fmt(qs, 9), r6))
        ctl &= r6
        out("controls all as expected: %s" % ctl)
        status &= ctl

    out("RESULT: %s" % ("PASS (Table 1 re-derived: Phi(z0, 1, kappa) > 0 at z0 = 0.09604060, so 1/z0 > 10.412263; "
                        "the mathematics of Sections 3 and 4 is not checked)" if status else "FAIL"))
    return 0 if status else 1


if __name__ == "__main__":
    sys.exit(main())
