"""42a: the explicit value C_42 >= 1/2 + 2^-4096* (R. Pinho, PR #310) for Biro's strict lower bound, replayed from
the written proof (LOWER_PROOF.md at rpinho/math-certificates @ 36325fa) and its parameter file alone (Quantavane
replays).  Standard library only; run with python3 -I.  The claimant's verify_lower_exact.py was neither read nor run.

What this checks.  The note proves, for every n and every complex tuple z_1 = 1, z_2, ..., z_n, that some power sum
S_l (1 <= l <= n) has |S_l| > q = 1/2 + delta, by Biro's argument made explicit (first escape of the b increments,
radial balance, the summed coefficients becoming parallel, transfer back to b).  We read the proof step by step and
re-derived each implication; this program checks, with exact arithmetic of our own, every numerical comparison and
algebraic identity those steps use:
  - every inequality in r = 1/sqrt(2) and c = cos(pi/8) is decided exactly in Q(sqrt 2) (numbers a + b sqrt 2 with
    rational a, b; the sign is exact), every square root by squaring;
  - the parameter chain (5)-(31) is evaluated from the parameter file (all powers of two), and every derived exponent
    the file states (t0 > 2^1021, k/t0 > 2^1022, b <= 2^-254, product < 2^-340, final error < 1/16) is shown to
    follow, with its slack;
  - the polynomial identities behind (11), (15), (23) and (29) are checked by exact expansion, and the coefficient
    identities (1) and (3) on random tuples of Gaussian rationals (with a perturbed tuple rejected).
--controls re-runs the same checks on altered parameter files, each of which must be refused.

What rests on mathematics (read, not computed): the geometric steps (projections, the first-escape angle, the law of
cosines), the Cauchy-Schwarz steps, the compactness argument of section 8.

usage: python3 -I biro_check.py CERTS_ROOT [--controls]
"""
import hashlib
import json
import os
import random
import sys
import time
from fractions import Fraction as F

CERT_SHA = "f6681b2f4bf138ad04b8cd62c585d9b13a50057caeda0bc676262a31a1a1c79b"
PROOF_SHA = "3acab5f81d45eb8d1db381355bc9333acf8bea4e251c04c06cd6cf518d92ba1a"


class Q2:
    """a + b sqrt(2), a and b rational; exact sign."""
    __slots__ = ("a", "b")

    def __init__(self, a, b=0):
        self.a, self.b = F(a), F(b)

    @staticmethod
    def of(x):
        return x if isinstance(x, Q2) else Q2(x)

    def __add__(self, o):
        o = Q2.of(o); return Q2(self.a + o.a, self.b + o.b)
    __radd__ = __add__

    def __sub__(self, o):
        o = Q2.of(o); return Q2(self.a - o.a, self.b - o.b)

    def __rsub__(self, o):
        return Q2.of(o) - self

    def __neg__(self):
        return Q2(-self.a, -self.b)

    def __mul__(self, o):
        o = Q2.of(o); return Q2(self.a * o.a + 2 * self.b * o.b, self.a * o.b + self.b * o.a)
    __rmul__ = __mul__

    def inv(self):
        n = self.a * self.a - 2 * self.b * self.b
        assert n != 0
        return Q2(self.a / n, -self.b / n)

    def __truediv__(self, o):
        return self * Q2.of(o).inv()

    def __rtruediv__(self, o):
        return Q2.of(o) * self.inv()

    def sign(self):
        a, b = self.a, self.b
        if a >= 0 and b >= 0:
            return 0 if (a == 0 and b == 0) else 1
        if a <= 0 and b <= 0:
            return -1
        if a > 0:  # b < 0: compare a with |b| sqrt 2
            return 1 if a * a > 2 * b * b else (-1 if a * a < 2 * b * b else 0)
        return 1 if 2 * b * b > a * a else (-1 if 2 * b * b < a * a else 0)

    def __lt__(self, o): return (self - o).sign() < 0
    def __le__(self, o): return (self - o).sign() <= 0
    def __gt__(self, o): return (self - o).sign() > 0
    def __ge__(self, o): return (self - o).sign() >= 0
    def eq(self, o): return (self - o).sign() == 0

    def approx(self):
        return float(self.a) + float(self.b) * 2 ** 0.5


R = Q2(0, F(1, 2))          # r = 1/sqrt 2
P2 = lambda e: F(2) ** e    # exact power of two (e may be negative)


class Run:
    def __init__(self, quiet=False):
        self.n = 0; self.fails = []; self.quiet = quiet

    def check(self, label, ok, detail=""):
        self.n += 1
        if not ok:
            self.fails.append(label)
        if not self.quiet:
            print(f"  check {label}: {'PASS' if ok else 'FAIL'}" + (f" ({detail})" if detail else ""))
        return ok


def lg2(x):
    """floor(log2 x) for a positive rational, exactly."""
    x = F(x); n, d = x.numerator, x.denominator
    e = n.bit_length() - d.bit_length()
    if F(2) ** e > x:
        e -= 1
    if F(2) ** (e + 1) <= x:
        e += 1
    return e


def poly_mul(p, q):
    out = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    return out


def poly_pow(p, k):
    out = [1]
    for _ in range(k):
        out = poly_mul(out, p)
    return out


def chain(c, run):
    """The parameter chain of sections 2-7, from the parameter file c.  Returns True if every check passes."""
    dE, eE, lE = c["delta_exponent"], c["eta_exponent"], c["lambda_exponent"]
    zE, epsE, sE = c["uniform_error_exponent"], c["parallel_deficit_exponent"], c["early_fraction_exponent"]
    bE, T0, K = c["radial_error_exponent"], c["t0_lower_exponent"], c["scale_ratio_lower_exponent"]
    PE, fin = c["contraction_product_upper_exponent"], F(c["final_error_upper"])
    d, eta, lam, Z, eps, sig = P2(-dE), P2(-eE), P2(-lE), P2(-zE), P2(-epsE), P2(-sE)
    q = F(1, 2) + d
    ok = run.check
    # (5)
    ok("(5) 7/10 < r < 3/4", Q2(F(7, 10)) < R < Q2(F(3, 4)))
    ok("(5) q < 9/16 < r and r - q > 1/8", Q2(q) < Q2(F(9, 16)) < R and R - q > F(1, 8))
    # (6), (9)
    e = d / R
    ok("(6) q/r - r = delta/r exactly (r^2 = 1/2), and e = delta/r <= 2 delta", (q / R - R).eq(e) and e <= 2 * d)
    ok("(9) 1 - r > 1/4, so A >= (1 - r)/e > 1/(8 delta)", 1 - R > F(1, 4) and (1 - R) / e > Q2(1 / (8 * d)),
       f"1/(8 delta) = 2^{lg2(1 / (8 * d))}")
    # (10), (11)
    ok("(10) 1 + q <= 2 and r - q > 1/8 (|b_t| <= 2A/t, >= A/(8t))", 1 + q <= 2 and R - q > F(1, 8))
    prod = F(1); good = True
    for t in range(1, 401):
        prod *= 1 + F(2, t)
        good &= prod == F((t + 1) * (t + 2), 2)
    ok("(11) prod_{j<=t} (1 + 2/j) = (t+1)(t+2)/2, t = 1..400 exactly (and the step (t+3)/(t+1) is the same identity)",
       good)
    # (13)
    lhs = Q2(q * q) - R * R * (1 - (R + eta) * (R + eta))
    rad = Q2(d + d * d) + R * eta + Q2(eta * eta / 2)
    ok("(13) q^2 - r^2 (1 - (r + eta)^2) = delta + delta^2 + r eta + eta^2/2 exactly", lhs.eq(rad))
    ok("(13) r + eta < 1 (so e A (r + eta) <= e A)", R + eta < 1)
    # (14)-(16)
    ok("(14) (5/2)*8 = 20 for the bad indices", F(5, 2) * 8 == 20)
    ecoef = (3 / R) * (R + F(5, 2) + 20 * R / eta)
    ok("(16) (3/r)(r + 5/2 + 20 r/eta) = 3 + 15/(2r) + 60/eta < 20 + 60/eta",
       ecoef.eq(3 + Q2(F(15, 2)) / R + 60 / eta) and ecoef < Q2(20 + 60 / eta))
    fcoef = (3 / R) * (F(5, 2) * R + 4)
    ok("(16) (3/r)((5/2) r + 4) = 15/2 + 12/r < 25, and 2(1 + r) <= 4", fcoef < 25 and 2 * (1 + R) <= 4,
       f"{fcoef.approx():.6f}")
    # (17)
    a17 = lg2(lam / (8 * d)); exact17 = P2(a17) == lam / (8 * d)
    ok("(17) A_{t0} >= lambda A > lambda/(8 delta) = 2^%d, and t0 >= 2 (A_0 = 1, A_1 <= 3)" % a17,
       exact17 and P2(a17) > 3)
    ok("(17) 4t^2 - (t+2)^2 = (3t+2)(t-2) >= 0 for t >= 2 (expanded)",
       poly_mul([2, 3], [-2, 1]) == [-4, -4, 3])
    ok("(17) 4 t0^2 > 2^%d gives t0 > 2^%d (the file's t0 bound)" % (a17, T0), 2 * T0 <= a17 - 2,
       f"slack 2^{(a17 - 2 - 2 * T0) / 2:.1f}")
    # (18)
    ok("(18) A_{t0} < 3 lambda A and A/A_{t0} <= 2 (k/t0)^2 give (k/t0)^2 > 1/(6 lambda) >= 2^(2*%d)" % K,
       6 * P2(2 * K) <= 1 / lam)
    ok("(18) 2k^2 (t0+1)(t0+2) - k(k+1) t0^2 = k(k t0^2 + 6 k t0 + 4k - t0^2) (expanded at 50 random points)",
       all(2 * k * k * (t + 1) * (t + 2) - k * (k + 1) * t * t == k * (k * t * t + 6 * k * t + 4 * k - t * t)
           for k, t in ((random.randrange(1, 10 ** 9), random.randrange(1, 10 ** 9)) for _ in range(50))))
    # (19)
    bb = P2(-(bE + 1))
    ok("(19) delta + delta^2 + r eta + eta^2/2 <= (2^-%d)^2 and r eta <= 2^-%d, so b <= 2^-%d" % (bE + 1, bE + 1, bE),
       rad <= Q2(bb * bb) and R * eta <= Q2(bb))
    ok("(19) 20 + 60/eta <= 2^%d" % (eE + 6), 20 + 60 / eta <= P2(eE + 6))
    tot19 = P2(-bE) + P2(eE + 6) * 2 * d / lam + 25 / P2(T0)
    ok("(19) b + 2^%d (2 delta)/lambda + 25/2^%d < Z = 2^-%d" % (eE + 6, T0, zE), tot19 < Z,
       f"total < 2^{lg2(tot19) + 1}")
    # (20), (21), (22)
    ok("(20) Z <= 1/10, so 7/5 <= R_t <= 8/5", Z <= F(1, 10))
    c2 = (1 + R) / 2                      # cos^2(pi/8)
    ok("(21) cos^2(pi/8) = (1 + r)/2 > (9/10)^2, 2(9/10) > q, and q < r = 2 cos(pi/8) sin(pi/8)",
       c2 > F(81, 100) and F(9, 5) > q and Q2(q) < R)
    ok("(22) 1/((7/5)(3/2)) = 10/21 < 1/2, 1/(4 (9/10)) = 5/18 < 1/3, 1/(16 (9/10)(7/5)(3/2)) < 1/20",
       1 / (F(7, 5) * F(3, 2)) < F(1, 2) and 1 / (4 * F(9, 10)) < F(1, 3)
       and 1 / (16 * F(9, 10) * F(7, 5) * F(3, 2)) < F(1, 20))
    ok("(22) (delta + delta^2)/(4 (9/10)(7/5)) <= delta and Z(1/2 + 1/3 + 1/20) + delta <= Z + delta",
       (d + d * d) / (4 * F(9, 10) * F(7, 5)) <= d and F(1, 2) + F(1, 3) + F(1, 20) <= 1)
    ok("(22) w_t >= 1/t when R_t >= 7/5 and t >= 4 ((7/5)(t-1) >= t iff t >= 7/2); t0 > 2^%d >= 4" % T0,
       F(7, 5) * 3 >= 4 and T0 >= 2)
    rs = [(F(random.randrange(1, 10 ** 6), 10 ** 6), F(random.randrange(1, 10 ** 6), 10 ** 6)) for _ in range(20)]
    a = 6 * (Z + d)
    ok("(23) a(1 - w/3) + 2w(Z + delta) = a for a = 6(Z + delta) (20 random w); F(h) - h = (1 - h^2)/(3h)",
       all(a * (1 - w / 3) + 2 * w * (Z + d) == a for w, _ in rs)
       and all((2 * h * h + 1) / (3 * h) - h == (1 - h * h) / (3 * h) for _, h in rs))
    p = poly_mul(poly_pow([-1, 3], 3), [1, 1])
    diff = [x - y for x, y in zip([0, 0, 0, 0, 27], p + [0] * (5 - len(p)))]
    ok("(23) 27 t^4 - (3t - 1)^3 (t + 1) = 18 t^2 - 8 t + 1, positive (discriminant 64 - 72 < 0)",
       diff == [1, -8, 18, 0, 0] and 8 * 8 - 4 * 18 * 1 < 0)
    ok("(23) the product < (t0/k)^(1/3) < 2^(-%d/3) <= 2^-%d" % (K, PE), 3 * PE <= K)
    ok("(23) 6(Z + delta) + 2^-%d < eps = 2^-%d" % (PE, epsE), 6 * (Z + d) + P2(-PE) < eps,
       f"left < 2^{lg2(6 * (Z + d) + P2(-PE)) + 1}")
    # (24)-(28)
    ok("(24) Z + delta <= eps/2 and eps <= 1/100 (so cos beta_k >= F(h) - 2(Z + delta) >= 1 - 2 eps)",
       Z + d <= eps / 2 and eps <= F(1, 100))
    ok("(25) (Z + (3/2) eps)/(1 - eps) <= 3 eps, (8/5)/(9/10) * 2 <= 4, eps <= sqrt(eps)",
       (Z + F(3, 2) * eps) / (1 - eps) <= 3 * eps and F(8, 5) / F(9, 10) * 2 <= 4 and eps <= 1)
    if epsE % 4 == 0:
        s = P2(-epsE // 2)                # sqrt(eps), exact
        ok("(26) sqrt(2 eps)/(1 - eps) + eps/(1 - eps) <= 3 sqrt(eps), i.e. sqrt 2 + s <= 3(1 - s^2), s = 2^-%d"
           % (epsE // 2), Q2(s, 0) + Q2(0, 1) <= Q2(3 * (1 - s * s)))
        ok("(27) 7 + 3q <= 16 (|S| <= q < 9/16)", 7 + 3 * q <= 16)
        ok("(28) delta, delta^2 <= sqrt(eps); 18 <= 25, so d <= 5 eps^(1/4) = 5 * 2^-%d" % (epsE // 4),
           d <= s and d * d <= s and 18 <= 25)
        e4 = P2(-epsE // 4)
    else:
        ok("(26)-(28) eps is a power 2^-4j (as the note's eps^(1/4) = 2^-60 needs)", False)
        e4 = None
    # (29)-(31)
    ok("(29) k > t0 * 2^%d > 2^%d and sigma k > 2 (so m = floor(sigma k) >= sigma k/2 >= 1, m < k)" % (K, T0 + K),
       P2(T0 + K - sE) >= 2)
    bern = poly_pow([1, 1], 8)            # (1 + x)^8: every coefficient >= 0 and the linear one is 8
    ok("(29) (1 + x)^8 >= 1 + 8x for x >= 0 (expanded: coefficients %s)" % bern,
       bern[0] == 1 and bern[1] == 8 and all(x >= 0 for x in bern))
    ok("(29) sigma^(1/8) = 2^-%d exactly" % (sE // 8), sE % 8 == 0)
    ok("(30) 2/r = 2 sqrt 2 < 3, and k/j <= k/m <= 2/sigma", 2 / R < 3)
    if e4 is not None and sE % 8 == 0:
        E = 2 * P2(-(sE // 8)) + (6 / sig) * 5 * e4
        ok("(31) q + 1/2 < 2 and 2 * 2^-%d + (6/sigma) * 5 * 2^-%d < %s" % (sE // 8, epsE // 4, fin),
           q + F(1, 2) < 2 and E < fin, f"E/A < {float(E):.6f}")
    else:
        ok("(31) the early-fraction and eps exponents are divisible by 8 and 4", False)
    ok("(32) r (r - 1/2) > %s (the contradiction; exactly 1/2 - r/2 = 0.1464...), and > 7/50" % fin,
       R * (R - F(1, 2)) > fin and R * (R - F(1, 2)) > F(7, 50))
    return not run.fails


class G:
    """Gaussian rationals."""
    __slots__ = ("x", "y")

    def __init__(self, x, y=0):
        self.x, self.y = F(x), F(y)

    def __add__(self, o): return G(self.x + o.x, self.y + o.y)
    def __sub__(self, o): return G(self.x - o.x, self.y - o.y)
    def __mul__(self, o): return G(self.x * o.x - self.y * o.y, self.x * o.y + self.y * o.x)
    def scale(self, k): return G(self.x * k, self.y * k)
    def __eq__(self, o): return self.x == o.x and self.y == o.y


def identities(run, rng, trials=40, perturb=False):
    """(1) and (3) on random tuples z_1 = 1, z_2..z_n of Gaussian rationals; (15) on random nonnegative a_j."""
    ok1 = ok3 = True
    for _ in range(trials):
        n = rng.randrange(2, 10)
        z = [G(1)] + [G(F(rng.randrange(-9, 10), rng.randrange(1, 7)), F(rng.randrange(-9, 10), rng.randrange(1, 7)))
                      for _ in range(n - 1)]
        T = n + 4
        S = [None] + [None] * T
        pw = [G(1)] * n
        for l in range(1, T + 1):
            pw = [pw[i] * z[i] for i in range(n)]
            acc = G(0)
            for v in pw:
                acc = acc + v
            S[l] = acc
        b = [G(1)]                                   # B(x) = prod_{i>=2} (1 - z_i x)
        for i in range(1, n):
            nb = b + [G(0)]
            for j in range(len(b)):
                nb[j + 1] = nb[j + 1] - b[j] * z[i]
            b = nb
        b = b + [G(0)] * (T + 1 - len(b))
        if perturb:
            b[1] = b[1] + G(F(1, 1000))
        C = []; acc = G(0)
        for j in range(T + 1):
            acc = acc + b[j]; C.append(acc)
        for t in range(1, T + 1):                    # (1): C_{t-1} - t b_t = sum_{j<t} b_j S_{t-j}
            rhs = G(0)
            for j in range(t):
                rhs = rhs + b[j] * S[t - j]
            ok1 &= (C[t - 1] - b[t].scale(t)) == rhs
        for t in range(1, T + 1):                    # (3): 2 U_t - t C_t = sum_{j<t} C_j S_{t-j}
            U = G(0)
            for j in range(t):
                U = U + C[j]
            rhs = G(0)
            for j in range(t):
                rhs = rhs + C[j] * S[t - j]
            ok3 &= (U.scale(2) - C[t].scale(t)) == rhs
    ok15 = True
    for _ in range(trials):
        a = [F(rng.randrange(0, 50), rng.randrange(1, 9)) for _ in range(rng.randrange(2, 12))]
        for t in range(1, len(a)):
            A = [sum(a[:j + 1]) for j in range(len(a))]
            W = sum(A[:t]); M = sum(j * a[j] for j in range(1, t + 1))
            ok15 &= t * A[t] == W + M
    return ok1, ok3, ok15


def main():
    args = [x for x in sys.argv[1:] if not x.startswith("--")]
    controls = "--controls" in sys.argv
    if len(args) != 1:
        print(__doc__); sys.exit(2)
    t_start = time.time()
    root = args[0]
    cpath = os.path.join(root, "biro", "lower_certificate.json")
    ppath = os.path.join(root, "biro", "LOWER_PROOF.md")
    raw = open(cpath, "rb").read()
    print("42a [Pin26c] C_42 >= 1/2 + 2^-4096 (PR #310), from LOWER_PROOF.md and lower_certificate.json")
    run = Run()
    run.check("parameter file sha256 = the manifest's f6681b2f...c79b", hashlib.sha256(raw).hexdigest() == CERT_SHA)
    run.check("proof sha256 = the manifest's 3acab5f8...ba1a",
              hashlib.sha256(open(ppath, "rb").read()).hexdigest() == PROOF_SHA)
    c = json.loads(raw)
    print("  parameters: delta = 2^-%(delta_exponent)d, eta = 2^-%(eta_exponent)d, lambda = 2^-%(lambda_exponent)d, "
          "Z = 2^-%(uniform_error_exponent)d, eps = 2^-%(parallel_deficit_exponent)d, "
          "sigma = 2^-%(early_fraction_exponent)d" % c)
    random.seed(20261108)
    chain(c, run)
    rng = random.Random(310)
    ok1, ok3, ok15 = identities(run, rng)
    run.check("(1) C_{t-1} - t b_t = sum_{j<t} b_j S_{t-j}, 40 random tuples, n = 2..9, t = 1..n+4, exactly", ok1)
    run.check("(3) 2 U_t - t C_t = sum_{j<t} C_j S_{t-j}, same tuples, exactly", ok3)
    run.check("(15) t A_t = W_t + M_t for 40 random nonnegative sequences", ok15)
    main_ok = not run.fails
    if controls:
        print("  -- red controls (each altered input must be refused)")
        alt = [("the claim strengthened to 1/2 + 2^-4094, auxiliaries unchanged (t0 > 2^1021 no longer follows)",
                {"delta_exponent": 4094}),
               ("eta = 2^-200 (b <= 2^-254 fails)", {"eta_exponent": 200}),
               ("lambda = 2^-4000 (t0 > 2^1021 fails)", {"lambda_exponent": 4000}),
               ("Z = 2^-256 (b + ... < Z fails)", {"uniform_error_exponent": 256}),
               ("eps = 2^-252 (6(Z + delta) + 2^-340 < eps fails)", {"parallel_deficit_exponent": 252}),
               ("sigma = 2^-8 (the early part of (31) is too large)", {"early_fraction_exponent": 8}),
               ("final error 1/32 instead of 1/16 ((31) fails)", {"final_error_upper": "1/32"}),
               ("product bound 2^-341 (3 * 341 > 1022)", {"contraction_product_upper_exponent": 341})]
        for label, change in alt:
            cc = dict(c); cc.update(change)
            sub = Run(quiet=True)
            chain(cc, sub)
            run.check("control: " + label + " is refused", bool(sub.fails), "first failing step: " + sub.fails[0]
                      if sub.fails else "nothing failed")
        b1, b3, _ = identities(run, random.Random(311), trials=5, perturb=True)
        run.check("control: b_1 perturbed by 1/1000 breaks identities (1) and (3)", not b1 and not b3)
    print(f"  checks: {run.n}, failed: {len(run.fails)}")
    print(f"  run time: {time.time() - t_start:.1f} s")
    print("RESULT: ALL CHECKS PASS" if not run.fails else "RESULT: FAIL " + "; ".join(run.fails))
    sys.exit(0 if not run.fails else 1)


if __name__ == "__main__":
    main()
