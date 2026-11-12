#!/usr/bin/env python3
"""Exact replay of two one-dimensional lower-bound rows [Pat2026] (teorth/optimizationproblems, main 87b4ebb):

  20a  C_20a >= 8, where C_20a = sup over n and isotropic log-concave X in R^n of Var(|X|^2)/n: n = 1, X = E - 1 with
       E standard exponential. Exactly: the moments of E are int_0^oo t^k e^-t dt = k!, so E X = 0, E X^2 = 1
       (isotropic), E X^4 = sum_j C(4,j) (-1)^(4-j) j! = 9, Var(X^2) = 9 - 1 = 8; the density e^-(x+1) on [-1, oo) is
       log-concave (its log is affine on a half-line).
  20c  C_20c >= sqrt 3, where psi_mu = 1/h_mu and h_mu = inf_A mu^+(A)/min(mu(A), 1 - mu(A)): n = 1, mu uniform on
       [-sqrt 3, sqrt 3] (density 1/(2 sqrt 3), variance 3/3 = 1, mean 0, log-concave), A = (-oo, 0]: mu(A) = 1/2 and
       mu^+(A) = 1/(2 sqrt 3) (the density at the boundary point 0), so h_mu <= 1/sqrt 3 and psi_mu >= sqrt 3. Exact in
       Q(sqrt 3) (pairs a + b sqrt 3).
The rows' remarks on products of independent copies (higher n) are read, not checked; not needed for the bounds.
Red controls (--controls): X = E (not centred) is refused (E X = 1); the uniform law on [-1, 1] is refused (variance
1/3, not isotropic); claims C_20a >= 9 and C_20c >= 1.7321 (above sqrt 3) from these examples are refused.
usage: python3 -I pat2026_check.py [--controls]      (Python 3.9+, standard library only)
"""
import sys
from fractions import Fraction as F
from math import comb, factorial

CONTROLS = "--controls" in sys.argv


def moment_shifted_exp(k, shift):
    """E (E - shift)^k for E ~ Exp(1), exactly (shift rational)."""
    return sum(comb(k, j) * F(-shift) ** (k - j) * factorial(j) for j in range(k + 1))


def example_20a(shift=1):
    m1, m2, m4 = moment_shifted_exp(1, shift), moment_shifted_exp(2, shift), moment_shifted_exp(4, shift)
    return m1, m2, m4, m4 - m2 * m2


# numbers a + b sqrt3 with rational a, b
def mul(x, y):
    return (x[0] * y[0] + 3 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def inv(x):
    n = x[0] * x[0] - 3 * x[1] * x[1]
    return (x[0] / n, -x[1] / n)


def gt(x, y):
    """x > y for a + b sqrt3 numbers, exactly."""
    a, b = x[0] - y[0], x[1] - y[1]
    if a >= 0 and b >= 0:
        return a > 0 or b > 0
    if a <= 0 and b <= 0:
        return False
    return (a * a > 3 * b * b) if a > 0 else (3 * b * b > a * a)


def example_20c(half_width):
    """Uniform on [-w, w] (w a + b sqrt3 number): variance w^2/3, density 1/(2w); A = (-oo, 0]."""
    var = mul(half_width, half_width)
    var = (var[0] / 3, var[1] / 3)
    dens = inv((2 * half_width[0], 2 * half_width[1]))
    muA = F(1, 2)
    h = (dens[0] / muA, dens[1] / muA)          # mu^+(A)/min(mu(A), 1 - mu(A))
    psi = inv(h)
    return var, dens, psi


def main():
    ok = True
    m1, m2, m4, v = example_20a()
    print("20a row [Pat2026]: X = E - 1, E ~ Exp(1): E X = %s, E X^2 = %s, E X^4 = %s, Var(X^2) = %s" % (m1, m2, m4, v))
    ok1 = m1 == 0 and m2 == 1 and v == 8
    print("RESULT [Pat2026] 20a: n = 1, isotropic log-concave X with Var(|X|^2)/n = 8, so C_20a >= 8: %s"
          % ("HOLDS" if ok1 else "FAILS"))
    ok &= ok1
    var, dens, psi = example_20c((F(0), F(1)))
    print("20c row [Pat2026]: uniform on [-sqrt3, sqrt3]: variance %s + %s sqrt3, density %s + %s sqrt3; "
          "A = (-oo, 0]: mu(A) = 1/2, mu^+(A) = density; psi_mu >= %s + %s sqrt3"
          % (var[0], var[1], dens[0], dens[1], psi[0], psi[1]))
    ok2 = var == (F(1), F(0)) and psi == (F(0), F(1)) and dens == (F(0), F(1, 6))
    print("RESULT [Pat2026] 20c: psi_mu >= sqrt 3 = 1.7320508..., so C_20c >= sqrt 3: %s" % ("HOLDS" if ok2 else "FAILS"))
    ok &= ok2
    status = 0 if ok else 1
    if CONTROLS:
        print("red controls (each must fail):")
        c_ok = True
        m1, m2, m4, v = example_20a(shift=0)
        print("  X = E (no shift): E X = %s: %s" % (m1, "accepted (UNEXPECTED)" if m1 == 0 else "refused, as expected"))
        c_ok &= m1 != 0
        var, dens, psi = example_20c((F(1), F(0)))
        print("  uniform on [-1, 1]: variance %s: %s" % (var[0], "isotropic (UNEXPECTED)" if var == (1, 0) else
                                                         "not isotropic, refused as expected"))
        c_ok &= var != (F(1), F(0))
        m1, m2, m4, v = example_20a()
        print("  a claimed C_20a >= 9 from E - 1: %s" % ("holds (UNEXPECTED)" if v >= 9 else "refused, as expected"))
        c_ok &= not v >= 9
        _, _, psi = example_20c((F(0), F(1)))
        claim = not gt((F(17321, 10000), F(0)), psi)     # 1.7321 <= sqrt3 ?
        print("  a claimed C_20c >= 1.7321 from the uniform law: %s" % ("holds (UNEXPECTED)" if claim else
                                                                       "refused, as expected (sqrt 3 < 1.7321)"))
        c_ok &= not claim
        print("controls all as expected" if c_ok else "CONTROLS FAILED")
        if not c_ok:
            status = 1
    return status


if __name__ == "__main__":
    sys.exit(main())
