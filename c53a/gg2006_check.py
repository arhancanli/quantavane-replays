#!/usr/bin/env python3
"""Exact replay of the 53a lower-bound row [GG2006] (teorth/optimizationproblems, main 87b4ebb): C_53 >= 3, where
C_53 = sup_n (D(C_n^3) - 1)/(n - 1) and D is the Davenport constant.

The row's construction: the sequence S = e1^(n-1) e2^(n-1) e3^(n-1) over C_n^3 (each standard generator n - 1 times) is
zero-sum free, so D(C_n^3) >= |S| + 1 = 3(n - 1) + 1. Checked exactly for n = 2..9 by enumerating every nonempty
subsequence (by multiplicities (a, b, c), 0 <= a, b, c <= n - 1, not all 0: its sum is (a, b, c) mod n, never 0),
and by a second route for n = 2..5: all 2^(3(n-1)) - 1 nonempty index subsets of S summed directly. So
(D(C_n^3) - 1)/(n - 1) >= 3 and C_53 >= 3.
Red controls (--controls): S with one more copy of e1 (e1^n e2^(n-1) e3^(n-1)) has a zero-sum subsequence (e1^n) and is
refused; at n = 3, e1^2 e2^2 e3 (e1 + e2) (one e3 replaced by e1 + e2) has the zero sum e1^2 e2^2 (e1 + e2) and is
refused.
usage: python3 -I gg2006_check.py [--controls]      (Python 3.9+, standard library only)
"""
import itertools
import sys

CONTROLS = "--controls" in sys.argv


def zero_sum_free_by_mult(seq, n):
    """seq: list of (element, multiplicity); check every nonempty sub-multiset."""
    ranges = [range(m + 1) for _, m in seq]
    for mult in itertools.product(*ranges):
        if not any(mult):
            continue
        s = [0, 0, 0]
        for (g, _), k in zip(seq, mult):
            for i in range(3):
                s[i] = (s[i] + k * g[i]) % n
        if s == [0, 0, 0]:
            return False, mult
    return True, None


def zero_sum_free_by_subsets(elems, n):
    L = len(elems)
    for mask in range(1, 1 << L):
        s = [0, 0, 0]
        for i in range(L):
            if mask >> i & 1:
                for k in range(3):
                    s[k] = (s[k] + elems[i][k]) % n
        if s == [0, 0, 0]:
            return False
    return True


E1, E2, E3 = (1, 0, 0), (0, 1, 0), (0, 0, 1)


def main():
    ok = True
    print("53a row [GG2006]: S = e1^(n-1) e2^(n-1) e3^(n-1) over C_n^3 is zero-sum free, so D(C_n^3) >= 3(n-1) + 1")
    for n in range(2, 10):
        g1, _ = zero_sum_free_by_mult([(E1, n - 1), (E2, n - 1), (E3, n - 1)], n)
        g2 = zero_sum_free_by_subsets([E1] * (n - 1) + [E2] * (n - 1) + [E3] * (n - 1), n) if n <= 5 else None
        print("  n = %d: |S| = %d; zero-sum free by multiplicities: %s; by all index subsets: %s; (D - 1)/(n - 1) >= 3"
              % (n, 3 * (n - 1), g1, "not run (n > 5)" if g2 is None else g2))
        ok &= g1 and (g2 is None or g2)
    print("RESULT [GG2006]: D(C_n^3) >= 3(n - 1) + 1, so C_53 >= 3: %s" % ("HOLDS" if ok else "FAILS"))
    status = 0 if ok else 1
    if CONTROLS:
        print("red controls (each must fail):")
        c_ok = True
        g, m = zero_sum_free_by_mult([(E1, 4), (E2, 3), (E3, 3)], 4)
        print("  n = 4, e1^4 e2^3 e3^3: %s (multiplicities %s sum to 0)" % ("zero-sum free (UNEXPECTED)" if g else
                                                                           "refused, as expected", m))
        c_ok &= not g
        g, m = zero_sum_free_by_mult([(E1, 2), (E2, 2), (E3, 1), ((1, 1, 0), 1)], 3)
        print("  n = 3, e1^2 e2^2 e3 (e1 + e2): %s (multiplicities %s)" % ("zero-sum free (UNEXPECTED)" if g else
                                                                          "refused, as expected", m))
        c_ok &= not g
        print("controls all as expected" if c_ok else "CONTROLS FAILED")
        if not c_ok:
            status = 1
    return status


if __name__ == "__main__":
    sys.exit(main())
