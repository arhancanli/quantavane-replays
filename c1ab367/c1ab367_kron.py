#!/usr/bin/env python3
"""Program A for the step-function certificates of PR #367 (1a) and PR #368 (1b): exact integers, Kronecker substitution.

  c1ab367_kron.py 1a FILE [--claim P/Q] [--prefix DIGITS] [--printed DEC] [--argmax K] [--bump I D]
  c1ab367_kron.py 1b FILE [--claim P/Q] [--printed DEC] [--argmax M] [--bump I D] [--move I J D]

1a (constants/1a.md): f >= 0 on [-1/4, 1/4] with n equal steps of heights a_i. f*f is linear between the knots, where it
equals h c_k with h = 1/(2n), c_k = sum_{i+j=k} a_i a_j, and int f = h S, S = sum a_i, so
    max_t (f*f)(t) / (int f)^2 = 2n max_k c_k / S^2,
an upper bound for C_1a. FILE: one height per line ('#' lines skipped); each height an integer or an exact decimal.
1b (constants/1b.md): h on [0, 2] with n equal steps, heights k_i / Q in [0, 1], mass exactly 1 (sum k_i = nQ/2). The
overlap x -> int h(t)(1 - h(t + x)) dt (the second factor zero off [0, 2]) is linear between multiples of 2/n, so
    sup_x = (2/n) max_m sum_{0 <= i, i+m < n} h_i (1 - h_{i+m}) = 2 max_m [Q P_m - R_m] / (n Q^2),
with P_m the sum of k_i over the valid i and R_m = sum_i k_i k_{i+m}; an upper bound for C_1b. FILE: JSON with integer
keys n, Q and the list heights_int (read as data).
The whole sequence (c_k, or R_m) comes from one big-integer product X^2 or X*Y with X = sum a_i 2^(b i), b larger than
any coefficient needs, so the base-2^b digits of the product are the coefficients exactly (no carries). Every value
at the maximum, and 40 others at random, is summed again term by term. A fingerprint sum_k v_k 3^k mod 2^61 - 1 of the
whole sequence is printed for comparison with program B. Program B reduces the heights to the same primitive integers.
--bump I D adds the integer D to height I (1a: in units of the primitive integer heights; 1b: in units of 1/Q), and
--move I J D moves D/Q from height J to height I (red controls).
Standard library only; run with python3 -I.
"""
import json
import math
import random
import sys
import time
from fractions import Fraction

PRIME = (1 << 61) - 1


def fingerprint(vals):
    f, p = 0, 1
    for v in vals:
        f = (f + (v % PRIME) * p) % PRIME
        p = (p * 3) % PRIME
    return f


def digits(fr, nd=40):
    """fr >= 0 truncated to nd decimals, as a string."""
    q = fr.numerator * 10 ** nd // fr.denominator
    s = str(q).rjust(nd + 1, "0")
    return s[:-nd] + "." + s[-nd:]


def read_heights_txt(path):
    vals = []
    with open(path) as fh:
        for line in fh:
            t = line.strip()
            if not t or t.startswith("#"):
                continue
            vals.append(Fraction(t))         # exact: an integer or a finite decimal string
    den = 1
    for v in vals:
        den = den * v.denominator // math.gcd(den, v.denominator)
    ints = [int(v * den) for v in vals]
    g = 0
    for v in ints:
        g = math.gcd(g, v)
    g = g or 1
    return [v // g for v in ints], Fraction(den, g)     # primitive integer heights (the ratio is scale-free)


def pack(vals, b):
    w = b // 8
    return int.from_bytes(b"".join(v.to_bytes(w, "little") for v in vals), "little")


def unpack(x, b, count):
    w = b // 8
    raw = x.to_bytes(w * count, "little")
    return [int.from_bytes(raw[i * w:(i + 1) * w], "little") for i in range(count)]


def main():
    args = sys.argv[1:]
    kind, path = args[0], args[1]
    opts = {}
    i = 2
    while i < len(args):
        if args[i] == "--bump":
            opts["bump"] = (int(args[i + 1]), Fraction(args[i + 2]))
            i += 3
        elif args[i] == "--move":
            opts["move"] = (int(args[i + 1]), int(args[i + 2]), int(args[i + 3]))
            i += 4
        else:
            opts[args[i][2:]] = args[i + 1]
            i += 2
    t0 = time.time()
    ok = True

    def check(desc, cond, detail=""):
        nonlocal ok
        ok &= bool(cond)
        print("A: check %s: %s%s" % (desc, "PASS" if cond else "FAIL", (" (%s)" % detail) if detail else ""))
        return cond

    if kind == "1a":
        a, den = read_heights_txt(path)
        if "bump" in opts:
            bi, bd = opts["bump"]
            a[bi] += int(bd)
            print("A: CONTROL height %d changed by %s (units of the primitive integer heights)" % (bi, bd))
        n = len(a)
        print("A: 1a, n = %d heights, scaled by %s to primitive integers of at most %d bits" % (n, den, max(a).bit_length()))
        if not check("every height >= 0", all(v >= 0 for v in a)):
            print("A: VERDICT REFUSED")
            return 1
        S = sum(a)
        check("S = sum a_i > 0", S > 0)
        b = 2 * max(a).bit_length() + n.bit_length() + 2
        b += (-b) % 8
        X = pack(a, b)
        c = unpack(X * X, b, 2 * n - 1)
        cmax = max(c)
        arg = [k for k, v in enumerate(c) if v == cmax]
        rng = random.Random(367)
        ks = arg + [rng.randrange(2 * n - 1) for _ in range(40)]
        direct = all(sum(a[j] * a[k - j] for j in range(max(0, k - n + 1), min(k, n - 1) + 1)) == c[k] for k in ks)
        check("c_k summed again term by term at every maximum and at 40 random k", direct, "%d values" % len(ks))
        value = Fraction(2 * n * cmax, S * S)
        print("A: S = %d" % S)
        print("A: max_k c_k = %d at k = %s" % (cmax, ",".join(map(str, arg))))
    elif kind == "1b":
        with open(path) as fh:
            d = json.load(fh)
        n, Q, k = d["n"], d["Q"], list(d["heights_int"])
        if not (isinstance(n, int) and isinstance(Q, int) and all(isinstance(v, int) for v in k)):
            check("integer data", False)
            return 1
        if "bump" in opts:
            bi, bd = opts["bump"]
            k[bi] += int(bd)
            print("A: CONTROL height %d changed by %s (units of 1/Q)" % (bi, bd))
        if "move" in opts:
            mi, mj, md = opts["move"]
            k[mi] += md
            k[mj] -= md
            print("A: CONTROL %d/Q moved from height %d to height %d (mass unchanged)" % (md, mj, mi))
        print("A: 1b, n = %d heights k_i / Q, Q = %d" % (n, Q))
        check("len(heights_int) = n", len(k) == n)
        check("0 <= k_i <= Q for every i", all(0 <= v <= Q for v in k))
        check("mass exactly 1: sum k_i = n Q / 2", 2 * sum(k) == n * Q, "sum = %d, nQ/2 = %s" % (sum(k), Fraction(n * Q, 2)))
        b = 2 * Q.bit_length() + n.bit_length() + 2
        b += (-b) % 8
        prod = unpack(pack(k, b) * pack(k[::-1], b), b, 2 * n - 1)
        pre = [0]
        for v in k:
            pre.append(pre[-1] + v)
        vals = []
        for m in range(-(n - 1), n):
            lo, hi = max(0, -m), min(n, n - m)
            vals.append(Q * (pre[hi] - pre[lo]) - prod[n - 1 - m])
        vmax = max(vals)
        arg = [m - (n - 1) for m, v in enumerate(vals) if v == vmax]
        rng = random.Random(368)
        ms = arg + [rng.randrange(-(n - 1), n) for _ in range(40)]
        direct = all(sum(k[i] * (Q - k[i + m]) for i in range(max(0, -m), min(n, n - m))) == vals[m + n - 1] for m in ms)
        check("Q P_m - R_m summed again term by term at every maximum and at 40 random m", direct, "%d values" % len(ms))
        c = vals
        value = Fraction(2 * vmax, n * Q * Q)
        print("A: max_m [Q P_m - R_m] = %d at m = %s" % (vmax, ",".join(map(str, arg))))
    else:
        print(__doc__)
        return 2
    print("A: fingerprint of the whole sequence (mod 2^61-1): %d" % fingerprint(c))
    print("A: exact value %d/%d" % (value.numerator, value.denominator))
    print("A: = %s..." % digits(value))
    if "claim" in opts:
        check("the exact value equals the stated fraction %s" % opts["claim"], value == Fraction(opts["claim"]))
    if "prefix" in opts:
        nd = len(opts["prefix"].split(".")[1])
        check("the exact value truncated to %d decimals is the stated %s" % (nd, opts["prefix"]),
              digits(value, nd) == opts["prefix"])
    if "argmax" in opts:
        check("the maximum is attained at the stated index %s" % opts["argmax"], int(opts["argmax"]) in arg)
    if "printed" in opts:
        pr = Fraction(opts["printed"])
        check("the printed upper bound %s is at least the exact value" % opts["printed"], pr >= value,
              "printed - exact = %.3e" % float(pr - value))
    print("A: run time %.1f s" % (time.time() - t0))
    print("A: VERDICT %s" % ("PASS" if ok else "REFUSED"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
