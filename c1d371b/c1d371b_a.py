#!/usr/bin/env python3
"""Program A for the 1d step function of PR #371 (its eighth certificate, 1.4546676818): exact integers, one Kronecker squaring with balanced digits.

  c1d371b_a.py FILE [--claim P/Q] [--argmax M] [--printed DEC] [--bump I D] [--negate I]

1d (constants/1d.md): f real-valued on [-1/4, 1/4], n equal cells of width h = 1/(2n) and signed heights v_0..v_{n-1}.
f*f is linear between the knots t = (m + 1) h - 1/2, where it equals h c_m with c_m = sum_{i+j=m} v_i v_j, so the
maximum of |f*f| is attained at a knot, and with int f = h S, S = sum v_i,
    max_t |f*f|(t) / (int f)^2 = 2n max_m |c_m| / S^2,
an upper bound for C_1d. FILE: one integer per line ('#' lines and blank lines skipped).

The sequence c_m comes from one squaring X^2, X = sum v_i 2^(b i) (a signed integer), with b chosen so that every
|c_m| < 2^(b-1): the coefficients are the balanced base-2^b digits of X^2 (digits in [-2^(b-1), 2^(b-1)), carried
upwards). The values at every maximizing m, and 40 other m at random, are summed again term by term. A fingerprint
sum_m c_m 3^m mod 2^61 - 1 of the whole sequence is printed for comparison with program B (which shares no code).
--bump I D adds D to height I; --negate I flips the sign of height I (red controls).
Standard library only; run with python3 -I.
"""
import random
import sys
import time
from fractions import Fraction

P61 = (1 << 61) - 1


def read_ints(path):
    out = []
    with open(path, "rb") as fh:
        for raw in fh:
            t = raw.strip()
            if not t or t.startswith(b"#"):
                continue
            out.append(int(t))
    return out


def args_of(argv):
    opts = {"claim": None, "argmax": None, "printed": None, "bump": [], "negate": []}
    i = 2
    while i < len(argv):
        a = argv[i]
        if a == "--claim":
            opts["claim"] = Fraction(argv[i + 1]); i += 2
        elif a == "--argmax":
            opts["argmax"] = int(argv[i + 1]); i += 2
        elif a == "--printed":
            opts["printed"] = argv[i + 1]; i += 2
        elif a == "--bump":
            opts["bump"].append((int(argv[i + 1]), int(argv[i + 2]))); i += 3
        elif a == "--negate":
            opts["negate"].append(int(argv[i + 1])); i += 2
        else:
            raise SystemExit("unknown option " + a)
    return opts


def balanced_digits(x, b, count):
    """The first `count` balanced base-2^b digits of x >= 0 (b a multiple of 8)."""
    nb = b // 8
    raw = x.to_bytes((x.bit_length() + 7) // 8 + nb, "little")
    half, full = 1 << (b - 1), 1 << b
    out, carry = [], 0
    for k in range(count):
        d = int.from_bytes(raw[k * nb:(k + 1) * nb], "little") + carry
        if d >= half:
            d -= full
            carry = 1
        else:
            carry = 0
        out.append(d)
    rest = int.from_bytes(raw[count * nb:], "little") + carry
    return out, rest


def main():
    t0 = time.time()
    path = sys.argv[1]
    opts = args_of(sys.argv)
    v = read_ints(path)
    n = len(v)
    for i, d in opts["bump"]:
        v[i] += d
        print("CONTROL height %d changed by %d" % (i, d))
    for i in opts["negate"]:
        v[i] = -v[i]
        print("CONTROL height %d negated" % i)
    S = sum(v)
    vmax = max(abs(x) for x in v)
    print("cells n = %d, negative heights %d, largest |v_i| = %d, S = %d" % (n, sum(1 for x in v if x < 0), vmax, S))
    if S == 0:
        print("VERDICT REFUSED (S = 0)"); return
    bound = n * vmax * vmax                      # |c_m| <= n max|v|^2
    b = bound.bit_length() + 2
    b += (-b) % 8
    nb = b // 8                                  # X = P - N, P and N packed from their base-2^b digits as bytes
    P = int.from_bytes(b"".join((x if x > 0 else 0).to_bytes(nb, "little") for x in v), "little")
    N = int.from_bytes(b"".join((-x if x < 0 else 0).to_bytes(nb, "little") for x in v), "little")
    X = P - N
    sq = X * X
    c, rest = balanced_digits(sq, b, 2 * n - 1)
    if rest != 0:
        print("FAIL: digits beyond 2n - 1"); print("VERDICT REFUSED"); return
    M = max(abs(x) for x in c)
    arg = [m for m, x in enumerate(c) if abs(x) == M]
    rnd = random.Random(3718)
    probe = sorted(set(arg) | set(rnd.randrange(2 * n - 1) for _ in range(40)))
    for m in probe:
        lo, hi = max(0, m - n + 1), min(m, n - 1)
        if sum(v[i] * v[m - i] for i in range(lo, hi + 1)) != c[m]:
            print("FAIL: the squaring and the direct sum differ at m = %d" % m); print("VERDICT REFUSED"); return
    fp, p3 = 0, 1
    for x in c:
        fp = (fp + (x % P61) * p3) % P61
        p3 = p3 * 3 % P61
    val = Fraction(2 * n * M, S * S)
    print("max |c_m| = %d at m = %s (sign %s); %d values summed again term by term, all equal"
          % (M, ",".join(map(str, arg)), ",".join("+" if c[m] > 0 else "-" for m in arg), len(probe)))
    print("fingerprint sum c_m 3^m mod 2^61-1 = %d" % fp)
    q = val.numerator * 10 ** 30 // val.denominator
    print("exact value %d/%d" % (val.numerator, val.denominator))
    print("decimal (truncated) %s.%s" % (str(q)[:-30], str(q)[-30:]))
    ok = True
    if opts["claim"] is not None:
        same = val == opts["claim"]
        print("stated fraction: %s" % ("equal" if same else "DIFFERENT"))
        ok &= same
    if opts["argmax"] is not None:
        same = arg == [opts["argmax"]]
        print("stated maximizing knot m = %d only: %s" % (opts["argmax"], "yes" if same else "NO"))
        ok &= same
    if opts["printed"] is not None:
        pr = Fraction(opts["printed"])
        holds = val <= pr
        gap = pr - val
        print("printed upper bound %s: %s (printed minus exact = %.3e)"
              % (opts["printed"], "holds" if holds else "FAILS", float(gap)))
        ok &= holds
    print("run time %.1f s" % (time.time() - t0))
    print("VERDICT %s" % ("PASS" if ok else "REFUSED"))


if __name__ == "__main__":
    main()
