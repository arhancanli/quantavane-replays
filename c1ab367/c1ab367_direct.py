#!/usr/bin/env python3
"""Program B for the step-function certificates of PR #367 (1a) and PR #368 (1b): its own reader and the plain double
loops of c1ab367_direct.c (compiled into bin/ on first use with the system C compiler), which take no product of big
integers. This file only reads the heights, hands them to the C program as integers, and turns its integer output
into the exact value.

  c1ab367_direct.py 1a FILE [--claim P/Q] [--prefix DIGITS] [--printed DEC] [--argmax K] [--bump I D]
  c1ab367_direct.py 1b FILE [--claim P/Q] [--printed DEC] [--argmax M] [--bump I D] [--move I J D]

1a: value = 2n max_k c_k / S^2, c_k = sum_{i+j=k} a_i a_j, S = sum a_i (heights scaled to integers by the power of
ten that clears every decimal; the ratio does not depend on the scale). 1b: value = 2 max_m v_m / (n Q^2),
v_m = sum_i k_i (Q - k_{i+m}). The red-control options mean the same as in program A.
Standard library and a C compiler; run with python3 -I.
"""
import json
import math
import os
import subprocess
import sys
import time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))


def build():
    exe = os.path.join(HERE, "bin", "c1ab367_direct")
    src = os.path.join(HERE, "c1ab367_direct.c")
    if not os.path.exists(exe) or os.path.getmtime(exe) < os.path.getmtime(src):
        os.makedirs(os.path.dirname(exe), exist_ok=True)
        subprocess.check_call(["cc", "-O2", "-o", exe, src])
    return exe


def decimal_lines(path):
    """The height strings of a text certificate, '#' lines and blank lines skipped."""
    out = []
    with open(path) as fh:
        for line in fh:
            t = line.strip()
            if t and not t.startswith("#"):
                out.append(t)
    return out


def scale_to_integers(strs):
    """(primitive integers, scale): every string is d, d.ddd, or either with an exponent e-N / e+N (digits only).
    Each is mantissa * 10^-(f - e) with f its decimals and e its exponent; F = the largest f - e (at least 0)."""
    parts = []
    for s in strs:
        mant, _, ex = s.lower().partition("e")
        e = int(ex) if ex else 0
        if mant.count(".") > 1 or not mant.replace(".", "").isdigit() or (ex and not ex.lstrip("+-").isdigit()):
            raise ValueError("not a plain non-negative decimal: %r" % s)
        whole, _, frac = mant.partition(".")
        parts.append((int(whole + frac), len(frac) - e))
    F = max(0, max(d for _, d in parts))
    ints = [m * 10 ** (F - d) for m, d in parts]
    g = 0
    for v in ints:
        g = math.gcd(g, v)
    g = g or 1
    return [v // g for v in ints], Fraction(10 ** F, g)   # divided by their gcd: the primitive integer heights


def main():
    args = sys.argv[1:]
    kind, path = args[0], args[1]
    opts = {}
    i = 2
    while i < len(args):
        if args[i] == "--bump":
            opts["bump"] = (int(args[i + 1]), args[i + 2]); i += 3
        elif args[i] == "--move":
            opts["move"] = tuple(int(x) for x in args[i + 1:i + 4]); i += 4
        else:
            opts[args[i][2:]] = args[i + 1]; i += 2
    t0 = time.time()
    ok = True

    def check(desc, cond, detail=""):
        nonlocal ok
        ok &= bool(cond)
        print("B: check %s: %s%s" % (desc, "PASS" if cond else "FAIL", (" (%s)" % detail) if detail else ""))

    if kind == "1a":
        vals, scale = scale_to_integers(decimal_lines(path))
        if "bump" in opts:
            bi, bd = opts["bump"]
            vals[bi] += int(bd)
            print("B: CONTROL height %d changed by %s (units of the primitive integer heights)" % (bi, bd))
        if min(vals) < 0:
            check("every height >= 0", False)
            print("B: VERDICT REFUSED")
            return 1
        header = "1a %d\n" % len(vals)
    else:
        with open(path) as fh:
            d = json.load(fh)
        n, Q, vals = d["n"], d["Q"], list(d["heights_int"])
        if "bump" in opts:
            bi, bd = opts["bump"]
            vals[bi] += int(bd)
            print("B: CONTROL height %d changed by %s (units of 1/Q)" % (bi, bd))
        if "move" in opts:
            mi, mj, md = opts["move"]
            vals[mi] += md
            vals[mj] -= md
            print("B: CONTROL %d/Q moved from height %d to height %d (mass unchanged)" % (md, mj, mi))
        check("n = number of heights", n == len(vals), "%d heights" % len(vals))
        check("every k_i in [0, Q]", all(0 <= v <= Q for v in vals))
        check("total mass exactly 1 (2 sum k_i = n Q)", 2 * sum(vals) == n * Q)
        if not ok:
            print("B: VERDICT REFUSED")
            return 1
        header = "1b %d %d\n" % (n, Q)
    exe = build()
    inp = header + "\n".join(str(v) for v in vals) + "\n"
    r = subprocess.run([exe], input=inp.encode(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    out = r.stdout.decode()
    for line in out.splitlines():
        if not line.startswith("B: max") and not line.startswith("B: argmax") and not line.startswith("B: S "):
            print(line)
    if r.returncode != 0 or "VERDICT REFUSED" in out:
        print("B: the C program refused the input: %s" % r.stderr.decode().strip())
        print("B: VERDICT REFUSED")
        return 1
    got = {ln.split()[1]: ln.split()[2:] for ln in out.splitlines() if ln.startswith("B: ")}
    vmax = int(got["max"][0], 16)
    arg = [int(x) for x in got["argmax"]]
    S = int(got["S"][0])
    n = len(vals)
    if kind == "1a":
        value = Fraction(2 * n * vmax, S * S)
        print("B: S = %d\nB: max_k c_k = %d at k = %s" % (S, vmax, ",".join(map(str, arg))))
    else:
        value = Fraction(2 * vmax, n * Q * Q)
        print("B: max_m v_m = %d at m = %s" % (vmax, ",".join(map(str, arg))))
    q = value.numerator * 10 ** 40 // value.denominator
    s = str(q).rjust(41, "0")
    print("B: exact value %d/%d\nB: = %s.%s..." % (value.numerator, value.denominator, s[:-40], s[-40:]))
    if "claim" in opts:
        check("the exact value equals the stated fraction %s" % opts["claim"], value == Fraction(opts["claim"]))
    if "prefix" in opts:
        nd = len(opts["prefix"].split(".")[1])
        tq = value.numerator * 10 ** nd // value.denominator
        whole, frac = divmod(tq, 10 ** nd)
        check("the exact value truncated to %d decimals is the stated %s" % (nd, opts["prefix"]),
              "%d.%s" % (whole, str(frac).rjust(nd, "0")) == opts["prefix"])
    if "argmax" in opts:
        check("the maximum is attained at the stated index %s" % opts["argmax"], int(opts["argmax"]) in arg)
    if "printed" in opts:
        pr = Fraction(opts["printed"])
        check("the printed upper bound %s is at least the exact value" % opts["printed"], pr >= value,
              "printed - exact = %.3e" % float(pr - value))
    print("B: run time %.1f s" % (time.time() - t0))
    print("B: VERDICT %s" % ("PASS" if ok else "REFUSED"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
