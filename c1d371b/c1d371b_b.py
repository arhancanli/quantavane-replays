#!/usr/bin/env python3
"""Program B for the 1d step function of PR #371 (its eighth certificate, 1.4546676818): the C core c1d371b_b.c (its own reader, plain loops over every knot
in signed 128-bit integers), then the ratio 2n max|c_m| / S^2 compared with the claims by integer cross-multiplication.

  c1d371b_b.py FILE [--claim P/Q] [--argmax M] [--printed DEC] [--bump I D]

Shares no code with program A (c1d371b_a.py). The binary is built from c1d371b_b.c into bin/ if missing.
Standard library only; run with python3 -I.
"""
import os
import subprocess
import sys
import time
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "bin", "c1d371b_b")


def build():
    if not os.path.exists(BIN):
        os.makedirs(os.path.dirname(BIN), exist_ok=True)
        subprocess.check_call(["cc", "-O2", "-o", BIN, os.path.join(HERE, "c1d371b_b.c")])


def dec_to_ratio(s):
    """An exact decimal string as (numerator, denominator)."""
    neg = s.startswith("-")
    s = s.lstrip("+-")
    if "." in s:
        a, b = s.split(".")
    else:
        a, b = s, ""
    num, den = int(a + b), 10 ** len(b)
    return (-num if neg else num), den


def main():
    t0 = time.time()
    path = sys.argv[1]
    claim = argmax = printed = None
    bump = []
    k = 2
    while k < len(sys.argv):
        a = sys.argv[k]
        if a == "--claim":
            p, q = sys.argv[k + 1].split("/"); claim = (int(p), int(q)); k += 2
        elif a == "--argmax":
            argmax = int(sys.argv[k + 1]); k += 2
        elif a == "--printed":
            printed = sys.argv[k + 1]; k += 2
        elif a == "--bump":
            bump = [sys.argv[k + 1], sys.argv[k + 2]]; k += 3
        else:
            raise SystemExit("unknown option " + a)
    build()
    out = subprocess.run([BIN, path] + bump, capture_output=True, text=True)
    if out.returncode != 0:
        print("core error: " + out.stderr.strip()); print("VERDICT REFUSED"); return
    kv = {}
    for line in out.stdout.splitlines():
        if line.startswith("control "):
            print("CONTROL " + line[len("control "):]); continue
        key, _, rest = line.partition(" ")
        kv[key] = rest
    n, S, M = int(kv["n"]), int(kv["S"]), int(kv["maxabs"])
    arg = [int(x) for x in kv["argmax"].split()]
    print("cells n = %d, negative heights %s, largest |v_i| = %s, S = %d" % (n, kv["negative"], kv["vmax"], S))
    if S == 0:
        print("VERDICT REFUSED (S = 0)"); return
    num, den = 2 * n * M, S * S
    g = gcd(num, den)
    num, den = num // g, den // g
    print("max |c_m| = %d at m = %s (c_m = %s at the first)" % (M, ",".join(map(str, arg)), kv["valueatfirstmax"]))
    print("fingerprint sum c_m 3^m mod 2^61-1 = %s" % kv["fingerprint"])
    print("exact value %d/%d" % (num, den))
    ip, r = divmod(num, den)
    ds = []
    for _ in range(30):                          # long division, truncated
        r *= 10
        d, r = divmod(r, den)
        ds.append(str(d))
    print("decimal (truncated) %d.%s" % (ip, "".join(ds)))
    ok = True
    if claim is not None:
        same = num * claim[1] == den * claim[0]
        print("stated fraction: %s" % ("equal" if same else "DIFFERENT"))
        ok = ok and same
    if argmax is not None:
        same = arg == [argmax]
        print("stated maximizing knot m = %d only: %s" % (argmax, "yes" if same else "NO"))
        ok = ok and same
    if printed is not None:
        pn, pd = dec_to_ratio(printed)
        holds = num * pd <= pn * den
        gap = (pn * den - num * pd) / (pd * den)
        print("printed upper bound %s: %s (printed minus exact = %.3e)" % (printed, "holds" if holds else "FAILS", gap))
        ok = ok and holds
    print("run time %.1f s" % (time.time() - t0))
    print("VERDICT %s" % ("PASS" if ok else "REFUSED"))


if __name__ == "__main__":
    main()
