#!/usr/bin/env python3
"""Exact check of a step-function certificate for C_1a or C_1d (upper bounds) or C_1c (lower bound), with red controls.

Definitions (teorth/optimizationproblems, constants/1c.md and 1d.md):
  C_1c: min_{0<=t<=1} int f(x) f(x+t) dx <= C_1c (int f)^2 for every non-negative f.  For the step function
        f = sum_i w_i 1_[i/K,(i+1)/K) with integer w_i >= 0 the autocorrelation is linear between multiples of 1/K,
        so  C_1c >= K * min_{0<=k<=K} A_k / S^2,  A_k = sum_i w_i w_{i+k},  S = sum_i w_i.
  C_1a: max_t f*f(t) >= C_1a (int f)^2 for every non-negative f supported on [-1/4,1/4]; the same knot formula as
        C_1d below, with every v_i >= 0 (so every c_m >= 0), gives C_1a <= 2n * max_m c_m / S^2.
  C_1d: max_t |f*f(t)| >= C_1d (int f)^2 for every real f supported on [-1/4,1/4].  For n equal cells (h = 1/(2n))
        with integer values v_i, f*f is linear between multiples of h with knot values h * c_m, c = v * v, so
        C_1d <= 2n * max_m |c_m| / S^2,  S = sum_i v_i.

Method (Python standard library only, every number an integer or a Fraction):
  1. all A_k (all c_m) at once from big-integer products (Kronecker substitution; signed values are split into
     positive and negative parts, c = pp - 2pq + qq);
  2. the whole sequence is confirmed by its generating-polynomial identity at random points modulo 2^61 - 1;
  3. the deciding lag (knot) is summed again directly;
  4. the bound is the reduced fraction; it must equal the fraction stored in the certificate, imply the stated
     decimal claim, and beat the record values given on the command line (strictly).
Red controls (--controls): one value of the certificate changed by one unit must be refused (its fraction differs
from the stated one), a value changed by 1% must be refused, and the claim moved one unit of its last digit in the
improving direction must be refused.

Usage: python3 -I check.py CERT.json [--record DECIMAL_OR_FRACTION ...] [--controls]
"""
import json
import random
import sys
import time
from fractions import Fraction

P61 = (1 << 61) - 1


def pack(seq, slot):
    nb = slot // 8
    return int.from_bytes(b"".join(int(x).to_bytes(nb, "little") for x in seq), "little")


def unpack(big, slot, count):
    nb = slot // 8
    raw = big.to_bytes(nb * (count + 1), "little")
    return [int.from_bytes(raw[i * nb:(i + 1) * nb], "little") for i in range(count)]


def slot_bits(bound):
    return (bound.bit_length() + 2 + 7) // 8 * 8


def autocorrelation(w):
    n = len(w)
    slot = slot_bits(n * max(w) ** 2)
    coef = unpack(pack(w, slot) * pack(w[::-1], slot), slot, 2 * n - 1)
    return coef[n - 1:]  # coefficient n-1+k of W(x) W_rev(x) is A_k


def autoconvolution(v):
    n = len(v)
    p = [x if x > 0 else 0 for x in v]
    q = [-x if x < 0 else 0 for x in v]
    m = max(max(p), max(q))
    slot = slot_bits(n * m * m)
    P, Q = pack(p, slot), pack(q, slot)
    pp, pq, qq = (unpack(x, slot, 2 * n - 1) for x in (P * P, P * Q, Q * Q))
    return [a - 2 * b + c for a, b, c in zip(pp, pq, qq)]


def horner(seq, r):
    acc = 0
    for x in reversed(seq):
        acc = (acc * r + x) % P61
    return acc


def identity_ok(kind, vals, coeffs, rounds=3):
    rng = random.SystemRandom()
    for _ in range(rounds):
        r = rng.randrange(2, P61 - 1)
        if kind == "1c":
            full = list(reversed(coeffs[1:])) + coeffs  # A_{-k} = A_k: the full product W(x) W_rev(x)
            if horner(full, r) != horner(vals, r) * horner(vals[::-1], r) % P61:
                return False
        elif horner(coeffs, r) != horner(vals, r) ** 2 % P61:
            return False
    return True


def evaluate(cert):
    """Exact bound of the certificate: (kind, Fraction, details)."""
    if "weights_int" in cert:
        K = int(cert["K"])
        w = [int(x) for x in cert["weights_int"]]
        if any(x < 0 for x in w) or not any(w):
            raise ValueError("1c weights must be non-negative integers, not all zero")
        n = len(w)
        if n < K + 1:
            w = w + [0] * (K + 1 - n)
            n = len(w)
        A = autocorrelation(w)
        if not identity_ok("1c", w, A):
            raise ValueError("polynomial identity failed")
        lags = A[:K + 1]
        amin = min(lags)
        k = lags.index(amin)
        if sum(w[i] * w[i + k] for i in range(n - k)) != amin:
            raise ValueError("direct sum at the deciding lag disagrees")
        S = sum(w)
        return "1c", Fraction(K * amin, S * S), {"K": K, "n": n, "S": S, "deciding lag k": k, "min A_k": amin,
                                                 "lags at the minimum": lags.count(amin)}
    v = [int(x) for x in cert["v_int"]]
    kind = "1a" if cert.get("constant") == "C_1a" else "1d"
    if kind == "1a" and any(x < 0 for x in v):
        raise ValueError("C_1a needs a non-negative f")
    n = len(v)
    S = sum(v)
    if S == 0:
        raise ValueError("integral is zero")
    c = autoconvolution(v)
    if not identity_ok("1d", v, c):
        raise ValueError("polynomial identity failed")
    cmax = max(abs(x) for x in c)
    m = next(i for i, x in enumerate(c) if abs(x) == cmax)
    if abs(sum(v[i] * v[m - i] for i in range(max(0, m - n + 1), min(m, n - 1) + 1))) != cmax:
        raise ValueError("direct sum at the deciding knot disagrees")
    return kind, Fraction(2 * n * cmax, S * S), {"n": n, "S": S, "deciding knot m": m, "max |c_m|": cmax,
                                                 "knots at the maximum": sum(1 for x in c if abs(x) == cmax)}


def dec(fr, digits, up=False):
    q, r = divmod(fr.numerator * 10 ** digits, fr.denominator)
    q += 1 if (up and r) else 0
    s = str(q).rjust(digits + 1, "0")
    return s[:-digits] + "." + s[-digits:]


def verdict(cert, records, quiet=False):
    """(passed, lines). Refuses on any mismatch; the claim must be implied and strictly beat every record."""
    out = []
    try:
        kind, B, det = evaluate(cert)
    except ValueError as e:
        return False, ["REFUSED: " + str(e)]
    out.append(f"constant C_{kind}; " + "; ".join(f"{k} = {v}" for k, v in det.items() if k != "S"))
    out.append(f"exact bound = {B.numerator}/{B.denominator}")
    out.append(f"            = {dec(B, 25)}... (floor to 25 decimals)")
    stated = Fraction(int(cert["score_num"]), int(cert["score_den"]))
    if stated != B:
        return False, out + ["REFUSED: the exact value differs from the fraction stated in the certificate"]
    out.append("the stated fraction: EQUAL")
    claim = Fraction(cert["claim"])
    ok = claim <= B if kind == "1c" else claim >= B
    if not ok:
        return False, out + [f"REFUSED: the claim {cert['claim']} is not implied"]
    out.append(f"claim C_{kind} {'>=' if kind == '1c' else '<='} {cert['claim']}: implied")
    for r in records:
        R = Fraction(r)
        better = B > R if kind == "1c" else B < R
        if not better:
            return False, out + [f"REFUSED: not strictly better than the record {r}"]
        cl = claim > R if kind == "1c" else claim < R
        out.append(f"record {r}: beaten strictly, margin {float(abs(B - R)):.6e} "
                   f"(the printed claim also beats it: {cl})")
    return True, out


def main():
    args = sys.argv[1:]
    path = args[0]
    records = []
    i = 1
    while i < len(args):
        if args[i] == "--record":
            records.append(args[i + 1])
            i += 2
        else:
            i += 1
    t0 = time.time()
    cert = json.load(open(path))
    ok, lines = verdict(cert, records)
    print(f"certificate {path}")
    for ln in lines:
        print("  " + ln)
    print(f"  RESULT: {'PASS' if ok else 'FAIL'} ({time.time() - t0:.0f} s)")
    if "--controls" in args:
        key = "weights_int" if "weights_int" in cert else "v_int"
        vals = cert[key]
        j = len(vals) // 3
        while vals[j] == 0:
            j += 1
        kind = "1c" if key == "weights_int" else ("1a" if cert.get("constant") == "C_1a" else "1d")
        print("red controls (each must be refused):")
        c1 = dict(cert, **{key: vals[:j] + [vals[j] + 1] + vals[j + 1:]})
        ok1, l1 = verdict(c1, records)
        print(f"  C1 value {j} raised by one unit: {'REFUSED' if not ok1 else 'ACCEPTED (control failed)'} - {l1[-1]}")
        big = max(range(len(vals)), key=lambda t: abs(vals[t]))
        c2 = dict(cert, **{key: vals[:big] + [vals[big] * 101 // 100] + vals[big + 1:]})
        c2["score_num"], c2["score_den"] = cert["score_num"], cert["score_den"]
        ok2, l2 = verdict(c2, records)
        print(f"  C2 largest value {big} raised by 1%: {'REFUSED' if not ok2 else 'ACCEPTED (control failed)'} - {l2[-1]}")
        # the same 1% change with its own (correct) fraction stated: the claim itself must then fail
        k2, B2, _ = evaluate(c2)
        c2b = dict(c2, score_num=str(B2.numerator), score_den=str(B2.denominator))
        ok2b, l2b = verdict(c2b, records)
        print(f"  C3 same change, its own fraction stated: {'REFUSED' if not ok2b else 'ACCEPTED (control failed)'}"
              f" - {l2b[-1]}")
        cs = cert["claim"]
        d = len(cs.split(".")[1])
        step = Fraction(1, 10 ** d)
        tighter = Fraction(cs) + step if kind == "1c" else Fraction(cs) - step
        c4 = dict(cert, claim=dec(tighter, d))
        ok4, l4 = verdict(c4, records)
        print(f"  C4 claim tightened by one unit of its last digit ({c4['claim']}): "
              f"{'REFUSED' if not ok4 else 'ACCEPTED (control failed)'} - {l4[-1]}")
        if not ok or ok1 or ok2 or ok2b or ok4:
            print("  CONTROLS: FAIL")
            sys.exit(1)
        print("  CONTROLS: all refused, as required")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
