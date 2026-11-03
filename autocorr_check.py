#!/usr/bin/env python3
"""Independent exact check of step-function certificates for C_1c and C_1d (teorth/optimizationproblems).

Written from the definitions on the constant pages, without reading or running the certificates' own verifiers.

  autocorr_check.py 1c FILE    lower bound C_1c >= K * min_{0<=k<=K} A_k / S^2,  A_k = sum_i w_i w_{i+k}
  autocorr_check.py 1d FILE    upper bound C_1d <= 2n * max_m |c_m| / S^2,      c = v * v (autoconvolution)

1c: f = sum_i w_i 1_[i/K,(i+1)/K) with non-negative integer w; the autocorrelation is linear between multiples of
1/K, so its minimum over t in [0,1] is (1/K) min_{0<=k<=K} A_k and int f = S/K.
1d: n equal cells on [-1/4,1/4] (h = 1/(2n)), signed values v; f*f is linear between multiples of h, its value at
the knot 2a + m h is h * sum_{i+j=m-1} v_i v_j, so max_t |f*f| / (int f)^2 = max |c| / (h S^2).

Method, chosen to share nothing with the usual limb-split numpy scorers:
1. every A_k (or c_m) at once from one big-integer product (Kronecker substitution; for signed v the positive and
   negative parts are packed separately, c = pp - 2pq + qq);
2. the whole sequence is confirmed by a polynomial identity at 4 random points modulo the prime 2^61 - 1
   (sum_k A_k r^k = W(r) W~(r), resp. sum_m c_m r^m = V(r)^2), which a wrong coefficient fails with probability
   at most 2n/2^61 per point;
3. the deciding lag (the minimising k, resp. the maximising m) is summed again directly in O(n);
4. the bound is the reduced fraction; it is compared with the claimed fraction or decimal when the file has one.
Values given as JSON floats are read exactly as binary64 (fractions.Fraction(float)), or with --decimal as the
exact decimals printed in the file.
"""
import json
import random
import sys
from fractions import Fraction
from math import gcd

P61 = (1 << 61) - 1


def as_ints(vals):
    """Integers from JSON ints, or floats read exactly (scaled to a common power-of-two denominator)."""
    if all(isinstance(v, int) and not isinstance(v, bool) for v in vals):
        return [int(v) for v in vals], 1
    fr = [v if isinstance(v, Fraction) else Fraction(v) if isinstance(v, (int, float)) else Fraction(str(v)) for v in vals]
    den = 1
    for f in fr:
        den = den * f.denominator // gcd(den, f.denominator)
    return [int(f * den) for f in fr], den


def pack(seq, slot):
    return int.from_bytes(b"".join(int(x).to_bytes(slot // 8, "little") for x in seq), "little")


def unpack(big, slot, count):
    raw = big.to_bytes(slot // 8 * count + slot // 8, "little")
    w = slot // 8
    return [int.from_bytes(raw[i * w:(i + 1) * w], "little") for i in range(count)]


def slot_bits(bound):
    b = bound.bit_length() + 2
    return (b + 7) // 8 * 8


def autocorrelation(w):
    """A_k = sum_i w_i w_{i+k} for k = 0..n-1 (w non-negative)."""
    n = len(w)
    slot = slot_bits(n * max(w) ** 2)
    prod = pack(w, slot) * pack(w[::-1], slot)  # coefficient n-1+k of the product is A_k
    coef = unpack(prod, slot, 2 * n - 1)
    return coef[n - 1:]


def autoconvolution(v):
    """c_m = sum_{i+j=m} v_i v_j for m = 0..2n-2 (v signed)."""
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


def identity_check(w, coeffs, kind, rounds=4):
    rng = random.SystemRandom()
    for _ in range(rounds):
        r = rng.randrange(2, P61 - 1)
        if kind == "1c":  # A(r) * r^(n-1) + lower half = W(r) * W_rev(r); check the full product instead
            full_left = horner(w, r) * horner(w[::-1], r) % P61
            n = len(w)
            # full product coefficients: low part (k < n-1) mirrors A by symmetry: coefficient n-1-k = A_k
            full = list(reversed(coeffs[1:])) + coeffs
            if horner(full, r) != full_left:
                return False
        else:
            if horner(coeffs, r) != horner(w, r) ** 2 % P61:
                return False
    return True


def check_1c(d):
    K = int(d["K"])
    w, den = as_ints(d["weights_int"] if "weights_int" in d else d["weights"])
    assert den == 1, "1c weights must be integers"
    assert all(x >= 0 for x in w) and any(w), "weights must be non-negative and not all zero"
    n = len(w)
    A = autocorrelation(w)
    A += [0] * max(0, K + 1 - n)  # A_k = 0 for k >= n
    assert identity_check(w, A[:n], "1c"), "polynomial identity failed"
    lags = A[:K + 1]
    amin = min(lags)
    kmin = lags.index(amin)
    direct = sum(w[i] * w[i + kmin] for i in range(n - kmin))
    assert direct == amin, "direct sum at the minimising lag disagrees"
    S = sum(w)
    bound = Fraction(K * amin, S * S)
    return {"n": n, "K": K, "S": S, "argmin_k": kmin, "min_A": amin, "bound": bound,
            "ties": sum(1 for x in lags if x == amin)}


def check_1d(d):
    vals = next(d[k] for k in ("values", "v_int", "v", "cells", "construction") if k in d)
    v, den = as_ints(vals)
    n = len(v)
    S = sum(v)
    assert S != 0, "integral is zero"
    c = autoconvolution(v)
    assert identity_check(v, c, "1d"), "polynomial identity failed"
    cmax = max(abs(x) for x in c)
    m = next(i for i, x in enumerate(c) if abs(x) == cmax)
    direct = sum(v[i] * v[m - i] for i in range(max(0, m - n + 1), min(m, n - 1) + 1))
    assert abs(direct) == cmax, "direct sum at the maximising knot disagrees"
    bound = Fraction(2 * n * cmax, S * S)  # the common denominator den cancels: (c/den^2) / (S/den)^2
    return {"n": n, "S": S, "argmax_m": m, "max_abs_c": cmax, "bound": bound, "scaled_by": den,
            "sign_at_max": "+" if c[m] > 0 else "-", "ties": sum(1 for x in c if abs(x) == cmax)}


def controls(kind, d):
    """Red controls: a corrupted coefficient must fail the identity, a changed weight must change the bound or the
    identity, and scaling every value by 3 must leave the bound unchanged."""
    key = "weights_int" if kind == "1c" else next(k for k in ("values", "v_int", "v", "cells") if k in d)
    vals, _ = as_ints(d[key])
    seq = autocorrelation(vals)[:len(vals)] if kind == "1c" else autoconvolution(vals)
    j = len(seq) // 3
    bad = seq[:j] + [seq[j] + 1] + seq[j + 1:]
    out = ["corrupted coefficient %d rejected: %s" % (j, not identity_check(vals, bad, kind))]
    i = len(vals) // 2
    moved = vals[:i] + [vals[i] + 1] + vals[i + 1:]
    out.append("old coefficients rejected for weight %d + 1: %s" % (i, not identity_check(moved, seq, kind)))
    base = (check_1c if kind == "1c" else check_1d)(d)["bound"]
    scaled = dict(d)
    scaled[key] = [3 * x for x in vals]
    out.append("bound unchanged when every value is scaled by 3: %s" %
               ((check_1c if kind == "1c" else check_1d)(scaled)["bound"] == base))
    return out


def dec(fr, digits, up=False):
    """fr rounded down (or up) to `digits` decimals, as a string."""
    q, r = divmod(fr.numerator * 10 ** digits, fr.denominator)
    q += 1 if up and r else 0
    s = str(q).rjust(digits + 1, "0")
    return s[:-digits] + "." + s[-digits:]


def main():
    kind, path = sys.argv[1], sys.argv[2]
    # --decimal reads JSON floats as the exact decimals printed in the file instead of as binary64
    d = json.load(open(path), parse_float=Fraction) if "--decimal" in sys.argv else json.load(open(path))
    if isinstance(d, list):
        d = {"values": d}
    r = check_1c(d) if kind == "1c" else check_1d(d)
    b = r["bound"]
    print(f"{kind} {path}")
    for k, val in r.items():
        if k != "bound":
            print(f"  {k}: {val}")
    print(f"  bound = {b.numerator}/{b.denominator}")
    print(f"  decimal: floor {dec(b, 22)}  ceil {dec(b, 22, up=True)}")
    if "score_num" in d and "score_den" in d:
        claimed = Fraction(int(d["score_num"]), int(d["score_den"]))
        print("  claimed fraction in the file:", "AGREES" if claimed == b else f"DIFFERS ({claimed})")
    if "--controls" in sys.argv:
        for line in controls(kind, d):
            print("  control: " + line)


if __name__ == "__main__":
    main()
