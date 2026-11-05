#!/usr/bin/env python3
"""Exact replay of the finite k = 2 controlled-carry bound for the Gyarmati-Hennecart-Ruzsa sum-difference constant
(3a page, [Syl2026]: C_3a > 1.19281555*), from the theorem note's data alone.

  1. c3a_hist.c builds M = <G> cap [0, B] and the min-cost histogram of all Q = q^2 residues (each residue is a min
     over at most four digit pairs, since B < q); a full scatter of D x D gives the same histogram on small cases.
  2. Exactly, in integers and rationals: Z_2(x) = sum_c m(c) x^c at the dyadic x, the type counts
     N_r = floor(l x^c(r) / Z_2(x)) for r not 0 or h = Q/2, N_0 = l - sum N_r, and C = sum_r N_r c(r).
  3. In interval arithmetic (decimal, directed rounding, correctly rounded logarithms widened by one unit):
     log of the multinomial l! / (N_0! prod N_r!) from below with Robbins' bounds on ln n!, log P_+(x) from above,
     and theta = 1 + (log M + C log x - 2 l log P_+(x)) / log(4 B q^(2l) / (q - 1) + 1) from below.

Usage: python3 -I c3a_check.py syl2 [CERTS] [--controls] [--workers=N]   Python 3.9+, standard library only, and cc.
(CERTS/syl2026/k2_cost_hist.tsv, the package's histogram fetched by fetch.sh, is compared byte for byte if present.)
"""
import decimal
import hashlib
import os
import subprocess
import sys
import time
from decimal import Decimal
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "bin")
EXE = os.path.join(BIN, "c3a_hist")
PREC = 90

# [Syl2026] k = 2, as its theorem note states (xSyllaa/c3a-controlled-carry-k2 @ 5ca7cb8, proof/THEOREM.md)
SYL2 = dict(
    G=[1971, 2016, 2100, 2628, 2688, 2800], B=26972, q=27022, xhex="0x1.ffde1bafccc06p-1", X53=9004870224628742,
    ell=10 ** 18,
    hist_sha="4d16010f1b41d875c2912f397e76e8d4010e49b635e139ad1f9bcbd741956c63",   # the note's data/k2_cost_hist.tsv
    stated=dict(M=5869, MM=31274, costs=26540, cmax=53762, ch=24367, n0=3362808316958,
                C=26096544438312721237750),
    stated_bounds=dict(logM="19353442719252255618.41914067403", logU="15418306021682217276.01171204408",
                       logden="20408813256039808818.33351550874", Z2="297948.172635433171794",
                       Pplus="76.30733448716651847590", theta="1.192815557092986248985"),
    claim="1.19281555",
)


def build():
    os.makedirs(BIN, exist_ok=True)
    src = os.path.join(HERE, "c3a_hist.c")
    if not os.path.exists(EXE) or os.path.getmtime(EXE) < os.path.getmtime(src):
        subprocess.run(["cc", "-O2", "-std=c99", "-o", EXE, src], check=True)


def histogram(B, q, G, workers):
    procs = [subprocess.Popen([EXE, str(B), str(q), str(i), str(workers), *map(str, G)], stdout=subprocess.PIPE,
                              text=True) for i in range(workers)]
    hist, info = {}, {}
    for p in procs:
        out, _ = p.communicate()
        if p.returncode:
            raise RuntimeError("c3a_hist failed")
        for line in out.split("\n"):
            f = line.split()
            if not f:
                continue
            if f[0] == "c":
                hist[int(f[1])] = hist.get(int(f[1]), 0) + int(f[2])
            elif f[0] == "unreached":
                info["unreached"] = info.get("unreached", 0) + int(f[1])
            else:
                info[f[0]] = int(f[1])
    return hist, info


def types(hist, mprime, cmax, X, ell):
    """Exact Z_2 numerator (Z_2 = Znum / 2^(53 cmax)), type counts n(c) = floor(l x^c / Z_2), N_0 and C."""
    Znum = 0
    for c in range(cmax, -1, -1):           # Horner in X; Znum = sum_c m(c) X^c 2^(53 (cmax - c))
        Znum = Znum * X + (hist.get(c, 0) << (53 * (cmax - c)))
    counts, Xc = {}, 1
    for c in range(0, cmax + 1):            # n(c) = floor(l X^c 2^(53 (cmax - c)) / Znum)
        if c in mprime:
            counts[c] = ((ell * Xc) << (53 * (cmax - c))) // Znum
        Xc *= X
    N0 = ell - sum(mprime[c] * counts[c] for c in mprime)
    Ctot = sum(mprime[c] * counts[c] * c for c in mprime)
    return Znum, counts, N0, Ctot


# ---- interval arithmetic on decimals ------------------------------------------------------------------------------
LO = decimal.Context(prec=PREC, rounding=decimal.ROUND_FLOOR)
HI = decimal.Context(prec=PREC, rounding=decimal.ROUND_CEILING)
NEAR = decimal.Context(prec=PREC, rounding=decimal.ROUND_HALF_EVEN)


class I:
    """A closed interval [lo, hi] of decimals; every operation rounds outward."""

    def __init__(self, lo, hi=None):
        self.lo, self.hi = Decimal(lo), Decimal(lo if hi is None else hi)
        assert self.lo <= self.hi

    def __add__(self, o):
        o = o if isinstance(o, I) else I(o)
        return I(LO.add(self.lo, o.lo), HI.add(self.hi, o.hi))

    def __sub__(self, o):
        o = o if isinstance(o, I) else I(o)
        return I(LO.subtract(self.lo, o.hi), HI.subtract(self.hi, o.lo))

    def __mul__(self, o):
        o = o if isinstance(o, I) else I(o)
        ps = [(a, b) for a in (self.lo, self.hi) for b in (o.lo, o.hi)]
        return I(min(LO.multiply(a, b) for a, b in ps), max(HI.multiply(a, b) for a, b in ps))

    def __truediv__(self, o):
        o = o if isinstance(o, I) else I(o)
        assert o.lo > 0 or o.hi < 0, "division by an interval containing 0"
        ps = [(a, b) for a in (self.lo, self.hi) for b in (o.lo, o.hi)]
        return I(min(LO.divide(a, b) for a, b in ps), max(HI.divide(a, b) for a, b in ps))


def iln(x):
    """ln of a positive interval: decimal's ln is correctly rounded, so one unit either way encloses the true value."""
    assert x.lo > 0
    a, b = NEAR.ln(x.lo), NEAR.ln(x.hi)
    return I(NEAR.next_minus(a), NEAR.next_plus(b))


def ratio(num, den, digits=PREC):
    """The positive rational num/den as an interval of width 10^-digits, by integer division (no float)."""
    s = 10 ** digits
    f = num * s // den
    return I(LO.divide(Decimal(f), Decimal(s)), HI.divide(Decimal(f + 1), Decimal(s)))


def pi_interval():
    """pi = 16 atan(1/5) - 4 atan(1/239) in integers scaled by 10^(PREC+10); the alternating tails bound the error."""
    s = 10 ** (PREC + 10)

    def atan_inv(n):
        total, term, k, sign = 0, s // n, 1, 1
        while term:
            total += sign * (term // k)
            term //= n * n
            k += 2
            sign = -sign
        return total

    p = 16 * atan_inv(5) - 4 * atan_inv(239)
    err = 200   # truncations: a few units per term, far below this
    return I(LO.divide(Decimal(p - err), Decimal(s)), HI.divide(Decimal(p + err), Decimal(s)))


PI = None


def lnfact(n):
    """ln n! for an integer n >= 0: exact sums for n <= 20, Robbins' bounds above."""
    if n <= 20:
        f = 1
        for k in range(2, n + 1):
            f *= k
        return iln(I(f))
    ln_n = iln(I(n))
    base = I(n) * ln_n - I(n) + (iln(I(2) * PI * I(n)) * I("0.5"))
    return I(LO.add(base.lo, LO.divide(1, 12 * n + 1)), HI.add(base.hi, HI.divide(1, 12 * n)))


def check_syl2(root, controls, workers):
    global PI
    PI = pi_interval()
    r = SYL2
    st = r["stated"]
    t0 = time.time()
    B, q, G, ell = r["B"], r["q"], r["G"], r["ell"]
    Q = q * q
    print(f"Row syl2 ([Syl2026] k = 2, PR #238): G = {G}, B = {B}, q = {q}, Q = q^2 = {Q}, l = 10^18, x = {r['xhex']}")
    hist, info = histogram(B, q, G, workers)
    cmax, ch = max(hist), info["h"]
    facts = [("|M|", info["M"], st["M"]), ("|M+M|", info["MM"], st["MM"]), ("distinct costs", len(hist), st["costs"]),
             ("max cost", cmax, st["cmax"]), ("c(h), h = Q/2", ch, st["ch"])]
    ok = True
    for name, ours, theirs in facts:
        ok &= ours == theirs
        print(f"  {name}: {ours} (stated {theirs}) {'agrees' if ours == theirs else 'DIFFERS'}")
    good = info["unreached"] == 0 and sum(hist.values()) == Q and hist.get(0) == 1
    ok &= good
    print(f"  every one of the {Q} residues is reached, and only r = 0 has cost 0: {good} ({time.time() - t0:.0f} s)")
    ours = "# cost\tmultiplicity\n" + "".join(f"{c}\t{hist[c]}\n" for c in sorted(hist))
    h_ours = hashlib.sha256(ours.encode()).hexdigest()
    same = h_ours == r["hist_sha"]
    print(f"  our histogram in the package's layout: sha256 {h_ours} "
          f"({'byte-identical to its data/k2_cost_hist.tsv, as the note states' if same else 'DIFFERENT from the stated sha256'})")
    ok &= same
    if root and os.path.exists(os.path.join(root, "syl2026", "k2_cost_hist.tsv")):
        theirs = open(os.path.join(root, "syl2026", "k2_cost_hist.tsv"), "rb").read()
        print(f"  the fetched data file (sha256 {hashlib.sha256(theirs).hexdigest()[:8]}...) equals ours byte for byte: "
              f"{theirs == ours.encode()}")
        ok &= theirs == ours.encode()
    # exact: x = X / 2^53
    xn, xd = float.fromhex(r["xhex"]).as_integer_ratio()
    assert 2 ** 53 % xd == 0, xd
    X = xn * (2 ** 53 // xd)                  # x = X / 2^53 exactly
    assert 0 < X < 2 ** 53 and X == r["X53"], X
    print(f"  x = {r['xhex']} = {X} / 2^53 (as the note states)")
    mprime = {c: (m - 1 if c == ch else m) for c, m in hist.items() if c != 0}
    Znum, counts, N0, Ctot = types(hist, mprime, cmax, X, ell)
    nmin = min(counts[c] for c in mprime if mprime[c])
    print(f"  exact type counts: N_0 = {N0} (stated {st['n0']}) {'agrees' if N0 == st['n0'] else 'DIFFERS'}; "
          f"smallest N_r = {nmin} (every N_r >= 1: {nmin >= 1})")
    print(f"  C = sum N_r c(r) = {Ctot} (stated {st['C']}) {'agrees' if Ctot == st['C'] else 'DIFFERS'}; "
          f"even: {Ctot % 2 == 0}")
    ok &= N0 == st["n0"] and Ctot == st["C"] and nmin >= 1 and N0 >= 0 and Ctot % 2 == 0
    # interval bounds
    Z2 = ratio(Znum, (1 << (53 * cmax)))
    inM = [False] * (B + 1)
    inM[0] = True
    for n in range(1, B + 1):
        inM[n] = any(g <= n and inM[n - g] for g in G)
    Mset = [n for n in range(B + 1) if inM[n]]
    sums = set()
    for a in Mset:
        for b in Mset:
            if b >= a:
                sums.add(a + b)
    MM = sorted(sums)
    smax = MM[-1]
    Pnum, sset = 0, set(MM)
    for s in range(smax, -1, -1):
        Pnum = Pnum * X + ((1 << (53 * (smax - s))) if s in sset else 0)
    Pp = ratio(Pnum, (1 << (53 * smax)))
    logM = lnfact(ell) - lnfact(N0)
    acc = I(0)
    for c, m in mprime.items():
        if m:
            acc = acc + lnfact(counts[c]) * I(m)
    logM = logM - acc
    lnx = iln(I(X)) - iln(I(2)) * I(53)
    logU = I(-Ctot) * lnx + I(2 * ell) * iln(Pp)          # log(x^-C P_+^(2l)), from above
    lnq = iln(I(q))
    logden = I(2 * ell) * lnq + iln(I(4 * B) / I(q - 1)) + I(0, Decimal("1e-100"))
    num = logM - logU
    theta = I(1) + num / logden
    sb = r["stated_bounds"]
    print("  interval bounds (ours, then the note's):")
    for name, ours, side, theirs in (("log M", logM, "lo", sb["logM"]), ("log(x^-C P+^2l)", logU, "hi", sb["logU"]),
                                     ("log(4Bq^2l/(q-1)+1)", logden, "hi", sb["logden"]), ("Z_2(x)", Z2, "lo", sb["Z2"]),
                                     ("P_+(x)", Pp, "hi", sb["Pplus"]), ("theta", theta, "lo", sb["theta"])):
        v = ours.lo if side == "lo" else ours.hi
        print(f"    {name:22s} {'>=' if side == 'lo' else '<='} {NEAR.create_decimal(v):.40}  | {theirs}")
    holds = theta.lo > Decimal(r["claim"])
    print(f"  theta >= {NEAR.create_decimal(theta.lo):.25}: C_3a > {r['claim']} {'holds' if holds else 'NOT shown'}")
    ok &= holds
    if controls:
        print("Red controls:")
        c1 = not theta.lo > Decimal("1.19281556")
        print(f"  C1 the claim one unit higher in the eighth decimal (1.19281556): {'not shown, as expected' if c1 else 'SHOWN (control failed)'}")
        h2, i2 = histogram(B, q, [g for g in G if g != 2800], workers)
        c2 = len(h2) != len(hist) or i2["M"] != info["M"]
        print(f"  C2 the generator 2800 dropped: |M| = {i2['M']}, {len(h2)} distinct costs: "
              f"{'the stated data no longer match, as expected' if c2 else 'SAME (control failed)'}")
        _, _, N0b, Cb = types(hist, mprime, cmax, X - 1, ell)
        c3 = N0b != st["n0"] and Cb != st["C"]
        print(f"  C3 x moved down by one unit in its last bit: N_0 = {N0b}, C = {Cb}: "
              f"{'neither matches the note, as expected' if c3 else 'MATCHES (control failed)'}")
        ok &= c1 and c2 and c3
    print(f"RESULT: syl2 {'PASS' if ok else 'FAIL'}")
    print(f"run time {time.time() - t0:.0f} s")
    return ok


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    workers = next((int(a.split("=")[1]) for a in sys.argv[1:] if a.startswith("--workers=")),
                   min(8, os.cpu_count() or 1))
    if not args or args[0] != "syl2":
        print(__doc__)
        sys.exit(2)
    build()
    root = args[1] if len(args) > 1 else os.path.join(HERE, "certs")
    sys.exit(0 if check_syl2(root, "--controls" in sys.argv, workers) else 1)


if __name__ == "__main__":
    main()
