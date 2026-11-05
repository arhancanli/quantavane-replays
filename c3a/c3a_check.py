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
       python3 -I c3a_check.py syl3 [CERTS] [--controls] [--workers=N] [--parts=DIR]   the k = 3 instance (1.19359680*):
       its q^3 histogram by c3a_hist3.c (about 6 CPU hours; --parts reads the 64 part outputs of an earlier run)
       python3 -I c3a_check.py m2026|m2026b|k2026b [CERTS] [--controls]   [M2026]'s finite bounds in the odd base
       q = 27995 (Theorems 1-4 of its README, rebuilt from the statements), and [K2026b]'s 1.19102809 as a control.
(CERTS/syl2026/k2_cost_hist.tsv, the package's histogram fetched by fetch.sh, is compared byte for byte if present.)
"""
import decimal
import hashlib
import json
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

SYL3 = dict(
    G=[1971, 2016, 2100, 2628, 2688, 2800], B=26972, q=27022, k=3, xhex="0x1.ffde3d0d5025dp-1", X53=9004879181054557,
    ell=10 ** 24,
    hist_sha="3186f6e220fd4a36369bf23ce841a9d54e6179cfaff00b5f3189c058f0d85c63",   # the note's k3/data/k3_cost_hist.tsv
    stated=dict(M=5869, MM=31274, costs=54099, cmax=80643, ch=24367, n0=5794476231064238,
                C=39290318538006294714264853344),
    stated_bounds=dict(logM="29091636209390427206933565.63284740439668", logU="23165014522649654717578745.16491695279380",
                       logden="30613219884059713225423555.68181075805065", Z="173196158.486920854972319981746",
                       Pplus="77.3059734705102016646138966528",
                       theta="1.19359680912973029490196359640645777794320100728768251031863"),
    claim="1.19359680",
)
EXE3 = os.path.join(BIN, "c3a_hist3")
NPARTS3 = 64

def build():
    os.makedirs(BIN, exist_ok=True)
    for src, exe in (("c3a_hist.c", EXE), ("c3a_hist3.c", EXE3)):
        src = os.path.join(HERE, src)
        if not os.path.exists(exe) or os.path.getmtime(exe) < os.path.getmtime(src):
            subprocess.run(["cc", "-O2", "-std=c99", "-o", exe, src], check=True)


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



def parse_hist(text, hist, info):
    for line in text.split("\n"):
        f = line.split()
        if not f:
            continue
        if f[0] == "c":
            hist[int(f[1])] = hist.get(int(f[1]), 0) + int(f[2])
        elif f[0] in ("unreached", "pairs"):
            info[f[0]] = info.get(f[0], 0) + int(f[1])
        else:
            info[f[0]] = int(f[1])


def histogram3(B, q, G, workers, parts_dir=None):
    """The k = 3 histogram as the sum of NPARTS3 row parts of c3a_hist3, run `workers` at a time, or read from
    parts_dir/part_i.txt if all parts are there (outputs of the same program and arguments)."""
    hist, info = {}, {}
    if parts_dir and all(os.path.exists(os.path.join(parts_dir, f"part_{i}.txt")) for i in range(NPARTS3)):
        for i in range(NPARTS3):
            parse_hist(open(os.path.join(parts_dir, f"part_{i}.txt")).read(), hist, info)
        return hist, info, "read from " + parts_dir
    pending, running = list(range(NPARTS3)), []
    while pending or running:
        while pending and len(running) < workers:
            i = pending.pop(0)
            running.append(subprocess.Popen([EXE3, str(B), str(q), str(i), str(NPARTS3), *map(str, G)],
                                            stdout=subprocess.PIPE, text=True))
        p = running.pop(0)
        out, _ = p.communicate()
        if p.returncode:
            raise RuntimeError("c3a_hist3 failed")
        parse_hist(out, hist, info)
    return hist, info, f"computed with {workers} workers"


def check_syl3(root, controls, workers, parts_dir):
    global PI
    PI = pi_interval()
    r = SYL3
    st = r["stated"]
    t0 = time.time()
    B, q, G, ell, k = r["B"], r["q"], r["G"], r["ell"], r["k"]
    Q = q ** k
    print(f"Row syl3 ([Syl2026] k = 3, PR #238): G = {G}, B = {B}, q = {q}, Q = q^3 = {Q}, l = 10^24, x = {r['xhex']}")
    ok = True
    # self-test on the real data: the folded rows equal the plain triple loop on a few rows
    for part in (0, 9007):
        a = subprocess.run([EXE3, str(B), str(q), str(part), str(q // 2 + 1), *map(str, G)], capture_output=True, text=True).stdout
        b = subprocess.run([EXE3, str(B), str(q), str(part), str(q // 2 + 1), *map(str, G), "--plain"], capture_output=True,
                           text=True).stdout
        ha, ia, hb, ib = {}, {}, {}, {}
        parse_hist(a, ha, ia)
        parse_hist(b, hb, ib)
        same = ha == hb and ia.get("unreached") == ib.get("unreached")
        ok &= same
        print(f"  self-test: rows {part} and {part + q // 2 + 1}, folded = plain triple loop: {same}")
    hist, info, how = histogram3(B, q, G, workers, parts_dir)
    cmax, ch = max(hist), info["h"]
    facts = [("|M|", info["M"], st["M"]), ("|M+M|", info["MM"], st["MM"]), ("distinct costs", len(hist), st["costs"]),
             ("max cost", cmax, st["cmax"]), ("c(h), h = Q/2", ch, st["ch"])]
    for name, ours, theirs in facts:
        ok &= ours == theirs
        print(f"  {name}: {ours} (stated {theirs}) {'agrees' if ours == theirs else 'DIFFERS'}")
    good = info["unreached"] == 0 and sum(hist.values()) == Q and hist.get(0) == 1
    ok &= good
    print(f"  histogram {how}: every one of the {Q} residues is reached, and only r = 0 has cost 0: {good} "
          f"({info['pairs']} folded (X_0, X_-1) pairs swept)")
    odd = sorted(c for c, m in hist.items() if m % 2)
    sym = odd == sorted({0, ch})
    ok &= sym
    print(f"  the odd cost classes are exactly those of 0 and h (c(-r) = c(r)): {sym} ({odd})")
    ours = "# cost\tmultiplicity\n" + "".join(f"{c}\t{hist[c]}\n" for c in sorted(hist))
    h_ours = hashlib.sha256(ours.encode()).hexdigest()
    same = h_ours == r["hist_sha"]
    print(f"  our histogram in the package's layout: sha256 {h_ours} "
          f"({'byte-identical to its k3/data/k3_cost_hist.tsv, as the note states' if same else 'DIFFERENT from the stated sha256'})")
    ok &= same
    fetched = os.path.join(root, "syl2026k3", "k3_cost_hist.tsv") if root else None
    if fetched and os.path.exists(fetched):
        theirs = open(fetched, "rb").read()
        print(f"  the fetched data file (sha256 {hashlib.sha256(theirs).hexdigest()[:8]}...) equals ours byte for byte: "
              f"{theirs == ours.encode()}")
        ok &= theirs == ours.encode()
    xn, xd = float.fromhex(r["xhex"]).as_integer_ratio()
    assert 2 ** 53 % xd == 0, xd
    X = xn * (2 ** 53 // xd)
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
    Zk = ratio(Znum, (1 << (53 * cmax)))
    Mset = semigroup_mask(G, B)
    MM = sorted({a + b for a in Mset for b in Mset})
    smax, sset, Pnum = MM[-1], set(MM), 0
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
    logU = I(-Ctot) * lnx + I(k * ell) * iln(Pp)          # log(x^-C P_+^(k l)), from above
    logden = I(k * ell) * iln(I(q)) + iln(I(4 * B) / I(q - 1)) + I(0, Decimal("1e-100"))
    theta = I(1) + (logM - logU) / logden
    sb = r["stated_bounds"]
    print("  interval bounds (ours, then the note's):")
    for name, val, side, theirs in (("log M", logM, "lo", sb["logM"]), ("log(x^-C P+^3l)", logU, "hi", sb["logU"]),
                                    ("log(4Bq^3l/(q-1)+1)", logden, "hi", sb["logden"]), ("Z_3(x)", Zk, "lo", sb["Z"]),
                                    ("P_+(x)", Pp, "hi", sb["Pplus"]), ("theta", theta, "lo", sb["theta"])):
        v = val.lo if side == "lo" else val.hi
        print(f"    {name:22s} {'>=' if side == 'lo' else '<='} {NEAR.create_decimal(v):.45}  | {theirs}")
    holds = theta.lo > Decimal(r["claim"])
    print(f"  theta >= {NEAR.create_decimal(theta.lo):.30}: C_3a > {r['claim']} {'holds' if holds else 'NOT shown'}")
    ok &= holds
    if controls:
        print("Red controls:")
        c1 = not theta.lo > Decimal("1.19359681")
        print(f"  C1 the claim one unit higher in the eighth decimal (1.19359681): {'not shown, as expected' if c1 else 'SHOWN (control failed)'}")
        _, _, N0b, Cb = types(hist, mprime, cmax, X - 1, ell)
        c2 = N0b != st["n0"] and Cb != st["C"]
        print(f"  C2 x moved down by one unit in its last bit: N_0 = {N0b}, C = {Cb}: "
              f"{'neither matches the note, as expected' if c2 else 'MATCHES (control failed)'}")
        h2 = dict(hist)
        h2[cmax] -= 1
        h2[cmax - 1] = h2.get(cmax - 1, 0) + 1
        c3 = sorted(c for c, m in h2.items() if m % 2) != sorted({0, ch})
        print(f"  C3 one residue moved from the top cost {cmax} to {cmax - 1}: the symmetry check "
              f"{'refuses it, as expected' if c3 else 'ACCEPTS it (control failed)'}")
        ok &= c1 and c2 and c3
    print(f"RESULT: syl3 {'PASS' if ok else 'FAIL'}")
    print(f"run time {time.time() - t0:.0f} s")
    return ok

# [M2026]: explicit depth-N masked-digit bounds in an odd base q with a pruned difference alphabet (the 3a row
# 1.19230254*, and its value on the [K2026] mask), and [K2026b]'s finite certificate 1.19102809 as a positive control.
GM = (1518, 1524, 1587, 2024, 2032, 2116)
M2026 = {
    "m2026": dict(B=17086, q=27995, factor=5, x=(9996683050229444, 10 ** 16), n0=17074745018, L=50031092045455301,
                  claim="1.1923025406835", stated="1.19230254068351294895", json="finite_certificate_17086_27995.json",
                  title="[M2026] primary: mask B = 17086, odd base q = 27995 (3a row 1.19230254*, comment 1.1923025406835)"),
    "m2026b": dict(B=17032, q=27995, factor=5, x=(9996683050229444, 10 ** 16), n0=17075622128, L=50029864065809249,
                   claim="1.1922990483828", stated="1.19229904838280822699", json="finite_certificate_17032_27995.json",
                   title="[M2026] on the [K2026] mask B = 17032, odd base q = 27995 (3a row comment 1.1922990483828)"),
    "k2026b": dict(B=17032, q=34065, factor=2, x=(9996789017186568, 10 ** 16), n0=14962802398, L=51500691976683641,
                   claim="1.19102809", stated="1.19102809424738956213", json="kleinwaks_finite_certificate_17032_34065.json",
                   title="positive control: [K2026b]'s finite certificate, B = 17032, carry-free base Q = 2B + 1 = 34065 "
                         "(3a row 1.19102809)"),
}


def semigroup_mask(G, B):
    """M = span_N0(G) cap [0, B], by reachability."""
    reach = bytearray(B + 1)
    reach[0] = 1
    for v in range(1, B + 1):
        reach[v] = 1 if any(v >= g and reach[v - g] for g in G) else 0
    return [v for v in range(B + 1) if reach[v]]


def kappa_table(M, B):
    """kappa(d) = min{a + b : a, b in M, a - b = d} for 0 <= d <= B (None if d is not in M - M), with a witness."""
    kap, wit = [None] * (B + 1), [None] * (B + 1)
    for i, a in enumerate(M):
        for b in M[:i + 1]:
            d, c = a - b, a + b
            if kap[d] is None or c < kap[d]:
                kap[d], wit[d] = c, (a, b)
    return kap, wit


def pruned_alphabet(kap, B, q):
    """One representative per complementary residue pair {r, q - r}, r = 1..(q-1)/2: the cheaper of r and q - r that
    lies in M - M (r on a tie). Returns [(r, d)]."""
    out = []
    for r in range(1, (q - 1) // 2 + 1):
        c = [d for d in (r, q - r) if d <= B and kap[d] is not None]
        if c:
            out.append((r, min(c, key=lambda d: (kap[d], d != r))))
    return out


def distinct_residues(alphabet, q):
    """True if {0} and +-d for d in the alphabet are pairwise distinct mod q (Theorem 1's condition)."""
    res = {0}
    for dd in alphabet:
        res.add(dd % q)
        res.add(-dd % q)
    return len(res) == 1 + 2 * len(alphabet)


def witness_ok(Mset, dd, a, b, kap):
    """(a, b) in M^2 with a - b = d and a + b = kappa(d)."""
    return a in Mset and b in Mset and a - b == dd and a + b == kap[dd]


def counts_ok(N, cnt):
    """Non-negative counts with n_0 = N - 2 sum n_d >= 0."""
    return all(v >= 0 for v in cnt.values()) and N - 2 * sum(cnt.values()) >= 0


def check_m2026(root, name, controls):
    global PI
    PI = pi_interval()
    t0 = time.time()
    row = M2026[name]
    B, q, N, fac = row["B"], row["q"], 10 ** 13, row["factor"]
    xn, xd = row["x"]
    print(f"Row {name}: {row['title']}; G = {GM}, N = 10^13, x = {xn}/{xd}")
    ok = True

    def check(what, cond, detail=""):
        nonlocal ok
        ok &= bool(cond)
        print(f"  check {what}: {'PASS' if cond else 'FAIL'}{' (' + detail + ')' if detail else ''}")

    d = os.path.join(root, "m2026") if root else None
    cert = None
    if d and os.path.exists(os.path.join(d, row["json"])):
        cert = json.load(open(os.path.join(d, row["json"])))   # compared only; nothing below is taken from it
    M = semigroup_mask(GM, B)
    S = sorted({a + b for a in M for b in M})
    check("hypotheses: q odd, B + 1 <= q, 0 in M, B in M", q % 2 == 1 and B + 1 <= q and M[0] == 0 and M[-1] == B,
          f"|M| = {len(M)}, |M + M| = {len(S)}, max(M + M) = {S[-1]}")
    if cert:
        check("|M|, |M + M| and max(M + M) as the certificate states",
              (len(M), len(S), S[-1]) == (cert["mask_size"], cert["sum_support_size"], cert["sum_support_max"]))
    if d and name == "m2026":
        theirs = [int(v) for v in open(os.path.join(d, "mask_17086.txt")).read().split()]
        check("our M equals the package's mask_17086.txt", theirs == M)
        theirs = [int(v) for v in open(os.path.join(d, "sum_support_17086.txt")).read().split()]
        check("our M + M equals the package's sum_support_17086.txt", theirs == S)
    kap, wit = kappa_table(M, B)
    Dp = pruned_alphabet(kap, B, q)
    check("the pruned alphabet {0, +-d} has pairwise distinct residues mod q", distinct_residues([e for _, e in Dp], q),
          f"|D+| = {len(Dp)}, |D| = {1 + 2 * len(Dp)}" + (", a complete residue system" if 1 + 2 * len(Dp) == q else ""))
    # counts: n_d = floor(N x^kappa(d) / Z(x)), Z(x) = 1 + 2 sum_{D+} x^kappa(d); n_0 = N - 2 sum n_d
    x = I(LO.divide(xn, xd), HI.divide(xn, xd))
    xp = [I(1)]
    for _ in range(2 * B):
        xp.append(xp[-1] * x)
    Z = I(1)
    for _, dd in Dp:
        Z = Z + xp[kap[dd]] * I(2)
    cnt, undecided = {}, 0
    for _, dd in Dp:
        v = I(N) * xp[kap[dd]] / Z
        lo, hi = int(v.lo.to_integral_value(rounding=decimal.ROUND_FLOOR)), int(v.hi.to_integral_value(rounding=decimal.ROUND_FLOOR))
        undecided += lo != hi
        cnt[dd] = lo
    n0 = N - 2 * sum(cnt.values())
    L = sum(cnt[dd] * kap[dd] for _, dd in Dp)
    check("every count n_d = floor(N x^kappa(d) / Z(x)) is decided by the enclosure", undecided == 0)
    check("n_0 = N - 2 sum n_d >= 0 and L = sum n_d kappa(d) >= B", counts_ok(N, cnt) and L >= B, f"n_0 = {n0}, L = {L}")
    check("n_0 and L as stated", (n0, L) == (row["n0"], row["L"]), f"stated n_0 = {row['n0']}, L = {row['L']}")
    if d and name == "m2026":
        lines = open(os.path.join(d, "pruned_differences_17086_27995.tsv")).read().split("\n")
        rows = [l.split("\t") for l in lines[1:] if l.strip()]
        same = len(rows) == len(Dp)
        bad_wit, Mset, tsv_wit = 0, set(M), []
        for (r, dd), f in zip(Dp, rows):
            rr, fd, fa, fb, fk, fp, fn, fc = (int(v) for v in f)
            same &= (rr, fd, fk, fp, fn, fc) == (r, dd, kap[dd], dd % q, -dd % q, cnt[dd])
            bad_wit += not witness_ok(Mset, dd, fa, fb, kap)
            tsv_wit.append((dd, fa, fb))
        check("the package's pruned_differences table: r, d, kappa(d), both residues and n_d all equal ours", same,
              f"{len(rows)} rows")
        check("every witness pair (alpha_d, beta_d) in the table lies in M^2 with difference d and sum kappa(d)",
              bad_wit == 0)
    # log D_N from below, log S_N(x) from above, log(fac q^N) from above
    from collections import Counter
    mult = Counter(cnt.values())
    lnD = lnfact(N) - lnfact(n0)
    for n, m in mult.items():
        lnD = lnD - lnfact(n) * I(2 * m)
    P = I(0)
    for s in S:
        P = P + xp[s]
    lnS = I(N) * iln(P) - I(2 * L) * iln(x)
    den = iln(I(fac)) + I(N) * iln(I(q))
    theta_lo = LO.add(1, LO.divide(LO.subtract(lnD.lo, lnS.hi), den.hi))
    print(f"  log D_N >= {LO.plus(lnD.lo):.25f}" + (f" (certificate: {cert['log_D_N_lower_bound']})" if cert else ""))
    print(f"  log S_N(x) <= {HI.plus(lnS.hi):.25f}" + (f" (certificate: {cert['log_S_N_upper_bound']})" if cert else ""))
    print(f"  log({fac} q^N) <= {HI.plus(den.hi):.25f}" + (f" (certificate: {cert['log_2maxV_plus_1_upper_bound']})" if cert else ""))
    print(f"  theta >= {theta_lo:.25f} (stated, pure-integer bounds: {row['stated']})")
    check(f"theta > {row['claim']}, the printed bound", theta_lo > Decimal(row["claim"]))
    if controls:
        print("Red controls (the same check functions, on altered data):")
        two = [(r, dd) for r, dd in Dp if q - r <= B and kap[q - r] is not None and kap[r] is not None]
        if two:
            r, dd = two[0]
            refused = not distinct_residues([e for _, e in Dp] + [q - r if dd == r else r], q)
            print(f"  R1 both members of the pair {{{r}, {q - r}}} kept in the alphabet: "
                  f"{'refused by the residue check, as expected' if refused else 'ACCEPTED (control failed)'}")
            ok &= refused
        else:
            print("  R1 (not applicable: in the carry-free base every pair has a single member in M - M)")
        if d and name == "m2026":
            dd, a, b = tsv_wit[len(tsv_wit) // 2]
            refused = not witness_ok(set(M), dd, a + 1, b + 1, kap)
            print(f"  R2 one witness of the table shifted to ({a + 1}, {b + 1}) for d = {dd} (same difference, cost "
                  f"kappa + 2): {'refused by the witness check, as expected' if refused else 'ACCEPTED (control failed)'}")
            ok &= refused
        bad = {k: v + v // 50 for k, v in cnt.items()}
        refused = not counts_ok(N, bad)
        print(f"  R3 every count raised by 2% (n_0 = {N - 2 * sum(bad.values())}): "
              f"{'refused by the count check, as expected' if refused else 'ACCEPTED (control failed)'}")
        ok &= refused
        refused = not (theta_lo > Decimal(row["claim"]) + Decimal("0.00001"))
        print(f"  R4 a claim 1e-5 above the printed one is not implied by our theta: "
              f"{'refused, as expected' if refused else 'ACCEPTED (control failed)'}")
        ok &= refused
    print(f"RESULT: {name} {'PASS' if ok else 'FAIL'}")
    print(f"run time {time.time() - t0:.0f} s")
    return ok


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    workers = next((int(a.split("=")[1]) for a in sys.argv[1:] if a.startswith("--workers=")),
                   min(8, os.cpu_count() or 1))
    if not args or args[0] not in ("syl2", "syl3", *M2026):
        print(__doc__)
        sys.exit(2)
    root = args[1] if len(args) > 1 else os.path.join(HERE, "certs")
    if args[0] in M2026:
        sys.exit(0 if check_m2026(root, args[0], "--controls" in sys.argv) else 1)
    build()
    if args[0] == "syl3":
        parts = next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--parts=")), None)
        sys.exit(0 if check_syl3(root, "--controls" in sys.argv, workers, parts) else 1)
    sys.exit(0 if check_syl2(root, "--controls" in sys.argv, workers) else 1)


if __name__ == "__main__":
    main()
