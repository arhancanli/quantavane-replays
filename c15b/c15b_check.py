#!/usr/bin/env python3
"""Rigorous replay of the dual matrix multiplication exponent bounds proposed for the 15b page in PR #337:
[OAI2026] alpha > 0.465 (the seven-step completion script (12.8)) and [Var2026] alpha > 0.46652506 (a 15-step script,
ArseniyVarlamov/dual-alpha-witness @ d2868fa), from Section 12.1 of the preprint and the certificate data alone.
Varlamov's verify.py is not read or run.

Section 12.1 of [OAI2026] (Lemma 12.1, (12.6), (12.7)): base two, distinguished side X and its colour B. A state is
(rho, b, h_A, h_B, h_C), initially (1/2, 1, 0, 0, 0). A completion of size m >= 2 centred on B gives
    rho' = rho^m,  h_B' = m h_B,  h_i' = m ((rho - rho')/(1 - rho') h_B + (1 - rho)/(1 - rho') h_i)   (i = A, C);
centred on i in {A, C}, with j the other inactive colour,
    rho' = 1 - (1 - rho)^m,  h_B' = m (rho/rho' h_B + (rho' - rho)/rho' h_i),  h_i' = m h_i,
    h_j' = log((e^h_i + e^h_j)^m - e^(m h_i));
and b' = m b. Terminating with B and the larger inactive colour gives D = (rho h_B + (1 - rho) max(h_A, h_C))/b and
alpha >= min{1, 2D/log 2}. With delta = |h_A - h_C| this reduces to (12.9): a step centred on B keeps D and multiplies
delta by m(1 - rho)/(1 - rho^m); a step centred on the larger inactive colour sets z = log((1 + e^-delta)^m - 1),
delta' = |z| and D' = D + (1 - rho') max{0, z}/(m b).

We check, in interval arithmetic of our own (fixed point with 448 fractional bits, outward rounding, exp and log by
series with explicit remainders, log 2 = 2 atanh(1/3)):
  1. the certificate's fields: the script, its base-factor count b = product of the sizes, and the claimed rational;
  2. D by the reduced recurrence (12.9), and the strict inequality 2D/log 2 > claim, by exact rational comparison of
     the interval ends (2 D_lo > claim log2_hi);
  3. D again by the full three-entropy recurrence of Lemma 12.1 (centres chosen by comparing h_A and h_C, ties to A),
     and that the two enclosures of D overlap;
  4. for the paper's script: its centres B, A, C, B, C, B, A, the certified information listed after (12.9) (rho
     exactly, as it is a dyadic rational), and the rational chain (12.10) ending at 1548637/3327360 > 0.465.
Red controls go through the same check function: the 15-step bound raised by 10^-8, the 15-step script without its
last step, a plus step re-centred on the smaller inactive colour, and the paper's script against the 15-step bound.

Also: python3 -I c15b_check.py longer   checks the longer scripts in longer_scripts.json (our own data: a beam search
over sizes 2 to 8 and every choice of centre found nothing above 0.4665256234; these are its best at 15 and 42 steps).

Usage: python3 -I c15b_check.py var2026 [CERTS] [--controls]   Python 3.9+, standard library only.
"""
import hashlib
import json
import os
import sys
import time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
VAR2026 = dict(rel="var2026/certificate.json", sha="8a81020fd6ac36a64b5e8b03e2dbf8c2860160a7e019e9c31bfe4dc50d5db345",
               claim=Fraction(23326253, 50000000), base_factors=147456, printed="0.466525061423431")
PAPER = dict(script=[-2, 2, 2, -3, 2, -2, 3], claim=Fraction(465, 1000), printed="0.466101554898",
             centres="BACBCBA")

# ---------- interval arithmetic: an interval is (lo, hi), integers with lo/2^PREC <= x <= hi/2^PREC ----------
PREC = 448
ONE = 1 << PREC
ONE_IV = (ONE, ONE)


def cst(q):
    q = Fraction(q)
    return ((q.numerator << PREC) // q.denominator, -((-q.numerator << PREC) // q.denominator))


def add(a, b):
    return (a[0] + b[0], a[1] + b[1])


def sub(a, b):
    return (a[0] - b[1], a[1] - b[0])


def neg(a):
    return (-a[1], -a[0])


def mul(a, b):
    ps = [x * y for x in a for y in b]
    return (min(ps) >> PREC, -((-max(ps)) >> PREC))


def div(a, b):
    assert b[0] > 0, "division by an interval that is not positive"
    return (min((x << PREC) // y for x in a for y in b), max(-((-x << PREC) // y) for x in a for y in b))


def powi(a, m):
    assert a[0] >= 0
    lo = hi = ONE
    for _ in range(m):
        lo, hi = (lo * a[0]) >> PREC, -((-(hi * a[1])) >> PREC)
    return (lo, hi)


def max0(a):
    return (max(0, a[0]), max(0, a[1]))


def iabs(a):
    if a[0] >= 0:
        return a
    if a[1] <= 0:
        return neg(a)
    return (0, max(-a[0], a[1]))


def order(a, b):
    """+1 if a > b, -1 if a < b, 0 if both are the same exact number; None if the enclosures overlap."""
    if a == b and a[0] == a[1]:
        return 0
    if a[0] > b[1]:
        return 1
    if a[1] < b[0]:
        return -1
    return None


NEXP, NATANH = 48, 150


def fmul(x, y, up):
    """x y / 2^PREC rounded down, or up."""
    return -((-(x * y)) >> PREC) if up else (x * y) >> PREC


def idiv(x, n, up):
    return -((-x) // n) if up else x // n


def _exp_nonneg(t, up):
    """Bound e^(t/2^PREC) for an integer t >= 0, from below or (up) from above."""
    k = max(0, t.bit_length() - (PREC - 12))                 # y = t/2^k < 2^-12 as a value
    y = idiv(t, 1 << k, up)
    s = term = ONE
    for n in range(1, NEXP):
        term = idiv(fmul(term, y, up), n, up)
        s += term
    if up:
        s += 1                                                 # tail sum_{n >= 48} y^n/n! < 2^-576 < one unit
    for _ in range(k):
        s = fmul(s, s, up)
    return s


def exp_iv(a):
    def bound(x, up):
        if x >= 0:
            return _exp_nonneg(x, up)
        e = _exp_nonneg(-x, not up)                            # e^x = 1/e^-x
        return -((-ONE << PREC) // e) if up else (ONE << PREC) // e
    return (bound(a[0], False), bound(a[1], True))


def _atanh(num, den, up):
    """Bound atanh(s) for a rational 0 <= s = num/den <= 1/3 by its series sum s^(2n+1)/(2n+1)."""
    s = idiv(num << PREC, den, up)
    s2 = fmul(s, s, up)
    total, p = 0, s
    for n in range(NATANH):
        total += idiv(p, 2 * n + 1, up)
        p = fmul(p, s2, up)
    if up:
        total += 1                                             # tail < 3^-301 (9/8)/301 < 2^-477 < one unit
    return total


LOG2 = (2 * _atanh(1, 3, False), 2 * _atanh(1, 3, True))


def _log_pos(v, up):
    assert v > 0, "log of a number that is not positive"
    L = v.bit_length() - 1                                     # v = 2^L w with 1 <= w < 2; s = (w-1)/(w+1) < 1/3
    e = L - PREC
    t = 2 * _atanh(v - (1 << L), v + (1 << L), up)
    if e >= 0:
        return e * LOG2[1 if up else 0] + t
    return e * LOG2[0 if up else 1] + t


def log_iv(a):
    return (_log_pos(a[0], False), _log_pos(a[1], True))


def digits(a, n=40):
    """The common decimal prefix of an enclosure of a number in (0, 1), at most n digits."""
    lo, hi = Fraction(a[0], ONE), Fraction(a[1], ONE)
    s1, s2 = str(lo.numerator * 10 ** (n + 2) // lo.denominator), str(hi.numerator * 10 ** (n + 2) // hi.denominator)
    k = 0
    while k < len(s1) and k < len(s2) and s1[k] == s2[k]:
        k += 1
    return "0." + s1.zfill(n + 2)[:min(max(k, 1), n)]


# ---------- the recurrences ----------

def steps_of(script):
    """Certificate signs: negative = centre B, positive = the inactive colour with larger h (ties A); 'S' = smaller."""
    out = []
    for s in script:
        if isinstance(s, str):
            out.append((s[0], int(s[1:])))
        else:
            out.append(("B" if s < 0 else "L", abs(s)))
    return out


def reduced(steps):
    """(12.9). Returns D, the final rho, b, and a trace of (b, rho, delta, D, z) after each step (z = None at B)."""
    rho, b, D, dl = cst(Fraction(1, 2)), 1, (0, 0), (0, 0)
    trace = []
    for kind, m in steps:
        z = None
        if kind == "B":
            r2 = powi(rho, m)
            den = p = ONE_IV
            for _ in range(m - 1):                             # m(1 - rho)/(1 - rho^m) = m/(1 + rho + ... + rho^(m-1))
                p = mul(p, rho)
                den = add(den, p)
            dl = mul(div(cst(m), den), dl)
        else:
            r2 = sub(ONE_IV, powi(sub(ONE_IV, rho), m))
            u = exp_iv(neg(dl))
            if kind == "L":
                z = log_iv(sub(powi(add(ONE_IV, u), m), ONE_IV))
                D = add(D, mul(div(sub(ONE_IV, r2), cst(m * b)), max0(z)))
                dl = iabs(z)
            else:                                              # centred on the smaller inactive entropy
                z = log_iv(sub(powi(add(ONE_IV, u), m), powi(u, m)))
                inc = sub(mul(sub(ONE_IV, r2), z), mul(mul(cst(m), sub(r2, rho)), dl))
                D = add(D, div(inc, cst(m * b)))
                dl = add(mul(cst(m), dl), z)
        rho, b = r2, b * m
        trace.append((b, rho, dl, D, z))
    return D, rho, b, trace


def full(steps):
    """Lemma 12.1 with the three entropies. Returns D and the centres used."""
    rho, b = cst(Fraction(1, 2)), 1
    h = {"A": (0, 0), "B": (0, 0), "C": (0, 0)}
    centres = ""
    for kind, m in steps:
        if kind == "B":
            r2 = powi(rho, m)
            q = sub(ONE_IV, r2)
            wB, wi = div(sub(rho, r2), q), div(sub(ONE_IV, rho), q)
            for i in "AC":
                h[i] = mul(cst(m), add(mul(wB, h["B"]), mul(wi, h[i])))
            h["B"] = mul(cst(m), h["B"])
            centres += "B"
        else:
            o = order(h["A"], h["C"])
            assert o is not None, "the order of h_A and h_C is not determined"
            big, small = ("A", "C") if o >= 0 else ("C", "A")
            i, j = (big, small) if kind == "L" else (small, big)
            r2 = sub(ONE_IV, powi(sub(ONE_IV, rho), m))
            hB = mul(cst(m), add(mul(div(rho, r2), h["B"]), mul(div(sub(r2, rho), r2), h[i])))
            if i == big:                                       # log((e^h_i + e^h_j)^m - e^(m h_i)) around max = h_i
                e = exp_iv(neg(sub(h[i], h[j])))
                hj = add(mul(cst(m), h[i]), log_iv(sub(powi(add(ONE_IV, e), m), ONE_IV)))
            else:                                              # around max = h_j
                e = exp_iv(neg(sub(h[j], h[i])))
                hj = add(mul(cst(m), h[j]), log_iv(sub(powi(add(ONE_IV, e), m), powi(e, m))))
            h["B"], h[i], h[j] = hB, mul(cst(m), h[i]), hj
            centres += i
        rho, b = r2, b * m
    o = order(h["A"], h["C"])
    assert o is not None, "the order of h_A and h_C is not determined at termination"
    hmax = h["A"] if o >= 0 else h["C"]
    D = div(add(mul(rho, h["B"]), mul(sub(ONE_IV, rho), hmax)), cst(b))
    return D, centres


def alpha_iv(D):
    return div(mul(cst(2), D), LOG2)


def check_script(script, claim, label, base_factors=None, verbose=True):
    """The check. Returns (ok, facts)."""
    say = print if verbose else (lambda *a, **k: None)
    steps = steps_of(script)
    ok, facts = True, {}
    if not all(m >= 2 for _, m in steps):
        say("  FAIL: a completion of size below 2")
        return False, facts
    b = 1
    for _, m in steps:
        b *= m
    facts["b"] = b
    if base_factors is not None and b != base_factors:
        say(f"  FAIL: the sizes multiply to b = {b}, the certificate says {base_factors}")
        ok = False
    D1, _, _, _ = reduced(steps)
    D2, centres = full(steps)
    a1, a2 = alpha_iv(D1), alpha_iv(D2)
    facts.update(alpha=a1, centres=centres)
    width = Fraction(a1[1] - a1[0], ONE)
    say(f"  b = {b}; centres {centres}")
    say(f"  2D/log 2 by (12.9):     {digits(a1, 60)}...  (enclosure width {float(width):.1e})")
    say(f"  2D/log 2 by Lemma 12.1: {digits(a2, 60)}...")
    if not (D1[0] <= D2[1] and D2[0] <= D1[1]):
        say("  FAIL: the two recurrences disagree")
        ok = False
    strict = Fraction(2 * D1[0], ONE) > claim * Fraction(LOG2[1], ONE)
    say(f"  2D/log 2 > {claim.numerator}/{claim.denominator} = {float(claim)}: {'yes' if strict else 'NO'}")
    ok = ok and strict
    say(f"  {label}: {'OK' if ok else 'FAIL'}")
    return ok, facts


def check_paper_text():
    """The certified information after (12.9) and the chain (12.10), for the script (12.8)."""
    steps = steps_of(PAPER["script"])
    _, _, _, tr = reduced(steps)
    fx = lambda a: (Fraction(a[0], ONE), Fraction(a[1], ONE))
    inside = lambda a, lo, hi: Fraction(lo) < fx(a)[0] and fx(a)[1] < Fraction(hi)
    meets = lambda a, c: not (fx(a)[1] < fx(c)[0] or fx(c)[1] < fx(a)[0])
    # rho exactly: it is a dyadic rational
    r = [Fraction(1, 2)]
    for kind, m in steps:
        r.append(r[-1] ** m if kind == "B" else 1 - (1 - r[-1]) ** m)
    by_b = {t[0]: t for t in tr}
    log3 = log_iv(cst(3))
    rows = [
        ("rho_8 = 1 - (3/4)^4", r[3] == 1 - Fraction(3, 4) ** 4),
        ("delta_8 = log(9/7)", meets(by_b[8][2], log_iv(cst(Fraction(9, 7))))),
        ("D_8 = (0.5625/4) log 3", meets(by_b[8][3], mul(cst(Fraction(5625, 40000)), log3))),
        ("plus steps at b = 4, 8: z = log 3, log(7/9)",
         meets(by_b[4][4], log3) and meets(by_b[8][4], log_iv(cst(Fraction(7, 9))))),
        ("rho_24 = rho_8^3, 0.3194 < rho_24 < 0.3195",
         r[4] == r[3] ** 3 and Fraction(3194, 10000) < r[4] < Fraction(3195, 10000)),
        ("0.348 < delta_24 < 0.354", inside(by_b[24][2], Fraction(348, 1000), Fraction(354, 1000))),
        ("1 - rho_48 = (1 - rho_24)^2 > 0.463, 0.536 < rho_48 < 0.537",
         1 - r[5] == (1 - r[4]) ** 2 and 1 - r[5] > Fraction(463, 1000)
         and Fraction(536, 1000) < r[5] < Fraction(537, 1000)),
        ("0.635 < z_48 < 0.651", inside(by_b[48][4], Fraction(635, 1000), Fraction(651, 1000))),
        ("rho_96 = rho_48^2, delta_96 = 2 delta_48/(1 + rho_48) < 0.86",
         r[6] == r[5] ** 2 and meets(by_b[96][2], div(mul(cst(2), by_b[48][2]), cst(1 + r[5])))
         and fx(by_b[96][2])[1] < Fraction(86, 100)),
        ("1 - rho_288 = (1 - rho_48^2)^3 > 0.360, z_288 > 0.628",
         1 - r[7] == (1 - r[5] ** 2) ** 3 and 1 - r[7] > Fraction(360, 1000)
         and fx(by_b[288][4])[0] > Fraction(628, 1000)),
        ("log 2 < 0.6932 and log 3 > 1.098",
         Fraction(LOG2[1], ONE) < Fraction(6932, 10000) and fx(log3)[0] > Fraction(1098, 1000)),
    ]
    chain = Fraction(2) / Fraction(6932, 10000) * (Fraction(5625, 10000) * Fraction(1098, 1000) / 4
                                                    + Fraction(463, 1000) * Fraction(635, 1000) / 48
                                                    + Fraction(360, 1000) * Fraction(628, 1000) / 288)
    rows.append(("(12.10): the chain equals 1548637/3327360 and exceeds 0.465",
                 chain == Fraction(1548637, 3327360) and chain > Fraction(465, 1000)))
    ok = True
    for text, good in rows:
        print(f"    {'ok  ' if good else 'FAIL'} {text}")
        ok = ok and good
    return ok


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    controls = "--controls" in sys.argv
    if args and args[0] == "longer":
        ok = True
        for entry in json.load(open(os.path.join(HERE, "longer_scripts.json"))):
            print(entry["name"])
            ok = check_script(entry["script"], Fraction(entry["bound"]), f"alpha > {entry['bound']}",
                              entry["base_factors"])[0] and ok
        print(f"RESULT: {'OK' if ok else 'FAIL'}")
        return 0 if ok else 1
    if not args or args[0] != "var2026":
        print(__doc__)
        return 2
    root = args[1] if len(args) > 1 else os.path.join(HERE, "certs")
    t0 = time.time()
    path = os.path.join(root, VAR2026["rel"])
    data = open(path, "rb").read()
    sha = hashlib.sha256(data).hexdigest()
    print(f"certificate {VAR2026['rel']}: sha256 {sha} {'(pinned)' if sha == VAR2026['sha'] else '(NOT the pinned file)'}")
    if sha != VAR2026["sha"]:
        return 1
    cert = json.loads(data)
    claim = Fraction(cert["bound"])
    print(f"script {cert['script']}, base factors {cert['base_factors']}, bound {cert['bound']}")
    ok = claim == VAR2026["claim"]
    print(f"interval arithmetic: {PREC} fractional bits; log 2 in [{digits(LOG2, 50)}...]")

    print("\n[Var2026] the 15-step script")
    good, facts = check_script(cert["script"], claim, "alpha > 23326253/50000000", cert["base_factors"])
    ok = ok and good
    starts = digits(facts["alpha"], 60).startswith(VAR2026["printed"])
    print(f"  agrees with the printed 2D/log 2 = {VAR2026['printed']}...: {'yes' if starts else 'NO'}")
    ok = ok and starts

    print("\n[OAI2026] the paper's seven-step script (12.8)")
    good, facts = check_script(PAPER["script"], PAPER["claim"], "alpha > 0.465")
    ok = ok and good
    starts = digits(facts["alpha"], 60).startswith(PAPER["printed"])
    cen = facts["centres"] == PAPER["centres"]
    print(f"  agrees with the printed 0.466101554898 (README of [Var2026] and the 15b row): {'yes' if starts else 'NO'}")
    print(f"  centres B, A, C, B, C, B, A as stated after (12.8): {'yes' if cen else 'NO'}")
    print("  the certified information after (12.9) and the chain (12.10):")
    good = check_paper_text()
    ok = ok and starts and cen and good

    if controls:
        print("\nRed controls (each must FAIL):")
        s15 = cert["script"]
        tests = [
            ("C1 the 15-step bound raised by 10^-8", s15, claim + Fraction(1, 10 ** 8), cert["base_factors"]),
            ("C2 the 15-step script without its last step", s15[:-1], claim, None),
            ("C3 its fifth step re-centred on the smaller inactive colour",
             s15[:4] + ["S2"] + s15[5:], claim, None),
            ("C4 the paper's seven-step script against the 15-step bound", PAPER["script"], claim, None),
        ]
        for text, sc, cl, bf in tests:
            good, facts = check_script(sc, cl, text, bf, verbose=False)
            print(f"  {text}: 2D/log 2 = {digits(facts['alpha'], 15)}...: {'passed (WRONG)' if good else 'fails, as it must'}")
            ok = ok and not good
    print(f"\nRESULT: {'OK' if ok else 'FAIL'}: alpha > 23326253/50000000 = 0.46652506 ([Var2026]) and "
          f"alpha > 0.465 ([OAI2026]), given Section 12.1 of [OAI2026]")
    print(f"run time {time.time() - t0:.1f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
