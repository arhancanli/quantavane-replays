#!/usr/bin/env python3
"""Exact replay of the two explicit lower-bound rows of the 37a page (teorth/optimizationproblems, main 87b4ebb):

  [HKP2011]            C_37 >= log_3 6 = 1.63093...  (Kushilevitz's function h on 6 bits: s(h) = 6, deg(h) = 3)
  [BdW2002], [T2013]   C_37 >= log_2 3 = 1.58496...  (not-all-equal g on 3 bits: s(g) = 3, deg(g) = 2)

C_37 is the supremum of log s(f) / log deg(f) over Boolean f with deg(f) > 1 (the page's definition), so one function
gives the bound.  Written from the page alone; the two polynomials are read from the page's text as data.

For each function: the printed polynomial is evaluated at every point of the cube (it must be 0 or 1 there); the
multilinear coefficients are recomputed from the truth table by the Moebius transform and must equal the printed
ones (the multilinear representation is unique), which gives the degree; the sensitivity is computed at every point.
The exponent is compared with the printed decimals (written with an approximately-equal sign, so rounded to nearest)
by exact integer powers; the bound itself is the exact log s / log deg.  Anchor: the composition f(f, ..., f) on
disjoint blocks (36 bits for h, 9 bits for g), expanded exactly as a polynomial: its degree is deg(f)^2 and its
sensitivity at 0 is s(f)^2, the same exponent.  Integers only.
usage: python3 -I kushilevitz_check.py CERTS_ROOT [--controls]   (CERTS_ROOT holds page/37a.md, from ./fetch_page.sh)
"""
import itertools
import random
import re
import sys
from math import log

ROOT = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("--") else "certs"
CONTROLS = "--controls" in sys.argv


def need(cond, msg):
    if not cond:
        raise AssertionError(msg)


def parse_page(path):
    """The two polynomials as {frozenset(variable indices): integer coefficient}, from the page's text."""
    text = open(path, encoding="utf-8").read()
    m = re.search(r"h\(z_1,\\dots,z_6\)=\\sum_\{i=1\}\^6 z_i - \\sum_\{1\\le i<j\\le 6\} z_i z_j \+ \(([^)]*)\)", text)
    need(m is not None, "the Kushilevitz polynomial was not found on the page")
    h = {frozenset([i]): 1 for i in range(1, 7)}
    for i, j in itertools.combinations(range(1, 7), 2):
        h[frozenset([i, j])] = -1
    triples = [t.strip() for t in m.group(1).split("+")]
    for t in triples:
        idx = [int(k) for k in re.findall(r"z_(\d)", t)]
        need(len(idx) == 3 and t.replace(" ", "") == "".join(f"z_{k}" for k in idx), f"unexpected term {t!r}")
        key = frozenset(idx)
        need(len(key) == 3 and key not in h, f"repeated or degenerate term {t!r}")
        h[key] = 1
    need(len(triples) == 10, "expected ten cubic terms")
    m = re.search(r"g\(x_1,x_2,x_3\)=([^,$]*),\$", text)
    need(m is not None, "the not-all-equal polynomial was not found on the page")
    g = {}
    for sign, mono in re.findall(r"([+-]?)((?:x_\d)+)", m.group(1)):
        key = frozenset(int(k) for k in re.findall(r"x_(\d)", mono))
        g[key] = g.get(key, 0) + (-1 if sign == "-" else 1)
    return h, g, len(triples)


def evaluate(poly, x):
    """x: dict variable -> 0/1."""
    return sum(c for s, c in poly.items() if all(x[v] for v in s))


def truth_table(poly, n):
    return [evaluate(poly, {v + 1: (bits >> v) & 1 for v in range(n)}) for bits in range(1 << n)]


def moebius(table, n):
    """Multilinear coefficients over {0,1}^n: a_S = sum_{T subset S} (-1)^{|S|-|T|} f(T)."""
    a = list(table)
    for v in range(n):
        bit = 1 << v
        for m in range(1 << n):
            if m & bit:
                a[m] -= a[m ^ bit]
    return {frozenset(v + 1 for v in range(n) if (m >> v) & 1): c for m, c in enumerate(a) if c}


def sensitivities(table, n):
    return [sum(1 for v in range(n) if table[x] != table[x ^ (1 << v)]) for x in range(1 << n)]


def rounds_to(s, d, num, den):
    """The printed decimal num/den is log s / log d rounded to nearest (the page writes it with an approximately-equal
    sign): (num - 1/2)/den <= log s / log d < (num + 1/2)/den, i.e. d^(2 num - 1) <= s^(2 den) < d^(2 num + 1)."""
    p = s ** (2 * den)
    return d ** (2 * num - 1) <= p < d ** (2 * num + 1)


def compose(poly, n):
    """The polynomial of f(f(block 1), ..., f(block n)) on n*n variables (block b uses variables n*(b-1)+1 .. n*b)."""
    inner = []
    for b in range(n):
        inner.append({frozenset(v + n * b for v in s): c for s, c in poly.items()})
    out = {}
    for s, c in poly.items():
        acc = {frozenset(): c}
        for b in sorted(s):
            nxt = {}
            for s1, c1 in acc.items():
                for s2, c2 in inner[b - 1].items():
                    k = s1 | s2           # disjoint blocks, so the product stays multilinear
                    nxt[k] = nxt.get(k, 0) + c1 * c2
            acc = nxt
        for k, c1 in acc.items():
            out[k] = out.get(k, 0) + c1
    return {k: c for k, c in out.items() if c}


def composed_value(poly, n, x):
    """f(f(x_block1), ..., f(x_blockn)) as functions, x a list of n*n bits."""
    outer = {b + 1: evaluate(poly, {v + 1: x[n * b + v] for v in range(n)}) for b in range(n)}
    return evaluate(poly, outer)


def check(name, poly, n, s_claim, d_claim, printed, quiet=False):
    """Return (s, deg, ok)."""
    table = truth_table(poly, n)
    boolean = all(v in (0, 1) for v in table)
    if not boolean:
        bad = next(x for x in range(1 << n) if table[x] not in (0, 1))
        if not quiet:
            print(f"{name}: NOT Boolean: value {table[bad]} at x = {format(bad, f'0{n}b')[::-1]} (x_1 first)")
        return None, None, False
    coeffs = moebius(table, n)
    need(coeffs == {k: c for k, c in poly.items() if c}, f"{name}: Moebius coefficients differ from the printed ones")
    deg = max(len(k) for k in coeffs)
    sens = sensitivities(table, n)
    s = max(sens)
    at = [format(x, f"0{n}b")[::-1] for x in range(1 << n) if sens[x] == s]
    ok = (s == s_claim and deg == d_claim)
    num, den = printed
    ge = rounds_to(s, deg, num, den)
    if not quiet:
        ones = sum(table)
        print(f"{name}: Boolean on all {1 << n} points ({ones} ones); the Moebius coefficients equal the printed "
              f"polynomial ({len(coeffs)} monomials); deg = {deg}; s = {s}, attained at {len(at)} point(s), "
              f"including {'0' * n}: {'0' * n in at}")
        print(f"  log s / log deg = {log(s) / log(deg):.15f}; printed {num / den:.{len(str(den)) - 1}f}: "
              f"{'the exact value rounded to nearest, as the page writes it (with an approximately-equal sign)' if ge else 'NOT the exact value rounded'} (exact integer powers)")
    return s, deg, ok and ge


def anchor(name, poly, n, s1, d1):
    comp = compose(poly, n)
    deg = max(len(k) for k in comp)
    top = sum(1 for k in comp if len(k) == deg)
    zero = [0] * (n * n)
    f0 = composed_value(poly, n, zero)
    s0 = 0
    for v in range(n * n):
        y = list(zero)
        y[v] = 1
        s0 += composed_value(poly, n, y) != f0
    # the expanded polynomial agrees with function composition at random points
    rng = random.Random(37)
    for _ in range(200):
        x = [rng.randint(0, 1) for _ in range(n * n)]
        need(evaluate(comp, {v + 1: x[v] for v in range(n * n)}) == composed_value(poly, n, x), "expansion mismatch")
    print(f"  anchor {name}({name}, ..., {name}) on {n * n} bits: {len(comp)} monomials, deg = {deg} = {d1}^2 "
          f"({top} top monomials), s at 0 = {s0} = {s1}^2; the same exponent")
    need(deg == d1 * d1 and s0 == s1 * s1, f"{name}: composition anchor failed")


def self_test():
    rng = random.Random(2026)
    for _ in range(40):
        n = rng.randint(1, 8)
        poly = {}
        for _ in range(rng.randint(1, 12)):
            k = frozenset(v + 1 for v in range(n) if rng.random() < 0.4)
            poly[k] = poly.get(k, 0) + rng.randint(-3, 3)
        poly = {k: c for k, c in poly.items() if c}
        need(moebius(truth_table(poly, n), n) == poly, "self-test: Moebius transform is not the inverse")
    print("self-test: the Moebius transform recovers 40 random integer multilinear polynomials (n <= 8)")


def main():
    h, g, nt = parse_page(f"{ROOT}/page/37a.md")
    print(f"page: page/37a.md (main 87b4ebb); Kushilevitz h with {len(h)} monomials ({nt} cubic), g with {len(g)}")
    self_test()
    ok = True
    s, d, good = check("h [HKP2011]", h, 6, 6, 3, (163093, 100000))
    ok &= good
    anchor("h", h, 6, s, d)
    print(f"RESULT [HKP2011]: s(h) = {s}, deg(h) = {d}, so C_37 >= log_3 6 = {log(6) / log(3):.12f}... (printed ~1.63093): "
          f"{'HOLDS' if good else 'FAILS'}")
    s, d, good = check("g [BdW2002]", g, 3, 3, 2, (158496, 100000))
    ok &= good
    anchor("g", g, 3, s, d)
    print(f"RESULT [BdW2002]/[T2013]: s(g) = {s}, deg(g) = {d}, so C_37 >= log_2 3 = {log(3) / log(2):.12f}... (printed "
          f"~1.58496): {'HOLDS' if good else 'FAILS'}")
    if CONTROLS:
        print("red controls (each must be refused or fail):")
        res = []
        h1 = dict(h)
        del h1[frozenset([1, 3, 4])]
        res.append(("h without z1z3z4", check("  h - z1z3z4", h1, 6, 6, 3, (163093, 100000))[2] is False))
        h2 = dict(h1)
        h2[frozenset([1, 2, 3])] = 1
        res.append(("h with z1z2z3 for z1z3z4", check("  h, z1z2z3 for z1z3z4", h2, 6, 6, 3, (163093, 100000))[2] is False))
        res.append(("h with claimed sensitivity 7", check("  h, claimed s = 7", h, 6, 7, 3, (163093, 100000), quiet=True)[2] is False))
        res.append(("h with claimed exponent 1.63094", check("  h, claimed 1.63094", h, 6, 6, 3, (163094, 100000), quiet=True)[2] is False))
        g1 = dict(g)
        g1[frozenset([1, 2])] = -2
        res.append(("g with -2 x1x2", check("  g, -2 x1x2", g1, 3, 3, 2, (158496, 100000))[2] is False))
        res.append(("g with claimed degree 1", check("  g, claimed deg = 1", g, 3, 3, 1, (158496, 100000), quiet=True)[2] is False))
        for label, refused in res:
            print(f"  {label}: {'refused, as expected' if refused else 'NOT REFUSED'}")
        ok &= all(r for _, r in res)
        print("controls all as expected" if all(r for _, r in res) else "CONTROLS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
