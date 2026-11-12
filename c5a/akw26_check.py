#!/usr/bin/env python3
"""Exact replay of the 5a row [Akw26] C_5a <= 0.9428095* (teorth/optimizationproblems; Erdos problem 30) through John
Akwei's own certificates, Zenodo 10.5281/zenodo.24138215 (v1.0, 11 Nov 2026), written from his report's definitions
and the deposit README's file format alone. No program of the deposit was read or run; the README and the report were
read as text.

Sources (SHA-256 pinned below and in ../SHA256SUMS_akw26; ./fetch_akw26.sh downloads them):
  J. Akwei, "A Data Science Analysis of Erdos Problem 30", Version 6 (26 Oct 2026; johnakwei/Science @ f357556):
  Appendix D (kernels and weights), Theorem 5.4 and the text after it, Table 10 (the printed values, transcribed in
  TABLE10 below); README_symmetric_certificates.md of the deposit, "Replaying the certificates from this README alone".

Definitions (Appendix D): n >= 1, m = 2n; for 1 <= l <= n the kernel p^(l) in R^m is uniform on the 2l central cells,
p^(l)_i = 1/(2l) for n - l <= i <= n + l - 1, with weights lambda_l = 4 l^3 / n^4 (l < n), lambda_n = (2n - 1)/n^2.
File format (README): every number is an integer over D = 10^12. _P.txt has m rows and n columns, column k is kernel k
(paired with line k of _lam.txt); _lam.txt has n lines; _W_end1.txt has Lm rows and n columns, the boundary weights
w_k[j] for 0 <= j < Lm, with w_k[j] = D for j >= Lm; the other end carries the same weights in mirror image.
Checks, in integers and Fractions:
  0. format: the stated shapes, ASCII decimal integers, the pinned SHA-256;
  1. kernels nonnegative, every column summing to D (mass 1), row i equal to row m - 1 - i (symmetric); weights
     nonnegative and summing to D (sum of lambda = 1);
  2. covering: T(q) = sum_k lambda_k sum_i p_k[i] w_k[q + i] >= D^3 for every q = 0..Lm, computed by two routes that
     share no code (prefix sums over the constant runs of each kernel; one Kronecker-substituted big-integer product
     per kernel) and required identical at every q, plus direct summation at sample shifts; the least slack
     T(q) - D^3 and where it occurs;
  3. a = m sum_k lambda_k sum_i p_k[i]^2 / D^3 and b = 1 + 2 (sum_k lambda_k sum_{j<Lm} w_k[j]^2 / (m D^3) - L), b > 0;
     sqrt(ab) floored to 30 digits;
  4. a, b and sqrt(ab) against Table 10 and the README (each printed to 12 places, rounded); ab < g^2 for the printed
     certified bound and for the printed ten-place bound; for the row's certificate ab <= (9428095/10^7)^2. Every
     comparison of sqrt(ab) is an exact cross-multiplication of integers.
Reported as provenance (not needed for the bound): each column is Appendix D's kernel rounded (column k is l = k + 1,
supported on its 2l central cells, value within 1/2 of D/(2l) except the central pair n - 1, n, which takes the
remainder so that the column sums to D); each weight is within 1/2 of its exact value except lambda_(n-1), which is D
minus the others; the cells no covering inequality reads (j < n - l) all hold the stated repair increment.
If 1-3 hold and b > 0, Hou and Zhao's Theorem 2.1 (arXiv 2607.01169) in the symmetric case gives
h(N) <= N^(1/2) + sqrt(ab) N^(1/4) + O(1) (V6, Theorem 5.4). Not checked: that theorem (mathematics, read).

Where Akwei's construction differs from Varlamov's (PR #335, c5a_check.py), followed here as Akwei defines it:
  kernels and weights: the same Appendix D family, rounded to integers over 10^12 (Varlamov: exact 1/(2l), lambda_l);
  boundary weights: listed explicitly, Lm x n integers over 10^12 (the boundary program re-solved for the rounded
    kernels, rounded, plus a uniform increment); the multiplier is not in the files. Varlamov's are generated from
    stored integers V_k on the even sites (window sums, floor, repair r, over 10^18). Both rest on the discrete renewal
    multiplier of V6 (after Theorem 5.4): Varlamov's V_k / 10^18 are half of it, and the multiplier implied by Akwei's
    weights matches it closely (--multiplier prints how closely, in floating point, as information only);
  scale: one denominator D = 10^12 for p, lambda and w, covering threshold D^3 (Varlamov: 2 n^4 D with D = 10^18);
  depth: L = 6 at m = 1024, 6,145 inequalities (Varlamov: L = 12, 12,289); repair 11 x 10^-12 (Varlamov: r = 43922
    x 10^-18).
Red controls (--controls, on the row's certificate m = 1024, L = 6; each must be refused):
  C1  one kernel integer one unit lower (p_512[0]); C1b every boundary integer one unit lower (the uniform repair one
      unit smaller); C2 one unit of weight moved from lambda_1 to lambda_512 (the sum stays D); C2b lambda_512 one unit
      lower; C3 a claim just below the computed sqrt(ab) (its 20-digit floor); C4 one unit moved from p_512[0] to
      p_512[1] (mass kept, symmetry broken). Also reported: no single boundary integer, lowered by one unit, breaks
      any covering inequality.
Usage: python3 -I akw26_check.py [CERTS_ROOT] [--controls] [--multiplier]   (default ./certs; files in CERTS_ROOT/akw26/)
Python 3.9+, standard library only. About a minute and 400 MB with both flags.
"""
import hashlib
import os
import re
import sys
import time
from fractions import Fraction
from math import isqrt

HERE = os.path.dirname(os.path.abspath(__file__))
D = 10 ** 12
D3 = D ** 3
SHA = {
    "sym_polya_m1024_L6_P.txt": "27f217fb0f285cd06d5998568a995667daebdefcc90ef92a7ff255c2a311c987",
    "sym_polya_m1024_L6_W_end1.txt": "73e04925e142b46225a194e085f43c9603dddc33a79478c4a85267a484dea27a",
    "sym_polya_m1024_L6_lam.txt": "10b00a21a91f54c31fb1acbee9a8b044c6f97ba28d44296d4753555fbe3fa135",
    "sym_polya_m128_L8_P.txt": "147a085ee07cee665f2020ee207ac977a7160b54431312e1369e9295bb0c9416",
    "sym_polya_m128_L8_W_end1.txt": "73b7fa219999abb8eebea7416b9b99361a9bfe46d975b8848f14350409d3097e",
    "sym_polya_m128_L8_lam.txt": "28408a48c4fe1a9d9dac2f11872d066a9848fa98dcbb8eaeb3a88906d8cdd153",
    "sym_polya_m256_L8_P.txt": "f28570caabb8a15f92cf16905d6cc2f079ec5df7aa1ae98d0325095e46b0a113",
    "sym_polya_m256_L8_W_end1.txt": "7c1c2356876f51863ba27aa717036d1d3eed61e549a39ebd98994bb74e2fea8e",
    "sym_polya_m256_L8_lam.txt": "c69b2cc73aee1a69dc2e4ee0c5a086ff2d8f4ef1527a61bd81de7f0c7871a5b2",
    "sym_polya_m512_L8_P.txt": "bf9a8853a4ac6ca2049c0697a4827b635ea5cee4411df97ebe1eb8485e423627",
    "sym_polya_m512_L8_W_end1.txt": "6668bd18476dd3d59c416afb5a6168ae988349e4bf382f65d6af046714426233",
    "sym_polya_m512_L8_lam.txt": "6ed4f279e5e8a0a26c3d7f553f724f8fcd66d5432426fccce4599b4883be2a00",
    "README_symmetric_certificates.md": "c77659e6ff63e5f109797ecdbc1d9841d9c57dc76ad41597dc67a5dd30babdc9",
}
# Table 10 of the V6 report (p. 25), as printed: m, L, kernels, covering inequalities, a, b, gamma, certified bound;
# then the ten-place bounds of its caption and the repair increments (x 10^-12) it states. The row's certificate first.
TABLE10 = [
    (1024, 6, 512, 6145, "1.333332061750", "0.666667938506", "0.942809491335", "0.942809493", "0.9428094915", 11),
    (128, 8, 64, 1025, "1.333251953159", "0.666748046890", "0.942837809902", "0.942837811", "0.9428378101", 1),
    (256, 8, 128, 2049, "1.333312988234", "0.666687011774", "0.942816234419", "0.942816236", "0.9428162346", 1),
    (512, 8, 256, 4097, "1.333328247092", "0.666671753640", "0.942810840342", "0.942810842", "0.9428108405", 23),
]
ROW = (9428095, 10 ** 7)       # the 5a row: C_5a <= 0.9428095


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_ints(path, nrows, ncols):
    """The file as a list of nrows rows of ncols ints, or (None, reason) if it is not in the stated format."""
    with open(path, "rb") as f:
        data = f.read()
    lines = data.split(b"\n")
    if lines and lines[-1].strip() == b"":
        lines.pop()
    if len(lines) != nrows:
        return None, f"{len(lines)} lines, not {nrows}"
    out = []
    for ln, raw in enumerate(lines):
        try:
            toks = raw.decode("ascii").split()
        except UnicodeDecodeError:
            return None, f"line {ln + 1} is not ASCII"
        if len(toks) != ncols:
            return None, f"line {ln + 1} has {len(toks)} fields, not {ncols}"
        for t in toks:
            if not (t.isdigit() or (t[0] == "-" and t[1:].isdigit())):
                return None, f"line {ln + 1}: {t!r} is not a decimal integer"
        out.append([int(t) for t in toks])
    return out, ""


def load(root, m, L):
    n, base = m // 2, f"sym_polya_m{m}_L{L}"
    hashes, data = {}, {}
    for suf, shape in (("_P.txt", (m, n)), ("_lam.txt", (n, 1)), ("_W_end1.txt", (L * m, n))):
        name = base + suf
        path = os.path.join(root, "akw26", name)
        hashes[name] = sha256(path)
        rows, why = read_ints(path, *shape)
        if rows is None:
            raise SystemExit(f"{name}: not in the README's format: {why}")
        data[suf] = rows
    P, lam, W = data["_P.txt"], [r[0] for r in data["_lam.txt"]], data["_W_end1.txt"]
    pc = [list(c) for c in zip(*P)]            # pc[k][i] = p_k[i]
    sym_rows = all(P[i] == P[m - 1 - i] for i in range(m))
    del P, data
    wc = [list(c) for c in zip(*W)]            # wc[k][j] = w_k[j], 0 <= j < Lm
    del W
    return hashes, pc, lam, wc, sym_rows


def runs_of(p):
    """Maximal runs of equal nonzero values: (start, end, value), end exclusive."""
    out, i, m = [], 0, len(p)
    while i < m:
        j = i
        while j + 1 < m and p[j + 1] == p[i]:
            j += 1
        if p[i]:
            out.append((i, j + 1, p[i]))
        i = j + 1
    return out


def covering_runs(pc, lam, wc, m, M):
    """Route 1: T(q) for q = 0..M by prefix sums of each boundary vector (extended by D) over the kernel's runs."""
    T = [0] * (M + 1)
    for p, lk, w in zip(pc, lam, wc):
        S, acc = [0], 0
        for x in w:
            acc += x
            S.append(acc)
        for _ in range(m):
            acc += D
            S.append(acc)
        for s0, e, v in runs_of(p):
            c = lk * v
            T = [t + c * (hi - lo) for t, hi, lo in zip(T, S[e:e + M + 1], S[s0:s0 + M + 1])]
    return T


def covering_kronecker(pc, lam, wc, m, M):
    """Route 2: T(q) as the coefficients of sum_k W_k(X) P_k(X) at X = 256^nb, where W_k packs w_k (extended by D) and
    P_k packs lambda_k p_k reversed. Needs every entry >= 0; nb is chosen so that no coefficient can carry."""
    top = max(D, max(max(w) for w in wc))
    bound = sum(lk * sum(p) for p, lk in zip(pc, lam)) * top
    nb = bound.bit_length() // 8 + 2
    pad = D.to_bytes(nb, "little") * m
    acc = 0
    for p, lk, w in zip(pc, lam, wc):
        A = int.from_bytes(b"".join(x.to_bytes(nb, "little") for x in w) + pad, "little")
        B = int.from_bytes(b"".join((lk * x).to_bytes(nb, "little") for x in reversed(p)), "little")
        acc += A * B
    raw = acc.to_bytes(nb * (M + 2 * m), "little")
    return [int.from_bytes(raw[(q + m - 1) * nb:(q + m) * nb], "little") for q in range(M + 1)]


def covering_direct(pc, lam, wc, m, M, q):
    """Route 3, one shift at a time: the literal double sum."""
    t = 0
    for p, lk, w in zip(pc, lam, wc):
        s = 0
        for i in range(m):
            if p[i]:
                s += p[i] * (w[q + i] if q + i < M else D)
        t += lk * s
    return t


def floor_dec(fr, d):
    v = fr.numerator * 10 ** d // fr.denominator
    s = str(v).rjust(d + 1, "0")
    return s[:-d] + "." + s[-d:]


def floor_sqrt_dec(fr, d):
    v = isqrt(fr.numerator * 10 ** (2 * d) // fr.denominator)
    s = str(v).rjust(d + 1, "0")
    return s[:-d] + "." + s[-d:]


def rounds_to(x, printed):
    """x lies within half a unit of the printed decimal's last place."""
    h = Fraction(1, 2 * 10 ** len(printed.split(".")[1]))
    return Fraction(printed) - h <= x < Fraction(printed) + h


def sqrt_rounds_to(x2, printed):
    """sqrt(x2) lies within half a unit of the printed decimal's last place: squares compared exactly."""
    h = Fraction(1, 2 * 10 ** len(printed.split(".")[1]))
    lo, hi = Fraction(printed) - h, Fraction(printed) + h
    return lo * lo <= x2 < hi * hi


def below(x2, bound):
    """x2 < bound^2 for a printed decimal bound, by cross-multiplication of integers."""
    b = Fraction(bound)
    return x2.numerator * b.denominator ** 2 < b.numerator ** 2 * x2.denominator


def core(pc, lam, wc, m, L, routes=True):
    """The certificate checks 1-3. Returns a dict of facts; with routes=False only route 1 is used (controls)."""
    n, M = m // 2, L * m
    f = {}
    f["p_nonneg"] = all(x >= 0 for p in pc for x in p)
    f["p_mass"] = all(sum(p) == D for p in pc)
    f["p_sym"] = all(p == p[::-1] for p in pc)
    f["lam_nonneg"] = all(x >= 0 for x in lam)
    f["lam_sum"] = sum(lam) == D
    f["w_min"], f["w_max"] = min(min(w) for w in wc), max(max(w) for w in wc)
    T = covering_runs(pc, lam, wc, m, M)
    f["T"] = T
    if routes:
        if f["p_nonneg"] and f["lam_nonneg"] and f["w_min"] >= 0:
            f["routes_agree"] = covering_kronecker(pc, lam, wc, m, M) == T
        else:
            f["routes_agree"] = False
    f["fails"] = [q for q in range(M + 1) if T[q] < D3]
    A2 = sum(lk * sum(x * x for x in p) for p, lk in zip(pc, lam))
    B2 = sum(lk * sum(x * x for x in w) for w, lk in zip(wc, lam))
    f["a"] = Fraction(m * A2, D3)
    f["b"] = 1 + 2 * (Fraction(B2, m * D3) - L)
    f["ab"] = f["a"] * f["b"]
    f["ok"] = (f["p_nonneg"] and f["p_mass"] and f["p_sym"] and f["lam_nonneg"] and f["lam_sum"] and not f["fails"]
               and f["b"] > 0 and f.get("routes_agree", True))
    return f


def provenance(pc, lam, wc, m, repair):
    """Is each column Appendix D's kernel, rounded; each weight its exact value, rounded; the unread cells = repair?"""
    n = m // 2
    supp = vals = True
    dev_p = Fraction(0)
    for k, p in enumerate(pc):
        l = k + 1
        lo, hi = n - l, n + l
        supp &= all(x == 0 for x in p[:lo]) and all(x == 0 for x in p[hi:]) and all(x > 0 for x in p[lo:hi])
        exact = Fraction(D, 2 * l)
        c = p[lo]
        rest = D - (2 * l - 2) * c
        expect = [0] * lo + [c] * (2 * l) + [0] * (m - hi)
        if rest % 2 == 0:
            expect[n - 1] = expect[n] = rest // 2
        vals &= abs(c - exact) <= Fraction(1, 2) and p == expect
        dev_p = max(dev_p, max(abs(x - exact) for x in p[lo:hi]))
    exact_l = [Fraction(4 * l ** 3 * D, n ** 4) for l in range(1, n)] + [Fraction((2 * n - 1) * D, n * n)]
    near = all(abs(lam[k] - exact_l[k]) <= Fraction(1, 2) for k in range(n) if k != n - 2)
    adj = lam[n - 2] == D - sum(lam[k] for k in range(n) if k != n - 2)
    dev_adj = lam[n - 2] - exact_l[n - 2]
    unread = set()
    for k, w in enumerate(wc):
        unread.update(w[:n - (k + 1)])
    return dict(supp=supp, vals=vals, dev_p=dev_p, lam_near=near, lam_adj=adj, dev_adj=dev_adj, unread=unread,
                repair_ok=unread == {repair})


def multiplier_note(pc, wc, m, L, repair):
    """Informational, in floating point (--multiplier; not part of the check). If w_k = (1/2) mu * p_k for one
    multiplier mu on the shifts, the widest column (p_n = 1/m on every cell, exact in all four files) gives
    mu_j - mu_(j-m) = 2m (w_n[j] - w_n[j-1]). Returns: the largest |mu| on odd shifts; the largest relative deviation of
    mu_(2k) from 4n^2/(2n-1) U_k, twice the discrete renewal multiplier of V6 (after Theorem 5.4), U the renewal
    sequence of f_j = 2/(2n-1) (1 <= j < n), f_n = 1/(2n-1); the largest |(1/2) mu * p_k - w_k| / D over every column
    and every read cell (j >= n - l)."""
    n, M = m // 2, L * m
    wn = [x - repair for x in wc[n - 1]]
    mu = [0.0] * M
    for j in range(M):
        mu[j] = 2 * m * (wn[j] - (wn[j - 1] if j else 0)) / D + (mu[j - m] if j >= m else 0.0)
    f = [0.0] + [2 / (2 * n - 1)] * (n - 1) + [1 / (2 * n - 1)]
    U = [1.0]
    for k in range(1, M // 2):
        U.append(sum(f[j] * U[k - j] for j in range(1, min(k, n) + 1)))
    dev = max(abs(mu[2 * k] - 4 * n * n / (2 * n - 1) * U[k]) / (4 * n * n / (2 * n - 1) * U[k]) for k in range(M // 2))
    odd = max(abs(mu[q]) for q in range(1, M, 2))
    S = [0.0]
    for x in mu:
        S.append(S[-1] + x)
    worst = 0.0
    for k, (p, w) in enumerate(zip(pc, wc)):
        l = k + 1
        runs = [(s0, e, v / D) for s0, e, v in runs_of(p)]
        for j in range(n - l, M):      # (1/2) sum_i p[i] mu_(j - i) over 0 <= j - i
            pred = 0.5 * sum(v * (S[max(0, j - s0 + 1)] - S[max(0, j - e + 1)]) for s0, e, v in runs)
            worst = max(worst, abs(pred - (w[j] - repair) / D))
    return odd, dev, worst


def readme_table(root):
    """The README's printed values: the 'Results' tables (both) and the ten-place sentence."""
    path = os.path.join(root, "akw26", "README_symmetric_certificates.md")
    text = open(path, encoding="utf-8").read()
    first = {}
    for g in re.finditer(r"^\| m = (\d+), L = (\d+), K = (\d+) \| (\d+) \| (\d+) \| ([0-9.]+) \| ([0-9.]+) \| "
                         r"([0-9.]+) \| γ < ([0-9.]+) \|$", text, re.M):
        m, L, K, ineq, rep = (int(x) for x in g.groups()[:5])
        first[(m, L)] = (K, ineq, g.group(6), g.group(7), g.group(8), g.group(9), rep)
    second = {}
    for g in re.finditer(r"^\| (\d+) \| (\d+) \| (\d+) \| (\d+) \| ([0-9.]+) \| ([0-9.]+) \| ([0-9.]+) \| "
                         r"γ < ([0-9.]+) \|$", text, re.M):
        m, L, K, ineq = (int(x) for x in g.groups()[:4])
        second[(m, L)] = (K, ineq, g.group(5), g.group(6), g.group(7), g.group(8))
    ten = re.search(r"To ten places the bounds are ([0-9.]+), ([0-9.]+), ([0-9.]+) and ([0-9.]+)\.", text)
    return first, second, (ten.groups() if ten else None)


def check_one(root, row, readme, verbose_controls):
    m, L, K, ineq, a_pr, b_pr, g_pr, cert_pr, ten_pr, repair = row
    n, M = m // 2, L * m
    t0 = time.time()
    name = f"sym_polya_m{m}_L{L}"
    hashes, pc, lam, wc, sym_rows = load(root, m, L)
    pinned = all(hashes[h] == SHA[h] for h in hashes)
    is_row = (m, L) == (1024, 6)
    print(f"Certificate {name}{' (the [Akw26] row of 5a)' if is_row else ''}: m = {m}, L = {L}, n = {n} kernels, "
          f"Lm = {M}, D = 10^12; sha256 of the three files as pinned: {pinned}")
    print(f"  files: _P {m} x {n}, _lam {n}, _W_end1 {M} x {n} decimal integers, as the README states: True")
    f = core(pc, lam, wc, m, L)
    print(f"  kernels: nonnegative {f['p_nonneg']}; every column sums to D (mass 1) {f['p_mass']}; "
          f"row i = row m - 1 - i (symmetric) {f['p_sym'] and sym_rows}")
    print(f"  mixing weights: nonnegative {f['lam_nonneg']}; sum = D (the lambdas sum to 1) {f['lam_sum']}")
    print(f"  boundary weights: min {f['w_min']}, max {f['w_max']} (= {floor_dec(Fraction(f['w_max'], D), 6)} D)")
    pv = provenance(pc, lam, wc, m, repair)
    print(f"  Appendix D family, rounded: column k supported on the 2l central cells of l = k + 1: {pv['supp']}; "
          f"value within 1/2 of D/(2l) except the central pair (cells n - 1, n), which takes the remainder: "
          f"{pv['vals']} (largest |p - D/(2l)| = {float(pv['dev_p']):.1f} units)")
    print(f"  weights within 1/2 of 4l^3/n^4 D and (2n - 1)/n^2 D except lambda_(n-1): {pv['lam_near']}; "
          f"lambda_(n-1) = D minus the others: {pv['lam_adj']} ({float(pv['dev_adj']):+.4f} units from its exact value)")
    print(f"  the {n * (n - 1) // 2} boundary cells no inequality reads (j < n - l) all hold {sorted(pv['unread'])}; "
          f"Table 10's repair increment {repair} x 10^-12: {pv['repair_ok']}")
    T = f["T"]
    sl = [t - D3 for t in T]
    least = min(sl)
    at = [q for q in range(M + 1) if sl[q] == least]
    pos = [(s, q) for q, s in enumerate(sl) if s > 0]
    s1, q1 = min(pos)
    inner = min((sl[q], q) for q in range(M - m + 1))
    sample = sorted({0, 1, 2, m - 1, q1, inner[1], at[0], M - m, M - 2, M - 1, M} | {k * (M // 8) for k in range(9)})
    direct = all(covering_direct(pc, lam, wc, m, M, q) == T[q] for q in sample)
    print(f"  covering q = 0..{M}: {M + 1 - len(f['fails'])} of {M + 1} hold (Table 10: {ineq}); both routes "
          f"identical at every q: {f['routes_agree']}; direct sums at {len(sample)} shifts equal: {direct}")
    print(f"    least slack T(q) - D^3: {'0 (the covering sum exactly 1)' if least == 0 else least} at q = {at}; "
          f"least positive slack {s1}/D^3 = {float(Fraction(s1, D3)):.6e} at q = {q1}; least slack with the window "
          f"inside the box (q <= {M - m}): {inner[0]}/D^3 = {float(Fraction(inner[0], D3)):.6e} at q = {inner[1]}")
    if f["fails"]:
        print(f"    FAILING at q = {f['fails'][:10]}")
    a, b, ab = f["a"], f["b"], f["ab"]
    g20 = floor_sqrt_dec(ab, 30)
    print(f"  a = {floor_dec(a, 25)}... = {a.numerator}/{a.denominator}")
    print(f"  b = {floor_dec(b, 25)}... (exact fraction, {b.numerator.bit_length()}/{b.denominator.bit_length()} "
          f"bits); b > 0: {b > 0}")
    print(f"  ab = {floor_dec(ab, 25)}...; sqrt(ab) = {g20}... (floor, 30 digits)")
    cmp_a, cmp_b, cmp_g = rounds_to(a, a_pr), rounds_to(b, b_pr), sqrt_rounds_to(ab, g_pr)
    print(f"  printed (Table 10): a {a_pr} {cmp_a}; b {b_pr} {cmp_b}; gamma {g_pr} {cmp_g} (each the value rounded "
          f"to 12 places); kernels {K} = n: {K == n}; inequalities {ineq} = Lm + 1: {ineq == M + 1}")
    first, second, ten = readme
    r1, r2 = first.get((m, L)), second.get((m, L))
    same_readme = (r1 == (K, ineq, a_pr, b_pr, g_pr, cert_pr, repair) and r2 == (K, ineq, a_pr, b_pr, g_pr, cert_pr)
                   and ten is not None and ten_pr in ten)
    print(f"  the README prints the same values (both tables, repair {repair}, ten-place bound {ten_pr}): {same_readme}")
    c_cert, c_ten = below(ab, cert_pr), below(ab, ten_pr)
    ten_least = Fraction(ten_pr) - Fraction(1, 10 ** 10)
    print(f"  ab < {cert_pr}^2 (the certified bound): {c_cert}; ab < {ten_pr}^2 (the ten-place bound): {c_ten}; "
          f"also ab < {floor_dec(ten_least, 10)}^2: {below(ab, floor_dec(ten_least, 10))}")
    valid = f["ok"] and pinned and direct and c_cert
    agree = (cmp_a and cmp_b and cmp_g and K == n and ineq == M + 1 and same_readme and c_ten and least == 0
             and pv["supp"] and pv["vals"] and pv["lam_near"] and pv["lam_adj"] and pv["repair_ok"])
    if is_row:
        rn, rd = ROW
        c_row = ab.numerator * rd * rd <= rn * rn * ab.denominator
        print(f"  the row: ab * 10^14 <= 9428095^2 (sqrt(ab) <= 0.9428095), integers cross-multiplied: {c_row}; "
              f"0.9428095 - sqrt(ab) = {float(Fraction(rn, rd) - Fraction(g20)):.6e}")
        valid &= c_row
    if "--multiplier" in sys.argv:
        odd, dev, worst = multiplier_note(pc, wc, m, L, repair)
        print(f"  multiplier (floating point, informational, not part of the check): the widest column implies one mu "
              f"with w_k = (1/2) mu * p_k; |mu| on odd shifts <= {odd:.1e}; on even shifts 2k, mu = 4n^2/(2n-1) U_k "
              f"(twice V6's discrete renewal multiplier) to a relative {dev:.1e}; every column rebuilt from this mu "
              f"agrees with the file to {worst:.1e} D")
    print(f"RESULT {name}: {'PASS' if valid else 'FAIL'}; agreement with Table 10, the README and Appendix D: "
          f"{agree} ({time.time() - t0:.0f} s)")
    ok = valid and agree
    if verbose_controls and is_row:
        ok &= controls(pc, lam, wc, m, L, f, sl)
    return ok, f


def controls(pc, lam, wc, m, L, f, sl):
    n, M = m // 2, L * m
    print("Red controls on sym_polya_m1024_L6 (altered copies, the same checks; each must be refused):")
    # how far a single boundary integer can move a covering sum: at q, the largest lambda_k p_k[i] with q + i < Lm
    top = [max(p) for p in pc]
    pref_max = []
    for p in pc:
        cur, run = 0, []
        for x in p:
            cur = max(cur, x)
            run.append(cur)
        pref_max.append(run)
    full = max(lk * t for lk, t in zip(lam, top))
    worst = None
    for q in range(M + 1):
        r = M - q
        e = full if r >= m else (max(lk * pm[r - 1] for lk, pm in zip(lam, pref_max)) if r > 0 else 0)
        if e and (worst is None or Fraction(sl[q], e) < worst[0]):
            worst = (Fraction(sl[q], e), q, e)
    print(f"  fact: no boundary integer lowered by one unit breaks an inequality: the least ratio of slack to the "
          f"largest single-entry effect is {float(worst[0]):.4f} >= 1 at q = {worst[1]}: {worst[0] >= 1}")
    res = []

    def run(label, pc2=None, lam2=None, wc2=None):
        g = core(pc2 or pc, lam2 or lam, wc2 or wc, m, L, routes=False)
        why = [w for w, bad in (("mass", not g["p_mass"]), ("symmetry", not g["p_sym"]),
                                ("sum of weights", not g["lam_sum"]), ("negative entry", not (g["p_nonneg"] and
                                                                                               g["lam_nonneg"])))
               if bad]
        fl = g["fails"]
        if fl:
            why.append(f"covering fails at {len(fl)} q ({fl[:6]}{' ...' if len(fl) > 6 else ''})")
        refused = not g["ok"]
        print(f"  {label}: {'; '.join(why) if why else 'nothing fails'}: "
              f"{'refused, as expected' if refused else 'ACCEPTED (control failed)'}")
        res.append(refused)

    p2 = [list(p) for p in pc]
    p2[n - 1][0] -= 1
    run(f"C1 p_{n}[0] one unit lower ({pc[n - 1][0]} -> {p2[n - 1][0]})", pc2=p2)
    w2 = [[x - 1 for x in w] for w in wc]
    run("C1b every boundary integer one unit lower (the uniform repair 11 -> 10)", wc2=w2)
    del w2
    l2 = list(lam)
    l2[0] -= 1
    l2[n - 1] += 1
    run(f"C2 one unit of weight moved from lambda_1 to lambda_{n} ({lam[0]} -> {l2[0]}, {lam[n - 1]} -> {l2[n - 1]}; "
        f"sum still D)", lam2=l2)
    l3 = list(lam)
    l3[n - 1] -= 1
    run(f"C2b lambda_{n} one unit lower", lam2=l3)
    claim = Fraction(floor_sqrt_dec(f["ab"], 20))
    c3 = not (f["ab"].numerator * claim.denominator ** 2 <= claim.numerator ** 2 * f["ab"].denominator)
    print(f"  C3 a claimed sqrt(ab) <= {floor_dec(claim, 20)} (just below the computed value): "
          f"{'not implied (ab > claim^2), refused, as expected' if c3 else 'IMPLIED (control failed)'}")
    res.append(c3)
    p4 = [list(p) for p in pc]
    p4[n - 1][0] -= 1
    p4[n - 1][1] += 1
    run(f"C4 one unit moved from p_{n}[0] to p_{n}[1] (mass kept)", pc2=p4)
    allok = all(res)
    print("controls all as expected" if allok else "CONTROLS FAILED")
    return allok


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    root = args[0] if args else os.path.join(HERE, "certs")
    t0 = time.time()
    rd = os.path.join(root, "akw26", "README_symmetric_certificates.md")
    rd_ok = sha256(rd) == SHA["README_symmetric_certificates.md"]
    print(f"[Akw26] John Akwei's symmetric certificates, Zenodo 10.5281/zenodo.24138215 v1.0; definitions from the V6 "
          f"report (Appendix D, Theorem 5.4, Table 10) and the README (sha256 as pinned: {rd_ok}); no deposit program "
          f"read or run")
    readme = readme_table(root)
    ok = rd_ok
    for row in TABLE10:
        r, _ = check_one(root, row, readme, "--controls" in sys.argv)
        ok &= r
    print(f"RESULT [Akw26]: all four Table 10 certificates {'PASS, as printed' if ok else 'NOT all as printed'}; by "
          f"Hou-Zhao Theorem 2.1 the m = 1024, L = 6 certificate gives C_5a <= sqrt(ab) < 0.9428095 "
          f"{'(shown)' if ok else '(see above)'}")
    print(f"run time {time.time() - t0:.0f} s")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
