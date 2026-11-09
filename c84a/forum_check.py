"""84a: Sawin's own example and the three finite-data rows of the 84a page that reuse Sawin's argument, replayed with
our own code (Quantavane replays).  Standard library only; run with python3 -I.

  [S2026]     1.014       Sawin's example in the proof of Theorem 1 (arXiv 2605.20579): T (13 primes), S (22 primes),
                          the weights k(p) and R = 72, as printed in the paper (the paper: numerator 3.8822...,
                          denominator 275.055..., delta = .014114...);

  [EPF90]     1.03184...  mlewko, Erdos Problems forum thread 90 (21 May 2026): T (74 primes), S (1082 primes) and
                          the weights k(p) read from the data file the post links (read as data; its inertness
                          witness table is compared with ours only after our own search);
  [MO511514]  1.03188...  spiderduckpig (MathOverflow answer 511531 and the same forum thread, 21 May 2026): mlewko's
                          T with 347, 353 replaced by 503, 601; S = all primes p <= 13441 not split in K, and all
                          split primes p < 3433; the same weights as mlewko's file for p <= 577, and 1 above (its
                          numerator and denominator are compared with the MathOverflow answer's 144.928398401812
                          and 4545.623074469570; the forum copy prints 144.928398401807 and 4545.623074469564, 5e-12
                          and 6e-12 below ours);
  [MO511514]  1.03190...  spiderduckpig (same thread, 'Edit: Slight improvement'): T = the first 74 odd primes with
                          337, 349, 373 replaced by 389, 443, 587; 'the re-optimized data' with |S| = 1088, U = 169
                          and the stated weight histogram, rebuilt here as the 919 smallest primes not split in K and
                          the 169 smallest split primes, with the weights of the rule below.

No program of the posts was read or run (none was linked except the data file).  The criterion is Sawin's (arXiv
2605.20579, proof of Theorem 1 with Lemmas 11 and 12 and Proposition 10), which we read and re-derived, as in our
replay of the 1.0152* row (sawin_check.py in this folder):
  (C1) T is a set of odd primes, and the number of q in T with q = 3 mod 4 is odd (so D = prod T = 3 mod 4 and
       K = Q(sqrt D) has discriminant 4D);
  (C2) every p in S is 1 mod 4 or inert in Q(sqrt q) for some q in T;
  (C3) the Golod-Shafarevich budget (9): #T + #S + #{p in S split in K} + 1 <= (#T - 1)^2/4;
  (C4) R > 1 and every k(p) is a positive integer;
then u(n) >= n^(1 + delta)/O(1) for arbitrarily large n, with (formula (11); e(p) = 2 for p = 2 or p in T, else 1)
  delta = [log(1 - 1/R) + log(2 pi/e)/2 + sum_p log(k(p)+1)/(4 e(p)) - log(4D)/8 - log(log sqrt(4D))/2]
          / log(1 + 2R prod_p p^(k(p)/(2 e(p)))).
Splitting and inertness are decided by our own Legendre symbols (Euler's criterion) and Kronecker symbols at 2;
delta is enclosed in 60-digit decimal arithmetic (ln, exp, sqrt correctly rounded), each value widened by 1e-50; the
product in the denominator is the exact integer X = prod_{e=1} p^(2k) prod_{e=2} p^k, taken to the power 1/4.
R is the decimal printed in each post, read as an exact rational.

The weights: for a fixed ratio lambda, the term of p in (numerator - lambda * denominator) is
log(k+1)/(4e) - lambda k log(p)/(2e), which is largest at the least k >= 1 with log((k+2)/(k+1)) <= 2 lambda log p
(e cancels).  We check that each row's weights are exactly these at lambda = its own delta.

What rests on mathematics (read, not computed): Sawin's Proposition 10 (Lemmas 2 to 9 behind it, Louboutin's class
number bound), Lemma 11 (the Golod-Shafarevich-Gaschuetz-Vinberg criterion and the Euler characteristic bound) and
Lemma 12 (the tower F, K = F(i), the root discriminant sqrt(4D), ramification and inertia).

usage: python3 -I forum_check.py CERTS_ROOT [--controls]
"""
import hashlib
import os
import re
import sys
import time
from decimal import Decimal as D, getcontext
from fractions import Fraction as F
from math import log

getcontext().prec = 60
W = D("1e-50")
DATA_SHA = "2783627ddb767c7e4c98d1d0da5a16d628b1652d83883398453b63db95b03f14"
LIMIT = 20000
SIEVE = bytearray([1]) * (LIMIT + 1)
SIEVE[0] = SIEVE[1] = 0
for _i in range(2, int(LIMIT ** 0.5) + 1):
    if SIEVE[_i]:
        SIEVE[_i * _i::_i] = bytearray(len(SIEVE[_i * _i::_i]))
PRIMES = [i for i in range(LIMIT + 1) if SIEVE[i]]


def is_prime(n):
    if n <= LIMIT:
        return n >= 0 and bool(SIEVE[n])
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


def legendre(a, p):
    """(a/p) for an odd prime p, by Euler's criterion."""
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def kron_disc(d, p):
    """Kronecker symbol (d/p) for a fundamental discriminant d and a prime p: +1 split, -1 inert, 0 ramified."""
    if p == 2:
        if d % 2 == 0:
            return 0
        return 1 if d % 8 in (1, 7) else -1
    return legendre(d, p)


def disc_of(m):
    """discriminant of Q(sqrt m) for squarefree m."""
    return m if m % 4 == 1 else 4 * m


def pi_dec():
    def arctan_inv(n):
        x = D(1) / n; x2 = x * x; term = x; s = x; k = 1
        while True:
            term *= -x2; k += 2
            t = term / k
            if abs(t) < D("1e-58"):
                return s
            s += t
    return 16 * arctan_inv(5) - 4 * arctan_inv(239)


PI = pi_dec()


def enclose_delta(T, S, k, R):
    """returns (lo, hi, numerator, denominator) for formula (11)."""
    Dp = 1
    for q in T:
        Dp *= q
    Tset = set(T)
    e = {p: (2 if (p == 2 or p in Tset) else 1) for p in S}
    Rd = D(R.numerator) / D(R.denominator)
    num = (1 - 1 / Rd).ln() + (2 * PI / D(1).exp()).ln() / 2
    for p in S:
        num += D(k[p] + 1).ln() / (4 * e[p])
    lam2 = D(4 * Dp)                                   # lambda^2 = 4D
    num -= lam2.ln() / 8 + (lam2.sqrt().ln()).ln() / 2
    X = 1
    for p in S:
        X *= p ** (2 * k[p]) if e[p] == 1 else p ** k[p]
    prod = D(X).sqrt().sqrt()                          # prod p^(k/(2e)) = X^(1/4)
    den = (2 * Rd * prod + 1).ln()
    nlo, nhi = num - W * (1 + abs(num)), num + W * (1 + abs(num))
    dlo, dhi = den - W * den, den + W * den
    lo = nlo / dhi if nlo > 0 else nlo / dlo
    hi = nhi / dlo if nhi > 0 else nhi / dhi
    return lo - W, hi + W, num, den


def best_k(p, lam):
    """the least k >= 1 maximizing log(k+1) - 2 lam k log p, in floating point (the weights affect only the value of
    delta, never the validity of the criterion)."""
    k = 1
    while log((k + 2) / (k + 1)) > 2 * lam * log(p):
        k += 1
    return k


def status_in(T, p):
    Dp = 1
    for q in T:
        Dp *= q
    return {1: "split", -1: "inert", 0: "ramified"}[kron_disc(disc_of(Dp), p)]


def check_cert(name, T, S, k, R, run, quiet=False):
    """the conditions (C1)-(C4) and the enclosure; returns (ok, lo, hi, details)."""
    ok = True
    n3 = sum(q % 4 == 3 for q in T)
    ok &= run("%s: T has %d distinct odd primes, %d of them = 3 mod 4 (odd)" % (name, len(T), n3),
              len(set(T)) == len(T) and all(q > 2 and is_prime(q) for q in T) and n3 % 2 == 1)
    ok &= run("%s: S has %d distinct primes, k(p) a positive integer for each and only these, R = %s > 1"
              % (name, len(S), float(R)),
              len(set(S)) == len(S) and all(is_prime(p) for p in S) and set(k) == set(S)
              and all(isinstance(k[p], int) and k[p] >= 1 for p in S) and R > 1)
    Dp = 1
    for q in T:
        Dp *= q
    status = {p: status_in(T, p) for p in S}
    nsplit = sum(s == "split" for s in status.values())
    ok &= run("%s: D = %d mod 4, so K = Q(sqrt D) has discriminant 4D; in S %d ramified, %d inert, %d split"
              % (name, Dp % 4, sum(s == "ramified" for s in status.values()),
                 sum(s == "inert" for s in status.values()), nsplit), Dp % 4 == 3)
    wit = {}
    for p in S:
        if p % 4 == 1:
            continue
        for q in sorted(T):
            if kron_disc(disc_of(q), p) == -1:
                wit[p] = q
                break
    need = [p for p in S if p % 4 != 1]
    ok &= run("%s: (C2) each of the %d primes of S that are not 1 mod 4 is inert in some Q(sqrt q), q in T"
              % (name, len(need)), all(p in wit for p in need))
    lhs = len(T) + len(S) + nsplit + 1; rhs = F((len(T) - 1) ** 2, 4)
    ok &= run("%s: (C3) budget #T + #S + #split + 1 = %d <= (#T - 1)^2/4 = %s" % (name, lhs, rhs), lhs <= rhs)
    lo, hi, num, den = enclose_delta(T, S, k, R) if R > 1 else (D(-1), D(-1), D(0), D(1))
    return ok, lo, hi, dict(status=status, wit=wit, num=num, den=den, nsplit=nsplit, n3=n3)


def hist(k):
    h = {}
    for v in k.values():
        h[v] = h.get(v, 0) + 1
    return dict(sorted(h.items(), reverse=True))


def parse_list(text, start):
    i = text.index("[", start)
    j = text.index("]", i)
    body = text[i + 1:j]
    if not re.fullmatch(r"[\d,\s]*", body):
        raise ValueError("unexpected characters in a list")
    return [int(x) for x in re.findall(r"\d+", body)], j


def parse_file(raw):
    text = raw.decode("ascii")
    T, _ = parse_list(text, text.index("T = ["))
    S, _ = parse_list(text, text.index("S = ["))
    k = {}
    for m in re.finditer(r"^k=(\d+):", text, re.M):
        lst, _ = parse_list(text, m.end())
        for p in lst:
            if p in k:
                raise ValueError("a prime has two weights")
            k[p] = int(m.group(1))
    wit = {}
    for m in re.finditer(r"^q=(\d+) witnesses \(count: (\d+)\):", text, re.M):
        lst, _ = parse_list(text, m.end())
        if len(lst) != int(m.group(2)):
            raise ValueError("a witness count differs from its list")
        for p in lst:
            wit[p] = int(m.group(1))
    return T, S, k, wit


class Run:
    def __init__(self, quiet=False):
        self.n = 0; self.fails = []; self.quiet = quiet

    def __call__(self, label, ok):
        self.n += 1
        if not ok:
            self.fails.append(label)
        if not self.quiet:
            print(f"  check {label}: {'PASS' if ok else 'FAIL'}")
        return ok


def first_odd_primes(n):
    return [p for p in PRIMES if p > 2][:n]


def main():
    args = [x for x in sys.argv[1:] if not x.startswith("--")]
    if len(args) != 1:
        print(__doc__); sys.exit(2)
    t0 = time.time()
    path = os.path.join(args[0], "epf90", "mlewko_sawin_data.txt")
    raw = open(path, "rb").read()
    run = Run()
    print("84a: Sawin's example [S2026] 1.014 and the finite-data rows [EPF90] 1.03184..., [MO511514] 1.03188..., "
          "1.03190... on Sawin's criterion")
    run("data file sha256 = 2783627d...3f14 (the file linked in mlewko's post, Google Drive 1ygqjOSa...)",
        hashlib.sha256(raw).hexdigest() == DATA_SHA)
    T0, S0, k0, wit_file = parse_file(raw)
    run("data file: |T| = %d, flat S has %d primes, %d with weights, %d witnesses listed"
        % (len(T0), len(S0), len(k0), len(wit_file)),
        len(T0) == 74 and len(S0) == 1082 and sorted(S0) == S0 and set(S0) == set(k0))
    rows = []

    # ---- [S2026]: Sawin's own example (arXiv:2605.20579, proof of Theorem 1), as printed in the paper
    TS = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]
    SS = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 47, 71, 79, 97, 101, 107, 109, 139, 151, 163, 167, 179]
    kS = dict(zip(SS, [50, 31, 21, 17, 14, 13, 12, 11, 10, 10, 8, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6]))
    rows.append(dict(name="[S2026]", T=TS, S=SS, k=kS, R=F(72), page=D("0.014"), claim="0.014114",
                     post=dict(nS=22, U=0, n3=7, num="3.8822", den="275.055", hist=hist(kS))))

    # ---- [EPF90]: the file as it is
    R0 = F("32.39719232352239")
    rows.append(dict(name="[EPF90]", T=T0, S=S0, k=k0, R=R0, page=D("0.03184"), claim="0.031849981670202",
                     post=dict(nS=1082, U=175, n3=39, num="143.25907766230074", den="4497.932813453685",
                               hist={22: 1, 13: 1, 9: 1, 7: 1, 6: 1, 5: 2, 4: 4, 3: 13, 2: 82, 1: 976})))

    # ---- [MO511514] 1.03188: mlewko's T with 347, 353 -> 503, 601; S by the post's thresholds
    T1 = sorted([q for q in T0 if q not in (347, 353)] + [503, 601])
    S1 = [p for p in PRIMES if (status_in(T1, p) == "split" and p < 3433)
          or (status_in(T1, p) != "split" and p <= 13441)]
    k1 = {p: (k0[p] if p <= 577 else 1) for p in S1}
    rows.append(dict(name="[MO511514] 1.03188", T=T1, S=S1, k=k1, R=F("32.36456381771092"),
                     page=D("0.03188"), claim="0.03188306553964",
                     post=dict(nS=1092, U=165, n3=39, num="144.928398401812", den="4545.623074469570",
                               hist={22: 1, 13: 1, 9: 1, 7: 1, 6: 1, 5: 2, 4: 4, 3: 13, 2: 82, 1: 986})))

    # ---- [MO511514] 1.03190: the first 74 odd primes with 337, 349, 373 -> 389, 443, 587
    T2 = sorted([q for q in first_odd_primes(74) if q not in (337, 349, 373)] + [389, 443, 587])
    nonsplit = [p for p in PRIMES if status_in(T2, p) != "split"][:1088 - 169]
    split = [p for p in PRIMES if status_in(T2, p) == "split"][:169]
    S2 = sorted(nonsplit + split)
    lam2 = 0.03190837173526785
    k2 = {p: best_k(p, lam2) for p in S2}
    rows.append(dict(name="[MO511514] 1.03190", T=T2, S=S2, k=k2, R=F("32.33973768065121"),
                     page=D("0.03190"), claim="0.03190837173526785",
                     post=dict(nS=1088, U=169, n3=41, num="144.15620697127656", den="4517.817711517471",
                               hist={22: 1, 13: 1, 9: 1, 7: 1, 6: 1, 5: 2, 4: 4, 3: 13, 2: 81, 1: 983})))
    run("[MO511514] 1.03188: T differs from mlewko's in 347, 353 -> 503, 601; S = the primes <= 13441 not split in "
        "K and the split primes < 3433 (largest %d and %d)" % (max(p for p in S1 if status_in(T1, p) != "split"),
                                                                max(p for p in S1 if status_in(T1, p) == "split")),
        len(T1) == 74 and set(T1) ^ set(T0) == {347, 353, 503, 601})
    run("[MO511514] 1.03190: T = the first 74 odd primes with 337, 349, 373 -> 389, 443, 587; S = the 919 smallest "
        "primes not split in K (to %d) and the 169 smallest split primes (to %d)" % (max(nonsplit), max(split)),
        len(T2) == 74 and set(T2) ^ set(first_odd_primes(74)) == {337, 349, 373, 389, 443, 587})

    results = {}
    for r in rows:
        name, T, S, k, R, post = r["name"], r["T"], r["S"], r["k"], r["R"], r["post"]
        print(f"  -- {name}")
        ok, lo, hi, det = check_cert(name, T, S, k, R, run)
        run("%s: |S| = %d, U = %d, %d primes of T are 3 mod 4, weight histogram %s: all as the post states"
            % (name, len(S), det["nsplit"], det["n3"], hist(k)),
            (len(S), det["nsplit"], det["n3"], hist(k)) == (post["nS"], post["U"], post["n3"], post["hist"]))
        lam = float(lo)
        kmin = min(k.values())
        run("%s: every weight is the best one for its prime at lambda = delta (k(2) = %d, k(3) = %d, k = %d from p = %d)"
            % (name, k[2], k[3], kmin, min(p for p in S if k[p] == kmin)), all(k[p] == best_k(p, lam) for p in S))
        nd = len(post["num"].split(".")[1]); dd = len(post["den"].split(".")[1])
        run("%s: numerator %s and denominator %s, the post's %s and %s to all %d and %d printed decimals"
            % (name, str(det["num"])[:20], str(det["den"])[:20], post["num"], post["den"], nd, dd),
            abs(det["num"] - D(post["num"])) <= D(10) ** -nd and abs(det["den"] - D(post["den"])) <= D(10) ** -dd)
        run("%s: delta in [%s, %s]" % (name, str(lo)[:26], str(hi)[:26]), ok and hi - lo < D("1e-40"))
        cd = len(r["claim"].split(".")[1])
        run("%s: the post's delta %s is our delta rounded to its %d printed decimals (%s)"
            % (name, r["claim"], cd, "rounded up" if D(r["claim"]) > hi else "rounded down"),
            abs(D(r["claim"]) - lo) <= D(10) ** -cd / 2)
        run("%s: delta > %s, so u(n) > n^(1 + %s)/O(1) for arbitrarily large n: the page's 1%s holds"
            % (name, r["page"], r["page"], str(r["page"])[1:]), lo > r["page"])
        results[name] = (lo, hi, det, r)

    # only now compare our inertness witnesses with the file's table (EPF90)
    det0 = results["[EPF90]"][2]
    fine = all(kron_disc(disc_of(q), p) == -1 for p, q in wit_file.items())
    run("[EPF90]: the file's witness table (%d entries, compared after our search) lists exactly the %d primes of S "
        "that are not 1 mod 4, each with a valid witness" % (len(wit_file), len(det0["wit"])),
        set(wit_file) == set(det0["wit"]) and fine)
    print("  the four rows: 1.014, 1.03184, 1.03188, 1.03190 hold (delta = %s..., %s..., %s..., %s...)"
          % tuple(str(results[n][0])[:16] for n in ("[S2026]", "[EPF90]", "[MO511514] 1.03188", "[MO511514] 1.03190")))

    if "--controls" in sys.argv:
        print("  -- red controls (each must be refused)")
        for name in results:
            lo, hi, det, r = results[name]
            tight = hi + D("1e-15")
            run("control: %s: delta >= %s (1e-15 above our enclosure) is not implied" % (name, str(tight)[:20]),
                not (hi >= tight))
        lo, hi, det, r = results["[EPF90]"]
        T, S, k = r["T"], list(r["S"]), dict(r["k"])
        p_split = next(p for p in PRIMES if p not in k and status_in(T, p) == "split")
        sub = Run(quiet=True); check_cert("x", T, S + [p_split], {**k, p_split: 1}, r["R"], sub)
        run("control: [EPF90]: the next split prime %d added breaks the budget (1332 + 2 > 1332.25)" % p_split,
            bool(sub.fails))
        p_ns = next(p for p in PRIMES if p not in k and status_in(T, p) != "split")
        sub = Run(quiet=True); check_cert("x", T, S + [p_ns], {**k, p_ns: 1}, r["R"], sub)
        run("control: [EPF90]: the next non-split prime %d added breaks the budget (1332 + 1 > 1332.25)" % p_ns,
            bool(sub.fails))
        Tb = [q for q in T if q != 631] + [641]
        sub = Run(quiet=True); check_cert("x", Tb, S, k, r["R"], sub)
        run("control: [EPF90]: T with 631 replaced by 641 (38 primes = 3 mod 4, an even number) is refused",
            bool(sub.fails))
        Tc = [q for q in T if q % 8 != 5]
        sub = Run(quiet=True); check_cert("x", Tc, S, k, r["R"], sub)
        run("control: [EPF90]: T without its %d primes = 5 mod 8 leaves p = 2 inert in no Q(sqrt q) (C2 fails)"
            % (len(T) - len(Tc)), bool(sub.fails) and not any(kron_disc(disc_of(q), 2) == -1 for q in Tc))
        kb = dict(k); kb[3] = 14
        lo_b, hi_b, num_b, den_b = enclose_delta(T, S, kb, r["R"])
        run("control: [EPF90]: k(3) = 14 instead of 13 lowers delta (%s < %s) and moves the numerator off the post's"
            % (str(hi_b)[:12], str(lo)[:12]), hi_b < lo and abs(num_b - D(r["post"]["num"])) > D("1e-6"))
        sub = Run(quiet=True)
        _, lo_r, hi_r, _ = check_cert("x", T, S, k, F(1), sub)
        run("control: [EPF90]: R = 1 is refused", bool(sub.fails))
        r1 = results["[MO511514] 1.03188"][3]
        S1b = sorted(set(r1["S"]) | {next(p for p in PRIMES if p > 13441 and status_in(r1["T"], p) != "split")})
        sub = Run(quiet=True); check_cert("x", r1["T"], S1b, {**r1["k"], S1b[-1]: 1}, r1["R"], sub)
        run("control: [MO511514] 1.03188: the threshold 13441 raised to the next non-split prime %d breaks the budget"
            % S1b[-1], bool(sub.fails))
    print(f"  checks: {run.n}, failed: {len(run.fails)}")
    print(f"  run time: {time.time() - t0:.1f} s")
    print("RESULT: ALL CHECKS PASS" if not run.fails else "RESULT: FAIL " + "; ".join(run.fails))
    sys.exit(0 if not run.fails else 1)


if __name__ == "__main__":
    main()
