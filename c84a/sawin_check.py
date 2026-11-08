"""84a: the row 1.0152* (u(n) > n^1.0152; Emmerich, arXiv 2606.03419v5, added in PR #315), replayed from the
certificate data alone with our own code (Quantavane replays).  Standard library only; run with python3 -I.  The
repository's verify_*.py and optimize_*.py were neither read nor run; only its data file is read, and its stored
results (the 'statuses', 'witnesses', 'checks', 'numerator', 'denominator', 'delta' fields) are compared with ours
only after our computation.

The criterion is Sawin's (arXiv 2605.20579, proof of Theorem 1 with Lemmas 11 and 12 and Proposition 10), which we
read and re-derived:
  (C1) T is a set of odd primes, and the number of q in T with q = 3 mod 4 is odd (so D = prod T = 3 mod 4 and the
       quadratic field Q = Q(sqrt D) has discriminant 4D);
  (C2) every p in S_Q is congruent to 1 mod 4 or inert in Q(sqrt q) for some q in T;
  (C3) the Golod-Shafarevich budget (9): #T + #S_Q + #{p in S_Q split in Q} + 1 <= (#T - 1)^2/4;
  (C4) R > 1 and every k(p) is a positive integer;
then u(n) >= n^(1+delta)/O(1) for arbitrarily large n, with (formula (11); e(p) = 2 for p = 2 or p in T, else 1)
  delta = [log(1 - 1/R) + log(2 pi/e)/2 + sum_p log(k(p)+1)/(4 e(p)) - log(4D)/8 - log(log sqrt(4D))/2]
          / log(2R prod_p p^(k(p)/(2 e(p))) + 1).
Splitting is decided by our own Legendre symbols (Euler's criterion) and Kronecker symbols at 2; delta is enclosed
in 60-digit decimal arithmetic (ln, exp, sqrt correctly rounded), each value widened by 1e-50; the product in the
denominator is the exact integer X = prod_{e=1} p^(2k) prod_{e=2} p^k, taken to the power 1/4.

What rests on mathematics (read, not computed): Sawin's Proposition 10 (Lemmas 2 to 9 behind it, Louboutin's class
number bound), Lemma 11 (the Golod-Shafarevich-Gaschuetz-Vinberg criterion and the Euler characteristic bound) and
Lemma 12 (the tower F, K = F(i), the root discriminant sqrt(4D), ramification and inertia).

usage: python3 -I sawin_check.py CERTS_ROOT [--controls]
"""
import hashlib
import json
import os
import sys
import time
from decimal import Decimal as D, getcontext
from fractions import Fraction as F

getcontext().prec = 60
W = D("1e-50")
DATA_SHA = "15fb28e33aae4918608e903a3d3b0a728d8ac567906d80146b5b1bdbd1f3ef86"


def is_prime(n):
    if n < 2:
        return False
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
    e = {p: (2 if (p == 2 or p in T) else 1) for p in S}
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


def check_cert(c, run, quiet=False):
    """returns (ok, lo, hi, details)."""
    T = [int(x) for x in c["T"]]
    S = [int(x) for x in c["S_Q"]]
    k = {int(p): int(v) for p, v in c["k"].items()}
    R = F(c["R"])
    ok = True
    ok &= run("T: %d distinct odd primes, %d of them = 3 mod 4 (odd)" % (len(T), sum(q % 4 == 3 for q in T)),
              len(set(T)) == len(T) and all(q > 2 and is_prime(q) for q in T) and sum(q % 4 == 3 for q in T) % 2 == 1)
    ok &= run("S_Q: %d distinct primes, k(p) a positive integer for each and only these, R = %s > 1" % (len(S), R),
              len(set(S)) == len(S) and all(is_prime(p) for p in S) and set(k) == set(S)
              and all(isinstance(k[p], int) and k[p] >= 1 for p in S) and R > 1)
    Dp = 1
    for q in T:
        Dp *= q
    dQ = disc_of(Dp)
    status = {p: {1: "split", -1: "inert", 0: "ramified"}[kron_disc(dQ, p)] for p in S}
    nsplit = sum(s == "split" for s in status.values())
    nram = sum(s == "ramified" for s in status.values()); nin = sum(s == "inert" for s in status.values())
    ok &= run("Q = Q(sqrt D), D = %d = 3 mod 4, discriminant 4D: in S_Q %d ramified, %d inert, %d split"
              % (Dp, nram, nin, nsplit), Dp % 4 == 3 and dQ == 4 * Dp)
    wit = {}
    for p in S:
        if p % 4 == 1:
            wit[p] = "1 mod 4"
        else:
            for q in sorted(T):
                if kron_disc(disc_of(q), p) == -1:
                    wit[p] = "inert in Q(sqrt %d)" % q
                    break
    ok &= run("(C2) every p in S_Q is 1 mod 4 or inert in some Q(sqrt q), q in T", len(wit) == len(S))
    lhs = len(T) + len(S) + nsplit + 1; rhs = F((len(T) - 1) ** 2, 4)
    ok &= run("(C3) budget #T + #S_Q + #split + 1 = %d <= (#T - 1)^2/4 = %s" % (lhs, rhs), lhs <= rhs)
    lo, hi, num, den = enclose_delta(T, S, k, R)
    return ok, lo, hi, dict(status=status, wit=wit, num=num, den=den)


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


def main():
    args = [x for x in sys.argv[1:] if not x.startswith("--")]
    if len(args) != 1:
        print(__doc__); sys.exit(2)
    t0 = time.time()
    path = os.path.join(args[0], "emm", "certificate_verification_all_v23.json")
    raw = open(path, "rb").read()
    run = Run()
    print("84a: Sawin's criterion on the four certificates of certificate_verification_all_v23.json (arXiv 2606.03419)")
    run("data file sha256 = 15fb28e3...ef86 (emmerichmtm/UnitDistanceProblemOptimizationOfSawinsLowerBound @ 2279d0a)",
        hashlib.sha256(raw).hexdigest() == DATA_SHA)
    certs = json.loads(raw, parse_float=str, parse_int=str)
    targets = {"sawin": D("0.014"), "rudolph_integer_es_discrete_recombination": D("0.0152"),
               "simple_rudolph_integer_es": D("0.0152"), "optimized": D("0.015")}
    results = {}
    for c in certs:
        key = c["key"]
        print(f"  -- {key}: {c['name']}")
        ok, lo, hi, det = check_cert(c, run)
        results[key] = (lo, hi, det, c)
        run("delta in [%s, %s] (numerator %s, denominator %s)" % (str(lo)[:24], str(hi)[:24], str(det["num"])[:20],
                                                                 str(det["den"])[:20]), ok and hi - lo < D("1e-40"))
        if key in targets:
            run("delta > %s, so u(n) > n^(1 + %s) for arbitrarily large n" % (targets[key], targets[key]),
                lo > targets[key])
        # only now compare with the file's stored results
        same_status = all(det["status"][int(p)] == s for p, s in c["statuses_in_Q_sqrt_D"].items())
        cl = D(c["claimed_delta"]); st = D(c["delta"])
        run("agrees with the file's stored splitting types, and its delta %s... and claimed %s lie in our enclosure "
            "to 1e-16" % (str(st)[:20], c["claimed_delta"]), same_status and abs(st - lo) < D("1e-40")
            and abs(cl - lo) < D("1e-16"))
    dr = results["rudolph_integer_es_discrete_recombination"]
    print(f"  the 84a row's certificate (discrete recombination, delta = 0.015263...): delta = {str(dr[0])[:22]}...")
    main_ok = not run.fails
    if "--controls" in sys.argv:
        print("  -- red controls (each must be refused)")
        base = dr[3]
        lo, hi = dr[0], dr[1]
        run("control: the exponent tightened to 1.0152629 (delta >= 0.0152629) is not implied", hi < D("0.0152629"))
        T = [int(x) for x in base["T"]]
        Dp = 1
        for q in T:
            Dp *= q
        p_split = next(p for p in range(3, 10 ** 5) if is_prime(p) and Dp % p and legendre(Dp, p) == 1
                       and str(p) not in base["k"])
        bad = json.loads(json.dumps(base)); bad["S_Q"] = list(bad["S_Q"]) + [str(p_split)]; bad["k"][str(p_split)] = "1"
        sub = Run(quiet=True); check_cert(bad, sub)
        run("control: adding %d (split in Q) breaks the budget (36 + 2 > 36)" % p_split, bool(sub.fails))
        p_bad = next(p for p in range(3, 10 ** 7, 4) if str(p) not in base["k"]
                     and all(legendre(q, p) != -1 for q in T if p % q) and is_prime(p))
        bad = json.loads(json.dumps(base)); i = bad["S_Q"].index("263"); bad["S_Q"][i] = str(p_bad)
        bad["k"][str(p_bad)] = bad["k"].pop("263")
        sub = Run(quiet=True); check_cert(bad, sub)
        run("control: 263 replaced by %d (3 mod 4, inert in no Q(sqrt q)) fails (C2)" % p_bad, bool(sub.fails))
        bad = json.loads(json.dumps(base)); bad["T"] = [q for q in bad["T"] if q != "43"]
        sub = Run(quiet=True); check_cert(bad, sub)
        run("control: T without 43 (six primes = 3 mod 4, an even number) is refused", bool(sub.fails))
        bad = json.loads(json.dumps(base)); bad["R"] = "1"
        sub = Run(quiet=True)
        try:
            check_cert(bad, sub)
        except Exception:
            sub.fails.append("R = 1: log(1 - 1/R) undefined")
        run("control: R = 1 is refused", bool(sub.fails))
    print(f"  checks: {run.n}, failed: {len(run.fails)}")
    print(f"  run time: {time.time() - t0:.1f} s")
    print("RESULT: ALL CHECKS PASS" if not run.fails else "RESULT: FAIL " + "; ".join(run.fails))
    sys.exit(0 if not run.fails else 1)


if __name__ == "__main__":
    main()
