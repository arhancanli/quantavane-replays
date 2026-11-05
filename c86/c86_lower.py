#!/usr/bin/env python3
"""Replay of the lower bound C86 >= 1.80220* (Schur-Siegel-Smyth trace constant; arXiv:2609.29298, D. Niedbala
Giraudin; the 86a row proposed in PR #315), from the certificate's data and the paper's Proposition 2 alone.

Proposition 2 (as the paper states it): if for every x >= 0
    x >= lambda + sum_{Q in A} lambda_Q log|Q(x)| + lambda_0 (2 U_eta(x) - I(eta)),
with lambda_Q >= 0, lambda_0 >= 0, Q integer polynomials and eta a probability measure of finite energy, then
lambda <= C86. We check:
  1. the data: integer monic Q, lambda_Q >= 0, lambda_0 >= 0, cells increasing from 0 to 8, masses >= 0; eta is the
     given measure scaled to total mass exactly 1 (the floats sum to 1 + 7.9e-17; the scaling is absorbed below);
  2. a rigorous lower bound for I(eta) (c86_lower.c energy);
  3. the inequality on [0, 14] by bisection (c86_lower.c check), with U_eta bounded cell by cell;
  4. the tail x >= 14: U_eta(x) <= log x, log|Q(x)| <= deg Q log x + log K_Q with K_Q = 1 + sum_i |a_i| 14^(i - deg),
     so g(x) >= h(x) = x - (sum lambda_Q deg Q + 2 lambda_0) log x - sum lambda_Q log K_Q + lambda_0 I, increasing for
     x >= 14; h(14) is computed in 60-digit decimal arithmetic.

Usage: python3 -I c86_lower.py ng [CERTS] [--controls] [--workers=N]   Python 3.9+, standard library only, and cc.
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
EXE = os.path.join(BIN, "c86_lower")
NG = dict(rel="ng/CERTIFICAT_C86.json", sha="3caa729730c5c52d755353ccf5e643372aadd0d792c37b14fcdb932ec6f1abe5",
          claim="1.80220", stated_I="3.1305148803e-5", X1=14)


def build():
    os.makedirs(BIN, exist_ok=True)
    src = os.path.join(HERE, "c86_lower.c")
    if not os.path.exists(EXE) or os.path.getmtime(EXE) < os.path.getmtime(src):
        subprocess.run(["cc", "-O2", "-std=c99", "-o", EXE, src, "-lm"], check=True)


def validate(c):
    P, L, lam0 = c["polynomes"], c["lambda_Q"], c["lambda_0"]
    e, m = c["aretes"], c["masses"]
    assert len(P) == len(L) and all(isinstance(a, int) for p in P for a in p), "integer polynomials"
    assert all(p[0] == 1 for p in P), "monic"
    assert all(x >= 0 for x in L) and lam0 >= 0, "lambda_Q >= 0, lambda_0 >= 0"
    assert len(e) == len(m) + 1 == c["n_cellules"] + 1, "cells"
    assert e[0] == 0 and e[-1] == 8 and all(a < b for a, b in zip(e, e[1:])), "edges increasing from 0 to 8"
    assert all(x >= 0 for x in m), "masses >= 0"
    return sum(Fraction(x) for x in m)


def data_text(c, lam0=None, masses=None):
    e, m = c["aretes"], masses if masses is not None else c["masses"]
    out = [str(len(m)), " ".join(float(x).hex() for x in e), " ".join(float(x).hex() for x in m), str(len(c["polynomes"]))]
    for p, l in zip(c["polynomes"], c["lambda_Q"]):
        out.append(f"{len(p) - 1} " + " ".join(float(a).hex() for a in p) + " " + float(l).hex())
    out.append(float(c["lambda_0"] if lam0 is None else lam0).hex())
    return "\n".join(out) + "\n"


def run(args, text):
    r = subprocess.run([EXE, *args], input=text, capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


def check_range(X1, Ilo, lam, text, workers):
    """certify the inequality on [0, X1], split into `workers` ranges run at once; returns (ok, summary)"""
    import tempfile
    fd, path = tempfile.mkstemp(suffix=".txt")
    with os.fdopen(fd, "w") as f:
        f.write(text)
    cuts = [X1 * i / workers for i in range(workers + 1)]
    procs = []
    for a, b in zip(cuts, cuts[1:]):
        h = open(path)
        procs.append((a, b, h, subprocess.Popen([EXE, "check", repr(a), repr(b), float(Ilo).hex(), repr(lam)],
                                                stdin=h, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)))
    ok, n, worst, worst_x, best, best_x, fails = True, 0, float("inf"), None, float("inf"), None, []
    for a, b, h, p in procs:
        out, _ = p.communicate()
        h.close()
        if p.returncode != 0 or not out.startswith("PASS"):
            ok = False
            fails.append(f"[{a:g}, {b:g}]: {out.strip()[:120]}")
            continue
        d = {}
        for line in out.strip().split("\n"):
            k, *v = line.split()
            d[k] = v
        n += int(d["PASS"][1])
        w = float(d["smallest_interval_bound"][0])
        if w < worst:
            worst, worst_x = w, float(d["smallest_interval_bound"][2])
        pv = float(d["smallest_point_value"][0])
        if pv < best:
            best, best_x = pv, float(d["smallest_point_value"][2])
    os.unlink(path)
    if not ok:
        return False, "FAIL " + "; ".join(fails)
    return True, (f"PASS on {n} intervals; smallest interval lower bound {worst:.10f} (at x = {worst_x:.6f}); "
                  f"smallest point value seen {best:.10f} (at x = {best_x:.6f})")


def tail_h14(c, Ilo):
    """h(14) in 60-digit decimal arithmetic (logs correctly rounded; one unit either way)."""
    ctx = decimal.Context(prec=60)
    lo = lambda v: ctx.next_minus(ctx.create_decimal(v)) if not isinstance(v, Decimal) else ctx.next_minus(v)
    ln_hi = lambda x: ctx.next_plus(ctx.ln(Decimal(x)))
    A = sum(Fraction(l) * (len(p) - 1) for p, l in zip(c["polynomes"], c["lambda_Q"]))
    B = 2 * Fraction(c["lambda_0"])
    logK = Decimal(0)
    for p, l in zip(c["polynomes"], c["lambda_Q"]):
        d = len(p) - 1
        K = 1 + sum(Fraction(abs(a), 14 ** i) for i, a in enumerate(p) if i > 0)    # p[i] multiplies x^(d - i)
        logK += Decimal(l) * ctx.next_plus(ctx.ln(Decimal(K.numerator) / Decimal(K.denominator)))
    h = (Decimal(14) - Decimal(A.numerator) / Decimal(A.denominator) * ln_hi(14) - Decimal(B.numerator) /
         Decimal(B.denominator) * ln_hi(14) - logK + Decimal(c["lambda_0"]) * Decimal(Ilo) - Decimal("1e-40"))
    slope = 1 - float(A + B) / 14
    return h, slope, float(A), float(B)


def check_ng(root, controls, workers):
    t0 = time.time()
    path = os.path.join(root, NG["rel"])
    raw = open(path, "rb").read()
    h = hashlib.sha256(raw).hexdigest()
    print(f"Row ng (arXiv:2609.29298 ancillary CERTIFICAT_C86.json, PR #315's 86a row): sha256 {h} "
          f"({'as fetched from arXiv v1' if h == NG['sha'] else 'NOT AS PINNED'})")
    assert h == NG["sha"]
    c = json.loads(raw)
    S = validate(c)
    P = c["polynomes"]
    print(f"  {len(P)} monic integer polynomials (degrees {[len(p) - 1 for p in P]}), "
          f"{sum(1 for l in c['lambda_Q'] if l > 0)} with lambda_Q > 0; lambda_0 = {c['lambda_0']}")
    print(f"  eta: {len(c['masses'])} cells on [0, 8], {sum(1 for x in c['masses'] if x > 0)} with positive mass; "
          f"total mass 1 + {float(S - 1):.2e} (scaled to exactly 1: U and I change by under 1e-15)")
    text = data_text(c)
    rc, out = run(["energy"], text)
    assert rc == 0, out
    vals = dict(line.split(None, 1) for line in out.strip().split("\n"))
    Ilo = float.fromhex(vals["lo"].strip())
    print(f"  I(eta) = {float(vals['I']):.12e} (error bound {vals['err'].strip()}); lower bound used {Ilo:.12e} "
          f"(the paper: {NG['stated_I']})")
    Ilo_used = Ilo - 1e-15                     # the scaling to mass 1
    h14, slope, A, B = tail_h14(c, Ilo_used)
    tail_ok = h14 > Decimal(NG["claim"]) and slope > 0
    print(f"  tail x >= 14: h(14) >= {float(h14):.6f}, h'(x) >= {slope:.4f} > 0 there "
          f"(sum lambda_Q deg Q = {A:.10f}, 2 lambda_0 = {B:.10f}): {'holds' if tail_ok else 'FAILS'}")
    lam = float(NG["claim"]) + 1e-12            # a little above the claim, for the scaling and the tail rounding
    good, msg = check_range(NG["X1"], Ilo_used, lam, text, workers)
    print("  [0, 14]: " + msg)
    ok = good and tail_ok
    print(f"  => for every x >= 0 the inequality holds with lambda = {NG['claim']}: C86 >= {NG['claim']} "
          f"{'(Proposition 2)' if ok else 'NOT shown'}")
    if controls:
        print("Red controls:")
        g1, out1 = check_range(NG["X1"], Ilo_used, c["certificat"] + 1e-6, text, workers)
        c1 = not g1
        print(f"  C1 lambda one millionth above the certificate's own value {c['certificat']}: "
              f"{'refused, as expected (' + out1[:110] + ')' if c1 else 'ACCEPTED (control failed)'}")
        g2, out2 = check_range(NG["X1"], 0.0, lam, data_text(c, lam0=0.0), workers)
        c2 = not g2
        print(f"  C2 the energy term dropped (lambda_0 = 0): {'refused at 1.80220, as expected' if c2 else 'ACCEPTED (control failed)'}")
        bad = dict(c)
        bad["lambda_Q"] = list(c["lambda_Q"])
        bad["lambda_Q"][3] = -bad["lambda_Q"][3]
        try:
            validate(bad)
            c3 = False
        except AssertionError as e:
            c3 = True
            print(f"  C3 one weight made negative: refused ({e})")
        if not c3:
            print("  C3 one weight made negative: ACCEPTED (control failed)")
        ok &= c1 and c2 and c3
    print(f"RESULT: ng {'PASS' if ok else 'FAIL'}")
    print(f"run time {time.time() - t0:.0f} s")
    return ok


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if not args or args[0] != "ng":
        print(__doc__)
        sys.exit(2)
    build()
    root = args[1] if len(args) > 1 else os.path.join(HERE, "certs")
    workers = next((int(a.split("=")[1]) for a in sys.argv[1:] if a.startswith("--workers=")), min(8, os.cpu_count() or 1))
    sys.exit(0 if check_ng(root, "--controls" in sys.argv, workers) else 1)


if __name__ == "__main__":
    main()
