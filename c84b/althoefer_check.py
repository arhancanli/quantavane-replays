"""84b: the row C84b <= 1.9993* (Althoefer's note of 28 May 2026 on [BSSZ2026, Section 5]; PR #275 restates the
ceiling of its chain), replayed with our own code (Quantavane replays).  Standard library only; run with python3 -I.

The note (improved_constant_052.tex, SHA-256 checked by fetch.sh; its numbers are read from the file and compared
with ours) keeps Section 5's construction A = G P (G the units of log-height <= Y, P = X + B^+(eps X)) and gets
    |A| >= (B X Y)^d,  |AA| <= (5 rho / Y)^d |A|^2,  |A+A| <= (M1 e^Y/(X Y^2) + M2/(X^2 Y^2))^d |A|^2,
    |A| <= e^(L d),  L > log(5Y + 1) + log(2 eps X + 1),
with eps = (phi - 1)/(phi + 1), rho^d bounding the regulator, B = eps/(rho sqrt C2), M2 = rho^2 C2/eps^2,
M1 = 4(1 + eps) M2 and C2 = 857.57; so max(|A+A|, |AA|) <= |A|^(2 - c) with c = -log(max(...))/L, up to (1+o(1))^d.
The regulator input is R_K <= (1+o(1))^d rho(s)^d for every s > 1, rho(s) = zeta(s) Gamma(s/2) C2^(s/2)/(2 pi^(s/2)).

Checks:
  1. the note's chain at its own parameters (s0 = 1.371966384, rho = 103, Y = 1415, X = floor(e^1423)):
     rho(s0) <= 103 (rigorous zeta and Gamma), eps = sqrt 5 - 2 and (1+eps)/(1-eps) = phi exactly, the printed
     B, M2, M1, 515/Y < 0.364, M1 e^Y/(X Y^2) < 0.136, the second additive term, L <= 1432, -log 0.364 > 0.0007 * 1432:
     so c >= 0.0007 and the exponent 2 - 0.0007 = 1.9993;
  2. the value of c at the note's choice (#275 and [T2026] say about 0.000706);
  3. the ceiling over every s > 1, X, Y > 0 (PR #275: c <= 0.0007150508): since -log q <= min(-log a, -log b) and
     L >= log(10 eps X Y), c <= g(Y, rho) = log(Y/(5 rho)) / (Y + log(8 (1+eps) rho C2 / eps)) for every X, which is
     decreasing in rho; so c <= max_Y g(Y, rho_lo) with rho_lo a rigorous lower bound of rho(s) over all s > 1
     (cells in s with monotone bounds: zeta decreasing, Gamma decreasing on (0, 1.46), the power increasing), and the
     maximum over Y bounded on cells as well.
zeta(s) by Euler-Maclaurin (N = 12, m = 9; f(x) = x^-s is completely monotone, so the remainder is bounded by the
first omitted term), log Gamma by Stirling's series at x + 12 (remainder bounded by the first omitted term) and the
recurrence; 50-digit decimal arithmetic (ln, exp correctly rounded), every enclosure widened by 1e-35.

What rests on mathematics (read, not computed): BSSZ2026 Section 5 (the tower with root discriminant <= C2, the unit
separation lemma, the lattice counts), the note's lattice-doubling lemma and its regulator inequality, and the
(1+o(1))^d absorption.

usage: python3 -I althoefer_check.py NOTE_DIR [--controls]
"""
import hashlib
import os
import re
import sys
from decimal import Decimal as D, getcontext, ROUND_FLOOR
from fractions import Fraction as Q

getcontext().prec = 50
W = D("1e-35")                      # widening of every enclosure
C2 = D("857.57")
BERN = [Q(1, 6), Q(-1, 30), Q(1, 42), Q(-1, 30), Q(5, 66), Q(-691, 2730), Q(7, 6), Q(-3617, 510), Q(43867, 798),
        Q(-174611, 330), Q(854513, 138)]        # B_2 .. B_22


def dq(q):
    return D(q.numerator) / D(q.denominator)


def machin_pi():
    def atan_inv(x, terms=60):
        s = Q(0)
        for k in range(terms):
            s += Q((-1) ** k, (2 * k + 1) * x ** (2 * k + 1))
        return s, Q(1, (2 * terms + 1) * x ** (2 * terms + 1))
    a, ea = atan_inv(5)
    b, eb = atan_inv(239)
    p = 16 * a - 4 * b
    err = 16 * ea + 4 * eb
    return dq(p), dq(err) + W


PI, PI_ERR = machin_pi()


def zeta(s, N=12, m=9):
    """(value, error bound) of zeta(s) for real s > 1, Euler-Maclaurin"""
    s = D(s)
    tot = sum(D(n) ** (-s) for n in range(1, N))
    tot += D(N) ** (1 - s) / (s - 1) + D(N) ** (-s) / 2
    rising = s                      # (s)_(2k-1)
    term = None
    for k in range(1, m + 2):
        term = dq(BERN[k - 1]) / D(factorial(2 * k)) * rising * D(N) ** (-s - 2 * k + 1)
        if k == m + 1:
            break
        tot += term
        rising *= (s + 2 * k - 1) * (s + 2 * k)
    return tot, abs(term) + W


def factorial(n):
    out = 1
    for i in range(2, n + 1):
        out *= i
    return out


def lgamma(x, shift=12, m=9):
    """(value, error bound) of log Gamma(x) for real x > 0, Stirling at z = x + shift and the recurrence"""
    x = D(x)
    z = x + shift
    val = (z - D("0.5")) * z.ln() - z + (2 * PI).ln() / 2
    term = None
    for k in range(1, m + 2):
        term = dq(BERN[k - 1]) / (D(2 * k) * D(2 * k - 1) * z ** (2 * k - 1))
        if k == m + 1:
            break
        val += term
    for j in range(shift):
        val -= (x + j).ln()
    return val, abs(term) + W + PI_ERR


def rho(s):
    """(lo, hi) of zeta(s) Gamma(s/2) C2^(s/2) / (2 pi^(s/2))"""
    s = D(s)
    z, ze = zeta(s)
    lg, lge = lgamma(s / 2)
    logrest = lg + (s / 2) * (C2.ln() - PI.ln()) - D(2).ln()
    err_rest = lge + (s / 2) * PI_ERR / PI + W
    lo = (z - ze) * (logrest - err_rest).exp() * (1 - W)
    hi = (z + ze) * (logrest + err_rest).exp() * (1 + W)
    return lo, hi


def read_note(note_dir):
    data = open(os.path.join(note_dir, "improved_constant_052.tex"), "rb").read()
    return hashlib.sha256(data).hexdigest(), data.decode("utf-8", "replace")


def chain_at(rho_v, Y, X_floor_exp, eps, verbose=True):
    """the note's chain at (rho, Y, X = floor(e^X_floor_exp)); returns c (a Decimal enclosure midpoint) and checks"""
    Y = D(Y)
    M2 = rho_v ** 2 * C2 / eps ** 2
    M1 = 4 * (1 + eps) * M2
    a = 5 * rho_v / Y
    # X = floor(e^t): e^Y / X <= e^(Y - t) / (1 - e^-t)
    t = D(X_floor_exp)
    eYoverX_hi = (Y - t).exp() / (1 - (-t).exp()) * (1 + W)
    b1 = M1 * eYoverX_hi / Y ** 2
    b2_log = M2.ln() - 2 * (t + (1 - (-t).exp()).ln()) - 2 * Y.ln()     # log of M2/(X^2 Y^2), X >= e^t - 1
    L_hi = (5 * Y + 1).ln() + t + (2 * eps + (-t).exp()).ln()             # log(2 eps X + 1) <= t + log(2 eps + e^-t)
    q = max(a, b1 + b2_log.exp())
    c = -q.ln() / L_hi
    return dict(M2=M2, M1=M1, B=eps / (rho_v * C2.sqrt()), a=a, b1=b1, b2_log=b2_log, L_hi=L_hi, c=c)


def ceiling(rho_lo, eps, Ymax=4000, step=D("0.05")):
    """max over Y of g(Y) = log(Y/(5 rho)) / (Y + K), K = log(8 (1+eps) rho C2 / eps), on cells [Ya, Yb]:
    g <= log(Yb/(5 rho)) / (Ya + K); g <= 0 for Y <= 5 rho; for Y >= Ymax, g <= log(Y/(5rho))/Y, decreasing."""
    K = (8 * (1 + eps) * rho_lo * C2 / eps).ln() * (1 - W)
    Y = 5 * rho_lo
    best, at = D(-1), None
    while Y < Ymax:
        Yb = Y + step
        g = (Yb / (5 * rho_lo)).ln() / (Y + K)
        if g > best:
            best, at = g, Y
        Y = Yb
    tail = (D(Ymax) / (5 * rho_lo)).ln() / D(Ymax)
    assert D(Ymax) > 5 * rho_lo * D(1).exp() and tail < best
    return best * (1 + W) + W, at, K


def record_ceiling(rho_f, eps, width=D("0.5"), Ymax=4000):
    """the looser ceiling of the record behind PR #275, as carlostoledo1891 describes it there (8 Nov): X only needs
    the additive factor below 1 (X > M1 e^Y / Y^2), so L > log(5Y + 1) + log(2 eps M1) + Y - 2 log Y; rho is the floor
    rho_f; each Y-cell [k w, (k+1) w] is bounded with the numerator log(Y/(5 rho)) at its right end and the denominator
    at its left end (both increase with Y for Y > 2). Cells from 5 rho_f (below it the numerator is negative) to Ymax;
    beyond Ymax the ratio is below log(Y/(5 rho))/(Y - 2 log Y), which decreases there. Returns (max, cell start)."""
    K = (8 * (1 + eps) * rho_f ** 2 * C2 / eps).ln()          # log(2 eps M1), M1 = 4 (1+eps) rho^2 C2 / eps^2
    k = int(5 * rho_f / width)
    best, at = D(-1), None
    while k * width < Ymax:
        Ya, Yb = k * width, (k + 1) * width
        if Yb > 5 * rho_f:
            den = ((5 * Ya + 1).ln() + K + Ya - 2 * Ya.ln()) * (1 - W)
            g = (Yb / (5 * rho_f)).ln() * (1 + W) / den
            if g > best:
                best, at = g, Ya
        k += 1
    Ym = D(Ymax)
    tail = (Ym / (5 * rho_f)).ln() / (Ym - 2 * Ym.ln())
    assert (Ym / (5 * rho_f)).ln() > D("1.01") and tail < best
    return best + W, at


def rho_lower_bound(verbose=False):
    """a lower bound of rho(s) over all s > 1, from cells; returns (rho_lo, where)"""
    lnratio = C2.ln() - PI.ln()          # log(C2/pi), up to PI_ERR/PI
    # s in (1, 1.02]: zeta(s) > 1/(s-1) >= 50, Gamma(s/2) >= Gamma(0.51) (Gamma decreasing on (0, 1.46)),
    # (C2/pi)^(s/2) >= (C2/pi)^(1/2)
    lg051, e051 = lgamma(D("0.51"))
    low1 = 50 * (lg051 - e051).exp() * (lnratio / 2).exp() / 2 * (1 - W)
    best, where = low1, "(1, 1.02]"
    # s in [1.02, 2.9]: cells; fine near the minimum
    s = D("1.02")
    edges = []
    while s < D("2.9"):
        h = D("0.00001") if D("1.36") <= s < D("1.385") else D("0.0002") if D("1.2") <= s < D("1.6") else D("0.002")
        edges.append(s)
        s += h
    edges.append(D("2.9"))
    cache = {}

    def zl(x):
        if x not in cache:
            z, ze = zeta(x)
            lg, lge = lgamma(x / 2)
            cache[x] = (z - ze, lg - lge)
        return cache[x]
    for sa, sb in zip(edges, edges[1:]):
        zb, lgb = zl(sb)
        low = zb * (lgb + (sa / 2) * (lnratio - PI_ERR / PI) - D(2).ln()).exp() * (1 - W)
        if low < best:
            best, where = low, f"[{sa}, {sb}]"
    # s >= 2.9: zeta >= 1, Gamma >= 0.885 (its minimum on (0, inf) is 0.8856...), (C2/pi)^(s/2) >= (C2/pi)^1.45
    low3 = D("0.885") * (D("1.45") * lnratio).exp() / 2
    assert low3 > best
    return best, where, len(edges) - 1


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if len(args) != 1:
        raise SystemExit(__doc__)
    h, tex = read_note(args[0])
    print(f"improved_constant_052.tex sha256 {h}")
    needed = ["C_2=857.57", "s=1.371966384", "101.9560759", "R_K\\le 103^d", "Y=1415", "X=\\lfloor e^{1423}\\rfloor",
              "=0.0000782645", "M_2=\\frac{103^2C_2}{\\varepsilon^2}=163256270.4", "M_1=4(1+\\varepsilon)M_2=807183391.9",
              "\\frac{515}{Y}=0.3639575971", "<0.136", "|A|\\le e^{1432d}", "-\\log(0.364)=1.010601",
              "0.0007\\cdot 1432=1.0024"]
    flat = re.sub(r"\s+", "", tex)
    missing = [n for n in needed if re.sub(r"\s+", "", n) not in flat]
    print(f"the note's numbers found in the file: {len(needed) - len(missing)} of {len(needed)}"
          + (f"; missing {missing}" if missing else ""))
    ok = not missing
    phi = (1 + D(5).sqrt()) / 2
    eps = (phi - 1) / (phi + 1)
    exact = abs(eps - (D(5).sqrt() - 2)) < W and abs((1 + eps) / (1 - eps) - phi) < W
    print(f"eps = (phi-1)/(phi+1) = sqrt 5 - 2 = {eps:.12f}, (1+eps)/(1-eps) = phi: {exact}")
    ok &= exact
    s0 = D("1.371966384")
    r0 = rho(s0)
    r_ok = r0[1] < 103
    print(f"rho({s0}) in [{r0[0]:.12f}, {r0[1]:.12f}] (the note: 101.9560759...); rho(s0) <= 103: {r_ok}")
    ok &= r_ok
    ch = chain_at(D(103), 1415, 1423, eps)
    checks = [
        ("B = eps/(103 sqrt C2) = 0.0000782645...", str(ch["B"]).startswith("0.0000782645")),
        ("M2 = 163256270.4... (the note's digits rounded to one decimal)", ch["M2"].quantize(D("0.1")) == D("163256270.4")),
        ("M1 = 807183391.9... (rounded to one decimal)", ch["M1"].quantize(D("0.1")) == D("807183391.9")),
        ("515/1415 < 0.364", ch["a"] < D("0.364")),
        ("M1 e^Y/(X Y^2) < 0.136", ch["b1"] < D("0.136")),
        ("M2/(X^2 Y^2) < 1e-1000", ch["b2_log"] < D(-1000) * D(10).ln()),
        ("L <= 1432 (log(5Y+1) + log(2 eps X + 1))", ch["L_hi"] <= 1432),
        ("-log(0.364) > 0.0007 * 1432", -D("0.364").ln() > D("0.0007") * 1432),
    ]
    for name, v in checks:
        print(f"  {name}: {v}")
        ok &= v
    print(f"  values: B = {ch['B']:.15f}, M2 = {ch['M2']:.6f}, M1 = {ch['M1']:.6f}, 515/Y = {ch['a']:.12f}, "
          f"M1 e^Y/(XY^2) <= {ch['b1']:.12f}, L <= {ch['L_hi']:.9f}")
    print(f"so c >= -log(0.364)/1432 = {-D('0.364').ln() / 1432:.12f} >= 0.0007: the exponent 2 - 0.0007 = 1.9993")
    print(f"c at the note's choice, without its roundings (rho = 103, Y = 1415, X = floor(e^1423)): {ch['c']:.12f}")
    # 3. the ceiling
    rho_lo, where, cells = rho_lower_bound()
    print(f"rho(s) >= {rho_lo:.9f} for every s > 1 ({cells} cells; the binding cell {where})")
    cmax, atY, K = ceiling(rho_lo, eps)
    print(f"ceiling: for every s > 1 and X, Y > 0, c <= {cmax:.12f} (cell at Y = {atY:.2f}; K = {K:.9f})")
    # how close the chain comes to the ceiling: rho(s0) (its upper end), Y = 1400.8, X = floor(e^t) at the crossing
    rr = r0[1]
    Yc = D("1400.8")
    M1c = 4 * (1 + eps) * rr ** 2 * C2 / eps ** 2
    tc = (M1c * Yc / (5 * rr * Yc ** 2)).ln() + Yc
    near = chain_at(rr, Yc, tc.quantize(D("0.001")), eps)["c"]
    print(f"  attained up to {cmax - near:.2E}: c = {near:.12f} at s = s0, Y = 1400.8, X = floor(e^{tc.quantize(D('0.001'))})")
    ceil_ok = cmax < D("0.0007150507")
    print(f"  so c <= 0.0007150508 (PR #275) holds, and even c < 0.0007150507: {ceil_ok}; 0.000719 is out of reach")
    ok &= ceil_ok and cmax < D("0.000719")
    # 4. the record's own ceiling, 0.00071505073..., from its three relaxations (carlostoledo1891 on #275, 8 Nov)
    rc, rat = record_ceiling(D("101.9"), eps)
    print(f"the record's looser chain (X with the additive factor below 1, rho = 101.9, Y-cells of width 1/2): "
          f"c <= {rc:.20f} (cell {rat} <= Y <= {rat + D('0.5')})")
    rec_ok = D("0.00071505073") <= rc < D("0.00071505074") and rc < D("0.0007150508") and rc > D("0.0007150507")
    print(f"  the record's 0.00071505073...: {rec_ok}; so #275's 0.0007150508 is it rounded up (0.0007150507 is below "
          f"it)")
    ok &= rec_ok
    if "--controls" in sys.argv:
        print("\nred controls (each must be refused):")
        c1 = not (rho(s0)[1] < D("101.95"))
        print(f"  rho(s0) <= 101.95 refused: {c1}")
        c2 = D(515) / D(1414) >= D("0.364")
        print(f"  Y = 1414 with the factor 0.364 kept: 515/1414 = {D(515) / D(1414):.10f} >= 0.364, refused: {c2}")
        ch3 = chain_at(D(103), 1415, 1421, eps)
        c3 = not (ch3["b1"] < D("0.136"))
        print(f"  X = floor(e^1421): M1 e^Y/(X Y^2) = {ch3['b1']:.6f} >= 0.136, refused: {c3}")
        c4 = not (-D("0.364").ln() > D("0.0007") * 1444)
        print(f"  L = 1444 in place of 1432: -log 0.364 > 0.0007 L refused: {c4}")
        c5 = not (ceiling(rho_lo, eps)[0] < D("0.0007138"))
        print(f"  a ceiling of 0.0007138, below the supremum, refused: {c5}")
        rq, _ = record_ceiling(D("101.953216"), eps)
        c6 = rq < D("0.00071505073")
        print(f"  the record's chain with rho = 101.953216 in place of 101.9 ({rq:.13f}) does not give 0.00071505073: {c6}")
        cok = c1 and c2 and c3 and c4 and c5 and c6
        print(f"red controls all as expected: {cok}")
        ok &= cok
    print(f"RESULT: {'PASS' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
