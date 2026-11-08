"""63a and 64a: the rows theta* = 0.3144831759741... (C63, Dirichlet divisor exponent, PR #313) and
2 theta* = 0.6289663519... (C64, Gauss circle exponent, PR #315), both from Li and Yang, arXiv 2308.14859v2,
Theorem 1.2 and Definition 1.1, replayed with our own code (Quantavane replays).  Standard library only; run with
python3 -I.

Definition 1.1 of the paper: -theta* is the unique solution on [-0.35, -0.3] of equation (1.4),
    -(8/25) x - (1/200) (sqrt(2(1-14x)) - 5 sqrt(-1-8x))^2 + 51/200 = -x.
The LaTeX of (1.4) is read from the paper's source (main.tex of the arXiv e-print, SHA-256 checked by fetch.sh)
and compared with our transcription before anything is computed.

Route A (exact algebra).  With u = sqrt(2(1-14x)), v = sqrt(-1-8x), (u - 5v)^2 = -23 - 228x - 10uv, so (1.4) is
364x + 74 + 10uv = 0.  Squaring gives 27524 x^2 + 13168 x + 1419 = 0, roots (-3292 +- 25 sqrt(1717))/13762; a root
of the squared equation solves (1.4) exactly when uv = -(364x + 74)/10 >= 0.  Every identity here is checked by
exact polynomial arithmetic, and the sign conditions by exact rational comparisons.
Route B (interval arithmetic, independent of route A).  Fixed point with 2^-400 units and outward rounding (integer
square roots rounded down for lower ends, up for upper ends): f(x) = LHS + x has f' > 0 on the whole interval (one
interval evaluation of f'), so (1.4) has at most one solution there; f changes sign at the ends; bisection encloses
the root.  Route B's enclosure must contain route A's root.

The printed values are then compared with the exact value: an upper bound printed in decimals must not be below it.

What rests on mathematics (read, not computed): the proof of Theorem 1.2 (R(X), Delta(X) = O(X^(theta* + eps))),
and, for 64a, that the page's radius normalization turns R(X) = O(X^(theta*+eps)) with X = t^2 into the exponent
2 theta* (N(t) - pi t^2 = R(t^2)).

usage: python3 -I theta_check.py PAPER_DIR [--controls]
"""
import hashlib
import math
import os
import re
import sys
from fractions import Fraction as Q

TRANSCRIPTION = r"-\frac{8}{25}x-\frac{1}{200}\Big(\sqrt{2(1-14x)}-5\sqrt{-1-8x}\Big)^2+\frac{51}{200}=-x"

# the equation as data: LHS = a1*x - c*(sqrt(p1 + p2 x) - k sqrt(q1 + q2 x))^2 + a0, solved for LHS = -x
EQ = dict(a1=Q(-8, 25), c=Q(1, 200), p1=2, p2=-28, k=5, q1=-1, q2=-8, a0=Q(51, 200))
INTERVAL = (Q(-35, 100), Q(-3, 10))

PRINTED = {  # row: (printed decimal, multiple of theta*, side)
    "63a (PR #313, README cell and row)": ("0.3144831759741", 1, "upper"),
    "64a (PR #315, README cell)": ("0.6289663519", 2, "upper"),
}


def read_equation(paper_dir):
    path = os.path.join(paper_dir, "main.tex")
    data = open(path, "rb").read()
    h = hashlib.sha256(data).hexdigest()
    text = data.decode("utf-8", "replace")
    m = re.search(r"\\label\{definition of theta\}\s*(.*?)\\end\{equation\}", text, re.S)
    if not m:
        raise SystemExit("FAIL: equation (1.4) not found in main.tex")
    tex = re.sub(r"\s+", "", m.group(1))
    return h, tex


# ---------------------------------------------------------------- route A: exact algebra in Q(sqrt 1717)
def route_a(eq, interval):
    """Return the exact solutions of (1.4) in the interval, as (r, s, d): x = r + s*sqrt(d)."""
    a1, c, p1, p2, k, q1, q2, a0 = (eq[n] for n in ("a1", "c", "p1", "p2", "k", "q1", "q2", "a0"))
    # (u - k v)^2 = u^2 + k^2 v^2 - 2k uv with u^2 = p1 + p2 x, v^2 = q1 + q2 x.
    # (1.4): a1 x - c(p1 + p2 x + k^2 (q1 + q2 x)) + 2ck uv + a0 + x = 0, i.e. uv = (A x + B) with
    A = -(a1 + 1 - c * (p2 + k * k * q2)) / (2 * c * k)
    B = -(a0 - c * (p1 + k * k * q1)) / (2 * c * k)
    # square: (p1 + p2 x)(q1 + q2 x) = (A x + B)^2  ->  qa x^2 + qb x + qc = 0
    qa = A * A - p2 * q2
    qb = 2 * A * B - (p1 * q2 + p2 * q1)
    qc = B * B - p1 * q1
    den = math.lcm(qa.denominator, qb.denominator, qc.denominator)
    ia, ib, ic = int(qa * den), int(qb * den), int(qc * den)
    g = math.gcd(math.gcd(ia, ib), ic)
    ia, ib, ic = ia // g, ib // g, ic // g
    disc = ib * ib - 4 * ia * ic
    assert disc > 0
    # disc = f^2 * d with d squarefree part (enough: pull out square factors by trial division)
    f, d = 1, disc
    p = 2
    while p * p <= d:
        while d % (p * p) == 0:
            d //= p * p
            f *= p
        p += 1
    sols, rejected = [], []
    for sign in (1, -1):
        # x = (-ib + sign f sqrt(d)) / (2 ia)
        r, s = Q(-ib, 2 * ia), Q(sign * f, 2 * ia)
        # exact checks, via comparisons in Q(sqrt d): t = r + s sqrt d
        lo_ok = cmp_qsqrt(r - interval[0], s, d) >= 0
        hi_ok = cmp_qsqrt(r - interval[1], s, d) <= 0
        uv = (A * r + B, A * s)                      # uv = A x + B
        uv_ok = cmp_qsqrt(uv[0], uv[1], d) >= 0
        dom = cmp_qsqrt(p1 + p2 * r, p2 * s, d) >= 0 and cmp_qsqrt(q1 + q2 * r, q2 * s, d) >= 0
        rec = dict(r=r, s=s, d=d, in_interval=lo_ok and hi_ok, uv_nonneg=uv_ok, domain=dom)
        (sols if (lo_ok and hi_ok and uv_ok and dom) else rejected).append(rec)
    return (ia, ib, ic), (A, B), sols, rejected


def cmp_qsqrt(r, s, d):
    """Sign of r + s*sqrt(d) (d > 0 not a square), exactly."""
    r, s = Q(r), Q(s)
    if s == 0:
        return (r > 0) - (r < 0)
    if r == 0:
        return (s > 0) - (s < 0)
    if (r > 0) == (s > 0):
        return 1 if r > 0 else -1
    # opposite signs: compare r^2 with s^2 d
    big = r * r - s * s * d
    if big == 0:
        return 0
    return (1 if r > 0 else -1) if big > 0 else (1 if s > 0 else -1)


def check_residual(eq, r, s, d):
    """(1.4) at x = r + s sqrt d, exactly: uv must equal A x + B and u^2 v^2 = (Ax+B)^2 (both in Q(sqrt d))."""
    p1, p2, q1, q2 = eq["p1"], eq["p2"], eq["q1"], eq["q2"]

    def mul(a, b):
        return (a[0] * b[0] + a[1] * b[1] * d, a[0] * b[1] + a[1] * b[0])
    u2 = (p1 + p2 * r, p2 * s)
    v2 = (q1 + q2 * r, q2 * s)
    return mul(u2, v2)


# ---------------------------------------------------------------- route B: fixed-point interval arithmetic
SH = 400
ONE = 1 << SH


def fx_floor(q):
    return (q.numerator << SH) // q.denominator


def fx_ceil(q):
    return -((-q.numerator << SH) // q.denominator)


def iv_const(q):
    return (fx_floor(Q(q)), fx_ceil(Q(q)))


def iv_add(a, b):
    return (a[0] + b[0], a[1] + b[1])


def iv_neg(a):
    return (-a[1], -a[0])


def iv_mul(a, b):
    ps = [a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1]]
    return (min(ps) >> SH, -((-max(ps)) >> SH))


def iv_div(a, b):
    assert b[0] > 0 or b[1] < 0, "division by an interval containing 0"
    qs = []
    for x in a:
        for y in b:
            qs.append(Q(x, y))
    return (fx_floor(min(qs)), fx_ceil(max(qs)))


def iv_sqrt(a):
    assert a[0] >= 0, "square root of a negative interval"
    return (math.isqrt(a[0] << SH), math.isqrt(a[1] << SH) + 1)


def f_and_df(eq, x):
    """Interval enclosures of f(x) = LHS(x) + x and of f'(x), for an interval x of fixed-point numbers."""
    c = lambda q: iv_const(q)
    u = iv_sqrt(iv_add(c(eq["p1"]), iv_mul(c(eq["p2"]), x)))
    v = iv_sqrt(iv_add(c(eq["q1"]), iv_mul(c(eq["q2"]), x)))
    w = iv_add(u, iv_neg(iv_mul(c(eq["k"]), v)))                    # u - k v
    f = iv_add(iv_add(iv_mul(c(eq["a1"] + 1), x), iv_neg(iv_mul(c(eq["c"]), iv_mul(w, w)))), c(eq["a0"]))
    du = iv_div(c(Q(eq["p2"], 2)), u)                               # u' = p2/(2u)
    dv = iv_div(c(Q(eq["q2"], 2)), v)                               # v' = q2/(2v)
    dw = iv_add(du, iv_neg(iv_mul(c(eq["k"]), dv)))
    df = iv_add(c(eq["a1"] + 1), iv_neg(iv_mul(c(2 * eq["c"]), iv_mul(w, dw))))
    return f, df


def route_b(eq, interval, steps=330):
    lo, hi = fx_floor(interval[0]), fx_ceil(interval[1])
    _, df = f_and_df(eq, (lo, hi))
    monotone = df[0] > 0
    flo, _ = f_and_df(eq, (lo, lo))
    fhi, _ = f_and_df(eq, (hi, hi))
    sign_change = flo[1] < 0 and fhi[0] > 0
    if not (monotone and sign_change):
        return dict(monotone=monotone, sign_change=sign_change, df=df, enclosure=None)
    a, b = lo, hi
    for _ in range(steps):
        m = (a + b) // 2
        fm, _ = f_and_df(eq, (m, m))
        if fm[1] < 0:
            a = m
        elif fm[0] > 0:
            b = m
        else:
            break
    return dict(monotone=monotone, sign_change=sign_change, df=df, enclosure=(a, b))


def fx_str(v, digits=40):
    sgn = "-" if v < 0 else ""
    v = abs(v)
    ip = v >> SH
    frac = ((v - (ip << SH)) * 10 ** digits) >> SH
    return f"{sgn}{ip}.{frac:0{digits}d}"


def qsqrt_dec(r, s, d, digits=50):
    """r + s sqrt(d), rounded down to `digits` decimals (exact integer arithmetic)."""
    scale = 10 ** digits
    # floor((r + s sqrt d) * scale) via isqrt on s^2 d scale^2
    num = Q(r) * scale
    S = Q(s)
    t = math.isqrt((S.numerator ** 2 * d * scale * scale) // (S.denominator ** 2))
    val = num + (t if S > 0 else -t - 1)
    v = math.floor(val)
    sgn = "-" if v < 0 else ""
    v = abs(v)
    return f"{sgn}{v // scale}.{v % scale:0{digits}d}"


def compare_printed(r, s, d, printed, mult):
    """Sign of printed - mult*theta*, where theta* = -(r + s sqrt d); exact."""
    p = Q(printed)
    return cmp_qsqrt(p + mult * r, mult * s, d)


def run(eq, interval, printed, label, verbose=True):
    out = []
    say = (lambda *a: out.append(" ".join(str(x) for x in a)))
    (ia, ib, ic), (A, B), sols, rejected = route_a(eq, interval)
    say(f"[{label}] route A: (1.4) <=> uv = {A} x + {B}; squared: {ia} x^2 + {ib} x + {ic} = 0")
    for rec in rejected:
        say(f"[{label}] route A: root {rec['r']} + {rec['s']}*sqrt({rec['d']}) rejected "
            f"(in interval {rec['in_interval']}, uv >= 0 {rec['uv_nonneg']}, radicands >= 0 {rec['domain']})")
    ok = len(sols) == 1
    if not ok:
        say(f"[{label}] route A: {len(sols)} admissible solutions in the interval; expected exactly one")
        return False, out, None
    x = sols[0]
    # exact residual: u^2 v^2 - (A x + B)^2 = 0 in Q(sqrt d), and uv >= 0 was checked
    u2v2 = check_residual(eq, x["r"], x["s"], x["d"])
    axb = (A * x["r"] + B, A * x["s"])
    sq = (axb[0] ** 2 + axb[1] ** 2 * x["d"], 2 * axb[0] * axb[1])
    resid_zero = (u2v2[0] - sq[0] == 0) and (u2v2[1] - sq[1] == 0)
    say(f"[{label}] route A: x = {x['r']} + ({x['s']})*sqrt({x['d']}), exact residual zero: {resid_zero}")
    theta = (-x["r"], -x["s"], x["d"])
    say(f"[{label}] theta*   = {qsqrt_dec(*theta)}...")
    say(f"[{label}] 2 theta* = {qsqrt_dec(2 * theta[0], 2 * theta[1], theta[2])}...")
    rb = route_b(eq, interval)
    say(f"[{label}] route B: f' on the whole interval in [{fx_str(rb['df'][0], 6)}, {fx_str(rb['df'][1], 6)}] "
        f"(> 0: {rb['monotone']}); sign change at the ends: {rb['sign_change']}")
    if rb["enclosure"] is None:
        say(f"[{label}] route B: no enclosure")
        return False, out, theta
    a, b = rb["enclosure"]
    contains = cmp_qsqrt(x["r"] - Q(a, ONE), x["s"], x["d"]) >= 0 and cmp_qsqrt(x["r"] - Q(b, ONE), x["s"], x["d"]) <= 0
    say(f"[{label}] route B: root in [{fx_str(a)}, {fx_str(b)}] (width 2^{(b - a).bit_length() - SH}); "
        f"contains route A's root: {contains}")
    ok = ok and resid_zero and rb["monotone"] and rb["sign_change"] and contains
    for row, (dec, mult, side) in printed.items():
        sgn = compare_printed(theta[0] * -1, theta[1] * -1, theta[2], dec, mult)  # printed - mult*theta*
        holds = sgn >= 0 if side == "upper" else sgn <= 0
        rel = {1: "above", 0: "equal to", -1: "below"}[sgn]
        say(f"[{label}] {row}: printed {dec} is {rel} the exact value {mult} theta*: "
            f"{'a valid' if holds else 'NOT a valid'} {side} bound as printed")
    return ok, out, theta


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if len(args) != 1:
        raise SystemExit(__doc__)
    h, tex = read_equation(args[0])
    print(f"main.tex sha256 {h}")
    same = tex == re.sub(r"\s+", "", TRANSCRIPTION)
    print(f"equation (1.4) in main.tex: {tex}")
    print(f"our transcription equals it: {same}")
    ok, out, theta = run(EQ, INTERVAL, PRINTED, "paper")
    print("\n".join(out))
    findings = [l for l in out if "NOT a valid" in l]
    status = ok and same
    if "--controls" in sys.argv:
        print("\nred controls (each must be refused or change the verdict):")
        ctl_ok = True
        # 1. the constant 51/200 changed to 52/200: the root moves, the printed 63a value no longer bounds it
        e = dict(EQ, a0=Q(52, 200))
        ok1, out1, th1 = run(e, INTERVAL, {"63a printed": ("0.3144831759741", 1, "upper")}, "ctl 51/200 -> 52/200")
        print("\n".join(out1))
        ctl_ok &= any("NOT a valid" in l for l in out1)
        # 2. the printed 63a value lowered by one unit in its last digit
        ok2, out2, _ = run(EQ, INTERVAL, {"63a lowered": ("0.3144831759740", 1, "upper")}, "ctl 63a - 1e-13")
        print("\n".join(out2[-1:]))
        ctl_ok &= "NOT a valid" in out2[-1]
        # 3. the interval moved to [-0.2, -0.1], which holds only the spurious root of the squared equation
        ok3, out3, _ = run(EQ, (Q(-2, 10), Q(-1, 10)), {}, "ctl interval [-0.2,-0.1]")
        print("\n".join(out3[:3]))
        ctl_ok &= not ok3
        # 4. the sign inside the square: sqrt(2(1-14x)) + 5 sqrt(-1-8x) (k = -5): route A and B must not
        #    reproduce the paper's root
        e4 = dict(EQ, k=-5)
        ok4, out4, th4 = run(e4, INTERVAL, {}, "ctl k = -5")
        print("\n".join(out4[:4]))
        ctl_ok &= (not ok4) or th4 != theta
        print(f"red controls all as expected: {ctl_ok}")
        status = status and ctl_ok
    print(f"RESULT: {'PASS' if status else 'FAIL'} (route A and route B agree; exactly one solution on [-0.35, -0.3])"
          + (f"; FINDING: {len(findings)} printed value(s) below the exact value" if findings else ""))
    sys.exit(0 if status else 1)


if __name__ == "__main__":
    main()
