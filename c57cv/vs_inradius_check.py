#!/usr/bin/env python3
"""Replay of the inradius half of the 57c row [vS2026], U <= 0.65628589* (W. van Staden, the univalent Bloch-Landau
constant; Preprints.org 202610.1936 and the repository WynandVStaden/bloch-landau-slit-discs @ a6f17b5), written from
the paper's LaTeX source (paper/bloch_slit_discs.tex: Section 2, the domain, and Proposition 3.1 with Lemmas 3.2 and 3.3)
and the data file certify/inradius_leaves_G.npz alone. inradius_certificate.py and every other program of the
repository were neither read nor run; the .npz file is read by our own parser (zip members in the .npy format: a fixed
header, then little-endian int64 or float64 values).

The claim: R(D) <= lambda = 1 + 10^-9, where D is the component containing 0 of {|z| < R} minus 24 slits,
R = 5.001117628756668, a slit with parameters (theta, t, b1, b2) being z = R e^s exp(i (theta + b1 s^2 + b2 s^4)),
-a <= s <= 0, a = log(R/t) (the paper's x = s/a), and the 24 slits the images of the five classes of Table 1 under
z -> e^{2 pi i/3} z and conjugation (theta -> -theta, b1 -> -b1, b2 -> -b2).

The argument (Proposition 3.1, re-derived): D is invariant under the two symmetries, so a disc of radius r in D can be
moved to one with centre p in S = {0 <= arg p <= pi/3, |p| < R} or p = 0; S lies in [0, 8]^2. If dyadic cells tile
[0, 8]^2 and each cell is (a) outside the open disc |p| < R, (b) inside arg p > pi/3, or (c) within lambda of one point q
outside D (a point of a slit or of the circle |z| = R) at its four corners (then at every point, by convexity), then
r <= |p - q| <= lambda.

What the data hold (decoded by us from the paper's description and checked below): cells (l, i, j, kind, type, k):
the cell [i 2^{3-l}, (i+1) 2^{3-l}] x [j 2^{3-l}, (j+1) 2^{3-l}], kind 0, 1, 2 for (a), (b), (c); for kind (c) the
witness is type 0: the point of slit k at s = pars (a double, read exactly), type 1: the tip of slit k (s = -a), or
type 2: the circle point R e^{i pars}; slits are numbered class by class (class 0: 0-2, class 1: 3-5, class 2: 6-11,
class 3: 12-17, class 4: 18-23), and the witnesses only use the slit of each class with r = 0, epsilon = +1 (k = 0, 3, 6,
12, 18) and k = 1 (class 0 rotated by 2 pi/3). Under this reading every check passes; it is the reading, not trust in
it, that the checks test (a wrong reading would put q far from its cells).

Checks, in exact rational arithmetic with outward-rounded series (method A) and again in 50-digit decimal arithmetic with
directed rounding (method B, a separate implementation of exp, log, sin and cos):
  1. the cells tile [0, 8]^2: all inside, pairwise distinct, no cell inside another (dyadic squares are nested or have
     disjoint interiors), areas summing to 64 exactly;
  2. kind (a): x0^2 + y0^2 >= R^2; kind (b): y0^2 > 3 x1^2, exactly;
  3. kind (c): -a < s <= 0 for slit points (as R e^s > t), the witness q enclosed in a box, and the largest distance
     from each corner to that box at most lambda, exactly.
Red controls (each must be refused): lambda = 1 - 10^-6 (the disc |z| < 1 lies in D, so no certificate can give less
than 1); one cell removed (the areas no longer add to 64); the slit-point witness of the smallest cell moved by 10^-2
along its slit (10^-6 is not enough: the witness is the nearest point, so the distance moves only to second order); a kind (a) cell claimed for a cell that meets the disc.

Usage: python3 vs_inradius_check.py NPZ [--controls]     (standard library only; about a minute)
"""
import ast
import decimal
import hashlib
import struct
import sys
import time
import zipfile
from fractions import Fraction as F

SHA_NPZ = "512a1c84dfcca9a5cdce985c927bf976781af892e9511077e186d0888aa1e6af"
R = F("5.001117628756668")
LAM = 1 + F(1, 10 ** 9)
# Table 1 of the paper: class -> (theta, t, b1, b2); theta of class 1 is pi/3 (marked None, handled with pi bounds)
CLASSES = [(F(0), F(1), F(0), F(0)),
           (None, F("1.9999999974092333"), F(0), F(0)),
           (F("0.029930188727091613"), F("2.789342409403121"), F("0.9026921410641155"), F("-0.3301459494660587")),
           (F("1.0157094972537686"), F("3.473381248030611"), F("-1.2656124700253997"), F("0.9794269477754911")),
           (F("0.5757776442359243"), F("4.342874184977691"), F("0.03225521024585182"), F("2.9360081049550737"))]


def need(c, msg):
    if not c:
        raise AssertionError(msg)


# ----------------------------------------------------------------------------------------------- the data file
def load_npz(path):
    z = zipfile.ZipFile(path)
    out = {}
    for n in z.namelist():
        b = z.read(n)
        need(b[:6] == b"\x93NUMPY", "not an .npy member: " + n)
        if b[6] == 1:
            hl, off = int.from_bytes(b[8:10], "little"), 10
        else:
            hl, off = int.from_bytes(b[8:12], "little"), 12
        hdr = ast.literal_eval(b[off:off + hl].decode("latin1"))     # a dict literal; literal_eval runs no code
        need(not hdr["fortran_order"], "fortran order")
        cnt = 1
        for s_ in hdr["shape"]:
            cnt *= s_
        data = b[off + hl:]
        need(len(data) == 8 * cnt, "size")
        fmt = {"<i8": "q", "<f8": "d"}[hdr["descr"]]
        out[n[:-4]] = (hdr["shape"], list(struct.unpack("<%d%s" % (cnt, fmt), data)))
    return out


# ----------------------------------------------------------------------------------------------- method A: rationals
P = 220
ONE = 1 << P


def rd(x):
    return F((x.numerator * ONE) // x.denominator, ONE)


def ru(x):
    return F(-((-x.numerator * ONE) // x.denominator), ONE)


def atan_inv(n, terms=90):
    s, prev = F(0), F(0)
    for k in range(terms + 1):
        prev = s
        s += F((-1) ** k, (2 * k + 1) * n ** (2 * k + 1))
    return min(s, prev), max(s, prev)


A5, A239 = atan_inv(5), atan_inv(239)
PI = (rd(16 * A5[0] - 4 * A239[1]), ru(16 * A5[1] - 4 * A239[0]))       # Machin


def exp_a(x, terms=60):
    """exp(x) for rational |x| <= 3: Taylor partial sum and Lagrange remainder |x|^(N+1) e^|x| / (N+1)! (e^|x| < 21)."""
    need(abs(x) <= 3, "exp range")
    s, term = F(0), F(1)
    for k in range(terms + 1):
        s += term
        term = term * x / (k + 1)
    rem = abs(term) * 21
    return rd(s - rem), ru(s + rem)


def sincos_a(phi_lo, phi_hi, terms=70):
    """sin and cos on [phi_lo, phi_hi] (|phi| <= 9): Taylor at phi_lo with Lagrange remainder, widened by the width."""
    x = phi_lo
    need(abs(x) <= 9 and phi_hi - phi_lo < F(1, 10 ** 20), "sincos range")
    s, c, term = F(0), F(0), F(1)
    for k in range(terms + 1):
        if k % 4 == 0:
            c += term
        elif k % 4 == 1:
            s += term
        elif k % 4 == 2:
            c -= term
        else:
            s -= term
        term = term * x / (k + 1)
    rem = abs(term) + (phi_hi - phi_lo)
    return (rd(s - rem), ru(s + rem)), (rd(c - rem), ru(c + rem))


def log_a(y, terms=120):
    """log(y) for rational 1 < y <= 6: 2 artanh(u), u = (y-1)/(y+1) <= 5/7, tail u^(2N+3)/((2N+3)(1-u^2))."""
    u = (y - 1) / (y + 1)
    s, p = F(0), u
    for k in range(terms + 1):
        s += p / (2 * k + 1)
        p *= u * u
    tail = p / ((2 * terms + 3) * (1 - u * u))
    return rd(2 * s), ru(2 * (s + tail))


def slit(k):
    """(theta interval, t, b1, b2) of slit k: class k // 3 for k < 6, else class 2 + (k - 6) // 6; within a class of
    six, slit 6 c' + m has r = m // 2 and epsilon = +1 for even m (only m = 0 is used by the data)."""
    if k < 6:
        cl, r, eps = k // 3, k % 3, 1
    else:
        cl, m = 2 + (k - 6) // 6, (k - 6) % 6
        r, eps = m // 2, (1 if m % 2 == 0 else -1)
    th, t, b1, b2 = CLASSES[cl]
    th_iv = (PI[0] / 3, PI[1] / 3) if th is None else (th, th)
    th_iv = (eps * th_iv[0], eps * th_iv[1]) if eps > 0 else (eps * th_iv[1], eps * th_iv[0])
    rot = (2 * PI[0] * r / 3, 2 * PI[1] * r / 3)
    return (th_iv[0] + rot[0], th_iv[1] + rot[1]), t, eps * b1, eps * b2


def witness_a(typ, k, par):
    """an enclosure (x_lo, x_hi, y_lo, y_hi) of the witness point q, and whether its parameter is admissible."""
    if typ == 2:                                   # circle point R e^{i par}
        (s_lo, s_hi), (c_lo, c_hi) = sincos_a(par, par)
        return (rd(R * c_lo), ru(R * c_hi), rd(R * s_lo), ru(R * s_hi)), True
    th, t, b1, b2 = slit(k)
    if typ == 0:                                   # slit point at s = par; admissible iff -a < s <= 0, i.e. R e^s > t
        s = par
        e_lo, e_hi = exp_a(s)
        ok = s <= 0 and R * e_lo > t
        rad = (R * e_lo, R * e_hi)
        phi = (th[0] + b1 * s * s + b2 * s ** 4, th[1] + b1 * s * s + b2 * s ** 4)
    else:                                          # the tip: s = -a, a = log(R/t), modulus t
        a_lo, a_hi = log_a(R / t)
        ok = True
        rad = (t, t)
        # phi(-a) = th + b1 a^2 + b2 a^4, enclosed over a in [a_lo, a_hi] (a > 0, so a^2 and a^4 are monotone)
        q2 = [b1 * a_lo ** 2, b1 * a_hi ** 2]
        q4 = [b2 * a_lo ** 4, b2 * a_hi ** 4]
        phi = (th[0] + min(q2) + min(q4), th[1] + max(q2) + max(q4))
    (s_lo, s_hi), (c_lo, c_hi) = sincos_a(phi[0], phi[1])
    xs = [rad[0] * c_lo, rad[0] * c_hi, rad[1] * c_lo, rad[1] * c_hi]
    ys = [rad[0] * s_lo, rad[0] * s_hi, rad[1] * s_lo, rad[1] * s_hi]
    return (rd(min(xs)), ru(max(xs)), rd(min(ys)), ru(max(ys))), ok


def corners_within(cell, box, lam):
    l, i, j = cell
    h = F(8, 2 ** l)
    lam2 = lam * lam
    for X in (i * h, (i + 1) * h):
        for Y in (j * h, (j + 1) * h):
            dx = max(abs(X - box[0]), abs(X - box[1]))
            dy = max(abs(Y - box[2]), abs(Y - box[3]))
            if dx * dx + dy * dy > lam2:
                return False
    return True


# ----------------------------------------------------------------------------------------------- method B: decimal
def method_b(rows, pars, lam_str="1.000000001"):
    D = decimal.Decimal
    prec = 50
    cf = decimal.Context(prec=prec, rounding=decimal.ROUND_FLOOR)
    cc = decimal.Context(prec=prec, rounding=decimal.ROUND_CEILING)
    ch = decimal.Context(prec=prec, rounding=decimal.ROUND_HALF_EVEN)
    tiny = D(10) ** -40
    # pi from Gauss's formula in integers at 10^60 (the formula: (18+i)^12 (57+i)^8 (239-i)^5 is a positive multiple of 1+i)
    sc = 10 ** 60

    def at(n, terms=60):
        lo = hi = 0
        sums = []
        for k in range(terms + 1):
            d_ = (2 * k + 1) * n ** (2 * k + 1)
            tl, th_ = sc // d_, -((-sc) // d_)
            if k % 2 == 0:
                lo, hi = lo + tl, hi + th_
            else:
                lo, hi = lo - th_, hi - tl
            sums.append((lo, hi))
        return min(sums[-1][0], sums[-2][0]), max(sums[-1][1], sums[-2][1])
    a18, a57, a239 = at(18), at(57), at(239)
    pi_lo = cf.scaleb(D(48 * a18[0] + 32 * a57[0] - 20 * a239[1]), -60)
    pi_hi = cc.scaleb(D(48 * a18[1] + 32 * a57[1] - 20 * a239[0] + 1), -60)
    Rd = D("5.001117628756668")
    lam = D(lam_str)
    classes = [(D(0), D(1), D(0), D(0)), (None, D("1.9999999974092333"), D(0), D(0)),
               (D("0.029930188727091613"), D("2.789342409403121"), D("0.9026921410641155"), D("-0.3301459494660587")),
               (D("1.0157094972537686"), D("3.473381248030611"), D("-1.2656124700253997"), D("0.9794269477754911")),
               (D("0.5757776442359243"), D("4.342874184977691"), D("0.03225521024585182"), D("2.9360081049550737"))]

    def sincos(x, w):
        """Taylor at x with |x| <= 9, Lagrange remainder; each partial sum accumulated twice (floor and ceiling)."""
        s_lo = s_hi = c_lo = c_hi = D(0)
        t_lo = t_hi = D(1)
        for k in range(80):
            m = k % 4
            if m == 0:
                c_lo, c_hi = cf.add(c_lo, t_lo), cc.add(c_hi, t_hi)
            elif m == 1:
                s_lo, s_hi = cf.add(s_lo, t_lo), cc.add(s_hi, t_hi)
            elif m == 2:
                c_lo, c_hi = cf.subtract(c_lo, t_hi), cc.subtract(c_hi, t_lo)
            else:
                s_lo, s_hi = cf.subtract(s_lo, t_hi), cc.subtract(s_hi, t_lo)
            # next term x^(k+1)/(k+1)!, enclosed (x may be negative: track the exact sign by products of both ends)
            cand = [cf.divide(cf.multiply(t_lo, x), D(k + 1)), cc.divide(cc.multiply(t_lo, x), D(k + 1)),
                    cf.divide(cf.multiply(t_hi, x), D(k + 1)), cc.divide(cc.multiply(t_hi, x), D(k + 1))]
            t_lo, t_hi = min(cand), max(cand)
        rem = cc.add(max(abs(t_lo), abs(t_hi)), w)
        return (cf.subtract(s_lo, rem), cc.add(s_hi, rem)), (cf.subtract(c_lo, rem), cc.add(c_hi, rem))

    def exp_iv(x):
        v = ch.exp(x)
        return cf.subtract(v, cc.multiply(tiny, v)), cc.add(v, cc.multiply(tiny, v))

    def ln_iv(x):
        v = ch.ln(x)
        e = cc.multiply(tiny, abs(v) + 1)
        return cf.subtract(v, e), cc.add(v, e)

    def slit_b(k):
        if k < 6:
            cl, r, eps = k // 3, k % 3, 1
        else:
            cl, m = 2 + (k - 6) // 6, (k - 6) % 6
            r, eps = m // 2, (1 if m % 2 == 0 else -1)
        th, t, b1, b2 = classes[cl]
        th_iv = (cf.divide(pi_lo, 3), cc.divide(pi_hi, 3)) if th is None else (th, th)
        if eps < 0:
            th_iv = (-th_iv[1], -th_iv[0])
        rot = (cf.divide(cf.multiply(2 * r, pi_lo), 3), cc.divide(cc.multiply(2 * r, pi_hi), 3))
        return (cf.add(th_iv[0], rot[0]), cc.add(th_iv[1], rot[1])), t, eps * b1, eps * b2

    def witness(typ, k, par):
        if typ == 2:
            (sl, sh), (cl_, chh) = sincos(par, D(0))
            return (cf.multiply(Rd, cl_), cc.multiply(Rd, chh), cf.multiply(Rd, sl), cc.multiply(Rd, sh)), True
        th, t, b1, b2 = slit_b(k)
        if typ == 0:
            s = par
            e_lo, e_hi = exp_iv(s)
            ok = s <= 0 and cf.multiply(Rd, e_lo) > t
            rad = (cf.multiply(Rd, e_lo), cc.multiply(Rd, e_hi))
            # b1 s^2 + b2 s^4 at 50 digits: every product is within one unit of the 50th digit, far below the
            # widening by tiny (10^-40) on both sides
            q = ch.add(ch.multiply(b1, ch.multiply(s, s)), ch.multiply(b2, ch.power(s, 4)))
            phi = (cf.subtract(cf.add(th[0], q), tiny), cc.add(cc.add(th[1], q), tiny))
        else:
            a_lo, a_hi = ln_iv(cc.divide(Rd, t))
            a_lo = cf.subtract(a_lo, tiny)
            ok = True
            rad = (t, t)
            q2 = [cf.multiply(b1, cf.multiply(a_lo, a_lo)), cc.multiply(b1, cc.multiply(a_hi, a_hi)),
                  cc.multiply(b1, cc.multiply(a_lo, a_lo)), cf.multiply(b1, cf.multiply(a_hi, a_hi))]
            q4 = [cf.multiply(b2, cf.power(a_lo, 4)), cc.multiply(b2, cc.power(a_hi, 4)),
                  cc.multiply(b2, cc.power(a_lo, 4)), cf.multiply(b2, cf.power(a_hi, 4))]
            phi = (cf.subtract(cf.add(cf.add(th[0], min(q2)), min(q4)), tiny),
                   cc.add(cc.add(cc.add(th[1], max(q2)), max(q4)), tiny))
        w = cc.subtract(phi[1], phi[0])
        (sl, sh), (cl_, chh) = sincos(phi[0], w)
        xs = [cf.multiply(rad[0], cl_), cc.multiply(rad[0], chh), cf.multiply(rad[1], cl_), cc.multiply(rad[1], chh),
              cc.multiply(rad[0], cl_), cf.multiply(rad[0], chh), cc.multiply(rad[1], cl_), cf.multiply(rad[1], chh)]
        ys = [cf.multiply(rad[0], sl), cc.multiply(rad[0], sh), cf.multiply(rad[1], sl), cc.multiply(rad[1], sh),
              cc.multiply(rad[0], sl), cf.multiply(rad[0], sh), cc.multiply(rad[1], sl), cf.multiply(rad[1], sh)]
        return (min(xs), max(xs), min(ys), max(ys)), ok

    cache = {}
    bad = 0
    lam2 = cf.multiply(lam, lam)
    for r_, p in zip(rows, pars):
        if r_[3] != 2:
            continue
        key = (r_[4], r_[5], p)
        if key not in cache:
            cache[key] = witness(r_[4], r_[5], D(p) if r_[4] != 1 else D(0))
        box, ok = cache[key]
        l, i, j = r_[0], r_[1], r_[2]
        h = ch.divide(D(8), ch.power(D(2), l))      # exact: 2^(3-l) has at most 23 significant digits here
        good = ok
        for X in (ch.multiply(D(i), h), ch.multiply(D(i + 1), h)):        # exact at 50 digits (at most 33 digits)
            for Y in (ch.multiply(D(j), h), ch.multiply(D(j + 1), h)):
                dx = max(abs(cc.subtract(X, box[0])), abs(cc.subtract(X, box[1])),
                         abs(cf.subtract(X, box[0])), abs(cf.subtract(X, box[1])))
                dy = max(abs(cc.subtract(Y, box[2])), abs(cc.subtract(Y, box[3])),
                         abs(cf.subtract(Y, box[2])), abs(cf.subtract(Y, box[3])))
                if cc.add(cc.multiply(dx, dx), cc.multiply(dy, dy)) > lam2:
                    good = False
        bad += not good
    return bad, len(cache)


# ----------------------------------------------------------------------------------------------- the checks
def check(rows, pars, lam, log=print, skip_tiling=False):
    t0 = time.time()
    ok = True
    if not skip_tiling:
        cells = set()
        inside = True
        for r_ in rows:
            l, i, j = r_[0], r_[1], r_[2]
            inside &= 0 <= i < 2 ** l and 0 <= j < 2 ** l and l >= 6
            cells.add((l, i, j))
        distinct = len(cells) == len(rows)
        nested = sum(1 for (l, i, j) in cells for l2 in range(0, l) if (l2, i >> (l - l2), j >> (l - l2)) in cells)
        area = sum(F(64, 4 ** l) for (l, i, j) in cells)
        tiles = inside and distinct and nested == 0 and area == 64
        log("  1. tiling: %d cells, inside [0,8]^2 with l >= 6: %s; pairwise distinct: %s; cells inside another: %d; "
            "total area %s: %s" % (len(rows), inside, distinct, nested, area, tiles))
        ok &= tiles
    na = nb = nc = fa = fb = fc = 0
    cache = {}
    for r_, p in zip(rows, pars):
        l, i, j, kind, typ, k = r_
        h = F(8, 2 ** l)
        x0, y0, x1 = i * h, j * h, (i + 1) * h
        if kind == 0:
            na += 1
            fa += not (x0 * x0 + y0 * y0 >= R * R)
        elif kind == 1:
            nb += 1
            fb += not (y0 * y0 > 3 * x1 * x1)
        else:
            nc += 1
            key = (typ, k, p)
            if key not in cache:
                cache[key] = witness_a(typ, k, F(p) if typ != 1 else F(0))
            box, adm = cache[key]
            fc += not (adm and corners_within((l, i, j), box, lam))
    log("  2. kind (a): %d cells, x0^2 + y0^2 >= R^2 fails on %d; kind (b): %d cells, y0^2 > 3 x1^2 fails on %d"
        % (na, fa, nb, fb))
    log("  3. kind (c): %d cells with %d distinct witnesses; parameter admissible and all four corners within "
        "lambda = %s of the witness box: fails on %d (%.0f s)" % (nc, len(cache), str(lam) if lam.denominator > 10 ** 12 else float(lam), fc, time.time() - t0))
    return ok and fa == 0 and fb == 0 and fc == 0, (na, nb, nc, len(cache))


def main():
    t0 = time.time()
    path = [a for a in sys.argv[1:] if not a.startswith("--")][0]
    controls = "--controls" in sys.argv
    sha = hashlib.sha256(open(path, "rb").read()).hexdigest()
    print("57c [vS2026], the inradius: R(D) <= 1 + 10^-9 (Proposition 3.1), from certify/inradius_leaves_G.npz")
    print("  sha256 %s: %s" % (sha, "as pinned" if sha == SHA_NPZ else "DIFFERENT"))
    d = load_npz(path)
    shape, cells = d["cells"]
    rows = [tuple(cells[6 * n:6 * n + 6]) for n in range(shape[0])]
    pars = d["pars"][1]
    eps = d["eps"][1][0]
    print("  %d cells; eps = %r (the file's; we use lambda = 1 + 10^-9 exactly)" % (len(rows), eps))
    print("  pi (Machin) enclosed to width %.1e" % float(PI[1] - PI[0]))
    print("Method A (exact rationals, outward-rounded series):")
    okA, counts = check(rows, pars, LAM)
    print("  counts (a, b, c, distinct witnesses) = %s; the paper states 2799, 464, 1730 and 610" % (counts,))
    okA &= counts == (2799, 464, 1730, 610)
    print("Method B (50-digit decimal, directed rounding, Gauss's pi, Decimal exp and ln widened):")
    badB, nw = method_b(rows, pars)
    print("  kind (c): %d witnesses; cells failing: %d" % (nw, badB))
    ok = okA and badB == 0 and sha == SHA_NPZ
    print("  so every disc in D has radius at most 1 + 10^-9: %s" % ok)
    if controls:
        print("Red controls (each must be refused):")
        okc, _ = check(rows, pars, 1 - F(1, 10 ** 6), log=lambda s: None, skip_tiling=True)
        print("  1. lambda = 1 - 10^-6: %s" % ("refused" if not okc else "ACCEPTED (bad)"))
        ok &= not okc
        okc, _ = check(rows[1:], pars[1:], LAM, log=lambda s: None)
        print("  2. the first cell removed: %s" % ("refused (the areas add to less than 64)" if not okc else "ACCEPTED (bad)"))
        ok &= not okc
        # the slit-point witness of the smallest cell (where the distance to the boundary is nearly 1), moved by 10^-2
        # along its slit (it is the nearest point, so the distance changes only to second order: 10^-6 is not enough)
        n = max((m for m, r_ in enumerate(rows) if r_[3] == 2 and r_[4] == 0), key=lambda m: rows[m][0])
        p2 = [p + 1e-2 if (rows[m][3] == 2 and rows[m][4] == 0 and p == pars[n] and rows[m][5] == rows[n][5]) else p
              for m, p in enumerate(pars)]
        okc, _ = check(rows, p2, LAM, log=lambda s: None, skip_tiling=True)
        print("  3. the slit-point witness of the smallest such cell (side 2^-%d) moved by 10^-2 along its slit: %s"
              % (rows[n][0] - 3, "refused" if not okc else "ACCEPTED (bad)"))
        ok &= not okc
        m = next(m for m, r_ in enumerate(rows) if r_[3] == 2)
        rows2 = list(rows)
        rows2[m] = rows[m][:3] + (0,) + rows[m][4:]
        okc, _ = check(rows2, pars, LAM, log=lambda s: None, skip_tiling=True)
        print("  4. a kind (c) cell relabelled kind (a): %s" % ("refused" if not okc else "ACCEPTED (bad)"))
        ok &= not okc
    print("RESULT: %s" % ("PASS" if ok else "FAIL"))
    print("run time %.0f s" % (time.time() - t0))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
