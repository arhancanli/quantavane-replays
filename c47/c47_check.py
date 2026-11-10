#!/usr/bin/env python3
"""Exact replay of periodic-measure lower bounds for C_47 = c_2 (47a page: the weak type (1,1) constant of the
centred maximal operator over axis-parallel squares in R^2).

A row is a measure mu, periodic under Px Z x Py Z with finitely many atoms (x, y, m) in the cell C = [0,Px) x [0,Py),
and a level lam (= 1 in every row). Written from the definitions and Tomar's Lemma 6.2 alone:

  Lemma 6.1/6.2. If S is a finite set of atoms and 0 < rho, 4 lam rho^2 < mu(S), every point x with |p - x|_inf <= rho
  for all p in S has M mu(x) >= lam; and if U is the union of the lattice translates of finitely many rectangles
  R_i, each inside R_{S_i}(rho_i) = intersection of the closed squares Q(p, rho_i), p in S_i, with such rho_i, then
  c_2 >= lam |U n C| / mu(C).  (Re-derived here: truncate mu to N x N cells, smear each atom m at p into
  m (2 eps)^-2 1_{Q(p,eps)}, apply the weak type inequality and let N -> infinity.)

The rectangles are those of Tomar's Section 6: for every box [x_lo,x_hi] x [y_lo,y_hi] whose sides lie on atom
coordinates and have length at most `reach`, S = the atoms in the box, rho_S = (k - 1) 2^-e with k the largest integer
having (k 2^-e)^2 <= mu(S)/(4 lam) (e = 40), and the rectangle [x_hi - rho, x_lo + rho] x [y_hi - rho, y_lo + rho].

  1. Every box with its lower-left corner in C is enumerated (masses by exact 2D prefix sums over the periodic grid of
     atom coordinates); each rectangle is cut by the lattice into pieces inside C.
  2. Method A: the pieces are subtracted one by one from C, exactly, leaving the uncovered pieces; the pieces that
     removed area form the contributing set K.
  3. Method B: |union K| by a sweep over x with the y-intervals of each slab merged (a different algorithm).
  4. Every uncovered piece is checked against every piece of the whole family (no overlap of positive area), so
     |union| <= |C| - |uncovered|; method B gives |union| >= |union K|. The two must be equal, which pins |U n C|.
  5. For every rectangle of K: its atom set is summed again directly from the atoms (not from prefix sums), and
     0 < rho, 4 lam rho^2 < mu(S), and piece - (lattice shift) inside R_S(rho) are checked exactly.
  6. lam |U n C| / mu(C) is compared with the stated fraction and the printed bound.

With --full the family is every box up to the side L* beyond which no box can give a rectangle of positive area
(mu(box) <= (L/Px + 2)(L/Py + 2) mu(C) <= L^2 for L >= L*), with Tomar's radii and with rigorous lower and upper
radii on a 2^-100 grid, which encloses |{M mu >= lam} n C| itself.

Usage:
  python3 -I c47_check.py selftest                  Lin's Phi, single atoms, Tomar's bare lattice, brute force
  python3 -I c47_check.py ROW [CERTS] [--controls] [--full]
  python3 -I c47_check.py measure FILE.json         evaluate any measure file (Tomar's input format)
Rows: tom2026, pr321, pr341; and lin2026 [--controls] (Lin's row, built from its closed forms, no file).
Python 3.9+, standard library only.
"""
import hashlib
import json
import os
import random
import sys
import time
from fractions import Fraction
from math import gcd, isqrt

HERE = os.path.dirname(os.path.abspath(__file__))
E_TOMAR = 40  # radii on the 2^-40 grid (Tomar, Section 6)

ROWS = {
    # [Tom2026], Theorem 1.3: c2lower_1p70272.json from Zenodo 10.5281/zenodo.23716078 (= GitHub v5 @ 3af09fc)
    "tom2026": dict(
        file=("tom2026/c2lower_1p70272.json", "e2ca1616be95afa6680e3ccf157cc16fdd39026ad0f5fa79ebd32e5214779888"),
        claimed=None, decimal="1.7027237335962937",  # the decimals PR #321 states for Tomar's exact ratio
        printed=["1.7027237335"], pr="Tomar, Zenodo 10.5281/zenodo.23716077, Theorem 1.3 (PRs #321 and #341)"),
    # PR #321 (rpinho): the measure half59-measure.json and the certificate half59.json @ 3d4cd55
    "pr321": dict(
        file=("pr321/half59-measure.json", "64585e1f89c428711dd1c36197de26326c23ac1c83e6aa132dbaec0a2dc54b16"),
        witnesses=("pr321/half59.json", "b4d88847ef0c41577213a34ef6fda3d0d9b52b672f6ce20e998b006a2031ba99"),
        claimed=Fraction(215956002657486428307163448091951082893881, 126800499112387565339223240212480000000000),
        printed=["1.70311", "1.70311634"], pr="#321"),
    # PR #350 (rpinho, stacked on #341): certificates/final.json @ 7273e97 (98 atoms and 1500 witness rectangles in one
    # file); its radii are on the 2^-48 grid (D = 2^48), so our own family uses that grid too
    "pr350": dict(
        file=("pr350/final.json", "bc342d18f9fe219828856bcaa9531669a1051c2685394e41dfbeb3e794e818eb"),
        witnesses=("pr350/final.json", "bc342d18f9fe219828856bcaa9531669a1051c2685394e41dfbeb3e794e818eb"),
        claimed=Fraction(1771220457083239070854869084751983500209800598673,
                         1039431510074799032049552767561236480000000000000),
        printed=["1.7040280575"], pr="#350", e=48, witness_only=True,
        finite=(214084879817, 2, Fraction(81179170702072811560461937774649547472268267876630428381100034109236137,
                                          47639573975778160894755636560089704255253920609263650078720000000000000))),
    # PR #350 (rpinho), update of 10 Nov 2026 (round 8, 14:13 UTC): certificates/final.json @ 85c4680 (440 atoms and
    # 7806 witness rectangles in one file, radii on the 2^-64 grid); the finite instance of release_validation/finite.json
    "pr350r8": dict(
        file=("pr350r8/final.json", "1afaf5d5c8968fffd3f98cf77de4d462dd1070b9dea4ec1456dd5a91087bfb44"),
        witnesses=("pr350r8/final.json", "1afaf5d5c8968fffd3f98cf77de4d462dd1070b9dea4ec1456dd5a91087bfb44"),
        claimed=Fraction(114388018195906255044427129545533845966664904690010061034899,
                         67069990589546886330812144279013633721316147200000000000000),
        printed=["1.7055022252"], pr="#350 (round 8)", e=64, witness_only=True,
        finite=(708844255899, 2,
                Fraction(2299016964454535926164008318906010429285009516613200633736490794747466283233740859,
                         1347999979410719285506487968385821088532379413340370461131840304447488000000000000))),
    # PR #341 (Doldrums, SCiEL): construction.json in Zenodo 10.5281/zenodo.24022791, which is restricted
    "pr341": dict(
        file=("pr341/construction.json", "5f7a8d9dfc88514c24995fd81b085b7c05ab4e2b9e43dcb923e48b3b6a7e310b"),
        claimed=Fraction(3132581668909806427917067970204393210052649, 1839223215401194968772738622942534413844480),
        printed=["1.7032090736"], pr="#341"),
}


# ----------------------------------------------------------------------------------------------------------- input

def q(v):
    """An exact rational from a JSON string or integer; JSON floats are refused."""
    if isinstance(v, bool) or not isinstance(v, (int, str)):
        raise ValueError(f"not an exact number: {v!r}")
    return Fraction(v)


class Measure:
    def __init__(self, Px, Py, atoms, lam=Fraction(1), reach=None):
        self.Px, self.Py, self.lam = Fraction(Px), Fraction(Py), Fraction(lam)
        self.atoms = [(Fraction(x), Fraction(y), Fraction(m)) for x, y, m in atoms]
        self.reach = None if reach is None else Fraction(reach)
        if self.Px <= 0 or self.Py <= 0 or self.lam <= 0:
            raise ValueError("periods and level must be positive")
        seen = set()
        for x, y, m in self.atoms:
            if not (0 <= x < self.Px and 0 <= y < self.Py):
                raise ValueError(f"atom ({x}, {y}) is outside the cell [0,Px) x [0,Py)")
            if m <= 0:
                raise ValueError(f"atom ({x}, {y}) has mass {m} <= 0")
            if (x, y) in seen:
                raise ValueError(f"atom ({x}, {y}) is listed twice")
            seen.add((x, y))
        self.mass = sum(m for _, _, m in self.atoms)

    @staticmethod
    def from_json(d):
        return Measure(q(d["Px"]), q(d["Py"]), [(q(a), q(b), q(c)) for a, b, c in d["pts"]], q(d.get("lam", 1)),
                       q(d["reach"]) if "reach" in d else None)

    def copy(self, **kw):
        m = Measure(kw.get("Px", self.Px), kw.get("Py", self.Py), kw.get("atoms", self.atoms),
                    kw.get("lam", self.lam), kw.get("reach", self.reach))
        return m


def load_json(path, sha):
    raw = open(path, "rb").read()
    got = hashlib.sha256(raw).hexdigest()
    if sha is not None and got != sha:
        raise SystemExit(f"SHA-256 of {path} is {got}, expected {sha}")
    return json.loads(raw.decode("utf-8")), got


# ---------------------------------------------------------------------------------------------- density cut-off

def sqrt_up(x, bits=64):
    """A rational >= sqrt(x) within 2^-bits (x a non-negative Fraction)."""
    return Fraction(isqrt(x.numerator * 2 ** (2 * bits) // x.denominator) + 1, 2 ** bits)


def cutoff(meas):
    """L* with: every box with longer side L >= L* has mass <= L^2, so 2 rho <= L and its rectangle has no interior.
    A closed interval of length L meets at most floor(L/P) + 2 cells, so mu(box) <= (L/Px + 2)(L/Py + 2) mu(C)."""
    Px, Py, mu = meas.Px, meas.Py, meas.mass / meas.lam
    d = mu / (Px * Py)
    if d >= 1:
        raise ValueError("density mu(C)/(lam Px Py) >= 1: no cut-off")
    A, B, C = 1 - d, 2 * d * (Px + Py), 4 * d * Px * Py
    L = (B + sqrt_up(B * B + 4 * A * C)) / (2 * A)
    assert A * L * L - B * L - C >= 0 and L >= B / (2 * A)  # the quadratic is >= 0 and increasing from L on
    return L


# --------------------------------------------------------------------------------------- the family of rectangles

def lcm_den(values):
    D = 1
    for v in values:
        D = D // gcd(D, v.denominator) * v.denominator
    return D


class Family:
    """All boxes with lower-left corner in C and sides <= Lmax, their radii rho = (k-1) 2^-e (or (k+1) 2^-e with
    radius='up'), and the pieces of their rectangles inside C. Integer coordinates in units 1/G, G = Dc 2^e."""

    def __init__(self, meas, Lmax, e=E_TOMAR, radius="tomar"):
        t0 = time.time()
        self.meas, self.e, self.Lmax, self.radius = meas, e, Fraction(Lmax), radius
        Px, Py, atoms = meas.Px, meas.Py, meas.atoms
        self.Dc = Dc = lcm_den([Px, Py] + [a[0] for a in atoms] + [a[1] for a in atoms])
        self.Dm = Dm = lcm_den([a[2] for a in atoms])
        self.G = Dc << e
        self.PX, self.PY = int(Px * Dc), int(Py * Dc)
        self.W, self.H = self.PX << e, self.PY << e
        ln, ld = meas.lam.numerator, meas.lam.denominator
        xs = sorted({int(a[0] * Dc) for a in atoms})
        ys = sorted({int(a[1] * Dc) for a in atoms})
        nx, ny = len(xs), len(ys)
        ix = {v: i for i, v in enumerate(xs)}
        iy = {v: i for i, v in enumerate(ys)}
        base = [[0] * ny for _ in range(nx)]
        for x, y, m in atoms:
            base[ix[int(x * Dc)]][iy[int(y * Dc)]] += int(m * Dm)
        Lm = int(self.Lmax * Dc)  # floor: integer sides <= Lm  <=>  sides <= Lmax
        ax = (self.PX + Lm) // self.PX + 2
        ay = (self.PY + Lm) // self.PY + 2
        xe = [xs[i % nx] + (i // nx) * self.PX for i in range(nx * ax)]
        ye = [ys[j % ny] + (j // ny) * self.PY for j in range(ny * ay)]
        self.xe, self.ye = xe, ye
        P = [[0] * (len(ye) + 1) for _ in range(len(xe) + 1)]  # P[i][j] = mass of grid points with index < (i, j)
        for i in range(len(xe)):
            Pi, Pn, bi = P[i], P[i + 1], base[i % nx]
            for j in range(len(ye)):
                Pn[j + 1] = Pn[j] + Pi[j + 1] - Pi[j] + bi[j % ny]
        self.P = P
        Kc, Km = ld * Dc * Dc, Dm * ln  # m/lam > L^2  <=>  M ld Dc^2 > L^2 Dm ln
        num_sh, den = 2 * e, 4 * Dm * ln
        up = radius == "up"
        rects = []
        nbox = 0
        for i0 in range(nx):
            x0 = xe[i0]
            P0 = P[i0]
            for i1 in range(i0, len(xe)):
                Lx = xe[i1] - x0
                if Lx > Lm:
                    break
                P1 = P[i1 + 1]
                xa, xb = xe[i1] << e, x0 << e
                for j0 in range(ny):
                    y0 = ye[j0]
                    c0 = P0[j0] - P1[j0]
                    for j1 in range(j0, len(ye)):
                        Ly = ye[j1] - y0
                        if Ly > Lm:
                            break
                        nbox += 1
                        M = P1[j1 + 1] - P0[j1 + 1] + c0
                        if M == 0:
                            continue
                        L = Lx if Lx > Ly else Ly
                        if M * Kc <= L * L * Km:
                            continue
                        k = isqrt((M * ld << num_sh) // den)
                        r = ((k + 1) if up else (k - 1)) * Dc
                        a, b = xa - r, xb + r
                        if a >= b:
                            continue
                        c, d = (ye[j1] << e) - r, (y0 << e) + r
                        if c >= d:
                            continue
                        rects.append((a, b, c, d, i0, i1, j0, j1, M, k))
        self.nbox, self.rects = nbox, rects
        # pieces inside C = [0,W] x [0,H]: (x0, x1, y0, y1, rect index, s, t) for the shift (s W, t H)
        W, H = self.W, self.H
        pieces = []
        for n, (a, b, c, d) in enumerate(r[:4] for r in rects):
            for s in range((-b) // W + 1, -((a - W) // W)):
                u0, u1 = max(a + s * W, 0), min(b + s * W, W)
                if u0 >= u1:
                    continue
                for t in range((-d) // H + 1, -((c - H) // H)):
                    v0, v1 = max(c + t * H, 0), min(d + t * H, H)
                    if v0 < v1:
                        pieces.append((u0, u1, v0, v1, n, s, t))
        self.pieces = pieces
        self.t_build = time.time() - t0


# ------------------------------------------------------------------------------------------------- union areas

def subtract_all(pieces, W, H, nb=48):
    """Method A: subtract the pieces from [0,W] x [0,H] in order of decreasing area; return the uncovered pieces and
    the indices of the pieces that removed area. A bucket grid only selects candidates; every test is exact."""
    order = sorted(range(len(pieces)), key=lambda i: -(pieces[i][1] - pieces[i][0]) * (pieces[i][3] - pieces[i][2]))
    bw, bh = -(-W // nb), -(-H // nb)
    live = {0: (0, W, 0, H)}
    buckets = {}

    def put(pid, u):
        for bx in range(u[0] // bw, min((u[1] - 1) // bw, nb - 1) + 1):
            for by in range(u[2] // bh, min((u[3] - 1) // bh, nb - 1) + 1):
                buckets.setdefault((bx, by), set()).add(pid)

    def drop(pid, u):
        for bx in range(u[0] // bw, min((u[1] - 1) // bw, nb - 1) + 1):
            for by in range(u[2] // bh, min((u[3] - 1) // bh, nb - 1) + 1):
                buckets[(bx, by)].discard(pid)

    put(0, live[0])
    nid = 1
    used = []
    for i in order:
        x0, x1, y0, y1 = pieces[i][:4]
        cand = set()
        for bx in range(x0 // bw, min((x1 - 1) // bw, nb - 1) + 1):
            for by in range(y0 // bh, min((y1 - 1) // bh, nb - 1) + 1):
                cand |= buckets.get((bx, by), set())
        hit = False
        for pid in cand:
            u0, u1, v0, v1 = u = live[pid]
            if u1 <= x0 or x1 <= u0 or v1 <= y0 or y1 <= v0:
                continue
            hit = True
            drop(pid, u)
            del live[pid]
            m0, m1 = max(u0, x0), min(u1, x1)
            for w in ((u0, x0, v0, v1) if u0 < x0 else None, (x1, u1, v0, v1) if x1 < u1 else None,
                      (m0, m1, v0, y0) if v0 < y0 else None, (m0, m1, y1, v1) if y1 < v1 else None):
                if w is not None:
                    live[nid] = w
                    put(nid, w)
                    nid += 1
        if hit:
            used.append(i)
    return list(live.values()), used


def union_area_sweep(rs):
    """Method B: area of a union of rectangles (x0, x1, y0, y1, ...) by slabs between consecutive x edges, the
    y-intervals of the rectangles spanning each slab sorted and merged."""
    xs = sorted({r[0] for r in rs} | {r[1] for r in rs})
    by_start = sorted(rs, key=lambda r: r[0])
    active, p, total = [], 0, 0
    for k in range(len(xs) - 1):
        xl, xr = xs[k], xs[k + 1]
        while p < len(by_start) and by_start[p][0] <= xl:
            active.append(by_start[p])
            p += 1
        active = [r for r in active if r[1] > xl]
        cov, end = 0, None
        for y0, y1 in sorted((r[2], r[3]) for r in active):
            if end is None or y0 > end:
                cov += y1 - y0
                end = y1
            elif y1 > end:
                cov += y1 - end
                end = y1
        total += (xr - xl) * cov
    return total


def overlap_free(unc, pieces, nb=48, W=None, H=None):
    """True when no uncovered piece meets any piece of the family in a set of positive area (bucketed, exact)."""
    bw, bh = -(-W // nb), -(-H // nb)
    buckets = {}
    for n, u in enumerate(unc):
        for bx in range(u[0] // bw, min((u[1] - 1) // bw, nb - 1) + 1):
            for by in range(u[2] // bh, min((u[3] - 1) // bh, nb - 1) + 1):
                buckets.setdefault((bx, by), []).append(n)
    for x0, x1, y0, y1 in (p[:4] for p in pieces):
        for bx in range(x0 // bw, min((x1 - 1) // bw, nb - 1) + 1):
            for by in range(y0 // bh, min((y1 - 1) // bh, nb - 1) + 1):
                for n in buckets.get((bx, by), ()):
                    u0, u1, v0, v1 = unc[n]
                    if u0 < x1 and x0 < u1 and v0 < y1 and y0 < v1:
                        return False
    return True


def pairwise_disjoint(unc):
    s = sorted(unc)
    for a in range(len(s)):
        u0, u1, v0, v1 = s[a]
        for b in range(a + 1, len(s)):
            w0, w1, z0, z1 = s[b]
            if w0 >= u1:
                break
            if v0 < z1 and z0 < v1:
                return False
    return True


# ------------------------------------------------------------------------------------- Lemma 6.2's hypotheses

def check_witness(fam, piece):
    """Exact check of one piece: S = atoms in its box, summed directly over the periodic copies; 0 < rho;
    4 lam rho^2 < mu(S); and piece - (s Px, t Py) inside R_S(rho) = [max x - rho, min x + rho] x [max y - rho, ...]."""
    meas, Dc, e = fam.meas, fam.Dc, fam.e
    a, b, c, d, i0, i1, j0, j1, M, k = fam.rects[piece[4]]
    bx0, bx1, by0, by1 = fam.xe[i0], fam.xe[i1], fam.ye[j0], fam.ye[j1]
    PX, PY = fam.PX, fam.PY
    mass, xs, ys = Fraction(0), [], []
    for x, y, m in meas.atoms:
        X, Y = int(x * Dc), int(y * Dc)
        for s in range(-((X - bx0) // PX), (bx1 - X) // PX + 1):
            for t in range(-((Y - by0) // PY), (by1 - Y) // PY + 1):
                mass += m
                xs.append(X + s * PX)
                ys.append(Y + t * PY)
    if mass * fam.Dm != M or not xs:
        return "mass differs from the prefix sum"
    rho = Fraction((k + 1) if fam.radius == "up" else (k - 1), 2 ** e)
    if rho <= 0:
        return "rho <= 0"
    if fam.radius != "up" and not 4 * meas.lam * rho * rho < mass:
        return "4 lam rho^2 >= mu(S)"
    r = int(rho * fam.G)
    sx, sy = piece[5] * fam.W, piece[6] * fam.H
    if not (piece[0] - sx >= (max(xs) << e) - r and piece[1] - sx <= (min(xs) << e) + r and
            piece[2] - sy >= (max(ys) << e) - r and piece[3] - sy <= (min(ys) << e) + r):
        return "piece not inside R_S(rho)"
    if max(xs) - min(xs) > int(fam.Lmax * Dc) or max(ys) - min(ys) > int(fam.Lmax * Dc):
        return "box longer than the cut-off"
    return None


# ------------------------------------------------------------------------------------------- one evaluation

class Result:
    pass


def evaluate(meas, Lmax, e=E_TOMAR, radius="tomar", log=print, witness_check=True):
    """Exact |U n C| for the family, by methods A and B with the two-sided check, and the exact ratio."""
    t0 = time.time()
    fam = Family(meas, Lmax, e, radius)
    t1 = time.time()
    unc, used = subtract_all(fam.pieces, fam.W, fam.H)
    t2 = time.time()
    K = [fam.pieces[i] for i in used]
    area_B = union_area_sweep(K)
    t3 = time.time()
    area_A = fam.W * fam.H - sum((u[1] - u[0]) * (u[3] - u[2]) for u in unc)
    disjoint = pairwise_disjoint(unc) and all(0 <= u[0] < u[1] <= fam.W and 0 <= u[2] < u[3] <= fam.H for u in unc)
    free = overlap_free(unc, fam.pieces, W=fam.W, H=fam.H)
    t4 = time.time()
    bad = []
    if witness_check:
        for p in K:
            why = check_witness(fam, p)
            if why:
                bad.append(why)
    t5 = time.time()
    R = Result()
    R.fam, R.unc, R.K = fam, unc, K
    R.area_A, R.area_B = area_A, area_B
    R.exact = area_A == area_B and disjoint and free
    R.area = Fraction(area_B, fam.G * fam.G)
    R.ratio = meas.lam * R.area / meas.mass
    R.bad = bad
    R.times = (t1 - t0, t2 - t1, t3 - t2, t4 - t3, t5 - t4)
    if log:
        log(f"  boxes {fam.nbox:,} (sides <= {float(Lmax):.6g}), rectangles of positive area {len(fam.rects):,}, "
            f"pieces in C {len(fam.pieces):,} ({R.times[0]:.0f} s)")
        log(f"  method A (subtraction): {len(unc):,} uncovered pieces, {len(K):,} pieces removed area "
            f"({R.times[1]:.0f} s); method B (slab sweep) on those {len(K):,} ({R.times[2]:.0f} s)")
        log(f"  |C| - uncovered = {'union K (B)' if area_A == area_B else 'DIFFERS FROM union K'}; uncovered pieces "
            f"pairwise disjoint: {disjoint}; no piece of the family meets them: {free} ({R.times[3]:.0f} s)")
        if witness_check:
            log(f"  Lemma 6.2 hypotheses, re-summed from the atoms, for all {len(K):,} rectangles used: "
                f"{'all hold' if not bad else 'FAIL: ' + bad[0]} ({R.times[4]:.0f} s)")
    return R


# ---------------------------------------------------------------------------------------- formatting helpers

def dec(x, n=20):
    """Decimal expansion of a non-negative rational, truncated to n places."""
    i = x.numerator * 10 ** n // x.denominator
    s = str(i).rjust(n + 1, "0")
    return s[:-n] + "." + s[-n:]


def printed_check(F, printed, log):
    ok = True
    for p in printed:
        b = Fraction(p)
        nd = len(p.split(".")[1])
        up = b + Fraction(1, 10 ** nd)
        imp = F >= b
        nimp = not (F >= up)
        log(f"  printed bound {p}: {'implied' if imp else 'NOT IMPLIED'} (margin {float(F - b):.3e}); "
            f"one last-digit unit higher, {dec(up, nd)}: {'not implied, as expected' if nimp else 'IMPLIED (bad)'}")
        ok &= imp and nimp
    return ok


# ------------------------------------------------------------------------------------------------- self-tests

class Iv:
    """Closed rational intervals for the closed forms (positive operands where it matters)."""

    def __init__(self, lo, hi=None):
        self.lo, self.hi = Fraction(lo), Fraction(lo if hi is None else hi)

    def __add__(s, o):
        o = o if isinstance(o, Iv) else Iv(o)
        return Iv(s.lo + o.lo, s.hi + o.hi)

    __radd__ = __add__

    def __sub__(s, o):
        o = o if isinstance(o, Iv) else Iv(o)
        return Iv(s.lo - o.hi, s.hi - o.lo)

    def __rsub__(s, o):
        return Iv(o) - s

    def __mul__(s, o):
        o = o if isinstance(o, Iv) else Iv(o)
        c = [s.lo * o.lo, s.lo * o.hi, s.hi * o.lo, s.hi * o.hi]
        return Iv(min(c), max(c))

    __rmul__ = __mul__

    def __truediv__(s, o):
        o = o if isinstance(o, Iv) else Iv(o)
        assert o.lo > 0
        return s * Iv(1 / o.hi, 1 / o.lo)

    def sqrt(s, bits=300):
        assert s.lo >= 0
        lo = Fraction(isqrt(s.lo.numerator * 4 ** bits // s.lo.denominator), 2 ** bits)
        hi = Fraction(isqrt(-(-s.hi.numerator * 4 ** bits // s.hi.denominator)) + 1, 2 ** bits)
        return Iv(lo, hi)


def lin_phi():
    """Phi from the closed form on the 47a page, and from the one in Lin's README, as rational intervals."""
    r22, r2 = Iv(22).sqrt(), Iv(2).sqrt()
    h = (5 + r22) / 6
    V = h + 1
    w = (17 + 4 * r22) / 9
    u = (2 + r22) / 3
    a = 2 * h - (2 * (2 + w)).sqrt() / 2 - u / 2
    b = V - r2 * u / 2 - w.sqrt() / 2
    page = (2 * h * V - 4 * a * b) / (1 + w)
    r11 = Iv(11).sqrt()
    readme = (Iv(77) / 2 + 8 * r22 - (8 + r22 - (70 + 8 * r22).sqrt()) *
              (11 + r22 - 2 * r2 - 2 * r11 - (17 + 4 * r22).sqrt())) / (26 + 4 * r22)
    return page, readme, 4 * a * b


def lin_measure(digits=60):
    """Lin's lattice with u = (2 + sqrt 22)/3 rounded down to 10^-digits and w, h, V formed from it exactly."""
    u = (2 + Fraction(isqrt(22 * 10 ** (2 * digits)), 10 ** digits)) / 3
    w = u * u - 1
    h = (1 + u) / 2
    V = h + 1
    return Measure(2 * h, V, [(0, 0, 1), (h, 0, w)])


def covered_at(meas, X, Y, e, reach, ideal=False):
    """The definition, pointwise: X, Y rationals. Sort the atoms by sup-distance d from (X, Y); with S(d) the atoms at
    distance <= d, the point is covered when rho_e(mu(S(d))) >= d (or 4 lam d^2 <= mu(S(d)) when ideal) for some d,
    the bounding box of S(d) having sides <= reach. Only d <= L*/2 can work: for 2d > L*, mu(Q(x, d)) < 4 d^2."""
    X, Y = Fraction(X), Fraction(Y)
    if getattr(meas, "_half", None) is None:
        meas._half = cutoff(meas) / 2
        meas._Dc = lcm_den([meas.Px, meas.Py] + [a[0] for a in meas.atoms] + [a[1] for a in meas.atoms])
        meas._Dm = lcm_den([a[2] for a in meas.atoms])
    Dm, ln, ld = meas._Dm, meas.lam.numerator, meas.lam.denominator
    S = lcm_den([Fraction(1, meas._Dc), X, Y])  # integer coordinates in units 1/S, exact
    Xi, Yi, PX, PY = int(X * S), int(Y * S), int(meas.Px * S), int(meas.Py * S)
    W = meas._half.numerator * S // meas._half.denominator  # d <= L*/2  <=>  d S <= floor(L* S / 2)
    Rch = None if reach is None else reach.numerator * S // reach.denominator
    pts = []
    for x, y, m in meas.atoms:
        x, y, M = int(x * S), int(y * S), int(m * Dm)
        for s in range(-((x - Xi + W) // PX), (Xi + W - x) // PX + 1):
            px = x + s * PX
            dx = abs(px - Xi)
            for t in range(-((y - Yi + W) // PY), (Yi + W - y) // PY + 1):
                py = y + t * PY
                dy = abs(py - Yi)
                pts.append((dx if dx > dy else dy, px, py, M))
    pts.sort()
    mass = 0
    x0 = x1 = y0 = y1 = None
    n = 0
    while n < len(pts):
        dd = pts[n][0]
        while n < len(pts) and pts[n][0] == dd:
            _, px, py, M = pts[n]
            mass += M
            x0 = px if x0 is None or px < x0 else x0
            x1 = px if x1 is None or px > x1 else x1
            y0 = py if y0 is None or py < y0 else y0
            y1 = py if y1 is None or py > y1 else y1
            n += 1
        if Rch is not None and (x1 - x0 > Rch or y1 - y0 > Rch):
            continue
        if ideal:
            if 4 * ln * dd * dd * Dm <= mass * ld * S * S:  # 4 lam d^2 <= mu(S(d))
                return True
        else:
            k = isqrt((mass * ld << (2 * e)) // (4 * Dm * ln))
            if (k - 1) * S >= dd << e:  # rho = (k - 1) 2^-e >= d
                return True
    return False


def brute_area(meas, Lmax, e):
    """|union| for a small measure from the definition: every box enumerated by direct loops with direct sums gives
    the candidate edges (mod the periods); each elementary cell of that grid is decided at its centre by covered_at."""
    Px, Py = meas.Px, meas.Py
    W = Lmax + 1
    copies = []
    for x, y, m in meas.atoms:
        for s in range(int(-W // Px) - 1, int((Px + W) // Px) + 2):
            for t in range(int(-W // Py) - 1, int((Py + W) // Py) + 2):
                copies.append((x + s * Px, y + t * Py, m))
    cx = sorted({c[0] for c in copies})
    cy = sorted({c[1] for c in copies})
    ex, ey = {Fraction(0), Px}, {Fraction(0), Py}
    for xl in [v for v in cx if 0 <= v < Px]:
        for xh in [v for v in cx if xl <= v <= xl + Lmax]:
            for yl in [v for v in cy if 0 <= v < Py]:
                for yh in [v for v in cy if yl <= v <= yl + Lmax]:
                    mass = sum(c[2] for c in copies if xl <= c[0] <= xh and yl <= c[1] <= yh)
                    if mass == 0:
                        continue
                    k = isqrt(int(mass / (4 * meas.lam) * 4 ** e))
                    rho = Fraction(k - 1, 2 ** e)
                    if xh - rho < xl + rho and yh - rho < yl + rho:
                        ex |= {(xh - rho) % Px, (xl + rho) % Px}
                        ey |= {(yh - rho) % Py, (yl + rho) % Py}
    ex, ey = sorted(ex), sorted(ey)
    area = Fraction(0)
    for i in range(len(ex) - 1):
        for j in range(len(ey) - 1):
            if covered_at(meas, (ex[i] + ex[i + 1]) / 2, (ey[j] + ey[j + 1]) / 2, e, Lmax):
                area += (ex[i + 1] - ex[i]) * (ey[j + 1] - ey[j])
    return area, (len(ex) - 1) * (len(ey) - 1)


def selftest(certs):
    ok = True
    log = print
    log("1. Lin's Phi (47a page): the closed form, and the evaluator on Lin's lattice")
    page, readme, slots = lin_phi()
    log(f"  47a closed form  Phi in [{dec(page.lo, 30)}, {dec(page.hi, 30)}]")
    log(f"  README form      Phi in [{dec(readme.lo, 30)}, {dec(readme.hi, 30)}]")
    good = (page.lo <= readme.hi and readme.lo <= page.hi and dec(page.lo, 10) == dec(page.hi, 10) == "1.6855099933"
            and dec(page.lo, 23) == "1.68550999335552518466528")
    log(f"  both agree, round down to 1.6855099933 and match Lin's 1.68550999335552518466528...: {good}")
    ok &= good
    lin = lin_measure(60)
    L = cutoff(lin)
    R = evaluate(lin, L, e=300, log=log)
    diff = max(abs(R.ratio - page.lo), abs(R.ratio - page.hi))
    good = R.exact and not R.bad and diff < Fraction(1, 10 ** 50)
    log(f"  u = (2 + sqrt 22)/3 to 60 digits, radii on 2^-300, every box up to L* = {float(L):.4f}: "
        f"{len(R.unc)} uncovered pieces, ratio = {dec(R.ratio, 30)}..., |ratio - Phi| < {float(diff):.1e}: {good}")
    ok &= good
    hole = sum(Fraction((u[1] - u[0]) * (u[3] - u[2]), R.fam.G ** 2) for u in R.unc)
    good = abs(hole - slots.lo) < Fraction(1, 10 ** 50)
    log(f"  uncovered area {dec(hole, 30)}... = the four slots 4ab = {dec(slots.lo, 30)}...: {good}")
    ok &= good
    big = [u for u in R.unc if (u[1] - u[0]) * (u[3] - u[2]) * 10 ** 20 > R.fam.G ** 2]
    inside = [covered_at(lin, Fraction(u[0] + u[1], 2 * R.fam.G), Fraction(u[2] + u[3], 2 * R.fam.G), 300, None,
                         ideal=True) for u in big]
    good = len(big) == 4 and not any(inside)
    log(f"  the uncovered pieces of area > 10^-20 are {len(big)} slots, and M mu < 1 at each centre "
        f"(from the definition, exact radii): {good}")
    ok &= good
    R40 = evaluate(lin, L, e=40, log=None)
    log(f"  same lattice with Tomar's radii on 2^-40: ratio = {dec(R40.ratio, 20)}..., below Phi by "
        f"{float(page.lo - R40.ratio):.2e}")
    ok &= R40.exact and R40.ratio < page.lo

    log("2. A single atom per cell (mass 1, period 10 x 10; only the atom's own square counts):")
    one = Measure(10, 10, [(3, 4, 1)])
    R = evaluate(one, cutoff(one), log=None)
    want = Fraction(4 * (2 ** 39 - 1) ** 2, 2 ** 80)
    good = R.exact and R.ratio == want
    log(f"  ratio = {R.ratio} = (1 - 2^-39)^2: {good}")
    ok &= good

    log("3. Tomar's measure without its 48 small atoms (the paper: 'the lattice alone gives only about 1.6447'):")
    path = os.path.join(certs, ROWS["tom2026"]["file"][0])
    if os.path.exists(path):
        d, _ = load_json(path, ROWS["tom2026"]["file"][1])
        m = Measure.from_json(d)
        bare = m.copy(atoms=m.atoms[:2])
        R = evaluate(bare, cutoff(bare), log=None)
        good = R.exact and abs(R.ratio - Fraction(16447, 10000)) < Fraction(5, 10 ** 5)
        log(f"  atoms {m.atoms[0][2]} at (0,0) and {float(m.atoms[1][2]):.10f} at (Px/2, 0): ratio = "
            f"{dec(R.ratio, 10)}: about 1.6447: {good}")
        ok &= good
    else:
        log("  (skipped: run fetch.sh first)")

    log("4. Brute force from the definition on small random measures (every grid cell decided by covered_at):")
    rng = random.Random(47)
    for case in range(8):
        n = rng.choice([1, 2, 2, 3])
        Px = Fraction(rng.randint(130, 260), 100)
        Py = Fraction(rng.randint(130, 260), 100)
        atoms, used = [], set()
        while len(atoms) < n:
            x, y = Fraction(rng.randrange(0, int(Px * 20)), 20), Fraction(rng.randrange(0, int(Py * 20)), 20)
            if (x, y) not in used:
                used.add((x, y))
                atoms.append((x, y, Fraction(rng.randint(20, 160), 100)))
        meas = Measure(Px, Py, atoms)
        if meas.mass >= Px * Py * Fraction(9, 10):
            atoms = [(x, y, m / 2) for x, y, m in atoms]
            meas = Measure(Px, Py, atoms)
        L = cutoff(meas)
        if case >= 6:  # a reach shorter than L*, as in the rows
            L = Fraction(int(L * 100) // 3, 100)
        e = 12
        R = evaluate(meas, L, e=e, log=None)
        ba, ncell = brute_area(meas, L, e)
        good = R.exact and not R.bad and R.area == ba
        log(f"  case {case}: {n} atom(s), period {Px} x {Py}, mass {float(meas.mass):.2f}, "
            f"{'L*' if case < 6 else 'reach'} = {float(L):.3f}, radii on 2^-{e}: fast {dec(R.ratio, 12)}, "
            f"brute ({ncell} cells) {dec(ba / meas.mass, 12)}: {'identical' if good else 'DIFFERENT'}")
        ok &= good
        pts = 0
        for _ in range(40):
            X, Y = Fraction(rng.randrange(0, 10 ** 6), 10 ** 6) * Px, Fraction(rng.randrange(0, 10 ** 6), 10 ** 6) * Py
            G = R.fam.G
            XG, YG = X * G, Y * G
            inside = any(p[0] <= XG <= p[1] and p[2] <= YG <= p[3] for p in R.fam.pieces)
            if inside != covered_at(meas, X, Y, e, L):
                ok = False
                log(f"    point ({X}, {Y}): union {inside}, definition {not inside}")
            pts += 1
        log(f"    and at {pts} random points the union of the pieces and the definition agree")
    log(f"RESULT: selftest {'PASS' if ok else 'FAIL'}")
    return ok


# ------------------------------------------------------------------------------------------------------- rows

def lin2026(controls):
    """The lower-bound row [Lin2026] of 47a (PR #192): c_2 >= Phi = 1.6855099933..., from the measure the row states.
    No file is read: Lin's lattice is built from the row's closed forms, with u = (2 + sqrt 22)/3 rounded down to
    10^-60 and w = u^2 - 1, h = (1 + u)/2, V = h + 1 formed from it exactly (the identities behind the row's w, h, V).
    That rational lattice is itself a periodic measure, so the exact ratio our evaluator returns for it is a lower
    bound for c_2 by Lemma 6.2 (as re-derived in the header), whatever the closed form says."""
    log = print
    t0 = time.time()
    ok = True
    log("47a row [Lin2026] (PR #192): c_2 >= Phi = (2hV - 4ab)/(1 + w), Lin's periodic measure")
    page, _, slots = lin_phi()
    log(f"  the row's closed form: Phi in [{dec(page.lo, 30)}, {dec(page.hi, 30)}]")
    good = dec(page.lo, 10) == dec(page.hi, 10) == "1.6855099933"
    log(f"  it rounds down to the row's 1.6855099933: {good}")
    ok &= good
    lin = lin_measure(60)
    L = cutoff(lin)
    log(f"  the measure: period {dec(lin.Px, 12)}... x {dec(lin.Py, 12)}..., masses 1 at (0, 0) and "
        f"w = {dec(lin.atoms[1][2], 12)}... at (h, 0); every box up to L* = {float(L):.4f}, radii on the 2^-300 grid")
    R = evaluate(lin, L, e=300, log=log)
    good = R.exact and not R.bad
    log(f"  ratio = lam |U n C| / mu(C) = {dec(R.ratio, 30)}... (exact rational, {len(str(R.ratio.denominator))}-digit "
        f"denominator): {'exact, both methods agree' if good else 'NOT EXACT'}")
    ok &= good
    diff = max(abs(R.ratio - page.lo), abs(R.ratio - page.hi))
    good = diff < Fraction(1, 10 ** 50)
    log(f"  |ratio - Phi| < {float(diff):.1e}: {good}")
    ok &= good
    ok &= printed_check(R.ratio, ["1.6855099933"], log)
    hole = sum(Fraction((u[1] - u[0]) * (u[3] - u[2]), R.fam.G ** 2) for u in R.unc)
    good = abs(hole - slots.lo) < Fraction(1, 10 ** 50)
    log(f"  uncovered area {dec(hole, 30)}... = the row's four slots 4ab = {dec(slots.lo, 30)}...: {good}")
    ok &= good
    big = [u for u in R.unc if (u[1] - u[0]) * (u[3] - u[2]) * 10 ** 20 > R.fam.G ** 2]
    inside = [covered_at(lin, Fraction(u[0] + u[1], 2 * R.fam.G), Fraction(u[2] + u[3], 2 * R.fam.G), 300, None,
                         ideal=True) for u in big]
    good = len(big) == 4 and not any(inside)
    log(f"  the uncovered pieces of area > 10^-20 are {len(big)} slots, and M mu < 1 at each centre (every closed "
        f"square about it, from the atoms): {good}")
    ok &= good
    if controls:
        log("Red controls (the same evaluator, other measures; each must fall below the row's 1.6855099933):")
        uu = (2 + Fraction(isqrt(22 * 10 ** 120), 10 ** 60)) / 3
        for what, (uh, ww) in (("u + 1/100 (h, V, w moved together)", (uu + Fraction(1, 100), None)),
                               ("w + 1/10, h and V kept", (uu, (uu * uu - 1) + Fraction(1, 10))),
                               ("w = 1 (equal masses), h and V kept", (uu, Fraction(1)))):
            h = (1 + uh) / 2
            w = uh * uh - 1 if ww is None else ww
            m = Measure(2 * h, h + 1, [(0, 0, 1), (h, 0, w)])
            Rc = evaluate(m, cutoff(m), e=300, log=None)
            good = Rc.exact and Rc.ratio < Fraction(16855099933, 10 ** 10)
            log(f"  {what}: ratio = {dec(Rc.ratio, 12)}..., {'below, as expected' if good else 'NOT BELOW'}")
            ok &= good
    log(f"RESULT: {'PASS' if ok else 'FAIL'} (c_2 >= {dec(R.ratio, 12)}... >= 1.6855099933)  ({time.time() - t0:.0f} s)")
    return ok


def row_header(key, meas, sha, fname, log):
    cfg = ROWS[key]
    log(f"Row {key} ({cfg['pr']}): {fname}, sha256 {sha}")
    log(f"  {len(meas.atoms)} atoms per cell, Px = {meas.Px} = {float(meas.Px):.12f}, Py = {meas.Py} = "
        f"{float(meas.Py):.12f}, lam = {meas.lam}, reach = {meas.reach}")
    log(f"  mu(C) = {meas.mass} = {float(meas.mass):.15f}")


def compare(key, F, log):
    cfg = ROWS[key]
    ok = True
    log(f"  exact ratio lam |U n C| / mu(C) = {F.numerator}/{F.denominator}")
    log(f"                                  = {dec(F, 25)}...")
    if cfg.get("claimed") is not None:
        same = F == cfg["claimed"]
        log(f"  stated fraction {cfg['claimed']}: {'IDENTICAL' if same else 'DIFFERENT (stated ' + dec(cfg['claimed'], 25) + ')'}")
        ok &= same
    if cfg.get("decimal"):
        same = dec(F, len(cfg["decimal"].split(".")[1])) == cfg["decimal"]
        log(f"  stated decimals {cfg['decimal']}...: {'agree' if same else 'DISAGREE'}")
        ok &= same
    ok &= printed_check(F, cfg["printed"], log)
    return ok


def perturbed_ratio(meas):
    R = evaluate(meas, meas.reach, log=None)
    return R.ratio, R.exact and not R.bad


def check_row(key, certs, controls, full):
    t0 = time.time()
    log = print
    cfg = ROWS[key]
    path = os.path.join(certs, cfg["file"][0])
    if not os.path.exists(path):
        log(f"Row {key} ({cfg['pr']}): {cfg['file'][0]} (sha256 {cfg['file'][1]}) is not in {certs}.")
        if key == "pr341":
            log("  Its only source, Zenodo record 10.5281/zenodo.24022791, is restricted (files: HTTP 403), so the")
            log("  63-atom measure could not be replayed. Arithmetic on the stated fraction alone:")
            printed_check(cfg["claimed"], cfg["printed"], log)
            log(f"  stated fraction = {dec(cfg['claimed'], 25)}...")
        log(f"RESULT: {key} NOT RUN (no data)")
        return None
    d, sha = load_json(path, cfg["file"][1])
    meas = Measure.from_json(d["measure"] if "measure" in d else d)
    row_header(key, meas, sha, cfg["file"][0], log)
    ok = True
    if cfg.get("witness_only"):
        # our own full family has 1.3 million rectangles here (over our memory cap), so the certificate's witnesses
        # are checked one by one instead, which is what the bound rests on
        Fw, ok1 = check_certificate(key, certs, meas, None, log)
        ok &= ok1 and compare(key, Fw, log)
        ok &= finite_check(key, certs, meas, log)
        ok &= sample_definition(key, certs, meas, 400, log)
        if controls:
            ok &= witness_controls(key, certs, meas, Fw, log)
        log(f"RESULT: {key} {'PASS' if ok else 'FAIL'}")
        log(f"run time {time.time() - t0:.0f} s")
        return ok
    e = cfg.get("e", E_TOMAR)
    log(f"Tomar's family (boxes with sides <= reach = {meas.reach}, radii (k-1) 2^-{e}):")
    R = evaluate(meas, meas.reach, e=e, log=log)
    ok &= R.exact and not R.bad
    F = R.ratio
    ok &= compare(key, F, log)
    holes = sum(Fraction((u[1] - u[0]) * (u[3] - u[2]), R.fam.G ** 2) for u in R.unc)
    log(f"  uncovered: {len(R.unc):,} pieces, area {float(holes):.10f} = {float(holes / (meas.Px * meas.Py)):.4%} of C")
    ok &= definition_checks(meas, R, 200, log)
    if "witnesses" in cfg:
        ok &= check_certificate(key, certs, meas, F, log)
    if "finite" in cfg:
        ok &= finite_check(key, certs, meas, log)
    if full:
        ok &= full_family(meas, F, log)
    if controls:
        ok &= red_controls(key, meas, F, log)
    log(f"RESULT: {key} {'PASS' if ok else 'FAIL'}")
    log(f"run time {time.time() - t0:.0f} s")
    return ok


def definition_checks(meas, R, n_random, log):
    """The pointwise definition (covered_at, every atom within L*/2, Tomar's radii and reach) against the union: no
    centre of an uncovered piece may be covered, and at random points of C the union and the definition must agree."""
    t0 = time.time()
    G, e = R.fam.G, R.fam.e
    bad = sum(covered_at(meas, Fraction(u[0] + u[1], 2 * G), Fraction(u[2] + u[3], 2 * G), e, meas.reach)
              for u in R.unc)
    rng = random.Random(4747)
    disagree = 0
    for _ in range(n_random):
        X = Fraction(rng.randrange(0, 10 ** 12), 10 ** 12) * meas.Px
        Y = Fraction(rng.randrange(0, 10 ** 12), 10 ** 12) * meas.Py
        XG, YG = X * G, Y * G
        inside = any(p[0] <= XG <= p[1] and p[2] <= YG <= p[3] for p in R.fam.pieces)
        disagree += inside != covered_at(meas, X, Y, e, meas.reach)
    log(f"  the definition, point by point: centres of the {len(R.unc):,} uncovered pieces covered: {bad}; "
        f"{n_random} random points where union and definition disagree: {disagree} ({time.time() - t0:.0f} s)")
    return bad == 0 and disagree == 0


def check_certificate(key, certs, meas, F, log):
    """PR #321's 559 witness rectangles, each checked exactly; the union of their lattice translates by method B."""
    path, sha = ROWS[key]["witnesses"]
    d, got = load_json(os.path.join(certs, path), sha)
    log(f"The certificate {path} (sha256 {got}): {len(d['witnesses'])} witnesses")
    mm = Measure.from_json(d["measure"])
    same = (mm.Px, mm.Py, mm.lam, mm.reach, sorted(mm.atoms)) == (meas.Px, meas.Py, meas.lam, meas.reach,
                                                                  sorted(meas.atoms))
    log(f"  its embedded measure equals the measure file: {same}")
    ok = same
    Px, Py = meas.Px, meas.Py
    rects, bad = [], []
    for w in d["witnesses"]:
        bx0, bx1, by0, by1 = (q(v) for v in w["box"])
        rho = q(w["radius"])
        x0, x1, y0, y1 = (q(v) for v in w["rectangle"])
        mass, xs, ys = Fraction(0), [], []
        for x, y, m in meas.atoms:
            for s in range(int((bx0 - x) // Px) - 1, int((bx1 - x) // Px) + 2):
                for t in range(int((by0 - y) // Py) - 1, int((by1 - y) // Py) + 2):
                    px, py = x + s * Px, y + t * Py
                    if bx0 <= px <= bx1 and by0 <= py <= by1:
                        mass += m
                        xs.append(px)
                        ys.append(py)
        if not xs or not rho > 0 or not 4 * meas.lam * rho * rho < mass:
            bad.append(("radius", w))
            continue
        if not (max(xs) - rho <= x0 < x1 <= min(xs) + rho and max(ys) - rho <= y0 < y1 <= min(ys) + rho):
            bad.append(("rectangle", w))
            continue
        if not (0 <= x0 and x1 <= Px and 0 <= y0 and y1 <= Py):
            bad.append(("outside C", w))
            continue
        rects.append((x0, x1, y0, y1))
    log(f"  every witness: S = atoms in its box (summed from the atoms), 0 < rho, 4 lam rho^2 < mu(S), rectangle "
        f"inside R_S(rho) and inside C: {'all hold' if not bad else 'FAIL (' + str(len(bad)) + ', first: ' + bad[0][0] + ')'}")
    ok &= not bad
    D = lcm_den([v for r in rects for v in r])
    area = Fraction(union_area_sweep([tuple(int(v * D) for v in r) for r in rects]), D * D)
    Fw = meas.lam * area / meas.mass
    if F is None:
        log(f"  lam |union of the witnesses| / mu(C) (method B, exact) = {dec(Fw, 25)}...")
        return Fw, ok
    log(f"  lam |union of the witnesses| / mu(C) (method B) = {dec(Fw, 25)}...: "
        f"{'equal to our ratio' if Fw == F else 'DIFFERENT from our ratio'}")
    ok &= Fw == F
    return ok


def finite_check(key, certs, meas, log):
    """The row's finite value: q ((N - 2h)/N)^2 with the PR's N and h equals the stated fraction, and h min(Px, Py)
    exceeds every witness radius (so every atom a witness of an inner cell needs is in the N x N block)."""
    N, h, stated = ROWS[key]["finite"]
    q_ = ROWS[key]["claimed"]
    val = q_ * Fraction(N - 2 * h, N) ** 2
    d, _ = load_json(os.path.join(certs, ROWS[key]["witnesses"][0]), ROWS[key]["witnesses"][1])
    rmax = max(q(w["radius"]) for w in d["witnesses"])
    good = val == stated and h * min(meas.Px, meas.Py) > rmax
    log(f"  the finite instance: q ((N - 2h)/N)^2 with N = {N}, h = {h} equals the stated fraction "
        f"{dec(stated, 25)}...: {val == stated}; h min(Px, Py) = {float(h * min(meas.Px, meas.Py)):.4f} > "
        f"the largest witness radius {float(rmax):.4f}: {h * min(meas.Px, meas.Py) > rmax}")
    ok = good and printed_check(stated, ROWS[key]["printed"], log)
    return ok


def sample_definition(key, certs, meas, n, log):
    """at random points of random witness rectangles, M mu >= lam from the definition (every closed square about the
    point, at each atom distance; exact rationals): a check of the witnesses against the definition, not a proof."""
    import random
    d, _ = load_json(os.path.join(certs, ROWS[key]["witnesses"][0]), ROWS[key]["witnesses"][1])
    rnd = random.Random(5)
    good = True
    for _ in range(n):
        w = rnd.choice(d["witnesses"])
        x0, x1, y0, y1 = (q(v) for v in w["rectangle"])
        X = x0 + (x1 - x0) * Fraction(rnd.randint(0, 1000), 1000)
        Y = y0 + (y1 - y0) * Fraction(rnd.randint(0, 1000), 1000)
        good &= covered_at(meas, X, Y, 300, None, ideal=True)
    log(f"  at {n} random points of the witness rectangles, M mu >= lam from the definition: {good}")
    return good


def witness_controls(key, certs, meas, Fw, log):
    """Each must be refused or not implied: a witness radius raised past sqrt(mu(S)/(4 lam)), a rectangle widened past
    R_S(rho), the heaviest atom lightened by 10% (some witnesses fail), and the bound one unit higher."""
    d, _ = load_json(os.path.join(certs, ROWS[key]["witnesses"][0]), ROWS[key]["witnesses"][1])
    ok = True
    log("Red controls:")

    def check_one(mm, w):
        bx0, bx1, by0, by1 = (q(v) for v in w["box"])
        rho = q(w["radius"])
        x0, x1, y0, y1 = (q(v) for v in w["rectangle"])
        mass, xs, ys = Fraction(0), [], []
        for x, y, m in mm.atoms:
            for s_ in range(int((bx0 - x) // mm.Px) - 1, int((bx1 - x) // mm.Px) + 2):
                for t_ in range(int((by0 - y) // mm.Py) - 1, int((by1 - y) // mm.Py) + 2):
                    px, py = x + s_ * mm.Px, y + t_ * mm.Py
                    if bx0 <= px <= bx1 and by0 <= py <= by1:
                        mass += m
                        xs.append(px)
                        ys.append(py)
        return bool(xs) and rho > 0 and 4 * mm.lam * rho * rho < mass and \
            max(xs) - rho <= x0 < x1 <= min(xs) + rho and max(ys) - rho <= y0 < y1 <= min(ys) + rho
    w0 = dict(d["witnesses"][0])
    bx0, bx1, by0, by1 = (q(v) for v in w0["box"])
    mass = sum(m for x, y, m in meas.atoms if bx0 <= x <= bx1 and by0 <= y <= by1)
    big = dict(w0, radius=str(Fraction(isqrt(int(mass * 10 ** 40 / (4 * meas.lam))) + 10 ** 6, 10 ** 20)))
    good = not check_one(meas, big)
    log(f"  C1 the first witness with its radius raised above sqrt(mu(S)/(4 lam)): {'refused, as expected' if good else 'ACCEPTED (bad)'}")
    ok &= good
    x0, x1, y0, y1 = (q(v) for v in w0["rectangle"])
    xs = [x + s_ * meas.Px for x, y, m in meas.atoms for s_ in range(-3, 4) for t_ in range(-3, 4)
          if bx0 <= x + s_ * meas.Px <= bx1 and by0 <= y + t_ * meas.Py <= by1]
    edge = max(xs) - q(w0["radius"])                   # the left edge of R_S(rho)
    wide = dict(w0, rectangle=[str(edge - Fraction(1, 10 ** 12)), str(x1), str(y0), str(y1)])
    good = not check_one(meas, wide)
    log(f"  C2 the first witness rectangle extended to 10^-12 past the left edge of R_S(rho): "
        f"{'refused, as expected' if good else 'ACCEPTED (bad)'}")
    ok &= good
    h = max(range(len(meas.atoms)), key=lambda n: meas.atoms[n][2])
    light = meas.copy(atoms=meas.atoms[:h] + [(meas.atoms[h][0], meas.atoms[h][1], meas.atoms[h][2] * Fraction(9, 10))] +
                      meas.atoms[h + 1:])
    fails = sum(1 for w in d["witnesses"] if not check_one(light, w))
    good = fails > 0
    log(f"  C3 the heaviest atom lightened by 10%: {fails} witnesses fail, as expected: {good}")
    ok &= good
    p = ROWS[key]["printed"][-1]
    nd = len(p.split(".")[1])
    up = Fraction(p) + Fraction(1, 10 ** nd)
    N, hh, stated = ROWS[key]["finite"]
    good = not (stated >= up) and not (Fw >= up)
    log(f"  C4 the bound {dec(up, nd)}, one unit higher in the last printed digit: "
        f"{'not implied, as expected' if good else 'IMPLIED (bad)'}")
    ok &= good
    return ok


def full_family(meas, F, log):
    """Every box up to L*: with Tomar's radii (is the reach binding?), and with radii (k -/+ 1) 2^-100, which enclose
    |{M mu >= lam} n C| (the union is monotone in every radius, and boxes longer than L* have no interior)."""
    ok = True
    L = cutoff(meas)
    log(f"Every box up to L* = {float(L):.6f} (beyond it mu(box) <= side^2, so no rectangle has interior):")
    R = evaluate(meas, L, log=log, witness_check=False)
    ok &= R.exact
    log(f"  ratio with Tomar's radii: {dec(R.ratio, 25)}...: "
        f"{'the same as with reach (the reach is not binding)' if R.ratio == F else 'LARGER than with reach'}")
    lo = evaluate(meas, L, e=100, log=None, witness_check=False)
    hi = evaluate(meas, L, e=100, radius="up", log=None, witness_check=False)
    ok &= lo.exact and hi.exact and lo.ratio <= hi.ratio
    log(f"  lam |{{M mu >= lam}} n C| / mu(C) in [{dec(lo.ratio, 25)}, {dec(hi.ratio, 25)}]")
    log(f"  (radii on 2^-100 below and above sqrt(mu(S)/4 lam); the 2^-40 radii give up "
        f"{float(lo.ratio - F):.3e} of it)")
    return ok


def red_controls(key, meas, F, log):
    """Each must fail: a moved atom, a changed mass, a wrong period (the fraction is no longer the stated one), a
    bound one unit higher in the last printed digit (not implied), and one oversized radius (refused by Lemma 6.2)."""
    ok = True
    target = ROWS[key].get("claimed") or F
    log("Red controls:")
    atoms = list(meas.atoms)
    i = len(atoms) - 1
    x, y, m = atoms[i]
    step = Fraction(1, 10 ** 6) if x + Fraction(1, 10 ** 6) < meas.Px else -Fraction(1, 10 ** 6)
    moved = meas.copy(atoms=atoms[:i] + [(x + step, y, m)])
    r, fine = perturbed_ratio(moved)
    good = fine and r != target
    log(f"  C1 last atom moved by {step} in x: ratio {dec(r, 16)}..., "
        f"{'differs from the stated fraction, as expected' if good else 'UNCHANGED (bad)'}")
    ok &= good
    h = max(range(len(atoms)), key=lambda n: atoms[n][2])
    changed = meas.copy(atoms=atoms[:h] + [(atoms[h][0], atoms[h][1], atoms[h][2] + Fraction(1, 10 ** 12))] +
                        atoms[h + 1:])
    r, fine = perturbed_ratio(changed)
    good = fine and r != target
    log(f"  C2 heaviest mass raised by 10^-12: ratio {dec(r, 16)}..., "
        f"{'differs from the stated fraction, as expected' if good else 'UNCHANGED (bad)'}")
    ok &= good
    wrong = meas.copy(Px=meas.Px + Fraction(1, 10 ** 9))
    r, fine = perturbed_ratio(wrong)
    good = fine and r != target
    log(f"  C3 period Px + 10^-9: ratio {dec(r, 16)}..., "
        f"{'differs from the stated fraction, as expected' if good else 'UNCHANGED (bad)'}")
    ok &= good
    p = ROWS[key]["printed"][-1]
    nd = len(p.split(".")[1])
    up = Fraction(p) + Fraction(1, 10 ** nd)
    good = not (F >= up)
    log(f"  C4 bound {dec(up, nd)} (one unit higher in the last printed digit): "
        f"{'not implied, as expected' if good else 'IMPLIED (bad)'}")
    ok &= good
    fam = Family(meas, meas.reach)
    big = max(range(len(fam.rects)), key=lambda n: (fam.rects[n][1] - fam.rects[n][0]))
    a, b, c, d, i0, i1, j0, j1, M, k = fam.rects[big]
    fam.rects[big] = (a, b, c, d, i0, i1, j0, j1, M, k + 2)  # rho = (k+1) 2^-40 > sqrt(mu(S)/4)
    piece = next(pc for pc in fam.pieces if pc[4] == big)
    why = check_witness(fam, piece)
    good = why is not None
    log(f"  C5 one radius raised by two grid units, to (k+1) 2^-40: {'refused (' + why + ')' if good else 'ACCEPTED (bad)'}")
    ok &= good
    return ok


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if not args:
        print(__doc__)
        return 2
    certs = args[1] if len(args) > 1 else os.path.join(HERE, "certs")
    if args[0] == "selftest":
        return 0 if selftest(certs) else 1
    if args[0] == "lin2026":
        return 0 if lin2026("--controls" in sys.argv) else 1
    if args[0] == "measure":
        d, sha = load_json(args[1], None)
        meas = Measure.from_json(d)
        print(f"{args[1]} (sha256 {sha}): {len(meas.atoms)} atoms, mu(C) = {float(meas.mass):.12f}")
        R = evaluate(meas, meas.reach if meas.reach is not None else cutoff(meas))
        print(f"  ratio = {R.ratio.numerator}/{R.ratio.denominator} = {dec(R.ratio, 25)}...; exact: {R.exact}")
        return 0 if R.exact and not R.bad else 1
    if args[0] not in ROWS:
        print(f"unknown row {args[0]}; rows: {', '.join(ROWS)}")
        return 2
    ok = check_row(args[0], certs, "--controls" in sys.argv, "--full" in sys.argv)
    return 0 if ok or ok is None else 1


if __name__ == "__main__":
    sys.exit(main())
