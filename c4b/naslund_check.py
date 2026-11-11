"""Replay of the 4b row [Nas2026b], C_4b >= 0.7580758 (exactly alpha = 473797394875551/625000000000000 =
0.7580758318008816; Naslund, hexagon:2610.00151v1, recorded in PR #219), from the paper's LaTeX source and the four
witness files of the archive embedded in it (bound.json, q215-depth3-product.json, q437-unequal32s.json,
binary-cert-0.7580758318008816.json; their SHA-256 as the paper's appendix states). No program of the archive or of
the repository was read or run: fetch.sh decodes the archive bytes from the TeX source and keeps only the four JSON
files.

The criterion (the paper's Theorem 3.3, read, not computed): ordered interval alphabets C_i in pairwise coprime square
bases b_i, widths in (0, 1), with Z_{f_i}(C_i) = sum_x w_x^{f_i} >= b_i^alpha and sum_i f_i > alpha, give
r(N) >= N^alpha / Q, so C_4b >= alpha. We check every finite hypothesis, exactly or with rigorous enclosures:

 1. Six Paley chains (Certificate 6.2): p prime, p = 3 mod 4, entries distinct mod p, and for every ordered pair
    s_j - s_i is a nonzero square mod p exactly when j > i; the contribution test (1 - f_p) log t > (2 alpha - 1) log p
    (Corollary 4.5) with the table's f_p.
 2. The (5,43) low code (Certificate 6.4): 4,913 distinct low points with depths (3, 3), widths n/10^14; the relation
    a ~> b of Definition 4.7 (at every prime the low words agree, or their first differing digit differs by a nonzero
    square, b minus a), computed for all ordered pairs; every related pair goes forward in the list (so it is acyclic),
    the heaviest directed path weighs at most 99999974386028/10^14 < 1 (Lemma 4.9 then supplies the intervals); and
    log sum w^f - (2 alpha - 1) log(215^3) > 0 (Corollary 4.10).
 3. The (19,23) low code (Certificate 6.6): 3,645 distinct low points with depths (3, 2), starts and widths over
    10^16, intervals in [0, 1], and sup I_a <= inf I_b for every related pair; the moment test with F = 19^3 23^2.
 4. The binary policy (Certificate 6.8, Definitions 5.2 and 5.4): 25 states, 94 branches; windows closed, of positive
    length, inside [0, 1]; at most one branch per residue r; the containment (5.2) and cyclic-order (5.3) conditions
    for the images J_{r+4 pi} of the child windows under reflection and z -> theta + u z; the root (state 1) has both
    windows [0, 1] and v_1 = 1; the 25 growth rows sum u^f v_j >= a v_s; c = min d_s/v_s; and the depth condition
    log c + (m - 1) log a - f log 2 - alpha m log 4 > 0 at m = 10^20 (Corollary 5.6).
 5. The budget: the nine contributions sum to more than alpha, exactly; the prime supports are disjoint.
Logarithms and powers: Python's decimal ln and exp (correctly rounded) at 70 digits, each result moved one unit in the
last place outwards, so every quantity is enclosed. The relation in 2 and 3 is computed twice: by prefix bit masks
for every pair, and by a direct per-pair function for 400 random targets (all sources), which must agree.

Usage: python3 -I naslund_check.py nas2026b|nas2026a CERTS_ROOT [--controls]
"""
import hashlib
import json
import os
import random
import sys
from decimal import Decimal as D, Context, ROUND_HALF_EVEN, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction as Fr

P = 70
CN = Context(prec=P, rounding=ROUND_HALF_EVEN)
CF = Context(prec=P, rounding=ROUND_FLOOR)
CC = Context(prec=P, rounding=ROUND_CEILING)

ROWS = {
    # Naslund, hexagon:2610.00151v1 (7 Oct 2026), repository commit 5a83f41; the row of PR #219
    "nas2026b": dict(
        pr=219, value="0.7580758", alpha=Fr(473797394875551, 625000000000000),
        tex_sha="4472c9a90c6282fb39e6768929c5a0a6ff28c4ec1f8289c89b4f1e16306881c6", commit="5a83f41",
        files={"bound.json": "080735a362d5fb736ae932d4c775e74c6c05ce3a25119e17cc0f1b928022bf84",
               "q215-depth3-product.json": "dae38a756de83c99d65973ff85f23a7068ce6e85a7a4f12400a923f24cf46574",
               "q437-unequal32s.json": "3ef1e78e59f6e90327bd739c8c114604a352bb4ae7f0faecc6d547b3cd0104bc",
               "binary-cert-0.7580758318008816.json": "1b388dd3cfda324a57d985fc96746a8d36b8e7d2e4ac25ecd264c854a8cb8af7"},
        odd215="q215-depth3-product.json", odd437="q437-unequal32s.json", binary="binary-cert-0.7580758318008816.json",
        F437=19 ** 3 * 23 ** 2, depths437=(3, 2),
        contrib={3: "0.1819189685063644313459576", 7: "0.0857695922261077385453455", 11: "0.1072043074619092931875010",
                 31: "0.0891366647927413831505830", 59: "0.0421426390460550181119771", 103: "0.0023654891206769647316637",
                 215: "0.0264983230432418544045862", 437: "0.0681150585658841330427119", 2: "0.1549247890379023302829202"},
        stated=dict(pairs215=6031254, maxpath215=99999974386028, m215="6.919e-18", pairs437=2702580,
                    maxright437=9999999218567850, m437="5.927e-18", rows="3.27e-19", growth="4.91e-19", depth="46.91",
                    c=">0.3551", surplus=Fr(15468032462, 10 ** 25), chains="2.4e-26")),
    # Naslund, September 2026 version (repository commit e5d6937); the row of PR #187
    "nas2026a": dict(
        pr=187, value="0.75806746", alpha=Fr(37903373, 50000000),
        tex_sha="d1ba46a136cd67b026e7b0fe9568c360a396947f81ed86efb85e9d63bb91af4b", commit="e5d6937",
        files={"bound.json": "836b231c32e771dfaf494a74de40b5213c8f8606865ac8309724e6a1fc5ec491",
               "q215-depth3-product.json": "dae38a756de83c99d65973ff85f23a7068ce6e85a7a4f12400a923f24cf46574",
               "q437-depth3-final.json": "7abd0c007c42fa1a5484ea8a4d02c56b04898209f74a0403c5503c03395e58d3",
               "certificate-alpha-0758067460.json": "4835df59b0cc8444abe274a2c0d923fa89d023687c52ec58a7c31a49a8b03c87"},
        odd215="q215-depth3-product.json", odd437="q437-depth3-final.json", binary="certificate-alpha-0758067460.json",
        F437=437 ** 3, depths437=(3, 3),
        contrib={3: "0.181945506487", 7: "0.085799249226", 11: "0.107233269134", 31: "0.089166212567",
                 59: "0.042173711274", 103: "0.002397851691", 215: "0.026537312662", 437: "0.067871875356",
                 2: "0.154942497274"},
        stated=dict(pairs215=6031254, maxpath215=99999974386028, m215="6.97840e-10", pairs437=77614486,
                    maxright437=9999999524038902, m437="6.95598e-10", rows="6.99293e-9", growth="3.05203e-9",
                    depth="28.32697028", c=Fr(285516532458, 803818168703), surplus=Fr(25671, 10 ** 12),
                    chains="7.60e-14")),
    # Jones, nk-lean @ 55acbf0 (Lean 4, Palomar PALOMAR-2026-09-21-000004); the row of PR #187
    "jon2026b": dict(
        pr=187, value="0.75806770413", alpha=Fr(75806770413, 10 ** 11), format="jones", commit="55acbf0",
        files={"certificate.json": "0e1a242044e62c8e4a240d7b3e00c379e963f9cc020a305c21985ed278fe29cc", "odd215.txt": "1ae65e9349a66e184c4376dd4196253d0abb32d3a0d0fc2dc265f7a273d4abd6",
               "odd437.txt": "fe026c70694a85892a749ea4417199325eb39ba79307c6351aa8de433585ffc0",
               "binary_policy.json": "0c7a47db2a19b619b42d731a9ab78013c2235e8f3ec90d2bfd7a1fe9114712c5"},
        F437=437 ** 3, depths437=(3, 3),
        stated=dict(surplus=Fr(38652, 10 ** 15), depth=10 ** 15, c=Fr(1, 4))),
}
ALPHA = None          # set per row in main()
CONTRIB = None


# ---------- rigorous decimal enclosures ----------
def dfrac_lo(q):
    return CF.divide(D(q.numerator), D(q.denominator))


def dfrac_hi(q):
    return CC.divide(D(q.numerator), D(q.denominator))


def ln_lo(x):
    """lower bound of ln(x) for a Decimal x > 0"""
    return CF.next_minus(CN.ln(x))


def ln_hi(x):
    return CC.next_plus(CN.ln(x))


def ln_q(q):
    """enclosure of ln(q) for a rational q > 0"""
    return ln_lo(dfrac_lo(q)), ln_hi(dfrac_hi(q))


def pow_lo(q, f):
    """lower bound of q^f for rational q > 0 and rational f >= 0"""
    lq = ln_q(q)
    fl, fh = dfrac_lo(f), dfrac_hi(f)
    # f * ln q: if ln q <= 0 use the larger f with the smaller ln, else the smaller f
    if lq[0] <= 0:
        t = CF.multiply(fh, lq[0])
    else:
        t = CF.multiply(fl, lq[0])
    return CF.next_minus(CN.exp(t))


def is_prime(n):
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


def squares_mod(p):
    return {(z * z) % p for z in range(1, p)} - {0}


# ---------- 1. chains ----------
def check_chain(p, s):
    Q = squares_mod(p)
    ok = is_prime(p) and p % 4 == 3 and len(set(x % p for x in s)) == len(s) and all(0 <= x < p for x in s)
    pairs = 0
    for i in range(len(s)):
        for j in range(len(s)):
            if i == j:
                continue
            pairs += 1
            sq = ((s[j] - s[i]) % p) in Q
            ok &= (sq == (j > i))
    return ok, pairs


def chain_margin(p, t, f):
    """lower bound of (1 - f) log t - (2 alpha - 1) log p"""
    lt = ln_q(Fr(t))
    lp = ln_q(Fr(p))
    one_minus_f = Fr(1) - f
    a = CF.multiply(dfrac_lo(one_minus_f), lt[0])
    b = CC.multiply(dfrac_hi(2 * ALPHA - 1), lp[1])
    return CF.subtract(a, b)


# ---------- 2, 3. low codes ----------
def rel_direct(a, b, primes, Q):
    """a ~> b (Definition 4.7), a != b, directly from the definition"""
    if a == b:
        return False
    for k, p in enumerate(primes):
        wa, wb = a[k], b[k]
        if wa == wb:
            continue
        j = next(i for i in range(len(wa)) if wa[i] != wb[i])
        if ((wb[j] - wa[j]) % p) not in Q[p]:
            return False
    return True


def related_masks(points, primes):
    """for every target b, the bit mask of the sources a with a ~> b, by prefix masks per prime"""
    n = len(points)
    Q = {p: squares_mod(p) for p in primes}
    masks = []                                   # per prime: dict prefix-tuple -> bit mask
    for k, p in enumerate(primes):
        pm = {}
        for idx, pt in enumerate(points):
            w = tuple(pt[k])
            for L in range(1, len(w) + 1):
                pm[w[:L]] = pm.get(w[:L], 0) | (1 << idx)
        masks.append(pm)
    cache = [dict() for _ in primes]
    out = []
    for b in points:
        allowed = (1 << n) - 1
        for k, p in enumerate(primes):
            wb = tuple(b[k])
            m = cache[k].get(wb)
            if m is None:
                m = masks[k].get(wb, 0)          # equal word
                for j in range(len(wb)):         # first difference at index j: equal prefix, square difference
                    pre = wb[:j]
                    for d in range(p):
                        if ((wb[j] - d) % p) in Q[p]:
                            m |= masks[k].get(pre + (d,), 0)
                cache[k][wb] = m
            allowed &= m
        out.append(allowed)
    for i in range(n):
        out[i] &= ~(1 << i)                     # a != b
    return out, Q


def bits(m):
    while m:
        low = m & -m
        yield low.bit_length() - 1
        m ^= low


def check_lowcode_215(d, controls_tweak=None):
    primes = d["primes"]
    pts = [tuple(tuple(w) for w in pt) for pt in d["points"]]
    den = d["width_denominator"]
    nums = list(d["width_numerators"])
    if controls_tweak:
        controls_tweak(pts, nums)
    n = len(pts)
    shape = all(len(pt) == 2 and len(pt[0]) == 3 and len(pt[1]) == 3 and all(0 <= x < primes[k] for k, w in enumerate(pt) for x in w) for pt in pts)
    distinct = len(set(pts)) == n
    pos = all(0 < x < den for x in nums)
    rel, Q = related_masks(pts, primes)
    npairs = sum(bin(m).count("1") for m in rel)
    forward = all((m >> i) == 0 for i, m in enumerate(rel))          # every source of b comes before b
    dist = [0] * n
    for i in range(n):
        best = 0
        for a in bits(rel[i]):
            if a < i and dist[a] > best:
                best = dist[a]
        dist[i] = best + nums[i]
    maxpath = max(dist)
    return dict(n=n, shape=shape, distinct=distinct, pos=pos, npairs=npairs, forward=forward, maxpath=maxpath,
                den=den, nums=nums, pts=pts, rel=rel, Q=Q, primes=primes)


def check_lowcode_437(d, controls_tweak=None, depths=(3, 2)):
    primes = d["primes"]
    pts = [tuple(tuple(w) for w in pt) for pt in d["points"]]
    den = d["width_denominator"]
    nums = list(d["width_numerators"])
    starts = list(d["start_numerators"])
    if controls_tweak:
        controls_tweak(pts, nums, starts)
    n = len(pts)
    ddep = tuple(d["depths"]) if "depths" in d else (d["depth"],) * len(d["primes"])
    depth_ok = ddep == tuple(depths)
    shape = depth_ok and all(len(pt) == 2 and len(pt[k]) == depths[k] and all(0 <= x < primes[k] for x in pt[k])
                             for pt in pts for k in range(2))
    distinct = len(set(pts)) == n
    inside = all(0 < w < den and 0 <= s and s + w <= den for s, w in zip(starts, nums))   # widths below one
    maxright = max(s + w for s, w in zip(starts, nums))
    rel, Q = related_masks(pts, primes)
    npairs = sum(bin(m).count("1") for m in rel)
    forward = all((m >> i) == 0 for i, m in enumerate(rel))
    bad = 0
    for i in range(n):
        for a in bits(rel[i]):
            if starts[a] + nums[a] > starts[i]:
                bad += 1
    return dict(n=n, shape=shape, distinct=distinct, inside=inside, maxright=maxright, npairs=npairs,
                forward=forward, bad=bad, den=den, nums=nums, pts=pts, rel=rel, Q=Q, primes=primes)


def moment_margin(nums, den, f, F):
    """lower bound of log(sum (n/den)^f) - (2 alpha - 1) log F"""
    s = D(0)
    for x in nums:
        s = CF.add(s, pow_lo(Fr(x, den), f))
    a = ln_lo(s)
    b = CC.multiply(dfrac_hi(2 * ALPHA - 1), ln_q(Fr(F))[1])
    return CF.subtract(a, b), s


def cross_check(res, rng, k=400):
    """direct per-pair relation for k random targets against every source: must equal the masks"""
    pts, rel, Q, primes = res["pts"], res["rel"], res["Q"], res["primes"]
    n = len(pts)
    ok = True
    for i in rng.sample(range(n), min(k, n)):
        m = 0
        for a in range(n):
            if rel_direct(pts[a], pts[i], primes, Q):
                m |= 1 << a
        ok &= (m == rel[i])
    return ok


# ---------- 4. binary policy ----------
def binary_check(d, f_override=None, scale_tweak=None, c_fixed=None):
    f = Fr(d["f"]) if f_override is None else f_override
    a = Fr(d["a"])
    m = d["depth"]
    root = d["root"]
    v = [Fr(x) for x in d["vector"]]
    states = d["states"]
    S = len(states)
    msgs = []
    ok = True

    def win(rec):
        lo1, hi1, lo2, hi2, q = rec
        E = None if (lo1, hi1) == (-1, -1) else (Fr(lo1, q), Fr(hi1, q))
        O = None if (lo2, hi2) == (-1, -1) else (Fr(lo2, q), Fr(hi2, q))
        return E, O

    W = [tuple(s["wins"]) if "wins" in s else win(s["windows"]) for s in states]
    for s, (E, O) in enumerate(W):
        present = [w for w in (E, O) if w is not None]
        if not present or any(not (0 <= w[0] < w[1] <= 1) for w in present):
            ok = False
            msgs.append(f"state {s}: windows not closed subintervals of [0,1] of positive length")
    n_cont = 0
    n_cyc = 0
    nbranches = 0
    for s, st in enumerate(states):
        E, O = W[s]
        rs = [br["r"] for br in st["branches"]]
        if len(set(rs)) != len(rs) or any(r not in (0, 1, 2, 3) for r in rs):
            ok = False
            msgs.append(f"state {s}: residues {rs} not distinct in 0..3")
        J = {}
        for bi, br in enumerate(st["branches"]):
            nbranches += 1
            r, j, tau, eta = br["r"], br["child"], br["swap"], br["reflect"]
            u = Fr(br["scale"][0], br["scale"][1])
            if scale_tweak and (s, bi) == scale_tweak[0]:
                u = u * scale_tweak[1]
            th = Fr(br["shift"][0], br["shift"][1])
            if not (0 <= j < S) or tau not in (0, 1) or eta not in (0, 1) or u <= 0:
                ok = False
                msgs.append(f"state {s} branch {bi}: bad fields")
                continue
            for pi in (0, 1):
                cw = W[j][pi ^ tau]                         # child's window of parity pi xor tau
                if cw is None:
                    continue
                lo, hi = cw
                if eta == 1:
                    lo, hi = 1 - hi, 1 - lo
                img = (th + u * lo, th + u * hi)
                k = r + 4 * pi
                J[k] = img
                parent = E if r % 2 == 0 else O
                n_cont += 1
                if parent is None or not (parent[0] <= img[0] and img[1] <= parent[1]):
                    ok = False
                    msgs.append(f"state {s} branch r={r} parity {pi}: containment fails")
        for k in range(8):
            k2 = (k + 1) % 8
            if k in J and k2 in J:
                n_cyc += 1
                if not (J[k][1] <= J[k2][0]):
                    ok = False
                    msgs.append(f"state {s}: cyclic order fails between classes {k} and {k2}")
    # root
    Er, Or = W[root]
    root_ok = Er == (0, 1) and Or == (0, 1) and v[root] == 1
    # growth rows: sum u^f v_j - a v_s, lower bounds
    row_lo = []
    for s, st in enumerate(states):
        acc = D(0)
        for bi, br in enumerate(st["branches"]):
            u = Fr(br["scale"][0], br["scale"][1])
            if scale_tweak and (s, bi) == scale_tweak[0]:
                u = u * scale_tweak[1]
            acc = CF.add(acc, CF.multiply(pow_lo(u, f), dfrac_lo(v[br["child"]])))
        row_lo.append(CF.subtract(acc, dfrac_hi(a * v[s])))
    rows_ok = all(x > 0 for x in row_lo) and all(x > 0 for x in v)
    # c = min_s d_s / v_s, d_s the length of the wider present window
    c_min = min(max((w[1] - w[0]) for w in W[s] if w is not None) / v[s] for s in range(S))
    c = c_min if c_fixed is None else c_fixed
    l2 = ln_q(Fr(2))
    la = ln_q(a)
    lc = ln_q(c)
    growth = CF.subtract(la[0], CC.multiply(dfrac_hi(ALPHA), CC.multiply(D(2), l2[1])))
    # log c + (m - 1) log a - f log 2 - alpha m log 4, lower bound
    depth = CF.add(lc[0], CF.multiply(D(m - 1), la[0]))
    depth = CF.subtract(depth, CC.multiply(dfrac_hi(f), l2[1]))
    depth = CF.subtract(depth, CC.multiply(CC.multiply(dfrac_hi(ALPHA), D(m)), CC.multiply(D(2), l2[1])))
    return dict(ok=ok, msgs=msgs, n_cont=n_cont, n_cyc=n_cyc, nbranches=nbranches, S=S, root_ok=root_ok,
                rows_ok=rows_ok, row_min=min(row_lo), c=c, c_min=c_min, growth=growth, depth=depth, f=f, m=m, a=a)


def load_jones_odd(path):
    """Jones's text format: a header 'p1 p2 count D', then per point three digits at p1, three at p2 (least
    significant first), the interval start and width numerators over D"""
    rows = open(path).read().split("\n")
    p1, p2, cnt, den = (int(x) for x in rows[0].split())
    pts, starts, nums = [], [], []
    for line in rows[1:]:
        if not line.strip():
            continue
        v = [int(x) for x in line.split()]
        assert len(v) == 8, "a row does not have 8 integers"
        pts.append((tuple(v[0:3]), tuple(v[3:6])))
        starts.append(v[6])
        nums.append(v[7])
    assert len(pts) == cnt, "the count in the header differs"
    return dict(primes=[p1, p2], depth=3, points=[[list(a), list(b)] for a, b in pts], width_denominator=den,
                width_numerators=nums, start_numerators=starts)


def jones_binary_to_naslund(pol, f, a, depth):
    """Jones's policy file in the field names of Naslund's witness: windows as [lo, hi] rationals (null = absent)"""
    states = []
    vec = []
    for st in pol["states"]:
        E, O = st["parity_windows"]
        ws = []
        for w in (E, O):
            ws.append(None if w is None else (Fr(w[0]), Fr(w[1])))
        brs = []
        for t in st["transitions"]:
            r, child, swap, refl, scale, shift = t
            u, th = Fr(scale), Fr(shift)
            brs.append({"r": r, "child": child, "swap": swap, "reflect": refl,
                        "scale": [u.numerator, u.denominator], "shift": [th.numerator, th.denominator]})
        states.append({"wins": ws, "branches": brs})
        vec.append(str(Fr(st["positive_weight"])))
    return {"f": str(f), "a": str(a), "depth": depth, "root": pol["root_state"], "vector": vec, "states": states,
            "alpha": str(ALPHA)}


def run_jones(key, R, root, controls):
    global ALPHA, CONTRIB
    ALPHA = R["alpha"]
    st = R["stated"]
    base = os.path.join(root, key)
    ok = True
    print(f"[Jon2026b] (PR #{R['pr']}): C_4b >= alpha = 75806770413/10^11 = 0.75806770413 (JD-Jones-ASES/nk-lean @ "
          f"{R['commit']}, certificates/square/)")
    for fn, h in R["files"].items():
        got = hashlib.sha256(open(os.path.join(base, fn), "rb").read()).hexdigest()
        note = "as the repository's README states" if fn == "odd437.txt" else "the file at the commit"
        print(f"  {fn} sha256 {got[:16]}...: {note}: {got == h}")
        ok &= got == h
    cert = json.load(open(os.path.join(base, "certificate.json")))
    ok &= Fr(cert["alpha"]) == ALPHA
    CONTRIB = {c["prime"]: Fr(c["f"]) for c in cert["chains"]}
    CONTRIB[215] = Fr(next(c["f"] for c in cert["odd"] if c["q"] == 215))
    CONTRIB[437] = Fr(next(c["f"] for c in cert["odd"] if c["q"] == 437))
    CONTRIB[2] = Fr(cert["binary"]["f"])
    print("\n1. the six Paley chains (the same chains as Naslund's; Jones's contributions)")
    tot_pairs = 0
    for c in cert["chains"]:
        p, svals = c["prime"], c["ordered_residues"]
        c_ok, pairs = check_chain(p, svals)
        tot_pairs += pairs
        mg = chain_margin(p, len(svals), CONTRIB[p])
        print(f"  p = {p:3d}, t = {len(svals):2d}: chain property on {pairs} ordered pairs: {c_ok}; f_p = {CONTRIB[p]}: "
              f"(1 - f_p) log t - (2 alpha - 1) log p >= {CN.plus(mg):.3E} > 0: {mg > 0}")
        ok &= c_ok and mg > 0
    ok &= tot_pairs == 244
    rng = random.Random(4)
    res = {}
    for q, F, label in ((215, 215 ** 3, "(5,43)"), (437, 437 ** 3, "(19,23)")):
        print(f"\n{2 if q == 215 else 3}. the {label} coupled-prime alphabet")
        d = load_jones_odd(os.path.join(base, f"odd{q}.txt"))
        r = check_lowcode_437(d, depths=(3, 3))
        res[q] = (r, d)
        print(f"  {r['n']} points, three digits per prime, digits in range: {r['shape']}; distinct: {r['distinct']}; "
              f"intervals in [0, 1] with positive widths: {r['inside']} (denominator {r['den']}; largest right endpoint "
              f"{r['maxright']})")
        print(f"  related pairs a ~> b: {r['npairs']:,}; sup I_a > inf I_b for {r['bad']} of them")
        cc = cross_check(r, rng)
        print(f"  the relation by a direct per-pair function for 400 random targets, all sources: equal: {cc}")
        mg, _ = moment_margin(r["nums"], r["den"], CONTRIB[q], F)
        print(f"  log sum w^f - (2 alpha - 1) log {q}^3 >= {CN.plus(mg):.4E} with f = {CONTRIB[q]}: {mg > 0}")
        ok &= r["shape"] and r["distinct"] and r["inside"] and r["bad"] == 0 and cc and mg > 0
    print("\n4. the binary policy (Naslund's geometry as Jones transcribed it; Jones's weights)")
    pol = json.load(open(os.path.join(base, "binary_policy.json")))
    f_b, a_b, m = Fr(pol["moment_power"]), Fr(pol["growth_lower"]), cert["binary"]["depth"]
    same = f_b == CONTRIB[2] and a_b == Fr(cert["binary"]["growth_lower"]) and m == st["depth"]
    db = jones_binary_to_naslund(pol, f_b, a_b, m)
    b = binary_check(db, c_fixed=st["c"])
    print(f"  {b['S']} states, {b['nbranches']} branches, depth m = {m}; f and a as in certificate.json: {same}")
    print(f"  windows, residues, containment and cyclic order: {b['ok']} ({b['n_cont']} containments, {b['n_cyc']} cyclic "
          f"comparisons)" + ("" if b["ok"] else f"; {b['msgs'][:3]}"))
    print(f"  root state {db['root']}: both windows [0, 1] and v = 1: {b['root_ok']}")
    print(f"  growth rows sum u^f v_j - a v_s: smallest >= {CN.plus(b['row_min']):.4E} > 0: {b['rows_ok']}")
    print(f"  the initialization factor c = 1/4 is at most min d_s/v_s = {float(b['c_min']):.6f}: {b['c_min'] >= st['c']}")
    print(f"  log a - alpha log 4 >= {CN.plus(b['growth']):.4E}: {b['growth'] > 0}")
    print(f"  log(1/4) + (m - 1) log a - f log 2 - alpha m log 4 >= {CN.plus(b['depth']):.6f} at m = 10^15: "
          f"{b['depth'] > 0}")
    ok &= same and b["ok"] and b["root_ok"] and b["rows_ok"] and b["c_min"] >= st["c"] and b["growth"] > 0
    ok &= b["depth"] > 0 and b["nbranches"] == 94 and b["S"] == 25
    print("\n5. the budget")
    total = sum(CONTRIB.values())
    surplus = total - ALPHA
    print(f"  sum f_i = alpha + {surplus} (the row: 38652/10^15): {surplus == st['surplus']}")
    ok &= surplus > 0 and surplus == st["surplus"]
    print("  so every hypothesis of the criterion holds: alpha = 0.75806770413 is attainable, C_4b >= 0.75806770413")
    cok = True
    if controls:
        print("\ncontrols (each must come out as stated):")
        r437, d437 = res[437]

        def shift(pts, nums, starts):
            for i in range(len(pts)):
                srcs = list(bits(r437["rel"][i]))
                if srcs:
                    starts[i] = starts[srcs[0]] + nums[srcs[0]] - 1
                    return
        qx = check_lowcode_437(d437, controls_tweak=shift, depths=(3, 3))
        c = qx["bad"] > 0
        print(f"[ctl] (19,23): one target start moved into its source's interval: {qx['bad']} violated pairs, refused: {c}")
        cok &= c
        mgy, _ = moment_margin(r437["nums"], r437["den"], CONTRIB[437] + Fr(1, 10 ** 12), 437 ** 3)
        c = mgy < 0
        print(f"[ctl] (19,23): contribution raised by 1e-12 fails the moment test (margin {CN.plus(mgy):.3E}): {c}")
        cok &= c
        by = binary_check(db, f_override=CONTRIB[2] + Fr(1, 10 ** 12), c_fixed=st["c"])
        c = not by["rows_ok"]
        print(f"[ctl] binary: f raised by 1e-12: smallest growth row {CN.plus(by['row_min']):.3E}, refused: {c}")
        cok &= c
        c = total - (ALPHA + 2 * surplus) < 0
        print(f"[ctl] the budget at alpha + 2 x the surplus fails: {c}")
        cok &= c
        print(f"controls all as expected: {cok}")
    resok = ok and cok
    print(f"\nRESULT: {'PASS' if resok else 'FAIL'} (every finite hypothesis of the criterion holds at alpha = "
          f"0.75806770413, so C_4b >= 0.75806770413)")
    return 0 if resok else 1


def main():
    global ALPHA, CONTRIB
    key = sys.argv[1]
    root = sys.argv[2]
    controls = "--controls" in sys.argv[3:]
    R = ROWS[key]
    if R.get("format") == "jones":
        return run_jones(key, R, root, controls)
    ALPHA, CONTRIB, st = R["alpha"], R["contrib"], R["stated"]
    base = os.path.join(root, key)
    ok = True
    tex = os.path.join(base, "paper.tex")
    th = hashlib.sha256(open(tex, "rb").read()).hexdigest()
    print(f"[{key[:3].capitalize()}{key[3:]}] (PR #{R['pr']}): C_4b >= alpha = {ALPHA.numerator}/{ALPHA.denominator} = "
          f"{float(ALPHA):.16g} (the row prints {R['value']})")
    print(f"paper TeX (enaslund/sarkozy-lower-bound-0.758 @ {R['commit']}) sha256 {th}: as pinned {th == R['tex_sha']}")
    ok &= th == R["tex_sha"]
    for fn, h in R["files"].items():
        got = hashlib.sha256(open(os.path.join(base, fn), "rb").read()).hexdigest()
        print(f"  {fn} sha256 {got[:16]}...: as the paper's appendix states {got == h}")
        ok &= got == h
    bound = json.load(open(os.path.join(base, "bound.json")))
    man = Fr(bound["alpha"]) == ALPHA and [Fr(c["power"]) for c in bound["odd_components"]] == \
        [Fr(CONTRIB[215]), Fr(CONTRIB[437])]
    print(f"  the manifest's alpha and odd contributions equal the paper's: {man}")
    ok &= man
    # 1. chains
    print("\n1. the six Paley chains")
    chains = [(p, s) for p, s in bound["prime_chains"]]
    tot_pairs = 0
    worst = None
    for p, s in chains:
        c_ok, pairs = check_chain(p, s)
        tot_pairs += pairs
        f = Fr(CONTRIB[p])
        mg = chain_margin(p, len(s), f)
        worst = mg if worst is None or mg < worst else worst
        print(f"  p = {p:3d}, t = {len(s):2d}: chain property on {pairs} ordered pairs: {c_ok}; f_p = {CONTRIB[p]}: "
              f"(1 - f_p) log t - (2 alpha - 1) log p >= {CN.plus(mg):.3E} > 0: {mg > 0}")
        ok &= c_ok and mg > 0
    print(f"  ordered pairs checked: {tot_pairs} (the paper: 244); smallest margin {CN.plus(worst):.3E} (the paper: "
          f"{st['chains']})")
    ok &= tot_pairs == 244
    # 2. (5,43)
    rng = random.Random(4)
    print("\n2. the (5,43) low code")
    d215 = json.load(open(os.path.join(base, R["odd215"])))
    r = check_lowcode_215(d215)
    print(f"  {r['n']} points, depths (3, 3), digits in range: {r['shape']}; distinct: {r['distinct']}; widths in (0,1): "
          f"{r['pos']}")
    print(f"  related pairs a ~> b: {r['npairs']:,} (the paper: {st['pairs215']:,}); every one goes forward in the list "
          f"(so ~> is acyclic): {r['forward']}")
    print(f"  heaviest directed path: {r['maxpath']}/10^14 (stated {st['maxpath215']}): "
          f"{r['maxpath'] == st['maxpath215']}; below 1: {r['maxpath'] < r['den']}")
    cc = cross_check(r, rng)
    print(f"  the relation by a direct per-pair function for 400 random targets, all sources: equal: {cc}")
    mg215, s215 = moment_margin(r["nums"], r["den"], Fr(CONTRIB[215]), 215 ** 3)
    print(f"  log sum w^f - (2 alpha - 1) log 215^3 >= {CN.plus(mg215):.4E} (the paper: > {st['m215']}): {mg215 > 0}")
    ok &= r["shape"] and r["distinct"] and r["pos"] and r["forward"] and r["maxpath"] < r["den"] and cc and mg215 > 0
    ok &= r["npairs"] == st["pairs215"] and r["maxpath"] == st["maxpath215"]
    # 3. (19,23)
    print("\n3. the (19,23) low code")
    d437 = json.load(open(os.path.join(base, R["odd437"])))
    q = check_lowcode_437(d437, depths=R["depths437"])
    print(f"  {q['n']} points, depths {R['depths437']}, digits in range: {q['shape']}; distinct: {q['distinct']}; "
          f"intervals in [0, 1] with positive widths: {q['inside']}; largest right endpoint {q['maxright']}/10^16 "
          f"(stated {st['maxright437']}): {q['maxright'] == st['maxright437']}")
    print(f"  related pairs a ~> b: {q['npairs']:,} (the paper: {st['pairs437']:,}); sup I_a > inf I_b for {q['bad']} of "
          f"them; all go forward in the list: {q['forward']}")
    cc = cross_check(q, rng)
    print(f"  the relation by a direct per-pair function for 400 random targets, all sources: equal: {cc}")
    mg437, s437 = moment_margin(q["nums"], q["den"], Fr(CONTRIB[437]), R["F437"])
    print(f"  log sum w^f - (2 alpha - 1) log F >= {CN.plus(mg437):.4E} with F = {R['F437']} (the paper: > {st['m437']}): "
          f"{mg437 > 0}")
    ok &= q["shape"] and q["distinct"] and q["inside"] and q["bad"] == 0 and cc and mg437 > 0
    ok &= q["npairs"] == st["pairs437"] and q["maxright"] == st["maxright437"]
    # 4. binary
    print("\n4. the binary policy")
    db = json.load(open(os.path.join(base, R["binary"])))
    b = binary_check(db)
    same_f = Fr(db["f"]) == Fr(CONTRIB[2]) and Fr(db["alpha"]) == ALPHA
    print(f"  {b['S']} states, {b['nbranches']} branches, depth m = {b['m']}; f_B and alpha as in the table and the "
          f"manifest: {same_f}")
    print(f"  windows, residues, containment (5.2) and cyclic order (5.3): {b['ok']} ({b['n_cont']} containments, each "
          f"two endpoint comparisons, and {b['n_cyc']} cyclic comparisons: {2 * b['n_cont'] + b['n_cyc']}, the paper's "
          f"304)" + ("" if b["ok"] else f"; {b['msgs'][:3]}"))
    print(f"  root state {db['root']}: both windows [0, 1] and v = 1: {b['root_ok']}")
    print(f"  growth rows sum u^f v_j - a v_s: smallest >= {CN.plus(b['row_min']):.4E} > 0 (the paper: > {st['rows']}): "
          f"{b['rows_ok']}")
    cst = st["c"]
    c_ok = (b["c"] == cst) if isinstance(cst, Fr) else (b["c"] > Fr(3551, 10000))
    print(f"  c = min d_s/v_s = {b['c'].numerator if b['c'].denominator < 10 ** 15 else ''}"
          f"{'/' + str(b['c'].denominator) if b['c'].denominator < 10 ** 15 else ''} = {float(b['c']):.6f} (the paper: "
          f"{cst}): {c_ok}")
    print(f"  log a - alpha log 4 >= {CN.plus(b['growth']):.4E} (the paper: > {st['growth']}): {b['growth'] > 0}")
    print(f"  log c + (m - 1) log a - f log 2 - alpha m log 4 >= {CN.plus(b['depth']):.6f} (the paper: > {st['depth']}): "
          f"{b['depth'] > 0}")
    ok &= b["ok"] and b["root_ok"] and b["rows_ok"] and b["growth"] > 0 and b["depth"] > 0 and same_f and c_ok
    ok &= 2 * b["n_cont"] + b["n_cyc"] == 304 and b["nbranches"] == 94 and b["S"] == 25
    # 5. budget
    print("\n5. the budget")
    total = sum(Fr(x) for x in CONTRIB.values())
    surplus = total - ALPHA
    supports = [{2}, {3}, {7}, {11}, {31}, {59}, {103}, {5, 43}, {19, 23}]
    disjoint = all(not (x & y) for i, x in enumerate(supports) for y in supports[i + 1:])
    print(f"  sum f_i = {total} = alpha + {surplus} (the paper: alpha + {st['surplus']}): {surplus == st['surplus']}; "
          f"prime supports pairwise disjoint: {disjoint}")
    ok &= surplus > 0 and disjoint and surplus == st["surplus"]
    print(f"  so every hypothesis of the criterion holds: alpha = {float(ALPHA):.16g} is attainable, C_4b >= {R['value']}")
    cok = True
    if controls:
        print("\ncontrols (each must come out as stated):")
        p, s = chains[3]
        s2 = list(s)
        s2[1], s2[2] = s2[2], s2[1]
        c = not check_chain(p, s2)[0]
        print(f"[ctl] the p = 31 chain with two entries swapped is refused: {c}")
        cok &= c
        mgx = chain_margin(103, 11, Fr(1, 2))
        c = mgx < 0
        print(f"[ctl] the p = 103 chain at contribution 1/2 fails the test (margin {CN.plus(mgx):.3E}): {c}")
        cok &= c
        slack = r["den"] - r["maxpath"]

        def widen(pts, nums):
            for i in range(len(nums)):
                nums[i] += slack + 1
        ry = check_lowcode_215(d215, controls_tweak=widen)
        c = ry["maxpath"] >= r["den"]
        print(f"[ctl] (5,43): every width raised by the slack {slack} + 1: heaviest path {ry['maxpath']} >= 10^14, "
              f"refused: {c}")
        cok &= c

        def reverse(pts, nums):
            pts.reverse()
            nums.reverse()
        rz = check_lowcode_215(d215, controls_tweak=reverse)
        c = not rz["forward"]
        print(f"[ctl] (5,43): the list reversed is not a topological order (refused): {c}")
        cok &= c
        bump = Fr(1, 10 ** 15) if key == "nas2026b" else Fr(1, 10 ** 8)
        mgy, _ = moment_margin(r["nums"], r["den"], Fr(CONTRIB[215]) + bump, 215 ** 3)
        c = mgy < 0
        print(f"[ctl] (5,43): contribution raised by {float(bump):.0e} fails the moment test (margin "
              f"{CN.plus(mgy):.3E}): {c}")
        cok &= c

        def shift437(pts, nums, starts):
            for i in range(len(pts)):
                srcs = list(bits(q["rel"][i]))
                if srcs:
                    a0 = srcs[0]
                    starts[i] = starts[a0] + nums[a0] - 1
                    return
        qx = check_lowcode_437(d437, controls_tweak=shift437, depths=R["depths437"])
        c = qx["bad"] > 0
        print(f"[ctl] (19,23): one target start moved 1/10^16 into its source's interval: {qx['bad']} violated pairs, "
              f"refused: {c}")
        cok &= c
        bx = binary_check(db, scale_tweak=((1, 2), Fr(10001, 10000)))
        c = not bx["ok"]
        print(f"[ctl] binary: one scale of state 1 raised by 0.01%: {bx['msgs'][0] if bx['msgs'] else 'no failure'}; "
              f"refused: {c}")
        cok &= c
        fb = Fr(1, 10 ** 17) if key == "nas2026b" else Fr(1, 10 ** 8)
        by = binary_check(db, f_override=Fr(CONTRIB[2]) + fb)
        c = not by["rows_ok"]
        print(f"[ctl] binary: f_B raised by {float(fb):.0e}: smallest growth row {CN.plus(by['row_min']):.3E}, "
              f"refused: {c}")
        cok &= c
        c = sum(Fr(x) for x in CONTRIB.values()) - (ALPHA + 2 * st["surplus"]) < 0
        print(f"[ctl] the budget at alpha + 2 x the surplus fails: {c}")
        cok &= c
        print(f"controls all as expected: {cok}")
    res = ok and cok
    print(f"\nRESULT: {'PASS' if res else 'FAIL'} (every finite hypothesis of the criterion holds at alpha = "
          f"{float(ALPHA):.16g}, so C_4b >= {R['value']})")
    return 0 if res else 1


if __name__ == "__main__":
    sys.exit(main())
