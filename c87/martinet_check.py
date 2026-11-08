"""87a: Martinet's row (Invent. Math. 44 (1978), Proposition 4.1 and Example 4.2, p. 69), as PR #332 restates it:
C87 <= 4 sqrt(70035) = 1058.5650664...  The field k' = Q(sqrt 2, sqrt(3.5.7.23.29)) is totally real, has an infinite
2-class field tower, and root discriminant 4 sqrt(70035).  Replayed with our own code (Quantavane replays);
standard library only, exact integer arithmetic; run with python3 -I.

Martinet's argument, which we re-derive step by step:
  k  = Q(sqrt p, sqrt 2) with p = 29 (totally real, degree n = 4);  m = q1 q2 q3 q4 = 3.5.7.23;
  K  = k(sqrt m) = Q(sqrt 2, sqrt 29, sqrt 2415), totally real of degree 8;
  k' = Q(sqrt 2, sqrt(pm)) = Q(sqrt 2, sqrt 70035), and K = k'(sqrt p) is unramified over k';
  section 4 (from section 2: genus theory and the Golod-Shafarevich inequality, p = 2, r1 = n, r2 = rho = 0,
  delta = 1): the 2-class field tower of K is infinite once the number r of primes of k ramified in K satisfies
  r >= n + 3 + 2 sqrt(2n + 1), i.e. r >= 13 for n = 4.  Since K/k' is unramified, the tower of k' is infinite too.

Checks (all exact):
  1. the hypotheses (i)-(iii) of Proposition 4.1 for p = 29, (q1..q4) = (3, 5, 7, 23), by our own Legendre symbols;
  2. every quadratic subfield's discriminant, and disc(k') as the product over the three quadratic subfields
     (conductor-discriminant formula for a biquadratic field): 2^8 (pm)^2, so rd(k') = 4 sqrt(pm);
  3. ramification and splitting of every prime l | 2pm in k, K and k', from the F_2-space of quadratic subfields
     (the inertia field is spanned by the subfields unramified at l, the decomposition field by those in which l
     splits); both sets are checked to be subgroups;
  4. K/k' unramified: e_l(K) = e_l(k') for every l | 2pm (all other primes are unramified in both), and all real;
  5. r = sum over l of the number of primes of k above l when e_l(K) = 2 e_l(k): must be >= 13;
     route 2 for the counts above odd l: Dedekind-Kummer, the factorization pattern of a defining polynomial of k
     modulo l, with the index of its order checked prime to l (discriminants by exact resultants);
  6. 4 sqrt(70035) enclosed by integer square roots; the printed decimals compared with it.

What rests on mathematics (read, not computed): section 2 of Martinet's paper (the 2-rank lower bound from genus
theory and the Golod-Shafarevich-Roquette criterion behind r >= n + 3 + 2 sqrt(2n+1)), the conductor-discriminant
formula, and that an unramified extension K of k' with an infinite tower forces an infinite tower over k'.

usage: python3 -I martinet_check.py [--controls]
"""
import math
import sys
from fractions import Fraction as Q
from itertools import product


def legendre(a, l):
    a %= l
    if a == 0:
        return 0
    return 1 if pow(a, (l - 1) // 2, l) == 1 else -1


def core(d):
    """squarefree part of a nonzero integer (sign kept)"""
    s, d = (-1 if d < 0 else 1), abs(d)
    out, f = 1, 2
    while f * f <= d:
        while d % (f * f) == 0:
            d //= f * f
        if d % f == 0:
            out *= f
            d //= f
        f += 1
    return s * out * d


def factor(n):
    fs, f = [], 2
    while f * f <= n:
        while n % f == 0:
            fs.append(f)
            n //= f
        f += 1
    if n > 1:
        fs.append(n)
    return fs


def quad_disc(d):
    d = core(d)
    return d if d % 4 == 1 else 4 * d


def quad_type(d, l):
    """'r' ramified, 's' split, 'i' inert, for the prime l in Q(sqrt d)"""
    D = quad_disc(d)
    if D % l == 0:
        return "r"
    if l == 2:
        return "s" if D % 8 == 1 else "i"
    return "s" if legendre(D, l) == 1 else "i"


def subfields(gens):
    """nonzero F_2 vectors v and the squarefree d_v of the quadratic subfield Q(sqrt(prod g_i^v_i))"""
    out = {}
    for v in product((0, 1), repeat=len(gens)):
        if any(v):
            x = 1
            for g, e in zip(gens, v):
                if e:
                    x *= g
            out[v] = core(x)
    return out


def is_subgroup(vs, s):
    z = tuple([0] * s)
    S = set(vs) | {z}
    return all(tuple(a ^ b for a, b in zip(x, y)) in S for x in S for y in S), len(S)


def efg(gens, l):
    """(e, f, g) of the prime l in the multiquadratic field Q(sqrt g : g in gens), from its quadratic subfields"""
    s = len(gens)
    sub = subfields(gens)
    unram = [v for v, d in sub.items() if quad_type(d, l) != "r"]
    split = [v for v, d in sub.items() if quad_type(d, l) == "s"]
    okU, nU = is_subgroup(unram, s)
    okS, nS = is_subgroup(split, s)
    assert okU and okS, f"subfield sets not subgroups at l = {l}"
    deg = 2 ** s
    e, g = deg // nU, nS
    f = nU // nS
    assert e * f * g == deg
    return e, f, g


# ---------------------------------------------------------------- route 2: Dedekind-Kummer for the counts in k
def pmod(a, m, l):
    a = [x % l for x in a]
    while len(a) >= len(m):
        if a[-1] == 0:
            a.pop()
            continue
        c = a[-1] * pow(m[-1], -1, l) % l
        sh = len(a) - len(m)
        for i, x in enumerate(m):
            a[sh + i] = (a[sh + i] - c * x) % l
        a.pop()
    while a and a[-1] == 0:
        a.pop()
    return a


def pmul(a, b, l):
    out = [0] * (len(a) + len(b) - 1) if a and b else []
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] = (out[i + j] + x * y) % l
    return out


def pgcd(a, b, l):
    while b:
        a, b = b, pmod(a, b, l)
    return a


def ppow_x(e, m, l):
    """x^e mod (m, l)"""
    result, base = [1], [0, 1]
    while e:
        if e & 1:
            result = pmod(pmul(result, base, l), m, l)
        base = pmod(pmul(base, base, l), m, l)
        e >>= 1
    return result


def factor_degrees(fpoly, l):
    """degrees of the irreducible factors of a squarefree monic polynomial mod l (distinct-degree factorization)"""
    f = [x % l for x in fpoly]
    degs, j = [], 1
    while len(f) > 1:
        xq = ppow_x(l ** j, f, l)
        h = xq + [0] * max(0, 2 - len(xq))
        h = h[:]
        h[1] = (h[1] - 1) % l
        while h and h[-1] == 0:
            h.pop()
        g = pgcd(f, h, l) if h else f
        dg = len(g) - 1
        if dg > 0:
            degs += [j] * (dg // j)
            # f = f / g
            q, r = pdivmod(f, g, l)
            assert not r
            f = q
        j += 1
        if j > 8:
            break
    return degs


def pdivmod(a, b, l):
    a = [x % l for x in a]
    q = [0] * max(0, len(a) - len(b) + 1)
    inv = pow(b[-1], -1, l)
    while len(a) >= len(b) and a:
        c = a[-1] * inv % l
        sh = len(a) - len(b)
        q[sh] = c
        for i, x in enumerate(b):
            a[sh + i] = (a[sh + i] - c * x) % l
        a.pop()
        while a and a[-1] == 0:
            a.pop()
    return q, a


def det(M):
    M = [[Q(x) for x in row] for row in M]
    n, d = len(M), Q(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None:
            return Q(0)
        if p != c:
            M[c], M[p] = M[p], M[c]
            d = -d
        d *= M[c][c]
        for r in range(c + 1, n):
            t = M[r][c] / M[c][c]
            for k in range(c, n):
                M[r][k] -= t * M[c][k]
    return d


def poly_disc(f):
    """discriminant of a monic polynomial (coefficients low to high) via the Sylvester resultant with f'"""
    n = len(f) - 1
    fp = [i * f[i] for i in range(1, n + 1)]
    hi_f, hi_fp = f[::-1], fp[::-1]
    size = 2 * n - 1
    rows = []
    for i in range(n - 1):
        rows.append([0] * i + hi_f + [0] * (size - i - len(hi_f)))
    for i in range(n):
        rows.append([0] * i + hi_fp + [0] * (size - i - len(hi_fp)))
    res = det(rows)
    return int((-1) ** (n * (n - 1) // 2) * res)


def main():
    controls = "--controls" in sys.argv
    status = run(29, (3, 5, 7, 23), verbose=True)
    if controls:
        print("\nred controls (each must make the criterion or a check fail):")
        c = []
        c.append(not run(29, (3, 5, 7, 19), verbose=False, label="q4 = 19 (19 = 3 mod 8)"))
        c.append(not run(29, (3, 5, 11, 23), verbose=False, label="q3 = 11 (11 = 3 mod 8)"))
        c.append(not run(13, (3, 5, 7, 23), verbose=False, label="p = 13"))
        c.append(not run(29, (3, 5, 7, 23), verbose=False, label="threshold r >= 14", extra=1))
        c.append(not run(29, (3, 7, 23, 31), verbose=False, label="(q) = (3, 7, 23, 31): m = 1 mod 4"))
        ok = check_printed("1058.565", upper=True, quiet=True) is False and \
            check_printed("1058.565066", upper=True, quiet=True) is False and \
            check_printed("1058.565067", upper=True, quiet=True) is True
        c.append(ok)
        print(f"  printed 1058.565 and 1058.565066 refused as upper bounds, 1058.565067 accepted: {ok}")
        print(f"red controls all as expected: {all(c)}")
        status = status and all(c)
    print(f"RESULT: {'PASS' if status else 'FAIL'}")
    sys.exit(0 if status else 1)


def check_printed(dec, upper, quiet=False):
    """exact comparison of a printed decimal with 4 sqrt(70035): returns True if it is a valid upper bound"""
    x = Q(dec)
    above = x > 0 and x * x > 16 * 70035
    if not quiet:
        print(f"  printed {dec}: {'above' if above else 'below'} 4 sqrt(70035), "
              f"{'a valid' if above else 'NOT a valid'} upper bound as printed")
    return above


def run(p, qs, verbose=True, label="Martinet's Example 4.2", extra=0):
    say = print if verbose else (lambda *a, **k: None)
    m = math.prod(qs)
    n = 4
    ok = True
    say(f"[{label}] p = {p}, (q1..q4) = {qs}, m = {m}, pm = {p * m} = {' . '.join(map(str, factor(p * m)))}")
    # 1. hypotheses of Proposition 4.1
    primes_ok = all(len(factor(x)) == 1 and x % 2 == 1 for x in (p,) + tuple(qs)) and len(set((p,) + tuple(qs))) == 5
    h1 = p % 4 == 1
    h2 = m % 4 == 3
    good = [q for q in qs if q % 8 in (1, 7) and legendre(q, p) == 1]
    h3 = len(good) >= 2
    say(f"[{label}] Prop. 4.1: five distinct odd primes {primes_ok}; (i) p = 1 mod 4 {h1}; (ii) m = -1 mod 4 {h2}; "
        f"(iii) q = +-1 mod 8 with (q/p) = +1: {good} (at least two: {h3})")
    hyp = primes_ok and h1 and h2 and h3
    # 2. discriminants
    kp_gens, k_gens, K_gens = [2, p * m], [2, p], [2, p, m]
    subs = subfields(kp_gens)
    discs = {d: quad_disc(d) for d in subs.values()}
    disc_kp = math.prod(abs(D) for D in discs.values())
    target = 2 ** 8 * (p * m) ** 2
    say(f"[{label}] quadratic subfields of k': " + ", ".join(f"Q(sqrt {d}) disc {D}" for d, D in sorted(discs.items()))
        + f"; product {disc_kp} = 2^8 (pm)^2: {disc_kp == target}")
    # rd(k') = disc^(1/4) = 4 sqrt(pm) exactly when disc = 2^8 (pm)^2
    ok &= disc_kp == target
    # 3-4. ramification
    ells = sorted(set([2] + factor(p * m)))
    unram = True
    r, rows = 0, []
    for l in ells:
        ek, fk, gk = efg(k_gens, l)
        eK, fK, gK = efg(K_gens, l)
        ekp, fkp, gkp = efg(kp_gens, l)
        unram &= (eK == ekp)
        ramifies = eK == 2 * ek
        if ramifies:
            r += gk
        rows.append((l, (ek, fk, gk), (eK, fK, gK), (ekp, fkp, gkp), ramifies))
    for l, a, b, c, rm in rows:
        say(f"[{label}]   l = {l:3d}: (e,f,g) in k {a}, in K {b}, in k' {c}; primes of k above l ramified in K/k: "
            f"{a[2] if rm else 0}")
    say(f"[{label}] K/k' unramified at every finite prime (e_l(K) = e_l(k') for all l | 2pm): {unram}; "
        f"all three fields totally real (generators positive): {all(g > 0 for g in K_gens)}")
    ok &= unram
    # route 2 for the counts in k above odd l: Dedekind-Kummer on a defining polynomial with index prime to l
    polys = {"alpha = sqrt 2 + sqrt p": alpha_poly(p), "beta = sqrt 2 + (1 + sqrt p)/2": beta_poly(p)}
    polys.update(small_elements(p))
    dk = math.prod(abs(quad_disc(d)) for d in subfields(k_gens).values())
    route2 = True
    for l in ells:
        if l == 2:
            continue
        ek, fk, gk = efg(k_gens, l)
        if ek != 1:
            continue
        used = None
        for name, f in polys.items():
            D = poly_disc(f)
            idx2 = Q(D, dk)
            if idx2.denominator == 1 and int(idx2) % l != 0:
                used = (name, f)
                break
        if used is None:
            say(f"[{label}]   l = {l}: no polynomial with index prime to l")
            route2 = False
            continue
        degs = factor_degrees(used[1], l)
        agree = len(degs) == gk and all(d == fk for d in degs)
        route2 &= agree
        say(f"[{label}]   route 2, l = {l}: {used[0]} mod {l} factors with degrees {degs}; agrees: {agree}")
    ok &= route2
    # 5. the criterion
    thr_ok = (r - n - 3 - extra) >= 0 and (r - n - 3 - extra) ** 2 >= 4 * (2 * n + 1)
    say(f"[{label}] r = {r}; criterion r >= n + 3 + 2 sqrt(2n+1) = {n + 3} + 2 sqrt({2 * n + 1})"
        + (f" (+{extra})" if extra else "") + f": {thr_ok}")
    ok &= thr_ok and hyp
    if verbose and p == 29 and tuple(qs) == (3, 5, 7, 23):
        # 6. the value
        N = 16 * 70035
        lo = math.isqrt(N * 10 ** 60)
        say(f"[{label}] 4 sqrt(70035) = sqrt({N}) in [{lo // 10**30}.{lo % 10**30:030d}, "
            f"{(lo + 1) // 10**30}.{(lo + 1) % 10**30:030d}]")
        say(f"[{label}] PR #332's check 1058.565^2 < 16.70035 < 1058.566^2: "
            f"{Q('1058.565') ** 2 < N < Q('1058.566') ** 2}")
        for dec in ("1058.565", "1058.565066", "1058.565067"):
            check_printed(dec, upper=True)
    if not verbose:
        print(f"  [{label}] hypotheses {hyp}, r = {r}, criterion {thr_ok}, disc {disc_kp == target}, "
              f"unramified {unram}: {'accepted' if ok else 'refused'}")
    return ok


def zw_mul(a, b, c):
    """(a0 + a1 w)(b0 + b1 w) with w^2 = w + c"""
    return (a[0] * b[0] + c * a[1] * b[1], a[0] * b[1] + a[1] * b[0] + a[1] * b[1])


def k_mul(u, v, c):
    """product in Z[sqrt 2, w], basis (1, s, w, s w), s^2 = 2, w = (1 + sqrt p)/2"""
    A, B, C, D = (u[0], u[2]), (u[1], u[3]), (v[0], v[2]), (v[1], v[3])
    ac, bd, ad, bc = zw_mul(A, C, c), zw_mul(B, D, c), zw_mul(A, D, c), zw_mul(B, C, c)
    re = (ac[0] + 2 * bd[0], ac[1] + 2 * bd[1])
    im = (ad[0] + bc[0], ad[1] + bc[1])
    return (re[0], im[0], re[1], im[1])


def charpoly(M):
    """characteristic polynomial (low to high, monic) by Faddeev-LeVerrier, exact"""
    n = len(M)
    M = [[Q(x) for x in r] for r in M]
    coeffs = [Q(1)]
    Mk = [[Q(0)] * n for _ in range(n)]
    I = [[Q(int(i == j)) for j in range(n)] for i in range(n)]
    ck = Q(1)
    for k in range(1, n + 1):
        Mk = [[sum(M[i][t] * (Mk[t][j] + ck * I[t][j]) for t in range(n)) for j in range(n)] for i in range(n)]
        ck = -sum(Mk[i][i] for i in range(n)) / k
        coeffs.append(ck)
    return [int(x) for x in coeffs[::-1]]


def small_elements(p):
    """defining polynomials of small elements theta of O_k = Z[sqrt 2, (1 + sqrt p)/2] (an integral basis, since
    the discriminant of this order is 8^2 p^2 = disc(k)), for route 2 at primes dividing the index of alpha"""
    c = (p - 1) // 4
    out = {}
    basis = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
    for t in product(range(-2, 3), repeat=3):
        th = (0,) + t
        cols = [k_mul(th, b, c) for b in basis]
        M = [[cols[j][i] for j in range(4)] for i in range(4)]
        out[f"theta = {t[0]} sqrt 2 + {t[1]} w + {t[2]} sqrt 2 w"] = charpoly(M)
    return out


def alpha_poly(p):
    # minimal polynomial of sqrt 2 + sqrt p: x^4 - 2(2+p) x^2 + (p-2)^2
    return [(p - 2) ** 2, 0, -2 * (2 + p), 0, 1]


def beta_poly(p):
    # w = (1 + sqrt p)/2, w^2 = w + (p-1)/4; beta = sqrt 2 + w satisfies, with c = (p-1)/4:
    # (b^2 + c - 2)^2 - (b^2 + c - 2)(2b - 1) - c (2b - 1)^2 = 0, expanded below
    c = (p - 1) // 4
    t = c - 2
    # (b^2 + t)^2 = b^4 + 2t b^2 + t^2
    # (b^2 + t)(2b - 1) = 2b^3 - b^2 + 2t b - t
    # c(2b - 1)^2 = 4c b^2 - 4c b + c
    return [t * t + t - c, -2 * t + 4 * c, 2 * t + 1 - 4 * c, -2, 1]


if __name__ == "__main__":
    main()
