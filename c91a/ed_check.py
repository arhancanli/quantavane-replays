"""91a: the lower bound C_91 >= 0.216276* of [J2026] (PR #273), the average edit distance of two random binary
strings, replayed from the argument in its README alone (Quantavane replays).  Standard library only.

The argument (re-derived here).  Every pair of strings has an optimal alignment (moves M, S, I, D) obeying five
local rules: (1) no I next to a D; (2) no S right after an I or a D; (3) the match ending a maximal run of
insertions (deletions) uses a letter different from every letter of the run; (4) no single match between a D and an
I (no D M I, no I M D); (5) a run of insertions (deletions) right after a substitution on (x_i, y_j) avoids the
letter x_i (y_j).  For random X, Y over k letters, the probability that a fixed move sequence obeying (1), (2), (4)
is valid and obeys (3) and (5) is a product of block factors: 1/k per match, (k-1)/k per substitution, and for a run
of r indels ended by a match (1/k)((k-1)/k)^r, or ((k-1)/k)^r (1/k) h(r) after a substitution, with
h(r) = 1/k + ((k-1)/k)((k-2)/(k-1))^r.  With x per letter and t per unit of cost, the expected number of such
alignments of cost <= a n is at most x^(-2n) t^(-a n) G(x, t), and G is finite when the transfer matrix A(x, t) of
an automaton with run lengths capped at R has spectral radius below 1 (capping only raises the weights).  So if
rho(A) < 1 and x^2 t^a > 1, then P(ed <= a n) -> 0 and C_91 >= a; the best a is -2 ln x / ln t.

What this program checks, with our own code (the repository's certify, count and canon programs were not read or
run; only the numbers x and t are taken from certify_rule5_R32.log):
  selftest   (i) for every binary pair up to the sizes listed, the least cost over alignments obeying (1)-(5) equals
             the edit distance (a dynamic programme over rule states); (ii) for every move sequence obeying (1), (2),
             (4) with at most 5 + 5 letters, the probability that it is valid and obeys (3) and (5), by enumeration
             of all letter assignments, equals the product of block factors; (iii) the automaton (uncapped), as a
             polynomial in x_X, x_Y and t, equals the sum of those products, and the brute-force expected counts.
  cert       A(x, t) for k = 2, R = 32 at the certificate's rational x and t: a rational positive v with
             (A v)_i < v_i for every state, checked exactly; independently, the run states eliminated in closed form
             (a 4 x 4 nonnegative matrix B, with rho(A) < 1 iff rho(B) < 1), certified the same way; then
             a* = -2 ln x / ln t in 60-digit decimal arithmetic with error bounds, compared with 0.216276.
  --controls x raised by 1% (rho >= 1, shown by a Collatz-Wielandt lower bound), rule 5 dropped (the automaton of
             rules 1-4 at the same x, t has rho >= 1), and the claim 0.216277 (above a*, not implied).

usage: python3 -I ed_check.py selftest | cert CERTS [--controls]
"""
import itertools
import os
import sys
import time
from decimal import Decimal, getcontext
from fractions import Fraction

K = 2


# ------------------------------------------------------------------------------------------------ the automaton
def automaton(k, R, rule5=True):
    """States and transitions.  Each transition: (src, dst, move, factor), factor a Fraction (the letter
    probability); move in M, S, I, D (I uses a Y letter, D an X letter, M and S one of each; S, I, D cost 1).
    Also end factors per state."""
    q = Fraction(k - 1, k)

    def h(r):
        return Fraction(1, k) + q * (Fraction(k - 2, k - 1) ** r if k > 1 else 0)

    st = ['N', 'Sx', 'MI', 'MD']
    for kind in ('I', 'Ip', 'D', 'Dp'):
        st += ['%s%d' % (kind, r) for r in range(1, R + 1)]
    idx = {s: i for i, s in enumerate(st)}
    tr = []
    for s in ('N', 'Sx', 'MI', 'MD'):
        tr.append((s, 'N', 'M', Fraction(1, k)))
        tr.append((s, 'Sx', 'S', q))
    tr += [('N', 'I1', 'I', 1), ('N', 'D1', 'D', 1), ('MI', 'I1', 'I', 1), ('MD', 'D1', 'D', 1)]
    if rule5:
        tr += [('Sx', 'Ip1', 'I', 1), ('Sx', 'Dp1', 'D', 1)]
    else:   # rules 1-4: a run after a substitution is an ordinary run
        tr += [('Sx', 'I1', 'I', 1), ('Sx', 'D1', 'D', 1)]
    end = {s: Fraction(1) for s in st}
    for r in range(1, R + 1):
        nxt = min(r + 1, R)
        for kind, mv, back in (('I', 'I', 'MI'), ('D', 'D', 'MD')):
            tr.append(('%s%d' % (kind, r), '%s%d' % (kind, nxt), mv, 1))
            tr.append(('%s%d' % (kind, r), back, 'M', Fraction(1, k) * q ** r))
            tr.append(('%sp%d' % (kind, r), '%sp%d' % (kind, nxt), mv, 1))
            tr.append(('%sp%d' % (kind, r), back, 'M', q ** r * Fraction(1, k) * h(r)))
            end['%sp%d' % (kind, r)] = q ** r
    return st, idx, tr, end


MOVE_LET = {'M': (1, 1, 0), 'S': (1, 1, 1), 'I': (0, 1, 1), 'D': (1, 0, 1)}   # (X letters, Y letters, cost)


def matrix(k, R, x, t, rule5=True):
    st, idx, tr, _ = automaton(k, R, rule5)
    n = len(st)
    A = [dict() for _ in range(n)]
    for s, d, mv, f in tr:
        lx, ly, c = MOVE_LET[mv]
        w = f * x ** (lx + ly) * t ** c
        A[idx[s]][idx[d]] = A[idx[s]].get(idx[d], 0) + w
    return st, A


# ------------------------------------------------------------------------------------------------- self-tests
def edit_distance(X, Y):
    n, m = len(X), len(Y)
    D = list(range(m + 1))
    for i in range(1, n + 1):
        prev, D[0] = D[0], i
        for j in range(1, m + 1):
            cur = min(D[j] + 1, D[j - 1] + 1, prev + (X[i - 1] != Y[j - 1]))
            prev, D[j] = D[j], cur
    return D[m]


def rule_min_cost(X, Y):
    """Least cost over alignments obeying rules 1-5, by a dynamic programme over (i, j, state).  State:
    ('N',) start or after an ordinary match; ('S', xi, yj) after a substitution; ('MI',) / ('MD',) after the match
    ending an insertion / deletion run; ('I', letters, avoid) / ('D', letters, avoid) inside a run, with the set of
    its letters and the letter it must avoid (rule 5; None if the run does not follow a substitution)."""
    n, m = len(X), len(Y)
    INF = 10 ** 9
    best = [[dict() for _ in range(m + 1)] for _ in range(n + 1)]
    best[0][0][('N',)] = 0

    def push(i, j, s, c):
        if c < best[i][j].get(s, INF):
            best[i][j][s] = c

    for i in range(n + 1):
        for j in range(m + 1):
            for s, c in list(best[i][j].items()):
                kind = s[0]
                inrun = kind in ('I', 'D')
                # M or S on (X[i], Y[j])
                if i < n and j < m:
                    a, b = X[i], Y[j]
                    if a == b:
                        if inrun:
                            if a not in s[1]:                          # rule 3
                                push(i + 1, j + 1, ('MI',) if kind == 'I' else ('MD',), c)
                        else:
                            push(i + 1, j + 1, ('N',), c)
                    elif not inrun:                                    # rule 2
                        push(i + 1, j + 1, ('S', a, b), c + 1)
                # I (uses Y[j])
                if j < m and kind != 'D' and kind != 'MD':             # rules 1, 4
                    y = Y[j]
                    if kind == 'I':
                        if s[2] is None or y != s[2]:
                            push(i, j + 1, ('I', s[1] | {y}, s[2]), c + 1)
                    else:
                        avoid = s[1] if kind == 'S' else None          # rule 5: avoid x_i
                        if avoid is None or y != avoid:
                            push(i, j + 1, ('I', frozenset([y]), avoid), c + 1)
                # D (uses X[i])
                if i < n and kind != 'I' and kind != 'MI':
                    xx = X[i]
                    if kind == 'D':
                        if s[2] is None or xx != s[2]:
                            push(i + 1, j, ('D', s[1] | {xx}, s[2]), c + 1)
                    else:
                        avoid = s[2] if kind == 'S' else None          # rule 5: avoid y_j
                        if avoid is None or xx != avoid:
                            push(i + 1, j, ('D', frozenset([xx]), avoid), c + 1)
    return min(best[n][m].values())


def seq_ok(seq):
    """rules 1, 2, 4 on a move sequence"""
    for p, c in zip(seq, seq[1:]):
        if (p, c) in (('I', 'D'), ('D', 'I'), ('I', 'S'), ('D', 'S')):
            return False
    for p, mm, c in zip(seq, seq[1:], seq[2:]):
        if mm == 'M' and ((p, c) == ('D', 'I') or (p, c) == ('I', 'D')):
            return False
    return True


def seq_probability(seq, k, n, m):
    """P over X in [k]^n, Y in [k]^m that seq is valid and obeys rules 3 and 5, by enumeration"""
    hits = 0
    for X in itertools.product(range(k), repeat=n):
        for Y in itertools.product(range(k), repeat=m):
            i = j = 0
            ok = True
            run, runlet, avoid, lastS = None, set(), None, None
            for mv in seq:
                if mv in ('M', 'S'):
                    a, b = X[i], Y[j]
                    if (mv == 'M') != (a == b):
                        ok = False
                        break
                    if mv == 'M' and run is not None and a in runlet:
                        ok = False
                        break
                    lastS = (a, b) if mv == 'S' else None
                    run, runlet, avoid = None, set(), None
                    i += 1
                    j += 1
                else:
                    letter = Y[j] if mv == 'I' else X[i]
                    if run is None:
                        run = mv
                        avoid = None if lastS is None else (lastS[0] if mv == 'I' else lastS[1])
                    if avoid is not None and letter == avoid:
                        ok = False
                        break
                    runlet.add(letter)
                    lastS = None
                    if mv == 'I':
                        j += 1
                    else:
                        i += 1
            hits += ok
    return Fraction(hits, k ** (n + m))


def block_product(seq, k):
    """the product of block factors for a sequence obeying rules 1, 2, 4"""
    q = Fraction(k - 1, k)
    h = lambda r: Fraction(1, k) + q * (Fraction(k - 2, k - 1) ** r if k > 1 else 0)
    p = Fraction(1)
    i = 0
    prevS = False
    while i < len(seq):
        mv = seq[i]
        if mv in ('M', 'S'):
            p *= Fraction(1, k) if mv == 'M' else q
            prevS = mv == 'S'
            i += 1
            continue
        r = 0
        while i < len(seq) and seq[i] == mv:
            r += 1
            i += 1
        if i < len(seq):                     # ended by a match (rules 1, 2 leave nothing else)
            assert seq[i] == 'M'
            p *= (q ** r * Fraction(1, k) * h(r)) if prevS else (Fraction(1, k) * q ** r)
            prevS = False
            i += 1
        else:
            p *= q ** r if prevS else 1
    return p


def automaton_poly(k, n, m):
    """sum over automaton paths (uncapped: R = n + m + 1) using n X letters and m Y letters: {cost: weight}"""
    st, idx, tr, end = automaton(k, n + m + 1)
    by_src = {}
    for s, d, mv, f in tr:
        by_src.setdefault(s, []).append((d, mv, f))
    out = {}
    frontier = {('N', 0, 0, 0): Fraction(1)}
    while frontier:
        nxt = {}
        for (s, a, b, c), w in frontier.items():
            if a == n and b == m:
                out[c] = out.get(c, 0) + w * end[s]
            for d, mv, f in by_src[s]:
                lx, ly, cc = MOVE_LET[mv]
                if a + lx <= n and b + ly <= m:
                    key = (d, a + lx, b + ly, c + cc)
                    nxt[key] = nxt.get(key, 0) + w * f
        frontier = nxt
    return out


def sequences(n, m):
    """all move sequences with n X letters and m Y letters obeying rules 1, 2, 4"""
    out = []

    def rec(seq, a, b):
        if a == n and b == m:
            out.append(tuple(seq))
            return
        for mv in 'MSID':
            lx, ly, _ = MOVE_LET[mv]
            if a + lx > n or b + ly > m:
                continue
            seq.append(mv)
            if seq_ok(seq[-3:]):
                rec(seq, a + lx, b + ly)
            seq.pop()

    rec([], 0, 0)
    return out


def selftest():
    k = K
    # (i) canonical form: the least rule-obeying cost equals the edit distance
    for n, m in ((8, 8), (7, 9), (6, 10), (4, 12)):
        bad = total = 0
        for X in itertools.product(range(k), repeat=n):
            for Y in itertools.product(range(k), repeat=m):
                total += 1
                if rule_min_cost(X, Y) != edit_distance(X, Y):
                    bad += 1
        print('canonical form (n, m) = (%d, %d): %d pairs, %d where the rule-obeying minimum differs from the edit '
              'distance' % (n, m, total, bad))
        if bad:
            raise SystemExit('selftest failed')
    # (ii) the probability of each move sequence, by enumeration of all letters, equals its block product
    for n, m in ((3, 3), (4, 4), (3, 5), (2, 6), (5, 2)):
        seqs = sequences(n, m)
        for seq in seqs:
            pr, bp = seq_probability(seq, k, n, m), block_product(seq, k)
            if pr != bp:
                raise SystemExit('selftest failed: %s: enumeration %s, product %s' % (''.join(seq), pr, bp))
        print('(n, m) = (%d, %d): %d move sequences obey rules 1, 2, 4; each probability (all %d letter assignments) '
              'equals its block product' % (n, m, len(seqs), k ** (n + m)))
    # (iii) the automaton, the sum of block products and the brute-force expected counts, by cost
    for n, m in ((5, 5), (4, 6), (6, 6), (7, 5), (3, 8)):
        seqs = sequences(n, m)
        poly = {}
        for seq in seqs:
            c = sum(MOVE_LET[v][2] for v in seq)
            poly[c] = poly.get(c, 0) + block_product(seq, k)
        auto = automaton_poly(k, n, m)
        brute = {}
        for X in itertools.product(range(k), repeat=n):
            for Y in itertools.product(range(k), repeat=m):
                for c, cnt in count_alignments(X, Y).items():
                    brute[c] = brute.get(c, 0) + Fraction(cnt, k ** (n + m))
        ok = poly == auto == brute
        print('(n, m) = (%d, %d): automaton, sum of block products (%d sequences) and brute-force expected counts by '
              'cost %s' % (n, m, len(seqs), 'agree' if ok else 'DISAGREE'))
        if not ok:
            raise SystemExit('selftest failed')


def count_alignments(X, Y):
    """number of alignments of X, Y obeying rules 1-5, by cost (dynamic programme with counts)"""
    n, m = len(X), len(Y)
    tab = [[dict() for _ in range(m + 1)] for _ in range(n + 1)]
    tab[0][0][('N',)] = {0: 1}

    def add(i, j, s, cmap, dc):
        d = tab[i][j].setdefault(s, {})
        for c, v in cmap.items():
            d[c + dc] = d.get(c + dc, 0) + v

    for i in range(n + 1):
        for j in range(m + 1):
            for s, cmap in list(tab[i][j].items()):
                kind = s[0]
                inrun = kind in ('I', 'D')
                if i < n and j < m:
                    a, b = X[i], Y[j]
                    if a == b:
                        if inrun:
                            if a not in s[1]:
                                add(i + 1, j + 1, ('MI',) if kind == 'I' else ('MD',), cmap, 0)
                        else:
                            add(i + 1, j + 1, ('N',), cmap, 0)
                    elif not inrun:
                        add(i + 1, j + 1, ('S', a, b), cmap, 1)
                if j < m and kind not in ('D', 'MD'):
                    y = Y[j]
                    if kind == 'I':
                        if s[2] is None or y != s[2]:
                            add(i, j + 1, ('I', s[1] | {y}, s[2]), cmap, 1)
                    else:
                        avoid = s[1] if kind == 'S' else None
                        if avoid is None or y != avoid:
                            add(i, j + 1, ('I', frozenset([y]), avoid), cmap, 1)
                if i < n and kind not in ('I', 'MI'):
                    xx = X[i]
                    if kind == 'D':
                        if s[2] is None or xx != s[2]:
                            add(i + 1, j, ('D', s[1] | {xx}, s[2]), cmap, 1)
                    else:
                        avoid = s[2] if kind == 'S' else None
                        if avoid is None or xx != avoid:
                            add(i + 1, j, ('D', frozenset([xx]), avoid), cmap, 1)
    out = {}
    for cmap in tab[n][m].values():
        for c, v in cmap.items():
            out[c] = out.get(c, 0) + v
    return out


# ---------------------------------------------------------------------------------------------------- the cert
def perron_vector(A, iters=20000):
    n = len(A)
    v = [1.0] * n
    for _ in range(iters):
        w = [sum(float(a) * v[j] for j, a in A[i].items()) for i in range(n)]
        s = max(w)
        w = [x / s + 1e-300 for x in w]
        if max(abs(a - b) for a, b in zip(v, w)) < 1e-16:
            v = w
            break
        v = w
    return v


def certify(A, label):
    """a rational v > 0 with (A v)_i < v_i for all i; returns the largest ratio (A v)_i / v_i (exact)"""
    n = len(A)
    # power iteration on A + I (same Perron vector, aperiodic)
    Ap = [dict(r) for r in A]
    for i in range(n):
        Ap[i][i] = Ap[i].get(i, 0) + 1
    v = perron_vector(Ap)
    vq = [Fraction(x) for x in v]
    worst = Fraction(0)
    for i in range(n):
        Av = sum(a * vq[j] for j, a in A[i].items())
        worst = max(worst, Av / vq[i])
    print('  %s: %d states, positive rational v; max_i (A v)_i / v_i = %.15f %s 1'
          % (label, n, float(worst), '<' if worst < 1 else '>='))
    return worst


def lower_rho(A):
    """Collatz-Wielandt: for v > 0, rho(A) >= min_i (A v)_i / v_i"""
    n = len(A)
    Ap = [dict(r) for r in A]
    for i in range(n):
        Ap[i][i] = Ap[i].get(i, 0) + 1
    v = perron_vector(Ap)
    vq = [Fraction(x) for x in v]
    return min(sum(a * vq[j] for j, a in A[i].items()) / vq[i] for i in range(n))


def reduced(k, R, x, t):
    """Eliminate the run states: B = A_dd + A_dr (I - A_rr)^(-1) A_rd on the four diagonal states, in closed form.
    An insertion run started from a diagonal state with weight x t continues with weight x t per letter; it is
    ended by a match with weight x^2 e(r) (e the end factor) and lands on MI.  With the cap at R the last state
    loops, so its tail is a geometric series."""
    q = Fraction(k - 1, k)
    h = lambda r: Fraction(1, k) + q * (Fraction(k - 2, k - 1) ** r if k > 1 else 0)
    u = x * t
    if u >= 1:
        raise SystemExit('x t >= 1')

    def run_sum(e):
        # sum_{r=1}^{R-1} u^r x^2 e(r) + u^R x^2 e(R) / (1 - u)
        s = sum(u ** r * x * x * e(r) for r in range(1, R))
        return s + u ** R * x * x * e(R) / (1 - u)

    plain = run_sum(lambda r: Fraction(1, k) * q ** r)
    afterS = run_sum(lambda r: q ** r * Fraction(1, k) * h(r))
    Mw, Sw = x * x * Fraction(1, k), x * x * t * q
    names = ['N', 'Sx', 'MI', 'MD']
    B = [dict() for _ in range(4)]
    for i in range(4):
        B[i][0] = Mw
        B[i][1] = Sw
    B[0][2] = B[0].get(2, 0) + plain      # N -> I-run -> MI
    B[0][3] = B[0].get(3, 0) + plain      # N -> D-run -> MD
    B[1][2] = B[1].get(2, 0) + afterS
    B[1][3] = B[1].get(3, 0) + afterS
    B[2][2] = B[2].get(2, 0) + plain      # MI -> I-run (I M I allowed), no D-run
    B[3][3] = B[3].get(3, 0) + plain
    return names, B


def a_star(x, t):
    getcontext().prec = 60
    lx = Decimal(x.numerator).ln() - Decimal(x.denominator).ln()
    lt = Decimal(t.numerator).ln() - Decimal(t.denominator).ln()
    err = Decimal(10) ** -55          # each ln correctly rounded at 60 digits; differences carry at most 2 ulps
    lo = (2 * (lx - err)) / (-(lt) + err) - Decimal(10) ** -50     # the division's own rounding, outward
    hi = (2 * (lx + err)) / (-(lt) - err) + Decimal(10) ** -50
    return lo, hi


def read_cert(path, k=2):
    for line in open(path):
        if line.startswith('k=%d ' % k):
            parts = dict(p.split('=', 1) for p in line.replace(':', ' ').split() if '=' in p)
            return Fraction(parts['x']), Fraction(parts['t'])
    raise SystemExit('k=%d line not found' % k)


def cert(root, controls):
    x, t = read_cert(os.path.join(root, 'j2026', 'certify_rule5_R32.log'))
    print('certificate (k = 2, R = 32): x = %s = %.12f, t = %s = %.12f, x t = %.6f' % (x, float(x), t, float(t), float(x * t)))
    st, A = matrix(K, 32, x, t)
    w1 = certify(A, 'rules 1-5, R = 32')
    names, B = reduced(K, 32, x, t)
    w2 = certify(B, 'run states eliminated (4 x 4)')
    lo, hi = a_star(x, t)
    claim = Decimal('0.216276')
    print('  a* = -2 ln x / ln t in [%s, %s]' % (str(lo)[:24], str(hi)[:24]))
    ok = w1 < 1 and w2 < 1 and x > 1 and 0 < t < 1 and lo > claim
    print('  x^2 t^0.216276 > 1 and rho < 1: C_91 >= 0.216276 %s (margin a* - 0.216276 >= %.6e)' % ('holds' if ok else 'NOT SHOWN', float(lo - claim)))
    if controls:
        x2 = x * Fraction(101, 100)
        _, A2 = matrix(K, 32, x2, t)
        lr = lower_rho(A2)
        print('  control: x raised by 1%%: rho(A) >= %.9f (Collatz-Wielandt), so no certificate exists %s'
              % (float(lr), '(as expected)' if lr >= 1 else 'UNEXPECTED'))
        _, A3 = matrix(K, 32, x, t, rule5=False)
        lr3 = lower_rho(A3)
        print('  control: rule 5 dropped (rules 1-4) at the same x, t: rho(A) >= %.9f %s'
              % (float(lr3), '(no certificate, as expected)' if lr3 >= 1 else 'UNEXPECTED'))
        bad = Decimal('0.216277')
        print('  control: the claim 0.216277 is %s' % ('not implied (above a*)' if bad > hi else 'IMPLIED'))
        ok = ok and lr >= 1 and lr3 >= 1 and bad > hi
    return ok


def main():
    t0 = time.time()
    if len(sys.argv) >= 2 and sys.argv[1] == 'selftest':
        selftest()
        print('RESULT selftest: canonical form, block products and automaton counts all agree')
    elif len(sys.argv) >= 3 and sys.argv[1] == 'cert':
        ok = cert(sys.argv[2], '--controls' in sys.argv)
        print('RESULT [J2026] C_91 >= 0.216276: %s' % ('PASS' if ok else 'FAIL'))
        if not ok:
            sys.exit(1)
    else:
        raise SystemExit(__doc__)
    print('run time %.0f s' % (time.time() - t0))


if __name__ == '__main__':
    main()
