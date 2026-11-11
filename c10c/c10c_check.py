"""Replay of the four certificate rows of the 10c page (Spencer's discrepancy constant), from the page alone.

The page defines C_10c as the least K with disc(A) <= K sqrt(n) for every n and every A in [-1,1]^{n x n}, where
disc(A) = min over x in {-1,1}^n of max_i |(Ax)_i|. So an n x n matrix A with entries in [-1,1] and disc(A) = d gives
C_10c >= d / sqrt(n). The page prints four such matrices:

  [G2026] (PR #74)  6 x 6 sign matrix, disc 4:      C_10c >= 4/sqrt(6)
  [X2026] (PR #86)  9 x 9 sign matrix, disc 5:      C_10c >= 5/3
  [L2026] (PR #129) 17 x 17 sign matrix, disc 7:    C_10c >= 7/sqrt(17)
  [H2026] (PR #170) 8 x 8 matrix over {-1,0,1}, disc 5 (column 8 is zero): C_10c >= 5/sqrt(8)

We read each matrix from the page's code block (as data), check its entries, and compute disc(A) exactly by two
exhaustive routes that share no code: route A walks all 2^n sign vectors in Gray-code order and updates Ax in exact
integers (one column change per step); route B enumerates every subset S of the columns (x_j = -1 on S) and evaluates
(Ax)_i = (|P_i| - |N_i|) - 2(|S & P_i| - |S & N_i|) with bit masks P_i, N_i of the +1 and -1 entries of row i. Both give
the full distribution of ||Ax||_inf. We also compare the page's T_i lists with the matrices, check its equality
witnesses, and compare each stated decimal with the exact bound (by squaring, in integers).

Red controls (each must come out as stated): every single-entry sign flip of each matrix (and, for [H2026], every
zero entry set to +1) lowers the discrepancy or leaves it, and the stated value plus one is never reached; [H2026] with
its zero column replaced by +1 (a sign matrix) has even discrepancy at most 4; a claimed discrepancy of d + 1 is refused.

Not checked: the page's remarks (the optimality of 5 at n = 8 by Kleitman's theorem, the covering-code framing and
K(n, R) values).

Usage: python3 -I c10c_check.py CERTS_ROOT [--controls]
"""
import hashlib
import os
import re
import sys
from fractions import Fraction

PAGE_SHA = "b8aa25a127d45239e25e78be9487f1b3bde4d62b87b05e24772f9b7c9e51365b"

# (key, PR, heading text on the page, n, claimed discrepancy, decimal printed in the table, exact form of the bound)
ROWS = [
    ("G2026", 74, "## Certificate for the $4/\\sqrt{6}$ lower bound", 6, 4, "1.632993", "4/sqrt(6)"),
    ("X2026", 86, "## Certificate for the $5/3$ lower bound", 9, 5, "1.666667", "5/3"),
    ("L2026", 129, "## Certificate for the $7/\\sqrt{17}$ lower bound", 17, 7, "1.697749", "7/sqrt(17)"),
    ("H2026", 170, "## Certificate for the $5/\\sqrt{8}$ lower bound", 8, 5, "1.767767", "5/sqrt(8)"),
]


def code_blocks(text):
    """all ``` blocks of the page, with their start offsets"""
    out = []
    for m in re.finditer(r"```[a-z]*\n(.*?)```", text, re.S):
        out.append((m.start(), m.group(1)))
    return out


def matrix_after(text, heading):
    """the first code block after the heading whose lines are all rows of integers of equal length"""
    pos = text.index(heading)
    for start, body in code_blocks(text):
        if start < pos:
            continue
        lines = [ln.split() for ln in body.strip().splitlines()]
        if lines and all(ln and all(re.fullmatch(r"-?\d+", t) for t in ln) for ln in lines):
            if len({len(ln) for ln in lines}) == 1:
                return [[int(t) for t in ln] for ln in lines], start
    raise ValueError("no matrix after " + heading)


def tsets_after(text, start, n):
    """the T_i lists printed after the matrix (T1 = {...}), as sets of column numbers 1..n"""
    for s, body in code_blocks(text):
        if s <= start:
            continue
        found = re.findall(r"T(\d+)\s*=\s*\{([0-9,\s]*)\}", body)
        if found:
            out = {}
            for i, members in found:
                out[int(i)] = {int(t) for t in members.replace(",", " ").split()}
            return out
    return None


def disc_route_a(A):
    """Gray-code walk over all x in {-1,1}^n; returns (min, histogram of ||Ax||_inf)"""
    n = len(A[0])
    m = len(A)
    cols = [[A[i][j] for i in range(m)] for j in range(n)]
    x = [1] * n
    y = [sum(row) for row in A]            # x = all ones
    hist = {}
    best = None
    for step in range(1 << n):
        v = max(abs(t) for t in y)
        hist[v] = hist.get(v, 0) + 1
        best = v if best is None or v < best else best
        if step == (1 << n) - 1:
            break
        j = ((step + 1) & -(step + 1)).bit_length() - 1   # the bit that changes in the Gray code
        if x[j] == 1:
            x[j] = -1
            c = cols[j]
            for i in range(m):
                y[i] -= 2 * c[i]
        else:
            x[j] = 1
            c = cols[j]
            for i in range(m):
                y[i] += 2 * c[i]
    return best, hist


def disc_route_b(A):
    """every subset S of columns (x_j = -1 on S), row values from bit masks and popcounts"""
    n = len(A[0])
    P = [sum(1 << j for j in range(n) if row[j] == 1) for row in A]
    N = [sum(1 << j for j in range(n) if row[j] == -1) for row in A]
    base = [bin(p).count("1") - bin(q).count("1") for p, q in zip(P, N)]
    hist = {}
    argmins = []
    best = None
    for S in range(1 << n):
        v = 0
        for b, p, q in zip(base, P, N):
            t = b - 2 * (bin(S & p).count("1") - bin(S & q).count("1"))
            if t < 0:
                t = -t
            if t > v:
                v = t
        hist[v] = hist.get(v, 0) + 1
        if best is None or v < best:
            best, argmins = v, [S]
        elif v == best and len(argmins) < 3:
            argmins.append(S)
    return best, hist, argmins


def bound_vs_decimal(d, n, dec):
    """compare the decimal string with d/sqrt(n): returns 'below', 'equal' or 'above' (exact, by squaring)"""
    q = Fraction(dec)
    lhs = q * q * n
    rhs = Fraction(d * d)
    return "below" if lhs < rhs else ("equal" if lhs == rhs else "above")


def digits(d, n, k=12):
    """d/sqrt(n) truncated to k decimals, exactly (integer square root)"""
    from math import isqrt
    # floor(d * 10^k / sqrt(n)) = floor(sqrt(d^2 * 10^(2k) / n))
    v = isqrt((d * d * 10 ** (2 * k)) // n)
    s = str(v).rjust(k + 1, "0")
    return s[:-k] + "." + s[-k:]


def main():
    root = sys.argv[1]
    controls = "--controls" in sys.argv[2:]
    path = os.path.join(root, "page", "10c.md")
    raw = open(path, "rb").read()
    sha = hashlib.sha256(raw).hexdigest()
    print("10c page (teorth/optimizationproblems @ 87b4ebb, constants/10c.md) sha256", sha)
    if sha != PAGE_SHA:
        print("RESULT: FAIL (the page is not the pinned file)")
        return 1
    text = raw.decode("utf-8")
    ok_all = True
    mats = {}
    for key, pr, heading, n, d, dec, exact in ROWS:
        A, start = matrix_after(text, heading)
        mats[key] = A
        print(f"\n[{key}] (PR #{pr}): C_10c >= {exact}, claimed disc(A) = {d} for an {n} x {n} matrix")
        shape = len(A) == n and all(len(r) == n for r in A)
        entries = sorted({v for r in A for v in r})
        in_range = all(-1 <= v <= 1 for v in entries)
        print(f"  matrix read from the page: {len(A)} x {len(A[0])}, entries {entries}; square {n} x {n}: {shape}; "
              f"all entries in [-1, 1]: {in_range}")
        zero_cols = [j + 1 for j in range(n) if all(A[i][j] == 0 for i in range(n))]
        dup = [(i + 1, k + 1) for i in range(n) for k in range(i + 1, n) if A[i] == A[k]]
        print(f"  identically zero columns: {zero_cols or 'none'}; repeated rows: {dup or 'none'}")
        T = tsets_after(text, start, n)
        if T is not None:
            want = {i + 1: {j + 1 for j in range(n) if A[i][j] == -1} for i in range(n)}
            tok = all(T.get(i) == want[i] for i in want) and len(T) == n
            print(f"  the page's T_i lists equal the -1 positions of the rows: {tok}")
            ok_all &= tok
        ra, ha = disc_route_a(A)
        rb, hb, arg = disc_route_b(A)
        same = (ra == rb) and (ha == hb)
        hs = ", ".join(f"{v}: {c}" for v, c in sorted(ha.items()))
        print(f"  route A (Gray code, 2^{n} = {1 << n} vectors): disc = {ra}; route B (bit masks): disc = {rb}; "
              f"histograms identical: {same}")
        print(f"  distribution of ||Ax||_inf over all sign vectors: {hs}")
        xs = []
        for S in arg[:1]:
            x = [(-1 if (S >> j) & 1 else 1) for j in range(n)]
            y = [sum(A[i][j] * x[j] for j in range(n)) for i in range(n)]
            xs.append((x, y))
            print(f"  a minimizer: x = {tuple(x)}, Ax = {tuple(y)}")
        good = shape and in_range and same and ra == d
        rel = bound_vs_decimal(d, n, dec)
        print(f"  disc(A) = {d} as claimed: {ra == d}; bound {exact} = {digits(d, n)}...; the table's decimal {dec} is "
              f"{rel} the exact value")
        ok_all &= good
    # the page's equality witnesses
    print("\nthe page's equality examples:")
    ex = [("G2026", (1, -1, 1, 1, 1, 1), (4, 2, 0, -2, 0, 2)),
          ("X2026", (-1, -1, 1, -1, 1, 1, 1, 1, 1), (3, 5, 5, 3, 5, 3, 3, -3, 3)),
          ("H2026", (-1, -1, -1, -1, -1, -1, 1, -1), None)]
    for key, x, want in ex:
        A = mats[key]
        y = tuple(sum(A[i][j] * x[j] for j in range(len(x))) for i in range(len(A)))
        v = max(abs(t) for t in y)
        print(f"  [{key}] x = {x}: Ax = {y}, ||Ax||_inf = {v}" + (f"; as printed: {y == want}" if want else ""))
        if want is not None:
            ok_all &= (y == want)
    # H2026's stated distribution: 5 (224 times) and 7 (32 times)
    _, hh, _ = disc_route_b(mats["H2026"])
    hd = hh == {5: 224, 7: 32}
    print(f"  [H2026] the stated distribution 5 (224 times) and 7 (32 times): {hd}")
    ok_all &= hd
    cok = True
    if controls:
        print("\ncontrols (each must come out as stated):")
        for key, pr, heading, n, d, dec, exact in ROWS:
            A = mats[key]
            # a corrupted certificate: the first single sign flip (row-major) whose discrepancy is below d
            low = None
            tried = 0
            for i in range(n):
                for j in range(n):
                    if A[i][j] == 0:
                        continue
                    B = [r[:] for r in A]
                    B[i][j] = -A[i][j]
                    tried += 1
                    va = disc_route_a(B)[0]
                    vb = disc_route_b(B)[0]
                    if va != vb:
                        cok = False
                        print(f"[ctl] [{key}] the two routes disagree on a flipped matrix: {va} vs {vb}")
                    if va < d:
                        low = (i + 1, j + 1, va)
                        break
                if low:
                    break
            c = low is not None
            print(f"[ctl] [{key}] the entry at row {low[0]}, column {low[1]} flipped: disc {low[2]} < {d}, so the "
                  f"corrupted certificate is refused (both routes; flip {tried} tried): {c}"
                  if c else f"[ctl] [{key}] no single flip lowers the discrepancy: {c}")
            cok &= c
            c = disc_route_a(A)[0] < d + 1
            print(f"[ctl] [{key}] a claimed discrepancy {d + 1} is refused: {c}")
            cok &= c
        # both routes against a naive triple loop on random small matrices over {-1, 0, 1}
        import random
        rng = random.Random(10)
        agree = True
        for _ in range(60):
            n = rng.randint(1, 7)
            B = [[rng.choice((-1, 0, 1)) for _ in range(n)] for _ in range(n)]
            naive = None
            for S in range(1 << n):
                x = [(-1 if (S >> j) & 1 else 1) for j in range(n)]
                v = max(abs(sum(B[i][j] * x[j] for j in range(n))) for i in range(n))
                naive = v if naive is None or v < naive else naive
            agree &= disc_route_a(B)[0] == naive == disc_route_b(B)[0]
        print(f"[ctl] routes A and B equal a naive triple loop on 60 random matrices over {{-1, 0, 1}}, n <= 7: {agree}")
        cok &= agree
        A = [r[:] for r in mats["H2026"]]
        for r in A:
            r[7] = 1
        v, h, _ = disc_route_b(A)
        c = v <= 4 and v % 2 == 0
        print(f"[ctl] [H2026] with column 8 set to +1 (a sign matrix at even n): disc {v}, even and at most 4: {c}")
        cok &= c
        print(f"controls all as expected: {cok}")
    ok = ok_all and cok
    print(f"\nRESULT: {'PASS' if ok else 'FAIL'} (disc = 4, 5, 7, 5 for the four matrices, as the page states; "
          f"C_10c >= 4/sqrt(6), 5/3, 7/sqrt(17), 5/sqrt(8))")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
