"""Turn a published Ramsey witness into the plain matrix the two checkers read, refusing anything malformed.

python3 -I to_matrix.py FORMAT IN OUT --n N [--complement] [--label NAME]

FORMAT is one of
  numpy      NaRT (google-research ramsey_number_bounds): a printed NumPy array, [[0 1 ...] [...]]; each inner
             list must hold exactly N integers 0/1.
  dimacs     Coniglio et al.: "p edge N M" and M lines "e u v" (vertices 1..N); every edge once, no loops.
  rows       N lines of N characters 0/1 (Cremer, Tatarevic, Quantavane).
  circulant  Yu: the connection set in braces after "NAME =" (--label), on Z_N; i ~ j when the circular distance
             min(|i-j|, N-|i-j|) is in the set.
OUT gets N lines of N characters, '1' for an edge of the input graph and '0' otherwise (swapped by --complement).
The input must be square, 0/1, symmetric, with a zero diagonal, on exactly N vertices; otherwise nothing is
written and the exit status is 1. Prints the SHA-256 of IN and OUT, N, the number of 1s, the degree range and
whether the matrix is circulant (row i is row 0 turned by i). Standard library only; nothing in IN is executed.
"""
import argparse
import hashlib
import re
import sys


def fail(msg):
    print('REFUSED: ' + msg)
    sys.exit(1)


def parse_numpy(text, n):
    if text.count('[') != n + 1 or text.count(']') != n + 1:
        fail('expected %d bracketed rows inside one outer list' % n)
    body = text.strip()
    if not (body.startswith('[[') and body.endswith(']]')):
        fail('not a printed 2-d array')
    rows = []
    for chunk in re.findall(r'\[([^\[\]]*)\]', body):
        vals = chunk.split()
        if len(vals) != n or any(v not in ('0', '1') for v in vals):
            fail('a row does not hold %d entries 0/1' % n)
        rows.append([int(v) for v in vals])
    if len(rows) != n:
        fail('%d rows, expected %d' % (len(rows), n))
    return rows


def parse_dimacs(text, n):
    a = [[0] * n for _ in range(n)]
    declared = None
    seen = 0
    for line in text.splitlines():
        f = line.split()
        if not f or f[0] == 'c':
            continue
        if f[0] == 'p':
            if declared is not None or len(f) != 4 or f[1] != 'edge':
                fail('bad problem line: %r' % line)
            if int(f[2]) != n:
                fail('problem line says %s vertices, expected %d' % (f[2], n))
            declared = int(f[3])
        elif f[0] == 'e':
            if declared is None or len(f) != 3:
                fail('edge before problem line or bad edge line: %r' % line)
            u, v = int(f[1]) - 1, int(f[2]) - 1
            if not (0 <= u < n and 0 <= v < n) or u == v:
                fail('bad edge %r' % line)
            if a[u][v]:
                fail('edge listed twice: %r' % line)
            a[u][v] = a[v][u] = 1
            seen += 1
        else:
            fail('unknown line: %r' % line)
    if declared is None or seen != declared:
        fail('problem line declares %s edges, file lists %d' % (declared, seen))
    return a


def parse_rows(text, n):
    lines = [l.strip() for l in text.splitlines() if l.strip()]
    if len(lines) != n:
        fail('%d non-empty lines, expected %d' % (len(lines), n))
    rows = []
    for l in lines:
        if len(l) != n or set(l) - {'0', '1'}:
            fail('a line is not %d characters 0/1' % n)
        rows.append([int(c) for c in l])
    return rows


def parse_circulant(text, n, label):
    m = re.search(re.escape(label) + r'\s*=\s*\{([0-9,\s]*)\}', text)
    if not m:
        fail('no "%s = {...}" in the file' % label)
    s = [int(x) for x in m.group(1).replace(',', ' ').split()]
    if len(set(s)) != len(s) or any(not (1 <= d <= n // 2) for d in s):
        fail('connection set has repeats or distances outside 1..N/2')
    ds = set(s)
    print('connection set (%s): %d distances, %s' % (label, len(s), ','.join(map(str, s))))
    return [[1 if i != j and min(abs(i - j), n - abs(i - j)) in ds else 0 for j in range(n)] for i in range(n)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('format', choices=['numpy', 'dimacs', 'rows', 'circulant'])
    ap.add_argument('inp')
    ap.add_argument('out')
    ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--complement', action='store_true')
    ap.add_argument('--label', default='S')
    a = ap.parse_args()
    raw = open(a.inp, 'rb').read()
    text = raw.decode('utf-8')  # Yu's certificate has prose with non-ASCII dashes; matrix formats are checked below
    n = a.n
    rows = {'numpy': lambda: parse_numpy(text, n), 'dimacs': lambda: parse_dimacs(text, n),
            'rows': lambda: parse_rows(text, n), 'circulant': lambda: parse_circulant(text, n, a.label)}[a.format]()
    for i in range(n):
        if rows[i][i]:
            fail('1 on the diagonal at %d' % i)
        for j in range(i + 1, n):
            if rows[i][j] != rows[j][i]:
                fail('not symmetric at %d,%d' % (i, j))
    if a.complement:
        rows = [[0 if i == j else 1 - rows[i][j] for j in range(n)] for i in range(n)]
    out = ''.join(''.join(map(str, r)) + '\n' for r in rows).encode('ascii')
    circulant = all(rows[i][j] == rows[0][(j - i) % n] for i in range(n) for j in range(n))
    deg = [sum(r) for r in rows]
    with open(a.out, 'wb') as f:
        f.write(out)
    print('in  %s sha256 %s (%d bytes, format %s)' % (a.inp, hashlib.sha256(raw).hexdigest(), len(raw), a.format))
    print('out %s sha256 %s' % (a.out, hashlib.sha256(out).hexdigest()))
    print('n %d ones %d degree %d..%d circulant %d complemented %d'
          % (n, sum(deg) // 2, min(deg), max(deg), int(circulant), int(a.complement)))


if __name__ == '__main__':
    main()
