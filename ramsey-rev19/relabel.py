"""An isomorphic copy of a colouring with the vertices renumbered by a seeded random permutation.

python3 -I relabel.py MATRIX OUT [--seed N]

Clique numbers do not change under renumbering, but a circulant matrix stops looking circulant, so
independent/ramsey_indep then searches every vertex set instead of using the rotations. This lets both checkers
run on a circulant witness with no symmetry used. Prints the permutation's SHA-256 and the new matrix's.
"""
import argparse
import hashlib
import random


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('matrix')
    ap.add_argument('out')
    ap.add_argument('--seed', type=int, default=20261104)
    a = ap.parse_args()
    rows = [l.strip() for l in open(a.matrix, encoding='ascii') if l.strip()]
    n = len(rows)
    assert all(len(r) == n for r in rows)
    p = list(range(n))
    random.Random(a.seed).shuffle(p)  # new vertex i is old vertex p[i]
    text = ''.join(''.join(rows[p[i]][p[j]] for j in range(n)) + '\n' for i in range(n))
    open(a.out, 'w', encoding='ascii').write(text)
    new = text.split()
    circ = all(new[i][j] == new[0][(j - i) % n] for i in range(n) for j in range(n))
    print('matrix %s n %d seed %d permutation sha256 %s; wrote %s sha256 %s circulant %d'
          % (a.matrix, n, a.seed, hashlib.sha256(','.join(map(str, p)).encode()).hexdigest(), a.out,
             hashlib.sha256(text.encode('ascii')).hexdigest(), int(circ)))


if __name__ == '__main__':
    main()
