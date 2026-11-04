"""The colour-c neighbourhood of vertex 0 in a circulant colouring, as a matrix the checkers read.

python3 -I nbhd.py MATRIX C OUT

Refuses unless MATRIX is circulant (row i equals row 0 turned by i), symmetric, with a zero diagonal. Translations
i -> i + 1 (mod n) then preserve both colours and act transitively, so a colour-c clique of size k exists iff one
contains vertex 0, iff the colour-c neighbourhood N_c(0) holds a colour-c clique of size k - 1. OUT is the
colouring induced on N_c(0) (vertices in increasing order). Standard library only.
"""
import hashlib
import sys


def main():
    path, c, out = sys.argv[1], sys.argv[2], sys.argv[3]
    rows = [l.strip() for l in open(path, encoding='ascii') if l.strip()]
    n = len(rows)
    if any(len(r) != n or set(r) - {'0', '1'} for r in rows):
        sys.exit('REFUSED: not a square 0/1 matrix')
    for i in range(n):
        if rows[i][i] != '0':
            sys.exit('REFUSED: diagonal')
        for j in range(n):
            if rows[i][j] != rows[j][i] or rows[i][j] != rows[0][(j - i) % n]:
                sys.exit('REFUSED: not symmetric and circulant at %d,%d' % (i, j))
    keep = [j for j in range(n) if j and rows[0][j] == c]
    text = ''.join(''.join(rows[a][b] for b in keep) + '\n' for a in keep)
    open(out, 'w', encoding='ascii').write(text)
    print('matrix %s n %d circulant 1; N_%s(0) has %d vertices; wrote %s sha256 %s'
          % (path, n, c, len(keep), out, hashlib.sha256(text.encode('ascii')).hexdigest()))


if __name__ == '__main__':
    main()
