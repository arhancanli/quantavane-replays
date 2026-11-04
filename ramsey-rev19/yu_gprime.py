"""Yu's second graph G' (arXiv:2608.18169v2, section 4): the image of G = C_251(S) under x -> 116x (mod 251).

python3 -I yu_gprime.py G_MATRIX OUT
Builds G' = C_251(116 S) from S as read off G's row 0, checks that 116 S lies in D0 u D3 (the paper's Table 1) and
that x -> 116x maps G onto G' edge for edge, and writes G' as n lines of n characters.
"""
import hashlib
import sys

P, MULT = 251, 116
D0 = {1, 2, 4, 5, 8, 10, 16, 20, 25, 32, 40, 47, 50, 51, 63, 64, 69, 80, 91, 94, 100, 102, 113, 123, 125}
D3 = {7, 14, 19, 27, 28, 29, 35, 38, 39, 53, 54, 56, 58, 61, 70, 76, 78, 95, 99, 106, 108, 111, 112, 116, 122}


def fold(d):
    d %= P
    return min(d, P - d)


def main():
    g = [l.strip() for l in open(sys.argv[1], encoding='ascii') if l.strip()]
    assert len(g) == P and all(len(r) == P for r in g)
    s = sorted(d for d in range(1, P // 2 + 1) if g[0][d] == '1')
    s2 = sorted({fold(MULT * d) for d in s})
    assert len(s2) == len(s) == 32 and set(s2) <= D0 | D3
    rows = [''.join('1' if i != j and fold(i - j) in set(s2) else '0' for j in range(P)) for i in range(P)]
    assert all(g[i][j] == rows[MULT * i % P][MULT * j % P] for i in range(P) for j in range(P))
    text = ''.join(r + '\n' for r in rows)
    open(sys.argv[2], 'w', encoding='ascii').write(text)
    print('S  = %s' % s)
    print("S' = %s (in D0 u D3: %d of D0, %d of D3); x -> 116x maps G onto G'" %
          (s2, len(set(s2) & D0), len(set(s2) & D3)))
    print("wrote %s sha256 %s" % (sys.argv[2], hashlib.sha256(text.encode('ascii')).hexdigest()))


if __name__ == '__main__':
    main()
