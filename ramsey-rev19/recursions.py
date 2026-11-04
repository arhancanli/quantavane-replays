"""Which cells of the survey's tables the new lower bounds move, directly or through the survey's own recursions.

python3 -I recursions.py

Starts from the lower bounds printed in Small Ramsey Numbers, revision 18 (24 April 2026), Tables Ia, IIa, IIb and
IIc (exact values count as lower bounds), closes them under the lower-bound recursions of its section 2.3, then adds
the new bounds (colourings since revision 18, and Tatarevic's inequality) and closes again. Prints every cell of the
four tables whose bound rises above both the printed value and the revision-18 closure, with how it is obtained,
and every printed cell that the revision-18 data already beat.

Rules (R symmetric; all for the stated ranges only):
  e   R(k,l) >= R(k,l-1) + 2k - 3                    k, l >= 3                 [BEFS], 2.3.e
  h2  R(k,l) >= R(k,l-1) + 2k - 2                    k >= 5, l >= 3            2.3.h with p = 2 and k - 2 [XSR2]
  h   R(k,p+q-1) >= R(k,p) + R(k,q) + c              2 <= p <= q, k >= 3, c = k - 3 (k - 2 if k >= 5), and c = p - 2
  i   R(2k-1,l) >= 4 R(k,l-1) - 3                    l >= 5, k >= 2            [XXER], 2.3.i
  c   R(3,3k+1) >= 4 R(3,k+1) - 3 (k >= 2);  R(3,4k+1) >= 6 R(3,k+1) - 5 (k >= 1)     2.3.c
  T   R(k+1,s+1) >= R(k,s) + 2k + 2s                 5 <= k <= s               Tatarevic, arXiv:2608.06531 (new)
"""

KMAX, LMAX = 24, 60

# Revision 18, as printed (rendered pages 6, 10, 11). Table Ia and IIa disagree on R(4,15): Ia 159 [NaRT], IIa 158.
REV18 = {
    3: {3: 6, 4: 9, 5: 14, 6: 18, 7: 23, 8: 28, 9: 36, 10: 40, 11: 47, 12: 53, 13: 61, 14: 67, 15: 74, 16: 82,
        17: 92, 18: 100, 19: 106, 20: 111, 21: 122, 22: 131, 23: 139, 24: 143, 25: 154, 26: 161, 27: 172, 28: 179,
        29: 190, 30: 197, 31: 208, 32: 217, 33: 227, 34: 234, 35: 248, 36: 255, 37: 267, 38: 278, 39: 290, 40: 298,
        41: 311, 42: 320, 43: 333, 44: 339, 45: 354, 46: 362, 47: 380, 48: 384, 49: 402},
    4: {4: 18, 5: 25, 6: 36, 7: 49, 8: 59, 9: 73, 10: 92, 11: 102, 12: 128, 13: 139, 14: 148, 15: 159, 16: 170,
        17: 200, 18: 205, 19: 213, 20: 234, 21: 242, 22: 314},
    5: {5: 43, 6: 59, 7: 80, 8: 101, 9: 133, 10: 149, 11: 183, 12: 203, 13: 233, 14: 267, 15: 275, 16: 293,
        17: 388, 18: 396, 19: 411, 20: 424, 21: 441, 22: 492, 23: 521},
    6: {6: 102, 7: 115, 8: 134, 9: 183, 10: 204, 11: 262, 12: 294, 13: 347, 15: 401, 16: 434, 17: 548, 18: 614,
        19: 710, 20: 878, 21: 888, 22: 1070},
    7: {7: 205, 8: 219, 9: 252, 10: 292, 11: 405, 12: 417, 13: 511, 16: 629, 17: 729, 18: 797, 19: 908, 21: 1214},
    8: {8: 282, 9: 329, 10: 343, 11: 457, 13: 817, 15: 873, 17: 1005, 18: 1049, 19: 1237, 21: 1617},
    9: {9: 565, 10: 581},
    10: {10: 798, 15: 1313},
}
for k, v in zip(range(11, 24), [1597, 1640, 2557, 2989, 5485, 5605, 8917, 11005, 17885, 21725, 30925, 39109, 49421]):
    REV18.setdefault(k, {})[k] = v

# Cells the four tables show (printed or left blank).
SHOWN = {(k, l) for k in range(3, 11) for l in range(max(k, 3), 16)}
SHOWN |= {(k, l) for k in range(3, 9) for l in range(15, 24)}
SHOWN |= {(3, l) for l in range(24, 51)}
SHOWN |= {(k, k) for k in range(11, 25)}

# New since revision 18: (k, l, bound, source). Colourings unless marked.
NEW = [
    (3, 17, 93, 'NaRT repository, 22 Oct 2026'),
    (4, 15, 163, 'NaRT repository, 22 Oct 2026'),
    (4, 16, 174, 'NaRT, arXiv:2603.09172v3-v5 and repository'),
    (4, 18, 209, 'NaRT, arXiv:2603.09172v3-v5 and repository'),
    (4, 19, 219, 'NaRT, arXiv:2603.09172v3-v5 and repository'),
    (4, 20, 237, 'NaRT, arXiv:2603.09172v3-v5 and repository'),
    (4, 21, 244, 'NaRT repository, 11 May 2026'),
    (4, 20, 252, 'Yu, arXiv:2608.18169'),
    (5, 15, 279, 'Canli, Zenodo 23857057'),
    (6, 8, 135, 'Cremer, arXiv:2610.12122'),
    (7, 13, 513, 'NaRT repository, 22 Oct 2026'),
    (8, 10, 345, 'Cremer, arXiv:2610.12122'),
    (8, 12, 511, 'NaRT repository, 22 Oct 2026'),
    (8, 12, 519, 'Canli, Zenodo 23920916'),
    (8, 18, 1052, 'NaRT repository, 22 Oct 2026'),
    (9, 11, 617, 'NaRT repository, 22 Oct 2026'),
    (10, 12, 856, 'NaRT repository, 22 Oct 2026'),
    (12, 12, 1641, 'Tatarevic, arXiv:2608.06531 (his inequality from R(11,11) >= 1597; colouring in his repository)'),
]
CONIGLIO = {24: 151, 25: 160, 26: 163, 28: 181, 29: 194, 30: 199, 31: 210, 32: 219, 33: 230, 34: 239, 35: 253,
            36: 261, 37: 272, 38: 281, 39: 294, 40: 305, 41: 317, 42: 327, 43: 340, 44: 348, 45: 364, 46: 371,
            47: 385, 48: 391, 49: 411}
NEW += [(3, l, b, 'Coniglio et al., arXiv:2608.18769') for l, b in CONIGLIO.items()]


def key(a, b):
    return (a, b) if a <= b else (b, a)


def close(base, with_t):
    best = {}  # (k, l) -> (value, how)
    for (k, l), (v, how) in base.items():
        if v > best.get((k, l), (0, ''))[0]:
            best[(k, l)] = (v, how)
    for k in range(2, KMAX + 1):
        best.setdefault((2, k), (k, 'R(2,%d) = %d' % (k, k)))

    def get(a, b):
        return best.get(key(a, b), (0, ''))[0]

    def offer(a, b, v, how):
        a, b = key(a, b)
        if a < 2 or a > KMAX or b > LMAX:
            return False
        if v > get(a, b):
            best[(a, b)] = (v, how)
            return True
        return False

    changed = True
    while changed:
        changed = False
        for (a, b) in list(best):
            v = best[(a, b)][0]
            for k, l in ((a, b), (b, a)):  # k is the fixed parameter, l grows
                if k >= 3 and l + 1 >= 3:
                    changed |= offer(k, l + 1, v + 2 * k - 3, 'e: R(%d,%d) + %d' % (k, l, 2 * k - 3))
                if k >= 5 and l + 1 >= 3:
                    changed |= offer(k, l + 1, v + 2 * k - 2, 'h2: R(%d,%d) + %d' % (k, l, 2 * k - 2))
                if k >= 2 and l + 1 >= 5:
                    changed |= offer(2 * k - 1, l + 1, 4 * v - 3, 'i: 4 R(%d,%d) - 3' % (k, l))
                if with_t and 5 <= k <= l:
                    changed |= offer(k + 1, l + 1, v + 2 * k + 2 * l, 'T: R(%d,%d) + %d' % (k, l, 2 * k + 2 * l))
        for k in range(3, KMAX + 1):
            for p in range(2, LMAX + 1):
                for q in range(p, LMAX + 2 - p):
                    vp, vq = get(k, p), get(k, q)
                    if not vp or not vq:
                        continue
                    c = max(k - 2 if k >= 5 else k - 3, p - 2)
                    changed |= offer(k, p + q - 1, vp + vq + c, 'h: R(%d,%d) + R(%d,%d) + %d' % (k, p, k, q, c))
        for k in range(2, LMAX):
            v = get(3, k + 1)
            if v and k >= 2:
                changed |= offer(3, 3 * k + 1, 4 * v - 3, 'c: 4 R(3,%d) - 3' % (k + 1))
            if v and k >= 1:
                changed |= offer(3, 4 * k + 1, 6 * v - 5, 'c: 6 R(3,%d) - 5' % (k + 1))
    return best


def main():
    printed = {(k, l): v for k, row in REV18.items() for l, v in row.items()}
    base18 = {c: (v, 'printed') for c, v in printed.items()}
    c18 = close(base18, with_t=False)
    basenew = dict(base18)
    for k, l, b, src in NEW:
        if b > basenew.get((k, l), (0, ''))[0]:
            basenew[(k, l)] = (b, 'colouring: ' + src)
    cnew = close(basenew, with_t=True)
    print('Printed cells that the revision-18 data already beat through these recursions:')
    for c in sorted(SHOWN):
        if c in printed and c18[c][0] > printed[c]:
            print('  R(%d,%d): printed %d, implied %d by %s' % (c + (printed[c], c18[c][0], c18[c][1])))
    print('\nCells raised by the new bounds (printed value or, for a blank cell, the revision-18 closure):')
    for c in sorted(SHOWN):
        old = max(printed.get(c, 0), c18.get(c, (0, ''))[0])
        v, how = cnew.get(c, (0, ''))
        if v > old:
            tag = 'blank' if c not in printed else 'printed %d' % printed[c]
            print('  R(%d,%d) >= %d   (%s; rev-18 closure %d)   %s' % (c + (v, tag, c18.get(c, (0, ''))[0], how)))


if __name__ == '__main__':
    main()
