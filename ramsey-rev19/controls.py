"""Two red controls for a checked colouring: each flips exactly one edge so that a forbidden clique must appear.

python3 -I controls.py MATRIX S T OUT_RED OUT_BLUE [--seed N]

MATRIX: n lines of n characters, '1' the colour that must avoid K_S, '0' the colour that must avoid K_T.
OUT_RED : one '0' entry (and its mirror) turned into '1' so that the '1' graph holds a K_S.
OUT_BLUE: one '1' entry (and its mirror) turned into '0' so that the '0' graph holds a K_T.
How: find a clique C of size S-2 in the colour whose common neighbourhood N in that colour has at least two vertices.
N holds no edge of that colour (C plus such an edge would already be a K_S), so any two vertices u, v of N are joined
in the other colour, and flipping {u, v} makes C + u + v a clique of size S. The new clique is printed and re-checked
directly, pair by pair, in the written file. A checker passes the control only if it reports the clique.
Standard library only; the randomised search uses a fixed seed.
"""
import argparse
import hashlib
import random
import sys


def read(path):
    rows = [l.strip() for l in open(path, encoding='ascii') if l.strip()]
    n = len(rows)
    assert all(len(r) == n for r in rows)
    return rows


def nbr_masks(rows, ch):
    n = len(rows)
    return [sum(1 << j for j in range(n) if j != i and rows[i][j] == ch) for i in range(n)]


def find_seed_clique(adj, size, rng, tries=20000):
    """A clique of `size` vertices whose common neighbourhood has >= 2 vertices (randomised greedy, then DFS)."""
    n = len(adj)
    full = (1 << n) - 1
    for _ in range(tries):
        clique, cand = [], full
        while len(clique) < size and cand:
            # among candidates prefer vertices keeping many candidates; random tie-break
            opts = []
            c = cand
            while c:
                b = c & -c
                v = b.bit_length() - 1
                c ^= b
                opts.append((bin(cand & adj[v]).count('1'), rng.random(), v))
            opts.sort(reverse=True)
            v = opts[min(len(opts) - 1, int(rng.random() ** 3 * min(len(opts), 8)))][2]
            clique.append(v)
            cand &= adj[v]
        if len(clique) == size and bin(cand).count('1') >= 2:
            return clique, cand
    return None, 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('matrix')
    ap.add_argument('s', type=int)
    ap.add_argument('t', type=int)
    ap.add_argument('out_red')
    ap.add_argument('out_blue')
    ap.add_argument('--seed', type=int, default=20261104)
    a = ap.parse_args()
    rows = read(a.matrix)
    n = len(rows)
    rng = random.Random(a.seed)
    print('matrix %s n %d sha256 %s' % (a.matrix, n, hashlib.sha256(open(a.matrix, 'rb').read()).hexdigest()))
    for ch, size, out, name in (('1', a.s, a.out_red, 'red'), ('0', a.t, a.out_blue, 'blue')):
        adj = nbr_masks(rows, ch)
        clique, common = find_seed_clique(adj, size - 2, rng)
        if clique is None:
            print('%s control: no seed clique found' % name)
            sys.exit(1)
        members = [v for v in range(n) if common >> v & 1]
        u, v = members[0], members[1]
        assert rows[u][v] != ch, 'the colouring already holds a forbidden clique'
        new = [list(r) for r in rows]
        new[u][v] = new[v][u] = ch
        big = sorted(clique + [u, v])
        assert all(new[x][y] == ch for x in big for y in big if x != y)
        text = ''.join(''.join(r) + '\n' for r in new)
        open(out, 'w', encoding='ascii').write(text)
        print('%s control: flip {%d,%d} from %s to %s; it completes a K_%d on %s; wrote %s sha256 %s'
              % (name, u, v, '0' if ch == '1' else '1', ch, size, ','.join(map(str, big)), out,
                 hashlib.sha256(text.encode('ascii')).hexdigest()))


if __name__ == '__main__':
    main()
