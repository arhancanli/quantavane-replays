#!/usr/bin/env python3
"""Finite-depth witness (method as described in words on page 4 of the note; own implementation):
chain of k = 2^43 + 1 copies of the 54-point block, v1 = (b+a, a+c), v_{k+1} = M v_k, M = [[a,b],[c,a]].
H0 = floor(D*M/q^52), H_{j+1} = floor(H_j^2 / D)  =>  H_j/D <= (M/q^52)^(2^j) entrywise.
Certificate: 1^T H43 v1 / (D q^54) > 9/8, hence |F_k| > q^(52k+2) for one finite family."""


def certify(a, b, c, qnum, qden, D=10 ** 80, sq=43):
    M = [[a, b], [c, a]]
    Q52n, Q52d = qnum ** 52, qden ** 52
    H = [[(D * M[i][j] * Q52d) // Q52n for j in range(2)] for i in range(2)]
    for _ in range(sq):
        H = [[(H[i][0] * H[0][j] + H[i][1] * H[1][j]) // D for j in range(2)] for i in range(2)]
    v1 = [b + a, a + c]
    s = sum(H[i][0] * v1[0] + H[i][1] * v1[1] for i in range(2))
    # s / (D q^54) > 9/8  <=>  8 s qden^54 > 9 D qnum^54
    return 8 * s * qden ** 54 > 9 * D * qnum ** 54


if __name__ == "__main__":
    cases = [
        ("Tilley, q=1.567866388204975", 7211663100, 4397311488, 11504366352, 1567866388204975, 10 ** 15),
        ("X* (one set), q=1.5678667956347", 7211751967, 4397440944, 11504366352, 15678667956347, 10 ** 13),
        ("X*+X' (two sets), q=1.567867203055880", 7211840834, 4397570400, 11504366352, 1567867203055880, 10 ** 15),
        ("control: one-set M with q=1.5678667956348 (above its rate; must fail)", 7211751967, 4397440944,
         11504366352, 15678667956348, 10 ** 13),
    ]
    for name, a, b, c, qn, qd in cases:
        print("%s: certificate holds = %s" % (name, certify(a, b, c, qn, qd)))
    print("points in the chain: 52*(2^43+1)+2 =", 52 * (2 ** 43 + 1) + 2)
    # two-set case: the 15-digit q needs a longer chain (2^50+1 blocks); a 13-digit q works at depth 2^43+1
    print("X*+X', q=1.5678672030558 (13 digits), 2^43+1 blocks:",
          certify(7211840834, 4397570400, 11504366352, 15678672030558, 10 ** 13))
    print("X*+X', q=1.567867203055880, 2^50+1 blocks:",
          certify(7211840834, 4397570400, 11504366352, 1567867203055880, 10 ** 15, sq=50))
