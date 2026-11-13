"""Program A: exact replay of a C_1c step-function certificate by one big-integer product (Kronecker substitution).
A_k = sum_i w_i w_{i+k}; bound = K * min_{0<=k<=K} A_k / S^2. Standard library only. Usage: c1c366_kron.py FILE [--bump J]"""
import json, sys
from fractions import Fraction
d = json.load(open(sys.argv[1])); K = d["K"]; w = [int(x) for x in d["weights_int"]]
assert all(x >= 0 for x in w) and len(w) > K
if "--bump" in sys.argv: w[int(sys.argv[sys.argv.index("--bump") + 1])] += 1
n = len(w); B = 2 * max(x.bit_length() for x in w) + n.bit_length() + 8
X = int.from_bytes(b"".join(x.to_bytes(B // 8 + 1, "little") for x in w), "little"); SL = 8 * (B // 8 + 1)
Y = int.from_bytes(b"".join(x.to_bytes(B // 8 + 1, "little") for x in reversed(w)), "little")
P = X * Y; mask = (1 << SL) - 1
A = [(P >> (SL * (n - 1 + k))) & mask for k in range(K + 1)]
S = sum(w); m = min(A); ks = [k for k, a in enumerate(A) if a == m]
assert A[0] == sum(x * x for x in w)
for k in ks + [1, K, K // 2]: assert A[k] == sum(w[i] * w[i + k] for i in range(n - k)), k   # direct re-sum
f = Fraction(K * m, S * S)
print("A: n", n, "K", K, "S", S, "min A_k", m, "at k", ks)
print("A: bound", f.numerator, "/", f.denominator, "=", float(f))
print("A: equals stated fraction:", f == Fraction(d["score_num"], d["score_den"]), "| stated argmin", d["argmin_k"], "in ours:", d["argmin_k"] in ks)
