# 52a, random 3-SAT: the first-moment bound [FP1983] 5.191 from the page's model (15 November 2026)

A value recomputed from its definition (the row's own argument), not a certificate. The page's model: a formula of
density $r = m/n$ is $m$ clauses chosen uniformly at random, with replacement, from the $8\binom n3$ clauses on three
distinct variables. For a fixed assignment exactly one of the 8 sign patterns on each triple is violated, so
$\mathbb E[\#\text{satisfying assignments}] = 2^n(7/8)^m$ and $\Pr[\text{satisfiable}] \le (2(7/8)^r)^n \to 0$ for
$r > r_1 = \log 2/\log(8/7)$. The row ("Direct first moment method") is $r_3 \le 5.191$.

- `c52a_a.py`: over ALL $(8\binom n3)^m$ formulas, the number of (formula, satisfying assignment) pairs equals
  $(8\binom n3)^m\, 2^n (7/8)^m$ exactly, for $(n, m) = (3, 1..3), (4, 1..3), (5, 1..2)$; then $r_1 < 5.191$ in integers:
  $2(7/8)^{5191/1000} < 1 \iff 8^{5191} > 2^{1000}\,7^{5191}$; and $r_1 > 5.19$.
- `c52a_b.py`: every assignment of $n \le 9$ variables satisfies exactly $7/8$ of the clauses, so
  $\sum_a (s_a/N)^m = 2^n(7/8)^m$ for every $m$; $r_1 \in [5.19089306968443160590833221439893720709\ldots]$ in decimal
  intervals.

| statement | both programs |
|---|---|
| $\mathbb E[\#\text{satisfying assignments}] = 2^n (7/8)^m$ in the page's model | yes (small cases, exhaustive) |
| $r_1 = \log 2/\log(8/7) = 5.1908930696\ldots < 5.191$ | yes |

**Red controls**, each refused by both programs: the bound claimed as 5.190; clauses on two variables (each satisfied
with probability 3/4, so the identity fails at once).

**Not checked:** nothing else is involved (Markov's inequality).

Reproduce: `./run_all.sh` (Python 3.9+, standard library only; under a minute). Output: `results/c52a.txt`.
