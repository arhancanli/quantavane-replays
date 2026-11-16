# 7a: Zeilberger and Zudilin's bound [ZZ2020] $\mu(\pi) \le 7.103205334137\ldots$ (16 November 2026)

From the paper's final formula (arXiv 1912.06345, Proposition 2, eq. (18), Lemma 3, eq. (10), section "World
record"): $N_1, N_2, N_3$ the roots of $108N^3 - 2359989N^2 + 138304N - 2048$, $c = \pi/(2\sqrt3) - \log(3\sqrt3/4)$,
$\mu \le 1 + (\log N_3 - \frac52\log 2 + 4 - c)/(-(\log\lvert N_1\rvert - \frac52\log2 + 4 - c))$. The arithmetic lemmas
and the saddle-point analysis are read, not checked. Two separately written programs:

- `c7azz_a.py`: exact rational intervals; the real root $N_3$ by bisection (the discriminant is negative), and
  $\lvert N_1\rvert^2 = 2048/(108 N_3)$ from the product of the roots; $c = \psi(2/3) - \psi(1/2)$ checked from Gauss's values.
- `c7azz_b.py`: 60-digit decimals; the saddle points $y_j$ of $g(y) = y(y^2+6y+25)^2/(y-25)^3$ (the stationary
  polynomial $2y^3 - 125y^2 - 500y - 625$, expanded and checked at random rationals), $N_j = g(y_j)$ checked against the cubic.

| quantity | value |
|---|---|
| $N_3$, $\lvert N_1\rvert$ | $21851.6913962195616968\ldots$, $0.029458495928116876\ldots$ |
| growth rates | $-1.90291648559998216823\ldots$ and $11.61389004533104820528\ldots$ |
| $1 + \sigma/\tau$ | $7.1032053341370017275057734\ldots$ |

Note: the row prints 7.103205334137, the truncation; the bound proved is $1.7\cdot10^{-15}$ above it (as an upper
bound, 7.103205334138 would be the safe rounding; the paper writes $7.103205334137\ldots$). Red controls, refused by
each: the claim 7.103205334136; 2047 for 2048 in the cubic and $2\log2$ for $\frac52\log 2$ (A); $(y-25)^2$ in $g$ and
$2\log 2$ (B).

Reproduce: `./run_all.sh` (Python 3.9+; seconds). Output: `results/c7azz.txt`.
