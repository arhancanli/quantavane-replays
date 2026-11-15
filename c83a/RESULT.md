# 83a, the Wirsing constant: [Schweiger1980] 0.38197 and [MR1987] 0.30366299, 0.30366327 from the definition (15 November 2026)

A value recomputed from its definition, not a replay of any source's computation. $C_{83} = |\lambda_2|$, the
subdominant eigenvalue of the Gauss-Kuzmin-Wirsing operator $\mathcal L f(x) = \sum_{k\ge1}(x+k)^{-2} f(1/(x+k))$
(the page's definition). The page's rows checked here: [Schweiger1980] $C_{83} \le 2/(3+\sqrt5) \approx 0.38197$ and
[MR1987] $0.30366299 \le C_{83} \le 0.30366327$ (eq. (5.20)).

**Reduction (Wirsing 1974; Knuth, TAOCP vol. 2, §4.5.3).** With $\varphi = (1+x) f$, $\mathcal L$ becomes
$T\varphi(x) = \sum_k p_k(x)\, \varphi(1/(x+k))$, $p_k = \frac{1+x}{(x+k)(x+k+1)}$, $\sum_k p_k = 1$. Abel summation with
$p_1 + \dots + p_k = k/(x+k+1)$ gives $(T\varphi)' = -V(\varphi')$ for the positive operator
$$V\psi(x) = \sum_{k\ge1}\Big[\frac{k}{(x+k+1)^2}\int_{1/(x+k+1)}^{1/(x+k)}\psi + \frac{1+x}{(x+k)^3(x+k+1)}\,\psi\Big(\frac1{x+k}\Big)\Big].$$
Its spectral radius is $|\lambda_2|$ (Wirsing's theorem; not re-proved here). If $\psi > 0$ on $[0,1]$ and
$a\psi \le V\psi \le b\psi$ there, then $a^n\psi \le V^n\psi \le b^n\psi$, so $a \le \rho(V) \le b$ (Collatz-Wielandt).
The identity $(T\varphi)' = -V(\varphi')$ is checked numerically in program B's self-test (the Gauss-Kuzmin operator
summed directly to $k = 400{,}000$ and differentiated: $0.6188707900$ against $0.6188707878$ at $x = 0.3$).

**The certificate** `psi.json`: a polynomial $\psi$ of degree 18 (coefficients multiples of $2^{-64}$), found by
`c83a_search.py` (collocation at 19 Chebyshev points and power iteration, 60 digits, no error control; it reproduces
`psi.json` byte for byte) and trusted for nothing. Claimed: $0.303663002898732\,\psi \le V\psi \le
0.303663002898733\,\psi$ on $[0, 1]$.

- `c83a_a.py`: $V\psi$ in closed form. With $y = x + k$, $X = 1 + x$ and $k = (y+1) - X$, every summand is a sum of
  $(\alpha + \beta X)\, y^{-a}(y+1)^{-b}$; partial fractions in $y$ and $y + 1$ and the sum over $k$ give
  $V\psi = \sum_{i\ge2} R_i(X)\,\zeta(i, X) + E(X)$ ($R_i$ linear, $E$ a Laurent polynomial; the $1/y$ parts cancel
  exactly, as checked). Hurwitz $\zeta(n, X_0)$ by Euler-Maclaurin (30 terms, 25 Bernoulli corrections, remainder
  $4\zeta(52)(n)_{51}/((2\pi)^{52}(X_0+30)^{n+51})$). $[0,1]$ in 32 cells; on each the Taylor polynomial of
  $V\psi - c\psi$ to order 25 with Lagrange remainders ($\zeta(i,\cdot)$ and $X^{-i}$ are completely monotone).
  Exact rationals rounded outward to $2^{-300}$. Self-tests: the partial fractions, $\zeta(2) = \pi^2/6$ and
  $\zeta(3)$ to 44 digits, $\zeta(5,a) - \zeta(5,a+1) = a^{-5}$, the closed form against direct summation.
- `c83a_b.py`: no partial fractions and no Hurwitz zeta at non-integer points. With $w = 1/(x + k + \frac12)$:
  $1/(x+k) = 2w/(2-w)$, $1/(x+k+1) = 2w/(2+w)$, $k = 1/w - \frac12 - x$; the summand is a power series in $w$ with
  coefficients linear in $x$, built by series arithmetic (Horner for $\psi(1/(x+k))$) to $w^{150}$, tail by a Cauchy
  bound on $|w| = 6/5$ (the summand is analytic in $|w| < 2$, and $|w| \le 2/3$). The sum over $k$ of $w^s$ is
  $\zeta(s, x + \frac32)$; about $x = \frac12$ this is one Taylor polynomial in $t = x - \frac12$ of degree 201
  needing only $\zeta(n, 2) = \zeta(n) - 1$, by the Cohen-Rodriguez Villegas-Zagier acceleration of
  $\eta(n) = (1 - 2^{1-n})\zeta(n)$ (error at most $\eta(n)/T_{100}(3)$) for $n \le 40$ and by direct summation above.
  Decimal intervals at 70 digits; the polynomial shifted to the centres of 32 cells. Self-tests: $\zeta(2)$,
  $\zeta(4)$ against $\pi$ by Machin, $\eta(1) = \log 2$, and the derivative identity above.

| statement | program A | program B |
|---|---|---|
| $\psi > 0$ on $[0,1]$ | least cell bound $0.99972$ | $0.99972$ |
| $V\psi - 0.303663002898732\,\psi \ge 0$ | least cell bound $5.27 \cdot 10^{-16}$ | $5.27 \cdot 10^{-16}$ |
| $0.303663002898733\,\psi - V\psi \ge 0$ | least cell bound $1.93 \cdot 10^{-16}$ | $1.93 \cdot 10^{-16}$ |
| so $0.303663002898732 \le C_{83} \le 0.303663002898733$ | yes | yes |
| [Schweiger1980] $C_{83} \le 2/(3+\sqrt5)$; [MR1987] $C_{83} \le 0.30366327$ and $\ge 0.30366299$ | hold | hold |

Information: the first 50 digits of the 175-digit centre value of [GKW2026] lie in our interval (15 digits agree);
the 175-digit row itself is not replayed.

**Red controls**, each refused by both programs (with the cells where the inequality fails, not merely unproved):
the lower bound claimed as $0.303663002898733$; the upper bound claimed as $0.303663002898732$; coefficient 9 of
$\psi$ changed by $2^{-30}$; the point term of $V$ with $(x+k)^2$ in place of $(x+k)^3$.

**Not checked:** Wirsing's theorem that $\rho(V) = |\lambda_2|$ (used as stated by Wirsing and by Knuth).

Reproduce: `./run_all.sh` (Python 3.9+, standard library only; about five minutes). Output: `results/c83a.txt`.
