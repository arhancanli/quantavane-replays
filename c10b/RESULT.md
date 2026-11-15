# 10b, the complex Grothendieck constant: Haagerup's upper bound and the lower bound of [GFL26] (15 November 2026)

Two replays on the 10b page, each by two separately written programs that agree on every number they both print
(`compare.py` checks that every pair of enclosures overlaps). Nothing is downloaded: every input is a definition or an
exact decimal taken from the papers below, read as text. No program of any source was read or run.

Sources (arXiv source files, SHA-256 of the TeX file read):

- [GFL26] S. Guo, E. X. Fang, J. Lu, *An Improved Lower Bound for the Complex Grothendieck Constant*, arXiv 2609.07000v1
  (the only version), `main.tex` `7fc49256…c8ed`: Sections 2 to 4 and Appendix A.1.
- [HJN26] S. Heilman, C. Jones, G. Malavolta, *Sharper Bounds for the Complex Grothendieck Constant*, arXiv 2609.25462v1,
  `mainv11.tex` `b10e9e91…a1fd`: Section 2 (Haagerup's equation (2.5) and $\gamma_H$ (2.6)).
- G. Pisier, *Grothendieck's theorem, past and present*, arXiv 1101.4195, `grothendieck.arxiv-march31.tex`
  `214bc299…3e48`: Section 4 (Haagerup's bound and conjecture). It states the same equation as [HJN26].

## 1. [H1987] $C_{10b} \le 1.40491$ (Haagerup), from its definition

$K_0 \in (0, 1)$ is the unique solution of $K \int_0^{\pi/2} \cos^2 t\,(1 - K^2 \sin^2 t)^{-1/2}\,dt = \frac{\pi}{8}(K + 1)$,
and the bound is $8/(\pi(K_0 + 1))$.

- `c10b_series.py`: exact rationals. The integral is $\frac{\pi}{4} \sum_j a_j K^{2j}$, $a_j = \binom{2j}{j}^2/(16^j(j+1))$
  (Wallis; the coefficients checked three ways), so $\pi$ cancels from the equation; $K_0$ bracketed by exact
  evaluation with a rigorous tail; $\pi$ by Machin's formula with alternating-series bounds.
- `c10b_agm.py`: decimal intervals, every operation rounded outward. $K(m)$ and $E(m)$ by the arithmetic-geometric mean
  (Gauss, Legendre) with a bounded tail; $\pi$ by Salamin-Brent; $K_0$ by interval bisection. Legendre's relation at
  $m = 1/2$ holds to $5 \cdot 10^{-42}$.

| quantity | both programs |
|---|---|
| $K_0$ | $0.812557858821472442318593494612\ldots$ |
| $8/(\pi(K_0+1)) = \gamma_H^{-1}$ | $1.404909132735795535525448150614\ldots$: the printed 1.40491 holds ($8.7 \cdot 10^{-7}$ above) |
| information: Haagerup's conjectured value $(\int_0^{\pi/2} \cos^2 t\,(1 + \sin^2 t)^{-1/2}\,dt)^{-1}$ | $1.404575934663742032773958471548\ldots$ |
| information: [P1978] $e^{1-\gamma}$ | $1.5262051115958638804748887150\ldots$: the printed 1.52621 holds as an upper bound |

Also, for the rows PR #196 adds: $\gamma_H = \pi(K_0+1)/8 = 0.7117898066849979291\ldots$, below the
$0.71178980668499814812$ that [HJN26] states, and $\gamma_H^{-1} - 10^{-5} = 1.40489913273579553\ldots$ exceeds
[HJN26]'s $1.404898554746$ by $5.78 \cdot 10^{-7}$, as its Theorem 1.1 displays. Not checked: Haagerup's theorem
itself, and [HJN26]'s own bounds 1.404898554746 and 1.4 (their certificates are supplied with their code).

## 2. [GFL26] $C_{10b} > 1.35584631827168$ (PR #196): the certificate of Proposition 4.1

The finite conditions of the paper's Theorem 3.1 for the exact parameters of its Section 4.1 ($\Lambda = \{1, \ldots,
10\}$, $\rho$, $\beta_1, \ldots, \beta_{10}$, and the weight $W(s) = \alpha + \eta s + \zeta_1 s/(\delta_1 + s) +
\zeta_2 s/(\delta_2 + s)$; 17 decimals, each read as an exact rational): every Gram block $I - \mathcal M^{p,q}(W)$
positive definite, $S = R - b^T(B^{-1} + G)^{-1} b \le 1$, the weight mean $\mu_W$ and the ratio
$(1-\rho)/(\rho+\mu_W)$.

- `c10b_gfl_a.py`: balls on exact rationals (centre and radius, rounded to $2^{-256}$). Radial factors from the
  definition of $h_{a,b}$; $1/W = \Delta/N$ with $N$ cubic, its roots by exact rational bisection, residues, and each
  integral $\int P e^{-s}/(s + \lambda)\,ds$ by synthetic division plus $P(-\lambda) e^\lambda E_1(\lambda)$, $E_1$ by its
  power series with Euler's $\gamma$ (Euler-Maclaurin with Bernoulli numbers) and logarithms by atanh series; LDL in
  increasing index order; $S$ by elimination.
- `c10b_gfl_b.py`: decimal intervals at 100 digits. Radial factors in the Laguerre form of the paper's appendix;
  moments $\int s^k e^{-s}/N(s)\,ds$ for $k \le 2$ from the roots and for $k \ge 3$ by the recurrence over the cubic
  ($\sum_i n_i U_{k+i} = k!$); $e^\lambda E_1(\lambda)$ by its continued fraction (A&S 5.1.22), enclosed by two successive
  convergents with an alternation check; $S$ as the last pivot of the bordered LDL of
  $\left[\begin{smallmatrix} B^{-1}+G & b \\ b^T & R \end{smallmatrix}\right]$.

| statement of the paper | both programs |
|---|---|
| the poles of $1/W$ and $1/(2\rho + W)$ are real, simple and negative | $-\lambda$ with $\lambda = 11.464962864239\ldots$, $1.034936591516\ldots$, $0.177072838159\ldots$, and $11.503086144878\ldots$, $1.055532023895\ldots$, $0.198743266566\ldots$ |
| 132 nonempty blocks, 570 pivots, all positive; least pivot $> 1.2259567912 \cdot 10^{-5}$, in block $(2,1)$ | 132 blocks, 570 pivots, all positive; least pivot $0.000012259567912\ldots$, block $(p,q) = (2,1)$ |
| $1 - S > 1.8374949151 \cdot 10^{-9}$ | $S = 0.9999999981625050848743609\ldots$, $1 - S = 1.83749491512\ldots \cdot 10^{-9}$ |
| $\mu_W < 0.716627735168304$ | $\mu_W = 0.7166277351683039056785778\ldots$ |
| $(1-\rho)/(\rho+\mu_W) > 1.35584631827168$ | $1.3558463182716895859699624\ldots$ |

So the certificate is accepted, and with Theorem 3.1 it gives $K_G^{\mathbb C} > 1.35584631827168$.

**Not checked:** Theorem 3.1 (the weighted criterion) and the reductions of Sections 2 and 3 (read), and the dual
bound $K_* < 1.35584697425050$ of Section 5 (not replayed).

## Red controls

Each refused by both programs of its pair:

- H1: [H1987] claimed as 1.4049091327357955 (the value cut after 16 decimals); H2: claimed as 1.40490; H3: [P1978]
  claimed as 1.52620.
- G1: $\rho = 0.012$ (lowered by $3.9 \cdot 10^{-5}$): the Gram blocks still hold, but $1 - S = -9.28 \cdot 10^{-5}$.
- G2: $\beta_{10} = 0.0148$ (raised by 0.7%): block $(2,1)$ fails (pivot $-0.00702$), while $S < 1$.
- G3: $\zeta_2 = 1.3819$ (the weight lowered, which would raise the ratio to 1.35586512): both conditions fail.
- G4: the ratio claimed above 1.35584631827169; G5: $\mu_W$ claimed below 0.716627735168303.

## Notes on method

- The first run of the two programs disagreed in the 23rd digit of $S$ while each enclosure was $10^{-25}$ wide. The
  cause was in program B: Python's unary minus on a `Decimal` rounds to the default 28-digit context, so the poles
  passed to $E_1$ were off by about $10^{-28}$. It now negates exactly (`copy_negate`), and the same fix was made in
  `c10b_agm.py` (also `copy_abs` for the AGM tail bound). The two programs agree since; the exact $S$ recomputed in
  rational arithmetic from program A's centres confirms the value above.

Reproduce: `./run_all.sh` (Python 3.9+, standard library only; about three minutes on two cores, under 200 MB). Output:
`results/c10b.txt`.
