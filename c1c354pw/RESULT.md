# 1c: the piecewise-polynomial certificate of PR #354, $C_{1c} \ge 0.410378195$\* (14 November 2026)

Pull request: [#354](https://github.com/teorth/optimizationproblems/pull/354) (SeanSagait), update of 13 November 2026
(commit `a6fae1a0` of the pull request). Page: [1c](https://teorth.github.io/optimizationproblems/constants/1c.html).

## Input

One data file, fetched read-only: SeanSagait/c1c-autocorrelation-certificate @ `308005a40a6df640b4f2124ef7d812522b792f42`,
`pwcert_p20_dps40_S120.json` (55,076 bytes), SHA-256
`8aa13b614c8a4751dca0a2dde7ad556b09d9b9f225e824b0774ee22d7518340f`, the hash stated in the pull request and in the
repository's `SHA256SUMS`. The repository's README was read for the format. Its programs (`verify_pwcert.py`,
`verify_pwcert_independent.py`, `construction/`) were neither read nor run.

## The statement

$f = s^{-1} 1_{[-s/2,\,s/2)} + h$ with $s = 2^{-120}$; $h$ is zero below $t_1 = 22545746012851/2^{47}$, a polynomial of
degree 20 with rational coefficients on each of nine intervals from $t_1$ to $1 - t_1$, the constant
$T = 13866842223167558038152245155365550010625/2^{132}$ on $[1 - t_1, 1 + s/2)$, and zero beyond. The certificate
states $\int f = 1 + M$ and $g(t) = \int f(x) f(x+t)\,dx \ge T - 10^{-60}$ for $0 \le t \le 1$, hence
$C_{1c} \ge (T - 10^{-60})/(1 + M)^2$, the rational `bound` of the file.

## What was checked, and how

Two separately written programs, both exact (`fractions.Fraction`), standard library only:

| step | program A (`c1c354pw_a.py`) | program B (`c1c354pw_b.py`) |
| --- | --- | --- |
| structure | pieces contiguous from $t_1$ to $1 + s/2$, spike left of $t_1$, last two pieces equal to $T$ | the same, coded separately |
| $f \ge 0$ | Bernstein coefficients of every piece | centred Taylor forms with bisection, after dividing out zeros at the left end of a piece |
| $\int f$ | $= 1 + M$ exactly | $= 1 + M$ exactly |
| $g$ on $[0, 1 + s]$ | for each ordered pair of pieces, the bivariate antiderivative $\int_0^y P_i(u) P_j(u + \tau)\,du$ expanded in $y$ and $\tau$, the limits substituted on each range of $\tau$ | the jump formula $g(t) = \sum (-1)^{r+1} D_{b,r} D_{b',r'} \frac{r!\,r'!}{(r+r'+1)!} (t - (b'-b))_+^{r+r'+1}$ over the breakpoints $b, b'$ and the jumps $D$ of the Taylor coefficients; no product of pieces is integrated |
| consistency of $g$ | continuous at all 40 interior mesh points; $g(0) = \int f^2$; $2\int_0^\infty g = (\int f)^2$; $g = 0$ at $1 + s$ | the closed form vanishes identically beyond the last lag; it equals a direct integral of $f(x) f(x+t)$ at 11 rational lags |
| $g \ge T - 10^{-60}$ on $[0, 1]$ | Bernstein coefficients on every cell (48 leaves, depth at most 4) | centred Taylor forms (128 intervals, depth at most 35) |
| the bound | `bound` $= (T - 10^{-60})/(1 + M)^2$ exactly, and $\ge 0.410378195$ | the same by cross-multiplication |

## Result

- Both programs accept the certificate.
- Both obtain **the same exact $g$**: 41 polynomial cells on $[0, 1 + s]$ (42 mesh points; degree at most 41), with
  identical coefficients. The SHA-256 of the canonical listing (cell ends and coefficients in powers of $t - u$) is
  `4ef298f0bc3bd6c7495ff9a916b57ca646a0bb4032e4e13e904c39a97c9e911d` from each.
- $\int f = 1 + M = 2.491247835698739976\ldots$
- $\min_{[0,1]} g \ge T - 10^{-60}$ holds. On the two cells covering $[1 - t_1, 1]$, $g = T$ identically, as the
  repository's README says.
- `bound` $= 0.410378195722515456620959705443384843306244844\ldots$ (a rational with a 728-bit denominator), equal to
  $(T - 10^{-60})/(1 + M)^2$; the printed $0.410378195$ is a valid rounding down, and the digits agree with the
  39 digits printed in the pull request.

Information (program A only): the same Bernstein test with $T$ in place of $T - 10^{-60}$ also passes, so
$g \ge T$ on $[0, 1]$ with equality on $[1 - t_1, 1]$; the margin is not used up. Program B cannot decide that
sharper statement, because a centred form cannot certify a bound that is attained.

## Red controls

All refused by both programs (`results/pr354pw.txt`):

1. The bound on $g$ raised to $T + 2^{-132}$: the positivity test fails in both.
2. The constant coefficient of piece 4 lowered by $10^{-20}$: $\int f \ne 1 + M$, the positivity test fails, and the
   digest of $g$ changes (to the same new value in both programs).
3. $h$ on $[1 - t_1, 1)$ lowered from $T$ by $10^{-30}$: the structure check, $\int f$ and the positivity test fail.

## What was NOT checked

- How the function was found (the collocation solve, the choice of $t_1$): not needed for validity, not looked at.
- The pull request's remark that the value matches the continuum limit of this family to 31 digits.
- The earlier step-function certificate with $K = 2^{26}$ cells ($0.41037819264\ldots$, kept in the same repository):
  not replayed here; 67 million weights of up to 58 bits need more memory than this machine's limit allows.
- Nothing about upper bounds for $C_{1c}$.

What rests on mathematics rather than computation: that a non-negative integrable $f$ with
$\min_{[0,1]} g \ge m$ gives $C_{1c} \ge m/(\int f)^2$ (the definition on the 1c page), and, for program B, the jump
formula above (tested in the program as described).

## Reproduce

`./fetch.sh && ./run_all.sh` (Python 3.9+, standard library only; about four minutes on two cores, each run about
45 seconds; under 200 MB).
