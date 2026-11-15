# 57a and 57b, Bloch's and Landau's constants: the closed forms of [AG1937] and [Rad1943] (15 November 2026)

Closed forms named in the rows, recomputed and checked against the page's decimals; the extremal-function analysis
behind them (Ahlfors-Grunsky, Rademacher) is not checked. No program of any source was read or run.

- 57a [AG1937]: $B \le \frac{1}{\sqrt{1+\sqrt3}}\frac{\Gamma(1/3)\Gamma(11/12)}{\Gamma(1/4)} \approx 0.4719$.
- 57b [Rad1943]: $L \le \frac{\Gamma(1/3)\Gamma(5/6)}{\Gamma(1/6)} \approx 0.5433$.

`c57ab_a.py`: $\Gamma(x) = \exp(\log\Gamma(x+40))/\prod_{j<40}(x+j)$ with Stirling's series (20 Bernoulli terms; for real
$z > 0$ the remainder is at most the first omitted term), integer fixed point at $2^{-200}$ with its own exp, log, sqrt
and $\pi$; the final division by the exact rational product. `c57ab_b.py`: for $0 < s \le 1$ and $N = 120$,
$\Gamma(s) = N^s e^{-N}\sum_{k\ge0} N^k/(s(s+1)\cdots(s+k)) + E$ with $0 \le E \le N^{s-1}e^{-N}$, decimal intervals at
110 digits. Self-tests: $\Gamma(1/2)^2 = \pi$, $\Gamma(1) = 1$, $\Gamma(1/6) = 2^{-1/3}(3/\pi)^{1/2}\Gamma(1/3)^2$ (A);
$\Gamma(1/12)\Gamma(11/12) = \pi/\sin(\pi/12)$ (B).

| value | both programs (first 20 digits identical, checked by `run_all.sh`) |
|---|---|
| $\Gamma(1/3)$, $\Gamma(1/4)$ | $2.67893853470774763365\ldots$, $3.62560990822190831193\ldots$ |
| $\Gamma(11/12)$, $\Gamma(5/6)$, $\Gamma(1/6)$ | $1.05554656481346630231\ldots$, $1.12878702990812596126\ldots$, $5.56631600178023520425\ldots$ |
| [AG1937] | $0.47186165345268178487\ldots \approx 0.4719$: holds |
| [Rad1943] | $0.54325896534297670695\ldots \approx 0.5433$: holds |

**Red controls**, each refused by both programs: the Ahlfors-Grunsky form without $1/\sqrt{1+\sqrt3}$; Rademacher's
form with $\Gamma(1/6)$ and $\Gamma(5/6)$ swapped; the Ahlfors-Grunsky value printed as 0.4718.

**Not checked:** that these are the Bloch radius and the Landau radius of the extremal functions; the lower bounds
[CG1996] and [Yan1995] (theorems).

Reproduce: `./run_all.sh` (Python 3.9+, standard library only; under a minute). Output: `results/c57ab.txt`.
