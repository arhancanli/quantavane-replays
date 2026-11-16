# 8a: the lower bound $C_8 \ge 2/\log\gamma_1 \approx 0.755106$, with a certified zero of $\zeta$ (16 November 2026)

The row: a zero $\frac12 + i\gamma$ of $\zeta$ with $\gamma \ge 2$ forces $\frac12 \le 1 - \frac{1}{R\log\gamma}$, so
$C_8 = R \ge 2/\log\gamma$; the page takes $\gamma_1 = 14.1347\ldots$ and prints $2/\log\gamma_1 \approx 0.755106$.
Both programs prove, by a sign change of Hardy's real function $Z(t) = e^{i\theta(t)}\zeta(\frac12+it)$, that $\zeta$
has a zero on the critical line with

$$14.13472514173469378 < \gamma < 14.13472514173469380,$$

so $R \ge 2/\log\gamma \ge 2/\log(14.13472514173469380) = 0.75510606143218443826\ldots$, and $2/\log\gamma$ lies in
$[0.7551060614321844382602636, 0.7551060614321844386636567]$: the page's $0.755106$ is both the truncation and the
correct rounding, and it is below the value (the safe side for a lower bound). That this zero is the lowest one is
classical and is not checked here; it is not needed for the bound.

| | program A (`c8a_a.py`) | program B (`c8a_b.py`) |
|---|---|---|
| $\zeta(\frac12+it)$ | Euler-Maclaurin, $N = 40$, 30 Bernoulli terms, remainder $\frac{\lvert B_{2p}\rvert}{(2p)!}\lvert s\cdots(s+2p-1)\rvert\frac{N^{1-\sigma-2p}}{\sigma+2p-1}$ | $\eta(s)/(1-2^{1-s})$, $\eta$ by Algorithm 1 of Cohen, Rodriguez Villegas and Zagier ($n = 90$, exact rational weights), error $\le 2\sqrt{\cosh \pi t}/(3+\sqrt8)^n$ (derived in the file from $\int\lvert d\mu\rvert = \Gamma(\sigma)/\lvert\Gamma(s)\rvert$) |
| $e^{i\theta(t)}$ | Stirling for $\operatorname{Im}\log\Gamma(\frac14+\frac{it}2)$ at shift 20, 20 terms | the same series at shift 30, 15 terms; Bernoulli numbers by Akiyama-Tanigawa |
| arithmetic | integer intervals at $2^{-320}$, own sqrt, ln, atan, cos, sin, Machin's $\pi$ | decimal intervals at 110 digits, own cos, sin, atan, $\pi$ by Bailey-Borwein-Plouffe |
| $Z(10)$ | $-1.549194546181022\ldots$ | $-1.549194546181022\ldots$ |
| $Z(14.13472514173469378)$ | $-8.294279\cdot10^{-18}$ (width $2\cdot10^{-55}$) | $-8.294279\cdot10^{-18}$ |
| $Z(14.13472514173469380)$ | $+7.568930\cdot10^{-18}$ | $+7.568930\cdot10^{-18}$ |

Self-test in both: $\operatorname{Im}(e^{i\theta}\zeta)$ encloses 0 at every point evaluated ($Z$ is real only if
$\theta$ and $\zeta$ are both right). B also checks its weights on $\sum(-1)^k/(k+1) = \log 2$.

**Red controls**, refused by both: $\theta$ with the wrong sign (the imaginary part no longer encloses 0); a bracket
$(14.13472514173469381, 14.13472514173469383)$ above $\gamma_1$ (no sign change).

Reproduce: `./run_all.sh` (Python 3.9+, standard library only; a few seconds). Output: `results/c8a.txt`.
