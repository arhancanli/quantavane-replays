# 8a: the finite inputs of PR #373's deduction $C_8 \le 2/\log\gamma_1$ (16 November 2026)

PR [#373](https://github.com/teorth/optimizationproblems/pull/373) deduces $C_8 \le 2/\log\gamma_1 \approx 0.755107$\*
from OpenAI's $7/8$ quasi-Riemann hypothesis [OAI2026] and Platt and Trudgian's verification of RH to $3\cdot10^{12}$
[PT2021]. Two separately written programs check every finite input of that deduction:

| statement | result |
|---|---|
| (1) a zero on the critical line near $\gamma_1$ | Hardy's $Z$ changes sign on $(14.13472514173469378, 14.13472514173469380)$ (as in folder `c8a`) |
| (2) it is the lowest zero | $\operatorname{Re}\zeta(\sigma + 17.8i) > 0$ for $1/2 \le \sigma \le 2$ (A: at least 1.348; B: at least 0.947), so $\lvert S(17.8)\rvert < 1/2$; $\theta(17.8)/\pi + 1 = 0.99243438671870\ldots$, so $N(17.8) = 1$: the zero of (1) is the only zero with $0 < \gamma \le 17.8$ |
| (3) the arithmetic | $2/\log\gamma_1 = 0.75510606143218443826\ldots \le 0.755107$; $\gamma_1^4 < 39916.23 < 4\cdot 10^4$; $8/\log(3\cdot10^{12}) \le 0.2784582 < 0.279$ |

So any zero with $2 \le \lvert\gamma\rvert \le 3\cdot10^{12}$ (on the line by [PT2021]) needs $R \ge 2/\log\lvert\gamma\rvert \le
2/\log\gamma_1$, and any zero above $3\cdot10^{12}$ (with $\beta \le 7/8$ by [OAI2026]) needs $R \ge 8/\log\lvert\gamma\rvert
< 0.279$. The identity $N(T) = \theta(T)/\pi + 1 + S(T)$ is read, not checked. **Not checked:** [OAI2026]'s theorem
and [PT2021]'s computation.

- `c8a373_a.py` (A): Euler-Maclaurin ($N = 40$, 30 terms) with $\sigma$ as an interval, Stirling at shift 30 for
  $\theta$; integer intervals at $2^{-320}$, our own ln, exp, cos, sin, atan.
- `c8a373_b.py` (B): the Cohen-Rodriguez Villegas-Zagier acceleration of $\eta$ ($n = 90$) with $\sigma$ as an interval;
  its error bound $2\sqrt{\cosh \pi t}/(3+\sqrt8)^n$ holds for every $\sigma \ge 1/2$ (derived in the file); Stirling at
  shift 40; decimal intervals at 110 digits.

**Red controls**, refused by both: the same count at height 22 claimed as $N = 1$ (both find $N(22) = 2$); 0.755106
claimed as the upper bound; a bracket just above $\gamma_1$ (no sign change).

Reproduce: `./run_all.sh` (Python 3.9+; seconds). Output: `results/c8a373.txt`.
