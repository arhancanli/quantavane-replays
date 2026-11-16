# 19a [E1956]: Esseen's lower bound $c_E = \frac{\sqrt{10}+3}{6\sqrt{2\pi}} \approx 0.4097321837$ from its definition (16 November 2026)

A value from its definition, as the row states it ("achieved asymptotically by (centered/normalized) Bernoulli
sums"). For the two-point law $P(X=1) = p$, standardized, Esseen's lattice expansion gives
$\sqrt n\,\sup_x\lvert F_n(x)-\Phi(x)\rvert \to (\lvert\mu_3\rvert/6 + h/2)\,\varphi(0)$ with span $h = 1/\sqrt{pq}$ and
$\mu_3 = (q-p)/\sqrt{pq}$; divided by $\beta_3 = (p^2+q^2)/\sqrt{pq}$ this is
$g(p) = \frac{\lvert 1-2p\rvert/6 + 1/2}{(p^2+q^2)\sqrt{2\pi}}$, and $c_E = \max_p g(p)$. The expansion itself is
Esseen's theorem (read, not re-proved).

- `c19a_a.py`: with $u = \lvert1-2p\rvert$, $g = \frac{u+3}{3(1+u^2)\sqrt{2\pi}}$; the derivative's numerator
  $1 - 6u - u^2$ vanishes once on $[0,1]$, at $u^* = \sqrt{10}-3$ ($p^* = (4-\sqrt{10})/2 = 0.41886116991581\ldots$),
  and $f(u^*) = (\sqrt{10}+3)/2$ exactly in $\mathbb Q(\sqrt{10})$ (larger than $f(0) = 3$ and $f(1) = 2$); decimal
  value in integer fixed point with Machin's $\pi$.
- `c19a_b.py`: no calculus: branch-and-bound over $p \in (0, \frac12]$ with mean-value bounds (derivative enclosed by
  monotonicity), exact Fractions; the maximum of $\sqrt{2\pi}\,g$ is enclosed in an interval of width $6\cdot10^{-31}$
  that contains $(\sqrt{10}+3)/6$; $\pi$ by Bailey-Borwein-Plouffe.

Both: $c_E = 0.409732183702396342990568823497\ldots$, so the page's $0.4097321837$ is the truncation and the correct
rounding (and below the value, the safe side for a lower bound).

**Red controls**, refused by both (and giving the same wrong values in both): the value printed as 0.4097321838; the
lattice term $h/2$ doubled ($0.80338747289769276035\ldots$); the factor $p^2+q^2$ of $\beta_3$ dropped ($0.26596152026762178529\ldots$).

Reproduce: `./run_all.sh` (Python 3.9+, standard library only; under a second). Output: `results/c19a.txt`.
