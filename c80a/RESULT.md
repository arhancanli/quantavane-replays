# 80a, the Ising perceptron: alpha_star of [H2024] and [DS2025] from its definition, and [AT2024] 0.847 (15 November 2026)

Values recomputed from their definitions, not certificates. Read as text: Huang, arXiv 2404.18902 (e-print
`7bf2d750…9bf2`), §3; Ding and Sun, arXiv 1809.07742v1 (`aa33edbe…7267`), eqs. (q.psi.recursion), (GG) and
Proposition 1.3; Altschuler and Tikhomirov, arXiv 2401.15092 (`c3a51d6a…5e34`), the Proposition and its proof. No
program of any source was read or run.

**The definition** (kappa = 0; $\varphi$, $\Phi$ the normal density and distribution function, $Z \sim N(0,1)$):
$P(\psi) = \mathbb E\,\mathrm{th}(\sqrt\psi Z)^2$, $R_\alpha(q) = \frac{\alpha}{1-q}\,\mathbb E\,m(cZ)^2$ with
$m = \varphi/\Phi$ and $c = \sqrt{q/(1-q)}$, and
$\mathscr G(\alpha, q, \psi) = -\frac{(1-q)\psi}{2} + \mathbb E \log(2\,\mathrm{ch}(\sqrt\psi Z)) + \alpha\,\mathbb E\log\Phi(cZ)$.
$q_\star(\alpha)$ is the fixed point of $q = P(R_\alpha(q))$ in $(q_{lb}, q_{ub}) = (0.56394907949, 0.56394908030)$,
$\psi_\star = R_\alpha(q_\star)$, and $\alpha_\star$ is the root of $\mathscr G_\star(\alpha) = \mathscr G(\alpha, q_\star, \psi_\star)$
on $(0.833078599, 0.833078600)$, where Ding and Sun's Proposition 1.3 states that $\mathscr G_\star$ is well defined and
strictly decreasing (Huang quotes it in his §3). Both rows of the page, [H2024] (upper) and [DS2025] (lower),
are $\alpha_\star \approx 0.833$, each conditional on a numerical maximisation hypothesis that is not checked here.

**Claim checked:** $\mathscr G_\star(0.833078599) > 0 > \mathscr G_\star(0.833078600)$ with $q_\star$ inside
$(q_{lb}, q_{ub})$; with Proposition 1.3's monotonicity, $\alpha_\star \in (0.833078599, 0.833078600)$, which is
Proposition 1.3's own interval, and the page's $\approx 0.833$ holds.

- `c80a_a.py`: decimal intervals (60 digits; exp, ln and sqrt of the decimal module, correctly rounded, widened by
  one unit in the last place; $\Phi$ by its everywhere convergent series at 110 digits). A floating-point iteration
  gives $\tilde q$; then $P(R(\tilde q \mp 10^{-13})) - (\tilde q \mp 10^{-13})$ is shown to be positive, resp.
  negative, so $q_\star \in Q = [\tilde q - 10^{-13}, \tilde q + 10^{-13}]$; $\mathscr G$ is then enclosed with $Q$ and
  $R_\alpha(Q)$ as interval parameters. Integrals by Taylor models: on cells of radius $1/32$ over $[-10, 10]$, the
  integrand's Taylor polynomial of order 14 at the cell centre is integrated exactly, plus the Lagrange remainder
  from the 15th coefficient enclosed over the cell; the series come from differential equations
  ($\mathrm{th}' = 1 - \mathrm{th}^2$, $m' = -m(t+m)$, $(\log\Phi)' = m$, $\varphi' = -z\varphi$); tails beyond
  $|z| = 10$ by explicit bounds (Gordon's inequality for the Mills ratio).
- `c80a_b.py`: integer fixed point at $2^{-200}$ with its own exp (halving and squaring), log (atanh series), sqrt and
  $\pi$. The fixed point in $\psi$ instead: $R(P(\tilde\psi \mp 10^{-12})) - (\tilde\psi \mp 10^{-12})$ positive, resp.
  negative, and $q_\star = P(\psi_\star)$. $\Psi(y) = 1 - \Phi(y)$ by the alternating erf series for $y \le 4$ and by
  Laplace's continued fraction $\Psi(y) = \varphi(y)/(y + 1/(y + 2/(y + \dots)))$ beyond (successive approximants
  bracket it), so that $m(t) = 1/K(-t)$ and $\log\Phi(t) = \log\varphi(t) + \log K(-t)$ on the far left with no
  subtraction. The integrals of $m^2$ and $\log\Phi$ in the variable $t = cz$ on $[-11, 11]$ with the weight
  $\varphi(t/c)/c$; cells of width $1/32$ expanded at their LEFT end to order 13.

| | program A | program B |
|---|---|---|
| $\alpha = 0.833078599$: $q_\star$ in | $[0.5639490794961034, 0.5639490794963034]$ | $[0.5639490794961342, 0.5639490794962798]$ |
| $\mathscr G_\star(0.833078599)$ in | $[8.0290, 8.0461] \cdot 10^{-10}$ | $[8.0281, 8.0471] \cdot 10^{-10}$ |
| $\alpha = 0.833078600$: $q_\star$ in | $[0.5639490802916135, 0.5639490802918135]$ | $[0.5639490802916501, 0.5639490802917958]$ |
| $\mathscr G_\star(0.833078600)$ in | $[-2.8387, -2.8217] \cdot 10^{-10}$ | $[-2.8397, -2.8207] \cdot 10^{-10}$ |

Information: both programs also give $\mathscr G_\star(0.8330785997) > 0 > \mathscr G_\star(0.8330785998)$, so
$\alpha_\star \in (0.8330785997, 0.8330785998)$; $\psi_\star \approx 2.5763513$.

**[AT2024] 0.847.** The capacity is at most $\alpha$ when the Gardner-Derrida value $GD(\alpha) = \min_q GD(\alpha, q)$
is below $-\log 2$ (the source's Theorem, with Markov's inequality); at $q = \frac12$,
$GD(\alpha, \frac12) = \alpha\,\mathbb E\log H(u) + \frac12(1 + \log\frac12)$ with $\mathbb E\log H(u) = -1$ ($H(u)$ is
uniform). Both programs compute $\mathbb E\log H(u)$ by their quadrature ($-1$ to 18 digits) and
$GD(0.847, \frac12) + \log 2 = -4.264097200273 \cdot 10^{-4} < 0$, so the bound 0.847 holds. (Equivalently: the
step works for every $\alpha > \frac12 + \frac12\log 2 = 0.846573590280$.)

**Red controls**, each refused by both programs: $\alpha_\star$ claimed in $(0.833078600, 0.833078601)$ and in
$(0.833078598, 0.833078599)$ (the sign of $\mathscr G_\star$ is wrong at one end, and $q_\star$ leaves
$(q_{lb}, q_{ub})$ at the other); $F_q$ without its factor $(1-q)^{-1/2}$ ($q_\star \approx 0.3504$, outside the
interval); $\mathscr G$ with the sign of $\alpha\,\mathbb E\log\Phi(cZ)$ reversed ($\mathscr G_\star \approx 1.81$ at both
ends: no root here); [AT2024]'s step at $\alpha = 0.8465$ ($GD + \log 2 = +7.36 \cdot 10^{-5}$).

**Not checked:** the conditions of Proposition 1.3 other than the values (the bound on $(P \circ R_\alpha)'$, the
uniqueness of $q_\star$, the monotonicity of $\mathscr G_\star$); the maximisation hypotheses of [H2024] and
[DS2025]; the spherical-perceptron theorem used by [AT2024].

Reproduce: `./run_all.sh` (Python 3.9+, standard library only; about five minutes). Output: `results/c80a.txt`.
