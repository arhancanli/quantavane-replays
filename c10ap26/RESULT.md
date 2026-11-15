# 10a [P2026] $K_G^{\mathbb R} \le 1.7813$\*: a replay of the computational content (15 November 2026)

Source: Bo Peng, *A computer-assisted upper bound of 1.7813 for the real Grothendieck constant*, arXiv 2609.20074v1
(e-print `7db3732a…a36c`), read as text: Sections 2 to 5 and the exact parameters of Appendix A. No program of the
source (its ancillary files) was read or run. Our code runs on Arb through python-flint 0.6.0 from PyPI
(`REQUIREMENTS.txt` records the wheel and its SHA-256).

**What the theorem needs from computation.** With $P$, $\rho$ and $\gamma = 5000\pi/17813$ as in the paper,
$A_{ab} = \mathbb E[\operatorname{sgn}(W + P(X))\,h_a(W)h_b(X)]$ and $H(t) = \frac\pi2\sum A_{ab}^2\rho(t)^a(-t)^b = \sum b_n t^n$,
Proposition 2.3 gives $K_G \le \pi/(2\gamma) = 1.7813$ once $\gamma + \sum_{n\ge3}|b_n| < b_1$. The paper splits the sum at
$n = 301$: the head from the 22,952 coefficients $A_{ab}$ with $a + b \le 301$, and the tail by Proposition 4.2,
$\sum_{n>301}|b_n| \le B_3/\sqrt{10\cdot301^5}$ with $B_3^2 \le C_a I$, $C_a = \frac12 + \frac{\pi}{2\sqrt a}\coth\frac{\pi}{\sqrt a}$,
$a = 1/1000$, and $I$ the integral of the integrand displayed in Section 5.2 over $x, y \in \mathbb R$, $s = \theta/\pi \in [0,1]$.
Lemma 5.1 states $b_1 > 0.881850816$, $\sum_{3\le n\le301}|b_n| < 3\cdot10^{-8}$, $B_3 < 107$.

## The head

| | program A (`p26_a.py`) | program B (`p26_b.py`) | the paper |
|---|---|---|---|
| method | composite Gauss-Legendre on $[0, 11/2]$ (43 panels of 24 nodes, Arb's nodes and weights); error by Trefethen's Bernstein-ellipse bound with $|h_n(x+iy)| \le 1.0865\,e^{x^2/4}e^{|y|\sqrt n}$ (Cramér) and a disk enclosure of $P$; Hermite values by the three-term recurrence; tails on $[11/2,\infty)$ by Cramér's inequality and $\int_c^\infty h_b\varphi = h_{b-1}(c)\varphi(c)/\sqrt b$ | composite Gauss-Legendre on $[0, 28/5]$ (91 panels of 30 nodes); error by the derivative form of the Gauss-Legendre remainder with Cauchy estimates; Hermite values by Arb's Hermite function; $P$ from the explicit coefficients of $He_j$; $b_n$ by Horner in $\rho$ | complex ball quadrature on $[-16, 16]$ |
| largest remainder bound on an $A_{ab}$ | $2.9\cdot10^{-23}$ | $1.6\cdot10^{-28}$ | |
| $b_1$ | $0.8818508164233383791955\ldots$ | $0.8818508164233383791955\ldots$ | $0.88185081642333837919\ldots$ |
| $\sum_{3\le n\le301}\lvert b_n\rvert$ | $\le 2.916517129\cdot10^{-8}$ | $\le 2.916516517\cdot10^{-8}$ | $< 2.916517\cdot10^{-8}$ |

Also $R = \sum|c_j| = 0.99997620388569307924355896 < 1$ and $\gamma = 0.88182581642334060474\ldots < 0.881825817$.

## The tail: the integral $I$ and $B_3$

Both programs build the integrand of Section 5.2 from the paper's definitions, each in its own way: program A with the
paper's recurrence for $F_{ij} = \varphi_r^{-1}\partial_u^i\partial_v^j\varphi_r$, Bell polynomials by Faà di Bruno from
enclosures of $P^{(m)}$ (the exact Taylor shift at the box centre plus remainder) and Table 1; program B with the closed
sum $F_{ij} = \sum_k\binom ik\binom jk k!\,(r/d)^k G_{i-k}(\ell_u, -1/d)G_{j-k}(\ell_v, -1/d)$ and Bell coefficients by series
composition. **Each program checks its integrand against its own head** through the identity $D^3H(t) = T_{-t}\Phi_3(t)$
of Section 4 at real $t$ (left side $\sum n^3 b_n t^n$, right side a Gaussian expectation by the trapezoidal rule):
A: $0.4409254082$ both ways at $t = 0.5$, $-0.2645552449$ at $t = -0.3$; B: $0.3527403266$ at $t = 0.4$,
$-0.5291104899$ at $t = -0.6$; and the derivative part $E_s[\partial_{xy}\Phi_3] = \frac{d}{ds}E_s[\Phi_3]$ (A: $0.08510424$
against $0.08510427$; B: $0.04646985$ against $0.04646986$). Table 1 was also re-derived by hand from
$H(t) = \mathbb E_{\mathrm{corr}=-t}[K_{\rho(t)}(P(X), P(Y))]$: all nine entries agree.

**The partition.** Program A partitions $[0,8]\times[-8,8]\times[0,1]$ (times 2 by the reflection $(x,y)\mapsto(-x,-y)$)
adaptively from the paper's initial $8\times16\times64$ boxes, each split chosen among the three bisections, with an
upper bound of the integrand on every box in ball arithmetic (the quadratic form bounded below by the larger of the split
form, both coefficients at least $1/(1+R)$, and $\min(A\pm B)(u^2+v^2)$; $P$ exact by monotonicity on $|x| \ge 11/2$).
After 1,000,000 splits: **1,008,192 boxes**. Program B checks that they tile the domain exactly (the bisections inside every
initial box reconstructed, dyadic volumes) and bounds every box with its own integrand. Exterior ($|x| > 8$ or $|y| > 8$),
each program with its own constants: $c_8 = 2.797909143\cdot10^{-8}$, $c_8 8^{11} = 240.3385654$, $C_0 \le 8.39\cdot10^{15}$,
$C_1 \le 5.24\cdot10^{22}$, bound $< 10^{-5823}$ (the paper: $c_8 > 2.7979091434\cdot10^{-8}$, $C_0 < 10^{16}$, $C_1 < 10^{23}$,
$< 10^{-5000}$).

| on the same 1,008,192 boxes | program A | program B |
|---|---|---|
| $I \le$ | $288.5461013$ | $214.7490041$ |
| $B_3 \le \sqrt{C_a I}$ ($C_a = 50.1729413\ldots$) | $120.32126$ | $103.80072$ |
| $B_3 < 107$ (Lemma 5.1 as printed) | not with A's bound | **yes** |
| tail $\le B_3/\sqrt{10\cdot301^5}$ | $2.4206174\cdot10^{-5}$ | $2.0882578\cdot10^{-5}$ |

Program B's box bounds are tighter than A's (a different organisation of the same integrand); both are rigorous upper
bounds. For orientation only (not rigorous): midpoint grids give $I \approx 1.93$, so the true $B_3$ is about 10.

## The closing inequality

| | program A | program B |
|---|---|---|
| $\gamma + \sum_{3\le n\le301}\lvert b_n\rvert + \text{tail}$ | $0.881850051762316\ldots$ | $0.881846728166\ldots$ |
| $b_1 >$ | $0.881850816423338$ | $0.881850816423338$ |
| margin | $7.6466\cdot10^{-7} > 0$ | $4.0883\cdot10^{-6} > 0$ |

So both programs certify $\gamma + \sum_{n\ge3}|b_n| < b_1$, hence $K_G^{\mathbb R} \le 1.7813$ given the paper's analytic lemmas;
Lemma 5.1 holds as printed with program B's $B_3 \le 103.81$, and program A's own bound $B_3 \le 120.33$ suffices for the theorem.

**Red controls**, each refused by both programs: the claim $K_G \le 1.7812$ (margin $-4.9\cdot10^{-5}$, A; $-4.5\cdot10^{-5}$, B);
$c_1$ lowered by $10^{-4}$ (heads recomputed: $b_1 = 0.88175187\ldots$, margin $-1.6\cdot10^{-4}$); the sign of $\kappa_{11}$ in
Table 1 flipped (the identity $D^3H = T_{-t}\Phi_3$ fails: $0.4439$ against $0.4409$, and $0.3535$ against $0.3527$).

**Not checked (read, not re-proved):** Lemmas 2.1, 2.2, Proposition 2.3, Lemma 4.1 (the weighted trace inequality),
Proposition 4.2, the identity $D^3H = T_{-t}\Phi_3$ of Section 4 in general (checked numerically at real $t$ only), and
the continuity argument on the unit circle.

Reproduce: `PYTHON=venv/bin/python ./run_all.sh` with python-flint 0.6.0 (about 60 minutes on a loaded 18-core Mac, two
processes; the partition alone about 47 minutes). Output: `results/c10ap26.txt`.
