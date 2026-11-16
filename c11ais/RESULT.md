# 11a: $\delta$ of [IS2024], $C_{11a} \le \pi/2 - \delta$ (16 November 2026)

Ivanisvili and Stone (arXiv 2407.04835, Section 2.4): $\delta = \int_{1/2}^{1} \frac{\min\{4p^3(1-p)^3,\ (2p-1)^4p(1-p)\}}{45 - 3(p(1-p))^3}\frac{dp}{\sqrt{p(1-p)}}$.
Their theorem is read, not checked. In the variable $s = \arcsin(2p-1)$ the integrand is a rational function of
$\cos s$, $\sin s$ with the crossover at $\tan^2 s = 1/2$; expanding $1/(45 - 3\cos^6 s/64)$ as a geometric series gives
$\delta$ exactly as $a\,s^* + b\,\pi + c\sqrt2$ plus a tail below $10^{-40}$. Two separately written exact routes:
`c11ais_a.py` (reduction formulas for $\int\cos^n$) and `c11ais_b.py` (cosine-polynomial expansion with Chebyshev
values of $\sin(2js^*)$); they give the same $a, b, c$.

Both: $\delta = 0.00013279130299166442911386578\ldots$, so $\pi/2 - \delta = 1.5706635354919049548\ldots$; the row's
$\delta \approx 0.00013$ is a correct rounding. Red controls, refused by each: 44 for 45; the crossover at
$\sin^2 = 1/4$; the claim $\delta \approx 0.000133$.

Reproduce: `./run_all.sh` (Python 3.9+; seconds). Output: `results/c11ais.txt`.
