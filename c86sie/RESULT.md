# 86a: Siegel's bound [Sie45] $\approx 1.73361$ from its definition (16 November 2026)

The bound $e(1+\nu^{-1})^{-\nu}$, $\nu$ the positive root of $(1+\nu)\log(1+\nu^{-1}) + \frac{\log\nu}{1+\nu} = 1$, as
stated in Orloski, Talebizadeh Sardari and Smith (arXiv 2401.03252v3, Section 1.1 and Corollary 1.4; Siegel's 1945
paper was not read). Two separately written programs:

- `c86sie_a.py`: exact rationals; $f < 0$ on $(0, 1/100]$ and $f > 0$ on $[1, \infty)$ by explicit bounds, the sign
  pattern on 2000 cells of $[1/100, 1]$ with a derivative bound and $f' > 0$ on the undecided cells (one root), bisection.
- `c86sie_b.py`: decimal intervals; $f$ increasing on $(0, 1/2]$ and positive on $[1/2, \infty)$, Newton and a certified bracket.

| quantity | value |
|---|---|
| $\nu_0$ | $0.3144808694076345901929969631\ldots$ (the only positive root) |
| $e(1+1/\nu_0)^{-\nu_0}$ | $1.73361051698464779155694826\ldots$ |

The page's $\approx 1.73361$ is a correct rounding; the source's 1.7336105 a truncation. Red controls, refused by
each: the claim 1.73362; the sign before $\log\nu/(1+\nu)$ flipped; $\nu = 1/3$.

Reproduce: `./run_all.sh` (Python 3.9+; under a minute). Output: `results/c86sie.txt`.
