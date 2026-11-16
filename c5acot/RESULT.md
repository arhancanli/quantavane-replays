# 5a: the WT-triple witnesses of [COT26], $C_{5a} \le 0.94304$ (16 November 2026)

Carter, O'Bryant and Tao (arXiv 2610.36207v1): a WT-triple $(\upsilon, \alpha, \beta)$ (Definition 2.2) gives
$\operatorname{diam}(\mathcal A) \ge k^2 - 2\sqrt{\alpha\beta}\,k^{3/2}$ for $k \ge 4(1+\upsilon\sqrt{\beta/\alpha})^2$
(Theorem 2.3, read, not checked), and the row records $\sqrt{\alpha\beta} \le 0.94304$ for the triple
$(10, 0.944814, 0.941259)$ of Lemma 2.2. Its witness (and that of $(3, 1.630563, 0.547610)$, Corollary 1.1) is a random
pair of step functions: $K$ atoms with probabilities $p_i$, window heights $g^{(i)}$ (length $N$) and mirrored table
heights $h^{(i)}$ (length $3N$), exact rationals in dcartermath/sidon-error @ `ff5b032` (`asymptotic/`, `explicit/`;
SHA-256 in `SHA256SUMS`). `verify.py` and the search code were neither read nor run.

Two separately written exact programs check every condition of Definition 2.2: $t = 1$ on $(0, \sigma-\upsilon)$,
$\int w = 1$, $w, t \ge 0$ (strong), and $\mathbb E(w*t) \ge 1$ on $[0, \sigma]$; `c5acot_a.py` evaluates $w*t$ at the
knots by the index formula, `c5acot_b.py` by a generic step-function engine, and also checks that $\mathbb E(w*t)$ is
linear between knots (at all midpoints), so the knot minimum is the minimum.

| witness | $\alpha = \mathbb E\lVert t\rVert^2 - \sigma$ | $\beta = \mathbb E\lVert w\rVert^2$ | $2\sqrt{\alpha\beta}$ | threshold | min $\mathbb E(w*t)$ |
|---|---|---|---|---|---|
| `asymptotic` ($\upsilon = 10$, $K = 6$, $N = 200$) | $0.94481304378915084085\ldots$ | $0.94125802624178782959\ldots$ (the fraction printed in Appendix A.3) | $1.886067719637\ldots \le 1.88607$ | $482.34 \le 483$ | $1$ |
| `explicit` ($\upsilon = 3$, $K = 5$, $N = 32$) | $1.63056262452157843046\ldots$ | $0.54760956513343842563\ldots$ | $1.889880091156\ldots \le 1.8898801$ | $29.9987 < 30$ | $1$ |

So $\sqrt{\alpha\beta} = 0.943033859818\ldots \le 0.94304$, the row's value. Red controls, refused by each: the claim
0.94303; a middle-third table height set to 0.999 and the mirroring ignored (A: $\min \mathbb E(w*t)$ falls to 0.82);
probability moved between atoms and a window height raised (B). (B's first height control changed a zero entry and was
not refused; replaced before any result was kept.)

Reproduce: `./fetch.sh && ./run_all.sh` (Python 3.9+; about half a minute). Output: `results/c5acot.txt`.
