# C1d <= 1.454691416385622 (signed third autocorrelation constant, upper bound)

A new upper bound for $C_{1d}$ of the
[optimization constants repository](https://teorth.github.io/optimizationproblems/constants/1d.html):

$$C_{1d} \leq \frac{1091712020077839913516641060217129861120}{750476704392981840387650570565448814041} = 1.45469141638562121644\ldots \leq 1.454691416385622.$$

| | value | our margin |
| --- | --- | --- |
| this certificate | 1.45469141638562121644... | |
| best open claim, #312 (Julsgaard) | 1.45472766452875064885... (1.454727664528751) | -3.625e-5 |
| best on the live page, [SCIEL2026] | 1.457854884921989 | -3.163e-3 |

## The certificate

`c1d_upper_n163840.json`: a step function with $n = 163840$ equal cells on $[-1/4, 1/4]$ (cell width $h = 1/327680$)
and signed integer cell values $v_i$ (largest $\lvert v_i \rvert = 2^{56}$). As on the 1d page, $f*f$ is linear between
the knots $mh$, where it equals $h \sum_{i+j=m-1} v_i v_j$, so

$$\max_t \lvert f*f(t) \rvert \Big/ \Big(\int f\Big)^2 = 2n \max_m \lvert c_m \rvert / S^2, \qquad c = v * v,\ S = \sum_i v_i.$$

The maximum is attained only at $m = 47352$ (0-based index of $c$), with
$\lvert c_m \rvert = 13326562745090819256794934817103636$ and $S = 54789659768718470858$. The file also states the
reduced fraction (`score_num`, `score_den`) and the printed claim (`claim`).

## Re-derive

Python 3.9 or later, standard library only (about 30 s):

    python3 -I check.py c1d_upper_n163840.json --record 1.454727664528751 \
      --record 9900630255081376929482866688/6805830738283663732311831361 \
      --record 1.457854884921989 --record 1.4630715826075594 --controls

`check.py` computes every $c_m$ exactly from big-integer products (Kronecker substitution, with the positive and
negative parts packed separately: $c = pp - 2pq + qq$), confirms the whole sequence by the identity
$\sum_m c_m r^m = V(r)^2$ at random points modulo $2^{61}-1$, sums the deciding knot again directly, and requires the
result to equal the stated fraction, to imply the claim, and to beat every `--record` strictly. The red controls must
all be refused: one value raised by one unit, the largest value raised by 1% (with the old fraction, and again with its
own fraction, where the claim fails), and the claim lowered by one unit of its last digit (1.454691416385621). The
output is in `results/check.txt`.

## How it was found, and credit

The warm start is Julsgaard's 81,920-cell certificate from
[#312](https://github.com/teorth/optimizationproblems/pull/312) (sha256
`04c9c75909b70c215c312e85ef14679d90bc7f98257749c75f422d33048f7237`, 1.454727664528751), itself refined from the
10,240-cell function of [#269](https://github.com/teorth/optimizationproblems/pull/269). Each cell was split in two
(163,840 cells, same function), then L-BFGS on $\log \lVert v*v \rVert_p - 2 \log S$ polished it at
$p = 10^8, \ldots, 10^{11}$, and one restart of the $p$-norm continuation that #312 describes, from $p = 32$ with $p$
doubled up to about $7 \cdot 10^{10}$ (77,439 L-BFGS iterations, float64 FFTs on a Kaggle T4 GPU, 8 Nov 2026), lowered
the value by $3.62 \cdot 10^{-5}$; the search code is in [`../search`](../search). Only the published values of #312
were used; no other contributor's checker or verification code was read or run.

Prepared with AI assistance (Claude, Anthropic).
