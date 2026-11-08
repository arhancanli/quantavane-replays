# C1d <= 1.454683534422819 (signed third autocorrelation constant, upper bound)

A new upper bound for $C_{1d}$ of the
[optimization constants repository](https://teorth.github.io/optimizationproblems/constants/1d.html):

$$C_{1d} \leq \frac{3563627596515371110610821706288771891200}{2449761416959550533000893662063862027289} = 1.45468353442281856471\ldots \leq 1.454683534422819.$$

| | value | our margin |
| --- | --- | --- |
| this certificate | 1.45468353442281856471... | |
| best open claim, #312 (Julsgaard) | 1.45472766452875064884... (1.454727664528751) | -4.413e-5 |
| open claim, #345 (our first certificate, commit 36b3a05) | 1.45469141638562121644... (1.454691416385622) | -7.882e-6 |
| best on the live page, [SCIEL2026] | 1.457854884921989 | -3.171e-3 |

## The certificate

`c1d_upper_n163840.json`: a step function with $n = 163840$ equal cells on $[-1/4, 1/4]$ (cell width $h = 1/327680$) and signed integer cell values $v_i$ (largest absolute value $2^{56} - 8$). As on the 1d page, $f*f$ is linear between the knots $mh$, where it equals $h \sum_{i+j=m-1} v_i v_j$, so

$$\max_t \lvert f*f(t) \rvert \Big/ \Big(\int f\Big)^2 = 2n \max_m \lvert c_m \rvert / S^2, \qquad c = v * v,\ S = \sum_i v_i.$$

The maximum is attained only at $m = 50729$ (0-based index of $c$), with $|c_m| = 10875328358506381563143376789211340$ and $S = 49495064571728265733$. The file also states the reduced fraction (`score_num`, `score_den`) and the printed claim (`claim`).

## Re-derive

Python 3.9 or later, standard library only:

    python3 -I check.py c1d_upper_n163840.json --record 1.454727664528751 --record 9900630255081376929482866688/6805830738283663732311831361 --record 1.454691416385622 --record 1091712020077839913516641060217129861120/750476704392981840387650570565448814041 --record 1.457854884921989 --controls

`check.py` computes every $c_m$ exactly from big-integer products (Kronecker substitution, with the positive and negative parts packed separately: $c = pp - 2pq + qq$), confirms the whole sequence by its generating-polynomial identity at random points modulo $2^{61}-1$, sums the deciding knot again directly, and requires the result to equal the stated fraction, to imply the claim, and to beat every `--record` strictly. The red controls must all be refused: one value raised by one unit (its fraction no longer matches), the largest value raised by 1% (with the old fraction, and again with its own fraction, where the claim fails), and the claim moved one unit of its last digit in the improving direction (1.454683534422818). The output is in `results/check.txt`.

## How it was found, and credit

The warm start is Julsgaard's 81,920-cell certificate from [#312](https://github.com/teorth/optimizationproblems/pull/312) (sha256 `04c9c75909b70c215c312e85ef14679d90bc7f98257749c75f422d33048f7237`, 1.454727664528751), itself refined from the 10,240-cell function of [#269](https://github.com/teorth/optimizationproblems/pull/269). Each cell was split in two (163,840 cells, the same function), L-BFGS on $\log \lVert v*v \rVert_p - 2 \log S$ polished it at $p = 10^8, \ldots, 10^{11}$, and then came restarts of the $p$-norm continuation that #312 describes, each from the best vector so far with $p$ doubled up to about $7 \cdot 10^{10}$ (76,000 to 85,000 L-BFGS iterations each): from $p = 32$ on a Kaggle T4 GPU (8 Nov 2026, about 1.3 h; 1.4546914164), then on one CPU core from $p = 32$ (1.4546868875), $p = 16$ (no gain) and $p = 64$ (this certificate), about 2.5 h. The search code is in [`../search`](../search). Only the published values of #312 were used; no other contributor's checker or verification code was read or run.

Prepared with AI assistance (Claude, Anthropic).
