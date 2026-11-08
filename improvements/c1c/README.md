# C1c >= 0.4103598711 (Barnard-Steinerberger autocorrelation constant, lower bound)

A new lower bound for $C_{1c}$ of the
[optimization constants repository](https://teorth.github.io/optimizationproblems/constants/1c.html):

$$C_{1c} \geq \frac{11005597691381428229825143056897540096}{26819380902505338613619516881482225625} = 0.41035987114651620459\ldots \geq 0.4103598711.$$

| | value | our margin |
| --- | --- | --- |
| this certificate | 0.41035987114651620459... | |
| best open claim, #270 (vidigoat, update of 5 Nov, `c1c_v3`) | 0.41035965124998351121... (0.4103596512) | +2.199e-7 |
| open claim, #318 (Julsgaard) | 0.41035389181157810991... (0.4103538918) | +5.979e-6 |
| open claim, #344 (our first certificate, commit 36b3a05) | 0.41035979932958291478... (0.4103597993) | +7.182e-8 |
| best on the live page, [Kos2026] | 0.4103333 | +2.657e-5 |

## The certificate

`c1c_lower_K98304.json`: $f = \sum_i w_i 1_{[i/K,(i+1)/K)}$ with $K = 98304$ and $98368$ non-negative integer weights $w_i$ (largest $2^{56} - 8$; support $[0, 98368/98304)$). Its autocorrelation is linear between multiples of $1/K$, so

$$\min_{0 \le t \le 1} \int f(x) f(x+t)\,dx \Big/ \Big(\int f\Big)^2 = K \min_{0 \le k \le K} A_k / S^2, \qquad A_k = \sum_i w_i w_{i+k},\ S = \sum_i w_i,$$

the rule the 1c page uses. The minimum is attained only at $k = 98304$, with $A_k = 111954729119684125059256419442724$ and $S = 5178743177886439525$. The file also states the reduced fraction (`score_num`, `score_den`) and the printed claim (`claim`).

## Re-derive

Python 3.9 or later, standard library only:

    python3 -I check.py c1c_lower_K98304.json --record 0.4103596512 --record 640823574552569898220412928/1561614482809354047911104009 --record 0.4103538918 --record 2582127846607416062611587072/6292441470965870058830439361 --record 0.4103597993 --record 3669060891957752860604242716590637056/8941082674160596251510892487009553075 --record 0.4103333 --controls

`check.py` computes every $A_k$ exactly from one big-integer product (Kronecker substitution), confirms the whole sequence by its generating-polynomial identity at random points modulo $2^{61}-1$, sums the deciding lag again directly, and requires the result to equal the stated fraction, to imply the claim, and to beat every `--record` strictly. The red controls must all be refused: one value raised by one unit (its fraction no longer matches), the largest value raised by 1% (with the old fraction, and again with its own fraction, where the claim fails), and the claim moved one unit of its last digit in the improving direction (0.4103598712). The output is in `results/check.txt`.

## How it was found, and credit

The warm start is vidigoat's $K = 98304$ certificate from [#270](https://github.com/teorth/optimizationproblems/pull/270) (`c1c_v3`, sha256 `e0786a59bd9d8598d9f4a86fcee40575ea108a97b184b332ef05bb6764caec66`, 0.4103596512), which refines Julsgaard's [#318](https://github.com/teorth/optimizationproblems/pull/318), which refines earlier certificates of #270, Koshta [Kos2026] and rozetyp [R2026]. From it, L-BFGS-B on a log-sum-exp soft minimum of $K A_k / S^2$ over all lags, in cycles of re-softening and re-sharpening as #270 describes (20,000 iterations per sharpness stage, float64 FFTs): on a Kaggle T4 GPU (8 Nov 2026, about 2.6 h) one and a half cycles with sharpness $3 \cdot 10^7, 10^8, \ldots, 10^{10}$ gave 0.4103597993, and on one CPU core two more cycles (sharpness $10^8, \ldots, 10^{10}$, then $3 \cdot 10^7, 10^8, \ldots, 10^{10}$; about 1.5 h) gave this certificate. The minimum now sits at the last lag, $k = K$. The search code is in [`../search`](../search). Only the published weights of #270 were used; no other contributor's checker or verification code was read or run.

Prepared with AI assistance (Claude, Anthropic).
