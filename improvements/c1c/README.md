# C1c >= 0.4103780000 (Barnard-Steinerberger autocorrelation constant, lower bound)

A new lower bound for $C_{1c}$ of the
[optimization constants repository](https://teorth.github.io/optimizationproblems/constants/1c.html) (10 Nov 2026,
evening):

$$C_{1c} \geq \frac{3385460305559548364109249137976803328}{8249614514518895716878560524318391809} = 0.41037800003882762970\ldots \geq 0.4103780000.$$

| | value | our margin |
| --- | --- | --- |
| this certificate (`c1c_lower_K1048576.json`, K = 2^20) | 0.41037800003882762970... | |
| our 10 Nov morning certificate (`c1c_lower_K262144.json`, #344) | 0.41037741298866168317... (0.4103774129) | +5.871e-7 |
| best other open claim, #354 (SeanSagait, K = 524288, 10 Nov) | 0.41037727125141974306... (0.4103772) | +7.288e-7 |
| open claim, #270 (vidigoat, `c1c_v5`, 9 Nov) | 0.41036083390473056108... (0.4103608339) | +1.717e-5 |
| best on the live page, [Kos2026] | 0.4103333 | +4.470e-5 |

## The certificate

`c1c_lower_K1048576.json`: $f = \sum_i w_i 1_{[i/K,(i+1)/K)}$ with $K = 1048576$ and $1048577$ non-negative integer
weights (support $[0, 1048577/1048576)$). Its autocorrelation is linear between multiples of $1/K$, so

$$\min_{0 \le t \le 1} \int f(x) f(x+t)\,dx \Big/ \Big(\int f\Big)^2 = K \min_{0 \le k \le K} A_k / S^2, \qquad A_k = \sum_i w_i w_{i+k},\ S = \sum_i w_i,$$

the rule the 1c page uses. The minimum is attained only at $k = 179183$, with $A_k = 3228626542624996532544373643853$
and $S = 2872214218076168447$; all 880,597 lags from $167980$ to $K$ lie within $1.7 \cdot 10^{-13}$ of the minimum (measured exactly from the NTT digits; 607,400 within $10^{-13}$), and every lag below $167980$ gives at least $0.41037828$. The shape: one spike cell $w_0 = 1152921504606846848$ ($2^{60} - 128$)
carrying 40.14% of the mass, $w_i = 0$ for $1 \le i < 167980$ ($t_1 = 0.160198$), and a positive profile on
$167980 \le i \le K$. The file also states the reduced fraction (`score_num`, `score_den`) and the printed claim
(`claim`).

## Two exact checks that share no code

    python3 -I check.py c1c_lower_K1048576.json --record 0.4103774129 --record 846362574176502890154236447416909824/2062400481577886062273091058271244329 --record 0.4103772 --record 151911829776891904/370176031712590489 --record 0.4103608339 --record 2556315663995378551725293568/6229433836731244234885554241 --record 0.4103333 --controls
    python3 ntt_check.py c1c_lower_K1048576.json --fraction=3385460305559548364109249137976803328/8249614514518895716878560524318391809 --claim=0.4103780000 --controls

- `check.py` (Python standard library only, about four minutes for the bound, about twenty with the red controls): every $A_k$ from one big-integer product (Kronecker
  substitution), the whole sequence confirmed by its generating-polynomial identity at random points modulo
  $2^{61}-1$, the deciding lag summed again directly; the result must equal the stated fraction, imply the claim and
  beat every `--record` strictly. Red controls, all refused: one value raised by one unit, the largest value raised by
  1% (with the old fraction, and with its own, where the claim fails), the claim moved to 0.4103780001.
  Output: `results/check_K1048576.txt`.
- `ntt_check.py` (numpy, about half a minute): every $A_k$ modulo six NTT primes below $2^{31}$ by number-theoretic
  transforms (length $2^{22}$), Garner's algorithm for the mixed-radix digits (the product of the primes exceeds
  $n \max w^2$, which is checked), the least $A_k$ found by comparing digit vectors, then summed again directly; a
  self-test against brute force on random inputs runs first. Red controls, all refused: a weight raised by 1 changes
  the fraction, and the claim moved one unit up. Output: `results/ntt_K1048576.txt`.

Both give the same fraction, the same deciding lag and the same $A_k$.

## How it was found, and credit

The structure is that of the open pull request [#354](https://github.com/teorth/optimizationproblems/pull/354)
(SeanSagait, 9-10 Nov 2026): a spike at the left end, $f = 0$ up to $t_1 \approx 0.16$, then a profile with every lag
in $[t_1, 1]$ active, solved by a damped fixed-point iteration. We followed its description only (its programs were
not read or run) and made two changes: the spike is a single cell rather than a geometric run of cells, and the
weights are scaled to about $2^{60}$ (big integers), so rounding costs almost nothing. Our solver
(`../search/c1c_spike/fp.py`, unchanged since the 10 Nov morning certificate): with $T = 1$ and spike weight $M$, the
active-lag equations are $M g_k + \sum_i g_i g_{i+k} = T$ for $k \ge t_1 K$, so $g_k = (T - c_k)/M$ with $c$ the
autocorrelation of the profile $g$, iterated with adaptive damping and negative values clipped ($t_1$ is where $c$
crosses $T$); $M$ by golden-section search. This certificate is `python ladder.py 1048576` (858 s wall on one core,
$M = 641.639118$, $M/\sqrt{KT} = 0.626601$) followed by `python mkcert.py 1048576 OUT.json`. The unrounded values of
the family at $K = 2^{14}, 2^{16}, 2^{18}, 2^{20}$ are $0.41036567, 0.41037506, 0.41037741, 0.41037800$, which
follows $0.4103782 - 0.205/K$ (#354 reports the same law for its family): the family's limit is near $0.4103782$, so
larger $K$ adds at most about $2 \cdot 10^{-7}$ more (extrapolation, not claimed).

Prepared with AI assistance (Claude, Anthropic).

## Our earlier certificates, superseded (both remain valid lower bounds)

10 Nov 2026, morning (in #344 until this update):

$$C_{1c} \geq \frac{846362574176502890154236447416909824}{2062400481577886062273091058271244329} = 0.41037741298866168317\ldots \geq 0.4103774129$$

`c1c_lower_K262144.json`: $K = 262144$, $262145$ weights, one spike cell $w_0 = 2^{60}$ with 40.14% of the mass, zero
to $t_1 = 0.160198$, minimum only at $k = 219901$ ($A_k = 12914467989753767244785102041884$,
$S = 2872212026698506854$); all 220,150 lags from $41995$ to $K$ lie within $4 \cdot 10^{-14}$ of the minimum.
Checks: `python3 -I check.py c1c_lower_K262144.json --record 0.4103768 --record 935577600963641344/2279800930147535089 --record 0.4103608339 --record 2556315663995378551725293568/6229433836731244234885554241 --record 0.4103598711 --record 11005597691381428229825143056897540096/26819380902505338613619516881482225625 --record 0.4103333 --controls`
and `python3 ntt_check.py c1c_lower_K262144.json --fraction=846362574176502890154236447416909824/2062400481577886062273091058271244329 --claim=0.4103774129 --controls`
(outputs `results/check_K262144.txt`, `results/ntt_K262144.txt`).

8 Nov 2026:

$$C_{1c} \geq \frac{11005597691381428229825143056897540096}{26819380902505338613619516881482225625} = 0.41035987114651620459\ldots \geq 0.4103598711$$

`c1c_lower_K98304.json`: $K = 98304$, $98368$ non-negative integer weights (largest $2^{56} - 8$), minimum only at
$k = 98304$ ($A_k = 111954729119684125059256419442724$, $S = 5178743177886439525$).
Check: `python3 -I check.py c1c_lower_K98304.json --record 0.4103596512 --record 0.4103538918 --record 0.4103597993 --record 0.4103333 --controls`
(output `results/check.txt`). Found from vidigoat's $K = 98304$ certificate of
[#270](https://github.com/teorth/optimizationproblems/pull/270) (`c1c_v3`, sha256 `e0786a59...caec66`), which refines
Julsgaard's [#318](https://github.com/teorth/optimizationproblems/pull/318), by L-BFGS-B on a log-sum-exp soft minimum
in cycles of re-softening and re-sharpening (search code in [`../search`](../search)).
