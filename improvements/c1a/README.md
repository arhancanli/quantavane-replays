# C1a <= 1.5027427280 (Sidon set autocorrelation constant, upper bound)

A new upper bound for $C_{1a}$ of the [optimization constants repository](https://teorth.github.io/optimizationproblems/constants/1a.html) (10 Nov 2026, evening):

$$C_{1a} \leq \frac{29131355358281158426614181216278902472704}{19385457547511099408999611434350080635225} = 1.50274272798999968357\ldots \leq 1.5027427280.$$

| | value | our margin |
| --- | --- | --- |
| this certificate (`c1a_upper_n262144_v2.json`) | 1.50274272798999968357... | |
| vidigoat, #270, 10 Nov (`c1a_v3`, refined from our #346) | 1.50274273405473036507... (1.5027427341) | -6.065e-9 |
| our 8 Nov certificate (#346, below) | 1.50274273747657840978... (1.5027427375) | -9.487e-9 |
| best on the live page, [A2026] | 1.5027434868 | -7.588e-7 |

## The certificate

`c1a_upper_n262144_v2.json`: $n = 262144$ equal cells on $[-1/4,1/4]$ with integer values $v_i$ (non-negative), the largest of size $2^{56}$. For a step function on $n$ equal cells, $f*f$ is linear between multiples of $1/(2n)$ with knot values proportional to $c = v * v$, so the bound is $2n \max_m |c_m| / S^2$ with $S = \sum_i v_i$ (the page's rule). The maximum is attained only at m = 303135. The file states the reduced fraction (`score_num`, `score_den`) and the printed claim (`claim`).

## Two exact checks that share no code

    python3 -I check.py c1a_upper_n262144_v2.json --record 1.5027427341 --record 116537598371962553887369519059399403896832/77549933019818026777245548252121251545225 --record 1.5027427375 --record 29134478977604335569142246661390387511296/19387536037290895570476550697799766901521 --controls
    python3 ntt_check_conv.py c1a_upper_n262144_v2.json --fraction=29131355358281158426614181216278902472704/19385457547511099408999611434350080635225 --claim=1.5027427280 --controls

- `check.py` (Python standard library only; the same file as in `../c1c`): every $c_m$ from big-integer products (Kronecker substitution, signed values split into positive and negative parts), the whole sequence confirmed by its generating-polynomial identity at random points modulo $2^{61}-1$, the deciding knot summed again; the value must equal the stated fraction, imply the claim and beat every `--record` strictly. Red controls, all refused. Output: `results/check_v2.txt`.
- `ntt_check_conv.py` (numpy, new on 10 Nov): every $c_m$ modulo six NTT primes below $2^{31}$, Garner's algorithm for the mixed-radix digits, the signed value recovered from the symmetric residue range ($2n \max|v|^2$ below the product of the primes is checked), $\max_m |c_m|$ over exact integers and the deciding knot summed again, after a self-test against brute force on signed and non-negative inputs. Red controls, all refused: one value moved by one unit where $|c_m|$ grows at the deciding knot raises the bound, and the claim moved one unit down. Output: `results/ntt_v2.txt`. It also reproduces vidigoat's C1d fraction exactly (positive control).

Both give the same fraction and the same deciding knot.

## How it was found, and credit

Warm start: vidigoat's 10 Nov certificate in [#270](https://github.com/teorth/optimizationproblems/pull/270) ([vidigoat/optimization-records @ a559680e](https://github.com/vidigoat/optimization-records/tree/a559680e416be0eae03955c58372efff4a551fbb/tao/certificates), `c1a_v3/c1a_n262144_v3.json`, sha256 `b32a20bfba61aabf7f35ae3887a371c1c22f200082b835467deed4e6d37f5c14`), which had refined our earlier certificate; only its published values were used, its programs were neither read nor run. Search: `../search/c1ad_modal/c1ad_retake2.py` on Modal, 16 independent one-core calls (8 per constant), each L-BFGS-B ($v \ge 0$) on $\log \lVert v*v \rVert_p - 2 \log S$ with $p$ doubled from its own restart value $p_0$ (from $3 \cdot 10^7$ to $10^9$) up to $1.1 \cdot 10^{12}$, 3000 iterations a stage, the best vector checkpointed to a Modal volume after every improving stage. This certificate is job `1a_3` ($p_0 = 3 \cdot 10^7$), rounded to integers with the largest value $2^{56}$; about 55 minutes on one core.

Prepared with AI assistance (Claude, Anthropic).

## Our earlier certificates (superseded; they remain valid bounds)

A new upper bound for $C_{1a}$ of the
[optimization constants repository](https://teorth.github.io/optimizationproblems/constants/1a.html):

$$C_{1a} \leq \frac{29134478977604335569142246661390387511296}{19387536037290895570476550697799766901521} = 1.50274273747657840978\ldots \leq 1.5027427375.$$

| | value | our margin |
| --- | --- | --- |
| this certificate | 1.50274273747657840978... | |
| best open claim, #270 (vidigoat, `c1a_n262144`) | 1.50274335217817824811... (1.5027433522) | -6.147e-7 |
| best on the live page, [A2026] | 1.5027434868 | -7.493e-7 |

## The certificate

`c1a_upper_n262144.json`: a step function with $n = 262144$ equal cells on $[-1/4, 1/4]$ (cell width $h = 1/524288$) and non-negative integer cell values $v_i$ (largest $2^{56}$). As on the 1a page, $f*f$ is linear between the knots $mh$, where it equals $h \sum_{i+j=m-1} v_i v_j$, so

$$\max_t f*f(t) \Big/ \Big(\int f\Big)^2 = 2n \max_m c_m / S^2, \qquad c = v * v,\ S = \sum_i v_i.$$

The maximum is attained only at $m = 198994$ (0-based index of $c$), with $c_m = 222278434582552609017503712931750393$ and $S = 278478265128831162578$. The file also states the reduced fraction (`score_num`, `score_den`) and the printed claim (`claim`).

## Re-derive

Python 3.9 or later, standard library only:

    python3 -I check.py c1a_upper_n262144.json --record 1.5027433522 --record 3756858380445445962687721679242267976826786797829957539923492864000000000000/2500000000000000227858304266031237715040682299126697825903032255969919203361 --record 1.5027434868 --controls

`check.py` computes every $c_m$ exactly from big-integer products (Kronecker substitution), confirms the whole sequence by its generating-polynomial identity at random points modulo $2^{61}-1$, sums the deciding knot again directly, and requires the result to equal the stated fraction, to imply the claim, and to beat every `--record` strictly. The red controls must all be refused: one value raised by one unit (its fraction no longer matches), the largest value raised by 1% (with the old fraction, and again with its own fraction, where the claim fails), and the claim moved one unit of its last digit in the improving direction (1.5027427374). The output is in `results/check.txt`.

## How it was found, and credit

The warm start is vidigoat's 262,144-cell certificate from [#270](https://github.com/teorth/optimizationproblems/pull/270) (`c1a_n262144.json`, sha256 `1848a92032f54a5336acc1a2f78d2f04bed45c0eea2e1802b929a02ff17bf3de`, 1.5027433522), refined there from Atris's 65,536-step function [A2026], itself refined from the AgentDiscover construction [EA2026]. From it, on one CPU core (8 Nov 2026, about 4.6 h), L-BFGS-B with $v \geq 0$ on $\log \lVert v*v \rVert_p - 2 \log S$ (float64 FFTs): a polish at $p = 10^{10}, 10^{11}, 10^{12}$, then four cycles of restarts of the $p$-norm continuation from the best vector at $p_0 = 10^8, 10^7, 10^9, 10^6$, with $p$ doubled up to about $5 \cdot 10^{11}$ to $8 \cdot 10^{11}$ (from about 300 to 29,000 L-BFGS-B iterations per restart). The gains came from $p_0 = 10^7$ and above all from $p_0 = 10^6$. The search code is in [`../search`](../search). Only the published values of #270 were used; no other contributor's checker or verification code was read or run.

Prepared with AI assistance (Claude, Anthropic).
