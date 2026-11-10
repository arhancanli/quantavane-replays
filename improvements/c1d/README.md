# C1d <= 1.4546733414 (signed third autocorrelation constant, upper bound)

A new upper bound for $C_{1d}$ of the [optimization constants repository](https://teorth.github.io/optimizationproblems/constants/1d.html) (10 Nov 2026, evening):

$$C_{1d} \leq \frac{838861227937471890002426104328885370880}{576666392427609677715360509915372739121} = 1.45467334138563686456\ldots \leq 1.4546733414.$$

| | value | our margin |
| --- | --- | --- |
| this certificate (`c1d_upper_n163840_v2.json`) | 1.45467334138563686456... | |
| vidigoat, #270, 10 Nov (`c1d`, refined from our #345) | 1.45467334576751329963... (1.4546733458) | -4.382e-9 |
| our 9 Nov certificate (#345, below) | 1.45467334742095245159... (1.454673347420953) | -6.035e-9 |
| best on the live page, [SCIEL2026] | 1.457854884921989 | -3.182e-3 |

## The certificate

`c1d_upper_n163840_v2.json`: $n = 163840$ equal cells on $[-1/4,1/4]$ with integer values $v_i$, the largest of size $2^{56}$. For a step function on $n$ equal cells, $f*f$ is linear between multiples of $1/(2n)$ with knot values proportional to $c = v * v$, so the bound is $2n \max_m |c_m| / S^2$ with $S = \sum_i v_i$ (the page's rule). The maximum is attained only at m = 52929. The file states the reduced fraction (`score_num`, `score_den`) and the printed claim (`claim`).

## Two exact checks that share no code

    python3 -I check.py c1d_upper_n163840_v2.json --record 1.4546733458 --record 646233848324886328286774558720/444246710235775351931155445761 --record 1.454673347420953 --record 167771889099128792433071411029129035776/115333032942810266163095027979519535125 --controls
    python3 ntt_check_conv.py c1d_upper_n163840_v2.json --fraction=838861227937471890002426104328885370880/576666392427609677715360509915372739121 --claim=1.4546733414 --controls

- `check.py` (Python standard library only; the same file as in `../c1c`): every $c_m$ from big-integer products (Kronecker substitution, signed values split into positive and negative parts), the whole sequence confirmed by its generating-polynomial identity at random points modulo $2^{61}-1$, the deciding knot summed again; the value must equal the stated fraction, imply the claim and beat every `--record` strictly. Red controls, all refused. Output: `results/check_v2.txt`.
- `ntt_check_conv.py` (numpy, new on 10 Nov): every $c_m$ modulo six NTT primes below $2^{31}$, Garner's algorithm for the mixed-radix digits, the signed value recovered from the symmetric residue range ($2n \max|v|^2$ below the product of the primes is checked), $\max_m |c_m|$ over exact integers and the deciding knot summed again, after a self-test against brute force on signed and non-negative inputs. Red controls, all refused: one value moved by one unit where $|c_m|$ grows at the deciding knot raises the bound, and the claim moved one unit down. Output: `results/ntt_v2.txt`. It also reproduces vidigoat's C1d fraction exactly (positive control).

Both give the same fraction and the same deciding knot.

## How it was found, and credit

Warm start: vidigoat's 10 Nov certificate in [#270](https://github.com/teorth/optimizationproblems/pull/270) ([vidigoat/optimization-records @ a559680e](https://github.com/vidigoat/optimization-records/tree/a559680e416be0eae03955c58372efff4a551fbb/tao/certificates), `c1d/c1d_n163840.json`, sha256 `b05de6abc8e3d7a2253a0887b6ec8d0eec58beb32e4eb76ceb98df66e26377b6`), which had refined our earlier certificate; only its published values were used, its programs were neither read nor run. Search: `../search/c1ad_modal/c1ad_retake2.py` on Modal, 16 independent one-core calls (8 per constant), each L-BFGS on $\log \lVert v*v \rVert_p - 2 \log S$ with $p$ doubled from its own restart value $p_0$ (from $3 \cdot 10^7$ to $10^9$) up to $3.3 \cdot 10^{12}$, 3000 iterations a stage, the best vector checkpointed to a Modal volume after every improving stage. This certificate is job `1d_3` ($p_0 = 3 \cdot 10^7$), rounded to integers with the largest value $2^{56}$; about 55 minutes on one core.

Prepared with AI assistance (Claude, Anthropic).

## Our earlier certificates (superseded; they remain valid bounds)

A new upper bound for $C_{1d}$ of the
[optimization constants repository](https://teorth.github.io/optimizationproblems/constants/1d.html):

$$C_{1d} \leq \frac{167771889099128792433071411029129035776}{115333032942810266163095027979519535125} = 1.45467334742095245159\ldots \leq 1.454673347420953.$$

This certificate (9 November 2026) supersedes the one of 8 November in this folder
([c40060f](https://github.com/arhancanli/quantavane-replays/tree/c40060ff2d03ed022d970ef25047e70c23a54ad3/improvements/c1d),
$1.454683534422819$), which remains a valid bound, as does the first one
([36b3a05](https://github.com/arhancanli/quantavane-replays/tree/36b3a05e7d966ef9288507ecb17065323adc217e/improvements/c1d),
$1.454691416385622$).

| | value | our margin |
| --- | --- | --- |
| this certificate | 1.45467334742095245159... | |
| our 8 Nov certificate (c40060f; #345 as of 8 Nov) | 1.45468353442281856471... (1.454683534422819) | -1.019e-5 |
| best open claim by others, #312 (Julsgaard) | 1.45472766452875064884... (1.454727664528751) | -5.432e-5 |
| best on the live page, [SCIEL2026] | 1.457854884921989 | -3.182e-3 |

## The certificate

`c1d_upper_n163840.json`: a step function with $n = 163840$ equal cells on $[-1/4, 1/4]$ (cell width $h = 1/327680$) and signed integer cell values $v_i$ (largest absolute value $2^{56}$). As on the 1d page, $f*f$ is linear between the knots $mh$, where it equals $h \sum_{i+j=m-1} v_i v_j$, so

$$\max_t \lvert f*f(t) \rvert \Big/ \Big(\int f\Big)^2 = 2n \max_m \lvert c_m \rvert / S^2, \qquad c = v * v,\ S = \sum_i v_i.$$

The maximum is attained only at $m = 47433$ (0-based index of $c$), with $|c_m| = 10239983465522997585026331239570864$ and $S = 48027707199659295450$. The file also states the reduced fraction (`score_num`, `score_den`) and the printed claim (`claim`).

## Re-derive

Python 3.9 or later, standard library only:

    python3 -I check.py c1d_upper_n163840.json --record 1.454683534422819 --record 3563627596515371110610821706288771891200/2449761416959550533000893662063862027289 --record 1.454727664528751 --record 9900630255081376929482866688/6805830738283663732311831361 --record 1.454691416385622 --record 1.457854884921989 --controls

`check.py` (unchanged) computes every $c_m$ exactly from big-integer products (Kronecker substitution, with the positive and negative parts packed separately: $c = pp - 2pq + qq$), confirms the whole sequence by its generating-polynomial identity at random points modulo $2^{61}-1$, sums the deciding knot again directly, and requires the result to equal the stated fraction, to imply the claim, and to beat every `--record` strictly. The red controls must all be refused: one value raised by one unit (its fraction no longer matches), the largest value raised by 1% (with the old fraction, and again with its own fraction, where the claim fails), and the claim moved one unit of its last digit in the improving direction (1.454673347420952). The output is in `results/check.txt`.

## How it was found, and credit

The warm start is Julsgaard's 81,920-cell certificate from [#312](https://github.com/teorth/optimizationproblems/pull/312) (sha256 `04c9c75909b70c215c312e85ef14679d90bc7f98257749c75f422d33048f7237`, 1.454727664528751), itself refined from the 10,240-cell function of [#269](https://github.com/teorth/optimizationproblems/pull/269). Each cell was split in two (163,840 cells, the same function), L-BFGS on $\log \lVert v*v \rVert_p - 2 \log S$ polished it at $p = 10^8, \ldots, 10^{11}$, and then came restarts of the $p$-norm continuation that #312 describes, each from the best vector so far with $p$ doubled up to about $5$ to $7 \cdot 10^{10}$: on 8 Nov from $p = 32$ on a Kaggle T4 GPU, then $p = 32, 16, 64$ on one CPU core (the 8 Nov certificate); on 9 Nov, fourteen more on one CPU core (`../search/worker_cont_1d_numpy.py`, about 7.9 h), with $p_0$ cycling through $8, 128, 24, 48, 12, 96$: gains at the restarts $p_0 = 128, 24, 48, 96$ and $128$, the last giving this certificate. Only the published values of #312 were used; no other contributor's checker or verification code was read or run.

Prepared with AI assistance (Claude, Anthropic).
