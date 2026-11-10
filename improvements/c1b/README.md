# C1b <= 0.3808585498 (Erdős minimum overlap constant, upper bound)

A new upper bound for $C_{1b}$ of the
[optimization constants repository](https://teorth.github.io/optimizationproblems/constants/1b.html) (10 Nov 2026):

$$C_{1b} \leq \frac{1730339319993979413608789703258449}{4543259751217974174964184288067584} = 0.38085854975166311714\ldots \leq 0.3808585498.$$

| | value | our margin |
| --- | --- | --- |
| this certificate (`c1b_upper_n114688.json`) | 0.38085854975166311714... | |
| best on the page, [A2026] (Atris Labs; replayed by [T2026b]) | 0.38085855031642856010... (0.3808585504) | -5.648e-10 |
| [EA2026] (CodexProLong, EinsteinArena, 15 Aug 2026) | 0.38085857485785830914... (0.3808585749) | -2.511e-8 |

The EinsteinArena leaderboard also lists four scores between 0.3808563 and 0.3808566 whose constructions are not
public; we could not check them (that leaderboard's scorer rescales the mass to $n/2$, which can push heights above 1),
so they are not compared here.

## The certificate

`c1b_upper_n114688.json`: $n = 114688$ equal cells of $[-1,1]$ (width $2/n$), $f = a_i/D$ on cell $i$ with
$D = 2^{48}$ and integers $0 \le a_i \le D$, $\sum_i a_i = nD/2$ exactly (so $\int f = 1$), and $g = 1 - f$ on $[-1,1]$.
For two step functions on the same grid the overlap $\int f(t) g(x+t)\,dt$ is linear in $x$ between multiples of
$2/n$, so

$$\sup_x \int f(t)\,g(x+t)\,dt = \frac{2}{n D^2} \max_{|k| < n} X_k, \qquad X_k = \sum_i a_i (D - a_{i+k}),$$

with $a = 0$ outside the cells (the rule of the [A2026] row). The maximum is attained at $k = \pm 11464$ (both shifts
cover the whole support, so $X_k = X_{-k}$ there), with $X_k = 1730339319993979413608789703258449$; since
$n D^2 / 2 = 7 \cdot 2^{109}$, the bound is $X_k / (7 \cdot 2^{109})$. The next largest shift is lower by
$1.3 \cdot 10^{-12}$. Heights: 33,115 cells at 0, 7,705 at 1,
73,868 strictly between, support cells 11759 to 102905. The file also states the reduced fraction (`score_num`,
`score_den`) and the printed claim (`claim`).

## Two exact checks that share no code

    python3 -I check_1b.py c1b_upper_n114688.json --record 0.3808585504 --record 10919976354672639675467789613614419/28672000000000000000000000000000000 --record 0.3808585749 --controls
    python3 ntt_check_1b.py c1b_upper_n114688.json --fraction=1730339319993979413608789703258449/4543259751217974174964184288067584 --claim=0.3808585498 --controls

- `check_1b.py` (Python standard library only, about 7 seconds, controls included about a minute): every $X_k$ from one big-integer product (Kronecker
  substitution of $A(x)$ and the reversed $B(x)$, $b_i = D - a_i$), the whole sequence confirmed by its
  generating-polynomial identity at random points modulo $2^{61}-1$, the deciding shift summed again directly; the
  mass and the box $0 \le a_i \le D$ checked; the result must equal the stated fraction, imply the claim and beat
  every `--record` strictly (the records given are the [A2026] row, printed and exact, and the [EA2026] row). Red
  controls, all refused: one unit of mass moved between two cells where it changes the bound (old fraction stated),
  1% of a cell moved where it raises the bound above the claim (its own fraction stated), the mass broken, a height
  above $D$, the claim moved to 0.3808585497. Output: `results/check_n114688.txt`.
- `ntt_check_1b.py` (numpy, a second): every $X_k$ modulo six NTT primes below $2^{31}$, Garner's algorithm for the
  mixed-radix digits (the product of the primes exceeds $n D^2$, which is checked), the largest $X_k$ found by
  comparing digit vectors and summed again directly; a self-test against brute force on random inputs runs first.
  Red controls, all refused: one unit of mass moved changes the bound, the mass broken, the claim moved one unit down.
  Output: `results/ntt_n114688.txt`.

Both give the same fraction and the same deciding shifts. Both checkers also reproduce the [A2026] construction's
exact value $10919976354672639675467789613614419/28672000000000000000000000000000000$ (positive control).

## How it was found, and credit

Warm start: the [A2026] construction (Atris Labs, [atrislabs/records @ 1fa857a0](https://github.com/atrislabs/records/tree/1fa857a021ddbb4064708daf336758b866f75cb4),
`1b-erdos-minimum-overlap/construction.json`, sha256 `62a976fa1366...`; data only, its scorer was not run), which
refines the CodexProLong construction [EA2026]. Our search (`../search/c1b/c1b_opt.py`, numpy and scipy, one core):
L-BFGS-B on $\log \lVert C \rVert_p$ over the shifts ($C_k = \sum_i h_i (1 - h_{i+k})$, box $0 \le h \le 1$), with a
quadratic penalty keeping $\sum h = n/2$ and the exact mass restored after each stage; $p = 3 \cdot 10^9, 10^{10},
3 \cdot 10^{10}, 10^{11}$, 4000 iterations a stage, about 10 minutes in all
(`python c1b_opt.py construction.json out.json --ps 3e9,1e10,3e10,1e11 --iters 4000`). `mkcert_1b.py` rounds the
heights to multiples of $2^{-48}$ and fixes the mass exactly in integers (rounding moved the value by about $10^{-16}$).
The gain is small ($5.6 \cdot 10^{-10}$): this is a refinement of the [A2026] construction in its own basin, not a
new construction.

Prepared with AI assistance (Claude, Anthropic).
