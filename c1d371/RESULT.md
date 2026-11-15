# 1d: the 327,680-cell step function of PR #371 (15 November 2026)

[Swap2026], $C_{1d} \le 1.4546733294$\* ([#371](https://github.com/teorth/optimizationproblems/pull/371)), replayed by
two separately written exact programs that agree on every number:

- `c1d371_a.py` (A): Python integers. The whole sequence $c_m = \sum_{i+j=m} v_i v_j$ ($0 \le m \le 2n-2$) from one
  squaring $X^2$, $X = \sum_i v_i 2^{bi}$ (a signed integer, packed as $P - N$ from the positive and negative
  heights), read back as balanced base-$2^b$ digits; the values at the maximum and at 40 other knots summed again
  term by term.
- `c1d371_b.py` with `c1d371_b.c` (B): its own reader in C, then plain loops over every knot in signed 128-bit
  integers (each pair $i < j$ once, doubled, plus the square), after checking that $n \max|v_i|^2 < 2^{125}$. No product
  of big integers and no transform.

Both print the fingerprint $\sum_m c_m 3^m \bmod 2^{61}-1$ of all 655,359 knot values: they agree. No program of the
source was read or run; only the certificate file was downloaded.

| item | result |
|---|---|
| certificate | intervalplace/frontier @ `316b76c`, `c1d-signed-autocorrelation/C1d_327680_1.4546733294.txt`, SHA-256 `66eec291…0df3`, as stated in #371 |
| cells | $n = 327{,}680$ on $[-1/4, 1/4]$; 87,672 negative heights; largest $\lvert v_i\rvert = 10^{12}$; $S = \sum v_i = 1332927688303720$ |
| maximum | $\max_m \lvert c_m\rvert = 3943653272890604779247222$, at $m = 99803$ only, where $c_m > 0$ |
| exact value | $2n \max_m \lvert c_m\rvert / S^2 = 8076601902879958587898310656/5552175694520934184034268245 = 1.454673329385129058565\ldots$, the stated fraction |
| printed | 1.4546733294 holds ($1.49 \cdot 10^{-11}$ above the exact value) |

The reduction (the page's definition): for a step function with $n$ equal cells of width $h = 1/(2n)$, $f*f$ is linear
between knots, where it equals $h c_m$, so $\max_t \lvert f*f\rvert$ is attained at a knot, and with $\int f = hS$ the
ratio is $2n \max_m \lvert c_m\rvert / S^2$.

**Red controls**, each refused by both programs (and both give the same changed fraction): height 50000 plus one
unit; the bound claimed as 1.4546733293; height 0 negated; the last height plus $10^6$ units ($S$ changes).

**Not checked:** that this step function is the limit of the authors' search (irrelevant to the bound).

Reproduce: `./fetch.sh && ./run_all.sh` (Python 3.9+ and a C compiler; about four minutes on two cores, under
200 MB). Output: `results/c1d371.txt`.
