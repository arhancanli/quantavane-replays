# 1d: the 51,200-cell step function of PR #372 (16 November 2026)

[REV2026], $C_{1d} \le 1.4540711000681048$\* ([#372](https://github.com/teorth/optimizationproblems/pull/372)), replayed
by two separately written exact programs that agree on every number:

- `c1d372_a.py` (A): Python integers. The whole sequence $c_m = \sum_{i+j=m} v_i v_j$ ($0 \le m \le 2n-2$) from one
  squaring $X^2$, $X = \sum_i v_i 2^{bi}$, read back as balanced base-$2^b$ digits; the values at the maximum and at 40
  other knots summed again term by term.
- `c1d372_b.py` with `c1d372_b.c` (B): its own JSON scanner in C, then plain loops over every knot in signed 128-bit
  integers, after checking $\sum_i v_i^2 < 2^{125}$ (by Cauchy-Schwarz this bounds every $\lvert c_m\rvert$ and every
  partial sum; the guard $n \max\lvert v_i\rvert^2 < 2^{125}$ of folder `c1d371` does not hold here, since
  $\max\lvert v_i\rvert = 2^{56}$).

Both print the fingerprint $\sum_m c_m 3^m \bmod 2^{61}-1$ of all 102,399 knot values: they agree. No program of the
source was read or run; only the certificate file was downloaded.

| item | result |
|---|---|
| certificate | ArthurAnache/Hackathon-AI-Reasearch @ `c0f8d18`, `records/1d/autonomous/c006_31bbf628_int56.json`, SHA-256 `dbca822b…02d3`, as stated in #372 |
| cells | $n = 51{,}200$ on $[-1/4, 1/4]$; 14,167 negative heights; largest $\lvert v_i\rvert = 2^{56}$; $S = \sum v_i = 19122158522613018945$ |
| maximum | $\max_m \lvert c_m\rvert = 5192296860720458560490724514234172$, at $m = 12511$ only, where $c_m > 0$ |
| exact value | $2n \max_m \lvert c_m\rvert / S^2 = 21267647941510998263770007610303168512/14626277862557668614837367632277156521 = 1.454071100068104835482988\ldots$, the stated fraction |
| printed | 1.4540711000681048 is $3.55 \cdot 10^{-17}$ **below** the exact value, so as an upper bound it should read 1.4540711000681049 (which holds) |

**Red controls**, each refused by both programs (both give the same changed fraction): height 25000 plus one unit;
the bound claimed as 1.4540711000681047; height 0 negated; the last height plus $10^6$ units ($S$ changes).

Reproduce: `./fetch.sh && ./run_all.sh` (Python 3.9+ and a C compiler; under ten seconds). Output: `results/c1d372.txt`.
