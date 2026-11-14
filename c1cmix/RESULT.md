# 1c: the mixture certificate behind PR #354's remark on test-measure upper bounds (14 November 2026)

Pull request: [#354](https://github.com/teorth/optimizationproblems/pull/354) (SeanSagait), commit `f089e846` of the
pull request (11 November 2026), which adds a remark to the comments of the 1c page: two step functions $f_1, f_2$
($K = 4096$) and rational $\lambda_1 + \lambda_2 = 1$ satisfy
$\min_{[0,1]} (\lambda_1 f_1 \star f_1 + \lambda_2 f_2 \star f_2) = 0.4106833868\ldots$ exactly (each $f_i$ normalised
to $\int f_i = 1$), so for every probability measure $\nu$ on $[0, 1]$, $\sup_f \int (f \star f)\,d\nu \ge 0.41068$.
It changes no bound of the page.

## Inputs

Fetched read-only from SeanSagait/c1c-dual-bound-limit @ `0b6cbddd0566387110ddd1e44a0b9c1f6e3ed765` (data only):

| file | SHA-256 (as in the repository's `SHA256SUMS`) |
| --- | --- |
| `mix_K4096_P2.json` (two rows of 4097 non-negative integer weights, $\lambda = 8502051/2^{24},\ 8275165/2^{24}$) | `4d1edd69f7230fae5d6ec7184c58fdfdba2c5d7d23a2110f9be842d35683137c` |
| `mix_K1024_P2.json` (the smaller certificate) | `cb16c841bd3bc87e160cea598fa88c2b369f426f7f7d501508dd6b516b033b21` |

The repository's README was read for the statement. `verify_mix.py` and `verify_kronecker.py` were neither read nor
run.

## What was checked, and how

For a step function on equal steps of width $1/K$ the autocorrelation is linear between multiples of $1/K$, and at
$t = k/K$ its normalised value is $K a_k / S^2$ with $a_k = \sum_j w_j w_{j+k}$, $S = \sum_j w_j$. So the minimum over
$[0, 1]$ of the mixture is $\min_{0 \le k \le K} \sum_i \lambda_i K a_k^{(i)} / S_i^2$.

- `c1cmix_direct.py` (A): every $a_k$ by the direct double sum in Python integers; the mixture in `Fraction`.
- `c1cmix_limbs.py` (B): every weight split into 16-bit limbs, every limb correlation by numpy in exact 64-bit
  integers (each below $2^{45}$), recombined in Python integers; lags 0, 1 and $K$ compared with direct sums; the
  mixture over one common denominator in integers only, reduced at the end.

Both also check that the weights are non-negative integers, that there are at least $K + 1$ cells (so no $a_k$ with
$k \le K$ is an empty sum), and that $\lambda_1, \lambda_2 > 0$ with $\lambda_1 + \lambda_2 = 1$.

## Result

| certificate | $S_1$, $S_2$ | exact minimum (both programs) | where |
| --- | --- | --- | --- |
| $K = 4096$ | 281474976710653, 281474976710639 | 5279542067372831214509056145021140361894623145945986164923459 / 12855504354070095327289033191535703182400661039125234952783872 = 0.410683386817204933762798330888... | only at $k = 780$ |
| $K = 1024$ | 17592186044395, 17592186044419 | 23008408314619471245228187581841051992114046895188148613 / 56045551208694960130369998216028958643061427698579763200 = 0.410530502750233008530285641491... | only at $k = 202$ |

The $K = 4096$ fraction is the one stated in the repository's README, with the minimum at $k = 780$ as stated; the
remark's $0.4106833868\ldots$ and $0.41068$ are truncations of it. The $K = 1024$ value is the README's
$0.4105305028$.

## Red controls (the value differs in both programs, and the two programs still agree with each other)

1. The weight at index 780 of the first component plus 1.
2. The two $\lambda_i$ exchanged: the minimum becomes $0.3997242\ldots$, at $k = 3994$.
3. The first weight of the second component plus 1.

Full output: `results/mix.txt`.

## What was NOT checked

- The inference from the exact minimum to the statement about every probability measure $\nu$ (three lines in the
  README: $\sup_f \int (f \star f)\,d\nu \ge \sum_i \lambda_i \int (f_i \star f_i)\,d\nu \ge \min_{[0,1]} \sum_i
  \lambda_i f_i \star f_i$) is mathematics, not computation; we read it and it is elementary, but no program checks it.
- The piecewise linearity of the autocorrelation of an equal-step function (Kravitz, arXiv:2004.06611, proof of
  Theorem 3.1) is used, not re-proved.
- The remark's reading of what this limit means for the literature's upper-bound arguments, and the suggested limit
  near 0.41073 as $K$ grows.

## Reproduce

`./fetch.sh && ./run_all.sh` (Python 3.9+; program B needs numpy; about ten seconds).
