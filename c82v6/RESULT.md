# C82: PR #280's update of 14 November 2026 (v6), $C_{82} < 0.2519694883$\*, replayed by two separately written programs

Certificate: vidigoat/optimization-records @ `e3bae2e`, `tao/certificates/c82/v6/certificate.json`, SHA-256
`29dcbf21a10ddb600a9866fe8ee7e7b43e4357926d74f374a8b6fbf41ada71b5` (the hash the pull request comment states). Six
monic denominators of degree 104 (those of #330's certificate of 13 November, with new multiplicities), 37 numerators,
$D = 1{,}000{,}000{,}000{,}064$ (as stated). Only this data file was downloaded; the folder's checkers, receipts and
search code were neither read nor run. With $T(z) = z(1 - z)$, $u = e^{2\pi i t}$,
$g(t) = \sum_i q_i \log\lvert S_i(T(u))\rvert - \sum_j m_j \log\lvert R_j(T(u))\rvert$ and
$B = [\sum_j m_j \log M(R_j \circ T) + \int_0^1 \max(0, g)\,dt]/D$.

| | program 1: `../c82/c82_check.py`, unchanged (`c82v6_first.py` only registers the row) | program 2: `c82v6_second.py` (no code shared) |
|---|---|---|
| H1, H2 | integer coefficients, monic denominators, $\deg A = 0.897093836719\,D < D$ | the same, $\deg A = 117583883373992923/131072$ |
| H3: all 222 pairs coprime | gcd modulo $2^{61} - 1$ | Euclid modulo $10^9 + 7$ |
| roots of every $F \circ T$ (43 factors, 3,412 roots) | Gershgorin discs of the Weierstrass matrix after Aberth's iteration; radius at most $2^{-303}$ | Newton inclusion: some root within $N\lvert P(z)/P'(z)\rvert$ of $z$, $P(z)$ and $P'(z)$ evaluated exactly as Gaussian integers, discs pairwise disjoint; radius at most $2^{-247}$ |
| $\sum_j m_j \log M(R_j \circ T)$ | $[251969407876.3573919301, \ldots 9302]$ | $[251969407876.3573919301, \ldots 9301]$ (decimal logarithms) |
| $\int_0^1 \max(0, g)\,dt$ | $[80412.46953410, 80412.46953411]$: per-root series on 3,260 panels (82 crude), depth 9 to 44 | $[80412.46953282, 80412.46953816]$: Horner's value of $g$ at each panel centre, a degree-28 Taylor polynomial in the panel variable from root power sums, Cauchy remainder; 3,342 panels (82 crude), depth 11 to 36 |
| self-tests | Vieta; Horner against the roots at random points; Jensen's formula for the denominators and for 11 numerators | closed forms ($\pi$, $\ln 2$, $e^{2\pi i/6}$, two Mahler measures, two Jensen integrals, a moved root); Jensen's formula for the denominators: $[-251969407876.357395, -251969407876.357388]$ against $-251969407876.357391$ |
| $B$ | $[0.2519694882727008, 0.2519694882727009]$, width $4.8 \cdot 10^{-22}$ | $[0.2519694882727008, 0.2519694882727009]$, width $5.3 \cdot 10^{-18}$ |
| $B < 0.2519694883$ | holds, margin $2.73 \cdot 10^{-11}$ | holds, margin $2.73 \cdot 10^{-11}$ |

So $B = 0.25196948827270087879\ldots$ (program 1's enclosure). The two values the pull request comment states for $B$
are reported, not used: the first run's $B \le 0.25196948827272$ is above our enclosure; the second run's
$B \le 0.25196948827270$ is $8.8 \cdot 10^{-16}$ below it (the fourteen digits of $B$, truncated). The printed bound
0.2519694883 does not depend on either.

**Red controls**, refused by both programs: the bound just below the enclosure of $B$ ($0.2519694882726$) is not
implied; the first denominator added as a numerator with $q = 1/1000$ (H3); the first denominator doubled, not monic
(H1); every $q_i$ scaled so that $\deg A = D + 1$ (H2). Program 2's self-tests also refuse a root approximation moved
by $2^{-100}$ (its disc grows to radius at least $2^{-101}$) and two equal approximations (the discs overlap).

**What rests on mathematics rather than computation:** that $C_{82} \le B$ for every certificate meeting H1 to H3
(the certificates' notes, read by hand, not re-derived); Gershgorin's theorem (program 1); the Newton inclusion and
Cauchy's estimate (program 2); the series bounds of both. Program 2 takes the decimal module's logarithms as correctly
rounded (its specification) and widens each by one unit in the last of 60 digits.

Reproduce: `./fetch.sh && ./run_all.sh` (Python 3.9+, standard library only; about 15 minutes on two cores; program 1
uses two worker processes, program 2's integral two processes at once; every process below 20 minutes of CPU, the
largest about 4 minutes). Output: `results/pr280f_first.txt`, `results/pr280f_second.txt`, `results/pr280f.txt`.
