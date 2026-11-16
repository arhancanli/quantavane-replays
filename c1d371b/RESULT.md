# 1d: the eighth certificate of PR #371, $C_{1d} \le 1.4546676818$\* (16 November 2026)

PR [#371](https://github.com/teorth/optimizationproblems/pull/371) now records $C_{1d} \le 1.4546676818$\* (commit
`75912a8`, 15 November 19:46 UTC), with the certificate `C1d_327680_1.4546676818.txt` at intervalplace/frontier @
`1ab20c7` (SHA-256 `a02e9cb8…027a`, as stated). Our folder `c1d371` replayed the first certificate (1.4546733294);
this folder replays the current one with the same two separately written programs (A: one Kronecker squaring in
Python integers; B: plain loops over every knot in C, signed 128-bit integers, guard $n\max\lvert v_i\rvert^2 < 2^{125}$).
`verify.py`, `verify_c1d.c` and the search code were neither read nor run.

| item | result |
|---|---|
| cells | $n = 327{,}680$; 89,427 negative heights; largest $\lvert v_i\rvert = 10^{12}$; $S = 1150626355561932$ |
| maximum | $\max_m \lvert c_m\rvert = 2938681335289198603141896$, at $m = 116123$ only, where $c_m > 0$ |
| exact value | $40122795831148524928230686720/27582104377369449232322657763 = 1.454667681704099986987\ldots$, the stated fraction |
| printed | 1.4546676818 holds ($9.59\cdot10^{-11}$ above the exact value) |

Both programs give the fingerprint $\sum_m c_m 3^m \bmod 2^{61}-1 = 394736550967477101$ of all 655,359 knot values.
Red controls, each refused by both: height 50000 plus one unit; the bound claimed as 1.4546676817; height 0 negated;
the last height plus $10^6$ units.

Reproduce: `./fetch.sh && ./run_all.sh` (Python 3.9+ and a C compiler; a few minutes on two cores). Output:
`results/c1d371b.txt`.
