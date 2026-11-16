# 82a: the certificate of [Gri26], $C_{82} \le 0.2536331090204145$ (16 November 2026)

[Gri26] (Grinsztajn, PR [#87](https://github.com/teorth/optimizationproblems/pull/87), merged): Doche's construction
with denominator $Q = Q_1 Q_2 R_0 R_2 P_7 P_9$ ($R_0 = Q_1 - P_1^5P_2^5P_4P_8$, $R_2 = Q_2 + P_1^5P_2^5P_4P_7$, degree
140) and numerators $P_1, \dots, P_6, P_8$ with rational exponents. Data: maaxgrin/zhang-zagier-c82-bound @ `cb88f92`,
the record `qstar_lp_grid32768_maxden1e7` (line 1 of `runs/qstar_targeted_refined.jsonl`), which is the candidate its
certified run `runs/qstar_cert_roots_262144_target253635.jsonl` names; `fetch.sh` extracts its `denominator` and
`numerator` fields into `certificate.json` (SHA-256 in `SHA256SUMS`). No program of the repository was read or run.

Two separately written programs, the same two as for PR #280's v6 certificate:
program 1 (`../c82/c82_check.py`, through `c82gri_first.py`, which only registers the row) and program 2
(`../c82v6/c82v6_second.py`). Hypotheses H1-H3 (monic denominators, $\deg A = 0.8946\,D < D$, coprimality of all 42
pairs), every root of every $F\circ T$ certified, $\log M(Q\circ T) = 35.4992299049\ldots$, the integral of
$\max(0, g)$ by adaptive panels with Jensen self-tests:

| program | $B$ |
|---|---|
| 1 | $[0.2536269201610361,\ 0.2536269202551648]$ |
| 2 | $[0.2536269162787913,\ 0.2536269293400451]$ |

Both below the row's $0.2536331090204145$ (margin $6.2\cdot10^{-6}$; the row is the repository's certified interval
bound, which is looser than $B$). For information: the README's numerical value 0.253626920959 is $7\cdot10^{-10}$ above
our enclosure. Red controls refused (program 1): a bound just below our enclosure; $P_7$ added as a numerator (caught
by H3); $P_7$ doubled (not monic, H1); exponents scaled so that $\deg A = D + 1$ (H2); program 2's own controls in
`results/gri26_second.txt`.

Reproduce: `./fetch.sh && ./run_all.sh` (Python 3.9+; about a minute, single-threaded). Outputs in `results/`.
