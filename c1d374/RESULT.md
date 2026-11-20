# 1d: the 102,400-cell step function of PR #374 (20 November 2026)

[WEILL2026], $C_{1d} \le 1.4513064006879361$\* ([#374](https://github.com/teorth/optimizationproblems/pull/374)),
replayed by two separately written exact programs that agree on every number:

- `c1d374_a.py` (A): Python integers. Each height is the binary64 number of the JSON file, read by Python's `json` and
  turned into an exact rational by `float.as_integer_ratio` (every height here is an integer). The whole sequence
  $c_m = \sum_{i+j=m} v_i v_j$ ($0 \le m \le 2n-2$) comes from one squaring $X^2$, $X = \sum_i v_i 2^{bi}$, read back as
  balanced base-$2^b$ digits; the values at the maximum and at 40 other knots are summed again term by term.
- `c1d374_b.py` with `c1d374_b.c` (B): its own JSON scanner in C, each number converted by the C library's `strtod`
  and required to be an integer, then plain loops over every knot in signed 128-bit integers, after checking
  $\sum_i v_i^2 < 2^{125}$ (here $\sum_i v_i^2 < 2^{118}$; by Cauchy-Schwarz this bounds every $\lvert c_m\rvert$ and every
  partial sum).

Both print the fingerprint $\sum_m c_m 3^m \bmod 2^{61}-1$ of all 204,799 knot values: they agree. No program of the
source was read or run (`verify_standalone.py`, `independent_check.py`, `sanity_checks.py`, `controls.py`); only the
certificate file was downloaded.

| item | result |
|---|---|
| certificate | cweill/c1d-upper-bound @ `7889e93`, `c1d_n102400_1.4513064007.json`, SHA-256 `391cc726…1e2d`, as stated in #374 |
| cells | $n = 102{,}400$ on $[-1/4, 1/4]$; every height an integer-valued binary64 number; 32,150 negative, none zero; largest $\lvert v_i\rvert = 2^{56}$; $S = \sum v_i = 37625858860555748526$ |
| maximum | $\max_m \lvert c_m\rvert = 10032334463188519210813332564771440$, at $m = 30289$ only, where $c_m > 0$ |
| exact value | $2n \max_m \lvert c_m\rvert / S^2 = 513655524515252183593642627316297728000/353926313748615382700541734286037793169 = 1.451306400687936098797\ldots$, the stated fraction |
| printed | 1.4513064006879361 holds ($1.2 \cdot 10^{-18}$ above the exact value) |

Why the knots suffice: $f*f$ is linear between the knots $t = (m+1)h - 1/2$ ($h = 1/(2n)$), where it equals $h\,c_m$, so
$\max_t \lvert f*f(t)\rvert/(\int f)^2 = 2n \max_m \lvert c_m\rvert / S^2$ (the absolute value inside the maximum, as on
the 1d page).

**Red controls**, each refused by both programs (both give the same changed fraction): height 50000 plus one unit;
the bound claimed as 1.4513064006879360; height 0 negated; the last height plus $10^6$ units ($S$ changes).

Reproduce: `./fetch.sh && ./run_all.sh` (Python 3.9+ and a C compiler; about 15 seconds). Output: `results/c1d374.txt`.
