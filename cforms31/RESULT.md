# Closed forms on eight pages that had no replay folder: 6a, 11b, 35a, 88b, 23b, 46a, 56a, 61a (16 November 2026)

Closed forms named in the rows, against the decimals the pages print. Three are irrational (6a, 11b, 35a); the other
fifteen rows are rationals, checked as exact arithmetic. None of them is a certificate; the theorems behind the rows
are not checked. No program of any source was read or run.

`cforms_a.py`: integer fixed point at $2^{-200}$ with its own square root and logarithm (atanh series with a geometric
tail), exact Fractions for rationals. `cforms_b.py`: decimal intervals at 60 digits (correctly rounded sqrt and ln,
widened by one unit), rationals compared through integers. Both classify each printed decimal (truncation, correct
rounding, rounded up), and `run_all.sh` checks that the two programs print the same 30 digits for every row.

| page | row | closed form | value | page | relation |
|---|---|---|---|---|---|
| 6a | [AHS2022], [CL2022], [P2022], [S2022] | $(3-\sqrt5)/2$ | $0.38196601125010515179\ldots$ | $0.381966\ldots$ | truncation (lower bound: safe) |
| 11b | [KP2020] | $\log_2(3/2)$ | $0.58496250072115618145\ldots$ | $\approx 0.58496$ | correct rounding |
| 35a | [AP24] | $\log_2(1+\sqrt2)$ | $1.27155330316361197264\ldots$ | $\approx 1.271$ | truncation, not the rounding (lower bound: safe) |
| 88b | [Pol14b] | $4 - 28/157$ | $3.82165605095541401273\ldots$ | $3.821656\ldots$ | truncation |
| 23b | [AJM2025] | $17/92$ | $0.18478260869565217391\ldots$ | $\approx 0.18478$ | correct rounding |
| 46a | [Bo1991], [Wo1995]/[MVV1996], [TVV1998], [TV2000] | $4 - 2/15$, $4-2/11$, $4-2/9$, $4-2/7$ | | $3.86667$, $3.81818$, $3.77778$, $3.71429$ | correct roundings |
| 46a | [Tao2003]/[BG2011], [Gut2016], [Wan2022], [WW2022], [WW2024] | $10/3$, $13/4$, $3+3/13$, $3+3/14$, $22/7$ | | $3.33333$, $3.25$, $3.23077$, $3.21429$, $3.142857$ | correct roundings ($13/4$ exact) |
| 56a | [LRS1999], [KS2003] | $\frac12 - \frac1{2^2+1}$, $7/64$ | $0.3$, $0.109375$ | $0.3$, $0.109375$ | exact |
| 61a | [LRS1995], [KRS2003] | $\lambda_1 \ge \frac14 - \theta^2$ | $3/16 \leftrightarrow \frac14$, $975/4096 \leftrightarrow 7/64$ | | exact |

For information (the pages write these with $\approx$, as approximations): on the upper-bound pages 11b and 46a,
$0.58496$, $3.81818$, $3.33333$ and $3.142857$ round down (they are correct roundings, slightly below the bounds).

**Red controls**, refused by both: $(3+\sqrt5)/2$ for $(3-\sqrt5)/2$; $\log_2(3/2)$ printed as 0.58497; $4 - 2/13$ for
$4-2/15$; $\lambda_1 = \frac14 + \theta^2$ for $\frac14 - \theta^2$.

Reproduce: `./run_all.sh` (Python 3.9+, standard library only; a few seconds). Output: `results/cforms31.txt`.
