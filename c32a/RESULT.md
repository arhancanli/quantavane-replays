# 32a, one-shot channel simulation: the closed form of [Li25] (15 November 2026)

A closed form named in the row, recomputed and checked against the page's decimal; the coding theorem behind it is not
checked. No program of any source was read or run. Row: [Li25]
$\sum_{k\ge1}2^{-k-1}k\log_2 k - \log_2\log_2 e \approx 0.76$.

`c32a_a.py`: terms $k \le 400$ in integer fixed point ($2^{-200}$) with its own log; tail for $k > 400$ at most
$2\cdot401^2\,2^{-402}$ ($k\log_2 k \le k^2$, consecutive terms decrease by a factor below 0.51). `c32a_b.py`: terms
$k \le 300$ with the decimal module's ln (correctly rounded, widened by one unit in the last place) and directed
rounding; tail at most $301^2\,2^{-301}$.

| | program A (the total also from program B: first 20 digits identical, checked by `run_all.sh`) |
|---|---|
| $\sum_{k\ge1}2^{-k-1}k\log_2 k$ | $1.28853127577938852965\ldots$ |
| $\log_2\log_2 e$ | $0.52876637294489758911\ldots$ |
| [Li25] | $0.75976490283449090729\ldots \approx 0.76$: holds |

**Red controls**, each refused by both programs: $2^{-k}$ in place of $2^{-k-1}$ (2.048…); the value printed as 0.75.

Reproduce: `./run_all.sh` (Python 3.9+, standard library only; under a minute). Output: `results/c32a.txt`.
