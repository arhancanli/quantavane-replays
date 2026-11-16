# 71a: Hod's $\beta(1/2)$ [Hod2017] (16 November 2026)

Hod (arXiv 1711.00762, Theorem 4.4): any constant in the Fourier entropy/influence conjecture is at least
$\beta(1/2) > 6.4547837$, with $\beta(z) = 4z(1-z)\sum_{k\ge-2} h(q_{k+2})\pi_k$, $h$ the binary entropy,
$q_0 = z$, $q_{m+1} = 1/(1+q_m)$, $\pi_m = q_1\cdots q_m$, $\pi_{-1} = 1/z$, $\pi_{-2} = 1/(1-z)$ (Section 4, read through
the paper's HTML rendering; the theorem itself is read, not checked). Two separately written programs:
`c71hod_a.py` (the recursion in exact rationals, logarithms in fixed-point integers at $2^{-256}$) and `c71hod_b.py`
(the paper's Fibonacci closed forms $q_m = F_{m+2}/F_{m+3}$, $\pi_m = 2/F_{m+3}$ at $z = 1/2$, checked against the
recursion; 70-digit decimals). Tails bounded geometrically ($h \le 1$, $q_k \le 2/3$).

Both: $\beta(1/2) = 6.45478371656275855381683622859\ldots$, so the row's $> 6.4547837$ holds (margin $1.66\cdot10^{-8}$).
Red controls, refused by each: the claim 6.4547838; $z = 0.499$ ($\beta = 6.45449\ldots$) and $q_{m+1} = 1/(2+q_m)$ (A);
an off-by-one in the Fibonacci indices ($6.3839\ldots$) and the factor $4z(1-z)$ taken as 1/2 (B).

Reproduce: `./run_all.sh` (Python 3.9+; seconds). Output: `results/c71hod.txt`.
