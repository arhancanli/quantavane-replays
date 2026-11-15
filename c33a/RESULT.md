# 33a, the Ihara constant $A(2)$: the finite inputs of [DM2013] $A(2) \ge 0.316999\ldots$ (15 November 2026)

[DM2013] (I. Duursma, K.-H. Mak, *On lower bounds for the Ihara constants A(2) and A(3)*, Compositio Math. 149
(2013); arXiv 1102.4127v4, TeX source `iharaconstants.tex`, SHA-256 of the e-print `92602f0d…db46`, read as text)
proves $A(2) \ge 0.316999\ldots$ (Theorem 1.3 in the arXiv numbering) by the class field tower criterion of its
Theorem 1.1 applied to one explicit function field (Proposition 5.5):

- $H = \mathbb F_2(x, y)$, $y^2 + (x^3 + x + 1)y = x^2 + x$, stated: $(a_d(H)) = (6, 0, 1, 1, 6, \ldots)$, $g(H) = 2$;
- $k$ = the compositum of $v^2 + (x^5 + x^2 + 1)v = (x^2 + x)h$, $h$ in the span of $\{1, x, x^2, y, y^2\}$, an
  elementary abelian extension of degree 32 with conductor $2P_5 + 2P_5'$ in which the 6 rational places split; stated:
  $(a_d(k)) = (192, 0, 0, 0, 2, 16, 0, 16, 0, 64, \ldots)$, $g(k) = 343$;
- $S$ = 2 places of degree 5, 16 of degree 6, 15 of degree 8, 4 of degree 10, $\nu = 2$, $t = 192$; inequality (1.1)
  holds; with $f' = 2(2\cdot5(1-2^{-5}) + 16\cdot6(1-2^{-6}) + 15\cdot8(1-2^{-8}) + 4\cdot10(1-2^{-10}))$,
  $A(2) \ge 192/(343 - 1 + f'/2) \ge 0.316999\ldots$.

Two separately written programs check every finite statement; nothing is downloaded and no program of any source
was read or run.

- `c33a_a.py`: $\mathbb F_{2^d}$ by carry-less products modulo the smallest irreducible polynomial; points of $H$ over
  $\mathbb F_{2^d}$, $d \le 10$ (each affine point checked nonsingular; two rational places at infinity); for each of
  the 31 quadratic subfields $C_h$ the trace test, then $\#k = \sum_h \#C_h - 30\,\#H$ and places by Moebius inversion;
  the conductor of every character at the zeros of $A = x^5 + x^2 + 1$ from Laurent series in $t = x - x_0$ (the
  pole of $f_h = (x^2+x)h/A^2$ reduced by $\sqrt{c_{-2}}/t$: ramified iff $c_{-1} + \sqrt{c_{-2}} \ne 0$, exponent 2);
  $g(k)$ by the conductor-discriminant formula.
- `c33a_b.py`: $\mathbb F_{2^d}$ by logarithm tables over the largest primitive polynomial, the trace as a linear form,
  $w^2 + w = c$ by Gaussian elimination; closed points of $H$ of degree $e \le 10$ and the Frobenius vector
  $(\mathrm{Tr}\, f_h(P))_h$ of each (32 places of degree $e$ above $P$ if it is 0, 16 of degree $2e$ otherwise);
  conductors from Laurent series in the local parameter $A$ itself (series reversion); the genus independently from the
  zeta functions: $g(H) = 2$ is the only genus whose L-polynomial, solved from $\#H(\mathbb F_{2^d})$ with the
  functional equation, predicts the counts to $d = 14$, and for every one of the 31 characters the L-polynomial has
  degree 22 (not 20 or 24) by the same test, so $2g(k) = 2\cdot 2 + 31\cdot 22$.

| statement | both programs |
|---|---|
| $(a_d(H))_{d \le 10}$ | $(6, 0, 1, 1, 6, 13, 12, 33, 56, 105)$; $g(H) = 2$ |
| the 6 rational places of $H$ split completely in $k$ | yes ($a_1(k) = 192$) |
| conductor $2P_5 + 2P_5'$ for each of the 31 nontrivial characters | yes (310 of 310 point-character pairs in A, 62 of 62 place-character pairs in B; 31 of 31 L-polynomials of degree 22 in B) |
| $(a_d(k))_{d \le 10}$ | $(192, 0, 0, 0, 2, 16, 0, 16, 0, 64)$, as stated; $g(k) = 343$ |
| $S$ exists in $k$ | 2 of 2, 16 of 16, 15 of 16, 4 of 64 places |
| Theorem 1.1: $t \le \sum f_{\mathfrak p} = 266$; inequality (1.1) | $192 \le 266$; $75^2 - 2\cdot 2252 - 4\cdot 266 = 57 \ge 0$ |
| the bound | $192/(342 + f'/2) = 24576/77527 = 0.31699923897480877629\ldots$ (the paper's $0.316999\ldots$); Corollary 1.2's cruder form gives $6/19$, as stated |

**Red controls**, each refused by both programs: the bound claimed as $A(2) \ge 0.317$; $y^2$ replaced by $x^3$ in the
span of $h$ (the place counts change: $a_8, a_9, a_{10} = 48, 128, 256$); $A = x^5 + x^3 + 1$, also irreducible, in
place of $x^5 + x^2 + 1$ (its zeros are inert in $H$: no places of degree 5 for $S$); $S$ with 17 places of degree 6
(only 16 exist).

**Not checked** (read): Theorem 1.1 itself (Kuhnt's generator and relation rank bounds and the Golod-Shafarevich
step), that the compositum is the ray class field the paper names, and the average-conductor estimate behind $f'$
(each place of $S$ is ramified in at most a fraction $1 - 2^{-f}$ of the characters).

Reproduce: `./run_all.sh` (Python 3.9+, standard library only; under a minute). Output: `results/c33a.txt`.
