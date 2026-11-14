# 49a: the 252,303-set block of PR #370, $\mu_3 \ge 1.5678643$\* [Pin2026] (14 November 2026)

The row of [#370](https://github.com/teorth/optimizationproblems/pull/370): Koshta's 28-point two-pendant block
([Kos2026], replayed in `../c49` on 4 November) with its 21 extra sets replaced by 32 seven-point sets, 16 inside each
pendant copy of $H$. Certificate: rpinho/math-certificates @ `a50021b`, `c49/round1/certificates/best_1.5678643845728069793620.json`,
SHA-256 `53381b5f52264c43d56343ab27d698c598462d5abfd83c043365afbae9efb1e8` (the packet's `SHA256SUMS` gives this
value), read as data. No program of the packet (`checker.py`, the referee's C++ and the rest) was read or run; its
README was read as documentation.

Two separately written exact programs:

- **A** = our checker of 4 November, `../c49/sunflower_check.py`, unchanged, through `c49p370_a.py` (which only reads
  the JSON and registers the row's numbers): the box method, exhaustive over all pairs and triples, on two point
  partitions (core 8 + pendants 10 + 10, and 4 + 4 + 10 + 10); port counts; the Perron root as an exact surd; the
  rate enclosed by exact integer comparisons.
- **B** = `c49p370_b.py`, new, no shared code: the family as a quasi-reduced decision diagram over the 28 points
  (1,059 nodes), then memoised joins over the diagram for every pair (Sperner) and every triple of members
  (a triple of distinct sets is a 3-sunflower exactly when no point lies in exactly two of them); the rate as the
  larger root of $x^2 - \mathrm{tr}\,x + \det$, by integer comparisons. Its self-test compares both joins with the
  definitions on 3,000 random small families (both answers occur).

| check | A | B |
|---|---|---|
| 252,303 distinct members on 28 points, ports 0 and 1 | yes | yes |
| Sperner and 3-sunflower-free, exhaustively | yes (box method, two partitions) | yes (pair join 17,001 states; triple join 965,747 states) |
| port counts $N_{00}, N_{01}, N_{10}, N_{11}$ | 33980, 64728, 64728, 88867 | the same |
| Perron root | $64728 + \sqrt{3019700660}$, $3019700660 = 88867 \cdot 33980$, irrational | the same |
| $\rho^{1/26}$ | $[1.5678643845728069793620334, \ldots335)$ | $[1.5678643845728069793620, \ldots621)$ |
| printed 1.5678643 | holds, the truncation at 7 decimals | holds |

So the row's $(64728 + \sqrt{3019700660})^{1/26} = 1.56786438\ldots$ and the description's $1.567864384572806979\ldots$
are exact truncations, and the printed 1.5678643 is a valid lower bound. Program A also checks, as a second path, that
the certificate is Koshta's block minus the 21 sets and plus the 32 sets that `best_delta.json` lists (the 20 sets
inside the first copy of $H$ and the set with all of the second copy removed; 32 seven-point sets avoiding both ports
added, 16 per pendant), and scans every pair and triple through the 32 added sets against all members.

**A fact about the row text (no effect on the bound):** the row gives the certificate's hash as `53381b5f...e1e8`; the
file's SHA-256 (and the packet's `SHA256SUMS`) ends in `...b1e8`.

**Red controls**, all refused. A: a planted 3-sunflower comparable with no member, and a planted containment, each
caught by both box partitions; $N_{11} + 1$ moves the rate; $N_{11}$ lowered by 0.1% refutes 1.5678643; 1.5678644
refuted; without the added set 64768 (a valid family, $N_{00} = 33979$) the rate is 1.5678639771... and 1.5678643 is
refuted. B: a planted set that makes a 3-sunflower with two members and is comparable with none (only the triple
join fires); a member with one point removed; the added set 64768 removed; 1.5678644.

**Not checked (mathematics, read):** the chain construction and the transfer-matrix count of [IRSI2026] (the gluing
lemma of [DEGKM1997], which `../c49` tests on small instances), which turn a valid block into the bound
$\rho^{1/26}$ with $26 = 28 - 2$ new points per copy.

Reproduce: `./fetch.sh && ./run_all.sh` (Python 3.9+, standard library only; about two minutes, one process).
Output: `results/pin370.txt`.
