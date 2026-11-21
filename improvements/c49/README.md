# 49a: two more sets in Tilley's 54-point block, $C_{49} \ge 1.5678672030558807$

**The construction is Robert F. Tilley Jr.'s** ([PR #379](https://github.com/teorth/optimizationproblems/pull/379),
20 November 2026; note at [Zenodo 10.5281/zenodo.24405229](https://doi.org/10.5281/zenodo.24405229)), built on Dennis
Koshta's two-pendant block ([PR #260](https://github.com/teorth/optimizationproblems/pull/260)), Sylla's pendants, the
[Cla2026] seed families, IRSI's transfer count and the gluing lemma of Deuber, Erdős, Gunderson, Kostochka and Meyer.
Our replay of his bound is in [`../../frontier/tilley-c49-54point-pr379`](../../frontier/tilley-c49-54point-pr379).

**What we add.** Two sets, in *both* copies of the 28-point block $B$ before the join (labels as in the note: core
edges of $K_5$, ports $p = 01$, $q = 02$, pendants at $e_1 = 03$ and $e_2 = 04$):

- $X_2 = \{12, 14, 23\} \cup$ all ten points of the second pendant (this is the extra set of Koshta's block, which
  Tilley's $B$ omits);
- $X_1 = \{12, 14, 23\} \cup$ all ten points of the first pendant.

Neither contains a port, so each adds to $N_{00}$ of $B$: $B$ becomes 252,273 sets with port counts
33,950 / 64,728 / 64,728 / 88,867. The 54-point family (30,325,618,420 sets) $G \cup (P \times L(A_p)) \cup (L(A_q) \times P)$ stays Sperner
and three-sunflower-free, with

| | Tilley (PR #379) | with $X_1$, $X_2$ |
| --- | --- | --- |
| $N_{00}$ | 4,397,311,488 | 4,397,570,400 |
| $N_{01} = N_{10}$ | 7,211,663,100 | 7,211,840,834 |
| $N_{11}$ | 11,504,366,352 | 11,504,366,352 |
| $\rho = N_{01} + \sqrt{N_{00} N_{11}}$ | $7211663100 + \sqrt{50588282321810251776}$ | $7211840834 + \sqrt{50591260940311180800}$ |
| $\rho^{1/52}$ | 1.56786638820497592522... | **1.56786720305588072742...** |

Tilley's iteration (page 3 of his note) needs only a valid block with the same two ports, so
$C_{49} = \mu_3^S \ge \rho^{1/52} \ge 1.5678672030558807$, which is $8.1 \times 10^{-7}$ above his bound.

**How it was found.** A greedy search (`search/` in our records) over sets "core subset avoiding both ports, each
pendant empty or full", added to both copies and kept when program A found the 54-point family still valid. Of
Pinho's 32 and Koshta's 21 sets, only Koshta's extra set fits beside $P$; of the 255 sets of the second round, only
its mirror image $X_1$ also fits.

## Checks (two separately written exact programs, red controls; a fresh referee)

- Program A (`tilley_a.py`, Python, the two-level box method on our 4 November checker) and program B
  (`tilley_b.mjs`, JavaScript, brute-force pendant patterns, BigInt rate) both find the family with $X_1$ and $X_2$
  Sperner and three-sunflower-free, the counts above, and $\rho^{1/52} \in [1.567867203055880727424, +10^{-21})$.
- Red controls, refused: a third set (Pinho's 64768) added to both copies (A and B find 3-sunflowers); the decimal
  1.5678672030558808 (A: not a valid lower bound).
- Referee (a fresh agent with its own program, written from Tilley's note alone, 21 Nov 19:59 +04): VALID for one set
  ($X_2$, rate 1.56786679563474996619...) and for two sets ($X_1$, $X_2$, rate 1.56786720305588072742...), with every
  count above; its program agrees with brute force on 26,000 random small instances and re-checks every refusal's
  witness against a separate membership oracle. It also confirms that the iteration on page 3 needs only a valid block
  with the same ports and $N_{01} = N_{10}$, $0 < N_{00} \le N_{11}$, and that among all 1,024 sets 'any core subset,
  each pendant empty or full' only $X_1$ and $X_2$ fit (none more on top of both).
- Finite depth (the registry's criterion, the referee's own code after page 4 of the note): the 13-decimal value
  1.5678672030558 is certified by a chain of $2^{43}+1$ blocks; the 16-decimal 1.5678672030558807 needs a longer chain
  (the 15-decimal 1.567867203055880 passes at $2^{50}+1$).

```sh
./run_all.sh      # Python 3.9+ and node 18+; under a minute (the referee's program is in referee/)
```

## Authorship and the role of AI

This is the work of Arhan Canli (Quantavane). The search, the two programs, the runs and this README were written by AI
agents (Claude, by Anthropic) working for him. The construction we extend is Robert F. Tilley Jr.'s, and the set $X_2$
is Dennis Koshta's. Checked by two independent programs and a fresh referee; not peer reviewed. Questions and
corrections: arhan@quantavane.com.
