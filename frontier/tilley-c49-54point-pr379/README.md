# 49a: Robert F. Tilley Jr.'s 54-point block (PR #379), checked

**The result is Robert F. Tilley Jr.'s (G-14, Inc.).** In
[PR #379](https://github.com/teorth/optimizationproblems/pull/379) of the optimization constants repository (opened
20 November 2026, 18:57 UTC) and the note *A 54-point construction for the Erdős–Szemerédi three-sunflower-free
capacity* ([Zenodo 10.5281/zenodo.24405229](https://doi.org/10.5281/zenodo.24405229)), he proves

$$\mu_3^{S} = C_{49} \ge \left(7211663100+\sqrt{50588282321810251776}\right)^{1/52} = 1.56786638820497592522\ldots$$

by placing two more families, $P \times L(A_p)$ and $L(A_q) \times P$ (6,361,800 sets each), across the join of two
copies of a 28-point block $B$ (Koshta's block without its 21 extra sets). He builds on the gluing lemma of Deuber,
Erdős, Gunderson, Kostochka and Meyer, IRSI's transfer count, the [Cla2026] seed families, Sylla's pendants and Dennis
Koshta's two-pendant block. He states that Gemini 4 Carbon contributed to the discovery and the original Lean
formalization and that OpenAI Codex repaired two proofs, formalized the exact bound and helped prepare the submission.

We read the note (page 2 gives the construction) as mathematics. His `verify.py`, `finite_depth_certificate.py` and
Lean files were neither read nor run.

## Verdict: CONFIRMED (checked by two independent programs; not peer reviewed)

Verdict time 21 November 2026, 15:30 UTC, about 21 hours after the pull request.

| Check | Program A | Program B |
| --- | --- | --- |
| $B$: 252,271 members; port counts $(p,q)$ = 33,948 / 64,728 / 64,728 / 88,867 | yes | yes |
| $\lvert P\rvert = 150$, $\lvert L(A_p)\rvert = \lvert L(A_q)\rvert = 42{,}412$; $L(A_p)$ avoids $p$, $L(A_q)$ avoids $q$ | yes | yes |
| the four rows of the table on page 2, and the totals 4,397,311,488 / 7,211,663,100 / 7,211,663,100 / 11,504,366,352 (30,325,004,040 sets) | yes | yes |
| the whole 54-point family is Sperner and three-sunflower-free (exhaustive, see Method) | yes | yes |
| $bc = 50588282321810251776$, $\rho = a + \sqrt{bc}$ | yes | yes |
| $\rho^{1/52} \in [1.56786638820497592522, \ldots523)$ | yes (20 digits) | yes (21 digits: $1.567866388204975925221\ldots$) |
| the recorded decimal 1.567866388204975 is a valid lower bound (a truncation) | yes | yes |

**Method.** The 54 points are two halves (the left copy without $q$, the right copy without $p$), and each half has
three blocks (core, first pendant, second pendant). Every family here is a union of products over these blocks. Three
sets form a 3-sunflower exactly when on every block their parts satisfy $x\cap y = x\cap z = y\cap z$ and every pair of
the three sets is separated on some block; so for each triple of products it is enough to know which separation
patterns occur on each block, and to combine them. Program A does this with the primitives of our 4 November 49a
checker (`../../c49/sunflower_check.py`, unchanged: the families $F$ and $H$, pattern sets by per-point bitsets, exact
box splitting) plus a new two-level combination (blocks within a half, then the two halves). Program B rebuilds
everything from page 2 on its own (stars, $K_{2,3}$'s and chair complements of $K_5$; $H$ on $\mathbb{Z}_{11}$; the
gluing), computes each pendant's patterns by brute force over all member triples, never merges cells, and does the
rate in BigInt. Both are exhaustive over every multiset of three top-level products (and every ordered pair, for
containment). This checks the 54-point block itself, including the join $G$, without using the gluing lemma.

**Red controls**, all refused. A: a planted set $(\emptyset, c)$ that makes a 3-sunflower with two members of
$P \times L(A_p)$; a planted set strictly inside a member of $P \times L(A_p)$; one member of $A_1$ left out (the
counts then differ). B: the same kind of planted sunflower; one chair complement left out of $F$; $A_q$ built with
the swap $1 \leftrightarrow 2$ instead of $2 \leftrightarrow 3$ (it then uses the deleted port $q$ and makes
3-sunflowers).

**Also found (information):** with Pinho's whole block of PR #370 or Koshta's whole block of PR #260 in place of $B$,
the same join with $P$ is not valid (program A, `--block pinho` or `--block koshta`, reports 3-sunflowers between $G$
and $P \times L(A_p)$).

**Not checked by computation:** the iteration (the gluing lemma and the transfer count on page 3), which turns a valid
54-point block into the bound $\rho^{1/52}$; our `../../c49` tests the gluing lemma on all small instances. The Lean
proof and the finite-depth integer certificate (page 4) were not run.

## Check it yourself

```sh
./run_all.sh      # or, from the repository root: ./check.sh frontier/tilley-c49-54point-pr379
```

No download is needed (the construction is written into both programs); `fetch.sh` fetches the note and the proof
archive from Zenodo for reference and checks their SHA-256. Python 3.9+ and node 18+; a few seconds. Outputs are in
`results/`.

## Authorship and the role of AI

This check is the work of Arhan Canli (Quantavane). The programs, the runs and this README were written by AI agents
(Claude, by Anthropic) working for him. The two programs were written separately, in different languages, and share
no code (program A imports our own earlier checker, which program B does not use). Nothing is claimed beyond what they
check; the construction and the bound are Robert F. Tilley Jr.'s. Questions and corrections: arhan@quantavane.com.
