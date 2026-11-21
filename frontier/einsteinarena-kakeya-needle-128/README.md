# Discretized Kakeya needle, n = 128: the EinsteinArena submissions of Witness and Zlab-QiHuang, checked

**The constructions are the submitters'.** On [EinsteinArena](https://einsteinarena.com/problems/kakeya-needle-128),
where AI agents compete on open problems, the problem *Discretized Kakeya Needle (n = 128)* (Problem 6.9 of
[Mathematical exploration and discovery at scale](https://arxiv.org/abs/2511.02864); Section 4.4 of
[The Station](https://arxiv.org/abs/2608.23691)) asks for offsets $x_1, \ldots, x_{128}$ minimising the area of
$\bigcup_j T_j$, $T_j = \mathrm{conv}\{(x_j,0), (x_j + \tfrac{1}{128}, 0), (x_j + \tfrac{j}{128}, 1)\}$.

| Submission | Agent | Posted (UTC) | Stated score |
| --- | --- | --- | --- |
| 2816 | **Witness** | 19 Nov 2026, 15:48 | 0.10673456419616238 (first on the board) |
| 2849 | **Zlab-QiHuang** | 19 Nov 2026, 23:04 | 0.10674472375252743 (second) |
| 2772 | kgrey (previous best, for comparison) | 10 Nov 2026, 22:41 | 0.10675138995168142 |
| 2563 | the Station's published construction (the seed) | 8 Sep 2026 | 0.10706663656163481 |

The offsets are copied in `data/solutions.json` from the arena's public API (`/api/solutions/best?problem_id=26`,
read 21 November 2026, 15:33 UTC), with credit to the submitters. We read the problem statement, which includes the
arena's verifier source and describes its method (exact horizontal slices); program A follows that published method
in its own code, program B uses a different one. The arena's verifier was not run.

## Verdict: CONFIRMED for both (checked by two independent programs; not peer reviewed)

Verdict time 21 November 2026, 15:35 UTC: 47.8 hours after Witness's submission, 40.5 hours after Zlab-QiHuang's.

Both programs give the same exact rational area for each of the four submissions, and its nearest double is the
stated score:

- Witness: $\dfrac{120494307936022552266703253242502036091311510940726864328739}{1128915537750000000000000000000000000000000000000000000000000} = 0.106734564196162382270039983106\ldots$
- Zlab-QiHuang: $\dfrac{401685924056865660335698363413073200423137963891728095077}{3763051792500000000000000000000000000000000000000000000000} = 0.106744723752527426935673345084\ldots$
- kgrey $= 0.10675138995168142\ldots$ and the Station seed $= 0.10706663656163481\ldots$ as stated.

So Witness's construction lowers the best area on the board from kgrey's 0.1067513899... by $1.68 \times 10^{-5}$, and
Zlab-QiHuang's lies between the two.

- **Program A** (`kakeya_a.py`, Python integers and fractions): every height in $(0,1)$ where two of the 256 section
  endpoints $x_j + jy/128$ and $x_j + 1/128 + (j-1)y/128$ cross (21,187 heights for Witness's construction); the union
  length is affine between them, so the area is the exact trapezoid sum.
- **Program B** (`kakeya_b.mjs`, JavaScript BigInt rationals): Green's theorem. The area is the integral of $x\,dy$
  around the union's boundary; for each side of each triangle the heights where it lies in the interior of another
  triangle form an open interval per triangle, and the rest of the side is integrated exactly (479 boundary pieces
  for Witness's construction). Coinciding sides of neighbouring triangles have opposite directions and cancel.

**Red control**, refused by both: the 64th triangle of Witness's construction moved by 1/1000 gives a different area
(0.1068044479864619...), the same in both programs. A translation of all triangles leaves the area unchanged in both.

**Not checked:** whether these constructions were found by the agents or taken from elsewhere (the arena asks
readers to check); lower bounds for the problem.

## Check it yourself

```sh
./run_all.sh      # or, from the repository root: ./check.sh frontier/einsteinarena-kakeya-needle-128
```

Python 3 and node 18+; under a minute. `run_all.sh` checks `data/solutions.json` against `SHA256SUMS`, runs both
programs on all four submissions, compares the exact fractions, and runs the red control. Outputs are in `results/`.

## Authorship and the role of AI

This check is the work of Arhan Canli (Quantavane). The programs, the runs and this README were written by AI agents
(Claude, by Anthropic) working for him. The two programs were written separately, in different languages, with
different methods, and share no code. Nothing is claimed beyond what they check; the constructions are the
submitters'. Questions and corrections: arhan@quantavane.com.
