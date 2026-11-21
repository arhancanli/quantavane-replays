# Thirteen binary constant-weight codes (P. R. A. Iskandar, arXiv 2611.24958)

Checked by two independent programs; not peer reviewed.

## The claim (credit first)

Peter Rafat Adel Iskandar, "Thirteen Improved Binary Constant-Weight Codes with Explicit Cyclic Orbit Certificates",
[arXiv 2611.24958](https://arxiv.org/abs/2611.24958) (math.CO). It was submitted on 10 November 2026 and announced on
20 November 2026, and the codes were sent to A. E. Brouwer's table on 9 November (label Ad). The paper gives thirteen
codes: A(16,4,7) >= 841, A(20,4,5) >= 876, A(25,4,6) >= 7800, A(29,4,5) >= 4491, A(29,6,5) >= 333,
A(36,4,5) >= 11207, A(38,4,5) >= 13737, A(39,4,5) >= 15618, A(44,4,5) >= 25585, A(54,4,5) >= 60632,
A(55,4,5) >= 66144, A(60,4,5) >= 94105 and A(62,4,5) >= 107360. Together they have 408,729 words, generated from
11,386 orbit representatives under a permutation g that shifts consecutive coordinate blocks cyclically. The author
writes: "Anthropic Claude Code and OpenAI ChatGPT assisted with computational experimentation, certificate software,
structural analysis, literature checking, and manuscript preparation under the author's direction."

## One command

```sh
./fetch.sh && ./run_all.sh
```

or, from the repository root, `./check.sh frontier/iskandar-cw13-2611.24958` (from this folder:
`../../check.sh frontier/iskandar-cw13-2611.24958`). `fetch.sh` downloads `CLAIMS.csv`, `ORBIT_MANIFEST.csv`, the 13 word
lists and the 13 representative lists from the paper's ancillary files into `certs/`. It makes one request every
15 s, so it takes about 7 minutes, and it checks every SHA-256 against `SHA256SUMS`. The paper's three `.py` verifiers
are not downloaded. `run_all.sh` needs Python 3.9+ and a C compiler and takes under a minute. Output:
`results/cw13.txt`.

## Verdict: confirmed

All thirteen codes are valid. Each has the right length, every word has weight w, no word repeats, and the size is
exactly as claimed. Two words share at most w - d/2 coordinates, so the minimum distance is at least d; it is exactly d
for every code. Each list also equals the union of the orbits of its representatives under g, and the orbit sizes match
the manifest.

| code | largest intersection of two words (program B) | orbit sizes (both) |
|---|---|---|
| A(16,4,7) >= 841 | 5 | 1 x 1, 120 x 7 |
| A(20,4,5) >= 876 | 3 | 3 x 1, 291 x 3 |
| A(25,4,6) >= 7800 | 4 | 1560 x 5 |
| A(29,4,5) >= 4491 | 3 | 9 x 1, 1494 x 3 |
| A(29,6,5) >= 333 | 2 | 111 x 3 |
| A(36,4,5) >= 11207 | 3 | 1 x 7, 320 x 35 |
| A(38,4,5) >= 13737 | 3 | 723 x 19 |
| A(39,4,5) >= 15618 | 3 | 411 x 38 |
| A(44,4,5) >= 25585 | 3 | 595 x 43 |
| A(54,4,5) >= 60632 | 3 | 1144 x 53 |
| A(55,4,5) >= 66144 | 3 | 1248 x 53 |
| A(60,4,5) >= 94105 | 3 | 1595 x 59 |
| A(62,4,5) >= 107360 | 3 | 1760 x 61 |

Program A finds no t-subset (t = w - d/2 + 1) in two words of any code. Totals for both programs are 408,729 words and
11,386 representatives.

**The table.** We compared the sizes with two copies of Brouwer's table that we had saved. One is dated 3 November 2026
(Last-Modified header) and one 16 November 2026; we have no copy from 9 November. All thirteen sizes exceed the
3 November entries. For twelve cells, the paper's "previous" value equals the 3 November entry. For A(29,4,5) the paper
gives 4466 and the 3 November table 4423, and we cannot see the value on 9 November. Six of the improved entries were
ours (label Can). On 16 November the table lists ten of the codes with label Ad. For three cells it lists larger codes
with label Pan: A(29,6,5) 336 (main grid), A(36,4,5) 11247 and A(55,4,5) 66591. Each of the paper's thirteen sizes was
an improvement on the entry it cites.

## The two programs (written separately, no shared code)

- **A** = `cw_a.py` (Python). It checks length, alphabet, weight, repeats and size. For the distance, every t-subset of
  every word goes into one hash set; a t-subset met twice would mean two words at distance below d. It rebuilds each
  code from its representatives under g, as sets, and checks the SHA-256 of each list against `CLAIMS.csv`.
- **B** = `cw_b.c` (C, 64-bit bitsets). It checks popcount weights, sorts the words to find repeats, and checks the
  size. It computes the intersection |x & y| of every two words that share at least two coordinates, inside the bucket
  of a shared pair of coordinates. A pair at distance below d shares at least 3 coordinates, so none is missed. For the
  nine codes with at most 30,000 words it also intersects every pair directly, and the two maxima agree. It rebuilds the
  orbits by bit permutation and compares the sorted lists.

## Red controls (both programs refuse each)

In A(16,4,7), word 2 is replaced by:
1. word 1 with one 1 moved (distance 2 from word 1);
2. a copy of word 1;
3. word 1 with an extra 1 (weight 8).

## Not checked

The upper bounds U in the paper's Table 1 (Johnson recursion), the automorphism groups, and the search that found
the codes.

## Authorship and the role of AI

This check is the work of Arhan Canli (Quantavane). The programs, runs and this README were written by AI agents
(Claude, by Anthropic) working for him. The two programs were written separately. Nothing is claimed beyond what they
check, and the mathematics is the claimant's. Questions and corrections: arhan@quantavane.com.
