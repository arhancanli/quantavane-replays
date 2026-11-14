# 88a: the admissible tuples behind the upper-bound rows (14 November 2026)

Page: [88a, the bounded prime gap constant](https://teorth.github.io/optimizationproblems/constants/88a.html)
(main at `87b4ebb`). The page says: "Each upper bound is the diameter of an explicit admissible $k$-tuple, for a $k$
that the sieve analysis determines." Each row is therefore two statements: a sieve theorem, DHL[$k$, 2] (every
admissible $k$-tuple has infinitely many translates with at least two primes), and a finite one, that an admissible
$k$-tuple of the stated diameter exists. **This folder checks the finite statement only.**

## What was checked

For every row, with two separately written programs:

- `c88a_tuples.py` (A): the entries are distinct, there are exactly $k$, the diameter is as stated, and for every prime
  $p \le k$ the set of occupied residue classes mod $p$ is not everything (for $p > k$ nothing is needed).
- `c88a_crt.py` (B): the occupied classes as one bitmask per prime; a free class $a_p$ for each $p \le k$; an integer $n$
  with $n \equiv -a_p \pmod p$ by the Chinese remainder theorem; then $\gcd(n + h, \prod_{p \le k} p) = 1$ for every
  entry $h$, which says directly that no prime $p \le k$ divides any $n + h$.

| Row | $k$ | Tuple (source, SHA-256) | Result (both programs) |
| --- | --- | --- | --- |
| [Zha14] 70,000,000 | 3,500,000 | the $k_0$ primes that follow $k_0$ (the construction of Zhang's Theorem 1, as it is commonly stated; the paper itself was not downloaded) | $\pi(3{,}500{,}000) = 250150$, $\pi(7 \cdot 10^7) = 4118064$, so $\pi(7 \cdot 10^7) - \pi(k_0) = 3867914 > k_0$; the tuple runs from 3500017 to 63374611, diameter exactly 59,874,594 $< 7 \cdot 10^7$; every entry is a prime above $k_0$, so class 0 is free modulo every prime $p \le k_0$ (B strikes the multiples of all 250,150 such primes across the interval and finds no entry struck) |
| [Pol14a] 4680 | 632 | `admissible_632_4680.txt`, Sutherland's tables, `57beae96...7058` | 632 entries, diameter 4680, admissible (115 primes examined) |
| [May15] 600 | 105 | `admissible_105_600.txt`, Sutherland's tables, `d0cdbf5b...2277` | 105 entries, diameter 600, admissible (27 primes) |
| [Pol14b] 246 | 50 | `admissible_50_246.txt`, Sutherland's tables, `347bba78...6976` | 50 entries, diameter 246, admissible (15 primes) |
| [Zha26] 240 | 49 | the field `admissible_tuple.tuple` of `data/certificate.json`, hanxinzhang/prime-bound @ `2319570b`, `d99e5d7c...b6f5` | 49 entries, diameter 240, admissible (15 primes) |
| [Sta26] 240 | 49 | `admissible_49_240.txt`, Sutherland's tables, `c17f9b2c...a7f1` | 49 entries, diameter 240, admissible |
| [OAI26] 186 | 40 | the literal `admissibleTuple` of `Challenge.lean`, openai/PrimeGaps186 @ `61340d0b`, `b99a403b...e1f0` (read as text; it begins 0, 2, 6 and ends 182, 186, as the row says) | 40 entries, diameter 186, admissible (12 primes) |

Sutherland's 40-tuple of diameter 186 (`admissible_40_186.txt`, `d827864e...f50a`), a different tuple from the one
in [OAI26], is admissible as well (information only).

Sutherland's tables (math.mit.edu/~primegaps, linked from the page) are a web site without commits; the hashes pin
the files as retrieved on 14 November 2026. Full hashes are in `SHA256SUMS`.

## Red controls (all refused by both programs)

- [Zha14]: the first 3,500,000 primes instead of those after $k_0$: every class is occupied modulo each prime up to 47.
- For each of the six tuple rows: one entry replaced by an odd number (both classes mod 2 occupied), and one entry
  replaced by an even number that fills every class modulo an odd prime (3, or 11 for the 632-tuple).
- The 50-tuple offered as a 51-tuple, and with diameter 244.

`results/tuples.txt` holds the full output: 8 tuples accepted by both programs, 15 controls refused by both.

## What was NOT checked

- **No sieve theorem.** DHL[3500000, 2] (Zhang), DHL[632, 2] (Polymath8a), DHL[105, 2] (Maynard), DHL[50, 2]
  (Polymath8b), DHL[49, 2] ([Zha26], [Sta26]) and DHL[40, 2] ([OAI26]) are taken from the sources. In particular
  nothing here bears on the numerical certificates or the Lean development of [OAI26], or on the exact certificate
  of [Zha26] (only its tuple was read). The three 2026 rows rest on preprints.
- The values of $k$ are those of the sources as the page and the linked repositories report them; we did not derive
  them.
- That these diameters are the least possible for their $k$ (not needed for an upper bound).

## Reproduce

`./fetch.sh && ./run_all.sh` (Python 3.9+; program B's [Zha14] mode needs numpy; about 15 seconds, under 1 GB).
No downloaded file is executed: the Lean file is read as text, the JSON with the `json` module.
