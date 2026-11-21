# Erdős problem 985: a prime primitive root below every odd prime, the finite part (Howe, Huffman, Kohlmeyer, Kruer)

Checked by two independent programs; not peer reviewed.

## The claim (credit first)

Will Howe, Grant Huffman, Jensen Kohlmeyer and Liam Kruer, "A prime primitive root below every odd prime",
[H0W4RD-dev/erdos-985-prime-primitive-root](https://github.com/H0W4RD-dev/erdos-985-prime-primitive-root), first
commit 21 November 2026 00:05 UTC (checked at commit `99fd0a7`, 21 November 00:16 UTC). They prove that for every odd
prime p there is a prime q < p of multiplicative order p - 1 modulo p (Erdős problem 985). Above 2791 the proof uses
the theorem of McGown, Treviño and Trudgian under GRH together with Alpöge's announced proof of GRH; Lemma 2.1 settles
the 405 odd primes p <= 2791 by computation (`small-prime-witnesses.tsv`). A companion argument (`endpoint.md`) leaves
74,187 candidates, 18,802 of them prime moduli, certified in `x-equals-p-witnesses.zip`. The authors write: "The authors
acknowledge assistance from OpenAI Codex in preparing the manuscript and verifying the finite calculations."

## One command

```sh
./fetch.sh && ./run_all.sh
```

or, from the repository root, `./check.sh frontier/erdos985-prime-primitive-root` (from this folder:
`../../check.sh frontier/erdos985-prime-primitive-root`). `fetch.sh` downloads the two data files from the pinned commit
into `certs/`, unpacks the seven certificate strips and checks every SHA-256 against `SHA256SUMS`. Python 3.9+ and a C
compiler; about one minute. Output: `results/erdos985.txt`. The repository's `verify.py` was neither read nor run.

## Verdict: the finite statements are confirmed

- **(a) Lemma 2.1.** Every one of the 405 odd primes p <= 2791 has a prime primitive root below p. The table has a row
  for exactly these primes, and every witness q in it is prime, below p and of order p - 1; in all 405 rows it is the
  least such prime. The largest is 53, at p = 2791.
- **(b) Companion certificate.** The candidate moduli of the seven columns (odd n, 2792 <= n < u_k^2, p - 1 with
  exactly k = 3, ..., 9 distinct prime factors) number 4459, 17520, 26905, 19435, 5233, 625 and 10. That is 74,187 in
  all, equal to the number of words in each strip. Among them 944, 4146, 6835, 5215, 1462, 195 and 5 are prime, 18,802
  in all. Every one of these primes has a prime primitive root below it: the certificate's own witness (always correct)
  and our own least one, the largest being 173 at n = 16,456,441. No composite candidate has a primitive-root word.
  This is a direct check of the conclusion for those moduli; the character-sum argument that reduces the problem to
  them is not checked.
- **(c) Extension.** Every odd prime p <= 10^7 (664,578 primes) has a prime primitive root below p. The least ones sum
  to 3,911,151, and the largest is 211, first at p = 4,022,911. Program B alone also ran to 10^8 (5,761,454 primes; sum
  34,008,187; largest 223, first at p = 73,189,117) with `cc -O2 erdos985_b.c -o b && ./b certs/erdos985/small-prime-witnesses.tsv certs/erdos985/witnesses 100000000`.

Both programs print the same lines for (a), (b) and (c); `run_all.sh` compares them line by line.

## The two programs (written separately, no shared code)

- **A** = `erdos985_a.py` (Python, exact integers). A sieve with a smallest-prime-factor table factors p - 1. A prime q
  is a primitive root when q^((p-1)/r) is not 1 mod p for every prime r dividing p - 1. The certificate's candidates come
  from a depth-first search over the factorisations of n - 1, and n is tested by trial division. A also checks the
  table's factorisation and residue columns, and the certificate's compositeness words (a proper divisor, or a Fermat or
  Euler witness).
- **B** = `erdos985_b.c` (C, 64-bit integers with 128-bit products). For p <= 2791 it finds orders by walking the powers
  q, q^2, ... until 1, with no factorisation. For the certificate it runs one segmented sieve over n < 20443^2, which
  gives the number of distinct prime factors of n - 1 by dividing out every small prime, and the primality of n by
  Eratosthenes on the same segment. It matches candidates to words in increasing order. For the extension it uses a bit
  sieve and trial division.

## Red controls (both programs refuse each)

1. Table row p = 7 given q = 2, which has order 3 mod 7.
2. The first prime modulus of the rank-3 strip, n = 2797, has its witness 2 replaced by the prime 3, which is not a
   primitive root mod 2797.

## Not checked

The main theorem above 2791. It rests on McGown, Treviño and Trudgian's theorem under GRH and on Alpöge's announced proof
of GRH, and we did not check either. Also not checked: the analytic side of the companion argument (their Lemma 1 bounds,
c(p,p) < 50077/56784, inequality (1) and the column thresholds u and L).

## Authorship and the role of AI

This check is the work of Arhan Canli (Quantavane). The programs, runs and this README were written by AI agents
(Claude, by Anthropic) working for him. The two programs were written separately. Nothing is claimed beyond what they
check, and the mathematics is the claimants'. Questions and corrections: arhan@quantavane.com.
