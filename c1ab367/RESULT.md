# 1a and 1b: the step functions of PR #367, of [Cum2026] and of PR #368 (14 November 2026)

Four certificates, each replayed by two separately written exact programs that agree on every number:

- `c1ab367_kron.py` (A): Python integers. The whole sequence ($c_k$ for 1a, $Q P_m - R_m$ for 1b) from one
  big-integer product by Kronecker substitution (base $2^b$ with $b$ larger than any coefficient needs, so no digit
  carries); every value at the maximum and 40 others at random summed again term by term.
- `c1ab367_direct.py` with `c1ab367_direct.c` (B): its own reader, then plain double loops in C over every lag, in
  unsigned 128-bit integers (heights of more than 56 bits split into equal limbs, two to four, with every size bound
  checked before the loop; the partial sums recombined in 384 bits). No product of big integers.

Both reduce decimal heights to the same primitive integers (scaled to integers, then divided by their gcd; the ratio
does not depend on the scale), and both print the fingerprint $\sum_t v_t 3^t \bmod 2^{61}-1$ of the whole sequence:
the fingerprints agree for all four certificates. No program of any source was read or run; the READMEs of
intervalplace/frontier and the README of Cumes's ancillary files were read as documentation.

| row | certificate (pinned) | $n$ | exact value | maximum | printed |
|---|---|---|---|---|---|
| 1a [Esc2026], the row in #367's files | intervalplace/frontier @ `9edc297`, `C1_262144_1.5026648332.txt`, `4c257d45…0a1f` (as stated in the row) | 262,144 | 735013322210029338350039415555883008 / 489139897323017668625465543540606041 $= 1.50266483317479630415\ldots$, the stated fraction | $k = 344084$ only, as stated | 1.5026648332, valid ($2.5 \cdot 10^{-11}$ above) |
| 1a, the file #367's description and title name | intervalplace/frontier @ `8c9fa72`, `C1_262144_1.5026648336.txt`, `b33995fe…5495` (as stated in the description) | 262,144 | 11484584234108926456768515550896128 / 7642811608550585358028269349510225 $= 1.50266483361440736448\ldots$, the stated fraction | $k = 316694$ only, as stated | see the note |
| 1a [Cum2026], recorded by #367 | Cumes's ancillary file, served at qwntl.com on 14 Nov (no commit): `49c83eb6.xz` (`c751eeb3…5d13`), unpacked `49c83eb6…2f2d` (its name is its SHA-256) | 262,144 (decimals, some with an exponent) | 10326234156771825140409842478149306017059563804408018744901632000000 / 6871947673599997545676727668735312122788541593125740957532237529929 $= 1.50266484077609912906120468\ldots$; the stated 23 digits | $k = 275085$ only | 1.50266484078, valid ($3.9 \cdot 10^{-12}$ above) |
| 1b [Grad2026], #368 | intervalplace/frontier @ `6d6c4ba`, `C1b_15360_0.3808564712.json`, `7bfe24a3…5b7b` (as stated) | 15,360, heights $k_i/2^{50}$ | 3707849734734953936493216938911489 / 9735556609752801803494680617287680 $= 0.38085647111543024238\ldots$, the stated fraction; $0 \le k_i \le Q$ and $\sum k_i = nQ/2$ exactly | shift $m = 6384$ only, as stated | 0.3808564712, valid ($8.5 \cdot 10^{-11}$ above) |

**Note (information, not a red control).** The title and description of PR #367 still give $C_{1a} \le 1.5026648336$
for the certificate at `8c9fa72`, whose exact value $1.50266483361440\ldots$ is $1.44 \cdot 10^{-11}$ above
1.5026648336; both programs record this. The row and README cell in the pull request's files were changed to the
later certificate (1.5026648332), which holds as printed.

**Red controls**, each refused by both programs: #367's row with one height plus one unit (a different fraction),
claimed as $\le 1.5026648331$, and with one height made negative; [Cum2026] with one height plus $10^{15}$ units
(the stated 23 digits fail) and claimed as $\le 1.50266484077$; #368 with a height equal to $Q$ plus one (above $Q$,
and the mass is off), with $1000/Q$ moved between two cells that meet at the maximizing shift (mass unchanged, a
different fraction), and claimed as $\le 0.3808564711$.

**Not checked:** that these step functions are the limits of the authors' searches (irrelevant to the bounds), and the
reduction from a step function to the bound, which is the page's definition ($f*f$ linear between knots for 1a, the
overlap linear between grid shifts for 1b).

Reproduce: `./fetch.sh && ./run_all.sh` (Python 3.9+ and a C compiler; about five minutes on two cores, under
300 MB). Output: `results/c1ab367.txt`.
