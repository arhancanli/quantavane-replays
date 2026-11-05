# Quantavane replays

Independent exact replays of certificates submitted to the
[optimization constants repository](https://github.com/teorth/optimizationproblems). Each replay uses our own code,
written from the definition on the constant page. The certificates' own verifiers were neither read nor run, and we
have no connection with any of the authors.

An index of every asterisked bound in that repository, with the replays posted for each (ours and others'), is in
[ASTERISK_INDEX.md](ASTERISK_INDEX.md) (snapshot of 5 November 2026).

## C1c and C1d (3 to 5 November 2026)

| Constant | Pull request | Certificate (pinned commit, SHA-256) | Exact value | Printed bound |
| --- | --- | --- | --- | --- |
| $C_{1c} \geq$ | [#292](https://github.com/teorth/optimizationproblems/pull/292) (Koshta's row) | DennisKoshta/c1c-autocorrelation-certificate @ `a4ee94ef29ff`, `K24576_0.4103333.json`, `963cb727...c9e` | 134648548624804778248790016 / 328144294688458468928334649 = 0.41033335274849333953..., minimum at k = 21264 only | 0.4103333, below the exact value |
| $C_{1c} \geq$ | [#318](https://github.com/teorth/optimizationproblems/pull/318) | Julsgaard/c1c-certificate @ `91e47b1a`, `cert/certificate.json` (K = 98,304, 98,368 weights), `fe11336b...52d2` | 2582127846607416062611587072/6292441470965870058830439361 = 0.4103538918115781099..., minimum only at k = 47318: the printed 0.4103538918 is a valid rounding down |
| $C_{1c} \geq$ | [#270](https://github.com/teorth/optimizationproblems/pull/270) (update of 4 Nov) | vidigoat/optimization-records @ `9f82508a`, `tao/certificates/c1c_v2/c1c_K98304_v2.json` (K = 98,304), `3bba62b0...4729` (no hash stated; this is the file at that commit) | 645458795582745374892957696/1572931393892132517694378321 = 0.4103540676275733581..., minimum only at k = 93282: the printed 0.4103540676 is a valid rounding down |
| $C_{1c} \geq$ | [#270](https://github.com/teorth/optimizationproblems/pull/270) (update of 5 Nov) | vidigoat/optimization-records @ `42d6a89f`, `tao/certificates/c1c_v3/c1c_K98304_v3.json` (K = 98,304, 98,368 weights), `e0786a59...ec66` | 640823574552569898220412928/1561614482809354047911104009 = 0.4103596512499835112..., minimum only at k = 98301: the printed 0.4103596512 is a valid rounding down |
| $C_{1d} \leq$ | [#312](https://github.com/teorth/optimizationproblems/pull/312) | Julsgaard/c1d-certificate @ `cc455ca0`, `cert/certificate.json` (81,920 cells), `04c9c759...237` | 9900630255081376929482866688 / 6805830738283663732311831361 = 1.45472766452875064884..., maximum only at i + j = 29124 | 1.454727664528751, above the exact value |
| $C_{1d} \leq$ | [#271](https://github.com/teorth/optimizationproblems/pull/271) | Keremozdemirra/tao-1c-1d-certificates @ `408f8d0b`, `c1d_n40960_1.4556234448.json` (40,960 cells), `be35060a...2da` | 279676466330630562685427712000 / 192135175720687825901095191121 = 1.45562344470022455281..., maximum only at i + j = 15313 | 1.4556234448, above |
| $C_{1d} \leq$ | [#269](https://github.com/teorth/optimizationproblems/pull/269) | echul-seal/c1d-upper-bound @ `2be42dfd`, `ub/c_f10240.json` from `c1d_records.tgz` (10,240 cells), `30d501e1...973` | 76358168487433178963290071040 / 52455219864513498497744894329 = 1.45568293650581519488..., maximum only at i + j = 3168 | 1.455682936505816, above |
| $C_{1d} \leq$ | [#286](https://github.com/teorth/optimizationproblems/pull/286) | huiihao/c1d-certificate @ `302817dc`, `construction.json` (400 cells, JSON floats), `c489b46b...8df` | read as binary64: 4563650116209721766004955898613711667200 / 3127595284199802099438492053245424408273 = 1.45915622115964902113...; read as the printed decimals: 1.45915622115964903535... | 1.4591562211596491, above both |

Every hash above is the one stated in the pull request (for #269 also the archive's, `b526326d...`), except the 4 Nov #270 row, where none was stated.

## Method

- $C_{1c}$: $f = \sum_i w_i 1_{[i/K,(i+1)/K)}$ with non-negative integer weights. The autocorrelation is linear between
  multiples of $1/K$, so the bound is $K \min_{0 \le k \le K} A_k / S^2$ with $A_k = \sum_i w_i w_{i+k}$ and
  $S = \sum_i w_i$ (we check $n > K$, so no $A_k$ in the range is an empty sum).
- $C_{1d}$: $n$ equal cells on $[-1/4,1/4]$, so $h = 1/(2n)$. $f*f$ is linear between multiples of $h$ and equals
  $h \sum_{i+j=m-1} v_i v_j$ at the knots, so the bound is $2n \max_m \lvert c_m \rvert / S^2$ with $c = v * v$.
- All sums are exact integers: every $A_k$ (or $c_m$) comes from one big-integer product (Kronecker substitution;
  positive and negative parts packed separately for signed values). The whole sequence is then confirmed by a
  polynomial identity at four random points modulo $2^{61}-1$, and the deciding lag is summed again directly.
- Red controls, each run on every certificate: one corrupted coefficient is rejected by the identity; the old
  coefficients are rejected after one weight changes by 1; scaling every value by 3 leaves the bound unchanged.

## C49: the 3-sunflower-free capacity (4 November 2026)

`c49/sunflower_check.py`, our own code: every base family is checked Sperner and 3-sunflower-free exhaustively (a pivot
method and a pair method, plus a box method for the two 250,000-set families; all agree with brute force on 3,400
random and glued test families), the port counts and transfer matrix are recomputed, and the rate is enclosed by exact
integer comparisons. Every printed bound is the 7-decimal truncation of the exact value.

| Row | Certificate (SHA-256) | Exact value |
| --- | --- | --- |
| IRSI v1.0, 1.5617847 | `cert_c49_chain20.json` (`c30a74fe...715b`), 6,127 sets on 20 points | rho^(1/18) = 1.56178470895053727967... |
| IRSI v1.1, 1.5621220 | `cert_c49_chain20_s17.json` (`c12d0dd6...fc7b`), 6,150 sets | 1.56212205600681534404... |
| IRSI v1.2, 1.5647730 | `cert_c49_chain20_v12.json` (`7472112c...4e82e`), 6,341 sets | 1.56477309608283129921... |
| [Cla2026] F, 1.5664286 | F rebuilt from the page, equal to `k5_family_masks.txt` (`34cab736...7659`), 75 sets | rho = 20 + 2 sqrt(66), rho^(1/8) = 1.56642866225003666041... |
| [Cla2026] F and H, 1.5666969 | H rebuilt from the page, equal to `z11_family_masks.txt` (`09c02bec...a807`), 121 sets | rho^(1/17) = 1.56669697875431533092... |
| [Syl2026] unit chain, 1.5675517 | `F.txt`, `H.txt`, `B.txt` (hashes as in their SHA256SUMS); gadget G (250,659 sets on 28 points) rebuilt from the page | rho^(1/60) = 1.56755175575818821977... |
| [Syl2026] B', 1.5678015 | `Bprime.txt` (`2a723857...5cec`), 4,364 sets | rho = 1127 + sqrt(925536), rho^(1/17) = 1.56780152016874834047... |
| [Kos2026], 1.5678599 | `two_pendant_extra_1.5678599.txt` (`10ee18f9...4838`), 252,292 sets on 28 points | rho = 64728 + sqrt(3018723123), rho^(1/26) = 1.56785990242458988585... |

What rests on mathematics rather than computation: the gluing lemma (DEGKM 1997, Lemma 4.4, as stated in the
claimants' READMEs; we proved it by hand and tested it on 80,889 small gluings and on explicit chains up to 155,892
sets), the chain argument, and the page's limit step. The 62-point unit of the 1.5675517 row (about 1.09e12 members)
is valid by that lemma. Two notes on the page: the [Kos2026] row says 15 sets are added at each pendant, but the
certificate (correctly) adds 806 per pendant, 11 x 46 + 4 x 75, as #260 says; and the IRSI rows print matrices in the
IRSI README's convention, which differs from the later rows' (conjugate matrices, the same rho).

Reproduce: `cd c49 && ./fetch.sh ROOT && ./run_all.sh ROOT` (about 9 minutes on a loaded 18-core Mac).

## C86: the Schur-Siegel-Smyth trace constant (4 November 2026)

`c86/c86_check.py`, our own code, standard library only: the six finite conditions of Smith's criterion as the 86a
page states them. E1 (density, support, mass m), E3 and E6 in exact rationals; E2 (primitive, irreducible over Q by a
mod-p factor-degree sieve and by an exact root-subset test, totally real by exact Sturm counts); E4 (integral of
log|Q| d mu) from rigorous logarithms at scale 2^-128 and isolated roots; E5 (U_mu(x) - sum b_Q log|Q(x)| >= c on the
support) by per-cell Taylor models with explicit remainder bounds, exact handling of each cell's edge cusps and nearby
roots, and adaptive bisection. Floats only choose split and tangent points; they never decide anything.

#298 needs more: its density is piecewise linear on 5,375 unequally spaced rational nodes, and ten of its 28
auxiliaries have nonreal roots. There U(z) = Re sum_k j_k F(z - e_k), with j_k the slope jumps and
F(w) = w^2 log(w)/2 - 3w^2/4 (checked against brute-force quadrature at 22 real and 20 complex points). Every root is
enclosed in a Gershgorin disc of the Weierstrass matrix (the code of our C82 replay, radius at most 2^-305), the real
ones also by exact Sturm isolation, and the enclosures satisfy Vieta's formulas. E4 uses our own arctangent at the
nonreal roots; E5 a degree-24 Taylor model per support cell, with the nodes and roots within four half-widths of the
cell handled exactly on each piece.

| Row | Certificate (SHA-256) | Result |
| --- | --- | --- |
| [Kos2026], 1.8084990 | `certificate/c86_A19d7_robust_1.8084990.json` @90c57242 (`1a9e108e...845c`) | 19 auxiliary polynomials pass E2; E3 ratio 0.83542405; smallest E4 integral 1.99999997804e-6; E5: f >= 1.01159e-6 on all 633 cells (4,403 pieces), min f in [1.011592e-6, 1.011601e-6] near x = 2.869309284; T/m = 144679917840593158601/79999999999999881200 = 1.80849897300741716813... |
| [Pin2026], 1.8117251 | `cert_n1200_A11.json` @116a58f (`722a287b...e53d`) | 11 polynomials pass E2; E3 ratio 0.82977960; smallest E4 9.99999986901e-5; E5: f >= 5.03972e-5 on all 629 cells, min f near x = 0.709342054; T/m = 144938007612277893/79999999999994000 = 1.81172509515360954188... |
| #270, 1.8077686 | `c86_deg8_h400_1.807768536863.json` (`7ca02d70...1134`) | 19 polynomials (two of degree 8) pass E2; E3 ratio 0.84520107; smallest E4 9.98933032188e-11 (thin but rigorous); E5: f >= 4.86106e-7 on all 1,243 cells; T/m = 578485931795944337659/319999999999999303200 = 1.80776853686232999160... |
| #298, 1.80560940 | `c86/adaptive-cone/certificate.json` @0d82162 (`10d68ec2...0dad5`) | 28 polynomials (10 with nonreal roots, degree up to 21) pass E2; E3 ratio 0.8950913258; smallest E4 8.94715228390e-8 (x^2 - 4x + 2); E5: f >= 5.57254e-7 on all 3,791 support cells (118,757 pieces), min f in [5.5725403e-7, 5.5725405e-7] near x = 0.285640644; T/m = 288897505229428473/160000001279369600 = 1.80560939324616691129... |

Each printed bound is T/m rounded up, so each holds as printed. Red controls, all as expected: lowering the densest
cell by 1% makes E5 fail (and E4 for some Q); raising the largest weight by 0.1% makes E5 fail; reducible,
non-primitive and non-totally-real polynomials are refused; weights times 1.3 break E3; a bound 1e-15 low is refused
by E6. For #298 also: the weight of x - 2 raised to twice what breaks E5 at the minimiser makes f < 0 there and the
bisection refuses that cell; the highest node lowered by 1% makes E4 fail for two polynomials; a nonreal root moved by
1e-12 fails the Vieta test. What rests on theorems rather than computation: the minimum principle, Smith's Prop. 5.15
and Thm 1.5 (we check their hypotheses).

Reproduce: `cd c86 && ./fetch.sh ROOT && ./run_all.sh ROOT` (about 3.5 minutes; the #298 row about 1.5 minutes and 5
CPU-minutes on a loaded 18-core Mac).

## C44 and C70 (4 November 2026)

`c44/c44_check.py`: Pinho's C44 <= 4.3576549226 (rpinho/c44-certificate @f3fb25d). The paper's recurrence
U(b,d) = min(c_d, d 2^-d + max_{k<=d} U(min(b-1, B_k), k)), caps c_d = d/2 and c_3 = 5/4, with B_d from the 48
degree certificates (each checked exactly: integer nodes, moment identities, positivity) intersected with Wellens's
3b(b-1) <= 2(d^4-d^2) to degree 96, plus Wellens's tail 233045/2^95, computed exactly three ways: the bound is
172624496193406565421448533795/39614081257132168796771975168 = 4.35765492257445791338..., the page's rational, and
the printed value is it rounded up. As a function of the degree-3 cap c the bound is this value + clamp(c - 173/128,
0, 1/128), so the claim needs only D <= 173/128 at degree 3 (a cap of 11/8, from R_3 <= 11, returns Epi's bound).

`c70/c70_check.py`: Pinho's C70 >= 1.20194336 (rpinho/c70-certificate @680e074). The single Steiner term at
n = 10^12, k = 367879441172 in our own interval arithmetic at 2^-400: R in [1.201943368443150566994802779735524137849309813,
...814] (width 4.4e-117); the printed bound is R truncated (margin 8.44e-9). exp(1/(2e)) is enclosed too
(1.20194336847031446719...).

Reproduce: `cd c44 && ./fetch.sh ROOT && ./run_all.sh ROOT` (about 30 s), the same in `c70` (under a second).

## C82: the essential minimum of the Zhang-Zagier height (4 November 2026)

`c82/c82_check.py`, standard library only: the bound B = [sum_j m_j log M(R_j o T) + int_0^1 max(0, g) dt] / D with
T(z) = z(1 - z). The hypotheses (integer coefficients, monic denominators, deg A < D, coprimality of every pair) are
checked exactly; every root of every composed factor is enclosed in a disc of radius at most 2^-303 by Gershgorin's
theorem on the Weierstrass matrix; the integral is enclosed on an adaptive dyadic tiling with per-root series and
rigorous remainders; the comparison with the printed bound is exact.

| Row | Certificate (SHA-256) | Enclosure of B | Printed upper bound |
| --- | --- | --- | --- |
| [Pin26] | rpinho/c82-height-certificate @69d8964 `certificate.json` (`b33626f3...0fce`) | [0.2534446873318462, 0.2534446873318463] | 0.2534449 |
| [Epi26] | BakeLab/Epi-results @97fa6a4 `construction.json` (`e85e26fa...e559`) | [0.2535689302577047, 0.2535689302831411] | 0.2535744749068227 (valid, 5.5e-6 loose) |
| #280 | vidigoat/optimization-records @32898ac (`557da3e3...5fe0`) | [0.2534350401279890, 0.2534350401279891] | 0.2534350402 |
| #301 | rpinho @6d4e80c `degree44-denominators-v2` (`adfca075...7ac2`) | [0.2530918335514613, 0.2530918335514614] | 0.253091834 |
| #301 first filing | rpinho @6d4e80c `degree44-denominators` (`6a24d45c...f056`) | [0.2531232746206604, 0.2531232746206605] | 0.253123275 |
| control: Doche 2001 | [Pin26] `controls/doche2001.json` (`835a10b9...7fc2`) | [0.2544361553, 0.2544361561] | 0.25444 |

What rests on mathematics: that C82 <= B for any certificate meeting the hypotheses (the direct argument in the
certificates' notes, read step by step), Gershgorin's theorem and the series bounds.

Reproduce: `cd c82 && ./fetch.sh ROOT && ./run_all.sh ROOT --workers=14` (about 5.5 minutes).

## C16: the Brezis-Gallouet-Wainger remainder constant on the torus (5 November 2026)

`c16/c16_check.py`, standard library only, written from the definition on the 16a page ($\mathbb T^2 = [0, 2\pi]^2$,
unnormalized $L^2$ norms). Both sides of the update of 4 November in
[#270](https://github.com/teorth/optimizationproblems/pull/270):

- Upper bound, from the certificate's 353 $\delta$-intervals, each with its own $\mu$. Cauchy-Schwarz with weights
  $\lvert k\rvert^2(1 + \mu\lvert k\rvert^2)$ gives $4\pi\lvert u(0)\rvert^2/\lVert\nabla u\rVert^2 \le G(\mu)(1 + \mu\delta)/\pi$
  with $G(\mu) = \sum_{k \ne 0} 1/(\lvert k\rvert^2(1 + \mu\lvert k\rvert^2))$. Minus $\ln\delta + \ln(1 + \ln\delta)$ this is
  convex in $\delta$, so each interval is checked at its two ends. $G$ is summed by rows,
  $\sum_y 1/(y^2 + s^2) = (\pi/s)\coth(\pi s)$; rows past $N$ are bracketed by convexity with closed-form integrals; 50-digit
  decimal arithmetic, refined until $G$ is enclosed to $10^{-19}$ near the peak. For $\delta \ge 2^{18}$ we use our own
  bound: the lattice count $\#\{0 < \lvert k\rvert \le \rho\} \le \pi(\rho + 1/\sqrt 2)^2$ with $\mu = 1/(\delta\ln\delta)$.
- Lower bound: the untruncated conditional extremal $c_k = 1/(\lvert k\rvert^2(1 + \mu\lvert k\rvert^2))$ at the pull
  request's $\mu$ = `0x1.f42a1d5c1aabdp-4`, with $S_1 = G - H/\mu$ and $S_2 = H/\mu^2$,
  $H = \sum_{k \ne 0}(\lvert k\rvert^2 + 1/\mu)^{-2}$, enclosed the same way.

| Side | Certificate (SHA-256) | Our result | Printed bound |
| --- | --- | --- | --- |
| upper | vidigoat/optimization-records @5da2dc1 `tao/certificates/c16_pinned/certificate.json` (`2783dfb9...c88b`; no hash stated, this is the file at that commit) | largest interval-end value 2.15622552815421974 (2.6e-16 below the printed bound); tail $\delta \ge 2^{18}$ at most 2.0163 | 2.15622552815422 |
| lower | the trial function above | 2.15622552815421292 (enclosure width 2.7e-17; $\delta$ = 3.928836059681948) | 2.1562255281542128 |

So $2.15622552815421292 \le C_{16} \le 2.15622552815421974$ for the page's "Equivalently" definition, the maximum
over $\delta \ge 1$ (which reading is meant is asked in [#325](https://github.com/teorth/optimizationproblems/issues/325));
the upper bound holds for every $\delta \ge 1$, so it also covers a reading with $\delta$ large only. Red controls, all as expected: a target $10^{-17}$ below the
largest interval-end value is refused, one $\mu$ raised by 1% is refused, and one interval removed (a gap) is refused.

What rests on mathematics: Cauchy-Schwarz, translation invariance, the convexity in $\delta$, the row-sum identities,
the convexity of the row tails and the lattice count.

Reproduce: `cd c16 && ./fetch.sh && ./run_all.sh` (under a minute).

## C68: Korenblum's constant (5 November 2026)

`c68/c68_check.py`, standard library only: the quantitative hypotheses of the 68a row [Pin2026], $C_{68} \ge 54/125$
(added in [#239](https://github.com/teorth/optimizationproblems/pull/239)). Wikström's moment lemma (arXiv:2607.17748,
Lemma 3.1, as the certificate's NOTE.md recalls it) needs a probability measure $\mu$ on $[c^2, 1]$ with no atom at 1 and
$\int t^k\,d\mu \ge 1/(k+1)$ for $k \ge 0$, $\int \beta\, t^k\,d\mu \le 1/(k+1)$ for $k \ge 1$, where $\beta = 1 - g$ and
$g = 1/G$ is the envelope of the certificate's GEOMETRIC_LEMMA.md. We check, from the exact rationals of
`certificate-jung.json`: the measure (support, positive weights, total mass exactly 1); the admissible range of the
linear slope $l$ with $g(t) \ge l(1-t)$ on $[R, 1]$, re-derived ($y \ge a(1-t)^2$ from $-\log t \ge 1-t$ and
$\sin v \ge v - v^3/6$; $g$ increasing in $y$); $\beta$ at the twelve atoms; both moment families at every finite index in
320-bit fixed point with outward rounding; and both infinite tails (the atom at $R$ for the first family; the strip
$[t_0, t_1]$ for the second).

| Row | Certificate (SHA-256) | Our result | Printed bound |
| --- | --- | --- | --- |
| [Pin2026] | rpinho/korenblum-c68-certificate @9d48de6 `certificate-jung.json` (`80a056f2...002b`, as stated) | all hypotheses hold: slope $l$ = 0.02563746486913 below the analytic bound 0.025637464869136; first family $k \le 31278$ explicit (smallest margin $A_k/(R^{k+1}/(k+1))$ = 2.0961e-6, at $k$ = 5714), then the atom at $R$; second family $k \le 23999$ explicit (smallest margin $B_k/\sum_j w_j\beta_j t_j^k$ = 2.2512e-7, at $k$ = 8502), then the strip (gap 2.1063e-12) | 0.432 |

The two smallest margins were also recomputed in exact rational arithmetic (2.096105e-6 and 2.251249e-7). Red controls,
all as expected: half the last atom's mass moved to the one before, the slope raised by 1%, a false second-tail cutoff
at $k$ = 2000, and $c$ raised to 0.44 (first node moved to the new $c^2$) are each refused.

What rests on mathematics: Wikström's Lemma 3.1, the envelope of GEOMETRIC_LEMMA.md (two covering zeros, hyperbolic Jung,
the coefficient identity), and the elementary bounds above.

Reproduce: `cd c68 && ./fetch.sh && ./run_all.sh` (seconds).

## Reproduce

```
./fetch.sh      # downloads the certificates at the pinned commits and checks every SHA-256
./run_all.sh    # writes results/*.txt (Python 3.9+, standard library only; about two minutes)
```

Prepared with AI assistance (Claude, Anthropic), as the repository's AI use policy asks. Questions or corrections:
open an issue here. Arhan Canli, Quantavane (quantavane.com).
