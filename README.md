# Quantavane replays

Independent exact replays of certificates submitted to the
[optimization constants repository](https://github.com/teorth/optimizationproblems). Each replay uses our own code,
written from the definition on the constant page. The certificates' own verifiers were neither read nor run, and we
have no connection with any of the authors.

An index of every asterisked bound in that repository, with the replays posted for each (ours and others'), is in
[ASTERISK_INDEX.md](ASTERISK_INDEX.md) (snapshot of 5 November 2026).

## C1a, C1c and C1d (3 to 5 November 2026)

| Constant | Pull request | Certificate (pinned commit, SHA-256) | Exact value | Printed bound |
| --- | --- | --- | --- | --- |
| $C_{1c} \geq$ | [#292](https://github.com/teorth/optimizationproblems/pull/292) (Koshta's row) | DennisKoshta/c1c-autocorrelation-certificate @ `a4ee94ef29ff`, `K24576_0.4103333.json`, `963cb727...c9e` | 134648548624804778248790016 / 328144294688458468928334649 = 0.41033335274849333953..., minimum at k = 21264 only | 0.4103333, below the exact value |
| $C_{1c} \geq$ | [#318](https://github.com/teorth/optimizationproblems/pull/318) | Julsgaard/c1c-certificate @ `91e47b1a`, `cert/certificate.json` (K = 98,304, 98,368 weights), `fe11336b...52d2` | 2582127846607416062611587072/6292441470965870058830439361 = 0.4103538918115781099..., minimum only at k = 47318: the printed 0.4103538918 is a valid rounding down |
| $C_{1c} \geq$ | [#270](https://github.com/teorth/optimizationproblems/pull/270) (update of 4 Nov) | vidigoat/optimization-records @ `9f82508a`, `tao/certificates/c1c_v2/c1c_K98304_v2.json` (K = 98,304), `3bba62b0...4729` (no hash stated; this is the file at that commit) | 645458795582745374892957696/1572931393892132517694378321 = 0.4103540676275733581..., minimum only at k = 93282: the printed 0.4103540676 is a valid rounding down |
| $C_{1c} \geq$ | [#270](https://github.com/teorth/optimizationproblems/pull/270) (update of 5 Nov) | vidigoat/optimization-records @ `42d6a89f`, `tao/certificates/c1c_v3/c1c_K98304_v3.json` (K = 98,304, 98,368 weights), `e0786a59...ec66` | 640823574552569898220412928/1561614482809354047911104009 = 0.4103596512499835112..., minimum only at k = 98301: the printed 0.4103596512 is a valid rounding down |
| $C_{1a} \leq$ | [#270](https://github.com/teorth/optimizationproblems/pull/270) | vidigoat/optimization-records @ `9f82508a`, `tao/certificates/c1a/c1a_n262144.json` (262,144 cells, decimal strings read exactly), `1848a920...f3de` | 3756858380445445962687721679242267976826786797829957539923492864000000000000/2500000000000000227858304266031237715040682299126697825903032255969919203361 = 1.50274335217817824811..., maximum only at m = 100329 | 1.5027433522, above the exact value |
| $C_{1a} \leq$ | #270 (its first refinement) | same commit, `c1a_n131072.json` (131,072 cells), `ad1b7ff4...8bef` | 88396668012481333418540698869863464003678707923912272976102136217600000/58823529411764719617703895028066775569135640608648734518121194077814337 = 1.50274335621218231722..., maximum only at m = 131071 | 1.502743356212 (its README), above |
| $C_{1d} \leq$ | [#312](https://github.com/teorth/optimizationproblems/pull/312) | Julsgaard/c1d-certificate @ `cc455ca0`, `cert/certificate.json` (81,920 cells), `04c9c759...237` | 9900630255081376929482866688 / 6805830738283663732311831361 = 1.45472766452875064884..., maximum only at i + j = 29124 | 1.454727664528751, above the exact value |
| $C_{1d} \leq$ | [#271](https://github.com/teorth/optimizationproblems/pull/271) | Keremozdemirra/tao-1c-1d-certificates @ `408f8d0b`, `c1d_n40960_1.4556234448.json` (40,960 cells), `be35060a...2da` | 279676466330630562685427712000 / 192135175720687825901095191121 = 1.45562344470022455281..., maximum only at i + j = 15313 | 1.4556234448, above |
| $C_{1d} \leq$ | [#269](https://github.com/teorth/optimizationproblems/pull/269) | echul-seal/c1d-upper-bound @ `2be42dfd`, `ub/c_f10240.json` from `c1d_records.tgz` (10,240 cells), `30d501e1...973` | 76358168487433178963290071040 / 52455219864513498497744894329 = 1.45568293650581519488..., maximum only at i + j = 3168 | 1.455682936505816, above |
| $C_{1d} \leq$ | [#286](https://github.com/teorth/optimizationproblems/pull/286) | huiihao/c1d-certificate @ `302817dc`, `construction.json` (400 cells, JSON floats), `c489b46b...8df` | read as binary64: 4563650116209721766004955898613711667200 / 3127595284199802099438492053245424408273 = 1.45915622115964902113...; read as the printed decimals: 1.45915622115964903535... | 1.4591562211596491, above both |

Every hash above is the one stated in the pull request (for #269 also the archive's, `b526326d...`), except the 4 Nov
#270 row and the two $C_{1a}$ rows, where none was stated (these are the files at the commit the pull request links).

## Method

- $C_{1c}$: $f = \sum_i w_i 1_{[i/K,(i+1)/K)}$ with non-negative integer weights. The autocorrelation is linear between
  multiples of $1/K$, so the bound is $K \min_{0 \le k \le K} A_k / S^2$ with $A_k = \sum_i w_i w_{i+k}$ and
  $S = \sum_i w_i$ (we check $n > K$, so no $A_k$ in the range is an empty sum).
- $C_{1d}$: $n$ equal cells on $[-1/4,1/4]$, so $h = 1/(2n)$. $f*f$ is linear between multiples of $h$ and equals
  $h \sum_{i+j=m-1} v_i v_j$ at the knots, so the bound is $2n \max_m \lvert c_m \rvert / S^2$ with $c = v * v$.
- $C_{1a}$: the same knot formula for a non-negative step function on $[-1/4,1/4]$ (non-negativity is checked), since
  $\max_t (f*f)(t) \ge C_{1a}(\int f)^2$ for every such $f$ gives $C_{1a} \le 2n \max_m c_m / S^2$.
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


PR #270, update of 5 Nov (vidigoat c86_v2 @fd4373c, `certificate.json`, `09c05af0...8d1f` as stated): the format of #298 with
#298's 28 auxiliaries plus 3, the density re-optimised on 5,734 knots with dyadic edges (denominators up to 2^57; our
checker now takes the logarithms term by term when the edge grid is too fine to tabulate). Every condition E1-E6 holds:
T/m = 195153337643980465525925899810927/108086391056860991292157728915456 = 1.8055310732071364922... exactly, below the
printed 1.8055311 (the README's decimal 1.80553107320710 is 3.6e-14 low); smallest log integral 9.09e-8 (Q23); f certified
positive on all 5,732 support cells, min f in [1.0049371e-6, 1.0049372e-6]. `python3 -I c86_check.py vid2 certs --controls`.

PR #270, night update of 5 Nov (vidigoat c86_v3 @9e2b86c, `certificate.json`, `182e76bf...d397` as stated): the same
format with 33 auxiliaries and 6,999 knots. Every condition E1-E6 holds: T/m =
780597977251274792019360777060007/432345564227458232572250425393152 = 1.80549551525085145363... exactly, as the
folder's README states, 8.5e-8 below the printed 1.8054956; sum b deg Q / m = 0.8987818231; smallest log integral
[9.05960954785e-8, 9.05960954786e-8] (Q22); f certified positive on all 6,997 support cells (247,249 pieces), min f in
[1.0085284e-6, 1.0085286e-6]. `python3 -I c86_check.py vid3 certs --controls` (about 11 minutes on 12 cores).

The lower bound $C_{86} \ge 1.80220$\* (D. Niedbala Giraudin, arXiv:2609.29298, the row proposed in
[#315](https://github.com/teorth/optimizationproblems/pull/315)), from the certificate `CERTIFICAT_C86.json` (an
ancillary file of arXiv v1, `3caa7297...abe5`) and the paper's Proposition 2: if
$x \ge \lambda + \sum_Q \lambda_Q \log\lvert Q(x)\rvert + \lambda_0 (2U_\eta(x) - I(\eta))$ for every $x \ge 0$, with
$\lambda_Q, \lambda_0 \ge 0$, integer $Q$ and $\eta$ a probability measure of finite energy, then $C_{86} \ge \lambda$.
`c86/c86_lower.py` with `c86_lower.c`, written from the paper (its `prove_SSS_C86.py` and `verify_C86.py` were
neither read nor run):
- the data: 18 monic integer polynomials (14 with $\lambda_Q > 0$), $\lambda_0 \ge 0$, 6,932 cells on [0, 8], masses
  $\ge 0$ summing to $1 + 7.9 \cdot 10^{-17}$ (scaled to exactly 1);
- $I(\eta)$ as a double sum over the 3,891 cells with mass, by the closed-form double antiderivative for near cells
  and its series for far ones, with explicit error bounds: $3.130515030 \cdot 10^{-5}$ ($\pm 8.9 \cdot 10^{-9}$; the
  paper gives $3.1305148803 \cdot 10^{-5}$);
- the inequality on [0, 14] by bisection (819,821 intervals): $U_\eta$ bounded cell by cell (each cell's potential
  is unimodal with its minimum at the centre), $\max\lvert Q\rvert$ by the centred Taylor form with a fourth-order
  remainder; the smallest value of the left side minus the sum seen is 1.8022211884, against the certificate's
  1.8022211894912388;
- the tail $x \ge 14$: $U_\eta(x) \le \log x$ and $\log\lvert Q(x)\rvert \le \deg Q \log x + \log K_Q$, so the
  left side minus the sum is at least an increasing $h$ with $h(14) \ge 6.7073$.

Binary64 arithmetic; libm's log is taken to be within one ulp (8 ulp allowed), every sum carries its rounding
bound, and the margin is $2.1 \cdot 10^{-5}$. Red controls, all as expected: $\lambda$ one millionth above the
certificate's own value is refused (at $x = 0.0795$), dropping the energy term ($\lambda_0 = 0$) is refused at
1.80220, and a negative weight is refused. `python3 -I c86_lower.py ng certs --controls` (about 5 minutes on 12 cores
with the controls).

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
| #330 | rpinho/math-certificates @5461d49 `c82/degree104-cvp/certificate.json` (`27a4e611...d3f9`) | [0.2524048059550536, 0.2524048059550537] | 0.252404806 |
| #330, updated 5 Nov | rpinho/math-certificates @98208a4 `c82/degree104-joint/certificate.json` (`a4061273...e36e`) | [0.2520179247552093, 0.2520179247552094] | 0.2520179248 |
| #330, final certificate (5 Nov) | rpinho/math-certificates @3be6b2f `c82/degree104-joint/certificate.json` (`b71150bb...4171`) | [0.2519964889882011, 0.2519964889882012] | 0.2519964890 |
| #280 v4, 5 Nov | vidigoat/optimization-records @79ce0df `tao/certificates/c82/v4/certificate.json` (`c1981f1c...c72f`) | [0.2520079723284296, 0.2520079723284297] | 0.2520079724 |
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
| lower, [#319](https://github.com/teorth/optimizationproblems/pull/319) (Julsgaard) | the same family at $\mu = 128042/2^{20}$ | [2.15622552815420610, 2.15622552815420614] ($\delta$ = 3.928838101027620) | 2.156225528 (stated enclosure 2.1562255281542061 $\pm$ 6.1e-17) |
| upper, #270 night update (18 digits) | vidigoat/optimization-records @1e8e786 `tao/certificates/c16_pinned_v2/certificate.json` (`8c13abbb...05e8`, as in its SHA256SUMS) | all 383 intervals hold: the largest interval-end value is 2.15622552815421295831 (G bracket width 2.0e-21 there), 1.7e-18 below the printed bound | 2.15622552815421296 |
| lower, #270 night update | the same family at $\mu$ = `0x1.f42a1bd5af466p-4` | 2.15622552815421295612 (enclosure width 6.7e-21; 262,144 explicit rows) | 2.156225528154212956 |

So $2.15622552815421295612 \le C_{16} \le 2.15622552815421295831$ (from the night update; $2.15622552815421292 \le C_{16} \le 2.15622552815421974$ from the first) for the page's "Equivalently" definition, the maximum
over $\delta \ge 1$ (which reading is meant is asked in [#325](https://github.com/teorth/optimizationproblems/issues/325));
the upper bound holds for every $\delta \ge 1$, so it also covers a reading with $\delta$ large only. Red controls, all as expected: a target $10^{-17}$ below the
largest interval-end value is refused, one $\mu$ raised by 1% is refused, and one interval removed (a gap) is refused.

What rests on mathematics: Cauchy-Schwarz, translation invariance, the convexity in $\delta$, the row-sum identities,
the convexity of the row tails and the lattice count.

Reproduce: `cd c16 && ./fetch.sh && ./run_all.sh` (about 7 minutes; the night update refines G and H near the peak with
`--upper-width 1e-21 --max-rows 262144 --lower-factor 128`).

## C68: Korenblum's constant (5 November 2026)

`c68/c68_check.py`, standard library only: the quantitative hypotheses of the 68a row [Pin2026], $C_{68} \ge 54/125$
(added in [#239](https://github.com/teorth/optimizationproblems/pull/239)), and of the two certificates of
[#274](https://github.com/teorth/optimizationproblems/pull/274) (0.446 and 0.44605). Wikström's moment lemma (arXiv:2607.17748,
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
| [Pin2026] | rpinho/korenblum-c68-certificate @9d48de6 `certificate-jung.json` (`80a056f2...002b`, as stated) | all hypotheses hold: slope $l$ = 0.02563746486913 below the analytic bound 0.025637464869136; first family $k \le 31278$ explicit (smallest margin $A_k/(R^{k+1}/(k+1))$ = 2.0961e-6, at $k$ = 5714), then the atom at $R$; second family $k \le 23999$ explicit (smallest margin $B_k/(R^{k+1}/(k+1))$ = 2.2848e-7, at $k$ = 8502), then the strip (gap 2.1063e-12) | 0.432 |
| #274 | vidigoat/optimization-records @1a76f9d `tao/certificates/c68/certificate-c0446.json` (`2b04f03e...d019`, as in its SHA256SUMS) | all hypotheses hold: $A$ = 12.2233841499591, $\Psi$ = 25436.8054244252, $\mathrm{sech}^2(A/2) < 1/3$, $8\cosh^4(A/6) < \Psi$; first family $k \le 5806$ explicit (smallest margin 6.9901e-6, at $k$ = 407), then the atom at $R$; second family $k \le 1322$ explicit (smallest margin 6.9824e-6, at $k$ = 44), then one tail piece alone from $k$ = 1323 | 0.446 |
| #274, supplementary | same commit, `certificate-c044605.json` (`b8ea516b...054d`, as in its SHA256SUMS) | all hypotheses hold; margins 8.4733e-7 (at $k$ = 9) and 8.1559e-7 (at $k$ = 44); one tail piece alone from $k$ = 1324 | 0.44605 |

The two smallest margins were also recomputed in exact rational arithmetic (2.096105e-6 and 2.2847831e-7); every figure
agrees with the certificate's NOTE_JUNG.md. This certificate was replayed before by the two checkers of
[#274](https://github.com/teorth/optimizationproblems/pull/274) (vidigoat, 26 October); ours shares no code with them or
with the certificate's own. Red controls,
all as expected: half the last atom's mass moved to the one before, the slope raised by 1%, a false second-tail cutoff
at $k$ = 2000, and $c$ raised to 0.44 (first node moved to the new $c^2$) are each refused.

For #274 the envelope is $G = (1 + \sqrt{1 + \Psi/(2\sigma)})/2$ (its PROOF.md, Lemma 6), and on each of its 64 tail pieces
$[a, b]$ we derive our own slope: $-\log t \ge \delta(1 + \delta_{lo}/2)$ and $\sin x \ge x(1 - U^2/6)$ give $\sigma \ge (b_i\delta)^2$, and
$g$ is increasing in $\sigma$. Our second-tail start is found by the checker (the first $k$ at which one piece alone covers
all atoms), here 27 below #274's 1350; its smallest margins agree with those #274 states. Controls for #274: mass moved
between atoms 8 and 7, and between 2 and 3, and $c$ raised to 0.447 (its PROOF.md puts the framework's limit near
0.44607) are each refused.

What rests on mathematics: Wikström's Lemma 3.1, the envelopes (for [Pin2026] GEOMETRIC_LEMMA.md: two covering zeros,
hyperbolic Jung, the coefficient identity; for #274 PROOF.md: the flux and half-diameter lemmas), and the elementary
bounds above.

Reproduce: `cd c68 && ./fetch.sh && ./run_all.sh` (seconds).

## Ramsey lower bounds since Small Ramsey Numbers revision 18 (4 November 2026)

`ramsey-rev19/`: every new lower bound for a classical two-colour Ramsey number R(s,t) since revision 18 (24 April
2026) of Radziszowski's survey, each witness taken from its authors' own repository or data record and checked by two
exact clique programs (`ramsey_check`, Bron-Kerbosch with Tomita's pivot; `ramsey_indep`, MCQ branch and bound), with
one-edge red controls. Results, notes and the witnesses still too large for our budget are in
[`ramsey-rev19/README.md`](ramsey-rev19/README.md).

Reproduce: `cd ramsey-rev19 && ./fetch.sh ROOT && ./run_all.sh ROOT`.

## C3a: the Gyarmati-Hennecart-Ruzsa sum-difference constant (5 November 2026)

`c3a/c3a_check.py` with `c3a_hist.c`: the finite $k = 2$ controlled-carry bound of [Syl2026] ($C_{3a} > 1.19281555$\*,
[#238](https://github.com/teorth/optimizationproblems/pull/238)), rebuilt from its theorem note alone (mask
$\langle 1971, 2016, 2100, 2628, 2688, 2800\rangle \cap [0, 26972]$, $q = 27022$, $x$ = `0x1.ffde1bafccc06p-1`, $\ell = 10^{18}$):
- $M$ and $M + M$, and the min-cost histogram of all $q^2 = 730{,}188{,}484$ residues. Each residue is a min over at most
  four digit pairs, since $B < q$; a full scatter of $D \times D$ gives the same histogram on small cases. The histogram
  is byte-identical to the package's `k2_cost_hist.tsv`.
- $Z_2(x)$, the type counts and $C$ in exact integer arithmetic: $N_0 = 3362808316958$ and $C = 26096544438312721237750$,
  as stated, and every $N_r \ge 1$.
- The bound in 90-digit decimal interval arithmetic with Robbins' bounds on $\ln n!$.

| Row | Data | Our result | Printed bound |
| --- | --- | --- | --- |
| 3a [Syl2026] $k = 2$ | xSyllaa/c3a-controlled-carry-k2 @5ca7cb8 `proof/THEOREM.md` (the histogram file `4d16010f...6c63` only compared) | $\theta \ge 1.1928155570929862489856\ldots$; every certified quantity agrees with the note to the digits it prints | 1.19281555 |

Red controls, all as expected: the claim raised to 1.19281556 is not shown; dropping the generator 2800 changes the
data; $x$ one unit lower in its last bit changes $N_0$ and $C$.

What rests on mathematics: Lemmas A-D of the note (the encoding, injectivity, the sumset bound and [GHR2007] with
$V = 2U$), Kleinwaks's [K2026b] Proposition 3, and Robbins' bounds. The $k = 3$ instance (1.19359680\*) needs the
min-cost histogram of $q^3 \approx 2 \cdot 10^{13}$ residues and is not replayed here.

Reproduce: `cd c3a && ./fetch.sh && ./run_all.sh` (about 2 minutes).

## C89b and C89c: Ramsey multiplicity of $K_5$ and $K_6$ (5 November 2026)

The minimum density of monochromatic $K_t$ over 2-coloured graphons; an upper bound is the exact density of one step
graphon, $t(K_t, W) + t(K_t, 1 - W)$ with repeated blocks counted.

`c89/c89_check.py` (89c, $t = 6$): step graphons on the 1024 blocks of $\mathbb Z_2^{10}$ whose red density depends only on
$w = \lvert x \oplus y\rvert$ ($n_w / D$). With the first vertex at 0, each pair distance is the number of coordinates
whose column separates the pair; a red class ($n_w = D$) rules out blue, a blue class ($n_w = 0$) rules out red, and the
fractional classes give $f^e$ and $(1-f)^e$. Two C programs written apart count the $2^{50}$ tuples behind each
(colour, exponents) cell: `c89_multiset.c` over the $\binom{41}{10}$ multisets of column types with multinomial weights,
`c89_refine.c` by placing the vertices one at a time on cells of coordinates, pruning once both colours are ruled out.
The tables must agree exactly and each must sum to $2^{50}$; the density is then formed in exact rationals.

`c89/c89b_check.py` (89b, $t = 5$): the 256-block binary Cayley matrix of #307 is read in full and checked (equal
weights, entries $0$, $2648$, $4041$ over $4041$, symmetric, $R[x][y] = R[0][x \oplus y]$ for all 65,536 pairs, row 0 as
the two-mask definition says). `c89b_count.c` counts the 5-tuples behind each $p^k$ and $(1-p)^k$ with the first vertex
in block 0, and again over every ordered tuple with no symmetry assumed, which must give exactly 256 times as much.
For #277 (192 blocks, every block-pair density free), `c89b_general.c` assumes no symmetry: the numerators
$N(R)$ and $N(D - R)$ over all $192^5$ tuples, both as $\sum_{a,b} R_{ab}\operatorname{tr}(S^3)$ with $S = \operatorname{diag}(R_a \circ R_b) R$
and as $\sum_{c,d,e} R_{cd}R_{de}R_{ec}\, v^{\mathsf T} R v$ with $v_x = R_{xc}R_{xd}R_{xe}$, modulo ten primes below $2^{28}$; the two agree at
every prime, the numerators are rebuilt by the Chinese remainder theorem from nine and the tenth confirms.

| Row | Data (SHA-256) | Exact density | Printed upper bound |
| --- | --- | --- | --- |
| 89c [SK2026], [#279](https://github.com/teorth/optimizationproblems/pull/279) | JeetuSK0808/tao-constants-certificates @98b6ee3 `construction.json` (`f90cc190...ac49`, as in its SHA256SUMS.txt) | identical to the stated fraction: $0.742254861593625\ldots \cdot 2^{-14}$ | $0.7422548616 \cdot 2^{-14}$ |
| 89c [SCIEL2026], [#331](https://github.com/teorth/optimizationproblems/pull/331) | the eleven $n_w / 2^{32}$ stated in the pull request (the Zenodo record is restricted) | identical to the stated fraction: $0.742204211528050\ldots \cdot 2^{-14}$ | $0.7422042116 \cdot 2^{-14}$ |
| 89b round 4, [#307](https://github.com/teorth/optimizationproblems/pull/307) | rpinho/math-certificates @1d38853 `c89b/round4/certificates/c5_cayley256_fine.json.gz` (`78978a9f...0125`, as in its SHA256SUMS) | identical to the stated fraction 15058298592194934908767847/9131939205947984819238469632 $= 0.0016489705256017\ldots$; as a polynomial in $p$, exactly the stated $(8734 - 5530p + \ldots + 10p^6)/2^{22}$; at $p = 19/29$ the stated simpler alternative | 0.001648970526 |
| 89b [SK2026], [#277](https://github.com/teorth/optimizationproblems/pull/277) | JeetuSK0808/tao-constants-certificates @82a6e4d `c89b_stepgraphon_m192_D20.json.gz` (`19c43354...bc06`, as in its SHA256SUMS.txt) | identical to the stated fraction: $0.0016870051344094\ldots$ | 0.0016870052 |
| control: Franek 2002 | red iff $\lvert x \oplus y\rvert \in \{0, 1, 3, 4, 7, 8, 10\}$ | exactly $51162598917 / 2^{50}$, as published | |
| control: Parczyk-Pokutta-Spiegel-Szabo 2024 | their 192-vertex graph, `c5.graph6.txt` in Zenodo record 6602512 (`graphs.zip`, `6ff8a249...ec89`) | exactly $2320651 \cdot 2^{-24} \cdot 3^{-4}$, as their Theorem 1.2 states | |

Self-tests, all passing: brute force over every tuple on small cases against the programs ($t = 4, 5, 6$; Cayley and
arbitrary symmetric matrices), and $f \equiv 1/2$ giving exactly $2^{1 - \binom{t}{2}}$. Red controls, all as expected: a
bound one last-digit unit lower is not implied; a density above 1 is refused; one density raised by one unit changes
the fraction; for #307, one entry changed breaks the translation invariance, a second fractional value and an
asymmetric matrix are refused.

What rests on mathematics: that a step graphon's monochromatic $K_t$ density bounds the constant from above (random
colourings of blow-ups), translation invariance on $\mathbb Z_2^n$, and the counting identities above.

Reproduce: `cd c89 && ./fetch.sh && ./run_all.sh` (about 20 minutes).

## Reproduce

```
./fetch.sh      # downloads the certificates at the pinned commits and checks every SHA-256
./run_all.sh    # writes results/*.txt (Python 3.9+, standard library only; about eight minutes)
```

Prepared with AI assistance (Claude, Anthropic), as the repository's AI use policy asks. Questions or corrections:
open an issue here. Arhan Canli, Quantavane (quantavane.com).
