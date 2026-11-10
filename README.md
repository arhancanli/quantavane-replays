# Quantavane replays

Independent exact replays of certificates submitted to the
[optimization constants repository](https://github.com/teorth/optimizationproblems). Each replay uses our own code,
written from the definition on the constant page. The certificates' own verifiers were neither read nor run, and we
have no connection with any of the authors.

An index of every asterisked bound in that repository, with the replays posted for each (ours and others'), is in
[ASTERISK_INDEX.md](ASTERISK_INDEX.md) (snapshot of 5 November 2026).

## C1a, C1c and C1d (3 to 7 November 2026)

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
| $C_{1d} \leq$ | [#326](https://github.com/teorth/optimizationproblems/pull/326) (the STEPS2026 row) | steps-re/c1d-certificate @ `9a091494`, `construction.json` (600 cells, JSON floats), `3dd39fae...a28c` | read as binary64: 698746209786798359455483373677229695641600 / 478522078482569620433992386249161163190849 = 1.46021728402287398524...; read as the printed decimals: 1.46021728402287397507... | 1.460217284022874, above both |
| $C_{1d} \leq$ | [#326](https://github.com/teorth/optimizationproblems/pull/326) (the AER2026 row) | unclegareth-create/alphaevolve_repository_of_problems @ `4b3bbd2a`, `.../problem_4_signed_autoconvolution_bounds.py`, `1a290ad1...54b2`: the list `STEP_FUNCTION_C4_PRIME_N600` (600 cells), read from the file with Python's `ast` parser and never run | read as binary64: 1.45998904690538324315...; read as the printed decimals: 1.45998904690538314495... | 1.459989046905384, above both |
| $C_{1d} \leq$ | [#224](https://github.com/teorth/optimizationproblems/pull/224) and [#243](https://github.com/teorth/optimizationproblems/pull/243) (the two CERT2026 rows; 7 November) | extremal-values-registry @ `50520b9`, `certificates/third_autocorrelation_abs_inside/construction.json` (300 cells, JSON floats; its values equal those printed in #224), `51f864ee...0a0e` (as the 1d page states) | read as binary64: 2039224160683882663565441886346340189798400 / 1393796574908163689489385839036226787118881 = 1.46307158260755932746..., the fraction #224 states, maximum only at m = 414; read as the printed decimals: 1.46307158260755929963... | the starred 1.463071582607559 is below both, so it is not an upper bound as printed; #243's 1.4630715826075594 is above both |
| $C_{1c} \geq$ | [#243](https://github.com/teorth/optimizationproblems/pull/243) (the [R2026] row, whose asterisk came off on one replay; 7 November) | rozetyp/barnard-steinerberger-lower-bound @ `19f3971`, `certificate.json` (K = 6,144, 6,148 integer weights), `62ac49dd...fd0e` (as the 1c page states) | 21072342208353792 / 51354782902387969 = 0.41032871754918741449..., the page's fraction, minimum only at k = 4665 | 0.4103287, below the exact value |
| $C_{1a} \leq$ | [#232](https://github.com/teorth/optimizationproblems/pull/232) (the [A2026] row, whose asterisk came off on one replay; 7 November) | atrislabs/records @ `1fa857a0`, `1a-sidon-autocorrelation/construction.json` (65,536 non-negative rational heights), `7387452d...7abf` (the page states `7387452d1da7`) | 645423412969467179759478160871915520 / 429496729600000035389440000000000729 = 1.50274348670958337974..., the fraction of the repository's `score_output.json`, maximum only at m = 25072 | 1.5027434868, above the exact value |

Every hash above is the one stated in the pull request (for #269 also the archive's, `b526326d...`), except the 4 Nov
#270 row and #270's two $C_{1a}$ rows, where none was stated (these are the files at the commit the pull request links),
and CERT2026, whose hash is the one the 1d page states for the registry copy.

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

## C15b: the dual matrix multiplication exponent (6 November 2026)

`c15b/c15b_check.py`, standard library only, interval arithmetic of our own (fixed point with 448 fractional bits,
outward rounding, exp and log by series with explicit remainders): the two rows proposed for the 15b page in
[#337](https://github.com/teorth/optimizationproblems/pull/337), $\alpha > 0.465$\* ([OAI2026], the seven-step
completion script (12.8)) and $\alpha > 0.46652506$\* ([Var2026], a 15-step script, `dual-alpha-witness` @d2868fa,
`certificate.json`, `8a81020f...d345`), written from Section 12.1 of the preprint and the certificate data alone
(Varlamov's `verify.py` was neither read nor run). Base two, distinguished side $X$ and its colour $B$; each completion
of size $m$ updates $(\rho, b, h_A, h_B, h_C)$ by Lemma 12.1, and $\alpha \ge \min\{1, 2D/\log 2\}$ with
$D = (\rho h_B + (1-\rho)\max(h_A, h_C))/b$.

| Check | Result |
| --- | --- |
| the 15-step script | sizes multiply to $b = 147456$, as the certificate says |
| $2D/\log 2$, reduced recurrence (12.9) and full Lemma 12.1, which agree | $0.466525061423431426876062217\ldots$ (enclosure width $6 \times 10^{-131}$) |
| the claim | $2D/\log 2 > 23326253/50000000$, by exact comparison of the interval ends |
| the paper's script (12.8) | $2D/\log 2 = 0.466101554898471852177\ldots > 0.465$; centres B, A, C, B, C, B, A as stated |
| the certified information after (12.9) and the chain (12.10) | all eleven statements hold ($\rho$ exactly, as a dyadic rational), and the chain equals $1548637/3327360 > 0.465$ |

Red controls, all as expected, through the same check: the 15-step bound raised by $10^{-8}$ is not implied, the
script without its last step gives $0.4665195\ldots$, its fifth step re-centred on the smaller inactive colour gives
$0.4613601\ldots$, and the paper's script does not reach the 15-step bound. What rests on mathematics: the preprint's
completion and entropy arguments (Lemma 10.2, Theorem 9.3 and the output-size bound behind (12.7), and the proof of
Lemma 12.1), which we read but did not formalize.

Longer scripts (our own data, `c15b/longer_scripts.json`, checked by `c15b_check.py longer`): a beam search over every
script with completion sizes 2 to 8 and every choice of centre (B, or either inactive colour, which Lemma 12.1 also
allows) found nothing above $0.4665256234$. Its best 15-step script is the [Var2026] script with $-3$ as the 14th step,
$2D/\log 2 = 0.466525082245893813\ldots$, and its best overall needs 42 completions ($b \approx 8.0 \times 10^{15}$),
$2D/\log 2 = 0.466525623382331283\ldots$. So within this construction the 15-step script is within
$6 \times 10^{-7}$ of the best value found.

Reproduce: `cd c15b && ./fetch.sh && ./run_all.sh` (under a second).

## C5a: a Sidon set constant, Erdős problem 30 (6 November 2026)

`c5a/c5a_check.py`, standard library only, exact integer and rational arithmetic: the bound $C_{5a} \le 0.9428095$\*
([Akw26] on the 5a page) through the independent certificate of
[#335](https://github.com/teorth/optimizationproblems/pull/335) (A. Varlamov, `sidon-c5a-certificate` @091844e,
`certificate.json`, `fece8394...80ab`), written from its `PROOF.md` alone (its `verify.py`, `audit_checker.py`,
`build_certificate.py` and `reproduce.py` were neither read nor run). Hou and Zhao's vector-valued smoothing in the
symmetric case: centred uniform kernels $p^{(\ell)}$ on $[n-\ell, n+\ell)$ with weights $\lambda_\ell = 4\ell^3/n^4$
($\ell < n$) and $(2n-1)/n^2$, $n = 512$, $m = 1024$, boundary depth $L = 12$, and boundary weights
$W^{(\ell)}_j = \lfloor \ell^{-1}\sum V_k \rfloor + r$ rebuilt from the 6,144 integers $V_k$ and $r = 43922$ (our
window sums use prefix sums over the even sites).

| Check | Result |
| --- | --- |
| kernels and weights | nonnegative, symmetric, masses 1, $\sum_\ell \lambda_\ell = 1$ |
| covering inequalities, $0 \le q \le M = 12288$, over the denominator $2n^4D$ | all 12,289 hold; the least sum is exactly 1, at $q = M$ only |
| $a$, $b$ | $a = 349525/262144$, $b = 0.666667938238285\ldots > 0$, both equal to the note's exact fractions |
| bound | $\sqrt{ab} = 0.94280949115162893061\ldots$, so $ab \le 0.942809491151629^2 < 0.9428095^2$ by cross-multiplication |

Red controls, all as expected: the repair lowered to $r - 1$ breaks one covering inequality (so $r = 43922$ is the least
that works), $V_{100}$ halved breaks 1,224, and the claim lowered by one unit in its last digit is not implied. What
rests on mathematics: Hou and Zhao's Theorem 2.1 in the symmetric case, as `PROOF.md` restates it (the Sidon energy
bound, the covering sums at the boundary, weighted Cauchy-Schwarz and the choice $H \approx \sqrt{a/b}\,N^{3/4}$),
which we read step by step.

Reproduce: `cd c5a && ./fetch.sh && ./run_all.sh` (about 20 seconds).

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
| #330, update of 5 Nov evening | rpinho/math-certificates @a12c362 `c82/degree104-joint/certificate.json` (`ee9622b5...fdd5`, as in its SHA256SUMS) | [0.2519923409720083, 0.2519923409720084], inside the stated [0.25199234097200838698, 0.25199234097200839540] | 0.2519923410 |
| #280 v4, 5 Nov | vidigoat/optimization-records @79ce0df `tao/certificates/c82/v4/certificate.json` (`c1981f1c...c72f`) | [0.2520079723284296, 0.2520079723284297] | 0.2520079724 |
| #280 v5, 9 Nov | vidigoat/optimization-records @191064e `tao/certificates/c82/v5/certificate.json` (`8a5f363c...a909`, as the PR comment states); 8 denominators of degree 104, $D$ = 999,999,999,960 as stated | [0.2519889613870203, 0.2519889613870204], below the README's 0.25198896138704 | 0.2519889614 |
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
| lower, the 16a row of [#326](https://github.com/teorth/optimizationproblems/pull/326) (maximusrufus) | the trial function truncated at $0 < \lvert k\rvert \le 4000$, $\mu = 1221/10000$ (maximusrufus/c16-bgw-certificates @ `0c570700`, `lower_certificate.json`, `de314c10...6991`: only $R$ and $\mu$ used), `c16_trunc_check.py`: the three sums over the disc in integers at the scale $2^{256}$, every term rounded both ways | 2.15622477064423579034; the certificate's three sum enclosures overlap ours | 2.1562247706 |
| upper, the 16a row of [#246](https://github.com/teorth/optimizationproblems/pull/246) (maximusrufus; 7 November) | maximusrufus/c16-bgw-certificates @ `0c570700`, `upper_certificate.json` (`f46adc04...25a3`; no hash stated, this is the file at the commit #326 cites): 84 boxes tiling $[1, 4.9 \cdot 10^8]$ exactly, each with its own $\mu$ (read exactly from its decimal string). `c16_boxes_check.py` checks the 82 boxes that start below $2^{18}$ with their own $\mu$ (each endpoint rounded outward; the one that crosses $2^{18}$ up to $2^{18}$); for $\delta \ge 2^{18}$ our own tail bound replaces the certificate's | largest endpoint value 2.15622564574899300283, on the box $[3.928106, 3.931069]$ ($\mu$ = 0.12206997134993264), $8.5 \cdot 10^{-10}$ below the printed bound; tail at most 2.0163. Controls: the bound $10^{-15}$ below that value, one $\mu$ raised by 1%, and one box removed are each refused | 2.1562256466 |
| lower, the 16a row [Pat2026b] of #270 (10 November) | the trial function the row states, $c_k = 1/(\lvert k\rvert^2(\lvert k\rvert^2 + 41/5))$ on $0 < \lvert k\rvert \le 40000$ (5,026,547,528 lattice points; no file needed), `c16_trunc_big.py`: the three sums of `c16_trunc_check.py`, computed in C (`c16_trunc_fast.c`, every term floored and ceiled at the scale $2^{100}$ by unsigned 128-bit division; equal to the Python sums integer for integer at $R = 300$ and $1200$, and giving the [Pat2026] value at $R = 4000$ again); $\delta = 3.93181190026902397902$. Controls: the bound raised just above our enclosure, the cutoff halved ($2.1562254817$) and $\mu$ = 6/41 in place of 5/41 ($2.1558994050$) are each refused | [2.15622550442903480541, 2.15622550442903480542], the row's 2.15622550442903... ; 4.3e-10 above the printed bound | 2.156225504 |

So $2.15622552815421295612 \le C_{16} \le 2.15622552815421295831$ (from the night update; $2.15622552815421292 \le C_{16} \le 2.15622552815421974$ from the first) for the page's "Equivalently" definition, the maximum
over $\delta \ge 1$ (which reading is meant is asked in [#325](https://github.com/teorth/optimizationproblems/issues/325));
the upper bound holds for every $\delta \ge 1$, so it also covers a reading with $\delta$ large only. Red controls, all as expected: a target $10^{-17}$ below the
largest interval-end value is refused, one $\mu$ raised by 1% is refused, and one interval removed (a gap) is refused.

What rests on mathematics: Cauchy-Schwarz, translation invariance, the convexity in $\delta$, the row-sum identities,
the convexity of the row tails and the lattice count.

Reproduce: `cd c16 && ./fetch.sh && ./run_all.sh` (about 13 minutes, the #246 upper row about 3 and the [Pat2026b] row about 3 with a C compiler; the night update refines G and H near the peak with
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
| 3a [Syl2026] $k = 3$ | xSyllaa/c3a-controlled-carry-k2 @5ca7cb8 `k3/proof/THEOREM_k.md` (the histogram `3186f6e2...5c63` only compared) | the min-cost histogram of all $q^3 = 19{,}731{,}153{,}214{,}648$ residues by `c3a_hist3.c` (Lemma 0's carry automaton; equal to a plain triple loop and to a scatter of $D^3$ on seven small cases), byte-identical to the package's; $N_0 = 5794476231064238$ and $C = 39290318538006294714264853344$ exactly as stated, every $N_r \ge 1$; $\theta \ge 1.19359680912973029490196359640645\ldots$, every bound agreeing with the note to the digits it prints | 1.19359680 |

Red controls, all as expected: the claim raised to 1.19281556 is not shown; dropping the generator 2800 changes the
data; $x$ one unit lower in its last bit changes $N_0$ and $C$.

What rests on mathematics: Lemmas A-D of the note (the encoding, injectivity, the sumset bound and [GHR2007] with
$V = 2U$), Kleinwaks's [K2026b] Proposition 3, Lemma 0 of the general-$k$ note (the carry automaton), and Robbins'
bounds. The $k = 3$ histogram takes about 6 CPU hours (`python3 -I c3a_check.py syl3 certs --controls`, 64 row
parts; `--parts=DIR` reads the outputs of an earlier run). Its red controls: the claim raised to 1.19359681 is not
shown, $x$ one unit lower changes $N_0$ and $C$, and a histogram with one residue moved fails the symmetry check
(only the classes of 0 and $h$ may be odd).

[M2026] (merged in [#236](https://github.com/teorth/optimizationproblems/pull/236)), rebuilt from Theorems 1-4 of its
package's README alone (`c3a_check.py m2026|m2026b|k2026b`): $M = \langle G \rangle \cap [0, B]$ by reachability;
$\kappa(d) = \min\{a + b : a, b \in M, a - b = d\}$ over all pairs; one representative $d(r)$ for each complementary
residue pair $\{r, q - r\}$ (the cheaper one in $M - M$, $r$ on a tie), with the check that $0$ and the $\pm d(r)$ are
pairwise distinct mod $q$; the counts $n_d = \lfloor N x^{\kappa(d)} / Z(x) \rfloor$ with $Z(x) = 1 + 2\sum x^{\kappa(d)}$,
$n_0 = N - 2\sum n_d$ and $L = \sum n_d \kappa(d)$; then $\log D_N = \ln N! - \ln n_0! - 2\sum \ln n_d!$ from below
(Robbins), $\log S_N(x) = N \ln P_+(x) - 2L \ln x$ from above and $\log(5q^N)$ ($\log(2Q^N)$ in the carry-free base)
from above, in 90-digit decimal interval arithmetic. The package's data files are only compared.

| Row | Data | Our result | Printed bound |
| --- | --- | --- | --- |
| 3a [M2026], $B = 17086$, $q = 27995$, $N = 10^{13}$ | mizopl/ghr-c3a-odd-base-finite-bound @464ceda `object/` (`b1230dd4...1fe21` and the mask, sum-support and pruned-difference files, hashes as in its README) | $M$, $M + M$ and all 13,997 rows of the pruned table ($r$, $d$, $\kappa(d)$, both residues, $n_d$) equal ours; a complete residue system mod $q$; $n_0$ and $L$ as stated; $\theta \ge 1.19230254070568368329\ldots$ | 1.19230254 (row comment 1.1923025406835) |
| 3a [M2026] on the [K2026] mask, $B = 17032$, $q = 27995$ | `finite_certificate_17032_27995.json` (`db1517d0...ad04b`) only compared | $n_0$ and $L$ as stated; $\theta \ge 1.19229904840497896133\ldots$ | 1.1922990483828 |
| control: [K2026b]'s finite certificate, $B = 17032$, $Q = 2B + 1 = 34065$ | `kleinwaks_finite_certificate_17032_34065.json` (`5c85a63a...8d9b`) only compared | $n_0$ and $L$ as stated; $\theta \ge 1.19102809427254087548\ldots$ | 1.19102809 |

Red controls, all as expected (the same check functions on altered data): both members of a residue pair kept in the
alphabet are refused; a witness shifted to cost $\kappa(d) + 2$ is refused; counts raised by 2% (so $n_0 < 0$) are
refused; a claim $10^{-5}$ above the printed one is not implied. What rests on mathematics: Theorems 1-4 of the
package's README (base-$q$ injectivity on words with pairwise distinct residues, the multinomial count, the sumset bound
$\lvert U + U\rvert \le P_+(x)^N / x^{2L}$, and the [GHR2007] lemma for $V = 2U$ with $2\max V + 1 < 5q^N$), and
Robbins' bounds.

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

## C3e: a sums-and-differences constant (6 November 2026)

`c3e/c3e_check.py`, standard library only: mizopl's finite carry-free lower bounds $C_{3e} > 1.27255802$ ($B = 160$) and
$1.27255796$ ($B = 80$) from Theorems C1-C3 of their certificate folder's README alone (posted on
[#339](https://github.com/teorth/optimizationproblems/pull/339)). The counts $n_d$, $n_0$ and $L$ in exact integers,
$\log D_N$ from Robbins' bounds and $\log S_N(y)$ in decimal interval arithmetic at 90 digits: $1.27255802020997507810\ldots$
and $1.27255796845856979212\ldots$. Theorems C1-C3 are also checked by brute force on small cases ($B = 15$), with five
red controls. Reproduce: `cd c3e && ./fetch.sh && ./run_all.sh`.

## C59: the Bohr radius of the bidisc (6 November 2026)

`c59/c59_check.py`, standard library only, integers only: the certificate of
[#297](https://github.com/teorth/optimizationproblems/pull/297), $K_2 < 0.302825279490988$. $U^*U = UU^* = I$ exactly
for the $3 \times 3$ Gaussian-rational colligation, so the transfer function is in the Schur class; its 1,681
coefficients $f_{jk}$ ($j, k \le 40$) as exact Gaussian integers, checked against the closed form at two points; and the
finite Bohr sum minus 1 in $[2.2259 \cdot 10^{-42}, 2.2259 \cdot 10^{-42}]$ by integer square roots. Red controls: $r$
lowered to 0.30282527949097 (the sum falls below 1), one matrix entry changed by one unit (refused by unitarity), the
sum cut at $N = 20$ (below 1). Reproduce: `cd c59 && ./fetch.sh && ./run_all.sh` (about a minute).

## C42: Turán's pure power-sum constant (6 November 2026)

`c42/c42_check.py`, standard library only: the two rows whose first-block cut-off is $\tau \ge 1/3$, Griego's
$C_{42} \le 0.6906538$ [Gri26] (added in [#85](https://github.com/teorth/optimizationproblems/pull/85)) and Pinho's
$0.6888863712$ [Pin26] ([#251](https://github.com/teorth/optimizationproblems/pull/251)), from the criterion as Griego's
note (section 4) and Pinho's paper (section 2) state it: with $K = I_\alpha(\tau)$, $D = I_a(\tau)$,
$I_A(x) = \int_0^x v^{A-1}/(1-v)\,dv$ and the pair integrals $J(t, r)$,
$Y = 1 + sK - \sum_j \delta_j (I_\alpha(1 - t_j) - K) + \tfrac12 \sum_{i,j} \delta_i \delta_j J(t_i, t_j)$, and the bound
holds if $|s| < C$, every $|\eta_j| < C$ and $|Y| < CD$. The radius conditions in exact rationals; $I$ and $J$ by their
power series with explicit tails, in decimal interval arithmetic at 50 digits ($x^\alpha$ through $e^{a \ln x}$ and
Taylor series for sine and cosine with bounded remainders).

| Row | Certificate (SHA-256) | Our result | Printed bound |
| --- | --- | --- | --- |
| [Gri26] | sebastian-griego/turan-c42-certificate @ `c8ddce14` (v1.0.0), `certificate/turan42_certificate.json` (`a76e9132...3a81`, as its .sha256 file states) | $\lvert Y\rvert/D \le 0.69065369515163118444\ldots$, $C^2D^2 - \lvert Y\rvert^2 \ge 1.550 \cdot 10^{-7}$; our $K$, $D$ and $Y$ overlap the certificate's printed enclosures | 0.6906538 |
| [Pin26] | rpinho/c42-certificate @ `a8a9bc30`, `certificates/c42_256.json` (`705a7291...afcab`, as the page states; 256 blocks, 16,512 pairs) | $\lvert Y\rvert/D \le 0.68888637098346564582276780\ldots$, $C^2D^2 - \lvert Y\rvert^2 \ge 2.3747 \cdot 10^{-10}$, both as the paper prints; $D$ to 30 digits as printed | 0.6888863712 |

Red controls, for each row: $C$ lowered just below our $\lvert Y\rvert/D$ (refused by the radius condition, which binds as
tightly), $\eta_0$ replaced by $s$ (refused by $\lvert Y\rvert < CD$), and the cut-off moved to $\tau = 3/10$ (refused: three
corrections would then reach $b_n$). A floating-point look at the construction itself (`--finite`, not part of the
proof): $\lvert Y_n\rvert / \sum_F \lvert P_m\rvert$ falls toward $\lvert Y\rvert/D$ like $n^{-a}$ ($a = \mathrm{Re}\,\alpha$), from
0.70436 at $n = 500$ to 0.69447 at $n = 4000$ for [Gri26], and from 0.70268 to 0.69255 for [Pin26]. What rests on
mathematics: the asymptotic reduction (the limits of the normalized sums), which the rows mark as unrefereed.

`c42/kos_check.py`, standard library only: Koshta's $C_{42} \le 0.6873358063$ [Kos2026]
([#263](https://github.com/teorth/optimizationproblems/pull/263)), with first-block cut-off $\tau = 41/200$. It rests on
the theorem in section 4 of his note, which extends the criterion to any $0 < \tau < 1/2$: $Y$ gains the terms
$\frac{(-1)^d}{d!} \sum \delta_{j_1} \cdots \delta_{j_d} T_d(t_{j_1}, \ldots, t_{j_d})$ over ordered tuples for
$2 \le d \le d^* = \lceil 1/\tau \rceil - 1$ (here $d^* = 4$), where
$T_d(u) = \int_{x \ge u,\ \sum x_i \le 1} (1 - \sum x_i)^{\alpha - 1} / (x_1 \cdots x_d)\,dx$ and $T_2 = J$. The note's
checker, verifier and cross-check code were neither read nor run, and $T_d$ comes from a series of our own, not from the
note's sections 5 and 7. On the slice $\sum x_i = S$, $1/(x_1 \cdots x_d) = S^{-1} \sum_k \prod_{i \ne k} 1/x_i$. So the
integral of $1/(x_1 \cdots x_d)$ over a slice is a sum of volume integrals for the tuple's sub-tuples, and each volume
integral is the integral of the sub-tuple's own slice integral. Expanded about the face $\sum x_i = 1$, all these series
converge with ratio $L/(L + \min u)$, where $L = 1 - \sum u$. That ratio is below 1 for every tuple: at most 0.742, 0.653
and 0.468 for $T_2$, $T_3$ and $T_4$ here, including the 250 triples with $L \ge \min u$, where the note's series about
the vertex diverges. A majorant carried through the same recursion bounds every tail. The arithmetic is
`c42_check.py`'s, with the angle in $x^{ib}$ reduced by an enclosure of $2\pi$ from Machin's formula.

| Row | Certificate (SHA-256) | Our result | Printed bound |
| --- | --- | --- | --- |
| [Kos2026] | DennisKoshta/c42-turan-power-sum-certificate @ `8e1797e0`, `certificates/c42_tau41_200_M64_0.6873358063.json` (`be0c406f...cbef`, as the note states; 64 blocks; 1,056 pairs, 2,443 triples and 609 quadruples with sum below 1) | $\lvert Y\rvert/D = 0.687335806186516616249961526259099396\ldots$ (enclosure width $1.4 \cdot 10^{-44}$), $C^2D^2 - \lvert Y\rvert^2 = 4.348178358021 \cdot 10^{-11}$; the note's 30 decimals of $\lvert Y\rvert/D$, its margin $4.348178358 \cdot 10^{-11}$ and its four-fold sum $1.214 \cdot 10^{-5} - 3.45 \cdot 10^{-6} i$ all agree | 0.6873358063 |

Controls (`--controls`), all as expected. The same code at $\tau \ge 1/3$ ($d^* = 2$) reproduces both rows above, with
$Y$ within $1.2 \cdot 10^{-48}$ of `c42_check.py`'s. $C$ lowered to 0.687335806186515 and $\eta_0$ replaced by $s$ are both
refused by $\lvert Y\rvert < CD$. $T_2$ on all 1,056 pairs overlaps `c42_check.py`'s $J$. On ten tuples with
$L < \min u$ (up to $L/\min u = 0.878$), $T_3$ and $T_4$ overlap the note's series about the vertex, which we also
evaluate rigorously. $T_3$ with $L \ge \min u$ (up to 1.878) agrees with $\int_0^L H_2(\sigma) K(L - \sigma)\,d\sigma$ to
$10^{-45}$, and $T_4$ with the note's scaling identity to $10^{-49}$, both by tanh-sinh quadrature (not rigorous). The
four-fold term raises $\lvert Y\rvert/D$ by $8.5 \cdot 10^{-7}$, about 7,500 times the gap $C - \lvert Y\rvert/D$, so it
cannot be left out. What rests on mathematics: the general-cut-off reduction of section 4, which the row marks as
unrefereed; this replay checks the limiting inequality.

A fix in `c42_check.py` (6 November): `abs()` and unary minus on a decimal round in Python's default context (28
digits, half-even), not outward. `I.mag` and `I.__neg__` now use the exact `copy_abs` and `copy_negate`. The effect was
below $10^{-28}$ relative, and the outputs for [Gri26] and [Pin26] are unchanged in every printed digit.

Reproduce: `cd c42 && ./fetch.sh && ./run_all.sh` (all three rows with their controls, about five minutes; `run_all.sh`
then also runs `thmA_check.py`, below).

## C42: three bounds on Bakeas's Theorem A (7 November 2026)

Three pending rows rest on one limit theorem: Bakeas's $C_{42} \le 0.68637$ [Bak26]
([#302](https://github.com/teorth/optimizationproblems/pull/302)), Pinho's $0.6863607950$ [Pin26b]
([#305](https://github.com/teorth/optimizationproblems/pull/305)) and Vatopoulos's $0.686360751$ [Vat26]
([#327](https://github.com/teorth/optimizationproblems/pull/327)). With $T > 0$, a phase $\varphi$ piecewise linear
through knots $0 = u_0 < \dots < u_K = \tfrac12$, $g = Te^{i\varphi}$, $\gamma = g(0)$, $a = \mathrm{Re}\,\gamma$, and
$d = \gamma - g$ on $(0, \tfrac12]$, $d = \gamma$ on $(\tfrac12, 1]$, Theorem A of Bakeas's `NOTES.md` (v1.1.0) says: if
$\kappa = \sup \lvert d\rvert < 1$ and $a \ge 0$, then $\limsup R_n \le \max(T, F)$, where
$F = \lvert p(1)\rvert / \int_0^{1/2} \lvert p(t)\rvert\,dt/(1-t)$, $p(t) = t^{-\gamma} + \int_0^t (t-v)^{-\gamma}\rho(v)\,dv$
and $\rho$ is the bounded solution of $t\rho(t) = d(t) + \int_0^t d(u)\rho(t-u)\,du$. Each row takes $T = C$ and needs
$F < C$. ([#310](https://github.com/teorth/optimizationproblems/pull/310), the explicit lower bound, does not bear on
these rows.)

**The proof.** We re-derived every step of Lemma 1, Lemma 2 (Harcos's linearization), Theorem A (steps 1 to 5) and
Lemma 3 with (8), and found no gap. In Lemma 2, $\sigma_k = \lambda u_k$ gives $\sum_{k>K}\sigma_k P_{n-k}/k = P_n$,
$s_k(z) = \sigma_k$ for $k \le n$ (logarithms modulo $x^{n+1}$), and $2(K+1) > n$ kills every product of two high terms,
so $W(1) = 0$. In the theorem: the exact identity $\sum P_m x^m = (1-x)^{\gamma - 1}\exp(\sum_k d(k/n)x^k/k)$ (also for
$k > K$, since $k/n > \tfrac12$ there); the contraction $\lVert\mathcal K\rVert \le \kappa$, $\lvert d(u)\rvert \le T\Lambda u$ and the
regularity of $\rho$ (one jump, at $\tfrac12$); the discrete recursion, its uniform bound $B = M/(1-\kappa)$ and the
error recursion, whose quadrature error is $O(\omega_d + \omega_\rho)$ plus $O(1/n)$ from at most two cells; the
splitting of $\Pi_n(j)$, the domination (B2) and the convergence (B1); dominated convergence for the denominator. The
proof uses the knots only through $\lvert d(u)\rvert \le T\Lambda u$ and the continuity of $d$ and $d(t)/t$ on $(0, \tfrac12]$,
so it holds verbatim for the knots $j^2/4608$ of #305 and the 34 knots of #327, as #305 says. The hypothesis
$a \ge 0$ enters only through the monotonicity of $v^{-a}$ in the Riemann-sum bounds of step 4. For a constant phase,
$\rho = 0$ on $(0,\tfrac12]$ and $\rho = \gamma/t$ beyond, which gives Harcos's formula
$F = \lvert 1 + \gamma\int_0^{1/2} t^{-\gamma}dt/(1-t)\rvert / \int_0^{1/2} t^{-a}dt/(1-t)$; at $\gamma = 0.43246 + 0.54237i$ it
is $0.6936763$ (rigorous series) with $T = \lvert\gamma\rvert = 0.6936763$. Lemma 3:
$e = \rho - \tilde\rho$ satisfies $e = R/t + \mathcal K e$, so $\lVert e\rVert \le \lVert R/t\rVert/(1-\kappa)$, and
$\lvert p - \tilde p\rvert \le \lVert e\rVert\,t^{1-a}/(1-a)$ gives (8). We also re-derived the integrated (Gronwall) bound of #305
(correct; not used here). This check was made with AI assistance and is not a referee report.

**The numbers.** `c42/thmA_check.py`, standard library only (with `c42_check.py`'s decimal intervals):
- Hypotheses, rigorously: $0 < T \le C < 1$, $a = T\cos\psi_0 > 0$, and $\kappa = \max(T, 2T\sin(E/2)) < 1$ with
  $E = \max\lvert\psi_j - \psi_0\rvert$ ($\lvert d\rvert = 2T\lvert\sin((\varphi-\psi_0)/2)\rvert$ on $(0,\tfrac12]$).
- A trial $\tilde\rho$, a discontinuous polynomial of degree $P$ on each cell of the grid $h = 1/N$ (which contains
  every knot, so every sum of knots too): floating-point collocation, then exact dyadic coefficients.
- Its residual $R = d + d*\tilde\rho - t\tilde\rho$, exactly per cell: on cell $n$,
  $(f*g) = h\sum_c \frac{\theta^{c+1}}{(c+1)!}U_{n,c} + h\sum_c\frac{(1-\theta)^{c+1}}{(c+1)!}V_{n-1,c}$, where $U$ and $V$ are
  two-dimensional convolutions (cell, degree) of $a!\,f_{i,a}$ and $b!\,g_{j,b}$ and of the same for the reversed
  polynomials. Each is one product of decimal integers (the arrays packed by Kronecker substitution, every rounding
  trapped). $d$ is replaced by its Taylor polynomial of degree $Q$ on each cell, remainder bounded; $\sup\lvert R/t\rvert$
  on each cell from the moduli of the Chebyshev coefficients, computed exactly.
- From the residual to $e = \rho - \tilde\rho$: $\lvert e\rvert \le r + \lvert\mathcal K\rvert\,\lvert e\rvert$ with
  $\lvert\mathcal K\rvert f(t) = t^{-1}\int_0^t\lvert d(u)\rvert f(t-u)\,du$. Any bounded $E \ge 0$ with
  $E \ge r + \lvert\mathcal K\rvert E$ bounds $\lvert e\rvert$ (for $w = \lvert e\rvert - E$, $w^+ \le \lvert\mathcal K\rvert w^+$, and
  $\lVert\lvert\mathcal K\rvert\rVert \le \kappa < 1$). We build $E$ on blocks of cells by a forward recursion, and
  $EN = \int_0^1(1-v)^{-a}E \ge \lvert p(1) - \tilde p(1)\rvert$,
  $ED = \frac{2}{1-a}\int_0^{1/2}(\tfrac12 - v)^{1-a}E \ge \int_0^{1/2}\lvert p - \tilde p\rvert\,dt/(1-t)$.
- $\tilde p(1) = 1 + h^{1-\gamma}\sum_c\int_0^1(c-y)^{-\gamma}\tilde\rho_{N-c}(y)\,dy$, by
  $(c-y)^{-\gamma} = c^{-\gamma}\sum_j\frac{(\gamma)_j}{j!}(y/c)^j$ ($\lvert(\gamma)_j/j!\rvert \le 1$) and $B(k+1, 1-\gamma)$ for $c = 1$.
- The denominator from below. On cell $n \ge 1$, $\tilde p = h^{-\gamma}(A(\xi) + h\,x^{1-\gamma}J(x))$, $x = \xi + \tfrac12$:
  $A$ is a polynomial from the series of $(n + \tfrac12 + \xi)^{-\gamma}$, of
  $\sum_k\tilde r_{n-1,k}B(k+1,1-\gamma)(\tfrac32 + \xi)^{k+1-\gamma}$ (the previous cell's polynomial continued over cell
  $n$) and, for the cells two or more back, of a second exact two-dimensional convolution; $J = \sum_k\Delta_kB(k+1,1-\gamma)x^k$
  holds the jumps $\Delta_k$ of $\tilde\rho$'s Taylor coefficients at the node. Then
  $\sqrt X \ge \tfrac32Xu - \tfrac12X^2u^3$ for $X, u \ge 0$ (the difference is $s(s-1)^2(s+2)/(2u)$, $s = u\sqrt X$), with
  $X = \lvert A\rvert^2$ exact and $u \approx 1/\lvert A\rvert$; the jump term to first order,
  $\lvert A + S\rvert \ge \lvert A\rvert + \mathrm{Re}(S\bar A)u - \eta\lvert S\rvert$ ($\eta \ge \sup\lvert u^2X - 1\rvert$), integrated
  with the moments $1/(k + 2 - \gamma)$; at the cells that start at a knot, the same inequality with $X = \lvert A + S\rvert^2$
  in generalized monomials $x^{k + m(1-\gamma) + m'(1-\bar\gamma)}$; cell 0 exactly ($\tilde p = h^{-\gamma}x^{-\gamma}q_0(x)$).
- $F \le (\lvert\tilde p(1)\rvert + EN)/(\tilde D - ED)$, compared exactly with $C$.

| Row | Profile (SHA-256) | Grid, $\sup\lvert R/t\rvert$ | Certified $F \le$ (diagnostic $F \ge$) | Claim $C = T$ |
| --- | --- | --- | --- | --- |
| [Bak26] ([#302](https://github.com/teorth/optimizationproblems/pull/302)) | Bakeas/turan-c42-bound @ `dc9c491c` (v1.1.0), `c42_profile_0.68637.json` (`dcea4edb...d395`, as #302 and the note state; 25 knots $j/48$; $a = 0.1660999$, $E = 0.7753$) | $N = 384$, $P = 10$; $4.8 \cdot 10^{-12}$ | 0.6863600829166 (0.6863600828267) | 0.68637: $C - F \ge 9.9 \cdot 10^{-6}$ |
| [Pin26b] ([#305](https://github.com/teorth/optimizationproblems/pull/305)) | rpinho/c42-certificate @ `009120a8`, `volterra-quadratic-knots/certificates/C42_compact_best.json` (`ac8d052d...5765`, as its README states; 49 knots $j^2/4608$; $a = 0.1118266$, $E = 0.8561$) | $N = 9216$, $P = 10$; $2.4 \cdot 10^{-9}$ | 0.6863607440170 (0.6863607438702) | 0.6863607950: $C - F \ge 5.09 \cdot 10^{-8}$ |
| [Vat26] ([#327](https://github.com/teorth/optimizationproblems/pull/327)) | NikosVat/turan-c42 @ `9668a60e`, `certificate/c42_profile_0.686360751.json` (`d015a829...5652`, as #327 states; 34 knots, least common denominator 18432; $a = 0.1016540$, $E = 0.8709$) | $N = 18432$, $P = 14$; $2.1 \cdot 10^{-9}$ | 0.6863607495095 (0.6863607477240) | 0.686360751: $C - F \ge 1.49 \cdot 10^{-9}$ |

In every row $E < \pi/3$, so $\kappa = T$. The error terms are $EN \le 7.6 \cdot 10^{-14}$, $1.8 \cdot 10^{-11}$ and
$2.0 \cdot 10^{-11}$ and $ED \le 4.0 \cdot 10^{-14}$, $9.2 \cdot 10^{-12}$ and $1.0 \cdot 10^{-11}$; the cells treated with
the kink kept exact number 24, 53 and 43. Floating point, not part of the proof: the construction itself gives
$F_n = 0.6908134$, $0.6878859$, $0.6868722$ ([Bak26]), $0.6908334$, $0.6878780$, $0.6868586$ ([Pin26b]) and $0.6908328$,
$0.6878775$, $0.6868593$ ([Vat26]) at $n$ = 1000, 4000, 16000, falling slowly toward the limits.

The three rows' own numbers, read after our runs: [Bak26]'s two Arb checkers $F \le 0.6863660$ and $0.6863615$
($N = 122880$) and its float solve 0.68636008; [Pin26b]'s verifiers $F \le 0.6863607743$ and $0.6863607774$
($N = 14155776$) and Bakeas's float solve 0.686360744; [Vat26]'s enclosure $[0.6863607494206, 0.6863607494291]$.
Our upper bounds are consistent with all of these: below the certified bounds of [Bak26] and [Pin26b] and at their
float values, and $8 \cdot 10^{-11}$ above the top of [Vat26]'s enclosure (our $EN$ and $ED$). For [Vat26] the margin is
resolution-limited on our side: at $P = 12$ on the same grid our upper bound is 0.6863607522834, above the claim,
because the residual on the cell after the knot $1/18432$ is then $7.3 \cdot 10^{-8}$. (#305's README prints its two
bounds as 0.6863607742 and 0.6863607773, its PR text as 0.6863607743 and 0.6863607774.)

Self-tests (`--selftest`), all passing: the exact convolution against direct summation, identical on six random
signed arrays with entries up to $10^{60}$; Harcos's constant phase ($T = 0.69367634$, $\psi = 0.89767188$), whose
closed form $0.69367630664\ldots$ lies inside our enclosure at $N = 16$ and 64, below its top by $2.0 \cdot 10^{-13}$ and
$8.0 \cdot 10^{-14}$; and convergence: the [Bak26] upper bound is 0.6863600829238, 0.6863600829167, 0.6863600829166 and
0.6863600829167 at $N$ = 96, 192, 384 and 768, [Pin26b] at $N = 4608$, $P = 12$ gives 0.6863607446243, and [Vat26] at
$P = 12$ gives 0.6863607522834 (above).

Red controls (`--controls`), all as expected. For each row, $T = C = 0.68636$ (at or below $F$) and one phase raised
([Bak26] $\psi_0 + 0.05$, [Pin26b] $\psi_3 + 0.02$, [Vat26] $\psi_{20} + 0.005$) all FAIL, and in every case our
diagnostic lower bound on $F$ exceeds $C$ (by $1.2 \cdot 10^{-5}$, $2.6 \cdot 10^{-5}$, $1.7 \cdot 10^{-6}$,
$9.4 \cdot 10^{-7}$, $1.7 \cdot 10^{-6}$ and $2.7 \cdot 10^{-7}$), so these claims are false, not only unproved. Refused:
[Bak26]'s phases all lowered by 0.3 ($\psi_0 = -1.626$, $\mathrm{Re}\,\gamma = -0.038 < 0$, which Theorem A excludes),
$T = C = 1$ ($\kappa = 1$), and [Vat26] on a grid that misses its knot $1/18432$.

Not read or run: any checker, solver or notebook of the three repositories. Read: Bakeas's `NOTES.md`, Pinho's
`README.md`, `BAKEAS_ADAPTATION.md` and `VOLTERRA_VERIFIER_MATH.md`, Vatopoulos's `NOTES_CERTIFIER.md` (and the
headings and opening claims of his `README.md`, `NOTES.md` and `MATH.md`), the three PR texts and the profiles. Two of
our devices also appear in Vatopoulos's
`NOTES_CERTIFIER.md`, the $U$/$V$ identity and the inequality $\sqrt X \ge \tfrac32Xu - \tfrac12X^2u^3$; we re-derived
both. The rest differs from all three: exact decimal products instead of ball arithmetic, our grids and degrees, the
comparison majorant instead of Lemma 3's sup bound or the Gronwall bound, the jump terms and the kink cells. What rests
on mathematics: Theorem A itself (unrefereed; our check above).

Reproduce: `cd c42 && ./fetch.sh && python3 -I thmA_check.py --selftest --controls --write` (the three rows take about
6 minutes, and about 25 with the self-tests and controls, on a loaded Apple machine; `./run_all.sh` runs it after the
checks above).

## C12a: the Euclidean travelling salesman constant (6 November 2026)

`c12a/c12a_check.c`, C11 with 128-bit integers: the computation behind $\beta \le 0.8618$ [O2026]
([#222](https://github.com/teorth/optimizationproblems/pull/222)), written from `RECIPE.md` of
devinokeefe/travelling-sales-person-constant-bounds @ `69b864b9` alone (its three checkers were neither read nor run).
Every integer of Sections 2 to 8 is recomputed from the parameters of a spec file and compared with its values file
(the masses, the tail term, every table line, every cell sum and the total), then the test
$A_0 h_2 \theta_2 \le \theta_1 n^{k+1} F k h_1$. The cell sums (Section 7; for $k = 9$, $4^9 \cdot 6^{10}$ cells, each
the least cost over $8!$ orders) use the shared-paths identity, with $H$ built in nested loops over the middle heights
(Held-Karp path costs over subsets, then the two-path covers), the reversal and the reflection on $a_1$ only. Without
any of the three identities, all cell sums of `tiny` to `k6` were also recomputed by all orders in every cell, and seven
each of `k7`, `k8` and `k9` by Held-Karp in every cell.

| Spec | Numbers compared | Total $A_0$ | Bound (rounded up) | Threshold |
| --- | --- | --- | --- | --- |
| `k5.txt` | 2,280, all equal | 1761033937996622 | 0.878820870371 | 8789/10000 |
| `k6.txt` | 6,192, all equal | 12593568504391125 | 0.872868757138 | 8729/10000 |
| `k7.txt` | 19,680, all equal | 92054418307910785 | 0.868075379727 | 8681/10000 |
| `k8.txt` | 70,482, all equal | 628708228549677266 | 0.864607001369 | 8647/10000 |
| `k9.txt` | 269,290, all equal | 4430997795015293889 | 0.861715646305 | 8618/10000 |

The test cases `tiny`, `rounding`, `bhh` and `k4` (146, 971, 8 and 952 numbers) agree as well. Red controls: one cell
sum of `k5` increased by 1 (refused; that number alone differs), the `k7` threshold set $10^{-12}$ below its rounded-up
bound (refused by the test; accepted at the bound), and mesh 127 in place of 128 for `k5` (refused; 1,779 of the 2,280
numbers differ). Not checked here: that `RECIPE.md` defines the integers of Section 4.1 of the paper, and the proof that
this number bounds $\beta$. Reproduce: `cd c12a && ./run_all.sh` (a C11 compiler; `k9` took 73 minutes on 12 threads of a
loaded Apple M5 Pro, and `./run_all.sh quick` leaves it out and takes about 7 minutes).

## C47: the centred maximal constant over squares in the plane (7 November 2026)

The weak type (1,1) constant $c_2$ of the centred maximal operator over axis-parallel squares. A lower bound is
$\lambda \lvert U \cap C\rvert/\mu(C)$ for a measure $\mu$ with finitely many atoms in each cell $C$ of the lattice
$P_x\mathbb Z \times P_y\mathbb Z$ (Tomar's Lemma 6.2, after Lin): $U$ is the union of the lattice translates of
rectangles $R_i \subset \bigcap_{p \in S_i} Q(p, \rho_i)$, each with a finite set $S_i$ of atoms and
$0 < \rho_i < \sqrt{\mu(S_i)/(4\lambda)}$. We re-derived the lemma (truncate $\mu$ to $N \times N$ cells, spread each atom
over $Q(p, \varepsilon)$, apply the weak type inequality to that function, let $N \to \infty$) and checked each step of
its proof in the paper.

`c47/c47_check.py`, standard library only, exact integers throughout: the rectangles of Tomar's Section 6. For every
box whose sides lie on atom coordinates, with length at most the file's `reach` (1396/125), $S$ is the set of atoms in
the box and $\rho_S = (k - 1)\,2^{-40}$, $k$ the largest integer with $(k\,2^{-40})^2 \le \mu(S)/4$. Every box with its
lower-left corner in $C$ is enumerated, with masses from 2D prefix sums over the periodic grid, and each rectangle is cut
into its pieces inside $C$. Two different methods then fix the area:
- subtraction: the pieces are removed from $C$ one by one, which leaves the uncovered pieces;
- a slab sweep: $\lvert\bigcup K\rvert$ for the pieces $K$ that removed area.

Every uncovered piece is checked against every piece of the family, so the two give $\lvert U \cap C\rvert$ from both
sides. For every rectangle in $K$, $\mu(S)$ is summed again directly from the atoms, and $0 < \rho$,
$4\lambda\rho^2 < \mu(S)$ and the containment are checked. The definition is also applied point by point: atoms sorted
by sup-distance $d$, covered when $\rho(\mu(\text{atoms within } d)) \ge d$. No uncovered piece has its centre covered,
and union and definition agree at 200 random points. With `--full`, every box up to the length $L^*$ is enumerated:
beyond $L^*$, $\mu(\text{box}) \le (L/P_x + 2)(L/P_y + 2)\mu(C) \le L^2$, so no rectangle has interior. This is run
with Tomar's radii and with radii on $2^{-100}$ below and above $\sqrt{\mu(S)/4}$. The second pair encloses
$\lvert\{M\mu \ge 1\} \cap C\rvert$ itself.

| Row | Data (SHA-256) | Exact ratio | Printed lower bound |
| --- | --- | --- | --- |
| [Tom2026], Theorem 1.3 (row added by [#321](https://github.com/teorth/optimizationproblems/pull/321)) | tejasvi/square-maximal-weak-type @ `3af09fc` (v5; the same file is in Zenodo 10.5281/zenodo.23716078), `lowerbound/c2lower_1p70272.json` (`e2ca1616...9888`, 50 atoms) | 482167717377842991882283179477574801392641890014291050323/283174368139841919615520431313000857449005056000000000000 $= 1.70272373359629370268\ldots$, the 16 decimals #321 states; the level set itself gives $1.70272373360313336487\ldots$ | 1.7027237335 |
| [#321](https://github.com/teorth/optimizationproblems/pull/321) | rpinho/math-certificates @ `3d4cd55`, `c47/round1/certificates/half59-measure.json` (`64585e1f...4b16`, 59 atoms) and `half59.json` (`b4d88847...ba99`) | identical to the stated fraction, $1.70311634551278326393\ldots$; the certificate's 559 witnesses, each checked, give the same union; the level set itself gives $1.70311634551668029559\ldots$ | 1.70311 (the row), 1.70311634 (#341's table) |
| [#341](https://github.com/teorth/optimizationproblems/pull/341) | `construction.json` (`5f7a8d9d...310b`) in Zenodo 10.5281/zenodo.24022791, which is restricted (HTTP 403): not replayed | stated: 3132581668909806427917067970204393210052649/1839223215401194968772738622942534413844480 $= 1.70320907363410346949\ldots$ | 1.7032090736 (implied by the stated fraction) |
| control: Lin's lattice [Lin2026] | $u = (2 + \sqrt{22})/3$ to 60 digits, radii on $2^{-300}$, every box up to $L^*$ | $\Phi$ to $1.3 \cdot 10^{-62}$; the uncovered set is exactly the four $a \times b$ slots, with $M\mu < 1$ at their centres | 1.6855099933 |

For both replayed rows the reach does not bind: every box up to $L^*$ (19.58 and 19.61) gives the same union. The
$2^{-40}$ radii give up $6.8 \cdot 10^{-12}$ and $3.9 \cdot 10^{-12}$ of the level set's ratio. Our 33,356 pieces in $C$
for #321 are the number of rectangles #321 reports for Tomar's checker.

Self-tests, all passing:
- $\Phi$ from the 47a page's closed form and from Lin's README agree, round down to 1.6855099933 and match Lin's
  1.68550999335552518466528;
- a single atom per cell gives exactly $(1 - 2^{-39})^2$;
- Tomar's measure without its 48 small atoms gives 1.6447274428 ("about 1.6447", his Section 6);
- on eight small random measures (two with a reach shorter than $L^*$, radii on $2^{-12}$), brute force agrees exactly
  with the fast method. Every cell of the grid of all candidate edges is decided at its centre by the pointwise
  definition, and boxes are summed directly.

Red controls, all as expected, for each replayed row: one atom moved by $10^{-6}$, the heaviest mass raised by $10^{-12}$
and $P_x$ raised by $10^{-9}$ each change the fraction; a bound one unit higher in the last printed digit is not implied;
one radius raised to $(k+1)\,2^{-40}$ is refused by $4\lambda\rho^2 < \mu(S)$. What rests on mathematics: Lemma 6.2 (re-derived),
and, for the level-set values only, that the union of the box rectangles with exact radii is $\{M\mu \ge 1\}$ (the
supremum in $M\mu$ is reached at the sup-distance of an atom). Reproduce: `cd c47 && ./fetch.sh && ./run_all.sh` (about
90 seconds). If the #341 file becomes available, put it at `certs/pr341/construction.json`; its SHA-256 is checked first.

## C42: Röhrig's eight-step profile (7 November 2026)

`c42/c42_check.py` (the code of the C42 section above) now also replays Röhrig's $C_{42} \le 0.688983$\* [R2026a], the C42
row of [#184](https://github.com/teorth/optimizationproblems/pull/184), from the criterion as his `notes/turan.md` states it
(Griego's, with eight middle values on the breaks $t_j = (8+j)/24$, $\tau = 1/3$) and `certificates/turan.json` alone
(`certify_turan.py` and `audit.py` were neither read nor run).

| Row | Profile (SHA-256) | Our result | Printed bound |
| --- | --- | --- | --- |
| [R2026a] | AndreRoehrig/optimization-constants-certificates @ `ed638739` (v1.0.2), `certificates/turan.json` (`18c20906...30f4`, as its MANIFEST states) | $\lvert s\rvert, \lvert\eta_j\rvert < C$ exactly (least $C^2 - \lvert\cdot\rvert^2 = 1.242\cdot10^{-6}$); $D \in [0.888190650452511954178223193231, \ldots232]$, $\lvert Y\rvert/D \le 0.68898209850406074850512099758364\ldots$, $C^2D^2 - \lvert Y\rvert^2 = 9.7997362463618\cdot10^{-7}$: all three 30-decimal enclosures of the note agree | 0.688983 |

Red controls as for the other rows (refused: $C$ just below our $\lvert Y\rvert/D$, $\eta_0$ replaced by $s$, the cut-off moved
to 3/10). Floating point, not part of the proof: $\lvert Y_n\rvert/\sum_F\lvert P_m\rvert$ = 0.70271, 0.69781, 0.69465, 0.69262 at
$n$ = 500 to 4000. The displayed enclosures of `c42_check.py` are now rounded outward (they were rounded to nearest; the
checks themselves were always exact). Reproduce: `cd c42 && ./fetch.sh && python3 -I c42_check.py certs --rows=roe --controls`.

## C51: the Erdős maximum-term constant (7 November 2026)

`c51/c51_check.py`, standard library only: Röhrig's $C_{51} \ge 0.5850788196744$\* [R2026b], the C51 row of
[#184](https://github.com/teorth/optimizationproblems/pull/184), from He and Tang's reduction (arXiv:2602.12217, Theorem 2.8
and Lemma 3.1; we read Lemmas 2.1 to 2.7 and found no gap) and the two exact parameters alone (`certify_maximum_term.py`
and `audit.py` were neither read nor run). For $K > 1$ and $\varepsilon = e^{i\alpha}$, $\beta(f_{K,\varepsilon}) = 1/A$ with
$A = \max_\theta\lvert P(\theta)\rvert$, $P(\theta) = 2\sum_{n\ge0}K^{-T_n}e^{i\alpha T_n}\cos((2n+1)\theta)$, $T_n = n(n+1)/2$.
Our method is neither the note's (Bernstein subdivision in $x = \cos^2\theta$) nor He and Tang's (a mesh and a Lipschitz
bound):
- $P = P_{12} + R_{12}$ with $\lvert R_{12}\rvert \le 2K^{-91}/(1 - K^{-14})$ ($1.1\cdot10^{-50}$), an exact rational;
- $F = \lvert P_{12}\rvert^2$ on $[0, 1.5708] \supset [0, \pi/2]$ (enough, as $P(-\theta) = P(\theta)$ and $P(\pi - \theta) = -P(\theta)$;
  $\pi$ from Machin's formula). On a cell, $F(m+u) \le F(m) + F'(m)u + F''(m)u^2/2 + G_3\lvert u\rvert^3/6$ with
  $G_3 = 2(g_0g_3 + 3g_1g_2) \ge \sup\lvert F'''\rvert$ and $g_k = \sum_n 2K^{-T_n}(2n+1)^k$ (Cauchy-Schwarz on the derivative
  vectors). $F$, $F'$, $F''$ are taken in 60-digit decimal interval arithmetic with outward rounding, with $e^{i\alpha T_n}$
  as powers of $e^{i\alpha}$, so no reduction modulo $2\pi$. The quadratic's maximum and every comparison are exact
  rationals. Cells above the target are bisected;
- a lower bound at the maxima, located in floating point and evaluated in intervals.

| Row | Data (SHA-256) | Our result | Printed bound |
| --- | --- | --- | --- |
| [R2026b] | AndreRoehrig/optimization-constants-certificates @ `ed638739` (v1.0.2), `certificates/maximum_term.json` (`dcf9d1bc...3683`, as its MANIFEST states): $K$ = 3.5681854310698435, $\alpha$ = 3.9615431126097391 | $A \in [1.709171425067672530527823486722, \ldots824486723]$, so $\beta \in [0.585078819674513477028433803123, 0.585078819674513477028434145442]$. The maximum is interior ($\theta \approx 0.6455415545$), $1.16\cdot10^{-15}$ above $\lvert P(0)\rvert$; 1,050 cells | 0.5850788196744 (margin $1.13\cdot10^{-13}$) |
| control: [S2026]'s parameters | $K$ = 3.568182317714, $\alpha$ = 3.961543335688 (as the note states them) | $\beta \in [0.585078819653352950854632279223, \ldots621542]$, maximum at $\theta = 0$ | 0.585078819653 holds; so does the note's $\beta(f_{\mathrm{old}}) < 0.585078819653352950854632621542$ |
| control: [HT2026]'s parameters | $K_0$ = 7137/2000, $\alpha_0$ = 198074929/50000000 | $A \in [1.709176397414446527842885593498, \ldots886593499]$, $\beta$ = 0.5850771175595...; their mesh value at $j$ = 1794536 ($N = 5$) recomputed: 1.709176397414741464507875377..., equal to all 25 digits they print | 0.58507 |

The note's upper bound 1.709171425067816128872617735945 lies $1.4\cdot10^{-13}$ above our upper end, so it holds; the printed
0.5850788196744 is not tight (one unit higher, 0.5850788196745, is also implied). Self-tests, all passing: our sine, cosine
and $\pi$; $F'$ and $F''$ against central differences; the two symmetries; the cell bound above $F$ at 21 points of each of
300 random cells; the tail bound against the next 39 terms. Red controls, all as expected: the least 13-decimal claim
above our enclosure (0.5850788196746) is false; the claim fails with the previous parameters, with $\alpha$ raised by
$10^{-6}$, and with the truncation at $N = 2$; a cover of $[0, 3/2]$ is refused by the range check; and the cover with a
target $10^{-20}$ below our lower bound is refused, at the cell of the maximum. What rests on mathematics: He and Tang's
Theorem 2.8 (the scaling identity, the maximum term at $r = K^m$, the endpoint reduction and the negligible tail of
$k - f$) and Lemma 3.1. Reproduce: `cd c51 && ./fetch.sh && ./run_all.sh` (about 20 seconds).

## 3b: the Katz-Tao sum-difference constant (7 November 2026)

`c3b/c3b_check.py`, standard library only: Farnamfar's $C_{3b} \le 209357544569848128/114631769327023765 = 1.82634836\ldots$\*
[F2026] ([#329](https://github.com/teorth/optimizationproblems/pull/329)) and the paper's all-groups certificate
($10304335425987192362/5641833174538280329 = 1.82641618\ldots$), from the paper (`main.tex`, Sections 2 to 5) and the
certificate data alone. The package's `verify_group.py`, `stackcheck.py`, `ring_lattice.py`, `check630*.py` and
`negative_controls.py` were neither read nor run (`fetch.sh` extracts only the JSON certificates). The certificate is a
stacked conditional-copy entropy proof. Vertices are pairs with the law of a difference-injective $(X, Y)$, built by
conditional copies (Lemma 4). Entropy terms are named by rational closures of labelled linear forms: $x, y, x+y, x-y$ at
every vertex, plus the cross forms. Rows are Shannon (sub, mono), copy identities (ci) and applications of library
lemmas; equal laws identify the tuples on a source and their images on the copy. We rebuilt:
- each library lemma from its own conditional-copy proof (sub, equal-law, ci and apply rows in abstract variables; a lemma
  is applied only after it is verified): the rows sum exactly to the lemma. The Matúš family $M_1, \ldots, M_{62}$ and the
  DFZ inequalities are therefore not taken from the literature;
- the relation module step by step: at each copy, the forms supported on the source that lie in span(cl $V$) + Rel, with
  prefix closures over $\mathbb Q$ and the difference-injectivity rule;
- the name of every queried tuple. Membership in the quotient is decided modulo a Mersenne prime $p$ on integer image
  vectors, which is exact because $p$ exceeds Hadamard's bound for every minor (checked). We also re-derive which
  equal-law gens are determined by the glue at their prefix (kept) and which are moved to the copy (their image must be
  a labelled form);
- the identity $\sum_r m_r R_r = a\,h(x_0) + b\,h(y_0) + c\,h(x_0+y_0) - h(x_0, y_0)$: every class must cancel after the
  equal-law identifications.

| Certificate (SHA-256) | Library | Configuration | Rows | Result |
| --- | --- | --- | --- | --- |
| `record.json` (`12e2e356...6a70`, as the 3b row and the index state) in Zenodo 10.5281/zenodo.23714415, `source-and-certificates.zip` (`29c07851...2f158`, as #329 states) | 72 lemmas, 1,797 proof rows, all verified | 42 vertices, 9 copies, 78 cross forms; relation rank 71, quotient dimension 13, as the paper states | 519 sub, 9 ci, 12 lemma; 167 equal laws (480 gens kept, 707 moved); 1,643 tuples | all 519 classes cancel exactly: $a + b + c = 209357544569848128/114631769327023765 = 1.8263483657186587536748\ldots$ |
| `all-groups.json` (`79848549...03f1`) | 77 lemmas (five DFZ rule instances added), 1,877 rows | 42 vertices, 9 copies, 85 cross forms; rank 71, dimension 13 | 566 sub, 1 mono, 9 ci, 11 lemma; 179 equal laws (511 kept, 754 moved); 1,775 tuples | all 573 classes cancel: $10304335425987192362/5641833174538280329 = 1.8264161855212042214076\ldots$ |

Self-tests: the paper's two standard-form certificates (a chain of pages instead of a star; Shannon rows only) pass the
same code: `standard.json`, 108 vertices, rank 177, dimension 39 (here $p = 2^{521}-1$), $2762/1509$; and
`standard-five-page.json`, 30 vertices, rank 46, dimension 14, $152/83$. Every count in the paper's table of certificate sizes
agrees for all four. Red controls, all as expected for both records: refused when one sub multiplier is raised by 10%, $c$ is lowered by
$10^{-12}$, one equal-law entry is removed, one equal-law gen is changed, a coefficient in the proof of `dfz109` is changed,
a glue reading of a page step is dropped, a used transported cross form is changed, or a lemma row gets a negative
multiplier.

Rational closures are sound in torsion-free groups: clear the denominators, and multiplication by an integer is
injective. That covers $\mathbb R$, the setting of $C_{3b}$. For the paper's wider statements we also replay its domain
audits (Section 4, Lemma 9 and Proposition 10). The relation module is rebuilt over $R = \mathbb Z[1/2]$ for the record
(groups without 2-torsion) and over $\mathbb Z$ for the other certificates (every abelian group), with prefix
$R$-closures and integer echelon bases. Its rational span must be the rational module above; every equal-law gen kept
as determined by the glue must be $R$-determined at its prefix; and cl$_R(T)$ = cl$_{\mathbb Q}(T)$ must hold on every
queried tuple. All three hold for all four certificates (record: 480 kept gens and 1,643 tuples over $\mathbb Z[1/2]$;
all groups: 511 and 1,775 over $\mathbb Z$; standard: 1,172 and 2,453; five-page: 70 and 344). Audit control: the record over
$\mathbb Z$ instead of $\mathbb Z[1/2]$ is refused on 9 tuples, each missing a form $f$ whose double $2f$ is determined.
That agrees with the paper: the record needs 2 inverted. What rests on mathematics: Lemma 4 (conditional copies), Lemma 6 (removal of difference injectivity, which needs
$a, b, c \ge 0$: checked), Lemma 7 (shared and transported readings) and Corollary 12. We read them and found no gap. Reproduce:
`cd c3b && ./fetch.sh && ./run_all.sh` (about a minute).

## C3a: the finite inputs of the 1.21418 bound (7 November 2026)

`c3a/o2026_check.py`, standard library only: the finite numbers behind Olafson's $C_{3a} > 1.21418$ [O2026]
([#148](https://github.com/teorth/optimizationproblems/pull/148); the Lean proof of the finite bound is registered as
Palomar entry PALOMAR-2026-10-29-000038), from section 6 of `c3a_two_temperature_proof.md`
(carson-olaf/c3a-two-temperature-bound @ `033952dc`) alone. The package's Python and GMP verifiers, its small-model test
and its Lean files were neither read nor run. The alphabet is $M = \langle 4, 5\rangle \cap [0, 80]$, $D = M - M$,
$\kappa(d) = \min\{a + b : a - b = d\}$, with $N = 10^{13}$, $x = 8149/10000$ and $y = 42571/50000$; these are the same $M$, $N$, $x$
and $y$ as the $B = 80$ row of our C3e replay. The counts are $n_d = n_{-d} = \lfloor Nx^{\kappa(d)}/Z\rfloor$,
$n_0 = N - 2\sum n_d$ and $L = \sum_{d>0} n_d\kappa(d)$; $W$ is the type class with these counts and
$\mathcal V = \{a \in M^N : \sum a_i \le L\}$.

| Check | Result |
| --- | --- |
| $M$, $M \pm M$, $\kappa$ | $\lvert M\rvert = 75$, $\lvert M + M\rvert = 155$, $M - M = [-80, 80]$; $\kappa(d) = d$ on $M$ and 9, 18, 13, 14, 17, 19 at the gaps 1, 2, 3, 6, 7, 11, as stated |
| counts (exact integers) | $n_0 = 1917647564302$ and $2L = 76435323145630$, as stated; every $n_d \ge 148290$ |
| $\log\lvert W\rvert = \log N! - \sum_d \log n_d!$ (Robbins' bounds) | $\ge 32160395410643.4068 \ge a = 32160395000000$ |
| $\log\lvert\mathcal V + \mathcal V\rvert \le N\log P_+(y) - 2L\log y$ | $\le 25272243943117.8189 \le b = 25272245000000$ |
| $\log\lvert\mathcal V - \mathcal V\rvert \le N\log Z - 2L\log x$ | $\le 32160395412376.7043 \le r = 32160396000000$ |
| (6.1) $\alpha(r + 8) + b + 2 < a$, $\alpha = 10709/50000$ | slack exactly $227404476791/6250 = 36384716.28656$, as stated |
| our own enclosures in the same final step | $1 + (\log\lvert W\rvert - \log 4 - \log\lvert\mathcal V+\mathcal V\rvert)/(\log 1296 + \log\lvert\mathcal V-\mathcal V\rvert) \ge 1.214181180897$ |

The three finite statements the proof uses, $W \subset \mathcal V - \mathcal V$ (symmetric minimum witnesses: both budgets equal
$L$), $\lvert\mathcal V + \mathcal V\rvert y^{2L} \le P_+(y)^N$ and $\lvert\mathcal V - \mathcal V\rvert x^{2L} \le Z^N$, are also checked by
enumeration on 63 small cases ($M = \langle4,5\rangle \cap [0,10]$, $N \le 4$, $x \in \{1/2, 4/5\}$). Red controls, all as
expected: $\alpha = 10710/50000$ breaks (6.1); $x$ one unit lower changes $n_0$ and $L$; $a$ raised by $10^9$ is refused; and
in a small case, with the budget lowered by one, $W \not\subset \mathcal V - \mathcal V$. What rests on mathematics, proved in
Lean and read here: the tilted counting bounds and the type count; cyclic compression to a prime $p$ (Bertrand), with
$\lvert A - A\rvert > \lvert W\rvert/2$ and $\lvert A + A\rvert \le 2\lvert\mathcal V + \mathcal V\rvert$; the dilation $V = 2A$; and the GHR
finite-set lemma in Kleinwaks's form. Not checked: the limiting theorem of sections 3 to 5 (the display values
$1.2141808960977855\ldots$ and $1.2141811809505\ldots$). Reproduce: `cd c3a && python3 -I o2026_check.py --controls`
(about 20 seconds; `./run_all.sh` includes it).

## C89: the Ramsey multiplicity of $K_4$ (7 November 2026)

`c89a/k4_check.py` with `k4_block.c` and `k4_dense.c`: the four upper bounds on the 89a page and in
[#270](https://github.com/teorth/optimizationproblems/pull/270) that come with a public step graphon, [VSJ2026]
$0.0301388876$\* and [SK2026] $0.0301388794$\* (added in [#255](https://github.com/teorth/optimizationproblems/pull/255)),
[Kos2026] $0.0301388857$\* ([#264](https://github.com/teorth/optimizationproblems/pull/264)) and #270's [Pat2026]
$0.0301388585$\*, written from the definition on the 89a page. None of the repositories' programs was read or run
(lukevs/k4-ramsey's Lean and Python; `exact_eval.py`, `float_eval.py` and `refine.py`; `rebuild.py` and the `verifier/`
folder; `exact89.py` and `refine_lbfgs.py`). For an equal-weight step graphon $W = M/D$ on $n$ classes the density is
$(N_{\rm red} + N_{\rm blue})/(n^4D^6)$, where $N_X = \sum_{a,b,c,d} X_{ab}X_{ac}X_{ad}X_{bc}X_{bd}X_{cd} = \sum_a
\operatorname{tr}((XD_a)^3)$ runs over all ordered 4-tuples, repeats included ($X = M$ or $D - M$,
$D_a = \operatorname{diag}(X_{a\cdot})$).

- The coarse structure is checked, not assumed: the classes fall into 192 groups of $s$ consecutive classes ($s = 20$;
  40 for [Kos2026]). In all four graphons, 34,368 of the 36,864 group blocks are constant at 0 or $D$, the 192 diagonal
  ones among them, and the other 2,496 (13 per group) are general; they form 1,152 triangles.
- `k4_block.c` writes $X = PKP^{\mathsf T} + F$ ($K$ the constant blocks, $F$ the general ones). With
  $A = PKP^{\mathsf T}D_a$ and $B = FD_a$, $\operatorname{tr}((A+B)^3) = \operatorname{tr}A^3 + 3\operatorname{tr}A^2B +
  3\operatorname{tr}AB^2 + \operatorname{tr}B^3$. The first three are sums over the 192 groups. $\operatorname{tr}B^3$
  runs over the triangles of general blocks, each treated by how many of its groups meet the root's group in a general
  block. First-level sums are exact in 128-bit integers, the rest modulo $2^{64}$ and five primes; the Chinese
  remainder theorem uses four moduli and the other two check it.
- `k4_dense.c` is a second program with no coarse structure: $\operatorname{tr}((XD_a)^3)$ by dense products modulo
  $2^{61} - 1$ for two sampled roots per colour and row. It equals `k4_block`'s value at every root.
- Self-tests: both programs against brute force over all $n^4$ tuples, per root and in total, on eight random
  block-structured graphons ($n \le 20$, denominators up to $2^{30}$).
- [Kos2026]'s 7680-class witness is rebuilt from VSJ's table and its amplitudes as its `FORMAT.md` (section 6) states,
  $M'_{2i+s,\,2j+t} = M_{ij} + (-1)^{s+t}A_{ij}$. Its dense SHA-256 is the one the delta file states (`c7f0e7cb...802e`).

| Row | Data (SHA-256) | Exact density | Printed upper bound |
| --- | --- | --- | --- |
| [VSJ2026] | lukevs/k4-ramsey @ `2cf216b`, `data/constructions/final3840.json.gz` (`d7fd4c92...7370`; the JSON inside, `05302cbc...1a29`, as VSJ state) | the page's fraction exactly, $0.0301388875664972218\ldots$; $N_{\rm red} + N_{\rm blue} = 519196432373018212667371525712277942005760$, the integer of VSJ's Lean certificate | 0.0301388876 (margin $3.4 \cdot 10^{-11}$) |
| [SK2026] | JeetuSK0808/tao-constants-certificates @ `3d9431a`, `c89_stepgraphon_D20.json.gz` (`82c28f77...631a`, as its SHA256SUMS.txt states) | the page's fraction exactly, $0.0301388793132559809\ldots$ | 0.0301388794 ($8.7 \cdot 10^{-11}$) |
| [Kos2026] | DennisKoshta/c89-k4-ramsey-multiplicity-certificate @ `4fdbe0c`, the amplitudes `delta/c89_sign3_crude_0.03013888570.delta.json.gz` (`514adc3a...d7ad3`, as its SHA256SUMS states), with VSJ's table | the page's fraction exactly, $0.0301388856953920524\ldots$, $1.87 \cdot 10^{-9}$ below [VSJ2026]; $N_{\rm red} + N_{\rm blue} = 8307142402237973815060289107446603942297600$, the $S$ its README states (249,480 amplitude entries; the witness's dense SHA-256 as the delta states) | 0.0301388857 ($4.6 \cdot 10^{-12}$) |
| #270 [Pat2026] | vidigoat/optimization-records @ `fb80010`, `tao/certificates/c89/c89_stepgraphon_D30_round2.json.gz` (`60ebbced...e7d8`, as the row states) | the row's fraction exactly, $0.0301388584652867313\ldots$; the first round (`c89_stepgraphon_D30.json.gz`, `079b0d0f...8480`) gives the fraction of the folder's README exactly, $0.0301388617423003\ldots$, as the row's parenthesis says | 0.0301388585 ($3.5 \cdot 10^{-11}$) |

Red controls, for each row, all as expected: one general entry raised by $1/D$ (with its mirror) changes the fraction;
the same entry changed on one side only is refused (not symmetric); the bound one unit lower in its last printed digit
is below the exact value. What rests on mathematics: $C_{89} \le F(W)$ for every step graphon (random colourings of
blow-ups, the first-moment argument the page states), which we re-derived.

Reproduce: `cd c89a && ./fetch.sh && ./run_all.sh` (Python 3.9+ and a C compiler; about 45 minutes on a loaded 18-core
Mac, of which the 7680-class row with its controls takes about 29).

## C91: the average edit distance of two random binary strings (7 November 2026)

`c91a/ed_check.py`, standard library only: [J2026]'s $C_{91} \ge 0.216276$\*, the lower bound of the new 91a page in
[#273](https://github.com/teorth/optimizationproblems/pull/273), from the argument in its README alone and the two
numbers $x$, $t$ of `certify_rule5_R32.log` (leandrejack-afk/chvatal-sankoff-bound @ `574220e`,
`edit-distance-lower/`, `7bb076dd...3fc9`). Its `certify_*.py`, `count_check*.py`, `canon_check*.c` and the review
code were neither read nor run. Every pair of strings has an optimal alignment obeying five local rules: no I next to
a D; no S right after an indel; the match ending a run avoids the run's letters; no single match between a D and an I;
a run right after a substitution avoids the substituted letter. We re-derived the rewriting argument. The expected
number of rule-obeying valid alignments of cost at most $an$ is at most $x^{-2n}t^{-an}G(x, t)$, and $G$ is finite when
the transfer matrix $A(x, t)$ of an automaton with run lengths capped at $R$ has spectral radius below 1 (capping only
raises the weights). So $\rho(A) < 1$ and $x^2t^a > 1$ give $C_{91} \ge a$.

| Check | Result |
| --- | --- |
| canonical form: the least cost over rule-obeying alignments (our dynamic programme over rule states) against the edit distance | equal for all 262,144 binary pairs of sizes $(8, 8)$, $(7, 9)$, $(6, 10)$ and $(4, 12)$ |
| each move sequence obeying the three sequence rules: its probability, by enumeration of every letter assignment, against the product of block factors | equal for all 139 sequences of sizes $(3, 3)$, $(4, 4)$, $(3, 5)$, $(2, 6)$, $(5, 2)$ |
| our automaton (uncapped) against the sum of block products and the brute-force expected counts over all string pairs, by cost | equal at sizes $(5, 5)$, $(4, 6)$, $(6, 6)$, $(7, 5)$, $(3, 8)$ |
| the certificate: $k = 2$, $R = 32$, $x = 1178119585379/950178956526$, $t = 86841948745/634295171877$ | a rational $v > 0$ with $(Av)_i < v_i$ on all 132 states, exactly: $\max_i (Av)_i/v_i = 0.999999998115539$ (the log prints 0.999999998116). Second method: the run states eliminated in closed form, a $4 \times 4$ nonnegative matrix whose spectral radius is below 1 exactly when $\rho(A)$ is, certified the same way (0.999999997869369) |
| $a^* = -2\ln x/\ln t$, 60-digit decimal arithmetic with error bounds | $0.2162762712874192826\ldots$: $C_{91} \ge 0.216276$ holds with margin $2.7 \cdot 10^{-7}$ |

Red controls, all as expected: $x$ raised by 1% gives $\rho(A) \ge 1.0189$ (a Collatz-Wielandt lower bound, so no
certificate exists); without rule 5, the automaton at the same $x$ and $t$ has $\rho \ge 1.0062$; the claim $0.216277$
is above $a^*$. What rests on mathematics: the rewriting argument for the five rules (checked exhaustively at the sizes
above), the independence of the block factors (checked by enumeration) and the first-moment step.

Reproduce: `cd c91a && ./fetch.sh && ./run_all.sh` (about two minutes).

## C42: an explicit value of Biró's lower bound (8 November 2026)

`c42/biro_check.py`, standard library only: [Pin26c]'s $C_{42} \ge 1/2 + 2^{-4096}$\* of
[#310](https://github.com/teorth/optimizationproblems/pull/310), from `LOWER_PROOF.md` and `inputs/lower_certificate.json`
(rpinho/math-certificates @ `36325fa`, `c42/biro-explicit-lower`; the hashes of its `lower_manifest.json`). Its
`verify_lower_exact.py` was neither read nor run. We re-derived sections 1 to 8 of the proof (Biró's first escape of
the $b$ increments, radial balance, the summed coefficients becoming parallel, the transfer back to $b$) and found no
false step. The program checks every numerical step, (5) to (32), exactly: inequalities in $r = 1/\sqrt 2$ and
$\cos(\pi/8)$ in $\mathbb{Q}(\sqrt 2)$, square roots by squaring, all parameters the file's powers of two. Each stated
intermediate exponent is shown to follow ($t_0 > 2^{1021}$ with $2^{0.5}$ to spare, $k/t_0 > 2^{1022}$,
$b \le 2^{-254}$, the product below $2^{-340}$, $E/A < 0.0386 < 1/16$ against $r(r-1/2) = 0.1464\ldots$). The identities
behind (11), (15), (23) and (29) are expanded exactly, and (1) and (3) are checked on 40 random tuples of Gaussian
rationals. Red controls, all refused: the claim strengthened to $1/2 + 2^{-4094}$ with the same auxiliaries,
$\eta = 2^{-200}$, $\lambda = 2^{-4000}$, $Z = 2^{-256}$, $\varepsilon = 2^{-252}$, $\sigma = 2^{-8}$, final error
$1/32$, product bound $2^{-341}$, and a perturbed $b_1$ in the identities. What rests on mathematics: the geometric and
analytic steps (read, not formalized). Reproduce: `cd c42 && ./fetch.sh && python3 -I biro_check.py certs --controls`.

## C68: the contact-row certificate for 0.447 (8 November 2026)

`c68/row_check.py`, standard library only: $C_{68} \ge 447/1000$\* of
[#299](https://github.com/teorth/optimizationproblems/pull/299) (the contact-row model), from `PROOF.md`,
`CONTACT_ROW_LEMMA.md` and the two data files `certificate-0447.json` and `row-cover.json` (rpinho/math-certificates
@ `09c8423`, `c68/contact-row/certificates`; `d9dd308f...2afe`, `fda0aa14...fc27`). Its `verify_measure.py` and
`verify_row.py` were neither read nor run.

| Check | Result |
| --- | --- |
| section 1 at $c = 447/1000$ (110-digit decimal) | $\Psi = 26316.4966 > 6$; $\mathrm{sech}^2(A/2) < 1/3$; $16\sinh^2(A/4) = 1827.35 < 2(1-8\kappa)\Psi = 2105.32$; $\tanh^2(A/8) = 0.82955 > (1+8\kappa)\pi^2/64$; $\mathrm{sech}^2(A/8) = 0.17045 < (1-4\kappa)/2$ |
| the row cover for (5), $z = 23/200$ | outside the root box $S > z$ exactly; the tree tiles $[0,17/50]^2 \times [7/10,1]$ (287 leaves); 17 PSD leaves fail (6) at their low corner and 4 constant leaves have $\max(S,J) > z$ (exact rationals); 265 grid leaves (4,096 to 32,768 cells; outward-rounded binary64, `Decimal.ln`) satisfy $\exp(\mathrm{mean}\, l_j) > z(1+D^2)$, smallest ratio 0.11502033 (the cover's `proposal_min`) |
| the small leaf | 8,704 bins, masses bounded below by $1/(2\sqrt{X(1-X)})$: $\log(\ell/16) + \sum \mathrm{mass}\log\gamma_v = -2.147125 > \log z = -2.162823$ |
| section 5 (320-bit fixed point) | $A_k \ge 0$ for $k \le 5998$ (smallest relative margin $7.8656 \cdot 10^{-6}$ at $k = 172$), then the atom at $R$; $B_k \ge 0$ for $k \le 1322$ with the atomic $\varepsilon$ term at $k = 1$ ($B_1 \ge 3.854 \cdot 10^{-6}$), then one piece of $[R,1]$ alone; our own slope bounds $g \ge l_i(1-t)$ on 42 pieces |
| the contact-row lemma, sanity test | 20,000 random two- and three-contact configurations: $\max_i T(\theta_i) \ge J$, (6) to rounding, (7) equals $a(q_0+a)/(1+a)$ |

Red controls, all refused: $z = 0.11503$; $\kappa = 29/250$; a split point outside its box; a grid leaf coarsened to
64 cells; half of atom 2's mass moved to atom 3; $c = 0.448$; the $k = 1$ $\varepsilon$ term with 19/200 (the last
three fail at $B_1$). What rests on mathematics: the Green/flux and angular inputs, the minimal enclosing disk, the
contact-row inequality and $J$, the outer-function steps behind (8), the small-region outer factors and the norm
implication. Reproduce: `cd c68 && ./fetch.sh && python3 -I row_check.py certs --controls` (about 11 minutes).

## C84a: the unit distance exponent, Sawin's criterion (8 November 2026)

`c84a/sawin_check.py`, standard library only: the 84a row $C_{84} \ge 1.0152$\* of
[#315](https://github.com/teorth/optimizationproblems/pull/315) (Emmerich, arXiv:2606.03419v5), from the data file
`certificate_verification_all_v23.json` (emmerichmtm/UnitDistanceProblemOptimizationOfSawinsLowerBound @ `2279d0a`,
`15fb28e3...ef86`; its parameters are those printed in Appendix B of v5). The repository's verify and optimize scripts
were neither read nor run. The criterion is Sawin's (arXiv:2605.20579, proof of Theorem 1 with Lemmas 11 and 12 and
Proposition 10): $T$ odd primes with an odd number $\equiv 3 \bmod 4$; every $p \in S_{\mathbb Q}$ either $1 \bmod 4$ or
inert in some $\mathbb{Q}(\sqrt q)$, $q \in T$; the budget $\#T + \#S_{\mathbb Q} + \#\{\text{split}\} + 1 \le (\#T-1)^2/4$
in $Q = \mathbb{Q}(\sqrt{\prod T})$; then formula (11). Splitting is decided by our own Legendre and Kronecker symbols,
and $\delta$ is enclosed in 60-digit decimal arithmetic.

| Certificate | $\#S_{\mathbb Q}$ (ramified, inert, split) | budget | $\delta$ |
| --- | --- | --- | --- |
| Sawin's published example, $R = 72$ | 22 (10, 12, 0) | 36 = 36 | 0.01411442867849823889... (his .014114...) |
| greedy, $R = 66.72240803$ | 22 (2, 20, 0) | 36 = 36 | 0.01517180563721325854... |
| evolution strategy, $R = 66.72416$ | 22 (3, 19, 0) | 36 = 36 | 0.01526166106841929424... |
| with discrete recombination, $R = 66.72416$ (the row) | 22 (3, 19, 0) | 36 = 36 | 0.01526286881700719166... > 0.0152 |

All four agree with the file's stored values and splitting types, compared only after our run. Red controls, all
refused: the exponent tightened to 1.0152629; the split prime 53 added; 263 replaced by 951851 ($\equiv 3 \bmod 4$,
inert in no $\mathbb{Q}(\sqrt q)$); $T$ without 43; $R = 1$. What rests on mathematics: Sawin's Proposition 10,
Lemma 11 (Golod-Shafarevich) and Lemma 12. Reproduce: `cd c84a && ./fetch.sh && ./run_all.sh` (about a second).

## C63 and C64: the Li-Yang exponent $\theta^*$ (8 November 2026)

`c63/theta_check.py`, standard library only, exact: the rows $C_{63} \le 0.3144831759741\ldots$ of
[#313](https://github.com/teorth/optimizationproblems/pull/313) and $C_{64} \le 2\theta^* = 0.6289663519\ldots$ of
[#315](https://github.com/teorth/optimizationproblems/pull/315), both from Li and Yang (arXiv:2308.14859v2), whose
Definition 1.1 makes $-\theta^*$ the unique solution on $[-0.35, -0.3]$ of
$-\tfrac{8}{25}x - \tfrac{1}{200}\big(\sqrt{2(1-14x)} - 5\sqrt{-1-8x}\big)^2 + \tfrac{51}{200} = -x$. The equation is
read from the paper's LaTeX source (the arXiv e-print, `6fe6f9c0...93de`) and compared with our transcription first.

| Check | Result |
| --- | --- |
| route A, exact algebra | (1.4) is $uv = -\tfrac{182}{5}x - \tfrac{37}{5}$; squared, $27524x^2 + 13168x + 1419 = 0$; of its two roots only $x = -(3292 + 25\sqrt{1717})/13762$ lies in the interval with $uv \ge 0$, and the residual is exactly zero |
| route B, interval arithmetic ($2^{-400}$ fixed point, outward rounding) | $f' \in [0.942, 1.127]$ on the whole interval, so at most one root; a sign change at the ends; the bisection enclosure contains route A's root |
| $\theta^*$ | $(3292 + 25\sqrt{1717})/13762 = 0.31448317597406140673\ldots$ |
| 63a, printed $0.3144831759741$ | above $\theta^*$: a valid upper bound as printed (the paper's "$0.3144831759741\cdots$" is a rounding, not a truncation) |
| 64a, README cell $0.6289663519$ | below $2\theta^* = 0.62896635194812281347\ldots$: not an upper bound as printed; $0.6289663520$ is |

Red controls, all as expected: $51/200$ changed to $52/200$ (the root moves to $0.3193626\ldots$ and the printed 63a
value no longer bounds it); the 63a value lowered by $10^{-13}$; the interval moved to $[-0.2, -0.1]$, which holds
only the spurious root of the squared equation; the sign inside the square flipped. What rests on mathematics: the
proof of Theorem 1.2, and for 64a the change to the radius normalization ($N(t) - \pi t^2 = R(t^2)$). Reproduce:
`cd c63 && ./fetch.sh && ./run_all.sh` (under a second; `c64/run_all.sh` runs the same check).

## C87: Martinet's totally real tower, $4\sqrt{70035}$ (8 November 2026)

`c87/martinet_check.py`, standard library only, exact: the 87a row as
[#332](https://github.com/teorth/optimizationproblems/pull/332) restates it, from Martinet (Invent. Math. 44 (1978),
Proposition 4.1 and Example 4.2, p. 69): $k = \mathbb{Q}(\sqrt{29}, \sqrt 2)$, $K = k(\sqrt{3 \cdot 5 \cdot 7 \cdot 23})$
and $k' = \mathbb{Q}(\sqrt 2, \sqrt{70035})$. Every ramification and splitting number comes from the $\mathbb{F}_2$-space
of quadratic subfields (inertia and decomposition fields as spans of the subfields unramified at $\ell$, resp. split at
$\ell$), and the counts in $k$ above $3, 5, 7, 23$ again from Dedekind-Kummer, with a defining polynomial whose index is
checked prime to $\ell$.

| Check | Result |
| --- | --- |
| Proposition 4.1 (i)-(iii) | $29 \equiv 1 \bmod 4$; $m = 2415 \equiv -1 \bmod 4$; $7$ and $23$ are $\pm 1 \bmod 8$ with $(q/29) = +1$ |
| $\mathrm{disc}(k')$ | $8 \cdot 280140 \cdot 560280 = 2^8 \cdot 70035^2$, so $\mathrm{rd}(k') = 4\sqrt{70035}$ |
| $K/k'$ | unramified at every prime ($e_\ell(K) = e_\ell(k')$ for $\ell \mid 2 \cdot 70035$), all fields totally real |
| primes of $k$ ramified in $K$ | 1 above 2, 2 above 3, 2 above 5, 4 above 7, 4 above 23 (none above 29): $r = 13 \ge n + 3 + 2\sqrt{2n+1} = 13$ |
| $4\sqrt{70035}$ | $1058.565066493316696961\ldots$; the old printed $1058.565$ is below it, as #332 says |

Red controls, all refused: $q_4 = 19$ and $q_3 = 11$ ($r = 11$); $p = 13$; the threshold raised to 14; $m \equiv 1 \bmod 4$.
What rests on mathematics: Martinet's section 2 (genus theory and the Golod-Shafarevich criterion behind
$r \ge n + 3 + 2\sqrt{2n+1}$) and the conductor-discriminant formula. Reproduce: `cd c87 && ./run_all.sh` (no download).

## C1a and C1b: eleven step-function rows (8 November 2026)

`../autocorr_check.py` (a `1b` mode added: $(2/n)\max_m \sum_i h_i(1 - h_{i+m})$ over every grid shift, with $0 \le h \le 1$
and $\sum h = n/2$ checked exactly; the autocorrelation from one big-integer product), on the heights of eleven rows,
each from the source its row cites; `c1ab/extract_values.py` copies the numbers out as text (Python files and notebook
cells with `ast`, the LaTeX with a regular expression; nothing downloaded is run). Both readings of the decimals are
run: as binary64 and as printed.

| Row | Source (pinned) | $n$ | Exact value (binary64; printed decimals) | Printed |
| --- | --- | --- | --- | --- |
| 1a [GGSWT2025] May 2025 | google-deepmind/alphaevolve_results @ `4226acb`, notebook cell 53 | 600 | 1.50529396844016080772...; 1.50529396844016081649... | 1.5053, valid |
| 1a [GGSWT2025] Dec 2025 | alphaevolve_repository_of_problems @ `8f44745`, notebook cell 60 | 1319 | 1.50316355468156091584...; 1.50316355468156092083... | 1.503164, valid |
| 1a [YKLBMWKCZGS2026] | test-time-training/discover @ `6c40e82`, `ttt_ac1_sequence.json` | 30000 | 1.50286289825582656290...; 1.50286289825582656243... | 1.502870, valid |
| 1a [T2026] | togethercomputer/EinsteinArena-new-SOTA @ `c388c6f`, `ours_2026.json` | 30000 | 1.50286285870531082009...; 1.50286285870531082455... | **1.502862, below the exact value** |
| 1a [EA2026] | EinsteinArena solution 2603 (fetched 8 Nov) | 65536 | 1.50274353195395969900...; 1.50274353195395970066... | 1.5027435320, valid |
| 1b [GGSWT2025] | alphaevolve_results @ `4226acb`, cell 78: 48 numbers, mirrored to 95 steps | 95 | 0.38092303510845014962...; 0.38092303510845014897... (rescaled) | 0.380924, valid |
| 1b [YKLBMWKCZGS2026] | discover @ `6c40e82`, `ttt_erdos_sequence.json` | 600 | 0.38087532321771874649...; 0.38087532321771874864... (rescaled) | 0.380876, valid |
| 1b [T2026] | EinsteinArena-new-SOTA @ `c388c6f`, `together_ai_2026.py` | 600 | 0.38087031058621990349...; 0.38087031058621990370... (rescaled) | 0.380871, valid |
| 1b [EA2026] | EinsteinArena solution 2507 (fetched 8 Nov) | 3584 | 0.38085857485785830914...; 0.38085857485785830932... (rescaled) | 0.3808585749, valid |
| 1b [A2026] | atrislabs/records @ `1fa857a`, `construction.json` (`62a976fa...`) | 114688 | 0.38085855031642856010... (exact rationals) | 0.3808585504, valid |
| 1b [H2016] | Haugland, arXiv:1609.08000, the 51-step function of section 2 (printed decimals) | 51 | 0.38092685343308697464... (the paper: 0.3809268534330870) | **0.380926, below the exact value** |

"Rescaled": the heights as published sum to $n/2$ only up to about $10^{-13}$ or less, so they are scaled to
$\sum h = n/2$ exactly (the repair described in #232), staying in $[0, 1]$. The 48 AlphaEvolve 1b numbers alone
sum to 24.246, not 24; mirrored (the 48, then the first 47 reversed) they sum to $47.5$ up to $10^{-20}$. Every maximum
is attained at one shift, except the AlphaEvolve and Haugland 1b functions (two, mirror images). Red controls on every
row: a corrupted coefficient fails the polynomial identity; for 1a, the old coefficients fail after one height changes
and scaling by 3 leaves the bound unchanged; for 1b, a height of 2 is refused, a sum moved off $n/2$ is refused, and
moving mass $1/D$ between cells changes the bound. The 1a [T2026] row should read 1.502863 (the repository's README
gives 1.50286286); the 1b [H2016] row read 0.380927 until
[#158](https://github.com/teorth/optimizationproblems/pull/158) changed it to the abstract's truncated 0.380926.
Reproduce: `cd c1ab && ./fetch.sh && ./run_all.sh` (about 70 seconds). `run_all.sh` checks every run against a table of expected outcomes (exit 0, except the three runs whose printed value the row's own heights do not support, tog_c1a_ours_2026 in both readings and c1b_haugland_2016 read as decimals, which must exit 1 with the "NOT a valid ... bound as printed" line) and exits non-zero only on an unexpected outcome, naming the row and why.

## C24: the Komlós constant of Akbas and Sra (8 November 2026)

`c24/komlos_check.py`, standard library only: the row $C_{24} \le 6.9013$\* of
[#250](https://github.com/teorth/optimizationproblems/pull/250), from Theorem 2.4 of Akbas and Sra
(arXiv:2609.27172v2), whose constant $C_{\cos}$ is the unique solution of
$F(C) = \tfrac{3\sqrt{2\pi}}{C}\int_0^3 (1 - r/3)^2 e^{-\pi^2 r^2/(2C^2)}\,dr = 1$. The equation is read from the paper's
LaTeX source (`komlos-paper.tex`, `386c47a6...92a2`) and compared with our transcription first. $F$ is strictly
decreasing in $C$ (our own one-line argument: with $s = \pi/C$, $dF/ds$ is a positive multiple of
$\int_0^{3s} w(u) e^{-u^2/2}(1 - u^2)\,du$ with $w$ non-negative and decreasing, which is at least
$w(1) \cdot 3s\,e^{-9s^2/2} > 0$), so $C_{\cos} < c$ exactly when $F(c) < 1$. Two evaluations in $2^{-400}$ fixed-point
interval arithmetic, with $\pi$ from Gauss's arctangent formula: the term-by-term alternating series (route A) and the
error-function form with positive-term series (route B).

| Check | Result |
| --- | --- |
| $F(6) - 1$ | $+0.12191512\ldots$ (both routes) |
| $F(6.9012574) - 1$, $F(6.9012575) - 1$ | $+1.856 \cdot 10^{-10}$, $-1.201 \cdot 10^{-8}$: the paper's enclosure $(6.9012574, 6.9012575)$ holds |
| $F(6.9013) - 1$ | $-5.196 \cdot 10^{-6}$: $C_{\cos} < 6.9013$, the row |
| $C_{\cos}$ by bisection | $6.901257401521836873264072186907\ldots$ (width $8 \cdot 10^{-38}$; route B agrees at both ends) |

Red controls, all as expected: the bound tightened to 6.9012 is refused; with the weight $(1 - r/3)$ instead of its
square, $F(6.9013) > 1$; with the 2 dropped from $2C^2$ the enclosure is refused; route A cut at three terms (with its
remainder bound) cannot decide the sign at 6.9013. What rests on mathematics: the proof of Theorem 2.4. Reproduce:
`cd c24 && ./fetch.sh && ./run_all.sh` (under a second).

## C9: the 134,753-vertex independent set of [IRCR2026] (8 November 2026)

`c9/strong_check.c`, our own C program: the 9a row $134753^{1/10} \approx 3.258020$ (Itty, Rosin, Carstensen and
Reichman, arXiv:2607.21517; added in [#174](https://github.com/teorth/optimizationproblems/pull/174)), from the data file
`c7/CC_10_7_134753.txt` of nathanielitty/lower-bounds-for-shannon-capacity @ `fe89660` (`60fb7b0a...755a`); the
repository's `verifier.py` and construction scripts were neither read nor run. Two members $x \ne y$ of
$\mathcal{C}_7^{\boxtimes 10}$ are adjacent when every coordinate differs by $0$ or $\pm 1 \bmod 7$. Route 1 looks up all
$3^{10} - 1 = 59048$ neighbours of every member in a bitmap of the $7^{10}$ vertices; route 2 groups the members by their
first five coordinates and compares each member with the members of the $3^5$ adjacent or equal groups on the last five.

| Check | Result |
| --- | --- |
| the file | 134,753 lines of ten values in $\{0, \ldots, 6\}$, no duplicates |
| adjacent ordered pairs | 0 in route 1, 0 in route 2: the set is independent, so $\alpha(\mathcal{C}_7^{\boxtimes 10}) \ge 134753$ |
| the printed value | $3.258020^{10} \le 134753 < 3.258021^{10}$ exactly; $134753^{1/10} = 3.25802073729324535952\ldots$ |

Red controls, each refused: member 0 with its first coordinate moved by one appended (8 adjacent ordered pairs in both
routes); a copy of member 0 appended (refused as a duplicate); the first 500 members under the wrong graph "differ by at
most 2" (48,176 pairs in both routes). Reproduce: `cd c9 && ./fetch.sh && ./run_all.sh` (a C compiler; about 90 seconds
with the controls, 40 MB).

## C84b: Althoefer's chain for the sum-product exponent (8 November 2026)

`c84b/althoefer_check.py`, standard library only: the 84b row $C_{84b} \le 1.9993$\* ([Al26], Althoefer's note on
[BSSZ2026, §5]) and the ceiling that [#275](https://github.com/teorth/optimizationproblems/pull/275) restates, from the
note's TeX alone (`improved_constant_052.tex`, `2f6bbdcb...5568`; its numbers are read from the file and compared). The
chain: $c = -\log\max(5\rho/Y,\ M_1e^Y/(XY^2) + M_2/(X^2Y^2)) / L$ with $L > \log(5Y+1) + \log(2\varepsilon X + 1)$,
$\varepsilon = (\varphi-1)/(\varphi+1)$, $M_2 = \rho^2C_2/\varepsilon^2$, $M_1 = 4(1+\varepsilon)M_2$, $C_2 = 857.57$, and
$\rho(s) = \zeta(s)\Gamma(s/2)C_2^{s/2}/(2\pi^{s/2})$ bounding the regulator. $\zeta$ by Euler-Maclaurin and $\log\Gamma$ by
Stirling's series, both with the first omitted term as the remainder bound, in 50-digit decimals.

| Check | Result |
| --- | --- |
| the note's chain at $s_0 = 1.371966384$, $\rho = 103$, $Y = 1415$, $X = \lfloor e^{1423}\rfloor$ | $\rho(s_0) = 101.95607583\ldots \le 103$; $515/Y = 0.36395760 < 0.364$; $M_1e^Y/(XY^2) \le 0.13523948 < 0.136$; $M_2/(X^2Y^2) < 10^{-1000}$; $L \le 1431.1139758 \le 1432$; $-\log 0.364 = 1.010601 > 0.0007 \cdot 1432$: $c \ge 0.0007$, the exponent $1.9993$ |
| $c$ at the note's choice, unrounded | $0.000706245573\ldots$ (#275 and [T2026]: about 0.000706) |
| the ceiling over every $s > 1$, $X, Y > 0$ | $\rho(s) \ge 101.953216$ for all $s > 1$ (5,115 cells); $c \le \log(Y/(5\rho))/(Y + \log(8(1+\varepsilon)\rho C_2/\varepsilon))$ for every $X$; so $c \le 0.000713944$, which the chain approaches to $4.5 \cdot 10^{-8}$ ($Y = 1400.8$) |

So the row's ceiling ($c \le 0.0007150508$, or $\le 0.0007150507$ before #275) holds with room; we did not reproduce its
value $0.00071505073\ldots$, which lies above our supremum. Red controls, all refused: $\rho(s_0) \le 101.95$; $Y = 1414$
with the factor $0.364$ kept; $X = \lfloor e^{1421}\rfloor$ (the additive factor becomes 0.9993); $L = 1444$; a ceiling of
$0.0007138$. What rests on mathematics: BSSZ2026 §5, the note's lattice-doubling lemma and regulator inequality, and the
$(1+o(1))^d$ absorption. Reproduce: `cd c84b && ./fetch.sh && ./run_all.sh` (about 25 seconds).

9 November: the record's own value, $0.00071505073\ldots$, is the same chain bounded more loosely in three places, as
carlostoledo1891 explained on #275: $X$ only needs the additive factor below 1, so
$L > \log(5Y+1) + \log(2\varepsilon M_1) + Y - 2\log Y$; $\rho$ is the floor $101.9$; and each $Y$-cell of width $1/2$ is
bounded with the numerator at its right end and the denominator at its left end. With these, `record_ceiling` gives
$0.00071505073368073619\ldots$ at $1399 \le Y \le 1399.5$ (cells from $5\rho$ to 4000, a monotone tail beyond), so #275's
$0.0007150508$ is that value rounded up. Red control: the same chain with $\rho = 101.953216$ gives $0.0007146808$.

## 84a: Sawin's example and three finite-data rows (9 November 2026)

`c84a/forum_check.py`, standard library only, on the criterion of our 1.0152\* replay above (Sawin, arXiv:2605.20579,
Lemmas 11 and 12 and Proposition 10): Sawin's own example ([S2026], the data printed in the proof of Theorem 1), and the
three rows of the 84a page marked "Unverified": [EPF90] (mlewko's data file from Erdős Problems forum thread 90,
`2783627d...3f14`, read as data) and spiderduckpig's two [MO511514] rows, whose $S$ we rebuild from the posts' rules.

| Row | $\#T$, $\#S$, $U$ | Budget | $\delta$ |
| --- | --- | --- | --- |
| [S2026] $1.014$ | 13, 22, 0 | $36 = 36$ | $0.0141144286784982389\ldots$ |
| [EPF90] $1.03184$ | 74, 1082, 175 | $1332 \le 1332.25$ | $0.0318499816702021619\ldots$ |
| [MO511514] $1.03188$ | 74, 1092, 165 | $1332 \le 1332.25$ | $0.0318830655396396548\ldots$ |
| [MO511514] $1.03190$ | 74, 1088, 169 | $1332 \le 1332.25$ | $0.0319083717352678483\ldots$ |

The parity and inertness conditions hold (for [EPF90] the file's witness table, compared after our search, lists the
same 546 primes); the weight histograms are the posted ones, and each weight is the best for its prime at the row's
$\delta$; numerators and denominators equal the posts' to every printed decimal (for 1.03188 the MathOverflow answer's;
the forum copy prints both about $5 \cdot 10^{-12}$ lower). Red controls, all refused: $10^{-15}$ above each enclosure,
the next split and non-split primes (budget), an even parity $T$, $T$ without its primes $5 \bmod 8$ (no witness for
$p = 2$), $R = 1$, the 1.03188 threshold raised. Reproduce: `cd c84a && ./fetch.sh && ./run_all.sh` (seconds).

## 37a: no Boolean function of degree 4 has sensitivity 10 (9 November 2026)

`c37/c37_check.py`, standard library only: the statement of
[#283](https://github.com/teorth/optimizationproblems/pull/283) (DennisKoshta), from the definition and the three
certificate files alone (`d10_4.cert.xz`, `3be6fa8b...e7e9f4`, decompressed `7dc38b37...7e6c`; `d7_3.cert`,
`d4_2.cert`), read in the format of the packet's `WRITEUP.md` §6; its checker, generator, verifier and tools were
neither read nor run. We rebuild the 638 Möbius conditions $a_S = 0$ ($\lvert S\rvert \ge 5$) from $(n, d) = (10, 4)$,
walk the branch-and-bound tree (478 branchings, 479 leaves, depth 15) and check at every leaf that
$\max_{g \in B_\ell} (y^\top A) g < 0$ exactly ($y^\top A$ by a superset-sum transform). As checks of $A$, $c \cdot h = 0$ at
every leaf for Proskurin's symmetric point (which solves the LP relaxation exactly) and for $\mathrm{NAE}_3 \circ
\mathrm{NAE}_3$; that function has degree 4 and sensitivity 9. Red controls, all refused: a leaf multiplier negated, the
last leaf deleted, a branch on a fixed coordinate. Reproduce: `cd c37 && ./fetch.sh && ./run_all.sh` (about 3 seconds).

## 3a: two capped digit constructions recounted exactly (9 November 2026)

`c3a/g2026_check.py` and `c3a/mi2026_check.py` (with `diff2d.c`), standard library only: the rows [G2026] $1.1740744$
([#71](https://github.com/teorth/optimizationproblems/pull/71); $A = \{0,2,\ldots,10\}$, base 21, 80 digits, digit sum
at most 150) and [MI2026] $1.1835129324$ ([#95](https://github.com/teorth/optimizationproblems/pull/95);
$A = \{0,3,4,6,\ldots,16\}$, base 33, 420 digits, cap 1392; parameters from the counts file of Zenodo 20794135, read as
data). No digit operation carries, so $\lvert U+U\rvert$ counts digit vectors $s$ with a split of both digit sums at most
$T$: that holds iff $n + \mu \le 2T$ ($n = \sum s_i$, $\mu$ the least $\lvert\sum y_i\rvert$, $y_i \in \{2a - s_i\}$),
counted by a dynamic program over (digit sum, reachable sums) with a saturation lemma (once every value of the right
parity in $[-K, K]$ is reachable it stays so). $\lvert U-U\rvert$ is a box sum over $[0,T]^2$ of the coefficients of a
$d$-th power, by the coordinatewise least split (in C modulo 35 primes for [MI2026], put together by the Chinese
remainder theorem, two more primes agreeing). Both programs equal brute force on 13 small instances.

| Row | $\lvert U+U\rvert$, $\lvert U-U\rvert$, $\max U$ | Lemma value | Run |
| --- | --- | --- | --- |
| [G2026] | the pull request's three integers (77, 96, 106 digits) | $1.17407444769352116336\ldots$ | 26 s |
| [MI2026] | the certificate's three integers (474, 591, 638 digits; at most 311,877 unsaturated states) | $1.18351293242186151065\ldots$ | 3.2 h |

Red controls, all refused: one unit more in each row's last digit; the cap one lower (the counts change); a digit that
would carry. Reproduce: `cd c3a && ./run_all.sh` (the [MI2026] row about three hours).

## C1c, C1a and C1d: six step functions of 9 and 10 November 2026, and the first 1c certificate of #270

`c1ac/` replays the certificates of [#354](https://github.com/teorth/optimizationproblems/pull/354) (SeanSagait) and of
the 9 and 10 November updates of [#270](https://github.com/teorth/optimizationproblems/pull/270) (vidigoat) with
`autocorr_check.py`, the checker of the first table (its method section applies unchanged). The certificates' own
verifiers and search programs were neither read nor run. #354's file is text (K, then the weights); `fetch.sh` copies
its integers into JSON after checking that every token is a canonical non-negative integer.

| Constant | Pull request | Certificate (pinned commit, SHA-256) | Exact value | Printed bound |
| --- | --- | --- | --- | --- |
| $C_{1c} \geq$ | [#354](https://github.com/teorth/optimizationproblems/pull/354) | SeanSagait/c1c-autocorrelation-certificate @ `9e7f3d19`, `K262144_0.4103768.txt` (K = 262,144, 262,400 weights), `062a2869...6769` | 935577600963641344/2279800930147535089 = 0.41037688360934932768..., minimum only at k = 46008 (S = 1509900967, A = 3568945316176) | 0.4103768, below the exact value |
| $C_{1c} \geq$ | #270 (26 Oct, its first certificate: the [Pat2026] row, replayed 10 Nov) | vidigoat/optimization-records @ `08c6bfdf`, `tao/certificates/c1c/c1c_K98304.json` (K = 98,304, 98,368 weights), `3fcd6581...f747` | 2594500369901311586790309888/6322755819680014072464230569 = 0.41034328129923190079..., minimum only at k = 97663 | 0.4103432812, below |
| $C_{1c} \geq$ | #270 (9 Nov, `c1c_v4`) | vidigoat/optimization-records @ `358b200b`, `tao/certificates/c1c_v4/c1c_K98304_v4.json`, `58dacf40...16f6` | 5282641762424933775753216/12873173145647069582596129 = 0.41036049951765035306..., minimum only at k = 98288 | 0.4103604995, below |
| $C_{1c} \geq$ | #270 (9 Nov, `c1c_v5`) | vidigoat/optimization-records @ `1b39a28d`, `tao/certificates/c1c_v5/c1c_K98304_v5.json`, `7021f269...676b` | 2556315663995378551725293568/6229433836731244234885554241 = 0.41036083390473056108..., minimum only at k = 45081 | 0.4103608339, below |
| $C_{1a} \leq$ | #270 (9 Nov, `c1a_v2`) | vidigoat/optimization-records @ `1b39a28d`, `tao/certificates/c1a_v2/c1a_n262144_v2.json` (262,144 non-negative integer heights), `4d8f8216...fd8d` | 29134419116089766846043850668467284344832/19387496227028661445940960928888412813689 = 1.50274273557163313362..., maximum only at m = 96751 | 1.5027427356, above |
| $C_{1d} \leq$ | #270 (10 Nov, `c1d`) | vidigoat/optimization-records @ `a559680e`, `tao/certificates/c1d/c1d_n163840.json` (163,840 signed integer values), `b05de6ab...77b6` | 646233848324886328286774558720/444246710235775351931155445761 = 1.45467334576751329963..., maximum of $\lvert c_m\rvert$ only at m = 165229 | 1.4546733458, above |
| $C_{1a} \leq$ | #270 (10 Nov, `c1a_v3`) | vidigoat/optimization-records @ `a559680e`, `tao/certificates/c1a_v3/c1a_n262144_v3.json`, `b32a20bf...5c14` | 116537598371962553887369519059399403896832/77549933019818026777245548252121251545225 = 1.50274273405473036507..., maximum only at m = 262143 | 1.5027427341, above |

Every fraction, deciding lag and printed bound is the one the pull request, its rows or its thread state. The hashes of
#354, `c1c_v4` and the first 1c certificate (the row states `3fcd6581f2c1`) are the ones stated there; for `c1c_v5`, `c1a_v2`, `c1d` and `c1a_v3` none was stated, and these are
the files at the commit the thread links (their folders' SHA256SUMS list the same hashes). Red controls, on every certificate: a
corrupted coefficient and a weight changed by 1 are refused by the polynomial identity, scaling by 3 leaves the bound
unchanged, and the printed bound moved one unit in its last digit past the exact value (0.4103768837, 0.4103432813, 0.4103604996,
0.4103608340, 1.5027427355, 1.4546733457, 1.5027427340) is refused. Reproduce: `cd c1ac && ./fetch.sh && ./run_all.sh`
(about five minutes).

## 22a: the ropelength lower bound 2.124 (10 November 2026)

`c22a/c22a_check.py`, standard library only: the finite part of [V2026] $C_{22a} \geq 2.124$
([#356](https://github.com/teorth/optimizationproblems/pull/356), Varlamov), written from the paper
(`paper.pdf` @ `6d11254`, `cdf1c6d1...a6fcc`, Sections 4 to 6); its `verify.py` was neither read nor run. With
$k = 50657/10000$ and $A = k(3/4)^{1/3}$ the paper needs $H(\ell) = \int_2^\infty 2t^{-3}\min\{\ell, \tfrac43(t+1)^3\}\,dt
\le A\ell^{1/3}$ for every $\ell > 0$, and $B = (256\pi^3/(3k^3))^{1/4} > 531/250$. Our route differs from the paper's
completed square: our own primitive of $\tfrac83(t+1)^3/t^3$, checked as an identity of Laurent polynomials, gives
$H = \ell/4$ for $\ell \le 36$ and $h(R) = 4R + 8\log(R/2) + 3 - 4/R$ for $\ell = \tfrac43(R+1)^3$; $\ell \le 36$ needs
$1296 < 48k^3$; for $R \ge 2$, $k(R+1) > h(R)$ is proved on $[2, 100]$ by 2,875 dyadic intervals (both sides increase,
so $\overline{h}(b) < k(a+1)$ suffices, with rational upper bounds for $\log$ from the $\operatorname{artanh}$ series and a
geometric tail) and beyond 100 by the tangent line of $\log$; $\pi$ comes from Gauss's formula
$48\arctan\frac1{18} + 32\arctan\frac1{57} - 20\arctan\frac1{239}$, itself confirmed by Gaussian integers. Then
$(531/250)^4 < B^4$ and $2.12403897858 < B < 2.12403897859$. Red controls, all refused: $k = 50655/10000$ (fails at
$R = 7.9336$, as the paper's Remark 5.3 predicts), the printed bound moved to 2.1241, a wrong primitive. Not checked: the
geometric lemmas of Sections 2 to 4, which we read and re-derived by hand. Reproduce: `cd c22a && ./fetch.sh &&
./run_all.sh` (about ten seconds).

## 5b: exact values of f(n) for Erdős problem #757 (10 November 2026)

`c5b/c5b_check.py`, standard library only: the exact values of [SK2026] on the 5b page
([#253](https://github.com/teorth/optimizationproblems/pull/253), merged 26 Oct), $f(n) = 4, 5, 5, 6, 6, 7, 7, 8, 8, 9$ for
$6 \le n \le 15$ and $f(17) = 10$, $f(19) = 11$, $f(21) = 12$, from the thirteen data files
(JeetuSK0808/tao-constants-certificates @ `afa2e10`, `c5b-erdos757-small-n/certificates/n6.txt` to `n21.txt`, each
hash as the repository's SHA256SUMS states); its `verify.py` and `check_all.sh` were neither read nor run. Each file is
checked to be a $(4,5)$-set of $n$ distinct integers by two routes (every 4-subset and its six differences; and every
pair of pairs with equal difference, attributed to the 4-sets that contain it), and its largest Sidon subset $h(A)$ is
computed by two routes (branch and bound with the difference test; and brute force over all subsets of sizes $h$ and
$h+1$ with the sum test).

| $n$ | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 17 | 19 | 21 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| $(4,5)$-set, both routes | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes | yes |
| $h(A)$, both routes | 4 | 5 | 5 | 6 | 6 | 7 | 7 | 8 | 8 | 9 | 10 | 11 | 12 |
| $n - \lfloor (8n-6)/17 \rfloor$ | 4 | 5 | 5 | 6 | 6 | 7 | 7 | 8 | 8 | 9 | 10 | 11 | 12 |

So each set attains the lower bound, and every value in the table holds given that bound. Red controls, all as
expected: the largest element of the $n = 21$ set replaced by $a_1 + a_2 - a_0$ (a parallelogram, two coincidences in one
4-set) is refused by both $(4,5)$ routes; a Sidon set of 21 powers of two ($h = 21$) refutes the claim; a repeated
element is refused; every file with $h$ claimed one larger is refused; and the routes agree on 200 random small sets.
What rests on mathematics: the lower bound $f(n) \ge n - \lfloor (8n-6)/17 \rfloor$ (Ma and Tang, arXiv:2602.23282,
Lemma 2.4 and Theorem 5.3 of Henning and Yeo, with the $-6$ kept), which we did not check. Not checked: the searches at
$n = 16, 18, 20, 22, 23, 25, 30$. Reproduce: `cd c5b && ./fetch.sh && ./run_all.sh` (about ten seconds).

## Reproduce

```
./fetch.sh      # downloads the certificates at the pinned commits and checks every SHA-256
./run_all.sh    # writes results/*.txt (Python 3.9+, standard library only; about eight minutes)
```

Prepared with AI assistance (Claude, Anthropic), as the repository's AI use policy asks. Questions or corrections:
open an issue here. Arhan Canli, Quantavane (quantavane.com).
