# Quantavane replays

Independent exact replays of certificates submitted to the
[optimization constants repository](https://github.com/teorth/optimizationproblems). Each replay uses our own code,
written from the definition on the constant page. The certificates' own verifiers were neither read nor run, and we
have no connection with any of the authors.

## Results (3 November 2026)

| Constant | Pull request | Certificate (pinned commit, SHA-256) | Exact value | Printed bound |
| --- | --- | --- | --- | --- |
| $C_{1c} \geq$ | [#292](https://github.com/teorth/optimizationproblems/pull/292) (Koshta's row) | DennisKoshta/c1c-autocorrelation-certificate @ `a4ee94ef29ff`, `K24576_0.4103333.json`, `963cb727...c9e` | 134648548624804778248790016 / 328144294688458468928334649 = 0.41033335274849333953..., minimum at k = 21264 only | 0.4103333, below the exact value |
| $C_{1d} \leq$ | [#312](https://github.com/teorth/optimizationproblems/pull/312) | Julsgaard/c1d-certificate @ `cc455ca0`, `cert/certificate.json` (81,920 cells), `04c9c759...237` | 9900630255081376929482866688 / 6805830738283663732311831361 = 1.45472766452875064884..., one maximising knot (m = 29124) | 1.454727664528751, above the exact value |
| $C_{1d} \leq$ | [#271](https://github.com/teorth/optimizationproblems/pull/271) | Keremozdemirra/tao-1c-1d-certificates @ `408f8d0b`, `c1d_n40960_1.4556234448.json` (40,960 cells), `be35060a...2da` | 279676466330630562685427712000 / 192135175720687825901095191121 = 1.45562344470022455281..., one maximising knot (m = 15313) | 1.4556234448, above |
| $C_{1d} \leq$ | [#269](https://github.com/teorth/optimizationproblems/pull/269) | echul-seal/c1d-upper-bound @ `2be42dfd`, `ub/c_f10240.json` from `c1d_records.tgz` (10,240 cells), `30d501e1...973` | 76358168487433178963290071040 / 52455219864513498497744894329 = 1.45568293650581519488..., one maximising knot (m = 3168) | 1.455682936505816, above |
| $C_{1d} \leq$ | [#286](https://github.com/teorth/optimizationproblems/pull/286) | huiihao/c1d-certificate @ `302817dc`, `construction.json` (400 cells, JSON floats), `c489b46b...8df` | read as binary64: 4563650116209721766004955898613711667200 / 3127595284199802099438492053245424408273 = 1.45915622115964902113...; read as the printed decimals: 1.45915622115964903535... | 1.4591562211596491, above both |

Every hash above is the one stated in the pull request (for #269 also the archive's, `b526326d...`).

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

## Reproduce

```
./fetch.sh      # downloads the certificates at the pinned commits and checks every SHA-256
./run_all.sh    # writes results/*.txt (Python 3.9+, standard library only; under a minute)
```

Prepared with AI assistance (Claude, Anthropic), as the repository's AI use policy asks. Questions or corrections:
open an issue here. Arhan Canli, Quantavane (quantavane.com).
