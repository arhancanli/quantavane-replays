# 10a, the real Grothendieck constant: the Davie-Reeds value K_DR and four closed forms (15 November 2026)

Values recomputed from their definitions, not certificates. No program of any source was read or run; read as text:
Jones and Malavolta, arXiv 2603.30039 (e-print `2274c6d9…6da6`), Lemma unique-c and the Davie-Reeds lemma.

**[Dav1984], [Ree1991] $K_{DR} \approx 1.67696\ldots$** (lower bound). As stated there: $C^*$ is the unique root in $(0,1)$
of $H(C) = 4\varphi(C)^2 - 4\Phi(-C) + 1$ ($H' = 4\varphi(1 - 2C\varphi) > 0$, as $2C\varphi(C) < 2\varphi(0) < 1$),
$\lambda^* = 2C^*\varphi(C^*)$, and $K_{DR} = 1/R(C^*)$ with
$R(C) = \frac{4\varphi^2 - 2C\varphi(4\Phi(-C) - 1)}{1 - 2C\varphi}$.

- `c10a_a.py`: integer fixed point ($2^{-200}$) with its own exp, log, sqrt, $\pi$ and the erf series; $C^*$ by
  bisection with exact dyadic points; at the root $4\Phi(-C^*) - 1 = 4\varphi(C^*)^2$, so $K_{DR} = (\pi/2)e^{C^{*2}}$.
- `c10a_b.py`: decimal intervals (60 digits; $\Phi$ by its convergent series at 110 digits); $C^*$ by bisection;
  $K_{DR} = 1/R(C)$ over the final $C$-interval with the full formula (no use of the identity).

| quantity | both programs (first 20 digits identical, checked by `run_all.sh`) |
|---|---|
| $H(0)$, $H(1)$ | $-0.363380227632$, $0.599578310371$ (the lemma's $H(1) \approx 0.59958$) |
| $C^*$ | $0.255730213166210014361869933728\ldots$ |
| $\lambda^*$ | $0.197479090994\ldots$ |
| $K_{DR}$ | $1.676956674215576237077855078853\ldots$, so the page's $1.67696\ldots$ holds |

**Closed forms**, each enclosed and checked against the page's decimal (within half a unit of its last digit):
[G1953] $\sinh(\pi/2) = 2.30129890230729487346\ldots \approx 2.30130$ (upper) and $\pi/2 = 1.57079632679\ldots \approx 1.57080$
(lower); [K1979] $\pi/(2\log(1+\sqrt2)) = 1.78221397819136911177\ldots \approx 1.782214$; [SLXCKKM26]
$6\pi/11 = 1.71359599286715994825\ldots \approx 1.71360$.

**Red controls**, each refused by both programs: $K_{DR}$ printed as 1.67695; $\Phi(C)$ in place of $\Phi(-C)$ in the
equation for $C^*$ (no root in $(0,1)$); $\sinh(\pi/2)$ printed as 2.30131.

**Not checked:** the games behind the bounds (the Davie-Reeds operator's SDP value, the theorems of [Hei26], [JM26],
[SLXCKKM26], [BMMN2011], Krivine's and Grothendieck's arguments); the computer-assisted row [P2026] 1.7813\*.

Reproduce: `./run_all.sh` (Python 3.9+, standard library only; under a minute). Output: `results/c10a.txt`.
