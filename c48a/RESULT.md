# 48a, the one-dimensional convex sub-Gaussian comparison constant: $C_{48} = c_0^2$ (15 November 2026)

[DP2026] (D. Davis, S. Power, *The sharp one-dimensional convex sub-Gaussian comparison constant*, arXiv 2604.03170v1,
the only version; TeX source `paper_draft_normal_and_psi_1.tex`, SHA-256 `6ac7863e…7b07`, read as text) proves
$C_{48} = c_0^2$ with $c_0$ given by the system (3)-(6):

$t_0 = \sqrt{2\ln 2}$; $A = t_0 + 2\int_{t_0}^\infty e^{-t^2/2}\,dt$, $B = A/2$; $a > t_0$ the root of
$H(a) = 2a e^{-a^2/2} + 2\int_a^\infty e^{-t^2/2}\,dt = B$; $p_0 = 2e^{-a^2/2}$; $z$ the root of $\bar\Phi(z) = p_0$;
$c_0 = B/\varphi(z)$.

The page prints $C_{48} = c_0^2 \approx 5.33386$ in both its upper and its lower table. Two separately written
programs evaluate the system with rigorous enclosures; nothing is downloaded and no program of any source was read or
run.

- `c48a_a.py`: balls on exact rationals (centres rounded to $2^{-220}$, the rounding added to the radius);
  $Q(x) = \int_x^\infty e^{-t^2/2}\,dt = \sqrt{\pi/2} - \sum_n (-1)^n x^{2n+1}/(2^n n!(2n+1))$ (alternating);
  $\pi$ by Machin's formula; $\ln 2$ by $2\,\mathrm{atanh}(1/3)$; roots by bisection on rationals with both ends
  certified.
- `c48a_b.py`: decimal intervals at 80 digits, outward rounding;
  $Q(x) = \sqrt{\pi/2} - e^{-x^2/2}\sum_n x^{2n+1}/(1\cdot3\cdots(2n+1))$ (positive terms); $\pi$ by the
  Bailey-Borwein-Plouffe series; $\exp$, $\ln$ and $\sqrt{\ }$ from `decimal` (correctly rounded) widened by one
  unit in the last place; roots by interval bisection.

| quantity | both programs (enclosures overlap; digits truncated) | the paper's Remark 2 |
|---|---|---|
| $t_0$ | $1.177410022515474691011569326459\ldots$ | |
| $A$ | $1.776574120161349358934736714095\ldots$ | |
| $a$ | $1.803338389684026369643742273621\ldots$ | 1.80334 |
| $p_0$ | $0.393423950816028211896671736470\ldots$ | 0.39342 |
| $z$ | $0.270406029923698269232988940363\ldots$ | 0.27041 |
| $c_0$ | $2.309515867366166682601860767377\ldots$ | 2.30952 |
| $C_{48} = c_0^2$ | $5.333863541616097216005578611413\ldots$ | 5.33386 |

Each of the paper's five figures is its value rounded to five decimals, and so is the page's 5.33386.

**Red controls**, each refused by both programs (and both give the same changed value): $c_0^2$ claimed to round to
5.33387; $t_0 = \sqrt{\ln 2}$ ($C = 5.6275\ldots$); $p_0$ without its factor 2 ($C = 10.2706\ldots$); $H$ with
$\int_x^\infty e^{-t^2/2}\,dt$ once instead of twice ($C = 5.0247\ldots$).

**Not checked:** Theorem 1 itself (that $c_\star = c_0$), read.

Reproduce: `./run_all.sh` (Python 3.9+, standard library only; about a minute). Output: `results/c48a.txt`.
