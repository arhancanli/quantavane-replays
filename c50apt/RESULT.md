# 50a: Parekh and Thompson's 0.533 [PT21] (17 November 2026)

Parekh and Thompson (arXiv 2105.05698, Theorem 2; read through the HTML rendering): the algorithm keeps the better of a
matching-based state (ratio 3/8 on small edges, $x_e \le 5/9 = \alpha(2)$, and 3/4 on large ones; Lemma 17, read) and the
Gharibian-Parekh product state (Lemma 18: at least 0.557931 on small edges and 0.498766 on large ones), and the proof
closes with a linear program in the share $s$ of small edges. The product-state numbers are minima of
$R(v) = (1 - F_3(v))/(1 - 3v)$, $F_3(v) = \frac{8}{3\pi} v\,{}_2F_1(\frac12,\frac12;\frac52;v^2)$, over $[-1, 1/3)$ and over
$[-5/9, 1/3)$ (the sign convention $x = -v$ of their Definition 5).

- `c50apt_a.py`, `c50apt_b.py`: our two GP19 programs of folder `c50a` (hypergeometric series in exact rationals;
  closed forms through arcsin in decimal intervals), copied with a range option (and adaptive cells in B):
  $\min_{[-1,1/3)} R \ge 0.498766$ and $\min_{[-5/9,1/3)} R \ge 0.557931$, both certified by both. (Numerically the
  minima are $0.4987674\ldots$ at $v \approx -0.966$ and $R(-5/9) = 0.5579321\ldots$.)
- `c50apt_lp.py`: the crossing of $\frac38 s + \frac34(1-s)$ and $0.557931 s + 0.498766(1-s)$:
  $s^* = 251234/434165$, $r^* = 231411/434165 = 0.5330024\ldots \ge 0.533$ (exact rationals, and again in decimals;
  $r^*$ increases with both constants, so the certified lower bounds suffice).

Red controls refused: the minima claimed as 0.498768 and 0.557933 (each by both programs); the bound claimed as 0.5331.

Reproduce: `./run_all.sh` (Python 3.9+; about a minute). Output: `results/c50apt.txt`.
