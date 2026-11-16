# 50a: Lee and Parekh's 0.595 [LP24] (17 November 2026)

Lee and Parekh (arXiv 2401.03616, Theorem 2; read through the HTML rendering): the better of the product state and the
matching state has ratio at least $\max_p \min_{x\in[-1,1/3)} \rho_p(x)$ with
$\rho_p(x) = [p(1 - F_3(x)) + (1-p)(1 + \frac65\max(0, -(1+3x)))]/(1-3x)$ ($x = v_i\cdot v_j$; the $\frac65$ is Lemma 4's
$\frac85$ scaling of $h^+$ times $\frac34$), stated as 0.595 at $p = 0.674$. The theorem is read, not checked. We certify
$\min_x \rho_{0.674}(x) \ge 0.595$ with two separately written programs, our [GP19] programs of folder `c50a`
(hypergeometric series in exact rationals; closed forms through arcsin in decimal intervals) with the ratio replaced:
the numerator decreases in $x$ and the denominator is largest at a cell's left end, so each cell $[a,b]$ gives
$\rho \ge \mathrm{num}(b)/(1-3a)$; on $[0, 1/3)$, $\rho \ge \mathrm{num}(1/3) \ge 0.807$.

Both: every cell of $[-1, 0]$ clears 0.595 (the least lower bounds 0.59534 and 0.59553, near $x \approx -0.457$; the
minimum is $0.59572\ldots$), so $\max_p\min_x\rho_p \ge 0.595$. We use the $F_3$ of [GP19] and Briët, de Oliveira Filho
and Vallentin, $F_3(x) = \frac{8}{3\pi}x\,{}_2F_1(\frac12,\frac12;\frac52;x^2)$. Red controls refused: the claim 0.5958;
$p = 0.9$ and $p = 0.5$.

Reproduce: `./run_all.sh` (Python 3.9+; under a minute). Output: `results/c50alp.txt`.
