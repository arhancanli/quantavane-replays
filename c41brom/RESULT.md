# 41b: the area of Romik's ambidextrous sofa [Romik2016] (16 November 2026)

Romik (arXiv 1606.08111, abstract, eq. (5) and Theorem 5): the area is $X + \arctan Y$ with $x^2(x+3) = 8$ and
$x(4x^2+3) = 1$, equivalently $\sqrt[3]{3+2\sqrt2} + \sqrt[3]{3-2\sqrt2} - 1 + \arctan\bigl[\frac12(\sqrt[3]{\sqrt2+1} - \sqrt[3]{\sqrt2-1})\bigr]$.
That the sofa has this area is Romik's theorem (read, not checked). Two separately written programs:
`c41brom_a.py` (the cubics' real roots in exact rationals, arctan by halving and series) and `c41brom_b.py` (the
radical form in 60-digit decimals, arctan by Euler's series; the radicals checked against the cubics).

Both: area $= 1.64495521842544085166880934\ldots$; the printed 1.64495521 [Romik2016] (and 1.64495 [Gibbs2014]) are
truncations. Red controls, refused by each: the claim 1.64495522; a sign slip in the second cubic and $=9$ for $=8$
in the first (A); $3+2\sqrt3$ for $3+2\sqrt2$ and the arctan term dropped (B).

Reproduce: `./run_all.sh` (Python 3.9+; seconds). Output: `results/c41brom.txt`.
