# 1c: the upper bound [BS2020] $1/(2(1+\theta_0))$ from its definition (16 November 2026)

The row: $C_{1c} \le \frac{1}{2(1+\theta_0)} = 0.4107675\ldots$, $\theta_0 = -\min_x \frac{\sin x}{x} = 0.2172336\ldots$
(Barnard and Steinerberger; the exact form is [MR2020, (1.2)]). The bound is their theorem (read, not checked); here
$\theta_0$ and the bound are computed from the definition by two separately written programs:

- `c1cbs20_a.py`: the critical point $x_0 \in (\pi, 3\pi/2)$ of $\sin x/x$ ($\sin x = x\cos x$) by bisection in exact
  rational intervals, $\theta_0 = 1/\sqrt{1+x_0^2}$, with the argument that this is the global minimum.
- `c1cbs20_b.py`: branch and bound on $\sin x/x$ over $[\pi/2, 4\pi]$ with a second-order Taylor bound, no use of the
  critical-point equation (80-digit decimals moved outward).

| quantity | value |
|---|---|
| $x_0$ | $4.493409457909064175307880927280\ldots$ |
| $\theta_0$ | $0.21723362821122165740827932556\ldots$ (the page's $0.2172336\ldots$ is a truncation) |
| $1/(2(1+\theta_0))$ | $0.41076748818940533767047553314\ldots$ |

Note: the page prints the bound as $0.4107675\ldots$, but its leading digits are $0.41076748\ldots$; 0.4107675 is the
correct rounding upward (valid for an upper bound). Red controls, refused by each program: $\theta_0$ claimed as
0.2172337...; the bound claimed as 0.4107674; $\theta_0$ from the maximum (A) / the search on $[\pi/2, \pi]$ only (B).

Reproduce: `./run_all.sh` (Python 3.9+; under a minute). Output: `results/c1cbs20.txt`.
