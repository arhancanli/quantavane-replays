# c22tao: closed-form bounds, replay round of 22 November 2026

Source: constant pages of [teorth/optimizationproblems](https://github.com/teorth/optimizationproblems) at main 87b4ebb
(pages 1b, 11a, 11b, 13a, 13b, 15a, 22a, 22b, 24a, 35a, 47a, 82a; credit to the authors named in each row's reference).
Each row states a bound as a closed form with printed digits. `tao22_check.py` evaluates every closed form with two programs
written separately (A: decimal arithmetic, pi by Machin, logarithms by atanh series; B: integer fixed point at 10^250,
pi by Gauss-Legendre, logarithms by AGM, roots by isqrt), requires them to agree to 40 digits, requires the printed digits to be
within one unit of their last place of the exact value, and runs a red control per row (the printed value moved by two units must
fail). 22 rows, all PASS. Scope: this checks the arithmetic of the stated closed forms, not the papers' proofs that the closed form
is a valid bound. Where the page prints a rounding that lies on the other side of the exact value (the `printed_on_safe_side=False`
lines), the exact closed form is the bound; the printed decimal is only its "approximately" display.
