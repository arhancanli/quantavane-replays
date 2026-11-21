# Our own improvements

Records of our own (not replays) for constants of the
[optimization constants repository](https://github.com/teorth/optimizationproblems). Each folder holds the
certificate, our exact checker (Python standard library only), its output with red controls, and a README with the
re-derive command and the credit for the warm start.

| Constant | New bound | Previous best (open pull request) | Folder |
| --- | --- | --- | --- |
| $C_{1a}$ (upper) | 1.5027427280 (exactly 1.50274272798999968357...; 10 Nov evening, superseding 1.5027427375 of 8 Nov) | 1.5027427341 (#270, vidigoat, 10 Nov) | [`c1a`](c1a) |
| $C_{1b}$ (upper) | 0.3808585498 (exactly 0.38085854975166311714...; 10 Nov, 114,688-step function refined from [A2026]) | 0.3808585504 ([A2026] on the page; no open pull request) | [`c1b`](c1b) |
| $C_{1c}$ (lower) | 0.4103780000 (exactly 0.41037800003882762970...; 10 Nov evening, K = 2^20, superseding 0.4103774129 of 10 Nov morning and 0.4103598711 of 8 Nov) | 0.4103772 (#354, K = 524288) | [`c1c`](c1c) |
| $C_{1d}$ (upper) | 1.4513061918499 (exactly 1.45130619184989692...; 21 Nov, 102,400-cell function refined from C. Weill's #374, superseding our 1.4546733414 of 10 Nov) | 1.4513064006879361 (#374, C. Weill, 17 Nov) | [`c1d`](c1d) |
| $C_{12}$ (lower), the Euclidean TSP (BHH) constant | 0.64295 (exactly >= 0.6429567398; 21 Nov: Dong-Cao Theorem 3.11 with our sharper certified $C_3 \geq 0.0114906492$ and $C_4 \geq 0.0054506029$ and their $C_5$) | 0.6421913 (Z. Dong and J. Cao, arXiv:2610.02809; the page lists 0.6277) | [`c12_bhh_lower`](c12_bhh_lower) |

Pull requests [#344](https://github.com/teorth/optimizationproblems/pull/344) (C1c) and
[#345](https://github.com/teorth/optimizationproblems/pull/345) (C1d) were opened on 8 Nov with the first certificates
(commit 36b3a05: 0.4103597993 and 1.454691416385622); the certificates here replace them.

On 10 Nov vidigoat ([#270](https://github.com/teorth/optimizationproblems/pull/270)) refined our C1a and C1d certificates further, to C1a <= 1.5027427341 and C1d <= 1.4546733458 (crediting ours as the warm starts); the 10 Nov evening certificates above start from those and go 6.065e-9 and 4.382e-9 lower.

Search code: [`search`](search). Prepared with AI assistance (Claude, Anthropic).
