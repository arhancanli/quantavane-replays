# Our own improvements

Records of our own (not replays) for constants of the
[optimization constants repository](https://github.com/teorth/optimizationproblems). Each folder holds the
certificate, our exact checker (Python standard library only), its output with red controls, and a README with the
re-derive command and the credit for the warm start.

| Constant | New bound | Previous best (open pull request) | Folder |
| --- | --- | --- | --- |
| $C_{1a}$ (upper) | 1.5027427375 (exactly 1.50274273747657840978...) | 1.5027433522 (#270) | [`c1a`](c1a) |
| $C_{1c}$ (lower) | 0.4103780000 (exactly 0.41037800003882762970...; 10 Nov evening, K = 2^20, superseding 0.4103774129 of 10 Nov morning and 0.4103598711 of 8 Nov) | 0.4103772 (#354, K = 524288) | [`c1c`](c1c) |
| $C_{1d}$ (upper) | 1.454673347420953 (exactly 1.45467334742095245159...; 9 Nov, superseding 1.454683534422819 of 8 Nov) | 1.454727664528751 (#312) | [`c1d`](c1d) |

Pull requests [#344](https://github.com/teorth/optimizationproblems/pull/344) (C1c) and
[#345](https://github.com/teorth/optimizationproblems/pull/345) (C1d) were opened on 8 Nov with the first certificates
(commit 36b3a05: 0.4103597993 and 1.454691416385622); the certificates here replace them.

On 10 Nov vidigoat ([#270](https://github.com/teorth/optimizationproblems/pull/270)) refined our C1a and C1d certificates further, to C1a <= 1.5027427341 and C1d <= 1.4546733458 (crediting ours as the warm starts); those are now the best open claims for C1a and C1d, and our rows above stay as superseded valid bounds.

Search code: [`search`](search). Prepared with AI assistance (Claude, Anthropic).
