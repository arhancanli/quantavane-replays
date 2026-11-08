# Our own improvements

Records of our own (not replays) for constants of the
[optimization constants repository](https://github.com/teorth/optimizationproblems). Each folder holds the
certificate, our exact checker (Python standard library only), its output with red controls, and a README with the
re-derive command and the credit for the warm start.

| Constant | New bound | Previous best (open pull request) | Folder |
| --- | --- | --- | --- |
| $C_{1a}$ (upper) | 1.5027427375 (exactly 1.50274273747657840978...) | 1.5027433522 (#270) | [`c1a`](c1a) |
| $C_{1c}$ (lower) | 0.4103598711 (exactly 0.41035987114651620459...) | 0.4103596512 (#270) | [`c1c`](c1c) |
| $C_{1d}$ (upper) | 1.454683534422819 (exactly 1.45468353442281856471...) | 1.454727664528751 (#312) | [`c1d`](c1d) |

Pull requests [#344](https://github.com/teorth/optimizationproblems/pull/344) (C1c) and
[#345](https://github.com/teorth/optimizationproblems/pull/345) (C1d) were opened on 8 Nov with the first certificates
(commit 36b3a05: 0.4103597993 and 1.454691416385622); the certificates here replace them.

Search code: [`search`](search). Prepared with AI assistance (Claude, Anthropic).
