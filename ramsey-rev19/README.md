# Ramsey lower bounds since Small Ramsey Numbers revision 18 (4 November 2026)

Independent exact checks of the new lower bounds for classical two-colour Ramsey numbers R(s,t) published since
revision 18 (24 April 2026) of Radziszowski's dynamic survey *Small Ramsey Numbers* (Electronic Journal of
Combinatorics, DS1). Every witness comes from its authors' own channel (repository at a pinned commit, or the paper's
data record), its SHA-256 is pinned in `SHA256SUMS`, and it was checked by two exact programs written apart from the
searches that found the colourings:

- `checkers/ramsey_check.c`: Bron-Kerbosch with Tomita's pivot and a greedy-colouring bound. It decides whether a red
  K_s or a blue K_t exists, over all vertices, using no symmetry.
- `checkers/ramsey_indep.c`: branch and bound with greedy-colouring bounds (MCQ). It finds the largest red and the
  largest blue clique exactly. For a circulant matrix it searches from vertex 0 only (rotations preserve both
  colours); for any other matrix, over all vertices.

Both are copied unchanged from our engine (commit ed11dfd); above 1,024 vertices `ramsey_indep` is compiled with room
for 2,048 (two `#define`s, done by `run_all.sh`). A colouring passes when `ramsey_check` finds neither clique and
`ramsey_indep`'s largest cliques are below s and t. `to_matrix.py` turns each published format into n lines of n
characters, `1` for the colour that must avoid K_s, and refuses anything not square, 0/1, symmetric with a zero
diagonal on exactly n vertices. Every script that reads a downloaded file runs under `python3 -I`.

## Results

| Bound | Rev. 18 | Source | File SHA-256 | n | ramsey_check | ramsey_indep (largest red / blue) | Result |
| --- | --- | --- | --- | --- | --- | --- | --- |
| R(3,17) >= 93 | 92 | NaRT repository (22 Oct 2026) | `71810667...336b` | 92 | no red, no blue clique, 0.5 s | largest 2 / 16, 0.3 s | PASS |
| R(3,24) >= 151 | 143 | Coniglio et al., arXiv:2608.18769 (19 Aug 2026); repository | `95d80b57...100f` | 150 | from vertex 0: N_red(0) no red K2, 0.0 s; N_blue(0) no blue K23, 42.7 s; all vertices: stopped at the 300 s cap | largest 2 / 23 (from vertex 0), 29.4 s | PASS |
| R(3,25) >= 160 | 154 | Coniglio et al., arXiv:2608.18769 (19 Aug 2026); repository | `313ce54a...a67f` | 159 | from vertex 0: N_red(0) no red K2, 0.0 s; N_blue(0) no blue K24, 42.6 s | largest 2 / 24 (from vertex 0), 44.5 s | PASS |
| R(3,26) >= 163 | 161 | Coniglio et al., arXiv:2608.18769 (19 Aug 2026); repository | `2cac81e0...e9e8` | 162 | from vertex 0: N_red(0) no red K2, 0.0 s; N_blue(0) no blue K25, 10.1 s | largest 2 / 25 (from vertex 0), 5.6 s | PASS |
| R(3,28) >= 181 | 179 | Coniglio et al., arXiv:2608.18769 (19 Aug 2026); repository | `f6c9ff19...f6bb` | 180 | from vertex 0: N_red(0) no red K2, 0.0 s; N_blue(0) no blue K27, 27.2 s | largest 2 / 27 (from vertex 0), 144.2 s | PASS |
| R(3,29) >= 194 | 190 | Coniglio et al., arXiv:2608.18769 (19 Aug 2026); repository | `31dc46a6...eb7b` | 193 | from vertex 0: N_red(0) no red K2, 0.0 s; N_blue(0) no blue K28, 203.1 s | largest 2 / 28 (from vertex 0), 354.1 s | PASS |
| R(3,30) >= 199 | 197 | Coniglio et al., arXiv:2608.18769 (19 Aug 2026); repository | `a445a432...e08e` | 198 | from vertex 0: N_red(0) no red K2, 0.0 s; N_blue(0) no blue K29, 60.4 s | largest 2 / 29 (from vertex 0), 150.2 s | PASS |
| R(3,31) >= 210 | 208 | Coniglio et al., arXiv:2608.18769 (19 Aug 2026); repository | `0da33eab...1f81` | 209 | from vertex 0: N_red(0) no red K2, 0.0 s; N_blue(0) no blue K30, 231.2 s | largest 2 / 30 (from vertex 0), 134.0 s | PASS |
| R(3,32) >= 219 | 217 | Coniglio et al., arXiv:2608.18769 (19 Aug 2026); repository | `9ef8db1d...54ae` | 218 | from vertex 0: N_red(0) no red K2, 0.0 s; N_blue(0) no blue K31, 243.2 s | largest 2 / 31 (from vertex 0), 615.4 s | PASS |
| R(3,33) >= 230 | 227 | Coniglio et al., arXiv:2608.18769 (19 Aug 2026); repository | `a41c978a...8d76` | 229 | from vertex 0: N_red(0) no red K2, 0.0 s; N_blue(0) no blue K32, 838.2 s | stopped at the 2400 s cap | not finished |
| R(3,34) >= 239 | 234 | Coniglio et al., arXiv:2608.18769 (19 Aug 2026); repository | `21cc082c...dd0c` | 238 | from vertex 0: N_red(0) no red K2, 0.0 s; N_blue(0) stopped by us after 589 s (budget) | stopped by us after 589 s (budget) | not finished |
| R(3,35) >= 253 | 248 | Coniglio et al., arXiv:2608.18769 (19 Aug 2026); repository | `265c8af6...b03b` | 252 | not run | not run | not finished |
| R(3,36) >= 261 | 255 | Coniglio et al., arXiv:2608.18769 (19 Aug 2026); repository | `04ecd5b0...8bb9` | 260 | not run | not run | not finished |
| R(3,37) >= 272 | 267 | Coniglio et al., arXiv:2608.18769 (19 Aug 2026); repository | `ca6931f5...11b4` | 271 | not run | not run | not finished |
| R(3,38) >= 281 | 278 | Coniglio et al., arXiv:2608.18769 (19 Aug 2026); repository | `b2151ed3...3990` | 280 | not run | not run | not finished |
| R(3,39) >= 294 | 290 | Coniglio et al., arXiv:2608.18769 (19 Aug 2026); repository | `57f83e87...f02d` | 293 | not run | not run | not finished |
| R(3,40) >= 305 | 298 | Coniglio et al., arXiv:2608.18769 (19 Aug 2026); repository | `64f89a70...1343` | 304 | not run | not run | not finished |
| R(3,41) >= 317 | 311 | Coniglio et al., arXiv:2608.18769 (19 Aug 2026); repository | `5898ffb6...7f50` | 316 | not run | not run | not finished |
| R(3,42) >= 327 | 320 | Coniglio et al., arXiv:2608.18769 (19 Aug 2026); repository | `0721f315...a6df` | 326 | not run | not run | not finished |
| R(3,43) >= 340 | 333 | Coniglio et al., arXiv:2608.18769 (19 Aug 2026); repository | `a84e6032...5aec` | 339 | not run | not run | not finished |
| R(3,44) >= 348 | 339 | Coniglio et al., arXiv:2608.18769 (19 Aug 2026); repository | `33c10272...abed` | 347 | not run | not run | not finished |
| R(3,45) >= 364 | 354 | Coniglio et al., arXiv:2608.18769 (19 Aug 2026); repository | `3db3bf27...b88f` | 363 | not run | not run | not finished |
| R(3,46) >= 371 | 362 | Coniglio et al., arXiv:2608.18769 (19 Aug 2026); repository | `8d7cbebd...6495` | 370 | not run | not run | not finished |
| R(3,47) >= 385 | 380 | Coniglio et al., arXiv:2608.18769 (19 Aug 2026); repository | `7781ebdc...5797` | 384 | not run | not run | not finished |
| R(3,48) >= 391 | 384 | Coniglio et al., arXiv:2608.18769 (19 Aug 2026); repository | `0a7502d4...9e9c` | 390 | not run | not run | not finished |
| R(3,49) >= 411 | 402 | Coniglio et al., arXiv:2608.18769 (19 Aug 2026); repository | `4142c5ce...9aaf` | 410 | not run | not run | not finished |
| R(4,15) >= 163 | 159 | NaRT repository (22 Oct 2026) | `56f3cf76...ec23` | 162 | no red, no blue clique, 0.6 s | largest 3 / 14, 0.7 s | PASS |
| R(4,16) >= 174 | 170 | NaRT, arXiv:2603.09172 v3-v5; repository (25 Mar 2026) | `def3d53e...2f6f` | 173 | no red, no blue clique, 2.2 s | largest 3 / 15 (from vertex 0), 0.1 s; renumbered copy, all vertices: largest 3 / 15, 2.2 s | PASS |
| R(4,18) >= 209 | 205 | NaRT, arXiv:2603.09172 v3-v5; repository (2 Apr 2026) | `5317e46e...9945` | 208 | no red, no blue clique, 6.3 s | largest 3 / 17, 4.1 s | PASS |
| R(4,19) >= 219 | 213 | NaRT, arXiv:2603.09172 v3-v5; repository (19 Apr 2026) | `8d237842...0f2c` | 218 | no red, no blue clique, 9.6 s | largest 3 / 18 (from vertex 0), 0.4 s; renumbered copy, all vertices: largest 3 / 18, 8.0 s | PASS |
| R(4,20) >= 237 | 234 | NaRT, arXiv:2603.09172 v3-v5; repository (19 Apr 2026) | `d33838ca...9e3e` | 236 | no red, no blue clique, 6.0 s | largest 3 / 19 (from vertex 0), 0.3 s; renumbered copy, all vertices: largest 3 / 19, 5.5 s | PASS |
| R(4,20) >= 252 | 234 | Yu, arXiv:2608.18169 (17 Aug 2026); Zenodo 21941603 | `cfc65790...712a` | 251 | no red, no blue clique, 36.9 s | largest 3 / 19 (from vertex 0), 1.9 s; renumbered copy, all vertices: largest 3 / 19, 33.0 s | PASS |
| R(4,20) >= 252 (G') | 234 | Yu, same paper: the copy x -> 116x, built from the same certificate | `cfc65790...712a` | 251 | no red, no blue clique, 33.4 s | largest 3 / 19 (from vertex 0), 2.1 s | PASS |
| R(4,21) >= 244 | 242 | NaRT repository (11 May 2026) | `3cdf3083...599e` | 243 | no red, no blue clique, 15.6 s | largest 3 / 20, 9.9 s | PASS |
| R(5,15) >= 279 | 275 | Canli, Zenodo 23857057 (1 Nov 2026); quantavane-ramsey-records | `489cb832...4b10` | 278 | no red, no blue clique, 22.9 s | largest 4 / 14 (from vertex 0), 0.3 s; renumbered copy, all vertices: largest 4 / 14, 11.4 s | PASS |
| R(6,8) >= 135 | 134 | Cremer, arXiv:2610.12122 (8 Oct 2026); repository v1.1.0 | `8ef378ad...c982` | 134 | no red, no blue clique, 0.0 s | largest 5 / 7, 0.0 s | PASS |
| R(7,13) >= 513 | 511 | NaRT repository (22 Oct 2026) | `e8160ebe...d70c` | 512 | stopped at the 120 s cap | not run | not finished |
| R(8,10) >= 345 | 343 | Cremer, arXiv:2610.12122; repository v1.1.0, witness 1 | `17881b78...c21d` | 344 | no red, no blue clique, 3.1 s | largest 7 / 9, 2.2 s | PASS |
| R(8,10) >= 345 | 343 | Cremer, arXiv:2610.12122; repository v1.1.0, witness 2 | `fd3c5a55...ad89` | 344 | no red, no blue clique, 3.1 s | largest 7 / 9, 2.1 s | PASS |
| R(8,12) >= 511 | - (471 implied) | NaRT repository (22 Oct 2026) | `5db94f26...a2df` | 510 | stopped at the 120 s cap | not run | not finished |
| R(8,12) >= 519 | - (471 implied) | Canli, Zenodo 23920916 (3 Nov 2026); quantavane-ramsey-records | `d84b5768...2407` | 518 | 3 Nov 2026 run, same file: no red, no blue clique, 39836070922 nodes, 8455 s wall | 3 Nov 2026 run: largest 7 / 11 | PASS (not re-run) |
| R(8,18) >= 1052 | 1049 | NaRT repository (22 Oct 2026) | `c5a4706f...f1ad` | 1051 | stopped at the 120 s cap | not run | not finished |
| R(9,11) >= 617 | - (608 implied) | NaRT repository (22 Oct 2026) | `ce12b6cb...afa1` | 616 | stopped at the 120 s cap | not run | not finished |
| R(10,12) >= 856 | - (846 implied) | NaRT repository (22 Oct 2026) | `4d270464...79e2` | 855 | stopped at the 120 s cap | not run | not finished |
| R(12,12) >= 1641 | 1640 | Tatarevic, arXiv:2608.06531 (6 Aug 2026); repository | `de07986c...1088` | 1640 | stopped at the 120 s cap | not run | not finished |

## Notes

- Circulant witnesses. `ramsey_indep` uses the rotations. For the five circulants outside the R(3,n) set we also ran it
  on a renumbered copy (`relabel.py`, seeded permutation), which it searches over all vertices: same clique numbers.
  For the R(3,n) certificates of Coniglio et al. the plain `ramsey_check` did not finish R(3,24) in 300 CPU seconds,
  so it was run on the two vertex-0 neighbourhoods instead (`nbhd.py`): a red K_s exists iff N_red(0) holds a red
  K_(s-1), and likewise in blue, because rotations act transitively and preserve both colours. Both programs then
  rest on that symmetry for those rows; `to_matrix.py`, `nbhd.py` and `ramsey_indep` each confirm the matrix is
  circulant.
- Yu's paper names a second graph, the image of the first under x -> 116x (mod 251); we built it from the paper's
  multiplier, confirmed the map is an isomorphism and its connection set lies in D0 u D3, and checked it too.
- Not finished here, to stay within about 4 core-hours: the second program on R(3,33), the certificates R(3,34) to
  R(3,49) of Coniglio et al., NaRT's R(7,13), R(8,12) >= 511, R(8,18), R(9,11), R(10,12), and Tatarevic's R(12,12).
  For R(3,n) the cost of both programs together grew from 16-557 s (n = 24-31) to 858 s (n = 32); on R(3,33)
  `ramsey_check` passed in 838 s but `ramsey_indep` stopped at its 2,400 s cap, and R(3,34) was stopped after 589 s
  per program. A fit to these runs puts the seventeen unfinished certificates at roughly 1,000 core-hours, most of it
  in the largest few (the paper reports several days of stand-alone verification for the largest graphs). Each large
  NaRT or Tatarevic witness ran 120 CPU seconds of `ramsey_check` without finishing; the nearest case we have measured
  to the end, R(8,12) on 511 and 518 vertices, took about two hours per program, so each needs at least that, and the
  1,051- and 1,640-vertex ones far more.
- Our own R(8,12) >= 519 colouring was checked to the end by both programs on 3 November 2026 (same file, same
  SHA-256; `results/qv_R8_12_n518.2026-11-03.*.txt`) and was not run again.
- Bounds that follow from these by the survey's own recursions: R(4,21) >= 257 (2.3.e from R(4,20) >= 252) and
  R(7,16) >= 649 (2.3.i from R(4,15) >= 163). Revision 18 also gives R(4,15) as 159 in Table Ia but 158 in
  Table IIa, and its own R(4,15) >= 159 already gives R(7,16) >= 633 by 2.3.i, where Table IIa has 629.

## Red controls

`controls.py` makes two controls per witness, each flipping exactly one edge: it finds a clique C of size s-2 with at
least two common neighbours in that colour, which are then joined in the other colour, and flips that edge, so
C and the two neighbours form a K_s (likewise in blue with t). The planted clique is re-checked pair by pair. Both
programs must report it.

| Witness | Control | ramsey_check | ramsey_indep (largest red / blue) |
| --- | --- | --- | --- |
| cremer_R6_8_n134 | one red edge made blue | red none, blue found, exit 1: fails as it should | largest 5 / 8: fails as it should |
| cremer_R6_8_n134 | one blue edge made red | red found, blue none, exit 1: fails as it should | largest 6 / 7: fails as it should |
| cremer_R8_10_n344_1 | one red edge made blue | red none, blue found, exit 1: fails as it should | largest 7 / 10: fails as it should |
| cremer_R8_10_n344_1 | one blue edge made red | red found, blue none, exit 1: fails as it should | largest 8 / 9: fails as it should |
| nart_R3_17_n92 | one red edge made blue | red none, blue found, exit 1: fails as it should | largest 2 / 17: fails as it should |
| nart_R3_17_n92 | one blue edge made red | red found, blue none, exit 1: fails as it should | largest 3 / 16: fails as it should |
| nart_R4_15_n162 | one red edge made blue | red none, blue found, exit 1: fails as it should | largest 3 / 15: fails as it should |
| nart_R4_15_n162 | one blue edge made red | red found, blue none, exit 1: fails as it should | largest 4 / 14: fails as it should |
| nart_R4_16_n173 | one red edge made blue | red none, blue found, exit 1: fails as it should | largest 3 / 16: fails as it should |
| nart_R4_16_n173 | one blue edge made red | red found, blue none, exit 1: fails as it should | largest 4 / 15: fails as it should |
| qv_R5_15_n278 | one red edge made blue | red none, blue found, exit 1: fails as it should | largest 4 / 15: fails as it should |
| qv_R5_15_n278 | one blue edge made red | red found, blue none, exit 1: fails as it should | largest 5 / 14: fails as it should |
| yu_R4_20_n251 | one red edge made blue | red none, blue found, exit 1: fails as it should | largest 3 / 20: fails as it should |
| yu_R4_20_n251 | one blue edge made red | red found, blue none, exit 1: fails as it should | largest 4 / 19: fails as it should |

## Reproduce

```
./fetch.sh ROOT      # the witness files at the pinned commits, SHA-256 checked
./run_all.sh ROOT    # converts, checks (cheap tier: under an hour of CPU here), writes results/
./run_all.sh ROOT all   # also the large witnesses, uncapped: many core-days
```

Python 3.9+ (standard library) and a C compiler. Arhan Canli, Quantavane (quantavane.com).
