# Search code for the C1c and C1d improvements

Both certificates were found on Kaggle (two sessions on 2 x Tesla T4, 8 Nov 2026, about 2.6 h each). Each session
ran `launcher_template.py` with the worker module pasted into `SRC` (the job list set to `c,dA,dB,dC` for the first
session and `c1,c2,c3` for the second); the launcher writes the module to `/tmp` and starts one process per job.

| file | session | job that produced a certificate |
| --- | --- | --- |
| `worker_kernel1_1d.py` | 1 | `dB`: the C1d certificate `../c1d/c1d_upper_n163840.json` (split to 163,840 cells, polish at p = 1e8..1e11, restart of the p-norm continuation from p = 32) |
| `worker_kernel2_1c.py` | 2 | `c1`: the C1c certificate `../c1c/c1c_lower_K98304.json` (soft-min cycles, beta 3e7..1e10, 20,000 L-BFGS-B iterations per stage) |

The other jobs of the two sessions finished with worse values: C1c from #318's certificate reached 0.4103539040;
C1c from #270's certificate split to K = 196608 reached 0.4103597338, and at K = 98304 from beta = 1e9 reached
0.4103597183; C1d jobs `dA` (163,840 cells) and `dC` (81,920 cells) reached 1.4546929051 and 1.4546997322.

Warm starts (published certificate data only): `vidigoat_c1c_K98304_v3.json` (#270, sha256 `e0786a59...caec66`),
`julsgaard_certificate.json` (#312, sha256 `04c9c759...f7237`), `julsgaard_c1c_K98304.json` (#318, sha256
`fe11336b...52d2`). The worker writes each new best as integers (largest absolute value 2^56); exactness is the job of
`check.py` in each certificate folder, never of the search.
