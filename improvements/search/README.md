# Search code for the C1a, C1c and C1d improvements

The searches ran on 8 Nov 2026: two Kaggle sessions on 2 x Tesla T4 (about 2.6 h each), then one CPU core per job.
Each Kaggle session ran `launcher_template.py` with a worker module pasted into `SRC` (job list `c,dA,dB,dC` for the
first session, `c1,c2,c3` for the second); the launcher writes the module to `/tmp` and starts one process per job.
The CPU runs used the `worker_local_*.py` modules directly (`python3 worker_local_1c.py Lc1`, with `C1_TOTAL` the time
budget in seconds), each started from the previous best integer certificate.

| constant | step | file | job | result |
| --- | --- | --- | --- | --- |
| C1c | 1 (GPU) | `worker_kernel2_1c.py` | `c1`: soft-min cycles from #270's `c1c_v3`, beta 3e7..1e10, 20,000 L-BFGS-B iterations per stage | 0.4103597993 |
| C1c | 2 (CPU) | `worker_local_1c.py` | `Lc1`: two more cycles from step 1 | 0.4103598711, the certificate |
| C1d | 1 (GPU) | `worker_kernel1_1d.py` | `dB`: #312's certificate split to 163,840 cells, polish at p = 1e8..1e11, restart of the p-norm continuation from p = 32 | 1.4546914164 |
| C1d | 2 (CPU) | `worker_local_1d.py` | `LdB`: restarts from step 1 at p = 32, 16, 64 | 1.4546835344, the certificate |
| C1d | 3 (CPU, 9 Nov) | `worker_cont_1d_numpy.py` | 14 restarts from the 8 Nov certificate, p0 cycling through 8, 128, 24, 48, 12, 96 (numpy and scipy; L-BFGS on log ||v*v||_p - 2 log S) | 1.4546733474, the certificate |
| C1a | 1 (CPU) | `worker_local_1a.py` | `aA`: #270's `c1a_n262144` polished at p = 1e10..1e12, then four cycles of restarts at p0 = 1e8, 1e7, 1e9, 1e6 | 1.5027427375, the certificate |

The other jobs ended with worse values: C1c from #318's certificate reached 0.4103539040; C1c from `c1c_v3` split to
K = 196608 reached 0.4103597338, and at K = 98304 from beta = 1e9 reached 0.4103597183; C1d jobs `dA` (163,840 cells)
and `dC` (81,920 cells) reached 1.4546929051 and 1.4546997322.

Warm starts (published certificate data only): `vidigoat_c1c_K98304_v3.json` (#270, sha256 `e0786a59...caec66`),
`vidigoat_c1a_n262144.json` (#270, sha256 `1848a920...f3de`), `julsgaard_certificate.json` (#312, sha256
`04c9c759...f7237`), `julsgaard_c1c_K98304.json` (#318, sha256 `fe11336b...52d2`). The workers write each new best as
integers (largest absolute value about 2^56); exactness is the job of `check.py` in each certificate folder, never of
the search.

C1c, 10 Nov: `c1c_spike/` (numpy and scipy, one CPU core): `fp.py` solves the spike-plus-profile equations of PR #354's
structure (one spike cell of weight M, every lag from t1 K to K active) by a damped fixed-point iteration, `opt.py`
chooses M by golden-section search, `ladder.py` runs K = 2^14, 2^16, 2^18 (about 45 s in all), and `mkcert.py`
scales the weights to 2^60, rounds them and writes the certificate. Written from #354's description; its programs were
not read or run.
The K = 2^20 certificate of 10 Nov (evening) is `python ladder.py 1048576` (858 s on one core, M = 641.639118)
followed by `python mkcert.py 1048576 OUT.json`; nothing in these files changed.

C1b, 10 Nov: `c1b/` (numpy and scipy, one core): `c1b_opt.py` minimises log ||C||_p over the shifts
(C_k = sum_i h_i (1 - h_{i+k}), 0 <= h <= 1) by L-BFGS-B with a quadratic penalty on the mass, restoring sum h = n/2
exactly after each stage; `python c1b_opt.py construction.json out.json --ps 3e9,1e10,3e10,1e11 --iters 4000`
(about 10 minutes) from the Atris Labs construction of the 1b page ([A2026], data only). `mkcert_1b.py out.json
cert.json` rounds the heights to multiples of 2^-48 and fixes the mass exactly in integers.
