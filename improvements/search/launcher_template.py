#!/usr/bin/env python3
# Kaggle launcher: writes the worker module (SRC) to /tmp and runs one worker per job (1c on GPU 0, three 1d
# searches spread over both GPUs); every output is in /kaggle/working.
import os, subprocess, sys
SRC = r"""%s"""
os.makedirs("/tmp/c1", exist_ok=True)
open("/tmp/c1/c1job.py", "w").write(SRC)
env = dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
names = os.environ.get("C1_JOBS", "c,dA,dB,dC").split(",")
procs = [subprocess.Popen([sys.executable, "/tmp/c1/c1job.py", nm], env=env) for nm in names]
codes = [p.wait() for p in procs]
print("all done", codes, flush=True)
