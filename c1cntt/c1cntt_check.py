"""Replay of PR #354's two latest C_1c certificates (SeanSagait/c1c-autocorrelation-certificate), K = 2^22 and
K = 2^24 cells, with our own exact NTT code (c1c_ntt.c). The repository's verify_ntt.py, verify_limbs.py, verify.py,
verify_independent.py and construction/ were neither read nor run; only the certificate text files are used.

The 1c page: C_1c is the least constant with min_{0<=t<=1} int f(x) f(x+t) dx <= C_1c (int f)^2 for all non-negative
f in L^1(R). For f = sum_i w_i 1_[i/K,(i+1)/K) the autocorrelation is linear between multiples of 1/K, so
C_1c >= K min_{0<=k<=K} A_k / S^2 with A_k = sum_i w_i w_{i+k}, S = sum_i w_i.

Usage: python3 -I c1cntt_check.py CERTS_ROOT BINARY [--controls]
"""
import hashlib
import os
import subprocess
import sys
from fractions import Fraction

ROWS = [
    # key, file (relative to CERTS_ROOT), sha256 of the text, stated fraction, stated lag, stated S, printed bound,
    # PR commit that states it, certificate repository commit, whether the printed bound is implied (our finding: the
    # K = 2^24 row prints 0.4103782, which is above its exact value 0.41037818308...; the expected outcome is recorded
    # here so that run_all.sh exits 0 exactly when every outcome is as found)
    ("k22", "k22/K4194304_0.4103781.txt", "d0278d18a330108b926281368b66aa41a069d7429451548b340bdb6366ca3494",
     Fraction(679700908027579442984780824576, 1656279491328950142501183554731), 3513048, 11219502722536334,
     "0.4103781", "2643147e48", "867e0e6f2858", True),
    ("k24", "k24/K16777216_0.4103782.txt", "f2b2e7771aa7f803c32efa91aacfd59abcc352719abfeb4281749df75d7c893a",
     Fraction(51659297209703959230693685329920, 125882172442308406613641411453441), 2691006, 11219722476171521,
     "0.4103782", "e06c3123e2", "f2b97a3caa1f", False),
]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run(binary, path, *extra):
    p = subprocess.run([binary, path, *extra], capture_output=True, text=True)
    out = {}
    for line in p.stdout.splitlines():
        parts = line.split()
        if parts:
            out[parts[0]] = parts[1:]
    return p.returncode, out, p.stdout


def dec(q, digits=20):
    n = q.numerator * 10 ** digits // q.denominator
    s = str(n).rjust(digits + 1, "0")
    return s[:-digits] + "." + s[-digits:]


def main():
    root, binary = sys.argv[1], sys.argv[2]
    controls = "--controls" in sys.argv[3:]
    ok = True
    cok = True
    findings = []
    for key, rel, want_sha, frac, lag, S_stated, printed, prc, certc, printed_ok in ROWS:
        path = os.path.join(root, rel)
        got = sha256(path)
        print(f"[{key}] {rel} (certificate commit {certc}; stated in #354 at {prc}) sha256 {got}: "
              f"as stated {got == want_sha}")
        ok &= got == want_sha
        rc, out, raw = run(binary, path)
        for line in raw.splitlines():
            print("  " + line)
        if rc != 0:
            print(f"  the exact computation failed (exit {rc})")
            ok = False
            continue
        K = int(out["K"][0])
        S = int(out["S"][0])
        Amin = int(out["Amin"][0])
        args = [int(a) for a in out["argmin"]]
        q = Fraction(K * Amin, S * S)
        print(f"  exact value K Amin / S^2 = {q.numerator}/{q.denominator} = {dec(q)}...")
        same = q == frac
        print(f"  equal to the stated fraction: {same}; minimum only at k = {args} (stated {lag}): "
              f"{args == [lag]}; S as stated: {S == S_stated}")
        b = Fraction(printed)
        nd = len(printed) - 2
        trunc = Fraction(q.numerator * 10 ** nd // q.denominator, 10 ** nd)
        if b <= q:
            nxt = b + Fraction(1, 10 ** nd)
            print(f"  printed {printed} <= exact value: True; the next printed value {dec(nxt, nd)} is above it: "
                  f"{nxt > q} (margin {float(q - b):.3e})")
            row_ok = nxt > q
        else:
            print(f"  printed {printed} <= exact value: False: the printed bound is {float(b - q):.3e} ABOVE the exact "
                  f"value, so it is NOT a valid lower bound as printed; the {nd}-decimal truncation {dec(trunc, nd)} is")
            findings.append(f"K = {K}: printed {printed} is not implied (exact {dec(q, 11)}...; {dec(trunc, nd)} is)")
            row_ok = True
        as_expected = (b <= q) == printed_ok
        print(f"  outcome as recorded in ROWS (printed bound implied: {printed_ok}): {as_expected}")
        ok &= same and args == [lag] and S == S_stated and row_ok and as_expected
        if controls:
            rc2, out2, _ = run(binary, path, "--bump", str(lag))
            q2 = Fraction(K * int(out2["Amin"][0]), int(out2["S"][0]) ** 2) if rc2 == 0 else None
            c = rc2 == 0 and q2 != frac
            print(f"  [ctl] w_{lag} + 1: the exact value becomes {dec(q2) if q2 else '?'}..., not the stated fraction, "
                  f"so the changed weights are refused: {c}")
            cok &= c
            rc3, out3, _ = run(binary, path, "--scale", "3")
            q3 = Fraction(K * int(out3["Amin"][0]), int(out3["S"][0]) ** 2) if rc3 == 0 else None
            c = rc3 == 0 and q3 == frac and [int(a) for a in out3["argmin"]] == [lag]
            print(f"  [ctl] every weight times 3 (values 9 times larger, same ratio): the same fraction and lag: {c}")
            cok &= c
            rc4, out4, _ = run(binary, path, "--fault")
            fails = int(out4.get("identity_failures", ["0"])[0])
            fl = out4.get("fault_lag", ["?"])[0]
            c = rc4 == 1 and fails >= 1 and out4.get("first_identity_failure", [None])[0] == fl
            print(f"  [ctl] one residue mod p1 corrupted (lag {fl}): the mod-p3 identity fails there "
                  f"({fails} failure), exit {rc4}: {c}")
            cok &= c
    if controls:
        p = subprocess.run([binary, "--selftest"], capture_output=True, text=True)
        print(p.stdout.strip())
        c = p.returncode == 0
        cok &= c
        print(f"controls all as expected: {cok}")
    res = ok and cok
    for f in findings:
        print("FINDING: " + f)
    print(f"RESULT: {'PASS' if res else 'FAIL'} (K = 2^22 and K = 2^24: each the stated fraction, lag and S exactly; "
          f"0.4103781 holds for both; 0.4103782 as printed for K = 2^24 does not)")
    return 0 if res else 1


if __name__ == "__main__":
    sys.exit(main())
