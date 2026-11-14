#!/usr/bin/env python3
"""Program A for the 49a block of PR #370: our 4 November checker ../c49/sunflower_check.py, unchanged, applied to
the new certificate. This file only reads the certificate (JSON: n, ports, masks), registers the row's stated
numbers, and calls the checker's own functions: the exhaustive box method on two point partitions (core 8 + pendants
10 + 10, and 4 + 4 + 10 + 10), port counts, the transfer matrix, the Perron root as an exact surd, the rate enclosed
by exact integer comparisons, and the checker's red controls (a planted 3-sunflower comparable with no member, a
planted containment, a changed port count, a printed bound one unit higher).

Second path inside A (information): the certificate equals Koshta's block (replayed on 4 November, ../c49) minus the
21 sets and plus the 32 sets that best_delta.json lists, and every pair and triple through one of the 32 added sets
is scanned against all 252,303 members (the checker's pivot scan).

Usage: python3 -I c49p370_a.py CERTS_ROOT [--controls]   (one process)
"""
import importlib.util
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("sunflower_check", os.path.join(HERE, "..", "c49", "sunflower_check.py"))
sc = importlib.util.module_from_spec(spec)
sys.modules["sunflower_check"] = sc  # so that the checker's worker processes can be pickled by name
spec.loader.exec_module(sc)

CERT = "pr370/best_1.5678643845728069793620.json"
DELTA = "pr370/best_delta.json"
KOS = "kos/two_pendant_extra_1.5678599.txt"
STATED_SHA = "53381b5f52264c43d56343ab27d698c598462d5abfd83c043365afbae9efb1e8"  # the packet's SHA256SUMS


def main():
    root = sys.argv[1]
    do_controls = "--controls" in sys.argv
    rep = sc.Report("ROW pin370: [Pin2026] (PR #370) @a50021b, mu_3 >= 1.5678643*")
    h = sc.sha256(os.path.join(root, CERT))
    rep.info("sha256 " + CERT, h)
    rep.check("sha256 equals the one in the packet's SHA256SUMS", h == STATED_SHA)
    rep.info("the row's fragment 53381b5f...e1e8", "first 8 agree: %s; last 4 agree: %s (the file ends in ...%s)"
             % (h.startswith("53381b5f"), h.endswith("e1e8"), h[-4:]))
    with open(os.path.join(root, CERT), "rb") as fh:
        doc = json.loads(fh.read().decode("ascii"))
    n, ports, fam = doc["n"], doc["ports"], doc["masks"]
    rep.check("one block: 28 points, in-port 0, out-port 1, 252303 members",
              n == 28 and ports == [0, 1] and len(fam) == 252303 and all(0 <= m < 1 << 28 for m in fam))
    rep.info("set sizes", dict(sorted(sc.Counter(map(sc.popcount, fam)).items())))
    kblocks = {"box": [sc.mask(range(8)), sc.mask(range(8, 18)), sc.mask(range(18, 28))],
               "box2": [sc.mask(range(4)), sc.mask(range(4, 8)), sc.mask(range(8, 18)), sc.mask(range(18, 28))]}
    sc.exhaustive(rep, "block", fam, ["box", "box2"], kblocks)
    N = sc.port_counts(fam, ports[0], ports[1])
    T = sc.matrix_checks(rep, N, ((64728, 88867), (33980, 64728)), " (ports 0 -> 1)")
    rep.check("port counts N00 = 33980, N01 = N10 = 64728, N11 = 88867 as the row states",
              (N[0][0], N[0][1], N[1][0], N[1][1]) == (33980, 64728, 64728, 88867))
    # second path (information): the stated difference from Koshta's block, and a scan through the added sets
    with open(os.path.join(root, DELTA), "rb") as fh:
        delta = json.loads(fh.read().decode("ascii"))
    kos = sc.read_sunflower3(os.path.join(root, KOS))[0][3]
    rep.check("Koshta's block: 252292 members (the file replayed in ../c49)", len(kos) == 252292)
    removed, added = delta["removed"], delta["added"]
    rep.check("the certificate is Koshta's block minus the 21 listed sets plus the 32 listed sets",
              len(removed) == 21 and len(added) == 32 and set(removed) <= set(kos) and not set(added) & set(kos)
              and set(fam) == (set(kos) - set(removed)) | set(added))
    rep.check("the 21 removed sets are Koshta's 20 sets inside the first copy of H and his extra set",
              sum(1 for m in removed if m & ~sc.mask(range(8, 18)) == 0) == 20
              and sum(1 for m in removed if m & sc.mask(range(18, 28)) == sc.mask(range(18, 28))) == 1)
    inside = [sum(1 for m in added if m & ~sc.mask(range(a, a + 10)) == 0) for a in (8, 18)]
    rep.check("the 32 added sets: seven points each, 16 inside each pendant copy of H, avoiding both ports",
              all(sc.popcount(m) == 7 and m & 3 == 0 for m in added) and inside == [16, 16], "per pendant %s" % inside)
    index = {m: i for i, m in enumerate(fam)}
    t = time.time()
    c = s = None
    for i in sorted(index[m] for m in added):  # the checker's pivot scan, in this process (no worker processes)
        ci, si, _ = sc.pivot_range(fam, i, i + 1, 1, everyone=True)
        c, s = c or ci, s or si
    sc.verify_witnesses(fam, c, s)
    rep.check("second path: no containment and no 3-sunflower through any of the 32 added sets (each scanned "
              "against all 252303 members)", c is None and s is None, "%.1f s" % (time.time() - t))
    S = sc.Surd(T)
    rep.check("Perron root equals 64728 + sqrt(3019700660), 3019700660 = 88867 * 33980",
              S.t == 2 * 64728 and S.D == 4 * 3019700660 and 3019700660 == 88867 * 33980)
    rep.info("points per copy", "28 - 2 = 26 (chain families on 26k + 2 points)")
    sc.bound_checks(rep, S, 26, "1.5678643", [("1.5678643", "page"), ("1.56786438", "the row text of PR #370"),
                                              ("1.567864384572806979", "PR #370 description"),
                                              ("1.5678643845728069793620", "the certificate's file name")])
    if do_controls:
        sc.controls(rep, "block", fam, ["box", "box2"], kblocks)
        sc.perturbation_controls(rep, N, T, S, 26, "1.5678643")
        fam2 = [m for m in fam if m != 64768]
        N2 = sc.port_counts(fam2, 0, 1)
        S2 = sc.Surd(sc.T_of(N2))
        rep.check("control: without the added set 64768 (still a valid family) N00 = 33979 and the printed "
                  "1.5678643 is refuted", N2[0][0] == 33979 and not sc.check_printed(S2, 26, "1.5678643")[0],
                  "rate would be %s" % sc.dec_str(S2.rate_floor(26, 10), 10))
    ok = rep.done()
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
