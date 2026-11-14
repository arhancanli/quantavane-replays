#!/usr/bin/env python3
"""Program 1 for PR #280's update of 14 November 2026 (v6): our existing C_82 checker, ../c82/c82_check.py, unchanged,
run on the v6 certificate. This file only registers the row (title, pinned SHA-256, printed bound, stated D and the
stated upper bounds of the pull request comment) and then calls that checker's own run_row, exactly as
`c82_check.py ROW CERTS_ROOT --controls --workers=N` does for the earlier rows.

  c82v6_first.py CERTS_ROOT [--controls] [--workers=N]       Python 3.9+, standard library only; run with python3 -I.
"""
import importlib.util
import os
import sys

here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("c82_check", os.path.join(here, "..", "c82", "c82_check.py"))
c82 = importlib.util.module_from_spec(spec)
sys.modules["c82_check"] = c82          # so that the worker processes find the module's functions by name
spec.loader.exec_module(c82)

c82.ROWS["pr280f"] = {
    "title": "PR #280 [Pat26] updated 14 Nov: C_82 < 0.2519694883* (open; tao/certificates/c82/v6 at e3bae2e)",
    "file": "pr280f/certificate.json",
    "sha": "29dcbf21a10ddb600a9866fe8ee7e7b43e4357926d74f374a8b6fbf41ada71b5",
    "sha_src": "the PR comment of 14 Nov", "printed": "0.2519694883", "D": 1000000000064,
    # the two values the PR comment states for B are compared and reported ("noted"), not conditions of the bound
    "stated": [("PR comment: Arb checker, tolerance 1e-11, upper bound B <=", "0.25196948827272", "noted"),
               ("PR comment: second run, tolerance 1e-12 and 256 bits, B <=", "0.25196948827270", "noted")],
    "encl": None,
}


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if len(args) != 1:
        print(__doc__)
        sys.exit(2)
    workers = 1
    for a in sys.argv:
        if a.startswith("--workers="):
            workers = int(a.split("=")[1])
    ok = c82.run_row("pr280f", args[0], "--controls" in sys.argv, workers)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
