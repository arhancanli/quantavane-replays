#!/usr/bin/env python3
"""Program 1 for the 82a row [Gri26] C_82 <= 0.2536331090204145 (maaxgrin/zhang-zagier-c82-bound @ cb88f92): our existing C_82 checker, ../c82/c82_check.py, unchanged,
run on the [Gri26] certificate (denominators P7, P9, Q1, Q2, R0, R2 of Doche's polynomials and their perturbations,
numerators P1..P6, P8 with rational exponents). This file only registers the row (title, pinned SHA-256, printed bound, stated D and the
values its README and run record state) and then calls that checker's own run_row, exactly as
`c82_check.py ROW CERTS_ROOT --controls --workers=N` does for the earlier rows.

  c82gri_first.py CERTS_ROOT [--controls] [--workers=N]       Python 3.9+, standard library only; run with python3 -I.
"""
import importlib.util
import os
import sys

here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("c82_check", os.path.join(here, "..", "c82", "c82_check.py"))
c82 = importlib.util.module_from_spec(spec)
sys.modules["c82_check"] = c82          # so that the worker processes find the module's functions by name
spec.loader.exec_module(c82)

c82.ROWS["gri26"] = {
    "title": "82a [Gri26] C_82 <= 0.2536331090204145 (merged row; maaxgrin/zhang-zagier-c82-bound @ cb88f92, the record"
             " qstar_lp_grid32768_maxden1e7 of runs/qstar_targeted_refined.jsonl, the one its certified run names)",
    "file": "gri26/certificate.json",
    "sha": "ccbda5c371802b765e7a01f5a9a5a7b7b84c907719a28ff79a7ef44ded0f3aee",
    "sha_src": "our extraction (fetch.sh) of line 1 of runs/qstar_targeted_refined.jsonl, sha256 698a2ed5...b1a6",
    "printed": "0.2536331090204145", "D": None,
    "stated": [("README: the numerical value for this Q", "0.253626920959", "noted"),
               ("runs/qstar_cert_roots_262144_target253635.jsonl: bound_upper", "0.2536331090204145", "noted")],
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
    ok = c82.run_row("gri26", args[0], "--controls" in sys.argv, workers)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
