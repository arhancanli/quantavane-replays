#!/usr/bin/env python3
"""Copy the step heights used by the c1mv replays out of their published files, as text, into {"values": [...]} JSON
files for ../autocorr_check.py.  Nothing downloaded is run: the LaTeX is read with a regular expression (the
verbatim blocks of the appendix, located by the sentence before each), the notebook cell with Python's ast module.

  mv2009_n208.json     Matolcsi and Vinuesa, arXiv 0907.1379v2, appendix: the best non-negative step function
                       (n = 208; "||f*f||_inf = 1.50972...")          -> the 1a row [MV2009]
  mv2009_n150s.json    the same appendix: the signed example (n = 150; "||f*f||_inf = 1.45810...") -> the note
  ae_c3prime.json      alphaevolve_results @ 4226acb, mathematical_results.ipynb, section B.3, cell 65
                       (height_sequence_4: the signed variant C_3' <= 1.4688)  -> the 1d row [GGSWT2025]
Each number is kept exactly as printed (a JSON number literal; a bare "0." becomes "0.0"), so autocorr_check.py
can read it either as binary64 or, with --decimal, as the printed decimal.
usage: python3 -I extract_mv.py CERTS_ROOT
"""
import ast
import json
import re
import sys

root = sys.argv[1]
NUM = re.compile(r"^-?\d+\.\d*(?:e[-+]?\d+)?$|^-?\d+$")


def write(name, toks, source):
    toks = [t + "0" if t.endswith(".") else t for t in toks]   # "0." (printed in the n = 208 list) is 0.0, exactly
    for t in toks:
        assert NUM.match(t), f"{name}: not a plain number: {t!r}"
        json.loads(t)
    with open(f"{root}/{name}", "w") as f:
        f.write('{"source": ' + json.dumps(source) + ', "values": [' + ", ".join(toks) + "]}\n")
    print(f"{name}: {len(toks)} numbers from {source}")


tex = open(f"{root}/mv2009/sch0903b.tex", encoding="utf-8").read()
blocks = [(m.start(), m.group(1).split()) for m in re.finditer(r"\\begin\{verbatim\}(.*?)\\end\{verbatim\}", tex, re.S)]


def block_after(phrase):
    i = tex.index(phrase)
    pos, toks = next((p, t) for p, t in blocks if p > i)
    return toks


n208 = block_after("is\nattained at $n = 208$")
assert len(n208) == 208, len(n208)
write("mv2009_n208.json", n208, "arXiv 0907.1379v2, appendix, n = 208")
n150 = block_after("is attained at $n=150$")
assert len(n150) == 150, len(n150)
write("mv2009_n150s.json", n150, "arXiv 0907.1379v2, appendix, n = 150 (signed)")

nb = json.load(open(f"{root}/ae/dm_mathematical_results.ipynb", encoding="utf-8"))
src = "".join(nb["cells"][65]["source"])
src = "\n".join(l for l in src.split("\n") if not l.lstrip().startswith(("%", "!")))
node = next(n for n in ast.walk(ast.parse(src)) if isinstance(n, ast.Assign) and getattr(n.targets[0], "id", "") == "height_sequence_4")
v = node.value.args[0] if isinstance(node.value, ast.Call) else node.value
assert isinstance(v, ast.List)
toks = [ast.get_source_segment(src, e).strip() for e in v.elts]
write("ae_c3prime.json", toks, "alphaevolve_results @ 4226acb, mathematical_results.ipynb, cell 65, height_sequence_4")
