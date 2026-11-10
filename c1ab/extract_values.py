"""Copy the step heights of the 1a and 1b rows replayed here out of their published files, as text, into
{"values": [...]} JSON files for ../autocorr_check.py.  Nothing downloaded is ever run: Python sources and notebook
cells are read with the ast parser, the LaTeX with a regular expression, the leaderboard JSON with json.
usage: python3 -I extract_values.py MODE SOURCE ARGS... OUT
  py    FILE.py VAR [ae_reflect|haugland_reflect]   a list literal (or np.array([...])) assigned to VAR
  cell  FILE.ipynb CELL VAR [ae_reflect]           the same inside one notebook code cell
  tex51 FILE.tex                                   Haugland's 51-step function (arXiv 1609.08000, section 2)
  arena FILE.json ID                               the 'values' of solution ID in an EinsteinArena response
  list  FILE.json                                  a JSON list of numbers (SimpleTES first autocorrelation)
  tuple0 FILE.json                                 the first array of a SimpleTES tuple ({'items': [{'data': [...]}...]})
  named FILE.json NAME                             the 'list' of the entry named NAME (ThetaEvolve's data.json)
"""
import ast
import json
import re
import sys
from fractions import Fraction


def literal_list(src, var):
    node = next(n for n in ast.walk(ast.parse(src)) if isinstance(n, ast.Assign) and getattr(n.targets[0], "id", "") == var)
    v = node.value.args[0] if isinstance(node.value, ast.Call) else node.value
    assert isinstance(v, ast.List)
    toks = [ast.get_source_segment(src, e).strip() for e in v.elts]
    for t in toks:
        float(t)  # plain number literals only
    return toks


def reflect(toks, how):
    if how == "ae_reflect":          # the half, then the half reversed without its last entry repeated
        return toks[:-1] + toks[::-1]
    if how == "haugland_reflect":    # f(x) = f(2 - x): the half, then the half reversed without the middle step
        return toks + toks[-2::-1]
    return toks


def tex51(path):
    """The paper's 51-step function: f = v on [a/51, b/51) for x <= 1, and f(x) = f(2 - x) for x > 1."""
    text = open(path, encoding="utf-8").read()
    start = text.index("51 steps")
    rows = re.findall(r"f\(x\)=([0-9.]+) \\text\{ for \} (?:0|\{(\d+)\\over 51\})\\leq x (<|\\leq) "
                      r"(?:\{(\d+)\\over 51\}|1)", text[start:])
    cells = []  # heights of the 2/51-steps on [0, 1 + 1/51), half-step at the end
    pos = Fraction(0)
    for val, a, rel, b in rows:
        lo = Fraction(int(a), 51) if a else Fraction(0)
        hi = Fraction(int(b), 51) if b else Fraction(1)
        assert lo == pos, "gap in the step list"
        width = hi - lo
        k = width / Fraction(2, 51)
        if hi == 1:
            assert k == Fraction(1, 2), "the last piece is the left half of the middle step"
            cells.append(val)
        else:
            assert k.denominator == 1
            cells += [val] * int(k)
        pos = hi
        if hi == 1:
            break
    assert pos == 1 and len(cells) == 26, (pos, len(cells))
    return cells + cells[-2::-1]


def main():
    mode, args = sys.argv[1], sys.argv[2:]
    out = args[-1]
    if mode == "py":
        toks = reflect(literal_list(open(args[0]).read(), args[1]), args[2] if len(args) > 3 else "")
    elif mode == "cell":
        src = "".join(json.load(open(args[0]))["cells"][int(args[1])]["source"])
        toks = reflect(literal_list(src, args[2]), args[3] if len(args) > 4 else "")
    elif mode == "tex51":
        toks = tex51(args[0])
    elif mode == "arena":
        d = json.load(open(args[0]), parse_float=str)
        sol = [s for s in d if s["id"] == int(args[1])]
        assert len(sol) == 1, "solution not in the response"
        toks = [str(t) for t in sol[0]["data"]["values"]]
    elif mode == "list":
        toks = [str(t) for t in json.load(open(args[0]), parse_float=str)]
    elif mode == "tuple0":
        x = json.load(open(args[0]), parse_float=str)["items"][0]
        assert x["dtype"] == "float64" and len(x["shape"]) == 1 and len(x["data"]) == x["shape"][0]
        toks = [str(t) for t in x["data"]]
    elif mode == "named":
        ent = [e for e in json.load(open(args[0]), parse_float=str) if e["name"] == args[1]]
        assert len(ent) == 1, "entry not found"
        toks = [str(t) for t in ent[0]["list"]]
    else:
        raise SystemExit(__doc__)
    for t in toks:
        assert re.fullmatch(r"-?[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?", t), "not a plain number: " + t
    open(out, "w").write('{"values": [' + ", ".join(toks) + "]}\n")
    print(f"{out}: {len(toks)} heights")


if __name__ == "__main__":
    main()
