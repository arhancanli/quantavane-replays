"""Take the seven decimal strings of q from code/capaxis.py (YouhuaLi/tight-alternating-knot-constant @ debf7c92) by
the route our coordinator allowed on 12 Nov 2026: Python's ast parser, the single dict literal assigned to QC, and the
one string Constant argument of each Fr(...) call in it. Nothing is executed and no other node of the file is printed
or used. Usage: python3 -I take_q.py CAPAXIS_PY [OUT_JSON]"""
import ast, hashlib, json, sys
path = sys.argv[1]
raw = open(path, "rb").read()
print("sha256", hashlib.sha256(raw).hexdigest())
tree = ast.parse(raw)
found = []
for node in ast.walk(tree):
    if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name) and node.targets[0].id == "QC":
        found.append(node.value)
    if isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name) and node.target.id == "QC":
        found.append(node.value)
assert len(found) == 1, f"{len(found)} assignments to QC"
d = found[0]
assert isinstance(d, ast.Dict), type(d).__name__
out = {}
for k, v in zip(d.keys, d.values):
    assert isinstance(k, ast.Constant) and type(k.value) is int, "key not an int constant"
    assert isinstance(v, ast.Call) and isinstance(v.func, ast.Name) and v.func.id == "Fr", "value not a call of Fr"
    assert len(v.args) == 1 and not v.keywords, "not one positional argument"
    a = v.args[0]
    assert isinstance(a, ast.Constant) and type(a.value) is str, "argument not a string constant"
    out[k.value] = a.value
print(json.dumps(out, indent=1))
if len(sys.argv) > 2:
    json.dump(out, open(sys.argv[2], "w"))
