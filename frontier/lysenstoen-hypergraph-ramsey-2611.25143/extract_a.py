"""For program A only: copy the 0/1 strings out of the ancillary JSON files into plain text files (data, no logic)."""
import json, sys, os
lit, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
pairs = [("R453_37_points.json", "colouring", "c37.txt"), ("R463_70_points.json", "colouring", "c70.txt"),
         ("R553_98_points.json", "colouring", "c98.txt"), ("R665_83_points_short.json", "short_form", "s83.txt"),
         ("R564_200_points_rule.json", "table", "t199.txt"), ("R664_600_points_rule.json", "table", "t599.txt")]
for f, key, o in pairs:
    with open(os.path.join(lit, f)) as fh:
        s = json.load(fh)[key]
    with open(os.path.join(out, o), "w") as fh:
        fh.write(s)
    print(o, len(s))
