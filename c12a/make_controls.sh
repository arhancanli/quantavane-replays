#!/bin/sh
# Write the red controls into controls/: altered copies of the fetched spec and values files (no other input).
#   (a) k5_sum_plus_one:   the k5 values with the cell sum S_100 increased by 1
#   (b) k7_threshold_below: the k7 spec with the threshold one unit of 10^-12 below the bound printed in results/k7.txt
#       (rounded up there, so below B); k7_threshold_at: the threshold equal to that printed bound (at least B)
#   (c) k5_mesh127:        the k5 spec and a copy of its values file with mesh 127 in place of 128, all else as stored;
#       k5_mesh127_orig:   the same spec with the unaltered k5 values file
# Usage: ./make_controls.sh [CERTS_ROOT]   (default: ./certs; results/k7.txt must exist)
set -e
here=$(cd "$(dirname "$0")" && pwd)
root=${1:-$here/certs}
cd "$here"
mkdir -p controls
rel=$(python3 -I -c 'import os, sys; print(os.path.relpath(sys.argv[1], sys.argv[2]))' "$root/tsp/data" controls)
spec=$root/tsp/data/specs
vals=$root/tsp/data/values
# (a)
sed "s#^values .*#values k5_sum_plus_one_values.txt#" "$spec/k5.txt" > controls/k5_sum_plus_one.txt
awk -v want=100 '
  $1 == "sums" { s = NR }
  s && NR == s + 1 + want {
    x = $1; i = length(x)
    while (i > 0 && substr(x, i, 1) == "9") { x = substr(x, 1, i - 1) "0" substr(x, i + 1); i-- }
    if (i == 0) x = "1" x; else x = substr(x, 1, i - 1) (substr(x, i, 1) + 1) substr(x, i + 1)
    print x; next }
  { print }' "$vals/k5.txt" > controls/k5_sum_plus_one_values.txt
# (b)
b=$(sed -n 's/^bound \([0-9]*\)\.\([0-9]*\) (rounded up)$/\1\2/p' results/k7.txt)
[ -n "$b" ] || { echo "no bound line in results/k7.txt"; exit 2; }
b=$(echo "$b" | sed 's/^0*//')
below=$(echo "$b - 1" | bc)
sed -e "s#^threshold .*#threshold $below 1000000000000#" -e "s#^values .*#values $rel/values/k7.txt#" \
    "$spec/k7.txt" > controls/k7_threshold_below.txt
sed -e "s#^threshold .*#threshold $b 1000000000000#" -e "s#^values .*#values $rel/values/k7.txt#" \
    "$spec/k7.txt" > controls/k7_threshold_at.txt
# (c)
sed -e "s#^mesh 128\$#mesh 127#" -e "s#^values .*#values k5_mesh127_values.txt#" "$spec/k5.txt" > controls/k5_mesh127.txt
sed "s#^mesh 128\$#mesh 127#" "$vals/k5.txt" > controls/k5_mesh127_values.txt
sed -e "s#^mesh 128\$#mesh 127#" -e "s#^values .*#values $rel/values/k5.txt#" "$spec/k5.txt" > controls/k5_mesh127_orig.txt
