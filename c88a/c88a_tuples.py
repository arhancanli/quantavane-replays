#!/usr/bin/env python3
"""88a, program A: the admissible tuples behind the rows of the bounded prime gap page.

Each upper bound H_1 <= D on the page is "DHL[k, 2] holds" (a sieve theorem, NOT checked here) applied to an
explicit admissible k-tuple of diameter D. This program checks the tuple only:

  * the entries are distinct integers, there are exactly k of them, and max - min = D;
  * admissibility: for every prime p there is a residue class mod p that no entry occupies. For p > k this is
    automatic (k entries cannot fill p classes), so every prime p <= k is examined, by listing the occupied classes.

For the [Zha14] row the tuple is the k0 = 3,500,000 primes that follow k0. There the check is: a sieve of
Eratosthenes to 7 * 10^7, the entries are the k0 primes after k0, every entry is a prime greater than k0 (so no entry
is divisible by a prime p <= k0, and class 0 mod p is free), the exact diameter, and Zhang's inequality
pi(7 * 10^7) - pi(k0) > k0.

Standard library only; run with python3 -I. Nothing downloaded is executed: the Lean file is read as text, the
JSON file with the json module, the tuple files as lines of integers.

usage:  c88a_tuples.py list  FILE K D [--swap OLD NEW]      one integer per line
        c88a_tuples.py lean  FILE K D [--swap OLD NEW]      def admissibleTuple : Finset N := {...}
        c88a_tuples.py json  FILE K D [--swap OLD NEW]      {"admissible_tuple": {"tuple": [...]}}
        c88a_tuples.py zhang K0 LIMIT [--from-first]        the K0 primes after K0 (or, as a control, the first K0)
"""
import json
import re
import sys


def small_primes(n):
    flags = bytearray([1]) * (n + 1)
    flags[0:2] = b"\x00\x00"
    i = 2
    while i * i <= n:
        if flags[i]:
            flags[i * i :: i] = bytearray(len(range(i * i, n + 1, i)))
        i += 1
    return [i for i in range(n + 1) if flags[i]]


def read_tuple(kind, path):
    if kind == "list":
        with open(path, encoding="ascii") as fh:
            return [int(line) for line in fh.read().split()]
    if kind == "lean":
        with open(path, encoding="utf-8") as fh:
            text = fh.read()
        hits = re.findall(r"def admissibleTuple : Finset \S+ :=\s*\{([0-9,\s]*)\}", text)
        if len(hits) != 1:
            raise SystemExit("A: expected exactly one admissibleTuple literal, found %d" % len(hits))
        return [int(tok) for tok in hits[0].replace(",", " ").split()]
    if kind == "json":
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
        values = data["admissible_tuple"]["tuple"]
        if not all(type(v) is int for v in values):
            raise SystemExit("A: a tuple entry is not an integer")
        return list(values)
    raise SystemExit("A: unknown kind " + kind)


def check_tuple(entries, k, diam):
    """Return (ok, lines). Every prime p <= k is examined by its set of occupied residue classes."""
    lines = []
    ok = True
    ordered = sorted(entries)
    distinct = len(set(entries)) == len(entries)
    lines.append("entries %d (claimed k = %d), distinct: %s" % (len(entries), k, distinct))
    if not distinct or len(entries) != k:
        ok = False
    width = ordered[-1] - ordered[0]
    lines.append("min %d, max %d, diameter %d (claimed %d)" % (ordered[0], ordered[-1], width, diam))
    if width != diam:
        ok = False
    tightest = None
    blocked = []
    primes = small_primes(max(k, 2))
    for p in primes:
        occupied = set(h % p for h in ordered)
        free = p - len(occupied)
        if free == 0:
            blocked.append(p)
        if tightest is None or free < tightest[1]:
            tightest = (p, free)
    lines.append("primes p <= k examined: %d (largest %d)" % (len(primes), primes[-1]))
    lines.append("fewest free classes: %d, at p = %d" % (tightest[1], tightest[0]))
    if blocked:
        ok = False
        lines.append("every class occupied modulo: %s" % blocked)
    lines.append("admissible: %s" % (not blocked))
    return ok, lines


def zhang(k0, limit, from_first):
    lines = []
    flags = bytearray([1]) * (limit + 1)
    flags[0:2] = b"\x00\x00"
    i = 2
    while i * i <= limit:
        if flags[i]:
            flags[i * i :: i] = bytearray(len(range(i * i, limit + 1, i)))
        i += 1
    pi_k0 = sum(flags[: k0 + 1])
    pi_limit = sum(flags)
    lines.append("pi(%d) = %d, pi(%d) = %d, difference %d" % (k0, pi_k0, limit, pi_limit, pi_limit - pi_k0))
    start = 0 if from_first else k0 + 1
    entries = []
    n = flags.find(1, start)
    while len(entries) < k0 and n >= 0:
        entries.append(n)
        n = flags.find(1, n + 1)
    ok = len(entries) == k0
    lines.append("entries %d (k0 = %d): the primes from %d to %d" % (len(entries), k0, entries[0], entries[-1]))
    width = entries[-1] - entries[0]
    lines.append("diameter %d; below %d: %s" % (width, limit, width < limit))
    if not from_first:
        ineq = pi_limit - pi_k0 > k0
        lines.append("Zhang's inequality pi(%d) - pi(%d) > %d: %s" % (limit, k0, k0, ineq))
        ok = ok and ineq
    ok = ok and width < limit
    above = entries[0] > k0
    lines.append("every entry is a prime greater than k0: %s" % above)
    # direct look at the first primes, by occupied classes (also what refuses the control)
    blocked = []
    for p in small_primes(47):
        if len({h % p for h in entries}) == p:
            blocked.append(p)
    lines.append("primes p <= 47, every class occupied modulo: %s" % (blocked if blocked else "none"))
    admissible = above and not blocked
    lines.append("admissible: %s" % admissible)
    return ok and admissible, lines


def main(argv):
    if len(argv) >= 4 and argv[1] == "zhang":
        ok, lines = zhang(int(argv[2]), int(argv[3]), "--from-first" in argv)
    elif len(argv) >= 5 and argv[1] in ("list", "lean", "json"):
        entries = read_tuple(argv[1], argv[2])
        if "--swap" in argv:
            j = argv.index("--swap")
            old, new = int(argv[j + 1]), int(argv[j + 2])
            if old not in entries or new in entries:
                raise SystemExit("A: bad --swap")
            entries[entries.index(old)] = new
        ok, lines = check_tuple(entries, int(argv[3]), int(argv[4]))
    else:
        raise SystemExit(__doc__)
    for line in lines:
        print("A: " + line)
    print("A: VERDICT %s" % ("PASS" if ok else "REFUSED"))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
