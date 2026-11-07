#!/usr/bin/env python3
"""Exact replay of the 3b upper bound C_3b <= 209357544569848128/114631769327023765 = 1.82634836...* [F2026] (the
Katz-Tao sum-difference constant, constants/3b.md; the row of PR #329), and of the paper's all-groups bound
10304335425987192362/5641833174538280329 = 1.82641618..., from the certificate data and the paper (Farnamfar, Zenodo
10.5281/zenodo.23714415, main.tex) alone. The package's verify_group.py, stackcheck.py, ring_lattice.py,
check630*.py and negative_controls.py were not read or run (they are not even extracted: fetch.sh unpacks only the
JSON certificates).

What a certificate proves (paper, Sections 2-5). Start with a difference-injective pair Z = (X, Y) in an abelian group
and build vertices z_v (pairs with the law of Z) by conditional copies: step s copies the source vertices U over
readings V (linear forms in old coordinates), independently of everything old given V (Lemma 4). Every vertex has
coordinates (x_v, y_v); a form is an integer combination of coordinates; the labelled forms are x_v, y_v, s_v = x_v +
y_v, d_v = x_v - y_v at every vertex and the certificate's extra ("cross") forms. Relations: at step s, every form
alpha supported on U lying in span(cl V) + Rel gives alpha(U') - alpha(U) = 0. cl(T): the labelled forms in span(T) +
Rel, closed under "d_v known => x_v, y_v known" (difference injectivity). h(T) := H[cl T].
Rows: sub  h(K,i) + h(K,j) - h(K,i,j) - h(K) >= 0 (multiplier >= 0); mono  h(K,i) - h(K) >= 0 (>= 0);
      ci   h(V,U') + h(O) - h(O,U') - h(V) = 0, O the entire old prefix (any sign);
      lemma  a verified information inequality, each abstract variable replaced by a tuple of forms (>= 0);
      equal laws  h(T) = h(T') for gens T on the source of a step and their images T' (a gen determined by the glue
      V is kept; a gen supported on U is moved to the copy literally, and its image must be a labelled form).
The certificate claims sum_r m_r R_r = a h(x0) + b h(y0) + c h(s0) - h(x0, y0) with every entropy class's
coefficient zero after the equal-law identifications. Then H(X,Y) <= a H(X) + b H(Y) + c H(X+Y) for difference-
injective pairs; Lemma 6 (a random selector) removes injectivity and Corollary 12 gives |A - B| <= M^(a+b+c).

Our checker, written from that description:
  - the 72 (record) or 77 (all groups) library lemmas, each from its own conditional-copy proof in abstract variables
    (sub, eqlaw, ci, apply rows; every earlier lemma verified before it is applied): the rows must sum exactly to the
    lemma's coefficient vector, all in exact rationals;
  - the configuration: vertices, the relation module step by step (prefix closures over Q, intersections with the
    source-supported forms), compared with the paper's counts (rank 71, quotient dimension 13 for both records);
  - every queried tuple named by its rational closure; span membership in the 13-dimensional quotient decided modulo
    p = 2^127 - 1 on integer image vectors, which is exact because p exceeds Hadamard's bound for every minor (checked);
  - every row, equal law (with its fixed/moved gens re-derived at its prefix) and the target; every class must cancel.
Rational closures are sound in torsion-free groups (clear denominators: n f = known and n is injective), in particular
in R, which is what C_3b concerns. The paper's Z[1/2] (no 2-torsion) and Z (every abelian group) domain audits are not
replayed here; for those wider statements see its Section 4 (Lemma 9, Proposition 10).

  python3 -I c3b_check.py [record|allgroups|both|selftest] [CERTS] [--controls]
(selftest: the paper's two standard-form certificates, 2762/1509 and 152/83, a different configuration)
Python 3.9+, standard library only.
"""
import copy
import hashlib
import json
import math
import os
import sys
import time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
PRIMES = ((1 << 127) - 1, (1 << 521) - 1)        # Mersenne primes; the least one above Hadamard's bound is used
CERTS = {
    "record": ("f2026/record.json", "12e2e356274bad4e044c1558908671cc7c543707b82450f83c22ee7a88046a70"),
    "allgroups": ("f2026/all-groups.json", "798485494f517e1ec7cbf8916cbbf046fdcf233d705dfc92dcf2c0fff07e03f1"),
    "standard": ("f2026/standard.json", "f0cc88a42829a4fd515dcad2bfecd9bc5c586339b9386c11199115f241218fe3"),
    "fivepage": ("f2026/standard-five-page.json", "209f5044bc7c745d3c89e6a3ca0979d48bbacf0fc2dcab1554662ed6b79e9654"),
}
STATED = {
    "record": dict(C=Fraction(209357544569848128, 114631769327023765), rank=71, quot=13, sub=519, mono=0, ci=9,
                   lemma=12, eqlaw=167, lib=72, extra=78, vertices=42),
    "allgroups": dict(C=Fraction(10304335425987192362, 5641833174538280329), rank=71, quot=13, sub=566, mono=1, ci=9,
                      lemma=11, eqlaw=179, lib=77, extra=85, vertices=42),
    # the paper's Theorem 3 (standard vertex forms only, a chain of pages; Shannon rows only): self-tests of our
    # reading of the format on a different configuration
    "standard": dict(C=Fraction(2762, 1509), rank=177, quot=39, sub=758, mono=0, ci=20, lemma=0, eqlaw=568, lib=0,
                     extra=0, vertices=108),
    "fivepage": dict(C=Fraction(152, 83), rank=46, quot=14, sub=103, mono=0, ci=8, lemma=0, eqlaw=59, lib=0,
                     extra=0, vertices=30),
}
SLOPES = {(1, 0): "x", (0, 1): "y", (1, 1): "s", (1, -1): "d"}


def fr(v):
    if isinstance(v, bool):
        raise ValueError("boolean where a number was expected")
    return Fraction(str(v)) if not isinstance(v, int) else Fraction(v)


# ------------------------------------------------------------------------------------------- abstract lemmas
class Lin:
    """A formal linear combination of entropy terms (keys: frozensets), exact rationals."""

    def __init__(self):
        self.c = {}

    def add(self, key, m):
        key = frozenset(key)
        if not key:
            return                                    # h(empty) = 0
        v = self.c.get(key, 0) + m
        if v:
            self.c[key] = v
        else:
            self.c.pop(key, None)


def vset(s):
    return frozenset(s)


def verify_lemma(name, L, verified, out):
    """Check one library lemma from its own proof; returns its coefficient dict (frozenset -> Fraction)."""
    base, allv = L["base"], L["vars"]
    assert len(set(allv)) == len(allv) and allv.startswith(base), f"{name}: variables"
    existing = list(base)
    steps = []
    for st in L["steps"]:
        new, src, over = st["new"], st["src"], st["over"]
        assert len(new) == len(src) >= 1, f"{name}: copy step shape"
        assert all(c not in existing for c in new) and len(set(new)) == len(new), f"{name}: copied variables not fresh"
        assert all(set(s) <= set(existing) and s for s in src), f"{name}: copy source not existing"
        assert set(over) <= set(existing), f"{name}: conditioning set not existing"
        steps.append(dict(new=new, src=list(src), over=over, old=list(existing)))
        existing += list(new)
    assert set(existing) == set(allv), f"{name}: variables do not match the steps"
    target = {}
    for k, v in L["coef"].items():
        assert set(k) <= set(base), f"{name}: target term outside the base variables"
        target[vset(k)] = target.get(vset(k), 0) + fr(v)
    acc = Lin()
    for row in L["proof"]:
        kind = row[0]
        if kind == "sub":
            i, j, K, m = vset(row[1]), vset(row[2]), vset(row[3]), fr(row[4])
            assert m >= 0 and (i | j | K) <= set(existing), f"{name}: bad sub row {row}"
            acc.add(i | K, m), acc.add(j | K, m), acc.add(i | j | K, -m), acc.add(K, -m)
        elif kind == "eqlaw":
            s, T, m = row[1], row[2], fr(row[3])
            st = steps[s]
            allowed = set(st["new"]) | set(st["over"])
            assert set(T) <= allowed, f"{name}: equal law outside new + conditioning variables {row}"
            img = set(c for c in T if c not in st["new"])
            for c in T:
                if c in st["new"]:
                    img |= set(st["src"][st["new"].index(c)])
            acc.add(vset(T), m), acc.add(vset(img), -m)
        elif kind == "ci":
            s, m = row[1], fr(row[2])
            st = steps[s]
            O, V, Un = vset(st["old"]), vset(st["over"]), vset(st["new"])
            acc.add(V | Un, m), acc.add(O, m), acc.add(O | Un, -m), acc.add(V, -m)
        elif kind == "apply":
            lname, args, m = row[1], row[2], fr(row[3])
            assert lname in verified, f"{name}: applies {lname} before it is verified"
            assert m >= 0, f"{name}: negative multiplier on an applied lemma"
            lb = verified[lname]["base"]
            assert len(args) == len(lb) and all(set(a) <= set(existing) for a in args), f"{name}: bad apply {row}"
            mp = dict(zip(lb, args))
            for S, cf in verified[lname]["coef"].items():
                acc.add(frozenset(c for v in S for c in mp[v]), m * cf)
        else:
            raise AssertionError(f"{name}: unknown proof row {kind}")
    diff = dict(acc.c)
    for k, v in target.items():
        diff[k] = diff.get(k, 0) - v
    bad = {k: v for k, v in diff.items() if v}
    if bad:
        raise AssertionError(f"{name}: the proof rows do not sum to the lemma ({len(bad)} terms differ)")
    return dict(base=base, coef=target)


def verify_library(lemmas, out):
    verified = {}
    pending = dict(lemmas)
    while pending:                                    # in dependency order (each applied lemma verified first)
        progress = False
        for name in list(pending):
            deps = {r[1] for r in pending[name]["proof"] if r[0] == "apply"}
            if deps <= set(verified):
                verified[name] = verify_lemma(name, pending.pop(name), verified, out)
                progress = True
        if not progress:
            raise AssertionError(f"lemma dependencies are circular or missing: {sorted(pending)}")
    return verified


# ----------------------------------------------------------------------------------------- rational linear algebra
class RSpan:
    """A subspace of Q^n by reduced row echelon rows (dicts col -> Fraction)."""

    def __init__(self, n):
        self.n, self.rows = n, {}                    # pivot col -> row

    def reduce(self, v):
        v = {k: Fraction(x) for k, x in v.items()}            # exact (never int / int -> float)
        for c in sorted(v):
            if c in v and c in self.rows and v[c]:
                f = v[c]
                for k, x in self.rows[c].items():
                    nv = v.get(k, 0) - f * x
                    if nv:
                        v[k] = nv
                    else:
                        v.pop(k, None)
        return {k: x for k, x in v.items() if x}

    def add(self, v):
        r = self.reduce(v)
        if not r:
            return False
        p = min(r)
        piv = r[p]
        r = {k: x / piv for k, x in r.items()}
        for q, row in self.rows.items():
            if p in row:
                f = row[p]
                for k, x in r.items():
                    nv = row.get(k, 0) - f * x
                    if nv:
                        row[k] = nv
                    else:
                        row.pop(k, None)
        self.rows[p] = r
        return True

    def contains(self, v):
        return not self.reduce(v)

    def copy(self):
        s = RSpan(self.n)
        s.rows = {k: dict(r) for k, r in self.rows.items()}
        return s


# ---------------------------------------------------------------------------------------------- the configuration
class Config:
    def __init__(self, d, out):
        self.d = d
        self.steps = []
        nv = 1
        for st in d["steps"]:
            U = list(st["U"])
            assert all(0 <= u < nv for u in U) and len(set(U)) == len(U)
            new = list(range(nv, nv + len(U)))
            self.steps.append(dict(U=U, V=[self.fkey(g) for g in st["V"]], new=new, nprev=nv,
                                   phi=dict(zip(U, new))))
            nv += len(U)
        self.nv = nv
        self.n = 2 * nv
        self.vec = {}                                 # labelled form -> dict coord -> int
        for v in range(nv):
            for sl in SLOPES:
                self.vec[("v", v, sl)] = {k: x for k, x in ((2 * v, sl[0]), (2 * v + 1, sl[1])) if x}
        self.extra_top = {}
        for k, e in enumerate(d.get("extra", [])):
            e = [int(x) for x in e]
            assert len(e) % 2 == 0 and len(e) <= self.n, f"extra form {k}: length"
            vecd = {i: x for i, x in enumerate(e) if x}
            assert vecd, f"extra form {k} is zero"
            self.vec[("e", k)] = vecd
            self.extra_top[("e", k)] = max(vecd) // 2      # highest vertex it involves
        self.labels = list(self.vec)
        self.rel = RSpan(self.n)
        self.build_relations(out)

    @staticmethod
    def fkey(g):
        if g[0] == "amb":
            return ("e", int(g[1]))
        sl = tuple(g[1])
        assert sl in SLOPES, f"unexpected slope {sl}"
        return ("v", int(g[0]), sl)

    def top(self, f):
        return f[1] if f[0] == "v" else self.extra_top[f]

    def closure_prefix(self, T, rel, nvert):
        """cl_Q(T) with the relations rel, among labelled forms on vertices < nvert."""
        sp = rel.copy()
        for f in T:
            sp.add(self.vec[f])
        cand = [f for f in self.labels if self.top(f) < nvert]
        members = set()
        while True:
            members = {f for f in cand if sp.contains(self.vec[f])}
            grew = False
            for v in range(nvert):
                if ("v", v, (1, -1)) in members:
                    for sl in ((1, 0), (0, 1)):
                        if ("v", v, sl) not in members:
                            sp.add(self.vec[("v", v, sl)])
                            grew = True
            if not grew:
                return members, sp

    def build_relations(self, out):
        for s, st in enumerate(self.steps):
            nprev = st["nprev"]
            clV, sp = self.closure_prefix(st["V"], self.rel, nprev)
            st["clV"] = clV
            # W = span(cl V) + Rel; intersect with forms supported on the source coordinates U
            ucols = sorted(c for u in st["U"] for c in (2 * u, 2 * u + 1))
            other = [c for c in range(2 * nprev) if c not in set(ucols)]
            order = other + ucols                      # eliminate the non-source coordinates first
            pos = {c: i for i, c in enumerate(order)}
            W = RSpan(self.n)
            gens = [self.vec[f] for f in clV] + [dict(r) for r in self.rel.rows.values()]
            for g in gens:
                W.add({pos[c]: x for c, x in g.items()})
            nother = len(other)
            added = 0
            for p_, row in list(W.rows.items()):
                if p_ >= nother:                       # zero on every non-source coordinate
                    alpha = {order[k]: x for k, x in row.items()}
                    rel = dict(alpha)
                    for c, x in alpha.items():         # alpha(U') - alpha(U)
                        u, xy = divmod(c, 2)
                        cc = 2 * st["phi"][u] + xy
                        rel[cc] = rel.get(cc, 0) - x
                    rel = {k: -x for k, x in rel.items() if x}
                    if self.rel.add(rel):
                        added += 1
            st["added"] = added
        self.rank = len(self.rel.rows)

    # ---- the final quotient, integer images and exact-rank membership modulo p
    def quotient(self):
        piv = set(self.rel.rows)
        self.free = [c for c in range(self.n) if c not in piv]
        idx = {c: i for i, c in enumerate(self.free)}
        self.img = {}
        B = 0
        for f in self.labels:
            r = self.rel.reduce(self.vec[f])
            den = 1
            for x in r.values():
                den = den * x.denominator // math.gcd(den, x.denominator)
            iv = [0] * len(self.free)
            for c, x in r.items():
                iv[idx[c]] = int(x * den)
            g = 0
            for x in iv:
                g = math.gcd(g, x)
            if g > 1:
                iv = [x // g for x in iv]
            self.img[f] = iv
            B = max([B] + [abs(x) for x in iv])
        k = len(self.free)
        hadamard = (math.isqrt(k) + 1) ** k * B ** k     # >= (sqrt(k) B)^k, every minor's bound
        self.hadamard, self.B = hadamard, B
        self.P = next((p for p in PRIMES if p > hadamard), PRIMES[-1])
        return k, B, hadamard

    def cl_name(self, T):
        """The rational closure of the tuple T (labelled forms), as a frozenset; exact via ranks modulo p."""
        key = frozenset(T)
        hit = self.cache.get(key)
        if hit is not None:
            return hit
        k = len(self.free)
        P = self.P
        basis = {}                                     # pivot -> row (mod p), reduced echelon

        def red(v):
            v = list(v)
            for pc, row in basis.items():
                if v[pc]:
                    f = v[pc]
                    v = [(a - f * b) % P for a, b in zip(v, row)]
            return v

        def add(v):
            r = red([x % P for x in v])
            nz = [i for i, x in enumerate(r) if x]
            if not nz:
                return
            pc = nz[0]
            inv = pow(r[pc], P - 2, P)
            r = [(x * inv) % P for x in r]
            for q, row in list(basis.items()):
                if row[pc]:
                    f = row[pc]
                    basis[q] = [(a - f * b) % P for a, b in zip(row, r)]
            basis[pc] = r

        for f in key:
            add(self.img[f])
        while True:
            members = frozenset(f for f in self.labels if not any(red([x % P for x in self.img[f]])))
            grew = False
            for f in members:
                if f[0] == "v" and f[2] == (1, -1):
                    for sl in ((1, 0), (0, 1)):
                        g = ("v", f[1], sl)
                        if g not in members:
                            add(self.img[g])
                            grew = True
            if not grew:
                break
        self.cache[key] = members
        return members


def form_list(lst, cfg):
    return [cfg.fkey(g) for g in lst]


def check_certificate(which, certs, out=print, mutate=None):
    path = os.path.join(certs, CERTS[which][0])
    raw = open(path, "rb").read()
    digest = hashlib.sha256(raw).hexdigest()
    d = json.loads(raw)
    if mutate:
        d = copy.deepcopy(d)
        mutate(d)
    st = STATED[which]
    t0 = time.time()
    out(f"Certificate {which}: {CERTS[which][0]}, sha256 {digest} ({'as stated' if digest == CERTS[which][1] else 'DIFFERS from ' + CERTS[which][1]})")
    ok = digest == CERTS[which][1] or mutate is not None
    # 1. the lemma library
    try:
        lib = verify_library(d.get("lemmas", {}), out)
    except AssertionError as e:
        out(f"  lemma library: REFUSED ({e})")
        return False
    out(f"  lemma library: all {len(lib)} lemmas verified from their own conditional-copy proofs "
        f"({sum(len(L['proof']) for L in d.get('lemmas', {}).values())} proof rows) [{time.time() - t0:.1f} s]")
    ok &= len(lib) == st["lib"]
    # 2. the configuration and relations
    cfg = Config(d, out)
    k, B, had = cfg.quotient()
    cfg.cache = {}
    out(f"  configuration: {cfg.nv} vertices, {len(cfg.steps)} copy steps, {len(d.get('extra', []))} extra forms; relation "
        f"rank {cfg.rank} (per step {[s['added'] for s in cfg.steps]}), quotient dimension {k}")
    out(f"  integer images: largest entry {B}; Hadamard bound for minors {float(had):.3e} < p = 2^{cfg.P.bit_length()} - 1: "
        f"{had < cfg.P}")
    ok &= cfg.nv == st["vertices"] and cfg.rank == st["rank"] and k == st["quot"] and had < cfg.P
    # 3. rows
    acc = {}

    def term(T, m):
        name = cfg.cl_name(T)
        if name:
            acc[name] = acc.get(name, 0) + m

    counts = {"sub": 0, "mono": 0, "ci": 0, "lemma": 0}
    bad = []
    for r in d["rows"]:
        t = r["type"]
        m = fr(r["mult"])
        counts[t] = counts.get(t, 0) + 1
        if t == "sub":
            if m < 0:
                bad.append("negative sub multiplier")
            K = form_list(r["K"], cfg)
            i, j = cfg.fkey(r["i"]), cfg.fkey(r["j"])
            term(K + [i], m), term(K + [j], m), term(K + [i, j], -m), term(K, -m)
        elif t == "mono":
            if m < 0:
                bad.append("negative mono multiplier")
            K = form_list(r["K"], cfg)
            i = cfg.fkey(r["i"])
            term(K + [i], m), term(K, -m)
        elif t == "ci":
            s = cfg.steps[r["step"]]
            O = [("v", v, sl) for v in range(s["nprev"]) for sl in ((1, 0), (0, 1))]
            Un = [("v", v, sl) for v in s["new"] for sl in ((1, 0), (0, 1))]
            V = s["V"]
            term(V + Un, m), term(O, m), term(O + Un, -m), term(V, -m)
        elif t == "lemma":
            if m < 0:
                bad.append("negative lemma multiplier")
            L = lib.get(r["lemma"])
            if L is None:
                bad.append(f"unverified lemma {r['lemma']}")
                continue
            args = [form_list(a, cfg) for a in r["args"]]
            if len(args) != len(L["base"]):
                bad.append("lemma arity")
                continue
            mp = dict(zip(L["base"], args))
            for S, cf in L["coef"].items():
                term([f for v in S for f in mp[v]], m * cf)
        else:
            bad.append(f"unknown row type {t}")
    # 4. equal laws: union-find on closure names
    parent = {}

    def find(x):
        while parent.get(x, x) != x:
            parent[x] = parent.get(parent[x], parent[x])
            x = parent[x]
        return x

    nfixed = nmoved = 0
    for e in d["eqlaw"]:
        s = cfg.steps[e["step"]]
        src, img = [], []
        for g in e["gens"]:
            f = cfg.fkey(g)
            src.append(f)
            if f in s["clV"]:                          # determined by the glue at that prefix: kept
                img.append(f)
                nfixed += 1
                continue
            vec = cfg.vec[f]
            if not all(c // 2 in s["phi"] for c in vec):
                bad.append(f"equal law at step {e['step']}: {g} neither determined by the glue nor on the source")
                continue
            moved = {2 * s["phi"][c // 2] + c % 2: x for c, x in vec.items()}
            hit = [h for h in cfg.labels if cfg.vec[h] == moved]
            if not hit:
                bad.append(f"equal law at step {e['step']}: the image of {g} is not a labelled form")
                continue
            img.append(hit[0])
            nmoved += 1
        a, b = cfg.cl_name(src), cfg.cl_name(img)
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    # 5. the target
    lam = {tuple(int(x) for x in k.split(",")): fr(v) for k, v in d["lambda"].items()}
    a, b, c = lam[(1, 0)], lam[(0, 1)], lam[(1, 1)]
    C = a + b + c
    tgt = {}
    for T, m in (([("v", 0, (1, 0))], a), ([("v", 0, (0, 1))], b), ([("v", 0, (1, 1))], c),
                 ([("v", 0, (1, 0)), ("v", 0, (0, 1))], Fraction(-1))):
        nm = cfg.cl_name(T)
        tgt[nm] = tgt.get(nm, 0) + m
    classes = {}
    for nm, m in acc.items():
        r = find(nm)
        classes[r] = classes.get(r, 0) + m
    for nm, m in tgt.items():
        r = find(nm)
        classes[r] = classes.get(r, 0) - m
    resid = {r: m for r, m in classes.items() if m}
    out(f"  rows: {counts['sub']} sub, {counts['mono']} mono, {counts['ci']} ci, {counts['lemma']} lemma; "
        f"{len(d['eqlaw'])} equal laws ({nfixed} gens kept as determined by the glue, {nmoved} moved to the copy); "
        f"{len(cfg.cache)} distinct tuples named")
    if bad:
        out(f"  REFUSED: {bad[:5]}{' ...' if len(bad) > 5 else ''}")
        ok = False
    out(f"  identity sum_r m_r R_r = a h(x0) + b h(y0) + c h(s0) - h(x0,y0): "
        f"{'every class cancels exactly' if not resid else f'{len(resid)} classes do NOT cancel'} "
        f"({len(classes)} classes after the equal-law identifications)")
    ok &= not resid and not bad
    same = all(counts[t] == st[t] for t in ("sub", "mono", "ci", "lemma")) and len(d["eqlaw"]) == st["eqlaw"]
    out(f"  counts as the paper's table of certificate sizes states: {same}")
    ok &= same
    out(f"  a = {a}, b = {b}, c = {c}; all >= 0: {min(a, b, c) >= 0}")
    ok &= min(a, b, c) >= 0
    out(f"  C = a + b + c = {C} = {dec(C, 22)}...; as stated: {C == st['C'] and C == fr(d['claimed_C'])}")
    ok &= C == st["C"] and C == fr(d["claimed_C"])
    out(f"  [{time.time() - t0:.1f} s]")
    return ok


def dec(x, n):
    x = Fraction(x)
    q = (x.numerator * 10 ** n) // x.denominator
    s = str(q).rjust(n + 1, "0")
    return s[:-n] + "." + s[-n:]


def controls(which, certs):
    print(f"Red controls for {which} (the same checks on altered copies of the certificate):")
    quiet = lambda *a, **k: None
    ok = True

    def run(desc, mut):
        nonlocal ok
        res = check_certificate(which, certs, out=quiet, mutate=mut)
        print(f"  {desc}: {'refused, as expected' if not res else 'ACCEPTED (control failed)'}")
        ok &= not res

    def m1(d):
        r = next(r for r in d["rows"] if r["type"] == "sub")
        r["mult"] = str(fr(r["mult"]) * Fraction(11, 10))
    run("one sub multiplier raised by 10%", m1)

    def m2(d):
        d["lambda"]["1,1"] = str(fr(d["lambda"]["1,1"]) - Fraction(1, 10 ** 12))
    run("c lowered by 10^-12 (a better claimed bound)", m2)

    def m3(d):
        d["eqlaw"].pop(len(d["eqlaw"]) // 2)
    run("one equal-law entry removed", m3)

    def m4(d):
        e = next(e for e in d["eqlaw"] if any(g[0] != "amb" and e["step"] >= 3 for g in e["gens"]))
        g = next(g for g in e["gens"] if g[0] != "amb")
        g[1] = [1, 1] if g[1] != [1, 1] else [1, 0]
    run("one equal-law gen changed to another reading", m4)

    def m5(d):
        L = d["lemmas"]["dfz109"]
        r = next(r for r in L["proof"] if r[0] == "sub")
        r[4] = str(fr(r[4]) + 1)
    run("one coefficient in the proof of dfz109 changed", m5)

    def m6(d):
        st = d["steps"][4]
        st["V"] = st["V"][:-1]
    run("one glue reading dropped from page step 4", m6)

    def m7(d):
        used = json.dumps(d["rows"]).replace(" ", "")
        k = next(i for i, e in enumerate(d["extra"]) if len(e) > 12 and f'["amb",{i}]' in used)
        e = d["extra"][k]
        j = next(i for i, x in enumerate(e) if x)
        e[j] = int(e[j]) + 1
    run("one coefficient of a transported cross form that the rows use changed", m7)

    def m8(d):
        r = next(r for r in d["rows"] if r["type"] == "lemma")
        r["mult"] = "-" + r["mult"].lstrip("-")
    run("a lemma row with a negative multiplier", m8)
    return ok


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    which = args[0] if args else "both"
    certs = args[1] if len(args) > 1 else os.path.join(HERE, "certs")
    t0 = time.time()
    ok = True
    for w in ({"both": ["record", "allgroups"], "selftest": ["standard", "fivepage"]}.get(which, [which])):
        r = check_certificate(w, certs)
        print(f"RESULT {w}: {'PASS' if r else 'FAIL'}")
        ok &= r
        if "--controls" in sys.argv:
            c = controls(w, certs)
            print(f"RESULT {w} controls: {'PASS' if c else 'FAIL'}")
            ok &= c
    print(f"run time {time.time() - t0:.0f} s")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
