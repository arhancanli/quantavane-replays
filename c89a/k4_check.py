"""89a: exact monochromatic K4 densities of the step graphons behind four upper bounds on C_89 (Quantavane replays).

Our own code, from the definition on the 89a page: for an equal-weight step graphon W = M/den on n classes,

    F(W) = n^-4 sum over all ordered (a,b,c,d), repeats included, of prod_{6 pairs} W + prod_{6 pairs} (1 - W),

and C_89 <= F(W) (random colourings of blow-ups).  The certificates are read as data only; none of their programs
(lukevs/k4-ramsey's Lean and Python, JeetuSK0808's exact_eval.py, DennisKoshta's rebuild.py and verifier,
vidigoat's exact89.py) was read or run.

  k4_block.c  the count by a coarse block decomposition (see its header), modulo 2^64 and five primes;
  k4_dense.c  a second program: single roots by dense products modulo 2^61 - 1, no coarse structure;
  this file   reads the witnesses, rebuilds [Kos2026]'s 7680-class witness from VSJ's table and its amplitudes
              (FORMAT.md, sections 2 and 6 of that repository, and its stated dense SHA-256), runs both programs,
              reconstructs N_red and N_blue by the Chinese remainder theorem and compares the exact fraction with
              each claim.  Standard library only.

usage:
  python3 -I k4_check.py selftest                 brute force on random block-structured graphons, both programs
  python3 -I k4_check.py ROW CERTS [--roots=K] [--controls]
             ROW in vsj, sk, kos, vid, vid1; CERTS the directory fetch.sh fills; K roots per colour for k4_dense
"""
import array
import gzip
import hashlib
import json
import os
import random
import subprocess
import sys
import time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, 'bin')
CERTS = None
WORK = os.path.join(HERE, 'bin', 'tmp')        # scratch (ignored by git), removed after each row
MODS = [2 ** 64, 2 ** 61 - 1, 2 ** 62 - 57, 2 ** 60 - 93, 2 ** 59 - 55, 2 ** 63 - 25]
REC = 4          # the first four moduli reconstruct; the last two check

ROWS = {
    'vsj': dict(file='vsj/final3840.json.gz', block=20, label='[VSJ2026] 0.0301388876*',
                claim=Fraction(8450462766487926638466333426306607129, 280384030360880691940646801777885184000),
                printed='0.0301388876'),
    'sk': dict(file='sk/c89_stepgraphon_D20.json.gz', block=20, label='[SK2026] 0.0301388794*',
               claim=Fraction(27220838459434583306277336063848530765605851167,
                              903180180540490237232196590042258640992206848000),
               printed='0.0301388794'),
    'kos': dict(file=None, block=40, label='[Kos2026] 0.0301388857*',
                claim=Fraction(1126728298914655736634696330762614467, 37384537381450758925419573570384691200),
                printed='0.0301388857'),
    'vid1': dict(file='vid1/c89_stepgraphon_D30.json.gz', block=20, label='#270 first round 0.0301388617',
                 claim=Fraction(31383471736748006978678023758588478740355521988220578341177988039,
                                1041295852679825698551191435837073366434387539971249757275291648000),
                 printed='0.0301388618'),
    'vid': dict(file='vid/c89_stepgraphon_D30_round2.json.gz', block=20, label='#270 [Pat2026] 0.0301388585*',
                 claim=Fraction(31383468324407329838135970970371081393990727093313271428891812741,
                                1041295852679825698551191435837073366434387539971249757275291648000),
                 printed='0.0301388585'),
}
KOS_DELTA = 'kos/c89_sign3_crude_0.03013888570.delta.json.gz'


def crt(res):
    x, mod = 0, 1
    for r, m in zip(res, MODS[:REC]):
        t = ((r - x) * pow(mod, -1, m)) % m
        x += mod * t
        mod *= m
    for r, m in zip(res[REC:], MODS[REC:]):
        if x % m != r:
            raise SystemExit('check modulus %d disagrees: CRT value is wrong' % m)
    return x


def read_step_graphon(path):
    """rational-step-graphon-v1 (VSJ, SK, vidigoat): returns n, den, flat list of numerators (row-major)."""
    with gzip.open(path, 'rb') as f:
        doc = json.load(f)
    if doc.get('schema') != 'rational-step-graphon-v1':
        raise SystemExit('unexpected schema in %s' % path)
    den, w, R = doc['edge_probability_denominator'], doc['block_weights'], doc['red_probability_numerators']
    n = len(w)
    if len(set(w)) != 1 or type(w[0]) is not int or w[0] <= 0:
        raise SystemExit('block weights are not equal')
    if type(den) is not int or den <= 0 or len(R) != n or any(len(r) != n for r in R):
        raise SystemExit('bad shape or denominator')
    flat = array.array('I')
    for r in R:
        for v in r:
            if type(v) is not int or not 0 <= v <= den:
                raise SystemExit('entry out of range')
        flat.extend(r)
    return n, den, flat


def dense_sha256(flat):
    a = array.array('i', flat)        # FORMAT.md section 2: little-endian int32, row-major
    if sys.byteorder != 'little':
        a.byteswap()
    return hashlib.sha256(a.tobytes()).hexdigest()


def build_kos():
    """[Kos2026]: M'[2i+s, 2j+t] = M_ij + (-1)^(s+t) A_ij over VSJ's table M (FORMAT.md section 6), A from the
    delta's amplitude blocks [I, J, rows] (I < J, coarse groups of 20; A_ji = A_ij; A = 0 elsewhere)."""
    n0, den0, base = read_step_graphon(os.path.join(CERTS, ROWS['vsj']['file']))
    with gzip.open(os.path.join(CERTS, KOS_DELTA), 'rb') as f:
        delta = json.load(f)
    if delta.get('schema') != 'c89-sign-delta-v1':
        raise SystemExit('unexpected delta schema')
    b = delta['base']
    if b['n'] != n0 or b['den'] != den0 or dense_sha256(base) != b['dense_sha256']:
        raise SystemExit('the base table is not the one the delta names')
    bb = delta['base_block']
    A = {}
    for I, J, rows in delta['amplitudes']:
        if not (type(I) is int and type(J) is int and 0 <= I < J < n0 // bb) or len(rows) != bb:
            raise SystemExit('bad amplitude block')
        for r in range(bb):
            if len(rows[r]) != bb:
                raise SystemExit('bad amplitude row')
            for c in range(bb):
                v = rows[r][c]
                if type(v) is not int:
                    raise SystemExit('non-integer amplitude')
                if v:
                    i, j = I * bb + r, J * bb + c
                    if (i, j) in A:
                        raise SystemExit('amplitude given twice')
                    A[(i, j)] = v
                    A[(j, i)] = v
    n = 2 * n0
    out = array.array('I', bytes(4 * n * n))
    for i in range(n0):
        rowi = base[i * n0:(i + 1) * n0]
        for si in (0, 1):
            o = (2 * i + si) * n
            for j in range(n0):
                v = rowi[j]
                out[o + 2 * j] = v
                out[o + 2 * j + 1] = v
    for (i, j), v in A.items():
        mij = base[i * n0 + j]
        if abs(v) > min(mij, den0 - mij):
            raise SystemExit('amplitude above min(w, 1 - w)')
        for si in (0, 1):
            for tj in (0, 1):
                sign = 1 if (si + tj) % 2 == 0 else -1
                out[(2 * i + si) * n + 2 * j + tj] = mij + sign * v
    w = delta['witness']
    sha = dense_sha256(out)
    print('Kos witness: %d amplitude entries, n = %d, block %d; dense SHA-256 %s (%s)'
          % (len(A), n, w['block'], sha, 'equal to the delta' if sha == w['sha256_dense'] else 'DIFFERENT'))
    if sha != w['sha256_dense']:
        raise SystemExit('rebuilt witness differs from the stated dense hash')
    if w['claim'] != '%d/%d' % (ROWS['kos']['claim'].numerator, ROWS['kos']['claim'].denominator):
        raise SystemExit('the delta states a different claim')
    return n, den0, out


def write_bin(flat, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    a = array.array('I', flat)
    if sys.byteorder != 'little':
        a.byteswap()
    with open(path, 'wb') as f:
        a.tofile(f)


def run_block(path, n, s, den, roots=()):
    cmd = ['nice', '-n', '10', os.path.join(BIN, 'k4_block'), path, str(n), str(s), str(den)] + [str(r) for r in roots]
    out = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
    res = {'root': {}}
    for line in out.splitlines():
        t = line.split()
        if t[0] == 'TOTAL':
            res[int(t[2])] = [int(v) for v in t[3:]]
        elif t[0] == 'ROOT':
            res['root'][(int(t[2]), int(t[4]))] = [int(v) for v in t[5:]]
        elif t[0] == 'PART':
            res.setdefault('part', {})[(int(t[2]), t[3])] = [int(v) for v in t[4:]]
        elif t[0] == 'n':
            res['header'] = line
    return res


def run_dense(path, n, den, colour, roots):
    cmd = ['nice', '-n', '10', os.path.join(BIN, 'k4_dense'), path, str(n), str(den), str(colour)] + [str(r) for r in roots]
    out = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
    return {int(t[4]): int(t[5]) for t in (l.split() for l in out.splitlines()) if t and t[0] == 'DENSE'}


def density(res, n, den):
    nr, nb = crt(res[0]), crt(res[1])
    return nr, nb, Fraction(nr + nb, n ** 4 * den ** 6)


# ---------------------------------------------------------------------------------------------------- self-test
def brute(M, n, den, colour):
    X = [[(M[i * n + j] if colour == 0 else den - M[i * n + j]) for j in range(n)] for i in range(n)]
    per = []
    for a in range(n):
        Xa = X[a]
        t = 0
        for b in range(n):
            ab = Xa[b]
            if not ab:
                continue
            Xb = X[b]
            for c in range(n):
                abc = ab * Xa[c] * Xb[c]
                if not abc:
                    continue
                Xc = X[c]
                t += abc * sum(Xa[d] * Xb[d] * Xc[d] for d in range(n))
        per.append(t)
    return per


def random_instance(rng, m, s, den, pf):
    n = m * s
    frac = [[False] * m for _ in range(m)]
    for I in range(m):
        for J in range(I + 1, m):
            frac[I][J] = frac[J][I] = rng.random() < pf
    M = [0] * (n * n)
    for I in range(m):
        for J in range(I, m):
            if frac[I][J]:
                vals = [[rng.randrange(den + 1) for _ in range(s)] for _ in range(s)]
                if len({v for r in vals for v in r}) == 1:
                    vals[0][0] = (vals[0][0] + 1) % (den + 1)
            else:
                c = rng.choice((0, den))
                vals = [[c] * s for _ in range(s)]
            for r in range(s):
                for c in range(s):
                    M[(I * s + r) * n + J * s + c] = vals[r][c]
                    M[(J * s + c) * n + I * s + r] = vals[r][c]
    return n, M


def selftest():
    rng = random.Random(20261107)
    path = os.path.join(WORK, 'selftest.bin')
    cases = [(6, 3, 65536, 0.7), (7, 2, 2 ** 30, 0.8), (5, 4, 2 ** 20, 0.9), (6, 2, 17, 0.6), (5, 3, 2 ** 30, 1.0),
             (8, 2, 65536, 0.5), (4, 5, 1000, 0.75), (6, 3, 2 ** 30, 0.85)]
    for k, (m, s, den, pf) in enumerate(cases):
        n, M = random_instance(rng, m, s, den, pf)
        write_bin(M, path)
        roots = sorted(rng.sample(range(n), 3))
        res = run_block(path, n, s, den, roots)
        ok = True
        for colour in (0, 1):
            per = brute(M, n, den, colour)
            tot = sum(per)
            ok &= all(tot % q == r for q, r in zip(MODS, res[colour]))
            ok &= crt(res[colour]) == tot
            for a in roots:
                ok &= all(per[a] % q == r for q, r in zip(MODS, res['root'][(colour, a)]))
            dn = run_dense(path, n, den, colour, roots)
            ok &= all(dn[a] == per[a] % MODS[1] for a in roots)
        print('selftest %d: m=%d s=%d den=%d n=%d %s; %s' % (k, m, s, den, n, res['header'].split('den ')[1].split(' ', 1)[1],
                                                           'agree with brute force' if ok else 'DISAGREE'))
        if not ok:
            raise SystemExit('selftest failed')
    # a matrix the block method must refuse: a constant block that is neither 0 nor den
    n, M = random_instance(rng, 4, 2, 100, 0.0)
    for r in range(2):
        for c in range(2):
            M[(0 * 2 + r) * n + 2 + c] = M[(2 + c) * n + r] = 50
    write_bin(M, path)
    p = subprocess.run([os.path.join(BIN, 'k4_block'), path, str(n), '2', '100'], capture_output=True, text=True)
    print('selftest: constant block of 50/100 %s' % ('refused' if p.returncode else 'ACCEPTED'))
    if not p.returncode:
        raise SystemExit('selftest failed')
    os.remove(path)


# ---------------------------------------------------------------------------------------------------------- rows
def load_row(key):
    if key == 'kos':
        return build_kos()
    return read_step_graphon(os.path.join(CERTS, ROWS[key]['file']))


def decimal(fr, digits=25):
    q, r = divmod(fr.numerator * 10 ** digits, fr.denominator)
    s = str(q).rjust(digits + 1, '0')
    return s[:-digits] + '.' + s[-digits:]


def check_row(key, nroots, controls):
    row = ROWS[key]
    n, den, flat = load_row(key)
    s = row['block']
    path = os.path.join(WORK, key + '.bin')
    write_bin(flat, path)
    rng = random.Random('roots-' + key)
    roots = sorted(rng.sample(range(n), nroots))
    res = run_block(path, n, s, den, roots)
    nr, nb, F = density(res, n, den)
    print('%s: %s' % (row['label'], res['header']))
    print('  N_red  = %d' % nr)
    print('  N_blue = %d' % nb)
    print('  F = %d/%d = %s...' % (F.numerator, F.denominator, decimal(F)))
    same = F == row['claim']
    printed = Fraction(row['printed'])
    print('  stated fraction: %s; printed upper bound %s is %s the exact value'
          % ('EQUAL' if same else 'DIFFERENT', row['printed'], 'above' if printed > F else 'NOT above'))
    for colour in (0, 1):
        dn = run_dense(path, n, den, colour, roots)
        for a in roots:
            bl = res['root'][(colour, a)][1]
            print('  root %d colour %d: k4_block %d, k4_dense %d mod 2^61-1: %s'
                  % (a, colour, bl, dn[a], 'agree' if bl == dn[a] else 'DISAGREE'))
            if bl != dn[a]:
                raise SystemExit('the two programs disagree')
    if not same or printed <= F:
        raise SystemExit('claim not confirmed')
    if controls:
        # red controls: one fractional entry (and its mirror) raised by 1 must change the fraction; a claim one
        # unit lower in its last printed digit must not be implied
        flat2 = array.array('I', flat)
        I, J = None, None
        for idx in range(n * n):
            i, j = divmod(idx, n)
            if i < j and 0 < flat2[idx] < den and flat2[idx] != flat2[i * n + (j // s) * s + ((j % s) + 1) % s]:
                I, J = i, j
                break
        flat2[I * n + J] += 1
        flat2[J * n + I] += 1
        write_bin(flat2, path)
        r2 = run_block(path, n, s, den)
        _, _, F2 = density(r2, n, den)
        print('  control: entry (%d,%d) raised by 1/den: F changes by %.3e, %s the stated fraction'
              % (I, J, float(F2 - F), 'FAILS' if F2 != row['claim'] else 'still equals'))
        if F2 == row['claim']:
            raise SystemExit('control not caught')
        flat2[I * n + J] -= 1          # (J, I) keeps the +1
        write_bin(flat2, path)
        p = subprocess.run([os.path.join(BIN, 'k4_block'), path, str(n), str(s), str(den)], capture_output=True, text=True)
        print('  control: the same entry changed on one side only: %s' % ('refused (' + p.stderr.strip() + ')' if p.returncode else 'ACCEPTED'))
        if not p.returncode:
            raise SystemExit('control not caught')
        lower = printed - Fraction(1, 10 ** (len(row['printed']) - 2))
        print('  control: the claim one unit lower in its last digit, %s, is %s'
              % (decimal(lower, len(row['printed']) - 2), 'not implied (below F)' if lower < F else 'IMPLIED'))
        if lower >= F:
            raise SystemExit('control not caught')
    os.remove(path)
    return F


def main():
    global CERTS
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    os.makedirs(WORK, exist_ok=True)
    t0 = time.time()
    if sys.argv[1] == 'selftest':
        selftest()
        print('RESULT selftest: all cases agree with brute force; the 50/100 constant block is refused')
    elif sys.argv[1] in ROWS and len(sys.argv) >= 3:
        CERTS = sys.argv[2]
        nroots = 2
        for a in sys.argv[3:]:
            if a.startswith('--roots='):
                nroots = int(a.split('=')[1])
        F = check_row(sys.argv[1], nroots, '--controls' in sys.argv)
        print('RESULT %s: PASS, F = %s... exactly the stated fraction; printed bound %s holds'
              % (ROWS[sys.argv[1]]['label'], decimal(F, 16), ROWS[sys.argv[1]]['printed']))
    else:
        raise SystemExit(__doc__)
    print('run time %.0f s' % (time.time() - t0))


if __name__ == '__main__':
    main()
