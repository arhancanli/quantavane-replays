"""Tabulate the check logs: one row per witness, both programs' verdicts and CPU seconds.

python3 -I summarize.py LOG_DIR [CONTROL_LOG_DIR]
A row passes when ramsey_check reports neither clique (exit 0) and ramsey_indep's largest cliques are below S and T.
Exit status 152 is the CPU cap (SIGXCPU) and 143 a stop by hand (SIGTERM): the run ended unfinished, which\ndecides nothing.
"""
import glob
import os
import re
import sys

# name: (bound, revision-18 value, source, s, t)
META = {
    'nart_R3_17_n92': ('R(3,17) >= 93', '92', 'NaRT repository (22 Oct 2026)', 3, 17),
    'nart_R4_15_n162': ('R(4,15) >= 163', '159', 'NaRT repository (22 Oct 2026)', 4, 15),
    'nart_R4_16_n173': ('R(4,16) >= 174', '170', 'NaRT, arXiv:2603.09172 v3-v5; repository (25 Mar 2026)', 4, 16),
    'nart_R4_18_n208': ('R(4,18) >= 209', '205', 'NaRT, arXiv:2603.09172 v3-v5; repository (2 Apr 2026)', 4, 18),
    'nart_R4_19_n218': ('R(4,19) >= 219', '213', 'NaRT, arXiv:2603.09172 v3-v5; repository (19 Apr 2026)', 4, 19),
    'nart_R4_20_n236': ('R(4,20) >= 237', '234', 'NaRT, arXiv:2603.09172 v3-v5; repository (19 Apr 2026)', 4, 20),
    'nart_R4_21_n243': ('R(4,21) >= 244', '242', 'NaRT repository (11 May 2026)', 4, 21),
    'nart_R7_13_n512': ('R(7,13) >= 513', '511', 'NaRT repository (22 Oct 2026)', 7, 13),
    'nart_R8_12_n510': ('R(8,12) >= 511', '- (471 implied)', 'NaRT repository (22 Oct 2026)', 8, 12),
    'nart_R8_18_n1051': ('R(8,18) >= 1052', '1049', 'NaRT repository (22 Oct 2026)', 8, 18),
    'nart_R9_11_n616': ('R(9,11) >= 617', '- (608 implied)', 'NaRT repository (22 Oct 2026)', 9, 11),
    'nart_R10_12_n855': ('R(10,12) >= 856', '- (846 implied)', 'NaRT repository (22 Oct 2026)', 10, 12),
    'yu_R4_20_n251': ('R(4,20) >= 252', '234', 'Yu, arXiv:2608.18169 (17 Aug 2026); Zenodo 21941603', 4, 20),
    'yu_R4_20_n251_Gprime': ('R(4,20) >= 252 (G\')', '234', 'Yu, same paper: the copy x -> 116x', 4, 20),
    'cremer_R6_8_n134': ('R(6,8) >= 135', '134', 'Cremer, arXiv:2610.12122 (8 Oct 2026); repository v1.1.0', 6, 8),
    'cremer_R8_10_n344_1': ('R(8,10) >= 345', '343', 'Cremer, arXiv:2610.12122; repository v1.1.0, witness 1', 8, 10),
    'cremer_R8_10_n344_2': ('R(8,10) >= 345', '343', 'Cremer, arXiv:2610.12122; repository v1.1.0, witness 2', 8, 10),
    'tatarevic_R12_12_n1640': ('R(12,12) >= 1641', '1640', 'Tatarevic, arXiv:2608.06531 (6 Aug 2026); repository', 12, 12),
    'qv_R5_15_n278': ('R(5,15) >= 279', '275', 'Canli, Zenodo 23857057 (1 Nov 2026); quantavane-ramsey-records', 5, 15),
    'qv_R8_12_n518': ('R(8,12) >= 519', '- (471 implied)', 'Canli, Zenodo 23920916 (3 Nov 2026); quantavane-ramsey-records', 8, 12),
}
CONIGLIO = {24: (151, 143), 25: (160, 154), 26: (163, 161), 28: (181, 179), 29: (194, 190), 30: (199, 197),
            31: (210, 208), 32: (219, 217), 33: (230, 227), 34: (239, 234), 35: (253, 248), 36: (261, 255),
            37: (272, 267), 38: (281, 278), 39: (294, 290), 40: (305, 298), 41: (317, 311), 42: (327, 320),
            43: (340, 333), 44: (348, 339), 45: (364, 354), 46: (371, 362), 47: (385, 380), 48: (391, 384),
            49: (411, 402)}
for l, (b, old) in CONIGLIO.items():
    META['coniglio_R3_%d_n%d' % (l, b - 1)] = ('R(3,%d) >= %d' % (l, b), str(old),
                                               'Coniglio et al., arXiv:2608.18769 (19 Aug 2026); repository', 3, l)


FILES = {
    'nart_R3_17_n92': 'nart/R3_17_93.txt', 'nart_R4_15_n162': 'nart/R4_15_163.txt',
    'nart_R4_16_n173': 'nart/R4_16_174.txt', 'nart_R4_18_n208': 'nart/R4_18_209.txt',
    'nart_R4_19_n218': 'nart/R4_19_219.txt', 'nart_R4_20_n236': 'nart/R4_20_237.txt',
    'nart_R4_21_n243': 'nart/R4_21_244.txt', 'nart_R7_13_n512': 'nart/R7_13_513.txt',
    'nart_R8_12_n510': 'nart/R8_12_511.txt', 'nart_R8_18_n1051': 'nart/R8_18_1052.txt',
    'nart_R9_11_n616': 'nart/R9_11_617.txt', 'nart_R10_12_n855': 'nart/R10_12_856.txt',
    'yu_R4_20_n251': 'yu/CERTIFICATE_R4_20.txt', 'yu_R4_20_n251_Gprime': 'yu/CERTIFICATE_R4_20.txt',
    'cremer_R6_8_n134': 'cremer/n134-gpt-6.1-sol.txt', 'cremer_R8_10_n344_1': 'cremer/n344-gpt-6.1-sol-1.txt',
    'cremer_R8_10_n344_2': 'cremer/n344-gpt-6.1-sol-2.txt', 'tatarevic_R12_12_n1640': 'tatarevic/r1212_1640.graph',
    'qv_R5_15_n278': 'quantavane/R5_15_K278.txt', 'qv_R8_12_n518': 'quantavane/R8_12_K518.txt',
}
for l, (b, old) in CONIGLIO.items():
    FILES['coniglio_R3_%d_n%d' % (l, b - 1)] = 'coniglio/R_3_%d_graphSize_%d.clq' % (l, b - 1)
SUMS = {}
_sums = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'SHA256SUMS')
if os.path.exists(_sums):
    for line in open(_sums):
        h, f = line.split()
        SUMS[f] = h


def parse(path):
    if not os.path.exists(path):
        return None
    t = open(path, encoding='utf-8', errors='replace').read()
    r = {'text': t}
    m = re.search(r'# summary \S+ \S+ exit (\d+) cpu_user ([\d.]+) cpu_sys ([\d.]+) wall ([\d.]+)', t)
    if not m:
        return {'running': True}
    r['exit'] = int(m.group(1))
    r['cpu'] = float(m.group(2)) + float(m.group(3))
    m = re.search(r'# matrix \S+ sha256 (\w+)', t)
    r['sha'] = m.group(1) if m else ''
    m = re.search(r'^n (\d+) red_K(\d+) (none|found) blue_K(\d+) (none|found)', t, re.M)
    if m:
        r['n'] = int(m.group(1))
        r['red'], r['blue'] = m.group(3), m.group(5)
    m1 = re.search(r"entries '1': .*largest clique (\d+)", t)
    m0 = re.search(r"entries '0': .*largest clique (\d+)", t)
    if m1 and m0:
        r['w1'], r['w0'] = int(m1.group(1)), int(m0.group(1))
    m = re.search(r'^n (\d+) symmetric \d circulant (\d)', t, re.M)
    if m:
        r['n'] = int(m.group(1))
        r['circ'] = int(m.group(2))
    return r


def fmt_check(r):
    if r is None:
        return 'not run', None
    if r.get('running'):
        return 'running', None
    if r['exit'] == 152:
        return 'stopped at the %d s cap' % round(r['cpu']), False
    if r['exit'] == 143:
        return 'stopped by us after %d s (budget)' % round(r['cpu']), False
    if 'red' in r:
        ok = r['red'] == 'none' and r['blue'] == 'none'
        return ('no red, no blue clique' if ok else 'clique FOUND') + ', %.1f s' % r['cpu'], ok
    return 'exit %d' % r['exit'], False


def fmt_nbhd(r, colour, k):
    if r is None:
        return 'not run', None
    if r.get('running'):
        return 'running', None
    if r['exit'] == 152:
        return 'stopped at the %d s cap' % round(r['cpu']), False
    if r['exit'] == 143:
        return 'stopped by us after %d s (budget)' % round(r['cpu']), False
    ok = r.get('red') == 'none' and r.get('blue') == 'none'
    return ('no %s K%d' % (colour, k) if ok else '%s K%d FOUND' % (colour, k)) + ', %.1f s' % r['cpu'], ok


def fmt_indep(r, s, t):
    if r is None:
        return 'not run', None
    if r.get('running'):
        return 'running', None
    if r['exit'] == 152:
        return 'stopped at the %d s cap' % round(r['cpu']), False
    if r['exit'] == 143:
        return 'stopped by us after %d s (budget)' % round(r['cpu']), False
    if 'w1' in r:
        ok = r['w1'] < s and r['w0'] < t
        return 'largest %d / %d%s, %.1f s' % (r['w1'], r['w0'], ' (from vertex 0)' if r.get('circ') else '',
                                              r['cpu']), ok
    return 'exit %d' % r['exit'], False


def main():
    logs = sys.argv[1]
    total = 0.0
    for f in glob.glob(os.path.join(logs, '**', '*.log'), recursive=True):
        r = parse(f)
        if r and not r.get('running'):
            total += r['cpu']
    print('| Bound | Rev. 18 | Source | File SHA-256 | n | ramsey_check | ramsey_indep (largest red / blue) | Result |')
    print('| --- | --- | --- | --- | --- | --- | --- | --- |')
    order = sorted(META, key=lambda k: (META[k][3], META[k][4], k))
    for name in order:
        bound, old, src, s, t = META[name]
        t_ = t
        chk = parse(os.path.join(logs, name + '.check.log'))
        ind = parse(os.path.join(logs, name + '.indep.log'))
        prior_c = os.path.join(logs, name + '.2026-11-03.check.txt')
        if chk is None and os.path.exists(prior_c):  # our own R(8,12) >= 519: the 3 Nov runs, same file
            t = open(prior_c).read()
            m = re.search(r'^n (\d+) red_K\d+ (none|found) blue_K\d+ (none|found) nodes (\d+) seconds ([\d.]+)', t, re.M)
            u = open(os.path.join(logs, name + '.2026-11-03.indep.txt')).read()
            w = [int(x) for x in re.findall(r'largest clique (\d+)', u)]
            ok = m.group(2) == m.group(3) == 'none' and w[0] < s and w[1] < t_
            print('| %s | %s | %s | `%s...%s` | %s | 3 Nov 2026 run, same file: no red, no blue clique, %s nodes, %.0f s wall '
                  '| 3 Nov 2026 run: largest %d / %d | %s (not re-run) |'
                  % (bound, old, src, SUMS[FILES[name]][:8], SUMS[FILES[name]][-4:], m.group(1), m.group(4),
                     float(m.group(5)), w[0], w[1], 'PASS' if ok else 'FAIL'))
            continue
        n1 = parse(os.path.join(logs, name + '.nbhd1.check.log'))
        n0 = parse(os.path.join(logs, name + '.nbhd0.check.log'))
        perm = parse(os.path.join(logs, name + '.perm.indep.log'))
        cells = []
        oks = []
        if n0 is not None or n1 is not None:
            a, oa = fmt_nbhd(n1, 'red', s - 1)
            b, ob = fmt_nbhd(n0, 'blue', t - 1)
            cells.append('from vertex 0: N_red(0) %s; N_blue(0) %s' % (a, b))
            oks += [oa, ob]
            if chk is not None:
                c, _ = fmt_check(chk)
                cells[-1] += '; all vertices: %s' % c
        else:
            c, oc = fmt_check(chk)
            cells.append(c)
            oks.append(oc)
        d, od = fmt_indep(ind, s, t)
        if perm is not None:
            e, oe = fmt_indep(perm, s, t)
            d += '; renumbered copy, all vertices: ' + e
            oks.append(oe)
        oks.append(od)
        n = int(re.search(r'_n(\d+)', name).group(1))
        result = 'PASS' if all(o is True for o in oks) else ('FAIL' if any(o is False and 'FOUND' in str(cells)
                                                                          for o in oks) else 'not finished')
        h = SUMS.get(FILES.get(name, ''), '')
        print('| %s | %s | %s | `%s...%s` | %s | %s | %s | %s |' % (bound, old, src, h[:8], h[-4:], n, cells[0], d,
                                                                   result))
    print()
    print('CPU seconds in these logs: %.0f (%.2f core-hours)' % (total, total / 3600))
    if len(sys.argv) > 2:
        print('\nRed controls (one edge flipped; each program must report the forbidden clique):')
        for f in sorted(glob.glob(os.path.join(sys.argv[2], '*.ctrl_*.log'))):
            r = parse(f)
            base = os.path.basename(f)[:-4]
            name, ctl, prog = base.rsplit('.', 2)[0], base.rsplit('.', 2)[1], base.rsplit('.', 2)[2]
            s, t = META[name][3], META[name][4]
            if prog == 'check':
                seen = 'found' in (r.get('red'), r.get('blue'))
                print('  %s %s ramsey_check: red %s, blue %s (exit %d) -> %s'
                      % (name, ctl, r.get('red'), r.get('blue'), r['exit'], 'fails as it should' if seen else 'MISSED'))
            else:
                seen = r.get('w1', 0) >= s or r.get('w0', 0) >= t
                print('  %s %s ramsey_indep: largest %s / %s -> %s'
                      % (name, ctl, r.get('w1'), r.get('w0'), 'fails as it should' if seen else 'MISSED'))


if __name__ == '__main__':
    main()
