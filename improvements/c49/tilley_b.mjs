// Program B (JavaScript, node >= 18): exact check of R. F. Tilley Jr.'s 54-point block (PR #379, 49a), written
// separately from program A (no shared code, nothing imported from our 4 November checker). Everything is rebuilt from
// page 2 of the note: K5 edge families, the Z_11 family H, the pendant gluing, the block B, the two halves and the
// three families G, P x L(Ap), L(Aq) x P.
//
// Method: every family here is a union of explicit 'cells' (core mask, side of pendant 1, side of pendant 2), each
// side one of: HA (the 46 members of H avoiding 0), HC (Y - {0} for the 75 members containing 0), E (only the empty
// set). No cells are merged. Separation patterns of member triples are computed by brute force on each pendant side
// triple (all |S1||S2||S3| triples) and directly on core masks; a half's pattern set is the OR-closure over all cell
// triples; the 54-point check runs over every multiset of three top boxes. Containment likewise.
// Rate: BigInt only. rho = a + sqrt(b c); a printed decimal m/10^k is a valid lower bound iff (m/10^k)^52 <= rho.
// Usage: node tilley_b.mjs [control]   controls: plant-sunflower, drop-chair, swap12 (exit 1 expected)
const control = process.argv[2] || '';
const say = (s) => console.log(s);

// ---- K5 edges, lexicographic: bit i <-> edge i
const EDGES = []; for (let a = 0; a < 5; a++) for (let b = a + 1; b < 5; b++) EDGES.push([a, b]);
const eid = (a, b) => EDGES.findIndex(([x, y]) => x === Math.min(a, b) && y === Math.max(a, b));
const emask = (pairs) => pairs.reduce((m, [a, b]) => m | (1 << eid(a, b)), 0);
const ALL = (1 << 10) - 1;
const P_EDGE = eid(0, 1), Q_EDGE = eid(0, 2), E1 = eid(0, 3), E2 = eid(0, 4);

// F: stars, K_{2,3}, complements of chairs
const F = new Set();
for (let v = 0; v < 5; v++) F.add(emask([0, 1, 2, 3, 4].filter((w) => w !== v).map((w) => [v, w])));
for (let a = 0; a < 5; a++) for (let b = a + 1; b < 5; b++) {
  const rest = [0, 1, 2, 3, 4].filter((w) => w !== a && w !== b);
  F.add(emask(rest.flatMap((w) => [[a, w], [b, w]])));
}
let chairs = 0;
for (let c = 0; c < 5; c++) for (let e = 0; e < 5; e++) if (e !== c) {
  const nb = [0, 1, 2, 3, 4].filter((w) => w !== c && w !== e);
  for (const a of nb) {
    if (control === 'drop-chair' && chairs === 0) { chairs++; continue; }
    F.add(ALL & ~emask([...nb.map((w) => [c, w]), [a, e]])); chairs++;
  }
}
if (F.size !== 75 && control !== 'drop-chair') throw new Error('F has ' + F.size);
// the chair with edges 01, 02, 03, 14 is one of them (the note's example)
const exampleChair = ALL & ~emask([[0, 1], [0, 2], [0, 3], [1, 4]]);

// H on Z_11: images of three base sets under x -> a x + b, a a nonzero square
const sq = [...new Set([1, 2, 3, 4, 5, 6, 7, 8, 9, 10].map((x) => (x * x) % 11))];
const Hs = new Set();
for (const base of [[0, 1, 2, 3, 4, 5, 6], [0, 1, 2, 3, 4, 6, 7], [0, 1, 2, 6, 9]])
  for (const a of sq) for (let b = 0; b < 11; b++) Hs.add(base.reduce((m, x) => m | (1 << ((a * x + b) % 11)), 0));
const H = [...Hs];
// pendant coordinates: the points 1..10 of Z_11 as bits 0..9
const HA = H.filter((m) => !(m & 1)).map((m) => m >> 1);
const HC = H.filter((m) => m & 1).map((m) => m >> 1);
const E = [0];
const FULL = [(1 << 10) - 1];
const SIDES = control === 'koshta-extra' || control === 'koshta-extra-pinho' ? { HA, HC, E, FULL } : { HA, HC, E };
// row 54 variation 'extras': env TILLEY_EXTRAS = JSON list of 28-bit masks in the coordinates of the PR #260/#370 files
// (bits 0..7 the core edges 01,02,12,13,14,23,24,34; bits 8..17 pendant 1 points 1..10; bits 18..27 pendant 2)
const KOSCORE = [[0, 1], [0, 2], [1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4]];
const EXTRA_CELLS = [];
if (control === 'extras') {
  for (const m of JSON.parse(process.env.TILLEY_EXTRAS || '[]')) {
    const core = KOSCORE.reduce((acc, [u, v], i) => (Math.floor(m / 2 ** i) % 2 ? acc | (1 << eid(u, v)) : acc), 0);
    const p1 = Math.floor(m / 2 ** 8) % 1024, p2 = Math.floor(m / 2 ** 18) % 1024;
    const k1 = p1 ? 'S' + p1 : 'E', k2 = p2 ? 'S' + p2 : 'E';
    if (p1) SIDES[k1] = [p1]; if (p2) SIDES[k2] = [p2];
    EXTRA_CELLS.push([core, k1, k2]);
  }
}

// cells: [core10, side1, side2], core with the glued edges e1, e2 removed
const cellsL = (X, glue1, glue2) => X.map((c) => [c & ~((1 << E1) | (1 << E2)),
  glue1 ? (c >> E1 & 1 ? 'HA' : 'HC') : 'E', glue2 ? (c >> E2 & 1 ? 'HA' : 'HC') : 'E']);
const A1 = [251, 425, 499, 506, 570, 707, 761, 875, 883, 889, 923, 947, 985, 986, 1002];
const A2 = [358, 471, 499, 501, 631, 707, 757, 758, 789, 883, 935, 947, 950, 982, 997];
const permute = (m, perm) => EDGES.reduce((acc, [a, b], i) => (m >> i & 1 ? acc | (1 << eid(perm[a], perm[b])) : acc), 0);
const Bcells = [...cellsL([...F], true, true), ...cellsL(A1, false, true), ...cellsL(A2, true, false)];
// row 54 variation: Koshta's extra set of PR #260 (core edges 12, 14, 23; pendant 1 empty; pendant 2 all ten points)
// added to the block B (both copies)
if (control === 'koshta-extra') Bcells.push([emask([[1, 2], [1, 4], [2, 3]]), 'E', 'FULL']);
for (const cell of EXTRA_CELLS) Bcells.push(cell);
const swapPerm = (u, v) => [0, 1, 2, 3, 4].map((x) => (x === u ? v : x === v ? u : x));
const Ap = A1.map((m) => permute(m, swapPerm(1, 3)));
const Aq = A1.map((m) => permute(m, control === 'swap12' ? swapPerm(1, 2) : swapPerm(2, 3)));
const LAp = cellsL(Ap, true, true), LAq = cellsL(Aq, true, true);
const P = [[0, 'HC', 'E'], [0, 'E', 'HC']];
const size = (cells) => cells.reduce((s, [, a, b]) => s + SIDES[a].length * SIDES[b].length, 0);
const portCount = (cells, bitA, bitB) => { // counts by (bitA in core, bitB in core)
  const N = [[0, 0], [0, 0]];
  for (const [c, a, b] of cells) N[c >> bitA & 1][c >> bitB & 1] += SIDES[a].length * SIDES[b].length;
  return N;
};
const NB = portCount(Bcells, P_EDGE, Q_EDGE);
say(`B: ${size(Bcells)} members, cells ${Bcells.length}; (p,q) counts 00 ${NB[0][0]} 01 ${NB[0][1]} 10 ${NB[1][0]} 11 ${NB[1][1]}`);
say(`F ${F.size} (the example chair's complement in F: ${F.has(exampleChair)}), H ${H.length} (HA ${HA.length}, HC ${HC.length}), |P| ${size(P)}, |L(Ap)| ${size(LAp)}, |L(Aq)| ${size(LAq)}`);
const avoidsP = LAp.every(([c]) => !(c >> P_EDGE & 1)), avoidsQ = LAq.every(([c]) => !(c >> Q_EDGE & 1));
say(`L(Ap) avoids p: ${avoidsP}; L(Aq) avoids q: ${avoidsQ}`);

// halves
const L0 = Bcells.filter(([c]) => !(c >> Q_EDGE & 1));
const L1 = Bcells.filter(([c]) => c >> Q_EDGE & 1).map(([c, a, b]) => [c & ~(1 << Q_EDGE), a, b]);
const R0 = Bcells.filter(([c]) => !(c >> P_EDGE & 1));
const R1 = Bcells.filter(([c]) => c >> P_EDGE & 1).map(([c, a, b]) => [c & ~(1 << P_EDGE), a, b]);
const TOP = [['G L0xR1', L0, R1], ['G L1xR0', L1, R0], ['PxL(Ap)', P, LAp], ['L(Aq)xP', LAq, P]];
if (control === 'plant-sunflower') {
  // the member (empty left, the first L(Ap) member's core and sides) with two members of P x L(Ap) whose left parts
  // lie in different pendants: a 3-sunflower (also a containment)
  TOP.push(['planted', [[0, 'E', 'E']], [LAp[0]]]);
}

// ---- patterns
const sep = (x, y, z) => (x !== y ? 1 : 0) | (x !== z ? 2 : 0) | (y !== z ? 4 : 0);
const sideKeys = Object.keys(SIDES);
const PEND = {};
for (const s1 of sideKeys) for (const s2 of sideKeys) for (const s3 of sideKeys) {
  let m = 0;
  for (const x of SIDES[s1]) for (const y of SIDES[s2]) { const k = x & y; for (const z of SIDES[s3]) if ((x & z) === k && (y & z) === k) m |= 1 << sep(x, y, z); }
  PEND[s1 + s2 + s3] = m;
}
const REL = {};
for (const s1 of sideKeys) for (const s2 of sideKeys) {
  let r = 0; for (const x of SIDES[s1]) for (const y of SIDES[s2]) if ((x & y) === x) r |= x === y ? 1 : 2;
  REL[s1 + s2] = r;
}
const orClose = (ms) => { let reach = 1; for (const m of ms) { let nr = 0; for (let r = 0; r < 8; r++) if (reach >> r & 1) for (let s = 0; s < 8; s++) if (m >> s & 1) nr |= 1 << (r | s); reach = nr; } return reach; };
const halfPatterns = (X, Y, Z) => {
  let reach = 0;
  for (const [cx, x1, x2] of X) for (const [cy, y1, y2] of Y) {
    const k = cx & cy;
    for (const [cz, z1, z2] of Z) {
      if ((cx & cz) !== k || (cy & cz) !== k) continue;
      const mc = 1 << sep(cx, cy, cz), m1 = PEND[x1 + y1 + z1], m2 = PEND[x2 + y2 + z2];
      if (!m1 || !m2) continue;
      reach |= orClose([mc, m1, m2]);
    }
  }
  return reach;
};
const halfRel = (X, Y) => {
  let res = 0;
  for (const [cx, x1, x2] of X) for (const [cy, y1, y2] of Y) {
    if ((cx & cy) !== cx) continue;
    const rc = cx === cy ? 1 : 2, r1 = REL[x1 + y1], r2 = REL[x2 + y2];
    if (!r1 || !r2) continue;
    if ((rc & 1) && (r1 & 1) && (r2 & 1)) res |= 1;
    if ((rc & 2) || (r1 & 2) || (r2 & 2)) res |= 2;
  }
  return res;
};
// note on halfRel: for one cell pair, a proper containment on the half needs every block to allow containment and one
// block to allow a proper one; r1, r2 nonzero and the core contained give the first; the second is (rc|r1|r2) & 2.
let sunflowers = 0, containments = 0;
for (let i = 0; i < TOP.length; i++) for (let j = i; j < TOP.length; j++) for (let k = j; k < TOP.length; k++) {
  const pl = halfPatterns(TOP[i][1], TOP[j][1], TOP[k][1]), pr = halfPatterns(TOP[i][2], TOP[j][2], TOP[k][2]);
  let hit = false;
  for (let a = 0; a < 8; a++) for (let b = 0; b < 8; b++) if ((pl >> a & 1) && (pr >> b & 1) && (a | b) === 7) hit = true;
  if (hit) { sunflowers++; say(`  3-sunflower with members from ${TOP[i][0]}, ${TOP[j][0]}, ${TOP[k][0]}`); }
}
for (let i = 0; i < TOP.length; i++) for (let j = 0; j < TOP.length; j++) {
  const rl = halfRel(TOP[i][1], TOP[j][1]), rr = halfRel(TOP[i][2], TOP[j][2]);
  if (rl && rr && ((rl & 2) || (rr & 2))) { containments++; say(`  containment: a member of ${TOP[i][0]} inside one of ${TOP[j][0]}`); }
}
say(`exhaustive over all multisets of three top boxes and all cell triples: 3-sunflowers ${sunflowers}, containments ${containments}`);

// ---- counts and rate (BigInt)
const M = [[0n, 0n], [0n, 0n]];
for (const [name, X, Y] of TOP) {
  const nx = [0, 0], ny = [0, 0];
  for (const [c, a, b] of X) nx[c >> P_EDGE & 1] += SIDES[a].length * SIDES[b].length;   // left: surviving port p
  for (const [c, a, b] of Y) ny[c >> Q_EDGE & 1] += SIDES[a].length * SIDES[b].length;   // right: surviving port q
  const col = [[nx[0] * ny[0], nx[0] * ny[1]], [nx[1] * ny[0], nx[1] * ny[1]]];
  say(`  ${name}: 00 ${col[0][0]} 01 ${col[0][1]} 10 ${col[1][0]} 11 ${col[1][1]}`);
  for (let x = 0; x < 2; x++) for (let y = 0; y < 2; y++) M[x][y] += BigInt(col[x][y]);
}
const total = M[0][0] + M[0][1] + M[1][0] + M[1][1];
say(`54-point block: ${total} members; 00 ${M[0][0]} 01 ${M[0][1]} 10 ${M[1][0]} 11 ${M[1][1]}`);
const a = M[0][1], b = M[0][0], c = M[1][1];
const sym = M[0][1] === M[1][0];
const bc = b * c;
say(`a = N01 = N10: ${sym}; rho = ${a} + sqrt(${bc})`);
const validLower = (str) => { // (m/10^k)^52 <= a + sqrt(bc)
  const [w, f] = str.split('.'); const k = BigInt(f.length); const m = BigInt(w + f); const D = 10n ** (52n * k);
  const lhs = m ** 52n - a * D; return lhs <= 0n || lhs * lhs <= bc * D * D;
};
// 20 decimals of rho^(1/52) by bisection on integers
let lo = 1567800000000000000000n, hi = 1567900000000000000000n; // scaled by 10^21
const scale = '000000000000000000000';
const okAt = (m) => { const s = m.toString(); return validLower(s.slice(0, 1) + '.' + s.slice(1)); };
while (hi - lo > 1n) { const mid = (lo + hi) / 2n; if (okAt(mid)) lo = mid; else hi = mid; }
const ls = lo.toString();
say(`rho^(1/52) in [${ls[0]}.${ls.slice(1)}, +1e-21)`);
const p1 = validLower('1.567866388204975'), p2 = validLower('1.56786638820497592522'), p3 = validLower('1.56786638820497592523');
say(`printed 1.567866388204975 valid: ${p1}; 1.56786638820497592522 valid: ${p2}; 1.56786638820497592523 valid: ${p3} (the note says beta < ...523)`);
const countsOK = size(Bcells) === 252271 && NB[0][0] === 33948 && NB[0][1] === 64728 && NB[1][0] === 64728 && NB[1][1] === 88867 &&
  size(P) === 150 && size(LAp) === 42412 && size(LAq) === 42412 && M[0][0] === 4397311488n && a === 7211663100n && sym &&
  c === 11504366352n && bc === 50588282321810251776n;
say(`counts as printed on page 2 and the radicand as in the Theorem: ${countsOK}`);
const pass = control === 'koshta-extra' || control === 'extras' ? (avoidsP && avoidsQ && !sunflowers && !containments)
  : countsOK && avoidsP && avoidsQ && !sunflowers && !containments && p1 && p2 && !p3;
if (control === 'extras') say(`VARIATION extras (${EXTRA_CELLS.length} cells): valid ${avoidsP && avoidsQ && !sunflowers && !containments}`);
if (control === 'koshta-extra') say(`VARIATION koshta-extra: valid ${avoidsP && avoidsQ && !sunflowers && !containments}; 1.5678667956347 valid lower bound: ${validLower('1.5678667956347')}; 1.5678667956348: ${validLower('1.5678667956348')}`);
say(`RESULT ${pass ? 'PASS' : 'FAIL'}`);
process.exit(pass ? 0 : 1);
