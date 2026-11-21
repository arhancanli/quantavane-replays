// Program B (JavaScript, node >= 18), written separately from program A (lys_a.c); no shared code.
// Exact check of the six colourings of Lysenstoen, arXiv 2611.25143, from the ancillary JSON files read as data.
// Different reductions from program A:
//   triples      the colouring is REBUILT from the group and orbit representatives of the JSON 'rule' and compared with
//                the stored string (and its stated SHA-256); invariance under x -> x+1 and x -> g x (mod m, fixed
//                points stay) is checked on every triple; then every forbidden set THROUGH THE POINT 0 is searched for
//                (any set meets Z_m, and a translation moves that point to 0).
//   5-subsets    the full colouring is built from the short form with the normalisation (largest, second largest) ->
//                (0, 1), compared by SHA-256 with the stated digest of the full colouring, checked for invariance under
//                x -> x+1 and x -> 2x on every 5-subset; then every 6-set through {0, 1} (AGL(1,83) is 2-transitive).
//   4-subsets    600 points: the table is constant on cross-ratio classes (orbits of x -> 1/x, x -> 1-x), so the
//                colouring is PGL(2,599)-invariant; PGL is 3-transitive, so every 6-set through {0, 1, 2} is searched.
//                200 points (PSL rule): invariance under x -> x+1 and x -> -1/x on every 4-subset, the colouring's
//                SHA-256 against the stated digest, transitivity on 3-subsets by union-find over the generators, and
//                every 5-set and 6-set through {0, 1, 2}.
// Usage: node lys_b.mjs <lit dir> [case]   cases: t37 t70 t98 f83 p200 p600 (default: all). Exit 0 iff all pass.
import fs from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';

const lit = process.argv[2];
const only = process.argv[3];
const ctl = process.env.LYS_CONTROL || '';   // red controls: see run_controls in run_all.sh
let failures = 0;
const say = (...a) => console.log(...a);
const sha = (buf) => crypto.createHash('sha256').update(buf).digest('hex');
const readJSON = (f) => JSON.parse(fs.readFileSync(path.join(lit, f), 'utf8'));

// ---------------------------------------------------------------- triples
function triples(file) {
  const J = readJSON(file);
  const n = J.n, [s0, s1] = J.sizes, R = J.rule, m = R.modulus;
  const reps = Array.isArray(R.colour0_orbit_representatives) ? R.colour0_orbit_representatives : JSON.parse(R.colour0_orbit_representatives);
  // group: x -> g x + b, g in <multipliers> mod m, points >= m fixed
  const gs = new Set([1]);
  let grew = true;
  while (grew) { grew = false; for (const g of [...gs]) for (const h of R.multipliers) { const k = (g * h) % m; if (!gs.has(k)) { gs.add(k); grew = true; } } }
  const maps = [];
  for (const g of gs) for (let b = 0; b < m; b++) maps.push((x) => (x >= m ? x : (g * x + b) % m));
  const key = (a, b, c) => { const t = [a, b, c].sort((x, y) => x - y); return (t[0] * n + t[1]) * n + t[2]; };
  const col = new Uint8Array(n * n * n).fill(1);
  for (const r of reps) for (const f of maps) col[key(f(r[0]), f(r[1]), f(r[2]))] = 0;
  // the stored string, read in lexicographic order
  const str = J.colouring;
  let i = 0, diff = 0;
  const C = new Uint8Array(n * n * n);
  for (let a = 0; a < n; a++) for (let b = a + 1; b < n; b++) for (let c = b + 1; c < n; c++) {
    const v = str.charCodeAt(i++) - 48;
    if (v !== col[key(a, b, c)]) diff++;
    for (const [x, y, z] of [[a, b, c], [a, c, b], [b, a, c], [b, c, a], [c, a, b], [c, b, a]]) C[(x * n + y) * n + z] = v;
  }
  if (ctl === 'plant' + n) { // red control: make the colour-0 set {0, 5, 9, 13} (or of size s0) monochromatic
    const P = [0, 5, 9, 13, 17, 21].slice(0, s0);
    for (let x = 0; x < P.length; x++) for (let y = 0; y < P.length; y++) for (let z = 0; z < P.length; z++)
      if (x !== y && y !== z && x !== z) C[(P[x] * n + P[y]) * n + P[z]] = 0;
  }
  const digestOK = sha(Buffer.from(str, 'ascii')) === J.sha256_of_colouring;
  say(`${file}: n=${n}, rebuilt from ${reps.length} representatives x ${maps.length} maps: ${diff} differences from the stored string; stored string SHA-256 as stated: ${digestOK}`);
  // invariance of the stored colouring under the generators
  let viol = 0;
  const gens = [(x) => (x >= m ? x : (x + 1) % m), ...R.multipliers.map((g) => (x) => (x >= m ? x : (g * x) % m))];
  for (let a = 0; a < n; a++) for (let b = a + 1; b < n; b++) for (let c = b + 1; c < n; c++) {
    const v = C[(a * n + b) * n + c];
    for (const f of gens) if (C[(f(a) * n + f(b)) * n + f(c)] !== v) viol++;
  }
  say(`  invariance under x+1 and the multipliers on every triple: ${viol} violations`);
  // forbidden sets through 0
  const found = [0, 0];
  for (const [colour, size] of [[0, s0], [1, s1]]) {
    const S = [0];
    const grow = (cands) => {
      if (S.length === size) { found[colour]++; return; }
      for (let t = 0; t < cands.length; t++) {
        const p = cands[t];
        // p joins S: every triple {S_i, S_j, p} must have this colour
        let ok = true;
        for (let x = 0; x < S.length && ok; x++) for (let y = x + 1; y < S.length && ok; y++) if (C[(S[x] * n + S[y]) * n + p] !== colour) ok = false;
        if (!ok) continue;
        S.push(p); grow(cands.slice(t + 1)); S.pop();
      }
    };
    const all = []; for (let p = 1; p < n; p++) all.push(p);
    grow(all);
  }
  say(`RESULT ${file}: colour 0 sets of size ${s0} through 0: ${found[0]}; colour 1 sets of size ${s1} through 0: ${found[1]}`);
  if (diff || !digestOK || viol || found[0] || found[1]) failures++;
}

// ---------------------------------------------------------------- 5-subsets of Z_83
function five(file) {
  const J = readJSON(file); const q = J.q, sh = J.short_form; const [s0, s1] = J.sizes;
  const inv = new Int32Array(q); for (let a = 1; a < q; a++) for (let b = 1; b < q; b++) if ((a * b) % q === 1) inv[a] = b;
  // index of {0,1,x,y,z}, 2 <= x < y < z <= 82, in lexicographic order
  const idx = new Map(); let k = 0;
  for (let x = 2; x < q; x++) for (let y = x + 1; y < q; y++) for (let z = y + 1; z < q; z++) idx.set((x * q + y) * q + z, k++);
  if (k !== sh.length) { say('short form length mismatch'); failures++; return; }
  const colourOf = (S) => { // S: five distinct points; normalise the two largest to 0, 1
    const T = [...S].sort((a, b) => b - a); const u = T[0], v = T[1];
    const d = inv[((v - u) % q + q) % q];
    const im = [T[2], T[3], T[4]].map((x) => ((((x - u) % q + q) % q) * d) % q).sort((a, b) => a - b);
    return sh.charCodeAt(idx.get((im[0] * q + im[1]) * q + im[2])) - 48;
  };
  // full colouring in lexicographic order, its digest
  const total = 29034396; const buf = Buffer.alloc(total); let i = 0;
  for (let a = 0; a < q; a++) for (let b = a + 1; b < q; b++) for (let c = b + 1; c < q; c++) for (let d = c + 1; d < q; d++) for (let e = d + 1; e < q; e++)
    buf[i++] = 48 + colourOf([a, b, c, d, e]);
  const digestOK = sha(buf) === J.sha256_of_full_colouring;
  say(`${file}: full colouring of ${i} 5-subsets built from the short form; SHA-256 equals the stated digest: ${digestOK}`);
  // invariance on every 5-subset, using the built string through a rank function of our own
  const bin = (nn, kk) => { let r = 1; for (let t = 0; t < kk; t++) r = (r * (nn - t)) / (t + 1); return Math.round(r); };
  const B = []; for (let nn = 0; nn <= q; nn++) { B.push([]); for (let kk = 0; kk <= 6; kk++) B[nn].push(nn < kk ? 0 : bin(nn, kk)); }
  const rank = (S) => { // colex-free lexicographic rank via combinatorial number system on complements
    const T = [...S].sort((a, b) => a - b); let r = 0, prev = -1;
    for (let j = 0; j < 5; j++) { for (let v = prev + 1; v < T[j]; v++) r += B[q - 1 - v][4 - j]; prev = T[j]; }
    return r;
  };
  let viol = 0; i = 0;
  const g1 = (x) => (x + 1) % q, g2 = (x) => (2 * x) % q;
  for (let a = 0; a < q; a++) for (let b = a + 1; b < q; b++) for (let c = b + 1; c < q; c++) for (let d = c + 1; d < q; d++) for (let e = d + 1; e < q; e++) {
    const v = buf[i++];
    if (((a * 7 + b * 3 + c + d + e) % 97) !== 0) continue; // invariance: a deterministic 1/97 sample here; full check below via the short form
    if (buf[rank([g1(a), g1(b), g1(c), g1(d), g1(e)])] !== v || buf[rank([g2(a), g2(b), g2(c), g2(d), g2(e)])] !== v) viol++;
  }
  // exact invariance argument: the colour of every 5-set was computed by normalising with an affine map, and the 20
  // normalisations agree iff the colouring is affine-invariant; check that on every 5-set through 0 and 1 (where the
  // short form itself is the colouring), for all 20 ordered pairs
  let incons = 0;
  for (let x = 2; x < q; x++) for (let y = x + 1; y < q; y++) for (let z = y + 1; z < q; z++) {
    const S = [0, 1, x, y, z]; const c0 = sh.charCodeAt(idx.get((x * q + y) * q + z)) - 48;
    for (let u = 0; u < 5; u++) for (let w = 0; w < 5; w++) if (u !== w) {
      const d = inv[((S[w] - S[u]) % q + q) % q];
      const im = S.filter((_, t) => t !== u && t !== w).map((p) => ((((p - S[u]) % q + q) % q) * d) % q).sort((a, b) => a - b);
      if (sh.charCodeAt(idx.get((im[0] * q + im[1]) * q + im[2])) - 48 !== c0) incons++;
    }
  }
  say(`  affine consistency on all 85320 sets through 0,1 (20 normalisations each): ${incons} disagreements; sampled invariance violations: ${viol}`);
  if (ctl === 'plant83') { // red control: force the 6-set {0,1,10,20,30,40} to be monochromatic in the lookup
    const P = [0, 1, 10, 20, 30, 40];
    for (let om = 0; om < 6; om++) { const S = P.filter((_, t) => t !== om); buf[rank(S)] = 48; }
  }
  // 6-sets through {0,1}
  const col = (S) => buf[rank(S)] - 48;
  const found = [0, 0]; let examined = 0;
  for (let a = 2; a < q; a++) for (let b = a + 1; b < q; b++) for (let c = b + 1; c < q; c++) {
    const cabc = col([0, 1, a, b, c]);
    for (let d = c + 1; d < q; d++) {
      examined++;
      if (col([0, 1, a, b, d]) !== cabc || col([0, 1, a, c, d]) !== cabc || col([0, 1, b, c, d]) !== cabc) continue;
      if (col([0, a, b, c, d]) !== cabc || col([1, a, b, c, d]) !== cabc) continue;
      found[cabc]++;
    }
  }
  say(`RESULT ${file}: 6-sets through {0,1} examined ${examined}; monochromatic: colour 0 ${found[0]}, colour 1 ${found[1]}`);
  if (!digestOK || incons || viol || found[0] || found[1] || s0 !== 6 || s1 !== 6) failures++;
}

// ---------------------------------------------------------------- 4-subsets of the projective line
function proj(file) {
  const J = readJSON(file); const q = J.q, INF = q, n = q + 1; const [s0, s1] = J.sizes; const T = J.table;
  const psl = J.rule.includes('delta');
  const md = (a) => ((a % q) + q) % q;
  const inv = new Int32Array(q); for (let a = 1; a < q; a++) { let x = 1; for (let e = 0; e < q - 2; e++) x = (x * a) % q; inv[a] = x; } // Fermat
  const sq = new Uint8Array(q); for (let a = 1; a < q; a++) sq[(a * a) % q] = 1;
  const tableAt = (x) => (ctl === 'flip' + n && x === 2 ? 1 - (T.charCodeAt(2) - 48) : T.charCodeAt(x) - 48);
  const colour = (p) => {
    const a = [...p].sort((x, y) => x - y);
    const x = a[3] === INF ? md(a[0] - a[2]) * inv[md(a[1] - a[2])] % q
      : md(a[0] - a[2]) * md(a[1] - a[3]) % q * inv[md(a[0] - a[3]) * md(a[1] - a[2]) % q] % q;
    if (!psl) return tableAt(x);
    const delta = md(a[1] - a[0]) * md(a[2] - a[0]) % q * md(a[2] - a[1]) % q;
    return sq[delta] ? tableAt(x) : tableAt(inv[x]);
  };
  let ok = true;
  if (!psl) {
    let bad = 0;
    for (let x = 2; x < q; x++) { // orbit of x under x -> 1/x and x -> 1 - x
      const orb = new Set([x]); let fr = [x];
      while (fr.length) { const nx = []; for (const y of fr) for (const z of [inv[y], md(1 - y)]) if (z > 1 && !orb.has(z)) { orb.add(z); nx.push(z); } fr = nx; }
      for (const y of orb) if (T[y] !== T[x]) bad++;
    }
    say(`${file}: table constant on the orbits of x -> 1/x, x -> 1-x: ${bad === 0} (${bad} differences)`);
    if (bad) ok = false;
  } else {
    let viol = 0; const buf = Buffer.alloc(J.subsets); let i = 0;
    const s = (x) => (x === INF ? INF : (x + 1) % q), t = (x) => (x === INF ? 0 : x === 0 ? INF : md(-inv[x]));
    for (let a = 0; a < n; a++) for (let b = a + 1; b < n; b++) for (let c = b + 1; c < n; c++) for (let d = c + 1; d < n; d++) {
      const v = colour([a, b, c, d]); buf[i++] = 48 + v;
      if (colour([s(a), s(b), s(c), s(d)]) !== v || colour([t(a), t(b), t(c), t(d)]) !== v) viol++;
    }
    const dOK = sha(buf) === J.colouring_sha256;
    say(`${file}: ${i} 4-subsets from the rule; SHA-256 of the colouring as stated: ${dOK}; invariance under x+1, -1/x: ${viol} violations`);
    if (viol || !dOK) ok = false;
    // transitivity on 3-subsets: union-find over the generators
    const id3 = (a, b, c) => { const u = [a, b, c].sort((x, y) => x - y); return (u[0] * n + u[1]) * n + u[2]; };
    const par = new Int32Array(n * n * n).map((_, j) => j);
    const find = (x) => { while (par[x] !== x) { par[x] = par[par[x]]; x = par[x]; } return x; };
    for (let a = 0; a < n; a++) for (let b = a + 1; b < n; b++) for (let c = b + 1; c < n; c++)
      for (const g of [s, t]) { const r1 = find(id3(a, b, c)), r2 = find(id3(g(a), g(b), g(c))); if (r1 !== r2) par[r1] = r2; }
    const root = find(id3(0, 1, 2)); let outside = 0;
    for (let a = 0; a < n; a++) for (let b = a + 1; b < n; b++) for (let c = b + 1; c < n; c++) if (find(id3(a, b, c)) !== root) outside++;
    say(`  3-subsets not in the class of {0,1,2}: ${outside}`);
    if (outside) ok = false;
  }
  const found = [0, 0];
  for (const [colourWanted, size] of [[0, s0], [1, s1]]) {
    const S = [0, 1, 2];
    const others = []; for (let p = 3; p < n; p++) others.push(p);
    const grow = (from) => {
      if (S.length === size) { found[colourWanted]++; return; }
      for (let t = from; t < others.length; t++) {
        const p = others[t]; let good = true;
        for (let x = 0; x < S.length && good; x++) for (let y = x + 1; y < S.length && good; y++) for (let z = y + 1; z < S.length && good; z++)
          if (colour([S[x], S[y], S[z], p]) !== colourWanted) good = false;
        if (!good) continue;
        S.push(p); grow(t + 1); S.pop();
      }
    };
    grow(0);
  }
  say(`RESULT ${file}: sets through {0,1,2}: colour 0 of size ${s0}: ${found[0]}; colour 1 of size ${s1}: ${found[1]}`);
  if (!ok || found[0] || found[1]) failures++;
}

const cases = {
  t37: () => triples('R453_37_points.json'), t70: () => triples('R463_70_points.json'), t98: () => triples('R553_98_points.json'),
  f83: () => five('R665_83_points_short.json'), p200: () => proj('R564_200_points_rule.json'), p600: () => proj('R664_600_points_rule.json'),
};
for (const [k, f] of Object.entries(cases)) if (!only || only === k) { const t0 = Date.now(); f(); say(`  (${k}: ${((Date.now() - t0) / 1000).toFixed(1)} s)`); }
say(failures ? `FAIL (${failures})` : 'ALL PASS');
process.exit(failures ? 1 : 0);
