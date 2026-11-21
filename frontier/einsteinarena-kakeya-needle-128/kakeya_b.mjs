// Program B (JavaScript BigInt), written separately from program A: the area of the union of the 128 triangles by
// Green's theorem, area = sum over the union's boundary of the integral of x dy (counter-clockwise). Each triangle
// T_j = conv{(x_j,0), (x_j+1/128,0), (x_j+j/128,1)} is traversed counter-clockwise: the base (dy = 0, contributes 0),
// the right side upwards, the left side downwards. A point of a side lies on the union's boundary iff it is in the
// interior of no other triangle; for each side and each other triangle that set of heights is an open interval
// (two affine inequalities), so the boundary part of a side is [0,1] minus a finite union of open intervals, and
// the integral of the affine x(y) over it is exact. Coinciding sides of two triangles (a left side of T_j on a right
// side of T_(j+1)) have opposite directions and cancel, as they should (the shared segment is inside the union).
// Rationals are pairs of BigInts. No floating point decides anything.
// Usage: node kakeya_b.mjs best_26.json SOLUTION_ID [control]   controls: shift (one triangle moved), translate
import fs from 'node:fs';

const [file, sid, control = ''] = process.argv.slice(2);
const sols = JSON.parse(fs.readFileSync(file, 'utf8'));
const sol = sols.find((s) => String(s.id) === sid);
const data = typeof sol.data === 'string' ? JSON.parse(sol.data) : sol.data;

const gcd = (a, b) => { a = a < 0n ? -a : a; b = b < 0n ? -b : b; while (b) [a, b] = [b, a % b]; return a; };
const R = (n, d = 1n) => { if (d < 0n) { n = -n; d = -d; } const g = gcd(n, d) || 1n; return [n / g, d / g]; };
const add = (x, y) => R(x[0] * y[1] + y[0] * x[1], x[1] * y[1]);
const sub = (x, y) => R(x[0] * y[1] - y[0] * x[1], x[1] * y[1]);
const mul = (x, y) => R(x[0] * y[0], x[1] * y[1]);
const div = (x, y) => R(x[0] * y[1], x[1] * y[0]);
const cmp = (x, y) => { const d = x[0] * y[1] - y[0] * x[1]; return d < 0n ? -1 : d > 0n ? 1 : 0; };
const ZERO = R(0n), ONE = R(1n), HALF = R(1n, 2n);
const dec = (s) => { // exact decimal string -> rational
  let t = s.trim(); let neg = false; if (t[0] === '-' || t[0] === '+') { neg = t[0] === '-'; t = t.slice(1); }
  const [ip, fp = ''] = t.split('.'); if (!/^\d+$/.test(ip + fp) || fp.length > 25) throw new Error('bad ' + s);
  const v = R(BigInt(ip + fp), 10n ** BigInt(fp.length)); return neg ? R(-v[0], v[1]) : v;
};
let xs = data.offsets.map(dec);
if (xs.length !== 128) throw new Error('need 128 offsets');
if (control === 'shift') xs[63] = add(xs[63], R(1n, 1000n));
if (control === 'translate') xs = xs.map((x) => add(x, R(7n, 3n)));
const N = 128n;
// sides as affine functions of y: x = c0 + c1 y
const left = xs.map((x, i) => [x, R(BigInt(i + 1), N)]);
const right = xs.map((x, i) => [add(x, R(1n, N)), R(BigInt(i), N)]);
// heights where f(y) < g(y) (f, g affine): an open half-line or all/nothing; returned as [lo, hi] open, clipped later
const lessSet = (f, g) => { // {y : f(y) < g(y)} = {y : (g0 - f0) + (g1 - f1) y > 0}
  const a = sub(g[0], f[0]), b = sub(g[1], f[1]);
  if (b[0] === 0n) return cmp(a, ZERO) > 0 ? [null, null] : 'empty';
  const root = div(R(-a[0], a[1]), b); // y where equal
  return b[0] > 0n ? [root, null] : [null, root];
};
const inter = (I, J) => {
  if (I === 'empty' || J === 'empty') return 'empty';
  const lo = I[0] === null ? J[0] : J[0] === null ? I[0] : (cmp(I[0], J[0]) > 0 ? I[0] : J[0]);
  const hi = I[1] === null ? J[1] : J[1] === null ? I[1] : (cmp(I[1], J[1]) < 0 ? I[1] : J[1]);
  return [lo, hi];
};
// integral of c0 + c1 y over [a, b]
const integ = (f, a, b) => add(mul(f[0], sub(b, a)), mul(mul(f[1], HALF), sub(mul(b, b), mul(a, a))));
let area = ZERO; let pieces = 0;
const sides = [];
for (let j = 0; j < 128; j++) { sides.push({ f: right[j], sign: 1, own: j }); sides.push({ f: left[j], sign: -1, own: j }); }
for (const { f, sign, own } of sides) {
  // covered heights: union over k != own of {0 < y < 1 : left_k(y) < f(y) < right_k(y)}
  const cov = [];
  for (let k = 0; k < 128; k++) {
    if (k === own) continue;
    let I = inter(lessSet(left[k], f), lessSet(f, right[k]));
    if (I === 'empty') continue;
    let lo = I[0] === null ? ZERO : I[0], hi = I[1] === null ? ONE : I[1];
    if (cmp(lo, ZERO) < 0) lo = ZERO; if (cmp(hi, ONE) > 0) hi = ONE;
    if (cmp(lo, hi) < 0) cov.push([lo, hi]);
  }
  cov.sort((p, q) => cmp(p[0], q[0]));
  // complement of the union of open intervals in [0,1]
  let pos = ZERO;
  for (const [lo, hi] of cov) {
    if (cmp(lo, pos) > 0) { area = sign > 0 ? add(area, integ(f, pos, lo)) : sub(area, integ(f, pos, lo)); pieces++; }
    if (cmp(hi, pos) > 0) pos = hi;
  }
  if (cmp(pos, ONE) < 0) { area = sign > 0 ? add(area, integ(f, pos, ONE)) : sub(area, integ(f, pos, ONE)); pieces++; }
}
console.log(`solution ${sol.id} by ${sol.agentName}, stated score ${sol.score}; control: ${control || 'none'}`);
console.log(`boundary pieces integrated: ${pieces}`);
console.log(`exact area = ${area[0]}/${area[1]}`);
const digits = (area[0] * 10n ** 30n / area[1]).toString().padStart(30, '0');
console.log(`area = 0.${digits}... (30 digits, truncated)`);
console.log(`Number(area) as the nearest double: ${Number(area[0] * 10n ** 25n / area[1]) / 1e25}`);
