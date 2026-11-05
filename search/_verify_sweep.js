/* 独立验证器 v2: 只用**最终 SVG 的真实几何**核对"悬挂环系链两侧夹角 = 90°+180°/n"(苯环 120°/120°)。
 *
 * v1 的缺陷(2026-11-05): 它按"键轴方向夹角"直接相减, 而键轴方向由"路径端点中最远两点连线"定,
 *   正负号是任意的 —— 同一个顶点两条键可能一正一反, 于是 120° 被量成 60°, 出现 a1+a2+inner=240
 *   的自相矛盾假失败(通过率被压到 0.4%)。本版改为:
 *     1) 每条键取其全部路径点中最远两点连线, 得到一条**无限长直线**;
 *     2) 相邻两条键的直线求交 → 得到**顶点真实坐标**;
 *     3) 用顶点向量算内角/系链角(向量有向, 无符号歧义), 系链方向用"背离环心"定向。
 *   这样量的是"图上看得见的东西", 不依赖任何内部预测。
 * 另核验: 每条块内键在最终 SVG 有对应路径(图元守恒)、边长一致性、坐标无 NaN、不越出画布。 */

const fs = require('fs'), path = require('path');
const DIR = __dirname;

globalThis.window = globalThis;
globalThis.addEventListener = () => {};
window.addEventListener = () => {};
window.document = { getElementById: () => null, createElement: () => ({ style: {} }), addEventListener: () => {} };
let src = fs.readFileSync(path.join(DIR, 'app.js'), 'utf8').replace(/'use strict';/, '');
src += '\n;globalThis.__setRDKit=(m)=>{RDKit=m;};';
(0, eval)(src);
if (typeof molToSVG !== 'function') { console.error('FATAL: molToSVG not loaded'); process.exit(1); }
/* 回归对照: NOSTEREO=1 时把"清双键 either 标记"关掉, 用于证明本轮修复没有影响环角度 */
if (process.env.NOSTEREO) { _prpClearBondStereo = (b) => b; console.log('[TEST] stereo clearing DISABLED'); }

/* ---------- 独立解析 molblock ---------- */
function parseBlk(blk) {
  const L = blk.split('\n');
  const c = L[3].trim().split(/\s+/).map(Number);
  const nA = c[0], nB = c[1];
  const atoms = [], adj = Array.from({ length: nA }, () => []);
  for (let i = 0; i < nA; i++) {
    const t = L[4 + i].trim().split(/\s+/);
    atoms.push({ x: +t[0], y: +t[1], el: t[3] });
  }
  const bonds = [];
  for (let i = 0; i < nB; i++) {
    const m = L[4 + nA + i].trim().split(/\s+/).map(Number);
    bonds.push([m[0] - 1, m[1] - 1]);
    adj[m[0] - 1].push(m[1] - 1); adj[m[1] - 1].push(m[0] - 1);
  }
  return { nAtoms: nA, atoms, bonds, adj };
}

/* ---------- 变换感知的最终 SVG 路径点 ---------- */
function parseTf(tf) {
  const ops = []; const re = /(translate|scale|rotate)\s*\(([^)]*)\)/g; let m;
  while ((m = re.exec(tf))) ops.push({ t: m[1], a: m[2].split(/[\s,]+/).map(Number) });
  return ops;
}
function applyOps(ops, p) {
  for (let i = ops.length - 1; i >= 0; i--) {
    const { t, a } = ops[i];
    if (t === 'translate') p = [p[0] + a[0], p[1] + (a[1] || 0)];
    else if (t === 'scale') p = [p[0] * a[0], p[1] * (a[1] !== undefined ? a[1] : a[0])];
    else if (t === 'rotate') {
      const r = a[0] * Math.PI / 180, ca = Math.cos(r), sa = Math.sin(r), cx = a[1] || 0, cy = a[2] || 0;
      const x = p[0] - cx, y = p[1] - cy;
      p = [cx + x * ca - y * sa, cy + x * sa + y * ca];
    }
  }
  return p;
}
function applyStack(st, p) { for (let i = st.length - 1; i >= 0; i--) p = applyOps(st[i], p); return p; }
/* 引号无关取属性 —— app.js 插入的 <g transform="..."> 是双引号, RDKit 原生路径是单引号, 两种都要认,
 * 否则工序二的 rotate 组与整体缩放组被静默忽略, "已被旋转修好"的环被量成旋转前的旧角度(假失败)。 */
function qattr(attrs, name) {
  let m = attrs.match(new RegExp(name + "='([^']*)'")); if (m) return m[1];
  m = attrs.match(new RegExp(name + '="([^"]*)"')); if (m) return m[1];
  return null;
}
function finalPaths(svg) {
  const stack = [], res = [];
  const re = /<(\/?)(path|g|text)\b([^>]*)>/g; let m;
  while ((m = re.exec(svg))) {
    const close = m[1] === '/', tag = m[2], attrs = m[3];
    if (tag === 'g') {
      if (close) stack.pop();
      else { const tf = qattr(attrs, 'transform'); if (tf) stack.push(parseTf(tf)); }
      continue;
    }
    if (tag === 'path') {
      const cls = qattr(attrs, 'class') || '';
      const d = qattr(attrs, 'd') || '';
      const pts = [], subs = [];
      for (const chunk of d.split(/[Mm]/).filter(s => s.trim())) {
        const sp = [];
        const reN = /(-?\d+\.?\d*)[,\s]+(-?\d+\.?\d*)/g; let q;
        while ((q = reN.exec(chunk))) sp.push(applyStack(stack, [+q[1], +q[2]]));
        if (sp.length) { subs.push(sp); for (const p of sp) pts.push(p); }
      }
      res.push({ cls, pts, subs });
    }
  }
  return res;
}

/* ---------- 键 → 无限长直线 ---------- */
function atomsInCls(cls) { const s = new Set(); let m; const re = /atom-(\d+)/g; while ((m = re.exec(cls))) s.add(+m[1]); return s; }
const _lineCache = new WeakMap();
function bondLine(fps, a, b) {
  let M = _lineCache.get(fps); if (!M) { M = new Map(); _lineCache.set(fps, M); }
  const key = Math.min(a, b) + '-' + Math.max(a, b);
  if (M.has(key)) return M.get(key);
  const pts = [], subLines = [];
  for (const p of fps) {
    if (p.cls.indexOf('bond-') < 0) continue;
    const at = atomsInCls(p.cls);
    if (!(at.has(a) && at.has(b))) continue;
    for (const q of p.pts) pts.push(q);
    for (const sp of (p.subs || [])) subLines.push(sp);
  }
  /* 关键(2026-11-05): 直线取"**最长子折线**"的首尾, 而不是"全部点里最远的一对"。
   * 双键/芳香键有内侧平行线(侧向偏移约 2.4px), 两条平行线四个端点里, **对角线**比对边更长
   * (√(L²+2.4²) > L), 于是"最远一对"会挑到对角线, 方向被扭 atan(2.4/L) ≈ 6°;
   * 两条这样的线再求交, 顶点被推偏好几个 px, 正六边形就被量成 109°/131° 的假相(v1 第二轮假失败)。
   * 最长子折线一定是主线, 首尾都在键轴上, 方向与键轴严格一致。 */
  let out = null;
  if (pts.length >= 2) {
    let best = null;
    for (const sp of subLines) {
      if (sp.length < 2) continue;
      const A = sp[0], B = sp[sp.length - 1];
      const d = Math.hypot(B[0] - A[0], B[1] - A[1]);
      if (!best || d > best.d) best = { d, A, B };
    }
    if (best && best.d > 1e-6) {
      const ux = (best.B[0] - best.A[0]) / best.d, uy = (best.B[1] - best.A[1]) / best.d;
      out = { p: [(best.A[0] + best.B[0]) / 2, (best.A[1] + best.B[1]) / 2], u: [ux, uy], len: best.d };
    }
  }
  M.set(key, out);
  return out;
}
function inter(l1, l2) {
  if (!l1 || !l2) return null;
  const cr = l1.u[0] * l2.u[1] - l1.u[1] * l2.u[0];
  if (Math.abs(cr) < 0.08) return null;              // 近平行, 无稳定交点
  const dx = l2.p[0] - l1.p[0], dy = l2.p[1] - l1.p[1];
  const t = (dx * l2.u[1] - dy * l2.u[0]) / cr;
  return [l1.p[0] + t * l1.u[0], l1.p[1] + t * l1.u[1]];
}
function sub(p, q) { return [p[0] - q[0], p[1] - q[1]]; }
function angDeg(u, v) {
  const nu = Math.hypot(u[0], u[1]), nv = Math.hypot(v[0], v[1]);
  if (nu < 1e-9 || nv < 1e-9) return null;
  let c = (u[0] * v[0] + u[1] * v[1]) / (nu * nv);
  c = Math.max(-1, Math.min(1, c));
  return Math.acos(c) * 180 / Math.PI;
}

/* ---------- 悬挂环检测(与 app.js 同门槛) ---------- */
function findRings(adj, n, t) {
  const rings = [];
  for (let s = 0; s < n; s++) {
    const ipath = [s], used = new Set([s]);
    (function dfs(v) {
      if (ipath.length === t) { if (adj[v].includes(s)) rings.push([...ipath]); return; }
      for (const nb of adj[v]) { if (nb === s || used.has(nb)) continue; used.add(nb); ipath.push(nb); dfs(nb); ipath.pop(); used.delete(nb); }
    })(s);
  }
  const u = [], seen = new Set();
  for (const r of rings) { const k = [...r].sort((x, y) => x - y).join(','); if (!seen.has(k)) { seen.add(k); u.push(r); } }
  return u;
}
function detectJobs(P) {
  const rings = [];
  for (const t of [3, 4, 5, 6, 7, 8]) for (const r of findRings(P.adj, P.nAtoms, t)) rings.push(r);
  const atomRings = Array.from({ length: P.nAtoms }, () => []);
  rings.forEach((r, i) => { for (const a of r) atomRings[a].push(i); });
  const isShared = (ri) => rings[ri].some(a => atomRings[a].length > 1);
  const reachExcl = (nb, ca) => {
    const seen = new Set([nb]); const st = [nb]; let cnt = 0;
    while (st.length) { const u = st.pop(); cnt++; for (const v of P.adj[u]) { if (seen.has(v) || ca.has(v)) continue; seen.add(v); st.push(v); } }
    return cnt;
  };
  const jobs = [];
  for (let ri = 0; ri < rings.length; ri++) {
    if (isShared(ri)) continue;
    const ringAtoms = new Set(rings[ri]);
    let tether = null, best = -1;
    for (const a of rings[ri]) for (const nb of P.adj[a]) { if (ringAtoms.has(nb)) continue; const sz = reachExcl(nb, ringAtoms); if (sz > best) { best = sz; tether = { ip: a, anchor: nb }; } }
    if (!tether) continue;
    const frag = new Set(ringAtoms); const q = [...ringAtoms];
    while (q.length) { const u = q.pop(); for (const nb of P.adj[u]) { if (frag.has(nb)) continue; if ((u === tether.ip && nb === tether.anchor) || (u === tether.anchor && nb === tether.ip)) continue; frag.add(nb); q.push(nb); } }
    if (frag.size > Math.max(8, Math.floor(P.nAtoms * 0.6))) continue;
    const ci = rings[ri].indexOf(tether.ip);
    const cyc = rings[ri].slice(ci).concat(rings[ri].slice(0, ci));
    jobs.push({ ringAtoms, tether, frag, cyc });
  }
  return jobs;
}

/* ---------- 单分子核验 ---------- */
function verifyOne(sm, label) {
  const svg = molToSVG(sm);
  if (!svg) return { err: 'molToSVG null' };
  if (/NaN|undefined/.test(svg)) return { err: 'NaN/undefined in svg' };
  let blk = null;
  { const mol = RDKit.get_mol(sm); if (!mol) return { err: 'get_mol' }; try { blk = mol.get_new_coords(); } catch (e) { } mol.delete(); }
  if (!blk) return { err: 'no blk' };
  const P = parseBlk(blk);
  if (P.nAtoms > 110) return { skip: 'too big ' + P.nAtoms };
  const fps = finalPaths(svg);
  const rotC = []; { const reR = /rotate\(\s*(-?[\d.]+)\s*,\s*(-?[\d.]+)\s*,\s*(-?[\d.]+)\s*\)/g; let mR; while ((mR = reR.exec(svg))) rotC.push([+mR[2], +mR[3]]); }
  const jobs = detectJobs(P);
  const results = [];
  for (const j of jobs) {
    const n = j.cyc.length, target = 90 + 180 / n, tInner = 180 - 360 / n;
    const { ip, anchor } = j.tether;
    const V = [];
    let bad = false;
    for (let k = 0; k < n; k++) {
      const L1 = bondLine(fps, j.cyc[k], j.cyc[(k + 1) % n]);
      const L2 = bondLine(fps, j.cyc[(k - 1 + n) % n], j.cyc[k]);
      const v = inter(L1, L2);
      if (!v) { bad = true; break; }
      V.push(v);
    }
    const Lt = bondLine(fps, ip, anchor);
    if (bad || !Lt) { results.push({ n, fail: 'no vertex/line', ip, anchor }); continue; }
    const ipP = V[0];
    const C = [V.reduce((s, p) => s + p[0], 0) / n, V.reduce((s, p) => s + p[1], 0) / n];
    let td = Lt.u;
    if ((ipP[0] - C[0]) * td[0] + (ipP[1] - C[1]) * td[1] < 0) td = [-td[0], -td[1]];
    const a1 = angDeg(td, sub(V[1], ipP));
    const a2 = angDeg(td, sub(V[n - 1], ipP));
    const inner = angDeg(sub(V[1], ipP), sub(V[n - 1], ipP));
    const innAll = [];
    for (let k = 0; k < n; k++) innAll.push(angDeg(sub(V[(k + 1) % n], V[k]), sub(V[(k - 1 + n) % n], V[k])));
    const lens = [];
    for (let k = 0; k < n; k++) lens.push(Math.hypot(V[(k + 1) % n][0] - V[k][0], V[(k + 1) % n][1] - V[k][1]));
    const allC = j.cyc.every(a => P.atoms[a].el === 'C');
    results.push({
      n, ip, anchor, tetherEl: P.atoms[anchor].el + P.atoms[ip].el,
      a1: +a1.toFixed(1), a2: +a2.toFixed(1), inner: +inner.toFixed(1),
      dev1: +(a1 - target).toFixed(1), dev2: +(a2 - target).toFixed(1),
      devInner: +(inner - tInner).toFixed(1),
      innMin: +Math.min(...innAll).toFixed(1), innMax: +Math.max(...innAll).toFixed(1),
      lenRatio: allC ? +(Math.max(...lens) / Math.min(...lens)).toFixed(2) : null,
      rotated: rotC.some(c => Math.hypot(c[0] - ipP[0], c[1] - ipP[1]) < 2),
      pass: Math.abs(a1 - target) <= 2.5 && Math.abs(a2 - target) <= 2.5 && Math.abs(inner - tInner) <= 3
    });
  }
  let missing = 0;
  for (const [a, b] of P.bonds) if (!bondLine(fps, a, b)) missing++;
  let oflow = 0;
  for (const p of fps) for (const q of p.pts) if (q[0] < -12 || q[0] > 262 || q[1] < -12 || q[1] > 207 || !isFinite(q[0]) || !isFinite(q[1])) oflow++;
  return { results, missing, oflow, nAtoms: P.nAtoms, label };
}

const initRDKitMod = require(path.join(DIR, 'lib/rdkit/RDKit_minimal.js')).default;
initRDKitMod({ locateFile: f => path.join(DIR, 'lib/rdkit', f) }).then((mod) => {
  globalThis.__setRDKit(mod);
  globalThis.RDKit = mod;
  mod.prefer_coordgen(true);
  sweep();
}).catch(e => { console.error('RDKit INIT FAIL', e); process.exit(1); });


/* 全库扫描 main: 逐分片抽样(省内存) + 联芳基类定向用例 */
function sweep() {
  const shardIds = fs.readdirSync(path.join(DIR, 'data')).filter(f => /^products\.\d+\.json$/.test(f)).sort((a, b) => +a.match(/\d+/)[0] - +b.match(/\d+/)[0]);
  const stats = { sweep: { rings: 0, pass: 0, missing: 0, oflow: 0, worst: [], failMiss: 0, nested: 0 } };
  const errs = [];
  let nDone = 0, t0 = Date.now();
  const PER_SHARD = 60;
  const patterns = [
    ['联芳醚 PhO', /ccccc1O[a-z]/, 80],
    ['联芳硫 PhS', /ccccc1S[a-z]/, 80],
    ['苄基 Bn', /O[Cc]c1ccccc1|CSc1ccccc1|NCCc1ccccc1/, 80],
  ];
  const patLeft = patterns.map(p => p[2]);
  const patSeen = patterns.map(() => 0);
  for (const f of shardIds) {
    let d = null;
    try { d = JSON.parse(fs.readFileSync(path.join(DIR, 'data', f), 'utf8')); } catch (e) { continue; }
    const step = Math.max(1, Math.floor(d.length / PER_SHARD));
    const picked = [];
    for (let i = 0; i < d.length && picked.length < PER_SHARD; i += step) {
      const sm = d[i].s || '';
      if (sm.length >= 20 && sm.length <= 200 && !sm.includes('.')) picked.push({ sm, name: d[i].n, tag: 'sweep' });
    }
    // 定向联芳基类
    for (let i = 0; i < d.length; i++) {
      const sm = d[i].s || '';
      for (let pi = 0; pi < patterns.length; pi++) {
        if (patLeft[pi] <= 0) continue;
        if (patterns[pi][1].test(sm)) {
          picked.push({ sm, name: d[i].n, tag: patterns[pi][0] });
          patLeft[pi]--; patSeen[pi]++; break;
        }
      }
      if (patLeft.every(x => x <= 0)) break;
    }
    for (const rec of picked) {
      let v;
      try { v = verifyOne(rec.sm, rec.tag + ' | ' + (rec.name || '').slice(0, 40)); }
      catch (e) { errs.push({ label: rec.tag, e: e.message }); continue; }
      nDone++;
      if (v.err) { errs.push({ label: rec.tag + ' | ' + (rec.name || ''), e: v.err, sm: rec.sm }); continue; }
      if (v.skip) continue;
      const s = stats[rec.tag] || (stats[rec.tag] = { rings: 0, pass: 0, missing: 0, oflow: 0, worst: [], failMiss: 0 });
      s.missing += v.missing; s.oflow += v.oflow;
      for (const r of v.results) {
        s.rings++;
        if (r.fail) { s.failMiss++; errs.push({ label: v.label, e: 'ring fail ' + JSON.stringify(r) }); continue; }
        if (r.pass) s.pass++;
        else if (s.worst.length < 8) s.worst.push({ label: v.label, r });
      }
    }
    d = null; picked.length = 0;
  }
  const secs = ((Date.now() - t0) / 1000).toFixed(1);
  console.log('\n===== 全库扫描 (' + nDone + ' 分子, ' + secs + 's, ' + (nDone ? (secs * 1000 / nDone).toFixed(0) : '-') + 'ms/分子含测量) =====');
  let totR = 0, totP = 0;
  for (const [tag, s] of Object.entries(stats)) {
    totR += s.rings; totP += s.pass;
    console.log('[' + tag + '] 环=' + s.rings + ' 通过=' + s.pass + (s.rings ? ' (' + (100 * s.pass / s.rings).toFixed(2) + '%)' : '')
      + ' | 缺失键=' + s.missing + ' 溢出点=' + s.oflow + ' 丢路径环=' + (s.failMiss || 0));
    for (const w of s.worst) console.log('  FAIL:', JSON.stringify(w.r), '|', w.label.slice(0, 65));
  }
  console.log('\n总计: 环=' + totR + ' 通过=' + totP + (totR ? ' (' + (100 * totP / totR).toFixed(2) + '%)' : ''));
  console.log('定向用例命中: PhO=' + patSeen[0] + ' PhS=' + patSeen[1] + ' Bn=' + patSeen[2]);
  console.log('错误分子数=' + errs.length);
  for (const e of errs.slice(0, 8)) console.log('  ERR:', e.e, '|', (e.label || '').slice(0, 55), '|', (e.sm || '').slice(0, 70));
  process.exit(errs.length || totP < totR ? 1 : 0);
}
