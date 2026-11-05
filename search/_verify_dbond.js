/* 验证"C=N / C=C 被画成交叉双键"是否修好。
 * 判据(纯几何, 只看最终 SVG): 对一个双键, 把它所有 bond 路径各看成一条直线(取首尾点),
 *   正常双键 = 两条**互相平行**的直线(与键轴夹角≈0);
 *   交叉双键 = 两条直线**与键轴成 ±10~30°**并在中点相交 → 直线与键轴夹角明显非 0。
 * 故指标 = maxDev("各直线方向与键轴方向的夹角")。正常应 < 3°, 交叉约 15~25°。
 * 同时统计 molblock 里 stereo=3/4 的双键个数(根因指标)。 */
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

/* ---- 解析 molblock 键块: 类型 + stereo ---- */
function parseBonds(blk) {
  const L = blk.split('\n');
  const c = (L[3] || '').trim().split(/\s+/).map(Number);
  const nA = c[0] | 0, nB = c[1] | 0;
  const b = [];
  for (let i = 0; i < nB; i++) {
    const s = L[4 + nA + i] || '';
    const a1 = parseInt(s.slice(0, 3), 10), a2 = parseInt(s.slice(3, 6), 10);
    const t = parseInt(s.slice(6, 9), 10), st = parseInt(s.slice(9, 12), 10);
    b.push({ a: a1 - 1, b: a2 - 1, t, st });
  }
  return b;
}

/* ---- 从 SVG 取每个 bond 的全部直线(忽略 <g> 变换: 组内旋转对各直线同样施加,
        本指标只比较**同一键内**直线之间的夹角, 变换无关) ---- */
function bondLines(svg) {
  const map = new Map();
  const re = /<path\b([^>]*)\/?>/g; let m;
  while ((m = re.exec(svg))) {
    const tag = m[1];
    const cls = (tag.match(/class='([^']*)'/) || [])[1] || '';
    if (cls.indexOf('bond-') < 0) continue;
    const ix = []; let q; const reA = /atom-(\d+)/g;
    while ((q = reA.exec(cls))) ix.push(+q[1]);
    if (ix.length !== 2) continue;
    const key = Math.min(ix[0], ix[1]) + '-' + Math.max(ix[0], ix[1]);
    const d = (tag.match(/d='([^']*)'/) || [])[1] || '';
    for (const chunk of d.split(/[Mm]/).filter(s => s.trim())) {
      const pts = []; const reN = /(-?\d+\.?\d*)[,\s]+(-?\d+\.?\d*)/g; let r;
      while ((r = reN.exec(chunk))) pts.push([+r[1], +r[2]]);
      if (pts.length < 2) continue;
      const A = pts[0], B = pts[pts.length - 1];
      const len = Math.hypot(B[0] - A[0], B[1] - A[1]);
      if (len < 2.5) continue;                       // 太短的残段不参与
      if (!map.has(key)) map.set(key, []);
      map.get(key).push({ A, B, len, ang: Math.atan2(B[1] - A[1], B[0] - A[0]) });
    }
  }
  return map;
}
const modPi = x => { let v = x % Math.PI; if (v < 0) v += Math.PI; return v; };
const devPi = (x, y) => { const d = Math.abs(modPi(x) - modPi(y)); return Math.min(d, Math.PI - d) * 180 / Math.PI; };

/* 某双键的 maxDev: 键轴 = 每条直线各自方向的最长直线方向;
   但"交叉"时两条线都不代表轴, 故取"使 maxDev 最小"的那条直线方向作为轴估计。 */
function maxDevOfLines(lines) {
  if (!lines || lines.length < 2) return null;
  let best = 1e9;
  for (const cand of lines) {
    let mx = 0;
    for (const l of lines) mx = Math.max(mx, devPi(l.ang, cand.ang));
    if (mx < best) best = mx;
  }
  return best;
}
/* 双键是否被画成交叉: 两条直线夹角 > 4°(正常平行应≈0) */
function linesSpread(lines) {
  if (!lines || lines.length < 2) return 0;
  let mx = 0;
  for (let i = 0; i < lines.length; i++) for (let j = i + 1; j < lines.length; j++)
    mx = Math.max(mx, devPi(lines[i].ang, lines[j].ang));
  return mx;
}

/* ---- 旧管线(不清 stereo), 用于给出"修复前"基线 ---- */
function renderOld(sm) {
  const mol = RDKit.get_mol(sm); if (!mol) return null;
  let out = null, blk = null;
  try { blk = mol.get_new_coords(); } catch (e) { }
  let m2 = null;
  try { m2 = blk ? RDKit.get_mol(blk) : null; } catch (e) { m2 = null; }
  try { out = (m2 || mol).get_svg(250, 195); } catch (e) { out = null; }
  if (m2) try { m2.delete(); } catch (e) { }
  mol.delete();
  return { svg: out, blk };
}

function audit(sm) {
  const old = renderOld(sm);
  const svgNew = molToSVG(sm);
  if (!old || !old.svg || !svgNew) return { err: 'render' };
  const bonds = parseBonds(old.blk);
  const dbl = bonds.filter(b => b.t === 2);
  const either = dbl.filter(b => b.st === 3 || b.st === 4);
  const Lold = bondLines(old.svg), Lnew = bondLines(svgNew);
  let xOld = 0, xNew = 0, worstOld = 0, worstNew = 0, overOld = [], overNew = [];
  for (const b of dbl) {
    const k = Math.min(b.a, b.b) + '-' + Math.max(b.a, b.b);
    const so = linesSpread(Lold.get(k)), sn = linesSpread(Lnew.get(k));
    if (so > 4) { xOld++; overOld.push({ k, so: +so.toFixed(1) }); }
    if (sn > 4) { xNew++; overNew.push({ k, sn: +sn.toFixed(1) }); }
    worstOld = Math.max(worstOld, so); worstNew = Math.max(worstNew, sn);
  }
  return { ndbl: dbl.length, either: either.length, xOld, xNew, worstOld: +worstOld.toFixed(1), worstNew: +worstNew.toFixed(1), overOld, overNew };
}

const CASES = [
  ['肟 C=N (DPPH醛肟酯)', 'O=P(ON=Cc1ccccc1)(c1ccccc1)c1ccccc1'],
  ['肟 C=N (对氟苯甲醛)', 'O=P(ON=Cc1ccc(F)cc1)(c1ccccc1)c1ccccc1'],
  ['腙 C=N', 'O=C(NN=Cc1ccc(=O)[nH]c1)c1ccc(-c2ccccc2)o1'],
  ['席夫碱 C=N', 'COc1ncncc1C=Nc1ccc(F)c(F)c1'],
  ['查尔酮 C=C', 'O=C(/C=C/c1ccccc1)c1ccc(N(C)C)cc1'],
  ['查尔酮 C=C (无E/Z标记)', 'O=C(C=Cc1ccccc1)c1ccc(N(C)C)cc1'],
  ['腈 C#N (对照, 非双键)', 'N#Cc1ccc(-c2ccc(C#N)cc2)cc1'],
  ['联苯 (全芳香对照)', 'c1ccc(-c2ccccc2)cc1'],
  ['阿司匹林', 'CC(=O)Oc1ccccc1C(=O)O'],
];

const initRDKitMod = require(path.join(DIR, 'lib/rdkit/RDKit_minimal.js')).default;
initRDKitMod({ locateFile: f => path.join(DIR, 'lib/rdkit', f) }).then(mod => {
  globalThis.__setRDKit(mod);
  globalThis.RDKit = mod;
  mod.prefer_coordgen(true);
  main();
}).catch(e => { console.error('RDKit INIT FAIL', e); process.exit(1); });

function main() {
  console.log('===== A. 定向用例 =====');
  console.log('名称'.padEnd(26) + '双键  either  |  旧:交叉/最差角  |  新:交叉/最差角');
  let badA = 0;
  for (const [name, sm] of CASES) {
    const r = audit(sm);
    if (r.err) { console.log(name.padEnd(26) + 'ERR ' + r.err + '  ' + sm); badA++; continue; }
    const flag = r.xNew > 0 ? '  <<< 仍有交叉' : '';
    if (r.xNew > 0) badA++;
    console.log(name.padEnd(26) + String(r.ndbl).padStart(3) + '   ' + String(r.either).padStart(3)
      + '    |  ' + String(r.xOld).padStart(2) + '/' + String(r.worstOld).padStart(5) + '°'
      + '      |  ' + String(r.xNew).padStart(2) + '/' + String(r.worstNew).padStart(5) + '°' + flag);
    if (r.overNew.length) console.log('     残留:', JSON.stringify(r.overNew.slice(0, 4)));
  }

  console.log('\n===== B. 全库抽样 =====');
  const shardIds = fs.readdirSync(path.join(DIR, 'data')).filter(f => /^products\.\d+\.json$/.test(f));
  const all = [];
  for (const f of shardIds) {
    try { all.push(...JSON.parse(fs.readFileSync(path.join(DIR, 'data', f), 'utf8'))); } catch (e) { }
  }
  console.log('总记录', all.length, '(分片', shardIds.length, ')');
  const N = 600, step = Math.max(1, Math.floor(all.length / N));
  let n = 0, totD = 0, totE = 0, totXO = 0, totXN = 0, wO = 0, wN = 0;
  const leftover = [], errs = [];
  for (let i = 0; i < all.length && n < N; i += step) {
    const sm = all[i].s || '';
    if (sm.length < 12 || sm.length > 150 || sm.includes('.')) continue;
    let r; try { r = audit(sm); } catch (e) { errs.push(sm + ' :: ' + e.message); continue; }
    if (r.err) { errs.push(sm + ' :: ' + r.err); continue; }
    n++;
    totD += r.ndbl; totE += r.either; totXO += r.xOld; totXN += r.xNew;
    wO = Math.max(wO, r.worstOld); wN = Math.max(wN, r.worstNew);
    if (r.xNew > 0 && leftover.length < 6) leftover.push({ sm, name: (all[i].n || '').slice(0, 30), x: r.xNew, w: r.worstNew });
  }
  console.log('抽样分子=' + n + ' 双键=' + totD + ' either标记=' + totE);
  console.log('旧管线交叉双键=' + totXO + ' (最差夹角 ' + wO.toFixed(1) + '°)');
  console.log('新管线交叉双键=' + totXN + ' (最差夹角 ' + wN.toFixed(1) + '°)');
  for (const r of leftover) console.log('  残留:', JSON.stringify(r));
  console.log('异常=' + errs.length);
  for (const e of errs.slice(0, 5)) console.log('  ', e);
  process.exit(badA || totXN ? 1 : 0);
}
