#!/usr/bin/env node
/* 产物库上线前的"画得出来"闸门 —— 用**线上同一个渲染器**(search/app.js::molToSVG) 实渲抽样,
 * 再用独立实现的几何度量判它画得对不对。任一项不过就退出码 1, 阻断 push。
 *
 * 起因(2026-11-26): 用户发现卡片上异丙基的键穿过苯并噻唑环 —— 结构数据没错, 是渲染器选坐标时
 *   只优化键长/角度、不管交叉。当时没有闸门, 只能事后全库抽样才发现(0.81%)。此后规定:
 *   **产物库上传前必须自己先跑本闸门, 不过不许推。**
 *
 * 用法:
 *   node render_gate.js --shards 79 80 81 82 83          # 指定分片(默认每片抽 1200 条)
 *   node render_gate.js --shards 84 --sample 3000         # 加大抽样
 *   node render_gate.js --shards 84 --full                # 全量逐条(慢: 约 30ms/条)
 *   node render_gate.js --calibrate --shards 79 80        # 只打印各指标分位, 用于定阈值
 *   node render_gate.js --appjs <path> ...                # 指定别的 app.js(回归对照用)
 *
 * 判据(阈值见 THRESH; 2026-11-26 起改为零容忍, 用户要求"全部挨个更正"):
 *   ① 键线交叉 = 0            两条不共端点的键线段真相交
 *   ② 键长比 0.70–2.60        逐键 d/中位; 1.65~2.60 记 WARN(拥挤分子外推取舍), >2.60 FAIL
 *   ③ 非键最近距离 ≥ 0.40×中位 原子叠在一起
 *   ④ 环内角平均偏差 ≤ 8°      5/6/7 元环相对 108/120/128.57
 *   ⑤ SVG 必须非空且含键路径
 *   → 任一片抽样里出现 1 条不合格或交叉即 FAIL, 列出 code 逐条修; 不再按缺陷率放行。
 */
'use strict';
const fs = require('fs'), path = require('path');

const ROOT = __dirname;
const SEARCH = path.join(ROOT, 'search');
const DATA = path.join(SEARCH, 'data');

/* ---------- 阈值(硬性) ---------- */
const THRESH = {
  cross: 0,          // 严格 0
  bondMin: 0.70,
  bondMax: 2.60,     // 硬上限: 超过就是"断线/游丝"级缺陷
  bondWarn: 1.65,    // 1.65~2.60 记 WARN: 拥挤分子"推远片断"必然拉长系链键, 属可接受的取舍
  minNB: 0.40,       // × 中位键长
  ringAngMean: 8.0,  // 度
  /* ★ 零容忍(2026-11-26 用户要求"全部挨个更正"): 不合格数、交叉数都必须精确为 0,
   *   不再用比例阈值放行。WARN(长键 1.65~2.6, 拥挤分子外推的合法取舍)仍不计入不合格。 */
};

function parseBlk(blk){
  const L = blk.split('\n');
  const c = (L[3] || '').trim().split(/\s+/).map(Number);
  const nA = c[0] | 0, nB = c[1] | 0;
  const xy = [];
  for (let i = 0; i < nA; i++){
    const t = (L[4+i] || '').trim().split(/\s+/);
    xy.push([+t[0], +t[1]]);
  }
  const bonds = [];
  for (let i = 0; i < nB; i++){
    const s = L[4+nA+i] || '';
    const a = parseInt(s.slice(0,3), 10) - 1, b = parseInt(s.slice(3,6), 10) - 1;
    if (a >= 0 && a < nA && b >= 0 && b < nA && a !== b) bonds.push([a, b]);
  }
  return { nA, xy, bonds };
}
function segCross(p1, p2, p3, p4){
  const d = (a, b, c) => (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]);
  const d1 = d(p3,p4,p1), d2 = d(p3,p4,p2), d3 = d(p1,p2,p3), d4 = d(p1,p2,p4);
  return ((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) && ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0));
}
function findRings(adj, n, t){
  const out = [], seen = new Set();
  for (let s = 0; s < n; s++){
    const path_ = [s], used = new Set([s]);
    (function dfs(v){
      if (path_.length === t){ if (adj[v].includes(s)) out.push([...path_]); return; }
      for (const nb of adj[v]){
        if (nb === s || used.has(nb)) continue;
        used.add(nb); path_.push(nb); dfs(nb); path_.pop(); used.delete(nb);
      }
    })(s);
  }
  const uniq = [];
  for (const r of out){
    const k = [...r].sort((a,b)=>a-b).join(',');
    if (!seen.has(k)){ seen.add(k); uniq.push(r); }
  }
  return uniq;
}
/* 度量: 全部基于"实际用于出图的那份坐标" */
function measure(blk, svg){
  const P = parseBlk(blk);
  const { nA, xy, bonds } = P;
  const m = { nA, nB: bonds.length, cross: 0, bondMin: 1, bondMax: 1, minNB: 9, ringAng: 0, svgOK: 0 };
  m.svgOK = (svg && svg.indexOf('bond-') >= 0) ? 1 : 0;
  if (!bonds.length) return m;
  const ds = bonds.map(([a,b]) => Math.hypot(xy[a][0]-xy[b][0], xy[a][1]-xy[b][1]));
  const med = ds.slice().sort((x,y)=>x-y)[(ds.length/2)|0];
  if (!(med > 0.2)) return m;
  for (const d of ds){
    const r = d / med;
    if (r < m.bondMin) m.bondMin = r;
    if (r > m.bondMax) m.bondMax = r;
  }
  // 交叉
  for (let i = 0; i < bonds.length; i++) for (let j = i+1; j < bonds.length; j++){
    const [a,b] = bonds[i], [e,f] = bonds[j];
    if (a === e || a === f || b === e || b === f) continue;
    if (segCross(xy[a], xy[b], xy[e], xy[f])) m.cross++;
  }
  // 非键最近距离
  const bonded = new Set();
  for (const [a,b] of bonds){ bonded.add(a+'-'+b); bonded.add(b+'-'+a); }
  for (let i = 0; i < nA; i++) for (let j = i+1; j < nA; j++){
    if (bonded.has(i+'-'+j)) continue;
    const d = Math.hypot(xy[i][0]-xy[j][0], xy[i][1]-xy[j][1]) / med;
    if (d < m.minNB) m.minNB = d;
  }
  // 环内角
  const adj = Array.from({length: nA}, () => []);
  for (const [a,b] of bonds){ adj[a].push(b); adj[b].push(a); }
  let sum = 0, cnt = 0;
  for (const t of [5,6,7]){
    const ideal = t === 6 ? 120 : (t === 5 ? 108 : 128.571);
    for (const r of findRings(adj, nA, t)){
      for (let k = 0; k < t; k++){
        const a = r[(k-1+t)%t], b = r[k], c = r[(k+1)%t];
        const v1 = [xy[a][0]-xy[b][0], xy[a][1]-xy[b][1]];
        const v2 = [xy[c][0]-xy[b][0], xy[c][1]-xy[b][1]];
        const n1 = Math.hypot(v1[0],v1[1]), n2 = Math.hypot(v2[0],v2[1]);
        if (n1 < 1e-6 || n2 < 1e-6) continue;
        let cs = (v1[0]*v2[0]+v1[1]*v2[1])/(n1*n2);
        cs = Math.max(-1, Math.min(1, cs));
        sum += Math.abs(Math.acos(cs)*180/Math.PI - ideal); cnt++;
      }
    }
  }
  m.ringAng = cnt ? sum/cnt : 0;
  return m;
}
function failReasons(m){
  const r = [];
  if (!m.svgOK) r.push('SVG无键路径');
  if (m.cross > THRESH.cross) r.push('交叉' + m.cross);
  if (m.bondMin < THRESH.bondMin) r.push('短键' + m.bondMin.toFixed(2));
  if (m.bondMax > THRESH.bondMax) r.push('长键' + m.bondMax.toFixed(2));
  if (m.minNB < THRESH.minNB) r.push('原子过近' + m.minNB.toFixed(2));
  if (m.ringAng > THRESH.ringAngMean) r.push('环角' + m.ringAng.toFixed(1) + '°');
  return r;
}
function warnReasons(m){
  const r = [];
  if (!failReasons(m).length && m.bondMax > THRESH.bondWarn) r.push('长键' + m.bondMax.toFixed(2) + '(系链拉伸)');
  return r;
}

/* ---------- 参数 ---------- */
const argv = process.argv.slice(2);
function argVal(name, def){
  const i = argv.indexOf(name);
  return i >= 0 ? argv[i+1] : def;
}
const FULL = argv.includes('--full');
const CAL = argv.includes('--calibrate');
const SAMPLE = +(argVal('--sample', 1200));
const APPJS = argVal('--appjs', path.join(SEARCH, 'app.js'));
let shards = [];
{
  const i = argv.indexOf('--shards');
  if (i >= 0){
    for (let k = i+1; k < argv.length && !argv[k].startsWith('--'); k++) shards.push(+argv[k]);
  }
  if (!shards.length){
    // 默认取最后 5 片(即最近新增的库)
    const nums = fs.readdirSync(DATA).map(f => {
      const m = f.match(/^products\.(\d+)\.json$/); return m ? +m[1] : -1;
    }).filter(x => x >= 0).sort((a,b)=>a-b);
    shards = nums.slice(-5);
  }
}

/* ---------- 加载渲染器 ---------- */
globalThis.window = globalThis;
globalThis.addEventListener = () => {};
window.addEventListener = () => {};
window.document = { getElementById: () => null, createElement: () => ({ style: {} }), addEventListener: () => {} };
let src = fs.readFileSync(APPJS, 'utf8').replace(/'use strict';/, '');
src += '\n;globalThis.__setRDKit=(m)=>{RDKit=m;};';
(0, eval)(src);
if (typeof molToSVG !== 'function'){ console.error('FATAL: molToSVG 未加载 (' + APPJS + ')'); process.exit(2); }

const init = require(path.join(SEARCH, 'lib/rdkit/RDKit_minimal.js')).default;
init({ locateFile: f => path.join(SEARCH, 'lib/rdkit', f) }).then(mod => {
  globalThis.__setRDKit(mod); mod.prefer_coordgen(true);
  const rows = [];
  let totN = 0, totR = 0, totFail = 0;
  const bad = [];
  for (const sh of shards){
    const p = path.join(DATA, 'products.' + sh + '.json');
    if (!fs.existsSync(p)){ console.log('  ! 分片不存在: ' + p); continue; }
    const arr = JSON.parse(fs.readFileSync(p, 'utf8'));
    const step = FULL ? 1 : Math.max(1, Math.ceil(arr.length / SAMPLE));
    let n = 0, fail = 0;
    const agg = { cross: 0, bondMin: 1, bondMax: 1, minNB: 9, ringAng: 0, warn: 0 };
    for (let i = 0; i < arr.length; i += step){
      const it = arr[i];
      if (!it || !it.s){ continue; }
      let svg = null;
      try { svg = molToSVG(it.s, 240, 190); } catch(e){ svg = null; }
      const blk = globalThis.__PRP_LASTBLK;
      n++;
      if (!svg || !blk){ fail++; bad.push({ sh, code: it.code, why: ['渲染失败'] }); continue; }
      const m = measure(blk, svg);
      agg.cross += m.cross;
      if (m.bondMin < agg.bondMin) agg.bondMin = m.bondMin;
      if (m.bondMax > agg.bondMax) agg.bondMax = m.bondMax;
      if (m.minNB < agg.minNB) agg.minNB = m.minNB;
      if (m.ringAng > agg.ringAng) agg.ringAng = m.ringAng;
      const rs = failReasons(m);
      if (rs.length){ fail++; if (bad.length < 40) bad.push({ sh, code: it.code, why: rs, s: it.s }); }
      else { const ws = warnReasons(m); if (ws.length) agg.warn++; }
    }
    rows.push({ sh, n, fail, ...agg });
    totN += n; totR += fail; totFail += fail;
  }

  if (CAL){
    console.log('\n===== 指标极值(本批样本, 用于定阈值) =====');
    console.log('分片   渲染数  最差键长比  非键最近  环角偏差  交叉合计');
    for (const r of rows){
      console.log(String(r.sh).padEnd(6), String(r.n).padStart(6), '   ',
        r.bondMin.toFixed(3) + '/' + r.bondMax.toFixed(3), '  ',
        r.minNB.toFixed(3), '  ', r.ringAng.toFixed(2) + '°', '  ', r.cross);
    }
    process.exit(0);
  }

  console.log('\n===== 渲染闸门 (' + (FULL ? '全量' : '抽样 ' + SAMPLE + '/片') + ', app.js=' + path.basename(APPJS) + ') =====');
  console.log('分片   渲染数   不合格   警告   交叉   最差键长比    非键最近   环角偏差');
  for (const r of rows){
    console.log(
      String(r.sh).padEnd(6),
      String(r.n).padStart(6),
      String(r.fail).padStart(7),
      String(r.warn || 0).padStart(6),
      String(r.cross).padStart(6),
      (r.bondMin.toFixed(2) + '~' + r.bondMax.toFixed(2)).padStart(12),
      r.minNB.toFixed(2).padStart(10),
      (r.ringAng.toFixed(1) + '°').padStart(9)
    );
  }
  if (bad.length){
    console.log('\n不合格样例(最多 12 条):');
    for (const b of bad.slice(0, 12)) console.log('  ' + b.code + ' [片' + b.sh + '] ' + b.why.join(', ') + '  ' + (b.s || '').slice(0, 60));
  }
  console.log('\n合计: 渲染 ' + totN + ' 条, 不合格 ' + totR + ' 条 (' + (totN ? (100*totR/totN).toFixed(2) : '0') + '%)');
  console.log('阈值: 交叉=' + THRESH.cross + ' 键长 ' + THRESH.bondMin + '~' + THRESH.bondMax
    + ' 非键最近≥' + THRESH.minNB + ' 环角≤' + THRESH.ringAngMean + '°'
    + ' 长键 WARN>' + THRESH.bondWarn + ' FAIL>' + THRESH.bondMax);
  const crossAll = rows.reduce((a, r) => a + r.cross, 0);
  /* ★ 2026-11-26 用户明确要求「全部挨个更正，打分不过关不要上传」→ 零容忍:
   *   不合格数必须精确为 0, 交叉必须精确为 0。旧的比例闸门(≤0.1%/≤0.3%)已废除 ——
   *   它会放行"每片 1~3 条带缺陷", 与用户要求冲突。抽样命中一条即 FAIL, 会列出 code 逐条修。 */
  const pass = totN > 0 && totR === 0 && crossAll === 0;
  console.log('零容忍: 不合格必须=0 (实 ' + totR + ') | 交叉必须=0 (实 ' + crossAll + ')');
  console.log(pass ? '\nRENDER GATE PASS' : '\nRENDER GATE FAIL —— 有 ' + totR + ' 条不合格/交叉, 不许上传, 先按 code 逐条修或剔除');
  process.exit(pass ? 0 : 1);
}).catch(e => { console.error('FATAL:', e && e.message); process.exit(2); });
