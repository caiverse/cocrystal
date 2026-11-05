/* 用线上同一个 molToSVG() 抽查 products.35.json 新增结构。
 * 注意: __PRP_LASTBLK 是 **_prpClearBondStereo 之前** 的坐标块(审计钩子),
 * 实际渲染用的是清理后的块, 所以必须显式调用线上同一个 _prpClearBondStereo 复核。 */
const fs = require('fs'), path = require('path');
const ROOT = 'C:/Users/bc/cocrystal/search';
const TMP = 'C:/Users/bc/AppData/Local/Temp';
globalThis.window = globalThis;
globalThis.addEventListener = () => {};
window.addEventListener = () => {};
window.document = { getElementById: () => null, createElement: () => ({ style: {} }), addEventListener: () => {} };
let src = fs.readFileSync(path.join(ROOT, 'app.js'), 'utf8').replace(/'use strict';/, '');
src += '\n;globalThis.__setRDKit=(m)=>{RDKit=m;};';
(0, eval)(src);

const arr = JSON.parse(fs.readFileSync(path.join(ROOT, 'data/products.35.json'), 'utf8'));
const N = 400, pick = [];
for (let i = 0; i < N; i++) pick.push(arr[Math.floor(i * (arr.length - 1) / (N - 1))]);
// 额外加入 4 个含孤立 C=C 的已知"either"样例
const known = ['CCCC=CCC(NNCCCC)P1(=O)Oc2ccccc2-c2ccccc21',
               'CCCC=CCC(NNc1ccc(N)cc1)P1(=O)Oc2ccccc2-c2ccccc21',
               'CCCC=CCC(NNC(=O)NC(C)C)P1(=O)Oc2ccccc2-c2ccccc21',
               'CCC=CCCC(NNC(=O)Nc1ccccc1C)P1(=O)Oc2ccccc2-c2ccccc21'];
known.forEach((s, i) => pick.push({ code: 'KNOWN' + (i + 1), s }));

function eitherCount(blk) {
  if (!blk) return -1;
  const L = blk.split('\n');
  const c = (L[3] || '').trim().split(/\s+/).map(Number);
  const nA = c[0] | 0, nB = c[1] | 0;
  let n = 0;
  for (let i = 4 + nA; i < 4 + nA + nB && i < L.length; i++) {
    const s = L[i] || '';
    if (parseInt(s.slice(6, 9), 10) !== 2) continue;
    const st = parseInt(s.slice(9, 12), 10) || 0;
    if (st === 3 || st === 4) n++;
  }
  return n;
}

const init = require(path.join(ROOT, 'lib/rdkit/RDKit_minimal.js')).default;
init({ locateFile: f => path.join(ROOT, 'lib/rdkit', f) }).then(mod => {
  globalThis.__setRDKit(mod); mod.prefer_coordgen(true);
  let ok = 0, bad = [];
  let rawEither = 0, afterEither = 0, dumped = 0;
  const t0 = Date.now();
  for (const r of pick) {
    let svg = null;
    try { svg = molToSVG(r.s, 300, 240); } catch (e) { bad.push(['EXCEPTION', r.code, String(e)]); continue; }
    if (!svg || /NaN|Infinity|undefined/.test(svg)) { bad.push(['BAD_SVG', r.code, r.s]); continue; }
    const blk = globalThis.__PRP_LASTBLK;
    const raw = eitherCount(blk);
    const cleaned = (typeof _prpClearBondStereo === 'function') ? _prpClearBondStereo(blk) : blk;
    const aft = eitherCount(cleaned);
    if (raw > 0) rawEither++;
    if (aft > 0) { afterEither++; bad.push(['EITHER_AFTER_CLEAR', r.code, r.s]); }
    // 把含孤立 C=C 的样例 SVG 落盘, 供转 PNG 目检
    if (r.code && r.code.startsWith('KNOWN') && dumped < 4) {
      fs.writeFileSync(path.join(TMP, '_new_' + r.code + '.svg'), svg);
      dumped++;
    }
    ok++;
  }
  const dt = Date.now() - t0;
  console.log('抽样数', pick.length, '通过', ok);
  console.log('清理前含 either 双键的分子数', rawEither, '(这些正是由线上 _prpClearBondStereo 修掉的)');
  console.log('清理后仍含 either 双键的分子数', afterEither, afterEither === 0 ? 'PASS' : 'FAIL');
  console.log('空/NaN/异常 SVG', bad.length, '耗时', dt, 'ms  平均', (dt / pick.length).toFixed(2), 'ms/分子');
  if (bad.length) console.log('异常样例', bad.slice(0, 8));
  console.log('已落盘目检 SVG', dumped);
  process.exit(bad.length ? 1 : 0);
});
