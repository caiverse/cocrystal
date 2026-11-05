/* 产物库检索站 — 前端逻辑
 * 数据: data/products.json (聚合后的统一 schema, 短键名)
 * 结构式: RDKit WASM 本地渲染
 */
'use strict';

let DATA = [];
let BLOBS = [];          // 每条记录的检索文本(小写)
let FILTERED = [];       // 当前结果在 DATA 中的下标(经排序后的显示顺序)
let MATCHED = [];        // 当前命中的下标(库顺序, 未经排序)
let PER_PAGE = 24;
let page = 0;
let RDKit = null;
let LIBRARY_ORDER = [];   // 产物库默认展示顺序(来自 stats.json), 供默认排序分组用
let ATOM_CNT = [];        // [结构检索] 每条记录的重原子数(来自分子式, 不含 H), 用于子结构/精确预筛
window.SIM_MAP = {};      // [结构检索] 记录 id -> 相似度(0..1), 仅相似度模式

const EXTRA_LABEL = {
  avail: '市售/文献', alkyl: '推荐烷化剂', ether_class: '类别', r_group: '烷氧基 R',
  derive_dir: '衍生化方向', note: '备注', tpsa: 'TPSA', hbd: 'HBD', hba: 'HBA',
  ccdc_raw: 'CCDC号', ccdc_doi: 'CCDC DOI', ref: '文献', inchikey: 'InChIKey',
  pubchem_cid: 'PubChem CID', sub_type: '取代基类型', product_type: '产物类型',
  iupac: 'IUPAC', price: '单价(元/g)'
};

function escapeHTML(s){
  return String(s == null ? '' : s).replace(/[&<>"']/g, c => (
    {'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]
  ));
}
// 分子式: 正体显示, 数字转下标 (C14H11NO2 -> C<sub>14</sub>H<sub>11</sub>NO<sub>2</sub>)
function fmtFormula(s){
  if (s == null) return '';
  return escapeHTML(s).replace(/(\d+)/g, '<sub>$1</sub>');
}
function debounce(fn, ms){
  let t; return function(){ clearTimeout(t); t = setTimeout(fn, ms); };
}

/* ---------- RDKit 初始化 ---------- */
async function initRDKit(){
  if (typeof initRDKitModule !== 'function'){
    console.warn('RDKit 脚本未加载');
    return;
  }
  try{
    RDKit = await initRDKitModule({ locateFile: (f) => 'lib/rdkit/' + f });
    RDKit.prefer_coordgen(true);   // 强制使用新版坐标生成器: 与首页侧栏结构图同一画法(取向规范、自适应居中、大小一致)
    console.log('RDKit 就绪');
    if (DATA.length) render();   // 数据已就绪则重绘当前页结构式
  }catch(e){
    console.warn('RDKit 初始化失败:', e);
    RDKit = null;
  }
}

/* ---------- 悬挂(末端)环"系链两侧夹角 = 90°+180°/n"(苯环即 120°/120°)规范化 ----------
 * 用户要求: 悬挂苯环系链左右两个角度必须是 120°。coordgen 有两类偏差, 分两道工序处理:
 * [工序一·重投影] 内部不规则(被压扁/拉伸)的环 —— 刚体旋转改变不了环顶点相对几何, 修不了;
 *   故在 SVG 上直接改写坐标, 把环重投影为正 n 边形(系链为对称轴, ipso 原位不动, 环心在系链反侧):
 *   两侧夹角精确 = 90°+180°/n。逐原子平移 Δ(环原子=新顶点-旧位置), 取代基随挂靠环原子整体平移;
 *   环内键(含双键内侧线)按(旧a,旧b)->(新a,新b)相似变换重写 d, 标签字形与其余片断键按 Δ 平移,
 *   无 class 小黑点三角按质心就近随 Δ 平移。安全阀: 已达标(误差<2°)不动; Δ>2.2 键长或新顶点
 *   距母核原子<0.55 键长则放弃(留待工序二); 防镜像(shoelace 方向核验)。
 * [工序二·刚体旋转] 内部已规则但取向偏离径向的环 —— 包 <g transform="rotate"> 绕 ipso 旋转到
 *   "系链沿环心->ipso 径向"。关键(2026-11-05): 锚点方向用键线方向(_prpBondDir)而非锚点估计坐标 ——
 *   RDKit 给 N/O 标签留位会把键线截短, 截短端点反推的锚点位置内缩 15~22px, 使方向偏差约 20°,
 *   曾把本已外向的环额外旋转 20° 造成 142°/98° 歪斜("与N相连的苯环两个角度不是120°");
 *   键线截短沿键轴故方向精确, 配合精确环原子位置即得 120°/120°。旋转时小黑点三角随质心归属一并转。
 * 环检测: "真正悬挂的末端环"(单连接点、不与它环融合、去环后其余连通、取代基≤12 原子)会被处理;
 *   取代基里再含悬挂环(苄氧基/苯硫基苯基等联芳基)由"嵌套作业"处理: 元素归最内层、rotate 组嵌套复合。 */
function _prpParseV2000(blk){ const L=blk.split('\n'); const c=L[3].trim().split(/\s+/).map(Number); const nA=c[0],nB=c[1]; return {nAtoms:nA,nBonds:nB,atomLines:L.slice(4,4+nA),bondLines:L.slice(4+nA,4+nA+nB)}; }
/* 端基键长度修正(**双向**): 过短拉到中位, 过长压回中位(只动度=1 的原子, 只沿原键方向)。
 * 起因(2026-11-05 用户: "C=O 双键一个被画成特别长, 一个被画成那么短"):
 *   坐标生成器偶发把某个**端基键**压到 0.8 Å(中位 1.5) 或拉到 2.5 倍, 画出来就是半截/一长条。
 * 抽样 300 分子统计: 短键(<0.72×中位) 165 条 **全是端基**; 过长键在 coordgen 坐标里 120 条则**全是非端基**
 *   (那是布局问题, 靠"逐分子选坐标"解决, 见 _prpPickBlk), 端基过长只是零星, 这里一并处理。
 * 端基只连一个原子, 沿原方向移动**不会移动其它任何原子、不改变任何键角**。 */
function _prpLengthenTerminalBonds(blk){
  const L = blk.split('\n');
  const c = (L[3]||'').trim().split(/\s+/).map(Number);
  const nA = c[0]|0, nB = c[1]|0;
  if (!(nA > 1 && nB > 0)) return blk;
  const xy = [];
  for (let i = 0; i < nA; i++){
    const t = (L[4+i]||'').trim().split(/\s+/);
    const x = +t[0], y = +t[1];
    if (!isFinite(x) || !isFinite(y)) return blk;
    xy.push([x, y]);
  }
  const nb = Array.from({length: nA}, () => []);
  const bIdx = new Int32Array(nA).fill(-1);   // 端基原子 -> 它那根键的序号(用于取键级)
  const ds = [];
  for (let i = 0; i < nB; i++){
    const s = L[4+nA+i] || '';
    const a = parseInt(s.slice(0,3), 10) - 1, b = parseInt(s.slice(3,6), 10) - 1;
    if (!(a >= 0 && a < nA && b >= 0 && b < nA) || a === b) continue;
    nb[a].push(b); nb[b].push(a);
    bIdx[a] = i; bIdx[b] = i;
    ds.push(Math.hypot(xy[a][0]-xy[b][0], xy[a][1]-xy[b][1]));
  }
  if (ds.length < 2) return blk;
  ds.sort((p, q) => p - q);
  const med = ds[(ds.length/2)|0];
  if (!(med > 0.2)) return blk;
  const LO = 0.85 * med, HI = 1.25 * med;    // 正常区间; 越界才修正
  const CLR = 0.50 * med;                    // 允许的最近非键合间隙
  const clearOf = (p, k) => {
    let m = 1e9;
    for (let q = 0; q < nA; q++){
      if (q === k) continue;
      const d = Math.hypot(xy[q][0]-p[0], xy[q][1]-p[1]);
      if (d < m) m = d;
    }
    return m;
  };
  let fixed = 0;
  for (let k = 0; k < nA; k++){
    if (nb[k].length !== 1) continue;                    // 只处理端基(度=1)
    const j = nb[k][0];
    const dx = xy[k][0]-xy[j][0], dy = xy[k][1]-xy[j][1];
    const d = Math.hypot(dx, dy);
    if (d < 1e-6 || (d >= LO && d <= HI)) continue;
    const cur = clearOf(xy[k], k);
    const f = med / d;
    const THR_HI = Math.min(cur, CLR);                       // 首选: 不比现在更挤, 且不低于半键长
    const THR_LO = Math.max(0.45 * med, 0.85 * cur);         // 次选: 允许略挤, 但不得跌破 0.45 键长
    /* 转向搜空位的三步:
     *   ① 原方向延长; ② 绕母原子小角度偏转(偏角由小到大, 优先小偏角);
     *   ③ 上面都不行就来全周细扫(3° 一格), 取最空的方向。
     * 单键端基(C/F/Cl/Br 这类可自由旋转的)才允许大角度; 双键端基(=O 之类)只在 ±30° 内偏,
     *   否则羰基会被摆到环的另一侧, 看着别扭(化学上没错, 但不像人画的)。 */
    const isMulti = (parseInt((L[4+nA+bIdx[k]]||'').slice(6,9), 10) || 1) > 1;
    const MAX_DEG = isMulti ? 30 : 180;
    let np = null, bestLoose = null, bestClr = -1;
    const tryAt = (pt) => {
      const cl = clearOf(pt, k);
      if (cl >= THR_HI) { np = pt; return true; }
      if (cl >= THR_LO && cl > bestClr){ bestClr = cl; bestLoose = pt; }
      return false;
    };
    {
      const base = Math.atan2(dy, dx);
      tryAt([xy[j][0] + dx * f, xy[j][1] + dy * f]);   // ① 原方向延长/压回
      if (!np){
        for (const deg of [6, 12, 18, 24, 30, 36, 45, 60]){
          if (deg > MAX_DEG) break;
          for (const sgn of [1, -1]){
            if (tryAt([xy[j][0] + med*Math.cos(base+sgn*deg*Math.PI/180),
                       xy[j][1] + med*Math.sin(base+sgn*deg*Math.PI/180)])) break;
          }
          if (np) break;
        }
      }
      if (!np && MAX_DEG > 60){                       // ③ 全周细扫
        for (let a = 0; a < 360; a += 3){
          if (tryAt([xy[j][0] + med*Math.cos(a*Math.PI/180), xy[j][1] + med*Math.sin(a*Math.PI/180)])) break;
        }
      }
      if (!np) np = bestLoose;
      if (!np){                                       // 最后退到"只改一半"(距离上仍改善)
        const f2 = 1 + 0.5 * (f - 1);
        tryAt([xy[j][0] + dx * f2, xy[j][1] + dy * f2]);
        if (!np) np = bestLoose;
      }
    }
    if (!np) continue;
    xy[k][0] = np[0]; xy[k][1] = np[1];
    L[4+k] = xy[k][0].toFixed(4).padStart(10) + xy[k][1].toFixed(4).padStart(10) + (L[4+k]||'').slice(20);
    fixed++;
  }
  if (fixed && globalThis.__PRP_DBG) console.log('[prp] terminal bonds fixed:', fixed);
  return L.join('\n');
}
/* 给一份坐标打"画出来有多难看"的分(越低越好), 用于在多个坐标生成器之间逐分子取优。
 * 分项:
 *   ① 键长离群: 取 |d/中位 - 1| 的最大值, 过长的罚 ×1.5(半截键只是短一截, 过长键会顶到别的基团);
 *   ② 环内角: 5/6/7 元环各顶点内角与理想值(108/120/128.57)之差的平均, 每度 2 分;
 *   ③ 原子挤在一起: 非键合最近距离不足 0.5 键长的部分, 每 0.1 键长 20 分。
 * 该分数只依赖几何, 与绘制样式无关。 */
function _prpScoreBlk(blk){
  if (!blk) return 1e9;
  let P; try { P = _prpParseV2000(blk); } catch(e){ return 1e9; }
  const nA = P.nAtoms;
  if (!(nA > 1)) return 1e9;
  const xy = [];
  for (let i = 0; i < nA; i++){
    const t = (P.atomLines[i]||'').trim().split(/\s+/);
    const x = +t[0], y = +t[1];
    if (!isFinite(x) || !isFinite(y)) return 1e9;
    xy.push([x, y]);
  }
  const adj = _prpBuildAdj(P.bondLines, nA);
  const ds = [];
  for (const b of P.bondLines){
    const m = b.match(/^\s*(\d+)\s+(\d+)/); if (!m) continue;
    const a = +m[1]-1, c = +m[2]-1;
    if (!(a >= 0 && a < nA && c >= 0 && c < nA)) continue;
    ds.push(Math.hypot(xy[a][0]-xy[c][0], xy[a][1]-xy[c][1]));
  }
  if (ds.length < 2) return 1e9;
  const sorted = ds.slice().sort((p, q) => p - q);
  const med = sorted[(sorted.length/2)|0];
  if (!(med > 0.2)) return 1e9;
  let bondPen = 0;
  for (const d of ds){
    const r = d / med;
    const dev = r < 1 ? (1 - r) : (r - 1) * 1.5;
    if (dev > bondPen) bondPen = dev;
  }
  let angSum = 0, angN = 0;
  for (const t of [5, 6, 7]){
    let rings; try { rings = _prpFindRings(adj, nA, t); } catch(e){ rings = []; }
    const ideal = t === 6 ? 120 : (t === 5 ? 108 : 128.571);
    for (const r of rings){
      for (let k = 0; k < t; k++){
        const a = r[(k-1+t)%t], b = r[k], cc = r[(k+1)%t];
        const v1x = xy[a][0]-xy[b][0], v1y = xy[a][1]-xy[b][1];
        const v2x = xy[cc][0]-xy[b][0], v2y = xy[cc][1]-xy[b][1];
        const n1 = Math.hypot(v1x, v1y), n2 = Math.hypot(v2x, v2y);
        if (n1 < 1e-6 || n2 < 1e-6) continue;
        let cs = (v1x*v2x + v1y*v2y) / (n1*n2);
        cs = Math.max(-1, Math.min(1, cs));
        angSum += Math.abs(Math.acos(cs) * 180 / Math.PI - ideal); angN++;
      }
    }
  }
  const angMean = angN ? angSum / angN : 0;
  let minNB = 1e9;
  if (nA <= 140){
    const bonded = new Set();
    adj.forEach((lst, i) => { for (const j of lst) bonded.add(i < j ? i+'-'+j : j+'-'+i); });
    for (let i = 0; i < nA; i++) for (let j = i+1; j < nA; j++){
      if (bonded.has(i+'-'+j)) continue;
      const d = Math.hypot(xy[i][0]-xy[j][0], xy[i][1]-xy[j][1]);
      if (d < minNB) minNB = d;
    }
  }
  const overPen = (minNB < 1e9) ? Math.max(0, 0.5 - minNB/med) / 0.1 : 0;
  return 100*bondPen + 2*angMean + 20*overPen;
}
/* 候选 C: 键长松弛 —— 把**非环键**的长度拉向中位值(两端各让一半), 环内键一律不动。
 * 用途: 两种坐标生成器都会偶尔留下"非端基键偏长/偏短"(如 coordgen 把某根 C-C 拉到 2.5 倍),
 *   而这种键没法靠移动端基修。松弛只改非环键长度, 不动环骨架, 因此不会把六边形搞歪;
 *   是否真的变好由 _prpScoreBlk 裁决(重叠会被罚分, 所以松弛出重叠不会被选中)。 */
function _prpRelaxBonds(blk){
  let P; try { P = _prpParseV2000(blk); } catch(e){ return blk; }
  const nA = P.nAtoms, L = blk.split('\n');
  if (!(nA > 1)) return blk;
  const xy = [];
  for (let i = 0; i < nA; i++){
    const t = (P.atomLines[i]||'').trim().split(/\s+/);
    const x = +t[0], y = +t[1];
    if (!isFinite(x) || !isFinite(y)) return blk;
    xy.push([x, y]);
  }
  const adj = _prpBuildAdj(P.bondLines, nA);
  const bonds = [];
  for (const b of P.bondLines){
    const m = b.match(/^\s*(\d+)\s+(\d+)/); if (!m) continue;
    const a = +m[1]-1, c = +m[2]-1;
    if (!(a >= 0 && a < nA && c >= 0 && c < nA) || a === c) continue;
    bonds.push([a, c]);
  }
  if (bonds.length < 2) return blk;
  /* 环内键判定: 去掉该键后两端是否仍连通(一次 BFS), 比枚举 3~8 元环快得多也稳得多。 */
  const isRing = (a, c) => {
    const seen = new Uint8Array(nA); seen[a] = 1;
    const st = [a];
    while (st.length){
      const u = st.pop();
      for (const v of adj[u]){
        if ((u === a && v === c) || (u === c && v === a)) continue;
        if (v === c) return true;
        if (seen[v]) continue;
        seen[v] = 1; st.push(v);
      }
    }
    return false;
  };
  const free = bonds.filter(([a, c]) => !isRing(a, c));
  if (!free.length) return blk;
  /* 环上原子一律钉死: 松弛只许挪"不在环里"的原子。
   *   起因(2026-11-05): 早先两端都对半挪, 会把环原子拖走 → 六边形被压歪、环内角罚分暴涨,
   *   于是"键长更匀"的松弛结果反而被 _prpScoreBlk 判负(实测 49.95 输给 42.97)。钉死环后
   *   松弛只在取代基内部动, 环保持正六边形, 才可能真的赢。 */
  const inRing = new Uint8Array(nA);
  for (const [a, c] of bonds) if (isRing(a, c)){ inRing[a] = 1; inRing[c] = 1; }
  const adjustable = free.filter(([a, c]) => !inRing[a] || !inRing[c]);
  if (!adjustable.length) return blk;
  const len = (a,c) => Math.hypot(xy[a][0]-xy[c][0], xy[a][1]-xy[c][1]);
  const ds = bonds.map(([a,c]) => len(a,c)).sort((p,q) => p-q);
  const med = ds[(ds.length/2)|0];
  if (!(med > 0.2)) return blk;
  const orig = xy.map(p => [p[0], p[1]]);
  for (let it = 0; it < 40; it++){
    let moved = 0;
    for (const [a,c] of adjustable){
      const d = len(a,c);
      if (d < 1e-6) continue;
      const diff = d - med;
      if (Math.abs(diff) < 0.01 * med) continue;
      const ux = (xy[c][0]-xy[a][0]) / d, uy = (xy[c][1]-xy[a][1]) / d;
      // 只有一端自由时由它独担全部位移; 两端都自由才各让一半。环上原子不动。
      const fa = !inRing[a], fc = !inRing[c];
      const wa = fa && fc ? 0.25 : (fa ? 0.5 : 0);
      const wc = fa && fc ? 0.25 : (fc ? 0.5 : 0);
      const sh = diff;
      xy[a][0] += ux * sh * wa; xy[a][1] += uy * sh * wa;
      xy[c][0] -= ux * sh * wc; xy[c][1] -= uy * sh * wc;
      moved++;
    }
    if (!moved) break;
  }
  /* 防挤压收尾: 松弛只该让键长变匀, 不该把两个原子推到几乎重合。
   *   实测该分子(三氯苯基+三甲氧基苯基)松弛把键长从 0.665 修到 0.99, 但把 9-18 两个原子
   *   挤到 0.25 键长, 打分反而变差 → 被弃用。这里在"原坐标 → 松弛坐标"之间二分,
   *   取"刚好不挤"的最大比例: 键长能修多少修多少, 绝不牺牲原子间距。 */
  if (nA <= 140){
    const bonded = new Set();
    adj.forEach((lst, i) => { for (const j of lst) bonded.add(i < j ? i+'-'+j : j+'-'+i); });
    const minNB = (pts) => {
      let m = 1e9;
      for (let i = 0; i < nA; i++) for (let j = i+1; j < nA; j++){
        if (bonded.has(i+'-'+j)) continue;
        const d = Math.hypot(pts[i][0]-pts[j][0], pts[i][1]-pts[j][1]);
        if (d < m) m = d;
      }
      return m;
    };
    const FLOOR = 0.50 * med;   // 与端基修正的 CLR 同口径: 半个键长以内算"挤"
    if (minNB(xy) < FLOOR){
      if (minNB(orig) < FLOOR) return blk;      // 本来就在一起, 松弛解决不了, 别动
      let lo = 0, hi = 1;
      for (let k = 0; k < 10; k++){
        const t = (lo + hi) / 2;
        const pts = orig.map((p, i) => [p[0] + t*(xy[i][0]-p[0]), p[1] + t*(xy[i][1]-p[1])]);
        if (minNB(pts) >= FLOOR) lo = t; else hi = t;
      }
      for (let i = 0; i < nA; i++){
        xy[i][0] = orig[i][0] + lo*(xy[i][0]-orig[i][0]);
        xy[i][1] = orig[i][1] + lo*(xy[i][1]-orig[i][1]);
      }
    }
  }
  for (let i = 0; i < nA; i++){
    L[4+i] = xy[i][0].toFixed(4).padStart(10) + xy[i][1].toFixed(4).padStart(10) + (L[4+i]||'').slice(20);
  }
  return L.join('\n');
}
/* 在若干候选坐标里挑分最低的一份; 无候选则返回 null。 */
function _prpPickBlk(list){
  let best = null, bestS = 1e9;
  for (const b of list){
    if (!b) continue;
    const s = _prpScoreBlk(b);
    if (s < bestS){ bestS = s; best = b; }
  }
  return best;
}
/* 清掉 molblock 键块里"双键立体化学=either"的标记(第 4 列 = 3 或 4)。
 * 关键修正(2026-11-05, C=N 被画成交叉): 主链路是
 *   mol.get_new_coords()  ->  molblock  ->  RDKit.get_mol(blk).get_svg()
 * get_new_coords() 写 molblock 时, 对**未指定 E/Z** 的双键会在键块第 4 列写 stereo=3(either);
 * 再解析回来画图, RDKit 按 either 画成**交叉双键**(两条短线在中点交叉成 X)——
 * 肟 C=N、腙 C=N、查尔酮 C=C、亚苄基 C=C 等全部中招(用户所见"这个 C=N 双键变成交叉了")。
 * 清 0 后恢复平行双线; 2D 构象本身仍在坐标里, E/Z 由坐标重新推断, 视觉不变。
 * 只动 **键类型=2(双键)且 stereo∈{3,4}** 的行 —— 单键楔形/虚楔(1/6)、三键、芳香键一律不碰,
 * 故原子手性与楔形键显示不受影响。V2000 定宽: 原子1[0,3) 原子2[3,6) 类型[6,9) stereo[9,12)。 */
function _prpClearBondStereo(blk){
  const L=blk.split('\n');
  const c=(L[3]||'').trim().split(/\s+/).map(Number);
  const nA=c[0]|0, nB=c[1]|0;
  if(!(nA>0 && nB>0)) return blk;
  let n=0;
  for(let i=4+nA; i<4+nA+nB && i<L.length; i++){
    const s=L[i]; if(!s || s.length<12) continue;
    if(parseInt(s.slice(6,9),10)!==2) continue;
    const st=parseInt(s.slice(9,12),10);
    if(st===3 || st===4){ L[i]=s.slice(0,9)+'  0'+s.slice(12); n++; }
  }
  if(n && globalThis.__PRP_DBG) console.log('[prp] double-bond either-stereo cleared:',n);
  return L.join('\n');
}
function _prpBuildAdj(bl,n){ const a=Array.from({length:n},()=>[]); for(const b of bl){const m=b.match(/^\s*(\d+)\s+(\d+)\s+(\d+)/); if(!m)continue; const x=+m[1]-1,y=+m[2]-1; a[x].push(y); a[y].push(x);} return a; }
function _prpFindRings(adj,n,t){ const rings=[]; for(let s=0;s<n;s++){ const path=[s],ip=new Set([s]); (function dfs(v){ if(path.length===t){ if(adj[v].includes(s)) rings.push([...path]); return;} for(const nb of adj[v]){ if(nb===s||ip.has(nb))continue; ip.add(nb); path.push(nb); dfs(nb); path.pop(); ip.delete(nb);} })(s);} const u=[],seen=new Set(); for(const r of rings){const k=[...r].sort((a,b)=>a-b).join(','); if(!seen.has(k)){seen.add(k);u.push(r);}} return u; }
function _prpGetPaths(svg){ const re=/<path\b([^>]*)\/?>/g; const res=[]; let m; while((m=re.exec(svg))){ const tag=m[0]; const cls=(tag.match(/class='([^']*)'/)||[])[1]||''; const d=(tag.match(/d='([^']*)'/)||[])[1]||''; res.push({tag,cls,d}); } return res; }
function _prpAtomsInClass(cls){ const s=new Set(); let m; const re=/atom-(\d+)/g; while((m=re.exec(cls))) s.add(+m[1]); return s; }
function _prpBondEndpoints(d){ const e=[]; let m; const re=/(-?\d+\.?\d*)[,\s]+(-?\d+\.?\d*)/g; while((m=re.exec(d))) e.push([parseFloat(m[1]),parseFloat(m[2])]); return e; }
function _prpPtSegDist(p,a,b){ const vx=b[0]-a[0], vy=b[1]-a[1]; const wx=p[0]-a[0], wy=p[1]-a[1]; const c1=vx*wx+vy*wy; if(c1<=0) return Math.hypot(p[0]-a[0],p[1]-a[1]); const c2=vx*vx+vy*vy; if(c2<=c1) return Math.hypot(p[0]-b[0],p[1]-b[1]); const t=c1/c2; return Math.hypot(p[0]-(a[0]+t*vx),p[1]-(a[1]+t*vy)); }
/* 顶点定位(碳原子用)。
 * 判据: 顶点 = 被"最多条不同键路径共享"的候选端点。
 * 关键修正(2026-11-05 三轮): 旧写法按"离多少条**线段**近(tol=2.5)"打分, 而芳香双键/三键的内侧平行线
 *   是**装饰线**, 其端点根本不落在顶点上(实测侧向偏移 2.43px —— 恰好小于 2.5, 于是被当成与顶点重合);
 *   结果内侧线端点拿到满分 3/3、真顶点只有 2/3, 顶点被内缩 1.4px, 正六边形被量成 111°/129° 的假相,
 *   旋转角(依赖环心与 ipso)也随之偏 1~2°。
 *   改判据后: 真顶点被 2~3 条**不同路径**共享, 内侧线端点只被它自己那 1 条路径共享 —— 天然胜出。
 *   (tol 收紧到 1.0px: 只接受"几乎重合"的端点, 从而免疫平行线的 2~2.5px 侧向偏移。) */
function _prpGetAtomPos(paths, k){
  // 候选点 = 该原子各键路径的线段端点, 每个端点带一条"键轴方向"(按 mod 180° 归族)。
  // 判据: 顶点处有 ≥2 条**不同直线**相交; RDKit 为给标签留位会把连向 N/O 的键在**中点拆成两段**,
  //   中点只是同一条直线的延续 —— 只算 1 族, 不会被误判为顶点。
  // (2026-11-05 实测: 旧"路径计数"判据让拆分的键中点以 2:2 平票取胜, 顶点被拉偏 9~12px。)
  const c=[];
  for(const p of paths){
    if(p.cls.indexOf('bond-')<0) continue;
    const ix=_prpAtomsInClass(p.cls); if(!ix.has(k)) continue;
    const e=_prpBondEndpoints(p.d); if(e.length<2) continue;
    for(let i=0;i+1<e.length;i+=2){
      const A=e[i], B=e[i+1];
      let f=Math.atan2(B[1]-A[1], B[0]-A[0]); f=((f%Math.PI)+Math.PI)%Math.PI;
      c.push({x:A[0], y:A[1], fam:f});
      c.push({x:B[0], y:B[1], fam:f});
    }
  }
  if(!c.length) return null;
  const TOL=1.0, FAMTOL=18*Math.PI/180;
  let best=null, bestScore=-1;
  for(const q of c){
    const fams=[];
    for(const r of c){
      if(Math.hypot(r.x-q.x, r.y-q.y)>TOL) continue;
      let dup=false;
      for(const f of fams){ let d=Math.abs(f-r.fam); d=Math.min(d, Math.PI-d); if(d<FAMTOL){ dup=true; break; } }
      if(!dup) fams.push(r.fam);
    }
    if(fams.length>bestScore){
      bestScore=fams.length;
      let sx=0, sy=0, n=0;
      for(const r of c){ if(Math.hypot(r.x-q.x, r.y-q.y)<=TOL){ sx+=r.x; sy+=r.y; n++; } }
      best=[sx/n, sy/n];
    }
  }
  return best;
}
// 块坐标(y 向上) -> SVG(y 向下) 的相似变换稳健拟合(IRLS/Cauchy 权重, 镜像+旋转+缩放+平移)。
// 实测: 400 个分子上残差中位数 0.05px —— 即 RDKit 的 get_svg 只是把块坐标做了相似变换(缩放/旋转/平移),
// 因此"块坐标经变换"就是每个原子在 SVG 里的**精确**位置, 完全免疫标签截短、双键内侧线、键被中点拆段等干扰。
function _prpFitSim(B, S){
  const n=B.length; if(n<3) return null;
  let w=new Array(n).fill(1), st=null;
  for(let it=0; it<10; it++){
    let sw=0, mbx=0, mby=0, msx=0, msy=0;
    for(let i=0;i<n;i++){ sw+=w[i]; mbx+=w[i]*B[i][0]; mby+=w[i]*B[i][1]; msx+=w[i]*S[i][0]; msy+=w[i]*S[i][1]; }
    if(!(sw>0)) return null;
    mbx/=sw; mby/=sw; msx/=sw; msy/=sw;
    let num=0, den=0;
    for(let i=0;i<n;i++){ const bx=B[i][0]-mbx, by=-(B[i][1]-mby), sx=S[i][0]-msx, sy=S[i][1]-msy;
      num+=w[i]*(bx*sy-by*sx); den+=w[i]*(bx*sx+by*sy); }
    if(Math.abs(den)<1e-12 && Math.abs(num)<1e-12) return null;
    const th=Math.atan2(num,den), ca=Math.cos(th), sa=Math.sin(th);
    let sN=0, sD=0;
    for(let i=0;i<n;i++){ const bx=B[i][0]-mbx, by=-(B[i][1]-mby);
      const rx=bx*ca-by*sa, ry=bx*sa+by*ca; const sx=S[i][0]-msx, sy=S[i][1]-msy;
      sN+=w[i]*(rx*sx+ry*sy); sD+=w[i]*(rx*rx+ry*ry); }
    if(!(sD>1e-12)) return null;
    const sc=sN/sD;
    if(!isFinite(sc)||sc<=0) return null;
    const res=new Array(n);
    for(let i=0;i<n;i++){ const bx=B[i][0]-mbx, by=-(B[i][1]-mby);
      const rx=bx*ca-by*sa, ry=bx*sa+by*ca;
      res[i]=Math.hypot(S[i][0]-(msx+sc*rx), S[i][1]-(msy+sc*ry)); }
    const sorted=res.slice().sort((a,b)=>a-b), med=sorted[Math.floor(n/2)];
    const sig=Math.max(0.3, 1.4826*med);
    for(let i=0;i<n;i++){ const z=res[i]/(2.5*sig); w[i]=1/(1+z*z); }
    st={mbx,mby,msx,msy,ca,sa,sc,med};
  }
  if(!st) return null;
  const T=(p)=>{ const bx=p[0]-st.mbx, by=-(p[1]-st.mby);
    return [st.msx+st.sc*(bx*st.ca-by*st.sa), st.msy+st.sc*(bx*st.sa+by*st.ca)]; };
  T.med=st.med; T.scale=st.sc; T.deg=Math.atan2(st.sa,st.ca)*180/Math.PI;
  return T;
}
// 每个原子的 SVG 精确位置: 先用键线交点(度>=2 的原子可靠)拟合相似变换, 再由块坐标预测全部原子。
function _prpPositions(paths, P){
  const pos=new Array(P.nAtoms).fill(null);
  const base=new Array(P.nAtoms).fill(null);
  const deg=new Array(P.nAtoms).fill(0);
  for(const b of P.bondLines){ const m=b.match(/^\s*(\d+)\s+(\d+)\s+(\d+)/); if(!m) continue; deg[+m[1]-1]++; deg[+m[2]-1]++; }
  const blkPos=(k)=>{ const L=(P.atomLines[k]||'').trim().split(/\s+/).map(Number); return [L[0],L[1]]; };
  const B=[], S=[];
  for(let k=0;k<P.nAtoms;k++){
    let s=null;
    try{ s=_prpAtomExact(paths,k); }catch(e){ s=null; }
    base[k]=s;
    if(s && deg[k]>=2 && isFinite(s[0]) && isFinite(s[1])){ B.push(blkPos(k)); S.push(s); }
  }
  if(B.length>=3){
    const T=_prpFitSim(B,S);
    if(T && T.med<3.5){
      for(let k=0;k<P.nAtoms;k++){ const q=blkPos(k); pos[k]=[T(q)[0], T(q)[1]]; }
    }
  }
  for(let k=0;k<P.nAtoms;k++) if(!pos[k]) pos[k]=base[k];
  return pos;
}
// ipso->锚点 键线方向(精确): RDKit 对连向 N/O 等标签原子的键线做了"截短"(沿键轴), 但方向不变。
// 取该键全部路径端点中"离 ipso 最远"的一簇(距最远点 3px 内)做平均 —— 既得键轴方向, 又抵消双键双线的侧向偏移。
function _prpBondDir(paths, a, b, ipP){
  const pts=[];
  for(const p of paths){
    if(p.cls.indexOf('bond-')<0) continue;
    const at=_prpAtomsInClass(p.cls);
    if(!(at.has(a)&&at.has(b))) continue;
    const e=_prpBondEndpoints(p.d);
    for(const q of e) pts.push(q);
  }
  if(pts.length<2) return null;
  let D2=-1;
  for(const q of pts){ const d2=(q[0]-ipP[0])*(q[0]-ipP[0])+(q[1]-ipP[1])*(q[1]-ipP[1]); if(d2>D2) D2=d2; }
  const D=Math.sqrt(D2), lim=Math.max(0,D-3)*Math.max(0,D-3);
  let sx=0,sy=0,n=0;
  for(const q of pts){ const d2=(q[0]-ipP[0])*(q[0]-ipP[0])+(q[1]-ipP[1])*(q[1]-ipP[1]); if(d2>=lim){ sx+=q[0]; sy+=q[1]; n++; } }
  if(!n) return null;
  return Math.atan2(sy/n-ipP[1], sx/n-ipP[0]);
}
// 杂原子精定位: 用"该原子各键的轴线的最小二乘交点"求中心。
//   动机(2026-11-05): 字形 bbox 中心对 N/O 有 2~4px 系统偏差(字母有上下伸展, bbox 不等于字面几何中心),
//   这会把含杂原子环的环心带偏、把重投影顶点带偏, 并把该处内部角量成 136°/104° 的假象。
//   键线在标签处被"沿键轴截短", 但轴线方向不变; 且每条键取**最长线段**(主线), 双键内侧平行线更短会被排除;
//   再要求该线段一端几乎落在另一原子的精确位置上(碳定位精确), 从而剔除偏移的内侧线。多键轴线交点即原子中心。
function _prpAxisPos(paths, k){
  const L=[];
  for(const p of paths){
    if(p.cls.indexOf('bond-')<0) continue;
    const at=[..._prpAtomsInClass(p.cls)];
    if(at.length!==2 || at.indexOf(k)<0) continue;
    const oth=(at[0]===k)?at[1]:at[0];
    const op=_prpGetAtomPos(paths, oth); if(!op) continue;
    const e=_prpBondEndpoints(p.d); if(e.length<2) continue;
    let bs=null, bl=-1;
    for(let i=0;i+1<e.length;i+=2){ const len=Math.hypot(e[i+1][0]-e[i][0], e[i+1][1]-e[i][1]); if(len>bl){ bl=len; bs=[e[i],e[i+1]]; } }
    if(!bs || bl<1e-6) continue;
    const d0=Math.hypot(bs[0][0]-op[0], bs[0][1]-op[1]), d1=Math.hypot(bs[1][0]-op[0], bs[1][1]-op[1]);
    let ref, dx, dy;
    if(d0<=d1 && d0<2.0){ ref=bs[0]; dx=bs[1][0]-bs[0][0]; dy=bs[1][1]-bs[0][1]; }
    else if(d1<d0 && d1<2.0){ ref=bs[1]; dx=bs[0][0]-bs[1][0]; dy=bs[0][1]-bs[1][1]; }
    else continue;                                  // 两端都不贴另一原子中心 -> 内侧平行线, 弃用
    const dl=Math.hypot(dx,dy); if(dl<1e-6) continue;
    L.push({x:ref[0], y:ref[1], nx:-dy/dl, ny:dx/dl});
  }
  if(L.length<2) return null;
  let a11=0,a12=0,a22=0,b1=0,b2=0;
  for(const l of L){
    a11+=l.nx*l.nx; a12+=l.nx*l.ny; a22+=l.ny*l.ny;
    const c=l.nx*l.x+l.ny*l.y; b1+=l.nx*c; b2+=l.ny*c;
  }
  const det=a11*a22-a12*a12;
  if(Math.abs(det)<1e-9) return null;
  const x=(a22*b1-a12*b2)/det, y=(a11*b2-a12*b1)/det;
  if(!isFinite(x)||!isFinite(y)) return null;
  return [x,y];
}
// 杂原子标签字形中心(兜底): RDKit 把 N/O/P 等标签画成 <path class='atom-k' d='...字形轮廓...'>(绝对坐标)。
// 碳原子无标签, 形状路径不存在, 返回 null。
function _prpLabelPos(paths, k){
  const key='atom-'+k;
  for(const p of paths){
    if(p.cls.indexOf('bond-')>=0) continue;
    const at=_prpAtomsInClass(p.cls);
    if(at.size!==1 || !at.has(k)) continue;
    if(p.cls.indexOf(key)<0) continue;
    const nums=[]; const re=/(-?\d+\.?\d*)[,\s]+(-?\d+\.?\d*)/g; let m;
    while((m=re.exec(p.d))) nums.push([parseFloat(m[1]),parseFloat(m[2])]);
    if(nums.length<2) continue;
    let mnx=1e9,mny=1e9,mxx=-1e9,mxy=-1e9;
    for(const q of nums){ if(q[0]<mnx)mnx=q[0]; if(q[0]>mxx)mxx=q[0]; if(q[1]<mny)mny=q[1]; if(q[1]>mxy)mxy=q[1]; }
    const gc=[(mnx+mxx)/2,(mny+mxy)/2];
    const ax=_prpAxisPos(paths, k);                 // 优先用键轴交点(更准); 与字形中心差太多则不信
    if(ax && Math.hypot(ax[0]-gc[0], ax[1]-gc[1])<6) return ax;
    return gc;
  }
  return null;
}
// 精确原子坐标: 带标签的杂原子 → 字形中心; 无标签的碳 → _prpGetAtomPos(碳原子无截短, 各键在顶点精确共点)。
// 说明: 不用"平行双线平均/直线求交"——双键的内侧平行线是自带侧向偏移的装饰线, 取平均会偏离真实顶点。
function _prpAtomExact(paths, k){
  const lp=_prpLabelPos(paths, k);
  if(lp) return lp;
  return _prpGetAtomPos(paths, k);
}

function postRotatePendantRings(svg, blk){
  if(!svg || !blk) return svg;
  const P=_prpParseV2000(blk); const adj=_prpBuildAdj(P.bondLines,P.nAtoms);
  const rings=[];
  for(const t of [3,4,5,6,7,8]){ for(const r of _prpFindRings(adj,P.nAtoms,t)) rings.push(r); }
  const atomRings=Array.from({length:P.nAtoms},()=>[]);
  rings.forEach((r,i)=>{ for(const a of r) atomRings[a].push(i); });
  const sharesAtom=(ri)=>{ for(const a of rings[ri]) if(atomRings[a].length>1) return true; return false; };
  // ★ 2026-11-05(嵌套悬挂环): 取代基上限 5 -> 12 —— "2-苄氧基/2-苯硫基苯基"等联芳基醚类的取代基
  //   片段含 7~8 个原子(比 5 大), 之前被叶环判定直接排除, 连工序都不进, 环歪 5~15° 无人处理
  //   (实测残留 4 环全部如此)。放宽到 12 后, 取代基里若还含另一悬挂环(如苄氧基的外苯环),
  //   由下方"嵌套作业"机制保证内外两环都正确处理: 元素归最内层作业, rotate 组嵌套复合。
  //   两端都是大片段的"桥连环"仍会被排除(次大分量远超 12)。
  const SUBST_MAX=12;
  const isLeaf=(ringAtoms)=>{ const rem=new Array(P.nAtoms).fill(true); for(const a of ringAtoms) rem[a]=false;
    const rest=[]; for(let i=0;i<P.nAtoms;i++) if(rem[i]) rest.push(i); if(!rest.length) return false;
    const visited=new Set(); const comps=[];
    for(const s0 of rest){ if(visited.has(s0)) continue; const comp=[]; const q=[s0]; visited.add(s0);
      while(q.length){ const u=q.pop(); comp.push(u); for(const nb of adj[u]){ if(rem[nb]&&!visited.has(nb)){ visited.add(nb); q.push(nb); } } }
      comps.push(comp); }
    if(!comps.length) return false;
    comps.sort((a,b)=>b.length-a.length);
    for(let k=1;k<comps.length;k++){ if(comps[k].length>SUBST_MAX) return false; }
    return true; };
  const reachExcl=(nb, ca)=>{ const seen=new Set([nb]); const st=[nb]; let cnt=0;
    while(st.length){ const u=st.pop(); cnt++; for(const v of adj[u]){ if(seen.has(v))continue; if(ca.has(v))continue; seen.add(v); st.push(v); } } return cnt; };
  const jobs=[];
  for(let ri=0; ri<rings.length; ri++){
    if(sharesAtom(ri)) continue;
    if(!isLeaf(rings[ri])) continue;
    const ringAtoms=new Set(rings[ri]);
    let tether=null, best=-1;
    for(const a of rings[ri]){ for(const nb of adj[a]){ if(ringAtoms.has(nb)) continue; const sz=reachExcl(nb, ringAtoms); if(sz>best){ best=sz; tether={ip:a, anchor:nb}; } } }
    if(!tether) continue;
    const frag=new Set(ringAtoms); const q=[...ringAtoms];
    while(q.length){ const u=q.pop(); for(const nb of adj[u]){
      if(frag.has(nb)) continue;
      if((u===tether.ip&&nb===tether.anchor)||(u===tether.anchor&&nb===tether.ip)) continue;
      frag.add(nb); q.push(nb); } }
    if(frag.size > Math.max(8, Math.floor(P.nAtoms*0.6))) continue;
    // 环的循环序(ipso 转到首位): rings[ri] 相邻元素成环(含首尾闭合)
    const cycIdx=rings[ri].indexOf(tether.ip);
    const cyc=rings[ri].slice(cycIdx).concat(rings[ri].slice(0,cycIdx));
    jobs.push({ringIdx:ri, ringAtoms, tether, frag, cyc});
  }
  if(!jobs.length) return svg;
  const paths0=_prpGetPaths(svg);
  const ap={}; const getPos=(a)=>{ if(ap[a]!==undefined) return ap[a]; const p=_prpGetAtomPos(paths0,a); ap[a]=p||null; return ap[a]; };
  // 全部原子的精确 SVG 位置(块坐标经稳健相似变换预测); 后续所有几何量都用它, 不再逐个"猜顶点"。
  const AP=_prpPositions(paths0, P);

  // ── 预处理: 内部不规则(被压扁/拉伸)的悬挂环 —— 刚体旋转不改变环顶点相对几何, 修不了这种歪斜,
  //    直接在 SVG 上把环重投影为正 n 边形(系链为对称轴, ipso 原位不动, 环心落在系链反侧):
  //    两侧夹角精确 = 90°+180°/n(苯环 120°/120°)。已达标(误差<2°)的环不动, 留给下方刚体旋转。──
  const dd=(x)=>{ let d=((x)%360+540)%360-180; return d; };
  const area2=(pts)=>{ let A=0; for(let k=0;k<pts.length;k++){ const p=pts[k], q2=pts[(k+1)%pts.length]; A+=p[0]*q2[1]-q2[0]*p[1]; } return A; };
  // ★ 嵌套作业(2026-11-05): 悬挂环的取代基本身可含另一悬挂环(2-苄氧基/苯硫基苯基等)。
  //   环内重投影位移(ringD, 作用于环顶点/相似变换)与外层作业对取代基的刚体携带位移(carryD)分开记账,
  //   totalD=两者之和 —— 内外作业谁先处理都不影响每原子总位移; d 改写单遍施加 totalD。
  const ringD=new Array(P.nAtoms).fill(null);
  const carryD=new Array(P.nAtoms).fill(null);
  const totalD=(a)=>{ const r=ringD[a], c=carryD[a]; if(r&&c) return [r[0]+c[0], r[1]+c[1]]; return r||c; };
  const doneJobs=[];
  for(const j of jobs){
    const n=j.cyc.length;
    const {ip, anchor}=j.tether;
    const ipP=AP[ip] || getPos(ip);
    if(!ipP) continue;
    const pos={}; let ok=true;
    for(const a of j.cyc){ const p=AP[a] || getPos(a); if(!p){ ok=false; break; } pos[a]=p; }
    if(!ok) continue;
    const o1=j.cyc[1], o2=j.cyc[n-1];
    let tA=_prpBondDir(paths0, ip, anchor, ipP);
    if(tA==null){ const anP=getPos(anchor); if(!anP) continue; tA=Math.atan2(anP[1]-ipP[1], anP[0]-ipP[0]); }
    // ★ 2026-11-05: 环键方向优先取"画线键轴"(截短沿键轴故方向精确, 免疫杂原子标签字形偏差),
    //   这才是用户看到的图形; 取不到(无该键路径)才退回预测位置向量。
    const _dirOf=(nb, fb)=>{ const t=_prpBondDir(paths0, ip, nb, ipP); if(t!=null) return t;
      return Math.atan2(fb[1]-ipP[1], fb[0]-ipP[0]); };
    const b1=_dirOf(o1, pos[o1]);
    const b2=_dirOf(o2, pos[o2]);
    const a1=Math.abs(dd((b1-tA)*180/Math.PI)), a2=Math.abs(dd((b2-tA)*180/Math.PI));
    const target=90+180/n;                       // 正 n 边形悬挂环: 系链两侧夹角
    if(Math.abs(a1-target)<2 && Math.abs(a2-target)<2) continue;   // 已达标, 不动
    // 环平均键长 / 环半径 / 理想顶点(环心在系链反侧 -> 环背向母核张开)
    let s=0; for(let k=0;k<n;k++){ const p=pos[j.cyc[k]], q2=pos[j.cyc[(k+1)%n]]; s+=Math.hypot(q2[0]-p[0],q2[1]-p[1]); } s/=n;
    if(s<1e-6) continue;
    const R=s/(2*Math.sin(Math.PI/n));
    const C=[ipP[0]-Math.cos(tA)*R, ipP[1]-Math.sin(tA)*R];
    let V=[]; for(let k=0;k<n;k++){ const a=tA+2*Math.PI*k/n; V.push([C[0]+R*Math.cos(a), C[1]+R*Math.sin(a)]); }
    // 防镜像: 旧环循环序与理想顶点序的 shoelace 符号若相反, 只翻 cyc 尾序(ipso 保持首位), V 不动。
    // ★ 坑(2026-11-05): cyc 与 V 若"同步翻转"则两次翻转相消 = 完全没翻, 环原子拿到镜像侧顶点
    //   (实测: 9-(2-甲基苯基)吖啶二酮 atom2 应去 V5(old2≈V5 实测重合), 同步翻转后去了 V1,
    //   位移 48px, 环被拉成自交形, ipso 夹角塌成 4°/116°)。仅翻 cyc 才是正确的镜像校正。
    //   (此前自交星形另有元凶: _sim 虚部符号反了(带镜像的相似变换), 已单独修正, 与本翻转无关。)
    const oldPts=j.cyc.map(a=>pos[a]);
    if(area2(oldPts)*area2(V)<0){
      j.cyc=[j.cyc[0]].concat(j.cyc.slice(1).reverse());
    }
    // 逐原子环内重投影增量 ringD(环原子=新顶点-旧位置; 不与外层携带位移混账)
    for(let k=0;k<n;k++){ const a=j.cyc[k]; ringD[a]=[V[k][0]-pos[a][0], V[k][1]-pos[a][1]]; }
    let maxd=0; for(const a of j.cyc) maxd=Math.max(maxd, Math.hypot(ringD[a][0],ringD[a][1]));
    if(maxd>2.2*s) continue;                     // 增量异常: 放弃重投影(留待旋转兜底)
    // 悬挂取代基随挂靠环原子整体平移(BFS 自环原子向片断传播):
    //   子原子刚体位移 = 挂靠原子总位移(ringD+carryD); 若取代基内还有被内层作业重投影过的原子,
    //   totalD 在施加时自动叠加, 顺序无关。
    { const seen=new Set(j.cyc); const q2=[...j.cyc];
      while(q2.length){ const x=q2.shift(); for(const nb of adj[x]){
        if(j.frag.has(nb)&&!seen.has(nb)){ seen.add(nb); const t=totalD(x); carryD[nb]= t? t.slice() : [0,0]; q2.push(nb); } } } }
    // 碰撞核验: 新顶点不得贴近母核原子
    let clash=false;
    outer: for(let k=0;k<n;k++){
      for(let ci=0; ci<P.nAtoms; ci++){
        if(j.frag.has(ci)) continue;
        const cp=getPos(ci); if(!cp) continue;
        if(Math.hypot(V[k][0]-cp[0], V[k][1]-cp[1])<0.55*s){ clash=true; break outer; }
      }
    }
    if(clash) continue;
    j.done=true; j.th=null; j.posX=pos;           // 已处理: 下方刚体旋转跳过
    doneJobs.push(j);
  }
  if(doneJobs.length){
    // 改写这些环的路径坐标: 环内键(含双键内侧线, 与主线同 class)按(旧a,旧b)->(新a,新b)相似变换;
    // 其余片断图元/杂原子标签字形按 Δ 平移; 无 class 小黑点三角按质心就近随 Δ 平移(否则残留旧顶点成游离黑点)。
    const _fmt=(v)=>{ const r2=Math.round(v*100)/100; return (Math.abs(r2)<1e-9?0:r2).toFixed(2); };
    const _sim=(oa,ob,na,nb)=>{ const vx=ob[0]-oa[0], vy=ob[1]-oa[1]; const L2=vx*vx+vy*vy;
      if(L2<1e-9) return null;
      const dx=nb[0]-na[0], dy=nb[1]-na[1];
      // 相似变换 w = d·conj(v)/|v|² : 实部=(d·v)/|v|², 虚部=(dy·vx - dx·vy)/|v|²。
      // ★ 坑(2026-11-05): 虚部原写作 (dx·vy - dy·vx)/|v|², 符号反了 → 该变换**带镜像**(定向反转),
      //   于是 ob 并不映射到 nb, 相邻两条键对同一顶点给出不同坐标(实测错开约 4px),
      //   环被拆成自交星形、内部角变成 90/90/120/120/90/90。
      const sc=(dx*vx+dy*vy)/L2;
      const ss=(dy*vx-dx*vy)/L2;
      return (p)=>{ const ex=p[0]-oa[0], ey=p[1]-oa[1];
        return [na[0]+sc*ex-ss*ey, na[1]+ss*ex+sc*ey]; }; };
    const _rwD=(d, xf)=>{ const nums=[]; d.replace(/-?\d+\.?\d*/g,t=>{nums.push(parseFloat(t)); return t;});
      if(nums.length<2 || nums.length%2!==0) return null;
      for(let k=0;k<nums.length;k+=2){ const p=xf([nums[k],nums[k+1]]); nums[k]=p[0]; nums[k+1]=p[1]; }
      let idx=0; return d.replace(/-?\d+\.?\d*/g, ()=>_fmt(nums[idx++])); };
    svg = svg.replace(/<path\b[^>]*\/?>/g, (tag)=>{
      const cls=(tag.match(/class='([^']*)'/)||[])[1]||'';
      const d=(tag.match(/d='([^']*)'/)||[])[1]||'';
      if(!d) return tag;
      const atoms=[]; { let mm; const reA=/atom-(\d+)/g; while((mm=reA.exec(cls))) atoms.push(+mm[1]); }
      const uniq=[...new Set(atoms)];
      let xf=null;
      if(!uniq.length){
        // 无 class 小黑点三角: 质心就近挂在某带 Δ 的片断原子上随其平移
        if(!/fill:none/.test((tag.match(/style='([^']*)'/)||[])[1]||'')) return tag;
        const pts=(d.match(/(-?\d+\.?\d*),(-?\d+\.?\d*)/g)||[]);
        if(pts.length<3) return tag;
        const nums0=pts.map(t=>t.split(',').map(Number));
        let cxm=0,cym=0; for(const q2 of nums0){ cxm+=q2[0]; cym+=q2[1]; } cxm/=nums0.length; cym/=nums0.length;
        let bestA=null, bestD=8;
        for(let a=0;a<P.nAtoms;a++){ if(!totalD(a)) continue; const p=getPos(a); if(!p) continue;
          const dd2=Math.hypot(p[0]-cxm,p[1]-cym); if(dd2<bestD){ bestD=dd2; bestA=a; } }
        if(bestA===null) return tag;
        const t0=totalD(bestA);
        xf=(p)=>[p[0]+t0[0], p[1]+t0[1]];
      } else if(uniq.length>2 || !uniq.some(a=>totalD(a))){
        return tag;                               // 不在任何重投影片断内
      } else if(uniq.length===2 && totalD(uniq[0]) && totalD(uniq[1])){
        // 两原子均带 Δ: 同属某 done 环的环内键 -> 相似变换(目标端点=旧端点+总位移, 外层携带已含其中);
        // 否则(环->取代基/取代基内部键, Δ相同)平移
        let jj=null;
        for(const j2 of doneJobs){ if(j2.ringAtoms.has(uniq[0]) && j2.ringAtoms.has(uniq[1])){ jj=j2; break; } }
        if(jj){
          const oa=jj.posX[uniq[0]], ob=jj.posX[uniq[1]];
          const ta=totalD(uniq[0]), tb=totalD(uniq[1]);
          if(oa&&ob&&ta&&tb){
            xf=_sim(oa, ob, [oa[0]+ta[0], oa[1]+ta[1]],
                         [ob[0]+tb[0], ob[1]+tb[1]]);
          }
        }
        if(!xf){ const a0=uniq[0]; const t1=totalD(a0); xf=(p)=>[p[0]+t1[0], p[1]+t1[1]]; }
      } else {
        const a0=uniq.find(x=>totalD(x));
        const t2=totalD(a0);
        xf=(p)=>[p[0]+t2[0], p[1]+t2[1]];
      }
      if(!xf) return tag;
      const nd=_rwD(d, xf);
      if(nd===null) return tag;
      return tag.replace(/d='[^']*'/, "d='"+nd+"'");
    });
  }

  // 关键修正(2026-11-05): 锚点方向改用 "ipso->锚点 键线方向"(精确) 计算, 不再用锚点估计坐标。
  //   根因: RDKit 为给 N/O 等杂原子标签留位, 会把连向它们的键线"截短"; 用截短端点反推的 N/O 位置系统性内缩
  //         (实测偏 15~22 px), 使 ancAng 偏差约 20°, 把 coordgen 本已"放射外向"的悬挂环额外旋转约 20°,
  //         造成 ipso 处两邻键夹角 142°/98° 的歪斜 —— 即用户所见"与N相连的苯环两个角度不是120°"。
  //         锚点为碳(无标签、不截短)的环不受影响 —— 故"有些库完全正常、有些库普遍歪斜"。
  //   键线被截短的方向沿键轴, 故 bond 方向本身精确; 配合环心/ipso(均为环碳, 不截短)即得 120°/120°。
  for(const j of jobs){
    if(j.done) continue;                          // 已在上方重投影处理
    const {ip, anchor}=j.tether;
    const ipP=AP[ip] || getPos(ip); if(!ipP){ j.th=null; continue; }
    // ★ 嵌套作业(2026-11-05): 若本环的 ipso 被外层 done 作业的取代基刚体携带位移(carryD)挪过,
    //   旋转轴心必须用挪动后的位置(元素坐标已在 d 改写时加了同一位移, 轴心不加就错位)。
    const dpc=totalD(ip);
    j.ipPos= dpc? [ipP[0]+dpc[0], ipP[1]+dpc[1]] : ipP;      // 旋转轴心(精确, 含外层携带)
    // 放射方向(指向环外) = ipso 处两条环键方向的"角平分线"的反向。
    // ★ 2026-11-05: 原来用"环心->ipso", 环心取环原子平均; 环上若有杂原子(N/O/S), 其标签字形中心
    //   有 2~4px 系统偏差, 把环心带偏约 1.6px → 旋转角偏 3.4°, 实测 tether=(123.4,116.5) 而不是 120/120。
    //   角平分线只依赖 ipso 与两个 ortho 三点(均为环碳, 顶点定位精确), 与杂原子位置、环是否正多边形都无关;
    //   正六边形下二者完全等价。取不到三点时才退回环心法。
    let radAng=null;
    { const n=j.cyc.length, o1=j.cyc[1], o2=j.cyc[n-1];
      const pO1=AP[o1]||getPos(o1), pO2=AP[o2]||getPos(o2);
      if(pO1&&pO2){
        const L1=Math.hypot(pO1[0]-ipP[0],pO1[1]-ipP[1]), L2=Math.hypot(pO2[0]-ipP[0],pO2[1]-ipP[1]);
        if(L1>1e-6&&L2>1e-6){
          const vx=-((pO1[0]-ipP[0])/L1+(pO2[0]-ipP[0])/L2);
          const vy=-((pO1[1]-ipP[1])/L1+(pO2[1]-ipP[1])/L2);
          if(Math.hypot(vx,vy)>1e-6) radAng=Math.atan2(vy,vx);
        }
      } }
    if(radAng==null){
      let cx=0,cy=0,nn=0;
      for(const a of j.ringAtoms){ const p=AP[a] || getPos(a); if(p){ cx+=p[0]; cy+=p[1]; nn++; } }
      if(!nn){ j.th=null; continue; } cx/=nn; cy/=nn;
      radAng=Math.atan2(ipP[1]-cy, ipP[0]-cx);                // SVG: 环心->ipso 向外(环原子精确, 含杂环)
    }
    let ancAng=_prpBondDir(paths0, ip, anchor, ipP);          // SVG: ipso->锚点 键线方向(精确, 免疫标签截短)
    if(ancAng==null){ const anP=getPos(anchor); if(!anP){ j.th=null; continue; } ancAng=Math.atan2(anP[1]-ipP[1], anP[0]-ipP[0]); }
    const th0svg=radAng-ancAng;
    // SVG 坐标(y 向下)与块坐标(y 向上)互为镜像; 此处角均在 SVG 空间, 直接用 -th0svg(旧式), 不再反向。
    j.th=-th0svg;
  }
  // 收集所有图元(<path> 与 <text>), 判定归属哪个 job 的片断
  const elems=[]; let m;
  const reP=/<path\b[^>]*>/g; while((m=reP.exec(svg))){ elems.push({start:m.index, end:m.index+m[0].length, full:m[0], cls:(m[0].match(/class='([^']*)'/)||[])[1]||'', owner:-1}); }
  const reT=/<text\b[^>]*>[\s\S]*?<\/text>/g; while((m=reT.exec(svg))){ elems.push({start:m.index, end:m.index+m[0].length, full:m[0], cls:(m[0].match(/class='([^']*)'/)||[])[1]||'', owner:-1}); }
  for(let ei=0; ei<elems.length; ei++){
    const at=_prpAtomsInClass(elems[ei].cls);
    // ★ 嵌套作业(2026-11-05): 元素归属"最内层"的旋转作业(frag 最小者), 外层作业的 rotate 组再嵌套包住内层组,
    //   两个旋转绕各自的 ipso(原始坐标)复合 = 先修内环取向、再整体随外环转, 数学上严格正确。
    //   "自己的环"图元(两原子同属某作业的 ringAtoms)优先, 防止环绕核的怪异 frag 包含把环错交给别的作业。
    //   done(重投影)作业不认领 —— 其元素已在 d 改写中原位更新, 但若被外层旋转作业的 frag 包含,
    //   会由外层认领、随外层组一起旋转(保持其重投影后的形状)。
    let bestRing=-1, bestFrag=-1;
    for(let ji=0; ji<jobs.length; ji++){ const j=jobs[ji]; if(j.th===null) continue;
      let inf=0; for(const a of at){ if(j.frag.has(a)) inf++; }
      if(!((at.size===2&&inf===2)||(at.size===1&&inf===1))) continue;
      const arr=[...at];
      const ringOwn = (at.size===2) ? (j.ringAtoms.has(arr[0])&&j.ringAtoms.has(arr[1])) : j.ringAtoms.has(arr[0]);
      if(ringOwn){ if(bestRing<0 || j.frag.size<jobs[bestRing].frag.size) bestRing=ji; }
      else { if(bestFrag<0 || j.frag.size<jobs[bestFrag].frag.size) bestFrag=ji; }
    }
    elems[ei].owner = bestRing>=0 ? bestRing : bestFrag;
  }
  // 携带 RDKit 在无 class 碳原子上画的"小黑点"三角(填充:none、黑色描边、仅 3 点):
  // 这类图元不带 atom- class, 上面按 class 的归属逻辑会漏掉它们; 若不随片断一起旋转,
  // 旋转后它们会残留在旧顶点位置, 形成"游离黑点/环变形"假象(小分子尤为明显, 因点较大未被 stripStrayDots 剔除)。
  // 故按质心位置判断: 若落在某 job 片断原子的容差范围内, 则归入该 job 一并旋转。
  for(const e of elems){ if(e.owner>=0) continue; if(e.cls!=='') continue;   // 仅未分类图元
    const style=(e.full.match(/style='([^']*)'/)||[])[1]||'';
    if(!/fill:none/.test(style)) continue;
    if(!/stroke:\s*#000000|stroke:\s*#000\b|stroke:\s*black/i.test(style)) continue;
    const d=(e.full.match(/d='([^']*)'/)||[])[1]||'';
    const pts=(d.match(/(-?\d+\.?\d*),(-?\d+\.?\d*)/g)||[]);
    if(pts.length!==3) continue;
    const nums=pts.map(s=>s.split(',').map(Number));
    const cxm=(nums[0][0]+nums[1][0]+nums[2][0])/3, cym=(nums[0][1]+nums[1][1]+nums[2][1])/3;
    let bestJi=-1;                                     // 小黑点也归最内层(与上面归属口径一致)
    for(let ji=0; ji<jobs.length; ji++){ const j=jobs[ji]; if(j.th===null) continue;
      let hit=false;
      for(const a of j.frag){ const p=getPos(a); if(p && Math.hypot(p[0]-cxm, p[1]-cym) < 7){ hit=true; break; } }
      if(hit && (bestJi<0 || j.frag.size<jobs[bestJi].frag.size)) bestJi=ji;
    }
    if(bestJi>=0) e.owner=bestJi;
  }
  // 删除片断图元(按 start 降序, 保持索引有效), 再按 job 包进 <g transform="rotate(...)">
  const removeSet=elems.filter(e=>e.owner>=0);
  removeSet.sort((a,b)=>b.start-a.start);
  for(const e of removeSet){ svg = svg.slice(0,e.start) + svg.slice(e.end); }
  const groups=[];
  { // ★ 嵌套作业(2026-11-05): 若作业 B 的 frag ⊂ 作业 A 的 frag(如苄氧基外苯环 ⊂ 内苯环作业),
    //   B 的 rotate 组嵌套进 A 的组内 —— 两个旋转都绕原始坐标轴心, 复合 = 先内后外, 严格正确。
    const rotIdx=[]; for(let ji=0; ji<jobs.length; ji++){ if(jobs[ji].th!==null) rotIdx.push(ji); }
    const fragSubset=(big,small)=>{ if(big.size<=small.size) return false; for(const a of small){ if(!big.has(a)) return false; } return true; };
    const kids=new Map();
    for(const a of rotIdx) for(const b of rotIdx){
      if(a!==b && fragSubset(jobs[a].frag, jobs[b].frag)){ if(!kids.has(a)) kids.set(a,[]); kids.get(a).push(b); } }
    const isKid=new Set(); for(const v of kids.values()) for(const b of v) isKid.add(b);
    const built=new Map();
    const buildGroup=(ji)=>{
      if(built.has(ji)) return built.get(ji);
      built.set(ji, null);                            // 防环
      const j=jobs[ji];
      const ipP=j.ipPos || getPos(j.tether.ip); if(!ipP) return null;
      const deg=j.th*180/Math.PI;
      let inner=removeSet.filter(e=>e.owner===ji).map(e=>e.full).join('');
      for(const b of (kids.get(ji)||[])){ const g=buildGroup(b); if(g) inner+=g; }
      if(!inner) return null;
      const g='<g transform="rotate('+deg.toFixed(3)+', '+ipP[0].toFixed(2)+', '+ipP[1].toFixed(2)+')">'+inner+'</g>';
      built.set(ji, g); return g;
    };
    for(const ji of rotIdx){ if(isKid.has(ji)) continue; const g=buildGroup(ji); if(g) groups.push(g); }
  }
  if(groups.length){
    const ci=svg.lastIndexOf('</svg>');
    if(ci>=0) svg = svg.slice(0,ci) + groups.join('') + svg.slice(ci);
  }
  // 整体缩放居中回原画布(必要时), 保证不裁切且显示尺寸稳定
  { const _ps=_prpGetPaths(svg); let _mnx=1e9,_mny=1e9,_mxx=-1e9,_myy=-1e9;
    for(const _p of _ps){ const _e=_prpBondEndpoints(_p.d); for(const _q of _e){ if(_q[0]<_mnx)_mnx=_q[0]; if(_q[0]>_mxx)_mxx=_q[0]; if(_q[1]<_mny)_mny=_q[1]; if(_q[1]>_myy)_myy=_q[1]; } }
    const _W=250,_H=195,_mg=10;
    if(!(_mnx>_mxx) && !(_mnx>=_mg && _mxx<=_W-_mg && _mny>=_mg && _myy<=_H-_mg)){
      const _bw=Math.max(_mxx-_mnx,1e-6), _bh=Math.max(_myy-_mny,1e-6);
      const _s=Math.min((_W-2*_mg)/_bw,(_H-2*_mg)/_bh);
      const _tx=(_W-_bw*_s)/2-_mnx*_s, _ty=(_H-_bh*_s)/2-_mny*_s;
      const oi=svg.indexOf('>', svg.indexOf('<svg'));
      const ci2=svg.lastIndexOf('</svg>');
      const inner2=svg.slice(oi+1, ci2);
      svg = svg.slice(0,oi+1) + '<g transform="translate('+_tx.toFixed(2)+','+_ty.toFixed(2)+') scale('+_s.toFixed(4)+')">' + inner2 + '</g>' + svg.slice(ci2);
    }
  }
  return svg;
}

/* 去除 RDKit 默认 SVG 残留在(中性)碳原子上的"小黑点"伪影:
 * 这些是无 class、fill:none、黑色描边、仅 3 个点的极小三角形路径;
 * 其坐标随分子增大而缩小, 但描边宽度固定 2px, 于是塌缩成实心黑点(苯环外散点)。
 * 苯等小分子同类三角形较大(>3.5px)看起来正常→保留; 带电杂原子的电荷标记(N+/O-)通常>3.5px→保留。
 * 仅去除塌缩成黑点的极小三角形, 不影响键线、彩色原子标签与合法电荷标记。 */
function stripStrayDots(svg){
  if(!svg) return svg;
  const re = /<path\b[^>]*\/?>/g; let m; const hits=[];
  while((m=re.exec(svg))){
    const tag = m[0];
    if(/class=/.test(tag)) continue;                       // 保留键线/原子标签(带 class)
    const style = (tag.match(/style='([^']*)'/)||[])[1]||'';
    const d = (tag.match(/d='([^']*)'/)||[])[1]||'';
    if(!/fill:none/.test(style)) continue;                 // 保留填充的彩色原子标签字形
    if(!/stroke:\s*#000000|stroke:\s*#000\b|stroke:\s*black/i.test(style)) continue;
    const pts = (d.match(/(-?\d+\.?\d*),(-?\d+\.?\d*)/g)||[]);
    if(pts.length!==3) continue;                           // 仅 3 点三角形
    const nums = pts.map(s=>s.split(',').map(Number));
    const xs=nums.map(p=>p[0]), ys=nums.map(p=>p[1]);
    const w=Math.max(...xs)-Math.min(...xs), h=Math.max(...ys)-Math.min(...ys);
    if(Math.max(w,h) > 3.5) continue;                      // 仅去除塌缩成黑点的极小三角形
    hits.push({idx:m.index, len:tag.length});
  }
  hits.sort((a,b)=>b.idx-a.idx);
  for(const h of hits){ svg = svg.slice(0,h.idx) + svg.slice(h.idx+h.len); }
  return svg;
}

/* ---------- 结构式绘制 ---------- */
const _SVG_CACHE = new Map();
/* W/H 默认 250×195(检索站卡片); 首页侧栏静态图要求 300×240, 由生成脚本传入。 */
function molToSVG(smiles, _w, _h){
  if (!RDKit || !smiles) return null;
  const W = _w || 250, H = _h || 195;
  const ck = W + 'x' + H + '|' + smiles;
  // 结果缓存: 现在一张图要试 A/E/C 三套坐标 + 三轮打分(约 30ms), 而同一个 SMILES 会被反复渲染
  //   (翻页、改筛选、来回切换分类), 缓存后重复渲染零成本。上限 400 条, 满了清最早的。
  const hit = _SVG_CACHE.get(ck);
  if (hit !== undefined){
    if (hit.blk) globalThis.__PRP_LASTBLK = hit.blk;   // 命中缓存也要把坐标还给审计钩子
    return hit.svg;
  }
  // 盐型 SMILES(含 '.') 优先用整体绘制; 若整体失败, 再按片段尝试, 但片段需按长度降序,
  // 优先绘制主分子而非反离子(如 [H+] / [Cl-]), 否则会只画出 H⁺ 等反离子。
  const parts = smiles.includes('.') ? smiles.split('.') : [smiles];
  const cands = [smiles].concat(parts.slice().sort((a, b) => b.length - a.length));
  for (const c of cands){
    let mol = null;
    try { mol = RDKit.get_mol(c); } catch(e){ mol = null; }
    if (!mol) continue;
    try {
      /* 坐标候选 A: RDKit 自带 depictor(get_new_coords) + 规范取向(straighten_depiction)。
       *   straighten 只是整体刚性旋转 —— 实测逐键长**一个数都不变**(0.797/1.5/1.5 前后完全一致), 环内角不变。 */
      let blkA = null;
      try { blkA = mol.get_new_coords(); } catch(e){ blkA = null; }
      if (blkA){
        try {
          const mt = RDKit.get_mol(blkA);
          if (mt){ try { mt.straighten_depiction(); blkA = mt.get_molblock() || blkA; } finally { mt.delete(); } }
        } catch(e){}
      }
      /* 坐标候选 E: coordgen 坐标(set_new_coords(true))。两套坐标各有胜负 ——
       *   RDKit depictor 偶发把端基键压成半截(如仲丁基吖啶二酮的双羰基一长一短);
       *   coordgen 则偶发把**非端基** C-C 拉到 2.5 倍(4-[2-烷氧基苯基亚甲基氨基]苯甲酸一族)。
       *   故两份都算、逐分子打分取优(_prpPickBlk), 不是固定用哪一个。 */
      let blkE = null;
      try { mol.set_new_coords(true); blkE = mol.get_molblock(); } catch(e){ blkE = null; }
      if (blkA) blkA = _prpLengthenTerminalBonds(blkA);
      if (blkE) blkE = _prpLengthenTerminalBonds(blkE);
      let blk = _prpPickBlk([blkA, blkE]);
      /* 候选 C: 在已选坐标上做"非环键键长松弛", 与 C 再比一次分 —— 专治非端基的过长/过短键。 */
      if (blk){
        let blkC = null;
        try { blkC = _prpRelaxBonds(blk); } catch(e){ blkC = null; }
        if (blkC) blk = _prpPickBlk([blk, blkC]);
      }
      // 调试/审计钩子: 让离线审计脚本能拿到"这一张图实际用的坐标"(浏览器控制台也可查)
      if (blk) globalThis.__PRP_LASTBLK = blk;
      // 关键修正(2026-11-05, C=N 交叉双键): RDKit 写 molblock 时会把"未指定 E/Z 的双键"
      //   标成 stereo=3(either), 再解析回来画图就会画成**交叉双键**(两条线在中点交叉成 X)——
      //   肟 C=N、腙 C=N、查尔酮 C=C 等库全部中招。清掉键的 stereo 字段即可恢复平行双线;
      //   原子手性(楔形键)不在键块这一列, 不受影响。
      if (blk) blk = _prpClearBondStereo(blk);
      let svg0 = null, molFromBlk = null;
      if (blk){ try { molFromBlk = RDKit.get_mol(blk); } catch(e){ molFromBlk = null; } }
      try { svg0 = (molFromBlk || mol).get_svg(W, H); } catch(e){ svg0 = null; }
      if (svg0 == null) svg0 = mol.get_svg(W, H);
      let svg = svg0;
      if (blk){ try { svg = postRotatePendantRings(svg0, blk); } catch(e){} }
      svg = stripStrayDots(svg);
      if (molFromBlk){ try{ molFromBlk.delete(); }catch(e){} }
      mol.delete();
      if (_SVG_CACHE.size >= 400) _SVG_CACHE.delete(_SVG_CACHE.keys().next().value);
      _SVG_CACHE.set(ck, {svg, blk});
      return svg;
    } catch(e){ if (mol) mol.delete(); }
  }
  return null;
}

/* ---------- 数据加载(支持分片) ---------- */
async function load(){
  const loadingText = document.getElementById('loading-text');

  // 优先读取分片清单(products.manifest.json), 回退到单一 products.json(向后兼容)
  let records = [];
  let manifest = null;
  try {
    const mres = await fetch('data/products.manifest.json');
    if (mres.ok) manifest = await mres.json();
  } catch(e){ manifest = null; }

  // 加载统计(含产物库默认展示顺序 library_order)
  let STATS = null;
  try {
    const sres = await fetch('data/stats.json');
    if (sres.ok) STATS = await sres.json();
  } catch(e){ STATS = null; }
  LIBRARY_ORDER = (STATS && Array.isArray(STATS.library_order)) ? STATS.library_order : [];

  if (manifest && Array.isArray(manifest.shards) && manifest.shards.length){
    if (loadingText) loadingText.textContent =
      `正在并行加载 ${manifest.shards.length} 个数据分片，共约 ${(manifest.count/10000).toFixed(1)} 万条…`;
    const parts = await Promise.all(manifest.shards.map(async (name) => {
      const r = await fetch('data/' + name);
      if (!r.ok) throw new Error('分片 ' + name + ' HTTP ' + r.status);
      return r.json();
    }));
    records = [].concat(...parts);
  } else {
    if (loadingText) loadingText.textContent = '正在加载产物库数据…';
    const res = await fetch('data/products.json');
    if (!res.ok) throw new Error('HTTP ' + res.status);
    records = await res.json();
  }

  DATA = records;
  BLOBS = DATA.map(r => [
    r.n, r.e, r.f, r.s, r.g, r.r, r.a, r.d, r.o, r.c,
    r.code,                                  // 产品代码(结构哈希短码)纳入检索: 支持按代码查
    r.x ? Object.values(r.x).join(' ') : ''
  ].filter(Boolean).join(' ').toLowerCase());

  // [结构检索] 重原子数(不含 H)预计算, 用于子结构/完全相同预筛
  ATOM_CNT = DATA.map(r => {
    const f = r.f || ''; let s = 0; let m; const re = /([A-Z][a-z]?)(\d*)/g;
    while ((m = re.exec(f)) !== null){ if (m[1] === 'H') continue; s += (m[2] ? parseInt(m[2], 10) : 1); }
    return s;
  });

  // 分类下拉: 优先按 stats.json 的 library_order 展示(用户指定重排顺序), 未列出的库追加到末尾
  const counts = {};
  DATA.forEach(r => { counts[r.g] = (counts[r.g] || 0) + 1; });
  let cats;
  if (STATS && Array.isArray(STATS.library_order) && STATS.library_order.length){
    const present = new Set(Object.keys(counts));
    cats = STATS.library_order.filter(c => present.has(c));
    Object.keys(counts).forEach(c => { if (!cats.includes(c)) cats.push(c); });
  } else {
    cats = Object.keys(counts).sort((a,b)=>counts[b]-counts[a]);
  }
  const sel = document.getElementById('catFilter');
  cats.forEach(c => {
    const o = document.createElement('option');
    o.value = c; o.textContent = `${c} (${counts[c].toLocaleString()})`;
    sel.appendChild(o);
  });
  document.getElementById('meta').textContent = '共 ' + DATA.length.toLocaleString() + ' 条';
}

/* ---------- 检索 ---------- */
function search(){
  const t0 = performance.now();
  window.SIM_MAP = {};   // 文本检索时清除结构检索遗留的相似度标注
  const q = document.getElementById('search').value.trim().toLowerCase();
  const cat = document.getElementById('catFilter').value;
  const terms = q.split(/\s+/).filter(Boolean);
  const isCodeQuery = q.startsWith('cc-');   // 以 CC- 开头视为产品代码精确/前缀检索
  const out = [];
  for (let i = 0; i < DATA.length; i++){
    if (cat && DATA[i].g !== cat) continue;
    let ok = true;
    if (isCodeQuery){
      const code = (DATA[i].code || '').toLowerCase();
      if (!code.startsWith(q)) ok = false;          // 完整码=精确命中, 部分码=前缀命中
    } else if (terms.length){
      const blob = BLOBS[i];
      for (const t of terms){ if (!blob.includes(t)){ ok = false; break; } }
    }
    if (ok) out.push(i);
  }
  MATCHED = out;
  applySort();
  page = 0;
  const dt = (performance.now() - t0).toFixed(0);
  const modeHint = isCodeQuery ? '（按产品代码精确检索）' : '';
  document.getElementById('status').innerHTML =
    '匹配 <b>' + FILTERED.length.toLocaleString() + '</b> 条 · 用时 ' + dt + ' ms' + modeHint;
  render();
}

/* ---------- 排序(前端, 仅影响显示顺序) ---------- */
function cmpNum(x, y){
  const xn = Number(x), yn = Number(y);
  const ax = isNaN(xn), ay = isNaN(yn);
  if (ax && ay) return 0;
  if (ax) return 1;   // 缺值排末尾
  if (ay) return -1;
  return xn - yn;
}
function applySort(){
  const field = document.getElementById('sortField').value;
  const dir = document.getElementById('sortDir').value === 'desc' ? -1 : 1;
  if (field === 'default'){
    // 默认: 按产物库顺序分组排列(同源 stats.json.library_order)
    if (LIBRARY_ORDER.length){
      const pos = c => { const i = LIBRARY_ORDER.indexOf(c); return i < 0 ? 1e9 : i; };
      FILTERED = MATCHED.slice().sort((a, b) => pos(DATA[a].g) - pos(DATA[b].g));
    } else {
      FILTERED = MATCHED.slice();
    }
    return;
  }
  FILTERED = MATCHED.slice().sort((a, b) => {
    const va = DATA[a], vb = DATA[b];
    let x, y;
    switch (field){
      case 'name':    x = va.n || ''; y = vb.n || ''; return dir * x.localeCompare(y, 'zh');
      case 'name_en': x = va.e || ''; y = vb.e || ''; return dir * x.localeCompare(y, 'zh');
      case 'formula': x = va.f || ''; y = vb.f || ''; return dir * x.localeCompare(y, 'zh');
      case 'cas':     x = va.c || ''; y = vb.c || ''; return dir * x.localeCompare(y, 'zh');
      case 'cat':     x = va.g || ''; y = vb.g || ''; return dir * x.localeCompare(y, 'zh');
      case 'src':     x = va.r || ''; y = vb.r || ''; return dir * x.localeCompare(y, 'zh');
      case 'mw':      return dir * cmpNum(va.m, vb.m);
      default:        return 0;
    }
  });
}
function onSortChange(){
  applySort();
  render();
}

/* ---------- 渲染 ---------- */
function cardHTML(r){
  const f = r.f || '';
  const realName = (r.n && r.n.trim() && r.n.trim() !== f.trim()) ? r.n.trim() : '';
  const cnameHTML = realName ? escapeHTML(realName) : (f ? fmtFormula(f) : '未命名产物');
  const formula = (f && realName) ? `<span class="formula">${fmtFormula(f)}</span>` : '';
  const mw = (r.m != null) ? `<span class="mw">MW ${r.m}</span>` : '';
  const cas = r.c ? `<span class="tag">CAS ${escapeHTML(r.c)}</span>` : '';
  const ename = (r.e && r.e.trim()) ? `<div class="ename">${escapeHTML(r.e.trim())}</div>` : '';
  const sim = (window.SIM_MAP && r.i != null && window.SIM_MAP[r.i] != null) ? window.SIM_MAP[r.i] : null;
  const simHTML = (sim != null) ? `<div class="sim">相似度 ${Math.round(sim * 100)}%</div>` : '';
  const tags = [`<span class="tag">${escapeHTML(r.g)}</span>`];
  if (r.o) tags.push(`<span class="tag">母体 ${escapeHTML(r.o)}</span>`);
  const codeHTML = r.code
    ? `<div class="pcode" title="点击复制产品代码" data-code="${escapeHTML(r.code)}">`
      + `<span class="pk">产品代码</span><code>${escapeHTML(r.code)}</code>`
      + `<span class="pcopy">复制</span></div>`
    : '';
  const reagents = r.a ? `<div class="kv"><span class="k">试剂</span><span>${escapeHTML(r.a)}</span></div>` : '';
  const ccdc = r.d ? `<div class="ccdc">CCDC：${escapeHTML(r.d)}</div>` : '';
  let extra = '';
  if (r.x && Object.keys(r.x).length){
    extra = '<div class="extra">' + Object.entries(r.x).map(([k,v])=>{
      const label = EXTRA_LABEL[k] || k;
      return `<div class="kv"><span class="k">${escapeHTML(label)}</span><span>${escapeHTML(v)}</span></div>`;
    }).join('') + '</div>';
  }
  const src = `<div class="kv"><span class="k">来源</span><span>${escapeHTML(r.r)}</span></div>`;
  return `<div class="card">
    <div class="struct"></div>
    <div class="cname">${cnameHTML}</div>
    ${ename}
    ${simHTML}
    ${codeHTML}
    <div class="crow">${formula}${mw}${cas}</div>
    <div class="tags">${tags.join('')}</div>
    ${reagents}
    ${ccdc}
    ${extra}
    ${src}
  </div>`;
}

function drawStructures(slice){
  const nodes = document.getElementById('grid').querySelectorAll('.struct');
  slice.forEach((idx, k) => {
    const node = nodes[k];
    if (!node) return;
    const r = DATA[idx];
    if (r.s && RDKit){
      const svg = molToSVG(r.s);
      node.innerHTML = svg ? svg : `<div class="noimg">结构渲染失败<br>${fmtFormula(r.f||'')}</div>`;
    } else if (r.s && !RDKit){
      node.innerHTML = '<div class="noimg">结构引擎加载中…</div>';
    } else {
      node.innerHTML = `<div class="noimg">源数据无 SMILES<br>分子式：${fmtFormula(r.f||'—')}</div>`;
    }
  });
}

function render(){
  const grid = document.getElementById('grid');
  const pager = document.getElementById('pager');
  const empty = document.getElementById('empty');
  const total = FILTERED.length;
  if (total === 0){
    grid.innerHTML = '';
    empty.style.display = 'block';
    pager.innerHTML = '';
    return;
  }
  empty.style.display = 'none';
  const pages = Math.ceil(total / PER_PAGE);
  if (page >= pages) page = pages - 1;
  const start = page * PER_PAGE;
  const end = Math.min(start + PER_PAGE, total);
  const slice = FILTERED.slice(start, end);
  grid.innerHTML = slice.map(idx => cardHTML(DATA[idx])).join('');
  drawStructures(slice);
  renderPager(pages, start, end, total);
}

function renderPager(pages, start, end, total){
  const pager = document.getElementById('pager');
  const sizes = [24, 50, 100, 200, 500];
  const selOpts = sizes.map(v =>
    `<option value="${v}"${v === PER_PAGE ? ' selected' : ''}>${v}</option>`).join('');
  let html = `<span class="pgsize-wrap">每页 <select id="pageSize" class="pgsize">${selOpts}</select> 条</span>`;
  html += `<button id="prev"${page<=0?' disabled':''}>上一页</button>`;
  html += `<span class="pginfo">第 ${page+1} / ${pages} 页 · 显示 ${start+1}-${end} / ${total.toLocaleString()}</span>`;
  html += `<button id="next"${page>=pages-1?' disabled':''}>下一页</button>`;
  pager.innerHTML = html;
  const ps = document.getElementById('pageSize');
  if (ps) ps.onchange = () => { PER_PAGE = parseInt(ps.value, 10) || 24; page = 0; render(); };
  const prev = document.getElementById('prev');
  const next = document.getElementById('next');
  if (prev && page > 0) prev.onclick = () => { page--; render(); window.scrollTo({top:0,behavior:'smooth'}); };
  if (next && page < pages-1) next.onclick = () => { page++; render(); window.scrollTo({top:0,behavior:'smooth'}); };
}

/* ---------- 启动 ---------- */
window.addEventListener('DOMContentLoaded', async () => {
  initRDKit();
  try{
    await load();
  }catch(e){
    document.getElementById('loading-text').textContent = '数据加载失败：' + e.message;
    console.error(e);
    return;
  }
  document.getElementById('loading').style.display = 'none';
  document.getElementById('search').addEventListener('input', debounce(search, 180));
  document.getElementById('catFilter').addEventListener('change', search);
  document.getElementById('sortField').addEventListener('change', onSortChange);
  document.getElementById('sortDir').addEventListener('change', onSortChange);
  // 产品代码: 点击卡片上的代码徽标即复制到剪贴板
  const grid = document.getElementById('grid');
  grid.addEventListener('click', (e) => {
    const pc = e.target.closest('.pcode');
    if (!pc) return;
    const code = pc.getAttribute('data-code') || '';
    if (!code) return;
    navigator.clipboard.writeText(code).then(() => {
      pc.classList.add('copied');
      const lbl = pc.querySelector('.pcopy');
      if (lbl) lbl.textContent = '已复制';
      setTimeout(() => { pc.classList.remove('copied'); if (lbl) lbl.textContent = '复制'; }, 1100);
    }).catch(() => {});
  });
  search();
});
