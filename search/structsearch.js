/* 结构检索（子结构 / 完全相同 / 相似度 Top-N）—— 基于站内已有 RDKit WASM, 纯前端、无服务器
 * 依赖 app.js 提供的全局: DATA, RDKit, MATCHED, FILTERED, page, PER_PAGE, applySort(), render()
 * 与 ATOM_CNT（每条记录重原子数）。
 *
 * 性能约定:
 *  - 所有长循环按时间切片: 每 ~24ms 让出主线程一次并刷新进度条, 页面不卡死;
 *  - 候选分子对象用完立即 delete(), 避免 WASM 堆内存无限增长(此前假死的根因之一);
 *  - 相似度用 2048 位查表(Uint8Array)代替逐位二分, 快一个数量级;
 *  - 运行中"运行结构检索"按钮变为"取消检索", 可随时中止。
 */
'use strict';

let FP_READY = false;
let FP_N = 0;
let FP_OFFSETS = null;     // Uint32Array(N+1)
let FP_BITS = null;        // Uint16Array(所有位, 全局拼接)
let KETCHER_READY = false;
let ketcherFrame = null;
let STRUCT_RUNNING = false;
let STRUCT_ABORT = false;
// window.SIM_MAP: 记录 id -> 相似度(0..1), 仅相似度模式用于卡片展示

const ABORT_ERR = '__STRUCT_ABORT__';

const FORMULA_RE = /([A-Z][a-z]?)(\d*)/g;
function parseFormulaElements(f){
  const m = {};
  if (!f) return m;
  FORMULA_RE.lastIndex = 0;
  let mm;
  while ((mm = FORMULA_RE.exec(f)) !== null){
    const el = mm[1];
    if (el === 'H') continue;
    const c = mm[2] ? parseInt(mm[2], 10) : 1;
    m[el] = (m[el] || 0) + c;
  }
  return m;
}
// 查询 SMILES/SMARTS -> {heavy, elem}（best-effort, 仅用于预筛; 正确性由 RDKit 保证）
// 通用解析元素符号, 覆盖 As/Se/Si 等杂元素(画板可能绘制含 As 的胂酸)
function parseStructElements(s){
  const elem = {};
  let heavy = 0;
  if (!s) return {heavy, elem};
  const n = s.length;
  let i = 0;
  while (i < n){
    const ch = s[i];
    if (ch === '['){
      const end = s.indexOf(']', i);
      if (end < 0) break;
      const inner = s.slice(i + 1, end);
      const em = inner.match(/^([A-Z][a-z]?)/);
      const el = em ? em[1] : '';
      if (el && el !== 'H'){ elem[el] = (elem[el] || 0) + 1; heavy++; }
      i = end + 1;
    } else if (ch === '#'){
      const em = s.slice(i).match(/^#(\d+)/);
      if (em){
        const an = parseInt(em[1], 10);
        const map = {5:'B',6:'C',7:'N',8:'O',9:'F',14:'Si',15:'P',16:'S',17:'Cl',33:'As',34:'Se',35:'Br',53:'I'};
        const el = map[an];
        if (el){ elem[el] = (elem[el] || 0) + 1; heavy++; }
        i += em[0].length; continue;
      }
      i++;
    } else {
      const two = s.slice(i, i + 2);
      let el = null;
      if (two === 'Cl' || two === 'Br' || two === 'As' || two === 'Se' || two === 'Si'){ el = two; i += 2; }
      else if ('BCNOFIPSbcnops'.includes(ch)){ el = ch.toUpperCase(); i += 1; }
      else { i += 1; continue; }
      if (el && el !== 'H'){ elem[el] = (elem[el] || 0) + 1; heavy++; }
    }
  }
  return {heavy, elem};
}
function elemSuperset(f, QE){
  if (!QE) return true;
  const m = parseFormulaElements(f);
  for (const e in QE){ if ((m[e] || 0) < QE[e]) return false; }
  return true;
}

/* ---------- fp.bin 加载（相似度模式按需懒加载） ---------- */
// 文件格式(与 aggregate.py 写出一致): Uint32 numRecords; 之后每条记录:
//   Uint16 count + count×Uint16 位索引(升序) —— 交错紧凑排列, 无独立 count 表
async function loadFP(){
  if (FP_READY) return;
  setStructStatus('正在加载指纹库 fp.bin（约 80MB，首次较慢）…');
  const t0 = performance.now();
  const res = await fetch('data/fp.bin', {cache: 'no-cache'});
  const buf = await res.arrayBuffer();
  const dv = new DataView(buf);
  FP_N = dv.getUint32(0, true);
  // 第一遍: 走一遍记录边界, 建每条记录的位区间 offsets, 并统计总位数
  FP_OFFSETS = new Uint32Array(FP_N + 1);
  let total = 0, p = 4;
  for (let i = 0; i < FP_N; i++){
    FP_OFFSETS[i] = total;
    const c = dv.getUint16(p, true);
    p += 2 + c * 2;
    total += c;
  }
  FP_OFFSETS[FP_N] = total;
  if (p !== buf.byteLength){
    console.warn('fp.bin 尾部校验不符: p=' + p + ' size=' + buf.byteLength + '（继续但结果可能不完整）');
  }
  // 第二遍: 填充位数据
  FP_BITS = new Uint16Array(total);
  p = 4;
  let last = performance.now();
  for (let i = 0; i < FP_N; i++){
    const c = dv.getUint16(p, true); p += 2;
    const base = FP_OFFSETS[i];
    for (let j = 0; j < c; j++){ FP_BITS[base + j] = dv.getUint16(p, true); p += 2; }
    if ((i & 65535) === 65535){
      if (STRUCT_ABORT) throw new Error(ABORT_ERR);
      const now = performance.now();
      if (now - last >= 24){ await new Promise(r => setTimeout(r, 0)); last = performance.now(); }
    }
  }
  FP_READY = true;
  console.log('fp.bin 加载', (performance.now() - t0).toFixed(0), 'ms, 记录', FP_N, ', 位总数', total);
}

/* ---------- SMILES 归一化 & SMARTS 识别 ---------- */
// 画板 Ketcher(Indigo) 输出的 SMILES 常把苯环写成 Kekulé 形式(大写 C + 显式双键),
// 而目标库 SMILES 由 RDKit 生成(芳香小写 c). 直接把 Kekulé SMILES 交给 RDKit.get_qmol
// 会保留非芳香标志, 导致子结构匹配对芳香环 0 命中. 这里统一经 RDKit 重新感知芳香性.
function normalizeSmiles(smi){
  if (!smi) return smi;
  let m = null;
  try { m = RDKit.get_mol(smi); } catch(e){ m = null; }
  if (!m) return smi;            // 解析失败则原样返回, 交由上层报错
  const out = m.get_smiles();
  m.delete();
  return (out && out.trim()) ? out.trim() : smi;
}
// SMARTS 含通配/逻辑/原子类约束, 不应被当作 SMILES 归一化(会破坏语义); 仅 SMILES 才归一化
function looksLikeSMARTS(s){
  return /[!&;*]/.test(s) || /\[#/.test(s) || /\$\(/.test(s);
}

/* ---------- UI 辅助 ---------- */
function SIM_MAP_RESET(){ if (window.SIM_MAP) window.SIM_MAP = {}; }
function setStructStatus(t){ const e = document.getElementById('structStatus'); if (e) e.textContent = t; }
function setStructProgress(f){
  const b = document.getElementById('structProgressBar');
  if (b) b.style.width = Math.max(0, Math.min(100, Math.round(f * 100))) + '%';
}
function setRunButtonRunning(on){
  const b = document.getElementById('btnStructRun');
  if (!b) return;
  b.textContent = on ? '取消检索' : '运行结构检索';
}
// 预筛候选(子结构: 重原子数>=QH; 完全相同: ===QH), 时间切片防卡顿
async function collectCands(mode, QH, QE){
  const cands = [];
  let last = performance.now();
  for (let i = 0; i < DATA.length; i++){
    if (!DATA[i].s) continue;
    if (mode === 'exact'){ if (ATOM_CNT[i] !== QH) continue; }
    else if (ATOM_CNT[i] < QH) continue;
    if (!elemSuperset(DATA[i].f, QE)) continue;
    cands.push(i);
    if ((i & 4095) === 4095){
      if (STRUCT_ABORT) throw new Error(ABORT_ERR);
      const now = performance.now();
      if (now - last >= 24){ setStructProgress(i / DATA.length); await new Promise(r => setTimeout(r, 0)); last = performance.now(); }
    }
  }
  return cands;
}

/* ---------- 主入口 ---------- */
async function runStructureSearch(){
  if (!RDKit){ alert('化学结构引擎尚未就绪，请稍候重试'); return; }
  if (STRUCT_RUNNING){ STRUCT_ABORT = true; setStructStatus('正在取消…'); return; }
  const mode = document.querySelector('input[name="structMode"]:checked').value;
  const raw = await getStructQuery();
  if (!raw){ alert('请先在画板绘制，或在文本框粘贴 SMILES / SMARTS'); return; }
  const isSmarts = looksLikeSMARTS(raw);
  const input = isSmarts ? raw : normalizeSmiles(raw);   // SMILES 经 RDKit 归一化(统一芳香性)
  const t0 = performance.now();
  setStructStatus('正在解析结构…');
  SIM_MAP_RESET();
  let qmol = null, qfpSet = null, QH = 0, QE = null;
  try {
    if (mode === 'similar'){
      if (isSmarts){ alert('相似度模式请使用 SMILES（不支持 SMARTS）'); return; }
      qfpSet = smilesToFpSet(input);
      if (!qfpSet){ alert('无法将该输入解析为分子（相似度请用 SMILES）'); return; }
    } else {
      qmol = RDKit.get_qmol(input);
      QH = qmol.get_num_atoms();
      QE = parseStructElements(input).elem;
      if (!QH){ alert('无法解析该结构'); return; }
    }
  } catch(e){ alert('结构解析失败：' + e.message); return; }

  STRUCT_RUNNING = true; STRUCT_ABORT = false; setRunButtonRunning(true);
  const matches = [];
  let aborted = false;
  try {
    if (mode === 'sub' || mode === 'exact'){
      setStructStatus('预筛中…');
      const cands = await collectCands(mode, QH, QE);
      const exact = (mode === 'exact');
      setStructStatus((exact ? '精确匹配' : '子结构匹配') + '中… 候选 ' + cands.length.toLocaleString() + ' 条（可点"取消检索"中止）');
      let last = performance.now();
      for (let k = 0; k < cands.length; k++){
        const i = cands[k];
        let mol = null; try { mol = RDKit.get_mol(DATA[i].s); } catch(e){ mol = null; }
        if (mol){
          let res = null;
          try { res = JSON.parse(mol.get_substruct_match(qmol)); } catch(e){ res = null; }
          const n = (res && res.atoms) ? res.atoms.length : 0;
          if (exact ? (n === QH) : (n > 0)) matches.push(i);
          try { mol.delete(); } catch(e){}          // 用完即删, 防 WASM 堆膨胀
        }
        if ((k & 255) === 255){
          if (STRUCT_ABORT) throw new Error(ABORT_ERR);
          const now = performance.now();
          if (now - last >= 24){
            setStructProgress(k / cands.length);
            setStructStatus((exact ? '精确匹配' : '子结构匹配') + '中… ' + (k + 1).toLocaleString() + '/' + cands.length.toLocaleString() + ' · 已匹配 ' + matches.length.toLocaleString());
            await new Promise(r => setTimeout(r, 0));
            last = performance.now();
          }
        }
      }
    } else if (mode === 'similar'){
      await loadFP();
      setStructStatus('相似度计算中…');
      const topN = Math.min(500, Math.max(1, parseInt(document.getElementById('structTopN').value, 10) || 50));
      // 查询指纹 -> 2048 位查表(比逐位二分快一个数量级)
      const QB = new Uint8Array(2048);
      let QBITS = 0;
      qfpSet.forEach(b => { if (b >= 0 && b < 2048){ QB[b] = 1; QBITS++; } });
      const idx = [], sims = [];
      const N = Math.min(FP_N, DATA.length);
      let last = performance.now();
      for (let i = 0; i < N; i++){
        const s = FP_OFFSETS[i], e = FP_OFFSETS[i + 1];
        if (e > s){
          let inter = 0;
          for (let p = s; p < e; p++) inter += QB[FP_BITS[p]];
          const denom = QBITS + (e - s) - inter;
          if (denom > 0){
            const t = inter / denom;
            if (t > 0){ idx.push(i); sims.push(t); }
          }
        }
        if ((i & 8191) === 8191){
          if (STRUCT_ABORT) throw new Error(ABORT_ERR);
          const now = performance.now();
          if (now - last >= 24){
            setStructProgress(i / N);
            await new Promise(r => setTimeout(r, 0));
            last = performance.now();
          }
        }
      }
      const order = idx.map((v, k) => k).sort((a, b) => sims[b] - sims[a]);
      for (let k = 0; k < order.length && k < topN; k++){
        const i = idx[order[k]];
        matches.push(i);
        if (window.SIM_MAP) window.SIM_MAP[DATA[i].i] = sims[order[k]];
      }
    }
  } catch(e){
    if (e && e.message === ABORT_ERR){ aborted = true; }
    else { alert('检索出错：' + e.message); }
  } finally {
    STRUCT_RUNNING = false; STRUCT_ABORT = false;
    setRunButtonRunning(false);
    setStructProgress(0);
  }

  if (qmol){ try { qmol.delete(); } catch(e){} }
  if (aborted){
    setStructStatus('已取消');
    document.getElementById('status').innerHTML = '结构检索已取消';
    return;
  }
  MATCHED = matches;
  if (mode === 'similar'){ FILTERED = MATCHED.slice(); }
  else { applySort(); }
  page = 0;
  const dt = (performance.now() - t0).toFixed(0);
  const topN = parseInt(document.getElementById('structTopN').value, 10) || 50;
  const modeLabel = mode === 'sub' ? '子结构' : mode === 'exact' ? '完全相同' : '相似度 Top-' + topN;
  document.getElementById('status').innerHTML =
    '结构检索（' + modeLabel + '）匹配 <b>' + matches.length.toLocaleString() + '</b> 条 · 用时 ' + dt + ' ms';
  render();
  setStructStatus('');
}

function smilesToFpSet(smiles){
  let mol = null;
  try { mol = RDKit.get_mol(smiles); } catch(e){ return null; }
  if (!mol) return null;
  const fp = mol.get_morgan_fp();
  const set = new Set();
  for (let i = 0; i < fp.length; i++) if (fp[i] === '1') set.add(i);
  mol.delete();
  return set;
}

/* ---------- 输入：优先 Ketcher 画板, 否则 SMILES/SMARTS 文本框 ---------- */
async function getStructQuery(){
  if (KETCHER_READY && ketcherFrame && ketcherFrame.contentWindow && ketcherFrame.contentWindow.ketcher){
    try {
      const k = ketcherFrame.contentWindow.ketcher;
      const smi = await k.getSmiles();
      if (smi && smi.trim()) return smi.trim();
    } catch(e){ /* 回退到文本框 */ }
  }
  const ta = document.getElementById('structInput');
  const v = ta ? ta.value.trim() : '';
  return v || '';
}

/* ---------- Ketcher 画板懒加载（自托管于 lib/ketcher/） ---------- */
function loadKetcher(){
  const wrap = document.getElementById('ketcherFrameWrap');
  const btn = document.getElementById('btnLoadKetcher');
  if (!wrap || wrap.dataset.loaded) return;
  wrap.dataset.loaded = '1';
  wrap.innerHTML = '<iframe id="ketcherIframe" src="lib/ketcher/index.html" ' +
    'style="width:100%;height:540px;border:1px solid #d3d8e0;border-radius:10px"></iframe>';
  ketcherFrame = document.getElementById('ketcherIframe');
  if (btn) btn.disabled = true;

  let pollTimer = null;
  const markReady = () => {
    if (KETCHER_READY) return;
    KETCHER_READY = true;
    const st = document.getElementById('ketcherStatus');
    if (st) st.textContent = '画板已就绪：可直接从画板结构检索';
    if (pollTimer){ clearInterval(pollTimer); pollTimer = null; }
  };
  // 官方 standalone 加载完成后向父窗口 postMessage({eventType:'init'}), 监听一次即可
  window.addEventListener('message', (ev) => {
    if (ketcherFrame && ev.source === ketcherFrame.contentWindow &&
        ev.data && ev.data.eventType === 'init'){
      markReady();
    }
  });
  // 兜底：若 'init' 消息因时序漏接, 轮询 contentWindow.ketcher 是否就绪
  pollTimer = setInterval(() => {
    try {
      if (ketcherFrame && ketcherFrame.contentWindow && ketcherFrame.contentWindow.ketcher){
        markReady();
      }
    } catch(e){ /* 尚未就绪或跨域, 忽略 */ }
  }, 800);
  setTimeout(() => {
    if (!KETCHER_READY){
      const st = document.getElementById('ketcherStatus');
      if (st) st.textContent = '画板加载较慢或不可用，可直接在下方粘贴 SMILES / SMARTS 检索';
      if (pollTimer){ clearInterval(pollTimer); pollTimer = null; }
    }
  }, 15000);
}

/* ---------- 初始化 ---------- */
window.addEventListener('DOMContentLoaded', () => {
  const btn = document.getElementById('btnStructRun');
  if (btn) btn.onclick = runStructureSearch;
  const lk = document.getElementById('btnLoadKetcher');
  if (lk) lk.onclick = loadKetcher;
  document.querySelectorAll('.exSmiles').forEach(b => {
    b.onclick = () => {
      const ta = document.getElementById('structInput');
      if (ta) ta.value = b.getAttribute('data-smiles');
    };
  });
  document.querySelectorAll('input[name="structMode"]').forEach(r => {
    r.onchange = () => {
      const tn = document.getElementById('structTopNWrap');
      if (tn) tn.style.display = (r.value === 'similar') ? 'inline-flex' : 'none';
    };
  });
});
