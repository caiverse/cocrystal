/* 结构检索（子结构 / 完全相同 / 相似度 Top-N）—— 基于站内已有 RDKit WASM, 纯前端、无服务器
 * 依赖 app.js 提供的全局: DATA, RDKit, MATCHED, FILTERED, page, PER_PAGE, applySort(), render()
 * 与 ATOM_CNT（每条记录重原子数）。
 */
'use strict';

let FP_READY = false;
let FP_N = 0;
let FP_OFFSETS = null;     // Uint32Array(N+1)
let FP_BITS = null;        // Uint16Array(所有位, 全局拼接)
let KETCHER_READY = false;
let ketcherFrame = null;
// window.SIM_MAP: 记录 id -> 相似度(0..1), 仅相似度模式用于卡片展示

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
async function loadFP(){
  if (FP_READY) return;
  const t0 = performance.now();
  const res = await fetch('data/fp.bin', {cache: 'no-cache'});
  const buf = await res.arrayBuffer();
  const dv = new DataView(buf);
  FP_N = dv.getUint32(0, true);
  FP_OFFSETS = new Uint32Array(FP_N + 1);
  let total = 0;
  for (let i = 0; i < FP_N; i++){
    FP_OFFSETS[i] = total;
    total += dv.getUint16(4 + i * 2, true);
  }
  FP_OFFSETS[FP_N] = total;
  FP_BITS = new Uint16Array(total);
  let p = 4 + FP_N * 2;
  for (let i = 0; i < FP_N; i++){
    const cnt = dv.getUint16(p, true); p += 2;
    for (let j = 0; j < cnt; j++){ FP_BITS[FP_OFFSETS[i] + j] = dv.getUint16(p, true); p += 2; }
  }
  FP_READY = true;
  console.log('fp.bin 加载', (performance.now() - t0).toFixed(0), 'ms, 位总数', total);
}
function binSearch(arr, s, e, b){
  let lo = s, hi = e - 1;
  while (lo <= hi){
    const mid = (lo + hi) >> 1;
    const v = arr[mid];
    if (v === b) return true;
    if (v < b) lo = mid + 1; else hi = mid - 1;
  }
  return false;
}
function tanimoto(qbits, idx){
  const s = FP_OFFSETS[idx], e = FP_OFFSETS[idx + 1];
  let inter = 0;
  for (let k = 0; k < qbits.length; k++){ if (binSearch(FP_BITS, s, e, qbits[k])) inter++; }
  const rcount = e - s;
  const denom = qbits.length + rcount - inter;
  return denom ? inter / denom : 0;
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

/* ---------- 主入口 ---------- */
async function runStructureSearch(){
  if (!RDKit){ alert('化学结构引擎尚未就绪，请稍候重试'); return; }
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

  const matches = [];
  if (mode === 'sub'){
    setStructStatus('预筛中…');
    const cands = [];
    for (let i = 0; i < DATA.length; i++){
      if (!DATA[i].s) continue;
      if (ATOM_CNT[i] < QH) continue;
      if (!elemSuperset(DATA[i].f, QE)) continue;
      cands.push(i);
    }
    setStructStatus('候选 ' + cands.length.toLocaleString() + ' 条，子结构匹配中…');
    let done = 0; const B = 2000; const molCache = new Map();
    for (let ci = 0; ci < cands.length; ci++){
      const i = cands[ci];
      let mol = molCache.get(i);
      if (mol === undefined){ try { mol = RDKit.get_mol(DATA[i].s); } catch(e){ mol = null; } molCache.set(i, mol); }
      if (!mol){ done++; if (done % B === 0) await yieldUI(done / cands.length); continue; }
      let res = null;
      try { res = JSON.parse(mol.get_substruct_match(qmol)); } catch(e){ res = null; }
      if (res && res.atoms && res.atoms.length > 0) matches.push(i);
      done++; if (done % B === 0) await yieldUI(done / cands.length);
    }
  } else if (mode === 'exact'){
    setStructStatus('精确匹配中…');
    for (let i = 0; i < DATA.length; i++){
      if (!DATA[i].s) continue;
      if (ATOM_CNT[i] !== QH) continue;
      if (!elemSuperset(DATA[i].f, QE)) continue;
      let mol = null; try { mol = RDKit.get_mol(DATA[i].s); } catch(e){ mol = null; }
      if (!mol) continue;
      let res = null;
      try { res = JSON.parse(mol.get_substruct_match(qmol)); } catch(e){ res = null; }
      if (res && res.atoms && res.atoms.length === QH) matches.push(i);
    }
  } else if (mode === 'similar'){
    await loadFP();
    setStructStatus('相似度计算中…');
    const topN = parseInt(document.getElementById('structTopN').value, 10) || 50;
    const sims = [];
    for (let i = 0; i < DATA.length; i++){
      const s = FP_OFFSETS[i], e = FP_OFFSETS[i + 1];
      if (e === s) continue;
      const t = tanimoto(qfpSet, i);
      if (t > 0) sims.push([i, t]);
    }
    sims.sort((a, b) => b[1] - a[1]);
    const top = sims.slice(0, topN);
    for (const [i, t] of top){ matches.push(i); if (window.SIM_MAP) window.SIM_MAP[DATA[i].i] = t; }
  }

  if (qmol){ try { qmol.delete(); } catch(e){} }
  MATCHED = matches;
  if (mode === 'similar'){ FILTERED = MATCHED.slice(); }
  else { applySort(); }
  page = 0;
  const dt = (performance.now() - t0).toFixed(0);
  const modeLabel = mode === 'sub' ? '子结构' : mode === 'exact' ? '完全相同' : '相似度 Top-' + (parseInt(document.getElementById('structTopN').value, 10) || 50);
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
function yieldUI(frac){
  setStructProgress(frac);
  return new Promise(r => setTimeout(r, 0));
}
function SIM_MAP_RESET(){ if (window.SIM_MAP) window.SIM_MAP = {}; }
function setStructStatus(t){ const e = document.getElementById('structStatus'); if (e) e.textContent = t; }
function setStructProgress(f){
  const b = document.getElementById('structProgressBar');
  if (b) b.style.width = Math.max(0, Math.min(100, Math.round(f * 100))) + '%';
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
