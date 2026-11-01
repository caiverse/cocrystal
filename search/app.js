/* 产物库检索站 — 前端逻辑
 * 数据: data/products.json (聚合后的统一 schema, 短键名)
 * 结构式: RDKit WASM 本地渲染
 */
'use strict';

let DATA = [];
let BLOBS = [];          // 每条记录的检索文本(小写)
let FILTERED = [];       // 当前结果在 DATA 中的下标(经排序后的显示顺序)
let MATCHED = [];        // 当前命中的下标(库顺序, 未经排序)
const PER_PAGE = 24;
let page = 0;
let RDKit = null;
let LIBRARY_ORDER = [];   // 产物库默认展示顺序(来自 stats.json), 供默认排序分组用

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
    console.log('RDKit 就绪');
    if (DATA.length) render();   // 数据已就绪则重绘当前页结构式
  }catch(e){
    console.warn('RDKit 初始化失败:', e);
    RDKit = null;
  }
}

/* ---------- 结构式绘制 ---------- */
function molToSVG(smiles){
  if (!RDKit || !smiles) return null;
  // 盐型 SMILES(含 '.') RDKit WASM 会整体返回 null, 改为逐片段尝试(丢弃反离子)
  const cands = smiles.includes('.') ? smiles.split('.').concat(smiles) : [smiles];
  for (const c of cands){
    let mol = null;
    try { mol = RDKit.get_mol(c); } catch(e){ mol = null; }
    if (!mol) continue;
    try {
      const svg = mol.get_svg(250, 195);
      mol.delete();
      return svg;
    } catch(e){ mol.delete(); }
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
    const mres = await fetch('data/products.manifest.json', { cache: 'no-cache' });
    if (mres.ok) manifest = await mres.json();
  } catch(e){ manifest = null; }

  // 加载统计(含产物库默认展示顺序 library_order)
  let STATS = null;
  try {
    const sres = await fetch('data/stats.json', { cache: 'no-cache' });
    if (sres.ok) STATS = await sres.json();
  } catch(e){ STATS = null; }
  LIBRARY_ORDER = (STATS && Array.isArray(STATS.library_order)) ? STATS.library_order : [];

  if (manifest && Array.isArray(manifest.shards) && manifest.shards.length){
    if (loadingText) loadingText.textContent =
      `正在并行加载 ${manifest.shards.length} 个数据分片，共约 ${(manifest.count/10000).toFixed(1)} 万条…`;
    const parts = await Promise.all(manifest.shards.map(async (name) => {
      const r = await fetch('data/' + name, { cache: 'no-cache' });
      if (!r.ok) throw new Error('分片 ' + name + ' HTTP ' + r.status);
      return r.json();
    }));
    records = [].concat(...parts);
  } else {
    if (loadingText) loadingText.textContent = '正在加载产物库数据…';
    const res = await fetch('data/products.json', { cache: 'no-cache' });
    if (!res.ok) throw new Error('HTTP ' + res.status);
    records = await res.json();
  }

  DATA = records;
  BLOBS = DATA.map(r => [
    r.n, r.e, r.f, r.s, r.g, r.r, r.a, r.d, r.o, r.c,
    r.x ? Object.values(r.x).join(' ') : ''
  ].filter(Boolean).join(' ').toLowerCase());

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
  const q = document.getElementById('search').value.trim().toLowerCase();
  const cat = document.getElementById('catFilter').value;
  const terms = q.split(/\s+/).filter(Boolean);
  const out = [];
  for (let i = 0; i < DATA.length; i++){
    if (cat && DATA[i].g !== cat) continue;
    let ok = true;
    if (terms.length){
      const blob = BLOBS[i];
      for (const t of terms){ if (!blob.includes(t)){ ok = false; break; } }
    }
    if (ok) out.push(i);
  }
  MATCHED = out;
  applySort();
  page = 0;
  const dt = (performance.now() - t0).toFixed(0);
  document.getElementById('status').innerHTML =
    '匹配 <b>' + FILTERED.length.toLocaleString() + '</b> 条 · 用时 ' + dt + ' ms';
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
  const tags = [`<span class="tag">${escapeHTML(r.g)}</span>`];
  if (r.o) tags.push(`<span class="tag">母体 ${escapeHTML(r.o)}</span>`);
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
  let html = `<button id="prev"${page<=0?' disabled':''}>上一页</button>`;
  html += `<span class="pginfo">第 ${page+1} / ${pages} 页 · 显示 ${start+1}-${end} / ${total.toLocaleString()}</span>`;
  html += `<button id="next"${page>=pages-1?' disabled':''}>下一页</button>`;
  pager.innerHTML = html;
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
  search();
});
