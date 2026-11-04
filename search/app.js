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

/* ---------- 悬挂(末端)六元环取向规范化: 让所有末端苯环水平(0°/180°), 视觉统一"正" ---------- */
/* 原理: RDKit/coordgen 不会对单点连接的悬挂环做取向规范化, 导致末端环随取代基位置歪斜。
   此处不改坐标(WASM 不接受外部 2D 坐标), 改为对 get_svg 产出的 SVG 做后处理: 识别"末端六元环"
   (仅一个母核连接点、且不与其它六元环融合), 将其绕 ipso 旋转到水平; 融合母核(吖啶/萘/蒽等)绝不旋转。 */
function _prpParseV2000(blk){ const L=blk.split('\n'); const c=L[3].trim().split(/\s+/).map(Number); const nA=c[0],nB=c[1]; return {nAtoms:nA,nBonds:nB,atomLines:L.slice(4,4+nA),bondLines:L.slice(4+nA,4+nA+nB)}; }
function _prpBuildAdj(bl,n){ const a=Array.from({length:n},()=>[]); for(const b of bl){const m=b.match(/^\s*(\d+)\s+(\d+)\s+(\d+)/); if(!m)continue; const x=+m[1]-1,y=+m[2]-1; a[x].push(y); a[y].push(x);} return a; }
function _prpFindRings(adj,n,t){ const rings=[]; for(let s=0;s<n;s++){ const path=[s],ip=new Set([s]); (function dfs(v){ if(path.length===t){ if(adj[v].includes(s)) rings.push([...path]); return;} for(const nb of adj[v]){ if(nb===s||ip.has(nb))continue; ip.add(nb); path.push(nb); dfs(nb); path.pop(); ip.delete(nb);} })(s);} const u=[],seen=new Set(); for(const r of rings){const k=[...r].sort((a,b)=>a-b).join(','); if(!seen.has(k)){seen.add(k);u.push(r);}} return u; }
function _prpBfsDist(adj,rs,start){ const d={[start]:0},q=[start]; while(q.length){const c=q.shift(); for(const nb of adj[c]) if(rs.has(nb)&&!(nb in d)){d[nb]=d[c]+1;q.push(nb);}} return d; }
function _prpDegree(adj){ return adj.map(a=>a.length); }
function _prpIsFused(ring, rings){ for(const o of rings){ if(o===ring) continue; const os=new Set(o); for(let i=0;i<ring.length;i++){ const a=ring[i], b=ring[(i+1)%ring.length]; if(os.has(a)&&os.has(b)) return true; } } return false; }
function _prpGetPaths(svg){ const re=/<path\b([^>]*)\/?>/g; const res=[]; let m; while((m=re.exec(svg))){ const tag=m[0]; const cls=(tag.match(/class='([^']*)'/)||[])[1]||''; const d=(tag.match(/d='([^']*)'/)||[])[1]||''; res.push({tag,cls,d}); } return res; }
function _prpAtomsInClass(cls){ const s=new Set(); let m; const re=/atom-(\d+)/g; while((m=re.exec(cls))) s.add(+m[1]); return s; }
function _prpBondEndpoints(d){ const e=[]; let m; const re=/(-?\d+\.?\d*),(-?\d+\.?\d*)/g; while((m=re.exec(d))) e.push([parseFloat(m[1]),parseFloat(m[2])]); return e; }
function _prpGetAtomPos(paths, k){ // 旋转不变: 取每条入射键"最靠近质心"的端点再平均, 即原子 k 的屏幕坐标
  const all=[]; const eps=[];
  for(const p of paths){ if(!p.cls.includes('bond-')) continue; const ix=_prpAtomsInClass(p.cls); if(!ix.has(k)) continue; const e=_prpBondEndpoints(p.d); if(e.length<2) continue; eps.push(e); for(const pt of e) all.push(pt); }
  if(!eps.length) return null;
  let cx=0,cy=0; for(const pt of all){cx+=pt[0];cy+=pt[1];} cx/=all.length; cy/=all.length;
  let sx=0,sy=0,n=0;
  for(const e of eps){ const d0=(e[0][0]-cx)**2+(e[0][1]-cy)**2, d1=(e[1][0]-cx)**2+(e[1][1]-cy)**2; const ep=d0<=d1?e[0]:e[1]; sx+=ep[0]; sy+=ep[1]; n++; }
  return n? [sx/n, sy/n] : null;
}
function _prpRotateD(d,piv,th){ return d.replace(/(-?\d+\.?\d*),(-?\d+\.?\d*)/g,(m,x,y)=>{ const dx=parseFloat(x)-piv[0],dy=parseFloat(y)-piv[1]; const nx=piv[0]+dx*Math.cos(th)-dy*Math.sin(th); const ny=piv[1]+dx*Math.sin(th)+dy*Math.cos(th); return nx.toFixed(2)+','+ny.toFixed(2); }); }
function postRotatePendantRings(svg, blk){
  if(!svg || !blk) return svg;
  const P=_prpParseV2000(blk); const adj=_prpBuildAdj(P.bondLines,P.nAtoms); const deg=_prpDegree(adj);
  const rings=_prpFindRings(adj,P.nAtoms,6);
  const out=[];
  for(const ring of rings){
    if(_prpIsFused(ring, rings)) continue;                       // 融合环(吖啶/萘/蒽母核)绝不旋转
    const rs=new Set(ring);
    let ext=[]; for(const a of ring) for(const nb of adj[a]) if(!rs.has(nb)) ext.push(nb);
    if(ext.length<1) continue;                                    // 完全内部的环
    let ipso=-1,best=-1; for(const a of ring){ let c=0; for(const nb of adj[a]) if(!rs.has(nb)) c++; if(c>best){best=c;ipso=a;} }
    const d=_prpBfsDist(adj,rs,ipso); let para=-1,pb=-1; for(const a of ring) if(d[a]>pb){pb=d[a];para=a;}
    if(para<0) continue;
    const tether = ext.find(x=>deg[x]>1); const t = (tether!==undefined)?tether:ext[0]; // 母核连接点(deg>1); 叶取代基(F/CH3…)随环一起转
    const S=new Set(ring); for(const a of ring) for(const nb of adj[a]) if(!rs.has(nb)&&nb!==t) S.add(nb);
    out.push({ipso,para,t,S});
  }
  const paths=_prpGetPaths(svg);
  for(const r of out){
    const ip=_prpGetAtomPos(paths,r.ipso); const pp=_prpGetAtomPos(paths,r.para);
    if(!ip||!pp) continue;
    const cur=Math.atan2(pp[1]-ip[1],pp[0]-ip[0])*180/Math.PI;
    let target=0; const ep=_prpGetAtomPos(paths,r.t); if(ep && ep[0]>ip[0]) target=180; else target=0;
    const th=(target-cur)*Math.PI/180;
    if(Math.abs(((th*180/Math.PI+180)%360)-180)<1e-6) continue;  // 已水平, 跳过
    for(const p of paths){ const ai=_prpAtomsInClass(p.cls); let rot=false; if(p.cls.includes('bond-')) rot=ai.size===2&&[...ai].every(x=>r.S.has(x)); else rot=ai.size===1&&r.S.has([...ai][0]); if(!rot) continue; const nd=_prpRotateD(p.d,ip,th); svg=svg.replace(p.tag, p.tag.replace(/d='[^']*'/, "d='"+nd+"'")); }
  }
  return svg;
}

/* ---------- 结构式绘制 ---------- */
function molToSVG(smiles){
  if (!RDKit || !smiles) return null;
  // 盐型 SMILES(含 '.') 优先用整体绘制; 若整体失败, 再按片段尝试, 但片段需按长度降序,
  // 优先绘制主分子而非反离子(如 [H+] / [Cl-]), 否则会只画出 H⁺ 等反离子。
  const parts = smiles.includes('.') ? smiles.split('.') : [smiles];
  const cands = [smiles].concat(parts.slice().sort((a, b) => b.length - a.length));
  for (const c of cands){
    let mol = null;
    try { mol = RDKit.get_mol(c); } catch(e){ mol = null; }
    if (!mol) continue;
    try {
      const svg0 = mol.get_svg(250, 195);
      let svg = svg0;
      try { const blk = mol.get_new_coords(); svg = postRotatePendantRings(svg0, blk); } catch(e){}
      mol.delete();
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
