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
function _prpBondEndpoints(d){ const e=[]; let m; const re=/(-?\d+\.?\d*)[,\s]+(-?\d+\.?\d*)/g; while((m=re.exec(d))) e.push([parseFloat(m[1]),parseFloat(m[2])]); return e; }
function _prpCentroid(d){ const e=_prpBondEndpoints(d); if(!e.length) return null; let sx=0,sy=0; for(const p of e){sx+=p[0];sy+=p[1];} return [sx/e.length, sy/e.length]; }
function _prpPtSegDist(p,a,b){ // 点到线段距离
  const vx=b[0]-a[0], vy=b[1]-a[1]; const wx=p[0]-a[0], wy=p[1]-a[1];
  const c1=vx*wx+vy*wy; if(c1<=0) return Math.hypot(p[0]-a[0],p[1]-a[1]);
  const c2=vx*vx+vy*vy; if(c2<=c1) return Math.hypot(p[0]-b[0],p[1]-b[1]);
  const t=c1/c2; return Math.hypot(p[0]-(a[0]+t*vx),p[1]-(a[1]+t*vy));
}
function _prpGetAtomPos(paths, k){ // 原子 k = 其所有入射键(线段)的汇聚点(对"键在 midpoint 断开为两段"与"芳香双线偏移"均稳健)
  const segs=[]; const cands=[];
  for(const p of paths){
    if(!p.cls.includes('bond-')) continue;
    const ix=_prpAtomsInClass(p.cls); if(!ix.has(k)) continue;
    const e=_prpBondEndpoints(p.d); if(e.length<2) continue;
    for(let i=0;i+1<e.length;i+=2){ const a=e[i], b=e[i+1]; segs.push([a,b]); cands.push(a); cands.push(b); }
  }
  if(!segs.length) return null;
  const tol=2.5;
  let best=null,bestScore=0;
  for(const c of cands){ let score=0; for(const s of segs){ if(_prpPtSegDist(c,s[0],s[1])<=tol) score++; } if(score>bestScore){ bestScore=score; best=c; } }
  if(!best) return null;
  let sx=0,sy=0,n=0; for(const q of cands){ if(Math.hypot(best[0]-q[0],best[1]-q[1])<=tol){ sx+=q[0]; sy+=q[1]; n++; } }
  return n?[sx/n, sy/n] : best;
}
function _prpRotateD(d,piv,th){ return d.replace(/(-?\d+\.?\d*)[,\s]+(-?\d+\.?\d*)/g,(m,x,y)=>{ const dx=parseFloat(x)-piv[0],dy=parseFloat(y)-piv[1]; const nx=piv[0]+dx*Math.cos(th)-dy*Math.sin(th); const ny=piv[1]+dx*Math.sin(th)+dy*Math.cos(th); return nx.toFixed(2)+','+ny.toFixed(2); }); }
function _prpAffineD(d, fn){ return d.replace(/(-?\d+\.?\d*)[,\s]+(-?\d+\.?\d*)/g,(m,x,y)=>{ const r=fn(parseFloat(x),parseFloat(y)); return r[0].toFixed(2)+','+r[1].toFixed(2); }); }
/* 悬挂(末端)六元环取向规范化: 让"无取代基的干净末端苯环"水平(0°/180°), 视觉统一"正"。
 * 重要: 仅对「非融合 且 仅有一个外部连接(即环上无任何取代基)」的末端苯环做旋转;
 *   旋转绕该环与母核相连的 ipso 原子(pivot, 固定不动), 整环刚体旋转, 数学上不可能断键。
 *   其余情形一律不旋转, 保留 RDKit 原布局:
 *     · 带取代基的苯环(乙酰/甲氧/甲基/F… )—— 旧逻辑会撕裂取代基支链, 现已排除;
 *     · 联苯/环-环相连(tether 落在另一个候选环内)—— 两环都不转, 避免撕裂环间键;
 *     · 稠合母核(萘/蒽/吖啶)—— 绝不旋转。
 *   这样可在消除"断键/错位"的前提下, 仅对最安全的一类结构做水平化。 */
/* 悬挂(末端)芳环取向规范化 + 系链加长(两阶段):
 * 只对「真正悬挂的末端芳环」处理, 判定严格:
 *   (1) 去掉该环全部原子后, 分子其余部分仍连通(末端/叶子环, 不是中心环);
 *   (2) 该环不与任何其它环共用原子(排除稠合/螺/桥联环——它们是刚性母核的一部分)。
 *   处理: 绕其与母核相连的 ipso 刚体旋转到正立(0°/180°); 若带取代基, 再沿系链轴向平移片断以"拉长系链"
 *         (加长放在缩放居中之后做, 否则被统一缩放抵消)。
 *   片断 = 环原子 + 沿非系链键 BFS 可达的取代基; 刚体操作, 数学上不断键。 */
function postRotatePendantRings(svg, blk){
  if(!svg || !blk) return svg;
  const P=_prpParseV2000(blk); const adj=_prpBuildAdj(P.bondLines,P.nAtoms);
  const rings=_prpFindRings(adj,P.nAtoms,6);
  // 原子 -> 所属环索引列表(用于判定"是否共用原子")
  const atomRings=Array.from({length:P.nAtoms},()=>[]);
  rings.forEach((r,i)=>{ for(const a of r) atomRings[a].push(i); });
  const sharesAtom=(ri)=>{ for(const a of rings[ri]) if(atomRings[a].length>1) return true; return false; };
  // 去掉该环原子后, 其余原子是否仍连通(叶子环判定)。
  // 容差: 环上取代基(Br/Cl/OH/CF3/NO2 等)被一并摘除后会成为孤立小碎片, 这些不算"骨架断裂"。
  const SUBST_MAX=5;
  const isLeaf=(ringAtoms)=>{ const rem=new Array(P.nAtoms).fill(true); for(const a of ringAtoms) rem[a]=false;
    const rest=[]; for(let i=0;i<P.nAtoms;i++) if(rem[i]) rest.push(i); if(!rest.length) return false;
    const visited=new Set(); const comps=[];
    for(const s0 of rest){ if(visited.has(s0)) continue; const comp=[]; const q=[s0]; visited.add(s0);
      while(q.length){ const u=q.pop(); comp.push(u); for(const nb of adj[u]){ if(rem[nb]&&!visited.has(nb)){ visited.add(nb); q.push(nb); } } }
      comps.push(comp); }
    if(!comps.length) return false;
    comps.sort((a,b)=>b.length-a.length);
    for(let k=1;k<comps.length;k++){ if(comps[k].length>SUBST_MAX) return false; } // 存在较大的第二分量 -> 骨架被切断 -> 非叶子环
    return true; };
  // 估算从 nb 出发、不进入 ca 的可达原子数(挑出连向主核的系链, 区分"连向母核"与"取代基")
  const reachExcl=(nb, ca)=>{ const seen=new Set([nb]); const st=[nb]; let cnt=0;
    while(st.length){ const u=st.pop(); cnt++; for(const v of adj[u]){ if(seen.has(v))continue; if(ca.has(v))continue; seen.add(v); st.push(v); } } return cnt; };
  // 收集可旋转的末端环
  const jobs=[];
  for(let ri=0; ri<rings.length; ri++){
    if(sharesAtom(ri)) continue;        // 稠合/螺/桥联 -> 刚性母核, 不旋转
    if(!isLeaf(rings[ri])) continue;     // 去掉后其余不连通(中心环) -> 不旋转
    const ringAtoms=new Set(rings[ri]);
    // 系链 = 环原子 -> 非环原子的外部键中, 锚点外连可达最多者(区分"连向母核"与"取代基")
    let tether=null, best=-1;
    for(const a of rings[ri]){ for(const nb of adj[a]){ if(ringAtoms.has(nb)) continue; const sz=reachExcl(nb, ringAtoms); if(sz>best){ best=sz; tether={ip:a, anchor:nb}; } } }
    if(!tether) continue;                // 无外部连接(孤立环)跳过
    // 片断 = 环原子 + 沿非系链键 BFS 可达的取代基
    const frag=new Set(ringAtoms); const q=[...ringAtoms];
    while(q.length){ const u=q.pop(); for(const nb of adj[u]){
      if(frag.has(nb)) continue;
      if((u===tether.ip&&nb===tether.anchor)||(u===tether.anchor&&nb===tether.ip)) continue; // 仅跳过系链
      frag.add(nb); q.push(nb); } }
    if(frag.size > Math.max(8, Math.floor(P.nAtoms*0.6))) continue; // 守卫: 片断过大则跳过, 避免破坏主结构
    const hasSubst = frag.size > ringAtoms.size;
    jobs.push({ringIdx:ri, ringAtoms, tether, frag, hasSubst});
  }
  const TETHER_FRAC=0.35, TETHER_MAXEXT=8;
  // 阶段A: 取向(绕 ipso 刚体旋转到水平; 系链旋转阶段不动因为 ipso 固定)
  { const paths=_prpGetPaths(svg); const ap={}; const getPos=(a)=>{ if(ap[a]!==undefined) return ap[a]; const p=_prpGetAtomPos(paths,a); ap[a]=p||null; return ap[a]; };
    for(const j of jobs){
      const {ip, anchor}=j.tether; const ipPos=getPos(ip), anPos=getPos(anchor); if(!ipPos||!anPos) continue;
      let np=-1; for(const nb of adj[ip]){ if(rings[j.ringIdx].includes(nb)){ np=nb; break; } }
      if(np<0) continue; const npPos=getPos(np); if(!npPos) continue;
      const cur=Math.atan2(npPos[1]-ipPos[1], npPos[0]-ipPos[0])*180/Math.PI;
      const norm=(x)=>{ while(x>180)x-=360; while(x<-180)x+=360; return x; };
      let target=Math.abs(norm(0-cur))<=Math.abs(norm(180-cur))?0:180;
      let th=(target-cur)*Math.PI/180;
      // 翻转判定: 旋转后环质心应位于 anchor 反方向; 否则翻转 180°
      let sx=0, sy=0, nn=0; for(const a of rings[j.ringIdx]){ const p=getPos(a); if(p){ sx+=p[0]; sy+=p[1]; nn++; } }
      if(nn){ const vx=sx/nn-ipPos[0], vy=sy/nn-ipPos[1];
        const rx=vx*Math.cos(th)-vy*Math.sin(th), ry=vx*Math.sin(th)+vy*Math.cos(th);
        const dA=[anPos[0]-ipPos[0], anPos[1]-ipPos[1]];
        if(rx*dA[0]+ry*dA[1]>0){ target=(target+180)%360; th=(target-cur)*Math.PI/180; } }
      const rot=(x,y)=>{ const rx=ipPos[0]+(x-ipPos[0])*Math.cos(th)-(y-ipPos[1])*Math.sin(th);
        const ry=ipPos[1]+(x-ipPos[0])*Math.sin(th)+(y-ipPos[1])*Math.cos(th); return [rx, ry]; };
      for(const pa of paths){ const ai=_prpAtomsInClass(pa.cls); if(ai.size===0) continue;
        const atoms=[...ai];
        if(ai.size===2 && j.frag.has(atoms[0])!==j.frag.has(atoms[1])) continue; // 系链: 一端在片断内、一端在外 -> ipso 固定不动
        if(!(atoms.some(a=>j.frag.has(a)))) continue; // 仅当路径涉及的任一原子属于片断才旋转(含单原子标签 N⁺/O⁻/Cl 等字形, 随环一起转)
        const nd=_prpAffineD(pa.d, rot); svg=svg.replace(pa.tag, pa.tag.replace(/d='[^']*'/, "d='"+nd+"'"));
      }
    }
  }
  // 整体缩放居中回原画布, 保证不裁切且显示尺寸稳定
  { const _ps=_prpGetPaths(svg); let _mnx=1e9,_mny=1e9,_mxx=-1e9,_myy=-1e9;
    for(const _p of _ps){ const _e=_prpBondEndpoints(_p.d); for(const _q of _e){ if(_q[0]<_mnx)_mnx=_q[0]; if(_q[0]>_mxx)_mxx=_q[0]; if(_q[1]<_mny)_mny=_q[1]; if(_q[1]>_myy)_myy=_q[1]; } }
    const _W=250,_H=195,_mg=10;
    if(!(_mnx>_mxx) && !(_mnx>=_mg && _mxx<=_W-_mg && _mny>=_mg && _myy<=_H-_mg)){
      const _bw=Math.max(_mxx-_mnx,1e-6), _bh=Math.max(_myy-_mny,1e-6);
      const _s=Math.min((_W-2*_mg)/_bw,(_H-2*_mg)/_bh);
      const _tx=(_W-_bw*_s)/2-_mnx*_s, _ty=(_H-_bh*_s)/2-_mny*_s;
      for(const _p of _ps){ const _nd=_prpAffineD(_p.d,(x,y)=>[_s*x+_tx,_s*y+_ty]); svg=svg.replace(_p.tag, _p.tag.replace(/d='[^']*'/, "d='"+_nd+"'")); }
    }
  }
  // 阶段B: 加长(缩放居中后, 仅带取代基者)—— 在缩放之后做, 否则被统一缩放抵消
  { const paths=_prpGetPaths(svg); const ap={}; const getPos=(a)=>{ if(ap[a]!==undefined) return ap[a]; const p=_prpGetAtomPos(paths,a); ap[a]=p||null; return ap[a]; };
    for(const j of jobs){ if(!j.hasSubst) continue;
      const {ip, anchor}=j.tether; const ipPos=getPos(ip), anPos=getPos(anchor); if(!ipPos||!anPos) continue;
      let ux=ipPos[0]-anPos[0], uy=ipPos[1]-anPos[1]; const ulen=Math.hypot(ux,uy)||1; ux/=ulen; uy/=ulen;
      const ext=Math.min(TETHER_MAXEXT, ulen*TETHER_FRAC); const dx=ux*ext, dy=uy*ext;
      const tr=(x,y)=>[x+dx, y+dy];
      for(const pa of paths){ const ai=_prpAtomsInClass(pa.cls); if(ai.size===0) continue;
        const atoms=[...ai];
        const isTether = ai.size===2 && (j.frag.has(atoms[0])!==j.frag.has(atoms[1]));
        if(isTether){ const a0=atoms[0], a1=atoms[1]; const p0=getPos(a0), p1=getPos(a1); if(!p0||!p1) continue;
          let mov,fixed; if(j.frag.has(a0)){ mov=tr(p0[0],p0[1]); fixed=p1; } else { mov=tr(p1[0],p1[1]); fixed=p0; }
          const nd='M '+mov[0].toFixed(2)+' '+mov[1].toFixed(2)+' L '+fixed[0].toFixed(2)+' '+fixed[1].toFixed(2);
          svg=svg.replace(pa.tag, pa.tag.replace(/d='[^']*'/, "d='"+nd+"'")); continue; }
        if(!(atoms.some(a=>j.frag.has(a)))) continue; // 含单原子标签 N⁺/O⁻/Cl 等字形, 随片断一起平移加长
        const nd=_prpAffineD(pa.d, tr); svg=svg.replace(pa.tag, pa.tag.replace(/d='[^']*'/, "d='"+nd+"'"));
      }
    }
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
