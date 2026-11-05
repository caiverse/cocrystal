/* 用**线上同一个渲染器**(search/app.js 的 molToSVG) 重新生成首页侧栏 14 张静态 SVG,
 * 保证线上线下完全一致。规格沿用原文件: 300×240。
 *
 * 另加一步「取景归一」(_reframe): RDKit get_svg 会留大片白边(墨迹最高只到 175/240),
 * 与旧版 Python 生成器"按最大边归一"的口径不一致; 这里按真实墨迹 bbox 重新拟合
 * transform, 让结构在 300×240 画布内居中并占满主要方向, 边距固定 10px。 */
const fs=require('fs'), path=require('path');
const ROOT=path.join(__dirname,'..');
globalThis.window=globalThis; globalThis.addEventListener=()=>{}; window.addEventListener=()=>{};
window.document={getElementById:()=>null,createElement:()=>({style:{}}),addEventListener:()=>{}};
let src=fs.readFileSync(path.join(__dirname,'app.js'),'utf8').replace(/'use strict';/,'');
src+='\n;globalThis.__setRDKit=(m)=>{RDKit=m;};';
(0,eval)(src);

const W=300, H=240, MARGIN=10;

/* 取墨迹 bbox(本 SVG 只有 M/L/Q, 数字两两成坐标), 再拟合 transform 居中铺满 */
function _reframe(svg){
  const ds=[...svg.matchAll(/<path[^>]*\bd='([^']+)'/g)].map(m=>m[1]);
  if(!ds.length) return svg;
  let x0=1e9,y0=1e9,x1=-1e9,y1=-1e9;
  for(const d of ds){
    const nums=d.match(/-?\d+(?:\.\d+)?/g)||[];
    for(let i=0;i+1<nums.length;i+=2){
      const x=+nums[i], y=+nums[i+1];
      if(x<x0)x0=x; if(x>x1)x1=x; if(y<y0)y0=y; if(y>y1)y1=y;
    }
  }
  const iw=x1-x0, ih=y1-y0;
  if(!(iw>0&&ih>0)) return svg;
  const k=Math.min((W-2*MARGIN)/iw, (H-2*MARGIN)/ih);
  const tx=W/2-k*(x0+x1)/2, ty=H/2-k*(y0+y1)/2;
  return svg.replace(/transform="translate\([^)]*\)\s*scale\([^)]*\)"/,
    `transform="translate(${tx.toFixed(2)},${ty.toFixed(2)}) scale(${k.toFixed(4)})"`);
}

// 从 rails.js 取 14 条(code/name/svg)
const rails=fs.readFileSync(path.join(ROOT,'rails.js'),'utf8');
const m=rails.match(/var RAILS = (\[.*?\]);/s);
if(!m){ console.error('rails.js 解析失败'); process.exit(1); }
const RAILS=JSON.parse(m[1]);

// code -> SMILES
const rec={};
for(const f of fs.readdirSync(path.join(__dirname,'data')).filter(f=>/^products\.\d+\.json$/.test(f))){
  let arr; try{arr=JSON.parse(fs.readFileSync(path.join(__dirname,'data',f),'utf8'));}catch(e){continue;}
  for(const it of arr) if(it.code) rec[it.code]=it;
}
const init=require(path.join(__dirname,'lib/rdkit/RDKit_minimal.js')).default;
init({locateFile:f=>path.join(__dirname,'lib/rdkit',f)}).then(mod=>{
  globalThis.__setRDKit(mod); mod.prefer_coordgen(true);
  let ok=0, miss=0;
  for(const d of RAILS){
    const r=rec[d.code];
    if(!r||!r.s){ console.log('缺 SMILES: '+d.code+' '+(d.name||'').slice(0,30)); miss++; continue; }
    let svg=molToSVG(r.s, W, H);
    if(!svg){ console.log('渲染失败: '+d.code); miss++; continue; }
    svg=_reframe(svg);
    fs.writeFileSync(path.join(ROOT, d.svg.replace(/\//g,path.sep)), svg);
    ok++;
    console.log('OK  '+d.svg+'  '+d.code+'  '+(d.name||'').slice(0,34));
  }
  console.log('\n生成 '+ok+' 张, 失败/缺失 '+miss+' 张');
});
