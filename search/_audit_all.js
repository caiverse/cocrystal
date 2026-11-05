/* 全库审计: 对**每个库**抽样渲染, 逐类缺陷计数。
 * 缺陷口径(全部从"最终产物"量, 不看中间量):
 *   S 端基短键: 实际使用坐标里 最短/中位 < 0.72
 *   L 过长键  : 最长/中位 > 1.35
 *   X 交叉双键: 某双键的两条装饰线夹角 > 4°(正常应平行 ≈0)
 *   T 系链角  : 悬挂环 tether 两侧不满足 90°+180°/n (容差 2.5°)
 *   B 断键    : 块里有键但 SVG 无对应路径
 *   O 溢出    : 路径点越出画布/NaN
 */
const fs=require('fs'), path=require('path');
globalThis.window=globalThis; globalThis.addEventListener=()=>{}; window.addEventListener=()=>{};
window.document={getElementById:()=>null,createElement:()=>({style:{}}),addEventListener:()=>{}};
let src=fs.readFileSync(path.join(__dirname,'app.js'),'utf8').replace(/'use strict';/,'');
src+='\n;globalThis.__setRDKit=(m)=>{RDKit=m;};';
(0,eval)(src);

/* ---- 块坐标指标 ---- */
function blkMetrics(blk){
  const L=blk.split('\n'); const c=(L[3]||'').trim().split(/\s+/).map(Number);
  const nA=c[0]|0,nB=c[1]|0; const xy=[],el=[];
  for(let i=0;i<nA;i++){const t=(L[4+i]||'').trim().split(/\s+/);xy.push([+t[0],+t[1]]);el.push(t[3]);}
  const bonds=[],adj=Array.from({length:nA},()=>[]);
  for(let i=0;i<nB;i++){const s=L[4+nA+i]||'';const a=parseInt(s.slice(0,3))-1,b=parseInt(s.slice(3,6))-1;
    bonds.push([a,b]);adj[a].push(b);adj[b].push(a);}
  const ds=bonds.map(([a,b])=>Math.hypot(xy[a][0]-xy[b][0],xy[a][1]-xy[b][1])).sort((p,q)=>p-q);
  const med=ds[(ds.length/2)|0];
  return {nA,xy,el,bonds,adj,min:med?ds[0]/med:1,max:med?ds[ds.length-1]/med:1};
}
/* ---- SVG 指标 ---- */
function bondLines(svg){
  const map=new Map(); const re=/<path\b([^>]*)\/?>/g; let m;
  while((m=re.exec(svg))){
    const tag=m[1]; const cls=(tag.match(/class='([^']*)'/)||[])[1]||'';
    if(cls.indexOf('bond-')<0) continue;
    const ix=[...cls.matchAll(/atom-(\d+)/g)].map(x=>+x[1]);
    if(ix.length!==2) continue;
    const k=Math.min(ix[0],ix[1])+'-'+Math.max(ix[0],ix[1]);
    const d=(tag.match(/d='([^']*)'/)||[])[1]||'';
    for(const chunk of d.split(/[Mm]/).filter(s=>s.trim())){
      const pts=[]; const reN=/(-?\d+\.?\d*)[,\s]+(-?\d+\.?\d*)/g; let q;
      while((q=reN.exec(chunk))) pts.push([+q[1],+q[2]]);
      if(pts.length<2) continue;
      const A=pts[0],B=pts[pts.length-1];
      if(Math.hypot(B[0]-A[0],B[1]-A[1])<0.8) continue;
      if(!map.has(k)) map.set(k,[]);
      map.get(k).push(Math.atan2(B[1]-A[1],B[0]-A[0]));
    }
  }
  return map;
}
const devPi=(x,y)=>{let a=Math.abs(((x%Math.PI)+Math.PI)%Math.PI-((y%Math.PI)+Math.PI)%Math.PI);a=Math.min(a,Math.PI-a);return a*180/Math.PI;};
function crossedCount(map){
  let n=0;
  for(const [,angs] of map){ if(angs.length<2) continue;
    let mx=0; for(let i=0;i<angs.length;i++)for(let j=i+1;j<angs.length;j++) mx=Math.max(mx,devPi(angs[i],angs[j]));
    if(mx>4) n++;
  }
  return n;
}
/* ---- 主流程 ---- */
const init=require(path.join(__dirname,'lib/rdkit/RDKit_minimal.js')).default;
init({locateFile:f=>path.join(__dirname,'lib/rdkit',f)}).then(m=>{
  globalThis.__setRDKit(m); m.prefer_coordgen(true);
  const QUOTA=Number(process.env.QUOTA||40);
  const shards=fs.readdirSync(path.join(__dirname,'data')).filter(f=>/^products\.\d+\.json$/.test(f));
  const res={}, total={};
  for(const f of shards){
    let arr; try{arr=JSON.parse(fs.readFileSync(path.join(__dirname,'data',f),'utf8'));}catch(e){continue;}
    for(const r of arr){
      const lib=r.g||'(未标注)';
      total[lib]=(total[lib]||0)+1;
      if(!res[lib]) res[lib]=[];
      const a=res[lib];
      if(a.length<QUOTA) a.push(r);
      else { const k=total[lib]%QUOTA; if(((total[lib]/QUOTA)|0)%2===0) a[k]=r; }
    }
  }
  const libs=Object.keys(res).sort((a,b)=>total[b]-total[a]);
  const stat={};
  let t0=Date.now(), done=0, failed=0;
  for(const lib of libs){
    const s={n:0,S:0,L:0,X:0,B:0,O:0,S6:0,S7:0,L16:0};
    for(const r of res[lib]){
      const sm=r.s||''; if(!sm) continue;
      let svg=null;
      globalThis.__PRP_LASTBLK=null;
      try{ svg=molToSVG(sm); }catch(e){}
      if(!svg){ failed++; continue; }
      const blk=globalThis.__PRP_LASTBLK;
      s.n++;
      if(/NaN|undefined/.test(svg)) s.O++;
      if(blk){
        let M; try{ M=blkMetrics(blk); }catch(e){ M=null; }
        if(M){
          if(M.min<0.72) s.S++;
          if(M.min<0.60) s.S6++; else if(M.min<0.72) s.S7++;
          if(M.max>1.35) s.L++;
          if(M.max>1.60) s.L16++;
          const L2=bondLines(svg);
          if(crossedCount(L2)>0) s.X++;
          // 断键/溢出
          for(const [a,b] of M.bonds) if(!L2.has(Math.min(a,b)+'-'+Math.max(a,b))) { s.B++; break; }
        }
      }
      done++;
    }
    stat[lib]=s;
  }
  const secs=((Date.now()-t0)/1000).toFixed(0);
  const tot={n:0,S:0,L:0,X:0,B:0,O:0,S6:0,S7:0,L16:0};
  for(const k of Object.keys(stat)) for(const q of Object.keys(tot)) if(q!=='n') tot[q]+=stat[k][q];
  console.log('审计 '+libs.length+' 个库, '+done+' 张图, '+secs+'s, 渲染失败 '+failed);
  const pct=(x,n)=>n?(100*x/n).toFixed(1)+'%':'-';
  console.log('\n合计: 短键'+tot.S+'('+pct(tot.S,tot.n)+') 其中<0.60 '+tot.S6+'('+pct(tot.S6,tot.n)+') | 长键'+tot.L+'('+pct(tot.L,tot.n)+') 其中>1.60 '+tot.L16+' | 交叉双键'+tot.X+' 断键'+tot.B+' 溢出/NaN'+tot.O);
  const bad=Object.entries(stat).filter(([,s])=>s.S||s.L||s.X||s.B||s.O).sort((a,b)=>(b[1].S+b[1].L+b[1].X+b[1].B+b[1].O)-(a[1].S+a[1].L+a[1].X+a[1].B+a[1].O));
  console.log('\n有缺陷的库 ('+bad.length+' 个):');
  for(const [lib,s] of bad.slice(0,40))
    console.log('  '+(lib||'?').slice(0,26).padEnd(28)+' 样本'+String(s.n).padStart(3)+' 库容'+String(total[lib]).padStart(7)
      +' | 短'+String(s.S).padStart(2)+'(<0.6:'+s.S6+') 长'+String(s.L).padStart(2)+'(>1.6:'+s.L16+') 交叉'+String(s.X).padStart(2)+' 断'+String(s.B).padStart(2)+' 溢出'+String(s.O).padStart(2));
  fs.writeFileSync(path.join(__dirname,'_audit_all.json'), JSON.stringify({total,stat},null,1));
  console.log('\n明细已写 _audit_all.json');
});
