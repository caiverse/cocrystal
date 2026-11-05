/* 对已知缺陷 SMILES 复测: 渲染后量"实际使用坐标"的最短/中位、最长/中位。 */
const fs=require('fs'), path=require('path');
globalThis.window=globalThis; globalThis.addEventListener=()=>{}; window.addEventListener=()=>{};
window.document={getElementById:()=>null,createElement:()=>({style:{}}),addEventListener:()=>{}};
let src=fs.readFileSync(path.join(__dirname,'app.js'),'utf8').replace(/'use strict';/,'');
src+='\n;globalThis.__setRDKit=(m)=>{RDKit=m;};';
(0,eval)(src);

function stats(blk){
  const L=blk.split('\n'); const c=(L[3]||'').trim().split(/\s+/).map(Number);
  const nA=c[0]|0,nB=c[1]|0; const xy=[];
  for(let i=0;i<nA;i++){const t=(L[4+i]||'').trim().split(/\s+/);xy.push([+t[0],+t[1]]);}
  const ds=[];
  for(let i=0;i<nB;i++){const s=L[4+nA+i]||'';const a=parseInt(s.slice(0,3))-1,b=parseInt(s.slice(3,6))-1;
    ds.push(Math.hypot(xy[a][0]-xy[b][0],xy[a][1]-xy[b][1]));}
  ds.sort((p,q)=>p-q); const med=ds[(ds.length/2)|0];
  return {min:+(ds[0]/med).toFixed(3), max:+(ds[ds.length-1]/med).toFixed(3), med:+med.toFixed(3), n:ds.length};
}

const list=fs.readFileSync(process.argv[2],'utf8').split('\n').map(s=>s.trim()).filter(Boolean);
const init=require(path.join(__dirname,'lib/rdkit/RDKit_minimal.js')).default;
init({locateFile:f=>path.join(__dirname,'lib/rdkit',f)}).then(m=>{
  globalThis.__setRDKit(m); m.prefer_coordgen(true);
  let bad=0;
  for(const sm of list){
    globalThis.__PRP_LASTBLK=null;
    let svg=null; try{ svg=molToSVG(sm); }catch(e){}
    const blk=globalThis.__PRP_LASTBLK;
    if(!svg||!blk){ console.log('失败 '+sm.slice(0,40)); bad++; continue; }
    const s=stats(blk);
    const ok=(s.min>=0.72 && s.max<=1.35);
    if(!ok) bad++;
    console.log((ok?'OK  ':'BAD ')+'min='+String(s.min).padEnd(6)+' max='+String(s.max).padEnd(6)+' med='+s.med+'  '+sm.slice(0,54));
  }
  console.log('\n共 '+list.length+' 个, 仍不达标 '+bad+' 个');
});
