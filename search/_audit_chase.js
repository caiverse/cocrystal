/* 追查残余短键: 只扫指定库, 高配额, 打印每个命中分子的 SMILES/键信息/环属性。 */
const fs=require('fs'), path=require('path');
globalThis.window=globalThis; globalThis.addEventListener=()=>{}; window.addEventListener=()=>{};
window.document={getElementById:()=>null,createElement:()=>({style:{}}),addEventListener:()=>{}};
let src=fs.readFileSync(path.join(__dirname,'app.js'),'utf8').replace(/'use strict';/,'');
src+='\n;globalThis.__setRDKit=(m)=>{RDKit=m;};';
(0,eval)(src);

function analyse(blk){
  const L=blk.split('\n'); const c=(L[3]||'').trim().split(/\s+/).map(Number);
  const nA=c[0]|0,nB=c[1]|0; const xy=[],el=[];
  for(let i=0;i<nA;i++){const t=(L[4+i]||'').trim().split(/\s+/);xy.push([+t[0],+t[1]]);el.push(t[3]);}
  const bonds=[],adj=Array.from({length:nA},()=>[]);
  for(let i=0;i<nB;i++){const s=L[4+nA+i]||'';const a=parseInt(s.slice(0,3))-1,b=parseInt(s.slice(3,6))-1;const ty=parseInt(s.slice(6,9));
    bonds.push({a,b,ty});adj[a].push(b);adj[b].push(a);}
  const ds=bonds.map(o=>Math.hypot(xy[o.a][0]-xy[o.b][0],xy[o.a][1]-xy[o.b][1])).sort((p,q)=>p-q);
  const med=ds[(ds.length/2)|0];
  const isRing=(a,c)=>{const seen=new Uint8Array(nA);seen[a]=1;const st=[a];
    while(st.length){const u=st.pop();for(const v of adj[u]){if((u===a&&v===c)||(u===c&&v===a))continue;
      if(v===c)return true; if(seen[v])continue; seen[v]=1; st.push(v);}} return false;};
  const out=[];
  for(const o of bonds){
    const d=Math.hypot(xy[o.a][0]-xy[o.b][0],xy[o.a][1]-xy[o.b][1]);
    const r=d/med;
    if(r<0.72||r>1.35) out.push({a:o.a,b:o.b,ty:o.ty,r:+r.toFixed(3),
      ea:el[o.a],eb:el[o.b],da:adj[o.a].length,db:adj[o.b].length,ring:isRing(o.a,o.b)});
  }
  return {nA,nB,med:+med.toFixed(3),bad:out};
}

const init=require(path.join(__dirname,'lib/rdkit/RDKit_minimal.js')).default;
init({locateFile:f=>path.join(__dirname,'lib/rdkit',f)}).then(m=>{
  globalThis.__setRDKit(m); m.prefer_coordgen(true);
  const QUOTA=Number(process.env.QUOTA||300);
  const ONLY=(process.env.LIBS||'').split('|').filter(Boolean);
  const shards=fs.readdirSync(path.join(__dirname,'data')).filter(f=>/^products\.\d+\.json$/.test(f));
  const bucket={};
  for(const f of shards){
    let arr; try{arr=JSON.parse(fs.readFileSync(path.join(__dirname,'data',f),'utf8'));}catch(e){continue;}
    for(const r of arr){
      const lib=r.g||'(未标注)';
      if(ONLY.length && !ONLY.some(k=>lib.includes(k))) continue;
      if(!bucket[lib]) bucket[lib]=[];
      if(bucket[lib].length<QUOTA) bucket[lib].push(r);
    }
  }
  let n=0, hits=0; const t0=Date.now();
  for(const lib of Object.keys(bucket)){
    for(const r of bucket[lib]){
      const sm=r.s||''; if(!sm) continue;
      globalThis.__PRP_LASTBLK=null;
      let svg=null; try{ svg=molToSVG(sm); }catch(e){}
      if(!svg) continue; n++;
      const blk=globalThis.__PRP_LASTBLK; if(!blk) continue;
      let A; try{ A=analyse(blk); }catch(e){ continue; }
      if(!A.bad.length) continue;
      hits++;
      console.log('--- '+lib.slice(0,20)+' / '+(r.code||'?')+' / nA='+A.nA+' med='+A.med);
      console.log('    '+(r.n||'').slice(0,48));
      for(const b of A.bad)
        console.log('    键 '+b.a+'-'+b.b+' ('+b.ea+b.eb+') ty='+b.ty+' r='+b.r+' deg='+b.da+'/'+b.db+(b.ring?' [环内]':' [链]'));
      console.log('    SMILES '+sm);
    }
  }
  console.log('\n扫 '+n+' 张, 有缺陷 '+hits+' 张, '+((Date.now()-t0)/1000).toFixed(0)+'s');
});
