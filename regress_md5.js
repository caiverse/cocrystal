'use strict';
/* 新旧 app.js 逐分子 SVG md5 回归: 用法 node regress_md5.js <appjs> <outfile> [step]
 * 与其它 node 脚本隔离: 每个进程只加载一个 app.js。 */
const fs = require('fs'), path = require('path'), crypto = require('crypto');
const APPJS = process.argv[2];
const OUT = process.argv[3];
const STEP = +(process.argv[4] || 400);
const ROOT = __dirname, SEARCH = path.join(ROOT, 'search'), DATA = path.join(SEARCH, 'data');
globalThis.window = globalThis; globalThis.addEventListener = () => {};
window.addEventListener = () => {};
window.document = { getElementById: () => null, createElement: () => ({ style: {} }), addEventListener: () => {} };
let src = fs.readFileSync(APPJS, 'utf8').replace(/'use strict';/, '');
src += '\n;globalThis.__setRDKit=(m)=>{RDKit=m;};';
(0, eval)(src);
const init = require(path.join(SEARCH, 'lib/rdkit/RDKit_minimal.js')).default;
init({ locateFile: f => path.join(SEARCH, 'lib/rdkit', f) }).then(mod => {
  globalThis.__setRDKit(mod); mod.prefer_coordgen(true);
  const files = fs.readdirSync(DATA).filter(f => /^products\.\d+\.json$/.test(f))
    .sort((a,b) => +a.match(/\d+/)[0] - +b.match(/\d+/)[0]);
  const out = [];
  for (const f of files){
    const arr = JSON.parse(fs.readFileSync(path.join(DATA, f), 'utf8'));
    for (let i = 0; i < arr.length; i += STEP){
      const it = arr[i]; if (!it || !it.s) continue;
      let svg = null; try { svg = molToSVG(it.s, 240, 190); } catch(e){}
      const md5 = svg ? crypto.createHash('md5').update(svg).digest('hex') : 'NULL';
      out.push(it.code + '\t' + md5);
    }
  }
  fs.writeFileSync(OUT, out.join('\n'));
  console.log('wrote', out.length, '->', OUT);
}).catch(e => { console.error('FATAL', e && e.stack); process.exit(2); });
