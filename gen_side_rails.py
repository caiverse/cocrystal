# -*- coding: utf-8 -*-
"""生成首页两侧滚动结构推送所需的 12 个产物库结构 SVG + rails.js。
数据源：search/data/products.*.json（真实产物库，含 SMILES）。
渲染：RDKit 2D 结构图（矢量 SVG，清晰可缩放）。
产物：side-structures/struct_NN.svg  +  rails.js（数据+构建逻辑，纯前端、无 fetch 依赖）。
"""
import os, re, json, glob
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D

ROOT = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(ROOT, "side-structures")
os.makedirs(OUTDIR, exist_ok=True)

W, H = 300, 240

def render(smiles, path):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return False
    try:
        Chem.SanitizeMol(m)
    except Exception:
        return False
    opts = Draw.MolDrawOptions()
    opts.bondLineWidth = 1.8
    opts.minFontSize = 15
    opts.maxFontSize = 24
    opts.useBWAtomPalette = False
    opts.scalingFactor = 1.05
    d = rdMolDraw2D.MolDraw2DSVG(W, H)
    d.SetDrawOptions(opts)
    d.DrawMolecule(m)
    d.FinishDrawing()
    svg = d.GetDrawingText()
    # 规范化：去掉固定宽高、加 viewBox 让其可被 CSS 自由缩放
    mopen = re.search(r"<svg([^>]*)>", svg)
    attrs = mopen.group(1)
    attrs = re.sub(r"\swidth\s*=\s*['\"][^'\"]*['\"]", "", attrs)
    attrs = re.sub(r"\sheight\s*=\s*['\"][^'\"]*['\"]", "", attrs)
    attrs = re.sub(r"\sviewBox\s*=\s*['\"][^'\"]*['\"]", "", attrs)
    newopen = '<svg%s viewBox="0 0 %d %d" width="%d" height="%d" preserveAspectRatio="xMidYMid meet">' % (attrs, W, H, W, H)
    svg = svg[:mopen.start()] + newopen + svg[mopen.end():]
    with open(path, "w", encoding="utf-8") as f:
        f.write(svg)
    return True

# ---- 分层抽样：每个类别取 1~2 个「中等大小、结构好看」的候选 ----
shards = sorted(glob.glob(os.path.join(ROOT, "search", "data", "products.*.json")))[:5]
cat_pool = {}   # cat -> list of (hatoms, smiles, name, code)
cat_count = {}
for sh in shards:
    try:
        data = json.load(open(sh, encoding="utf-8"))
    except Exception:
        continue
    for it in data:
        s = it.get("s")
        g = (it.get("g") or "其他").strip()
        if not s or len(s) > 160:
            continue
        cat_count[g] = cat_count.get(g, 0) + 1
        if g not in cat_pool:
            cat_pool[g] = []
        # 延迟解析：仅记录，渲染时再验证
        cat_pool[g].append((s, it.get("n", ""), it.get("code", "")))

# 选 12 个：优先多样类别；每类取重原子数在 [14,36] 之间的较好看样本
def hatoms(s):
    m = Chem.MolFromSmiles(s)
    return m.GetNumHeavyAtoms() if m else 0

picked = []
used_cats = []
# 第一轮：每类 1 个（取该类中重原子数最靠近 22 的）
for g in sorted(cat_pool, key=lambda k: -cat_count.get(k, 0)):
    if len(picked) >= 16:
        break
    best = None; bestscore = 1e9
    for (s, n, code) in cat_pool[g]:
        h = hatoms(s)
        if h < 12 or h > 38:
            continue
        sc = abs(h - 22)
        if sc < bestscore:
            bestscore = sc; best = (s, n, code)
    if best:
        picked.append((g,) + best); used_cats.append(g)

# 第二轮：若不足 12，再从已有类别补，凑满
if len(picked) < 16:
    for g in used_cats:
        if len(picked) >= 16:
            break
        for (s, n, code) in cat_pool[g]:
            if (g, s, n, code) in [(p[0], p[1], p[2], p[3]) for p in picked]:
                continue
            h = hatoms(s)
            if 12 <= h <= 38:
                picked.append((g, s, n, code)); break

# 渲染
rails = []
for i, (g, s, n, code) in enumerate(picked, 1):
    fn = "struct_%02d.svg" % i
    if render(s, os.path.join(OUTDIR, fn)):
        rails.append({"cat": g, "name": n or "", "code": code or "", "svg": "side-structures/" + fn})

# 写 rails.js（数据内联，无 fetch 依赖）
js = "// 首页两侧滚动结构推送（产物库真实结构）。由 gen_side_rails.py 生成。\n"
js += "(function(){\n  var RAILS = " + json.dumps(rails, ensure_ascii=False) + ";\n"
js += """
  function card(d){
    var a=document.createElement('a');
    a.className='rail-card'; a.href='search/';
    a.title=(d.name||'')+'  ('+(d.code||'')+')';
    a.innerHTML='<img loading="lazy" src="'+d.svg+'" alt="'+(d.name||'')+'">'
      +'<span class="rc-cat">'+(d.cat||'')+'</span>'
      +'<span class="rc-name">'+(d.name||'')+'</span>'
      +'<span class="rc-code">'+(d.code||'')+'</span>';
    return a;
  }
  var HALF=Math.ceil(RAILS.length/2);
  var SETS={left:RAILS.slice(0,HALF), right:RAILS.slice(HALF)};
  ['left','right'].forEach(function(side){
    var rail=document.querySelector('.side-rail.'+side); if(!rail) return;
    var set=SETS[side];
    var track=document.createElement('div'); track.className='rail-track';
    set.forEach(function(d){ track.appendChild(card(d)); });   // 第一份
    set.forEach(function(d){ track.appendChild(card(d)); });   // 第二份（无缝循环）
    rail.appendChild(track);
  });
})();
"""
with open(os.path.join(ROOT, "rails.js"), "w", encoding="utf-8") as f:
    f.write(js)

print("picked", len(picked), "rendered", len(rails))
print("categories:", [p[0] for p in picked])
for r in rails:
    print(" ", r["svg"], "|", r["cat"], "|", r["name"][:24], "|", r["code"])
