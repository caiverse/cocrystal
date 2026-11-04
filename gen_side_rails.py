# -*- coding: utf-8 -*-
"""生成首页两侧滚动结构推送所需的 12 个产物库结构 SVG + rails.js。
数据源：search/data/products.*.json（真实产物库，含 SMILES）。
渲染：RDKit 2D 结构图（矢量 SVG，清晰可缩放）。
产物：side-structures/struct_NN.svg  +  rails.js（数据+构建逻辑，纯前端、无 fetch 依赖）。
"""
import os, re, json, glob, math
from rdkit import Chem
from rdkit.Chem import Draw, rdCoordGen
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry import Point3D

ROOT = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(ROOT, "side-structures")
os.makedirs(OUTDIR, exist_ok=True)

W, H = 300, 240

# ---------- 悬挂(末端)六元环取向规范化: 让所有末端苯环水平(0°/180°), 视觉统一"正" ----------
# 与线上/离线检索站同一逻辑: RDKit/coordGen 不对单点连接的悬挂环做取向规范化, 导致末端环随取代基
# 歪斜。此处直接旋转坐标(比 SVG 后处理更干净): 识别"末端六元环"(仅一个母核连接点、且不与其它
# 六元环融合), 将其绕 ipso 旋转到水平; 融合母核(吖啶/萘/蒽等)绝不旋转。叶取代基(F/CH3…)随环一起转。
def _rings6(mol):
    return [tuple(r) for r in mol.GetRingInfo().AtomRings() if len(r) == 6]

def _is_fused(ring, rings):
    rs = set(ring); n = len(ring)
    for o in rings:
        if o == ring:
            continue
        os = set(o)
        for i in range(n):
            a, b = ring[i], ring[(i + 1) % n]
            if a in os and b in os:
                return True
    return False

def _bfs_dist(adj, rs, start):
    d = {start: 0}; q = [start]
    while q:
        c = q.pop(0)
        for nb in adj[c]:
            if nb in rs and nb not in d:
                d[nb] = d[c] + 1; q.append(nb)
    return d

def normalize_pendant_rings(mol):
    """让所有"末端(叶)六元环"(仅一个母核连接点、且不与其它六元环融合)取向水平(0°/180°)，
    并把系链(与母核相连的键)拉长以缓解取代基拥挤；融合母核(吖啶/萘/蒽等)绝不旋转。
    关键修正(对比旧版):
      1) 旋转集合 S_rot 只含本环 + 直接取代基, 绝不纳入"另一个环"的原子 —— 旧版会把嵌套的
         子悬挂环的连结点一起带走, 导致子环被拉歪(出现 114.9°/138.8°/-43.1° 等斜角与变形)。
      2) 系链拉长改为"整体平移"本悬挂片断(S_trans, 含嵌套子环, 作为刚体), 只拉伸与母核相连的
         那一根键, 片断内部键与子环取向完全不变。
      3) 处理顺序 root 先于 leaf, 保证嵌套子环随父环整体平移一次, 不会被重复平移。
    坐标层面旋转 => 原子标签(由 RDKit 按坐标绘制)始终随原子走, 不会脱附/漂浮。
    """
    rings = _rings6(mol)
    if not rings or mol.GetConformer() is None:
        return False
    n = mol.GetNumAtoms()
    adj = [[] for _ in range(n)]
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        adj[i].append(j); adj[j].append(i)
    deg = [len(adj[i]) for i in range(n)]
    ringof = {}
    for ri, r in enumerate(rings):
        for a in r:
            ringof.setdefault(a, []).append(ri)
    conf = mol.GetConformer()
    def _xy(i):
        p = conf.GetAtomPosition(i); return p.x, p.y
    pend = []
    for ri, ring in enumerate(rings):
        if _is_fused(ring, rings):
            continue
        rs = set(ring)
        ext = [nb for a in ring for nb in adj[a] if nb not in rs]
        if len(ext) < 1:
            continue
        ipso = max(ring, key=lambda a: sum(1 for nb in adj[a] if nb not in rs))
        tether = next((x for x in ext if deg[x] > 1), ext[0])
        pend.append((ri, ring, ipso, tether))
    def _is_leaf(t):
        return any(t[3] in o for o in rings)
    order = sorted(pend, key=lambda t: 1 if _is_leaf(t) else 0)  # 先处理 root(连母核), 后 leaf
    changed = False
    F = 1.25  # 系链拉长系数
    for ri, ring, ipso, tether in order:
        rs = set(ring)
        # S_trans / S_rot: 去掉系链键后的整个连通片断(本环 + 全部取代基支链 + 嵌套子环),
        # 作为刚体一起旋转/平移。这样多原子取代基(如 -CF3 / -COOMe / -OEt 的末端 O/F)不会
        # 被落在原地造成脱附/残线(旧版 S_rot 只含"直接邻原子", 漏掉支链末端原子 -> 浮空黑点)。
        visited = set(ring); stack = list(ring)
        while stack:
            a = stack.pop()
            for nb in adj[a]:
                if (a == ipso and nb == tether) or (a == tether and nb == ipso):
                    continue
                if nb not in visited:
                    visited.add(nb); stack.append(nb)
        Str = visited
        Srot = Str
        ix, iy = _xy(ipso)
        # 取向: 对齐"ipso -> 相邻环原子"这条环边到水平(0/180), 而非直径 -> 所有苯环统一平顶(正向),
        # 与线上检索站 JS 路径(postRotatePendantRings 阶段A)口径一致。对齐直径会时而平顶、时而尖顶
        # (尖顶向左/右看起来就"歪 30°"), 正是用户此前看到的"斜"。
        np_a = next((nb for nb in adj[ipso] if nb in rs), None)
        if np_a is not None:
            nx, ny = _xy(np_a); tx2, ty2 = _xy(tether)
            cur = math.atan2(ny - iy, nx - ix) * 180.0 / math.pi
            def _norm(x):
                while x > 180: x -= 360
                while x < -180: x += 360
                return x
            target = 0.0 if abs(_norm(0 - cur)) <= abs(_norm(180 - cur)) else 180.0
            th = (target - cur) * math.pi / 180.0
            # 翻转判定: 旋转后环质心应在 anchor(系链)反方向; 否则翻转 180°, 使环体不压到母核上。
            def _rp(i, ang):
                x, y = _xy(i); dx = x - ix; dy = y - iy
                return ix + dx * math.cos(ang) - dy * math.sin(ang), iy + dx * math.sin(ang) + dy * math.cos(ang)
            cxs = [_rp(a, th)[0] for a in ring]; cys = [_rp(a, th)[1] for a in ring]
            cxr = sum(cxs) / len(cxs) - ix; cyr = sum(cys) / len(cys) - iy
            adx = tx2 - ix; ady = ty2 - iy
            if cxr * adx + cyr * ady > 0:
                th += math.pi
            cos_t, sin_t = math.cos(th), math.sin(th)
            for i in Srot:
                x, y = _xy(i); dx = x - ix; dy = y - iy
                nx = ix + dx * cos_t - dy * sin_t
                ny = iy + dx * sin_t + dy * cos_t
                conf.SetAtomPosition(i, Point3D(nx, ny, 0.0))
            changed = True
        # 系链拉长: 把整个悬挂片断沿"远离母核"方向平移
        tx, ty = _xy(tether)
        L0 = math.hypot(ix - tx, iy - ty)
        if L0 > 1e-6:
            ux, uy = (ix - tx) / L0, (iy - ty) / L0
            delta = (F - 1.0) * L0
            for i in Str:
                x, y = _xy(i)
                conf.SetAtomPosition(i, Point3D(x + ux * delta, y + uy * delta, 0.0))
            changed = True
    return changed

def render(smiles, path):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return False
    try:
        Chem.SanitizeMol(m)
    except Exception:
        return False
    # 新版 2D 坐标生成器：布局更规范、取向更统一（避免旧版默认坐标的随机/歪斜取向）
    rdCoordGen.AddCoords(m)
    if m.GetConformer() is None:
        return False
    # 悬挂(末端)六元环取向规范化: 让所有末端苯环水平(0°/180°), 视觉统一"正"
    normalize_pendant_rings(m)
    # 归一化：居中对齐 + 统一缩放，使所有结构重心居中、视觉大小一致、不贴边
    conf = m.GetConformer()
    xs = [conf.GetAtomPosition(i).x for i in range(m.GetNumAtoms())]
    ys = [conf.GetAtomPosition(i).y for i in range(m.GetNumAtoms())]
    if not xs:
        return False
    minx, maxx, miny, maxy = min(xs), max(xs), min(ys), max(ys)
    w, h = maxx - minx, maxy - miny
    if max(w, h) < 1e-6:
        return False
    cx, cy = (minx + maxx) / 2.0, (miny + maxy) / 2.0
    scale = 200.0 / max(w, h)          # 统一归一化到约 200 坐标单位
    for i in range(m.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, Point3D((p.x - cx) * scale, (p.y - cy) * scale, 0.0))
    opts = Draw.MolDrawOptions()
    opts.bondLineWidth = 1.8
    opts.minFontSize = 15
    opts.maxFontSize = 24
    opts.useBWAtomPalette = False
    opts.scalingFactor = 1.0
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

# ---- 固定 14 条产物库结构（与既有 rails.js 完全一致），仅用修正后的渲染重画 ----
# 目的：用户已就这 14 个具体结构反馈"滚动栏结构乱"，故保持产物不变、只修正绘制。
FIXED_CODES = [
    "CC-2I8zxxJmA", "CC-ELLk2CRoi", "CC-bEuBVoQFS", "CC-mRx17fOXH", "CC-HR8BrrnY5",
    "CC-GF1TLL42I", "CC-HuU2g06cT", "CC-FVEwrOq5u", "CC-pFKxJ39xQ", "CC-qnIslewtq",
    "CC-vR1Vdc5SI", "CC-bGmr8GZBT", "CC-NqwQHOIFV", "CC-5EzjVe4QH",
]
rec_by_code = {}
for sh in glob.glob(os.path.join(ROOT, "search", "data", "products.*.json")):
    if sh.endswith("manifest.json"):
        continue
    try:
        data = json.load(open(sh, encoding="utf-8"))
    except Exception:
        continue
    if isinstance(data, list):
        for it in data:
            if isinstance(it, dict) and it.get("code") in FIXED_CODES:
                rec_by_code[it["code"]] = it

picked = []
for code in FIXED_CODES:
    r = rec_by_code.get(code)
    if r and r.get("s"):
        picked.append((r.get("g", "其他"), r["s"], r.get("n", ""), code))

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
