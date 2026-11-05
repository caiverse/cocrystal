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
try:
    from rdkit.Chem import rdDepictor   # 规范取向(整体不歪斜), 与检索站 JS 管线 straighten_depiction 对等
    _HAVE_STRAIGHTEN = True
except Exception:
    _HAVE_STRAIGHTEN = False

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

def _lengthen_terminal_bonds(mol):
    """端基(度=1)键保护: coordgen 偶发把端基键压到 0.8Å(中位 1.5), 画出来半截
    (检索站 JS 管线 _prpLengthenTerminalBonds 的同一根因与同一修法)。
    只放大不缩小(阈值 0.85×中位), 只沿原键方向 —— 不动其它任何原子、不改变任何键角。"""
    conf = mol.GetConformer()
    n = mol.GetNumAtoms()
    nb = [[] for _ in range(n)]
    ds = []
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        nb[i].append(j); nb[j].append(i)
        p1 = conf.GetAtomPosition(i); p2 = conf.GetAtomPosition(j)
        ds.append(math.hypot(p1.x - p2.x, p1.y - p2.y))
    if len(ds) < 2:
        return 0
    ds.sort()
    med = ds[len(ds) // 2]
    if med <= 0.2:
        return 0
    fixed = 0
    for k in range(n):
        if len(nb[k]) != 1:
            continue
        j = nb[k][0]
        pk = conf.GetAtomPosition(k); pj = conf.GetAtomPosition(j)
        d = math.hypot(pk.x - pj.x, pk.y - pj.y)
        if d < 1e-6 or d >= 0.85 * med:
            continue
        f = med / d
        conf.SetAtomPosition(k, Point3D(pj.x + (pk.x - pj.x) * f, pj.y + (pk.y - pj.y) * f, 0.0))
        fixed += 1
    return fixed

def normalize_pendant_rings(mol):
    """让所有"末端(叶)六元环"取向水平(0°/180°)、并把系链(连母核的键)拉长以缓解取代基拥挤；
    融合/螺/桥联母核与"中心环"绝不旋转。
    本实现与线上检索站 JS 路径 postRotatePendantRings 口径严格一致(同一套判定 + 同一套守卫),
    目的是让首页侧栏与检索站结构图行为统一、且杜绝 Z 形折弯/双旋转:
      - sharesAtom: 任一环原子属于 >=2 个环 -> 稠合/螺/桥联, 跳过(刚性母核)。
      - isLeaf(SUBST_MAX=5): 摘掉环原子后, 其余图若分裂出 >5 原子的大分量 -> 骨架被切断,
        非叶子环, 跳过(避免旋转中心环把骨架扯断)。
      - reachExcl 选系链: 锚点外连"可达原子数最大"者 = 真正连向母核的那根键; 取代基(CF3/NO2/OH)
        外连数小, 绝不会被误判为系链(旧版用 deg>1 启发式会误选取代基碳 -> 连链被纳入旋转 -> 双旋转 Z 弯)。
      - frag = 环原子 + 沿"非系链键" BFS 可达的取代基; 不含连链(连链不被旋转 -> 两个苯环各自
        绕自身 ipso 旋转, 连链端点(ipso)固定, 绝不会被二次旋转形成 Z 弯)。
      - 守卫: frag.size > max(8, floor(nAtoms*0.6)) -> 片断过大(连链把另一苯环也拉进来)则跳过,
        直接采用 coordgen 自然取向, 保持连通、不破坏主结构(与检索站完全一致)。
      - hasSubst: 仅当片断含取代基(片断 > 环)时才加长系链; 干净叶环不加长(避免产生断连空隙)。
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
    ringof = {}
    for ri, r in enumerate(rings):
        for a in r:
            ringof.setdefault(a, []).append(ri)
    conf = mol.GetConformer()
    def _xy(i):
        p = conf.GetAtomPosition(i); return p.x, p.y

    def shares_atom(ri):
        # 任一环原子属于 >=2 个环 -> 稠合/螺/桥联, 视为刚性母核不旋转
        return any(len(ringof.get(a, [])) > 1 for a in rings[ri])

    SUBST_MAX = 12
    # ★ 2026-11-05 对等升级(与检索站 app.js 同步): 5 -> 12 ——
    # "2-苄氧基/2-苯硫基/2-苯氧基苯基"等联芳基的取代基片段含 7~8 原子, 旧上限 5 会把这些悬挂环
    # 整个排除(连旋转都不做, 环歪 5~15°)。放宽到 12; 取代基里若还含另一悬挂环(嵌套),
    # 本实现按"外环先转、内环后转"的顺序逐个绕各自 ipso 旋转, 坐标级复合天然正确。
    def is_leaf(ring_atoms):
        # 摘掉环原子后, 其余图必须仍基本连通(只允许 <=SUBST_MAX 的小取代基碎片孤立)
        rem = [True] * n
        for a in ring_atoms:
            rem[a] = False
        rest = [i for i in range(n) if rem[i]]
        if not rest:
            return False
        visited = set(); comps = []
        for s0 in rest:
            if s0 in visited:
                continue
            comp = []; q = [s0]; visited.add(s0)
            while q:
                u = q.pop(); comp.append(u)
                for nb in adj[u]:
                    if rem[nb] and nb not in visited:
                        visited.add(nb); q.append(nb)
            comps.append(comp)
        if not comps:
            return False
        comps.sort(key=len, reverse=True)
        for k in range(1, len(comps)):
            if len(comps[k]) > SUBST_MAX:
                return False
        return True

    def reach_excl(nb, ca):
        # 从 nb 出发、不进入 ca(环原子)的可达原子数 -> 区分"连向母核"与"取代基"
        seen = {nb}; st = [nb]; cnt = 0
        while st:
            u = st.pop(); cnt += 1
            for v in adj[u]:
                if v in seen or v in ca:
                    continue
                seen.add(v); st.append(v)
        return cnt

    jobs = []
    for ri, ring in enumerate(rings):
        if shares_atom(ri):
            continue
        ring_atoms = set(ring)
        if not is_leaf(ring_atoms):
            continue
        # 系链 = 环原子 -> 非环原子的外部键中, 锚点外连可达最多者(连向母核)
        tether = None; best = -1
        for a in ring:
            for nb in adj[a]:
                if nb in ring_atoms:
                    continue
                sz = reach_excl(nb, ring_atoms)
                if sz > best:
                    best = sz; tether = (a, nb)
        if tether is None:
            continue
        ip, anchor = tether
        # 片断 = 环原子 + 沿非系链键 BFS 可达的取代基(不含连链)
        frag = set(ring_atoms); q = list(ring_atoms)
        while q:
            u = q.pop()
            for nb in adj[u]:
                if nb in frag:
                    continue
                if (u == ip and nb == anchor) or (u == anchor and nb == ip):
                    continue
                frag.add(nb); q.append(nb)
        if len(frag) > max(8, math.floor(n * 0.6)):
            continue  # 守卫: 片断过大(连链把另一苯环拉进来) -> 跳过, 不破坏主结构
        has_subst = len(frag) > len(ring_atoms)
        jobs.append((ri, ring_atoms, ip, anchor, frag, has_subst))

    changed = False
    F = 1.25  # 系链拉长系数(仅对含取代基者)
    non_frag_all = set(range(n))
    for ri, ring_atoms, ip, anchor, frag, has_subst in jobs:
        ix, iy = _xy(ip)
        ax, ay = _xy(anchor)
        # ★ 2026-11-05 对等升级(与检索站 JS 工序二同口径): 径向方向改用"ipso 处两条环键方向的角平分线"的反向。
        #   环心平均法在环上含杂原子(S/O/N)时会被其坐标带偏(检索站实测 123.4/116.5 而非 120/120);
        #   角平分线只依赖 ipso 与两个 ortho 三点, 与环是否正多边形、是否含杂原子都无关。
        #   正六边形下二者完全等价; 取不到三点才退回环心法。
        ring_l = list(rings[ri])
        ii = ring_l.index(ip)
        o1 = ring_l[(ii + 1) % len(ring_l)]; o2 = ring_l[(ii - 1) % len(ring_l)]
        o1x, o1y = _xy(o1); o2x, o2y = _xy(o2)
        l1 = math.hypot(o1x - ix, o1y - iy); l2 = math.hypot(o2x - ix, o2y - iy)
        cxr0 = sum(_xy(a)[0] for a in ring_atoms) / len(ring_atoms)
        cyr0 = sum(_xy(a)[1] for a in ring_atoms) / len(ring_atoms)
        if l1 > 1e-6 and l2 > 1e-6:
            vx = -((o1x - ix) / l1 + (o2x - ix) / l2)
            vy = -((o1y - iy) / l1 + (o2y - iy) / l2)
            if math.hypot(vx, vy) > 1e-6:
                rad_ang = math.atan2(vy, vx)           # 放射外向方向(角平分线反向)
            else:
                rad_ang = math.atan2(iy - cyr0, ix - cxr0)
        else:
            rad_ang = math.atan2(iy - cyr0, ix - cxr0)  # 环心法(退化兜底)
        anc_ang = math.atan2(ay - iy, ax - ix)         # 当前 anchor(系链)方向
        th0 = rad_ang - anc_ang
        # 在两种取向(th0=系链放射外向 与 th0+pi=系链指向环内)中选与母核/其余原子 clash 更少者 ->
        # 杜绝旋转把取代基(F/CF3/OH)甩到母核原子(N/碳)上造成重叠(典型: struct_13 的 F 与 N 重叠)。
        non_frag = non_frag_all - frag
        def _clash_for(th):
            # 只统计"两个都是杂原子(N/O/F/Cl/Br/S/P…)且中心距<0.5"的成对接近——这才是真正可见的
            # "标签互相压"(如 F 与 N/O 重叠)。碳不画标签, 碳-碳/杂原子-碳靠近不可见, 不计;
            # 阈值 0.5(约半根键长)对应标签文字实际重叠, 0.5~0.6 的"略近"不误判, 避免把环退回歪斜取向。
            c = 0
            sym_i_all = [mol.GetAtomWithIdx(i).GetSymbol() for i in frag]
            for idx_i, i in enumerate(frag):
                if sym_i_all[idx_i] == 'C':
                    continue
                x, y = _xy(i); dx = x - ix; dy = y - iy
                xi = ix + dx * math.cos(th) - dy * math.sin(th)
                yi = iy + dx * math.sin(th) + dy * math.cos(th)
                for j in non_frag:
                    if j in adj[i]:
                        continue
                    if mol.GetAtomWithIdx(j).GetSymbol() == 'C':
                        continue
                    pj = conf.GetAtomPosition(j)
                    if math.hypot(xi - pj.x, yi - pj.y) < 0.5:
                        c += 1
            return c
        c0 = _clash_for(th0); c1 = _clash_for(th0 + math.pi)
        # 若两种平顶取向都无法避免与母核 clash(如 ortho-F 直接连母核 N 的拓扑, 邻位取代基天然贴近母核),
        # 则放弃旋转、保留 coordgen 自然布局(零 clash)。clean 优先于 upright: 宁可略歪也不让标签重叠。
        if min(c0, c1) > 0:
            continue
        th = th0 if c0 <= c1 else th0 + math.pi
        cos_t, sin_t = math.cos(th), math.sin(th)
        for i in frag:
            x, y = _xy(i); dx = x - ix; dy = y - iy
            nxp = ix + dx * cos_t - dy * sin_t
            nyp = iy + dx * sin_t + dy * cos_t
            conf.SetAtomPosition(i, Point3D(nxp, nyp, 0.0))
        changed = True
        # 系链拉长: 仅含取代基者, 把片断沿"远离母核"方向平移(加长 ipso-anchor 那一根键, 同时进一步远离母核缓解拥挤)
        if has_subst:
            tx, ty = _xy(anchor)
            L0 = math.hypot(ix - tx, iy - ty)
            if L0 > 1e-6:
                ux, uy = (ix - tx) / L0, (iy - ty) / L0
                delta = (F - 1.0) * L0
                for i in frag:
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
    # ★ 2026-11-05 对等升级(与检索站 JS 管线同步, 用户反馈的同一批根因):
    # (1) 规范取向 —— coordgen 块坐标未做规范取向时整个分子斜着摆(检索站"整体歪斜"反馈);
    #     StraightenDepiction 只做整体刚性旋转, 实测逐键长一个数都不变, 不破坏环几何。
    if _HAVE_STRAIGHTEN:
        try:
            rdDepictor.StraightenDepiction(m)
        except Exception:
            pass
    # (2) 端基键保护 —— coordgen 偶发把端基(度=1)键压到 0.8Å 画成半截(检索站"C=O 一长一短"反馈)。
    try:
        _lengthen_terminal_bonds(m)
    except Exception:
        pass
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

def load_picked():
    """从 search/data/products.*.json 读取 FIXED_CODES 对应的 (cat, smiles, name, code)。"""
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
    return picked

def write_rails(picked):
    rails = []
    for i, (g, s, n, code) in enumerate(picked, 1):
        fn = "struct_%02d.svg" % i
        if render(s, os.path.join(OUTDIR, fn)):
            rails.append({"cat": g, "name": n or "", "code": code or "", "svg": "side-structures/" + fn})
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
    return rails

def main():
    picked = load_picked()
    rails = write_rails(picked)
    print("picked", len(picked), "rendered", len(rails))
    print("categories:", [p[0] for p in picked])
    for r in rails:
        print(" ", r["svg"], "|", r["cat"], "|", r["name"][:24], "|", r["code"])

if __name__ == "__main__":
    main()
