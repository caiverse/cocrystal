#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""上线前强制闸门：新分片产品代码 + 命名 + 索引一致性 + 渲染质量。FAIL 即禁止 push。

用法：
  python prepush_audit.py                    # 自动取 git 改动的分片（含未跟踪新分片）
  python prepush_audit.py --shards 52 53 54  # 指定分片
  python prepush_audit.py --all              # 全库语义核验（慢）
  python prepush_audit.py --full-render      # 渲染闸门逐条全量（默认每片抽 1200 条）

渲染闸门（2026-11-26 起强制，调 render_gate.js）：用线上同一个渲染器实渲抽样，判
键线交叉=0 / 键长比 0.70~1.65 / 非键最近距离≥0.40×中位 / 环内角平均偏差≤8°。
「数据对但画错」的库一律不许上传；没有跳过开关。

退出码：0=PASS，1=FAIL。
"""
import argparse, hashlib, glob, json, os, re, subprocess, sys

from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')

_B62 = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'
_CODE_SPACE = 62 ** 9


def make_code(inchikey, name='', formula=''):
    """与 product_search_build/aggregate.py::make_code() 完全一致。"""
    ik = (inchikey or '').strip() or ('NF:' + (name or '') + '|' + (formula or ''))
    h = int.from_bytes(hashlib.sha256(ik.encode('utf-8')).digest(), 'big') % _CODE_SPACE
    s = ''
    n = h
    if n == 0:
        return 'CC-0'
    while n:
        s = _B62[n % 62] + s
        n //= 62
    return 'CC-' + s.zfill(9)


def inchikey_from_smiles(smiles):
    try:
        m = Chem.MolFromSmiles(smiles)
        return Chem.MolToInchiKey(m) if m is not None else ''
    except Exception:
        return ''


REPO = 'C:/Users/bc/cocrystal'
DATA = REPO + '/search/data'
SHARD_RE = re.compile(r'products\.(\d+)\.json$')
CODE_RE = re.compile(r'^CC-[0-9A-Za-z]{9}(?:-\d+)?$')
# 名称中不允许出现原料/编号尾巴；命名必须是产物结构本身的名字
BAD_NAME_RE = re.compile(r'CAS\s*\d|\[aldehyde\s*:|\[anhydride\s*:|\(diphenyl(?:phosphinic|chloro)|from aromatic (?:amine|diamine)', re.I)


def changed_shards():
    out = subprocess.run(['git', '-C', REPO, 'diff', 'HEAD', '--name-only'],
                         capture_output=True, text=True, check=True).stdout.splitlines()
    out += subprocess.run(['git', '-C', REPO, 'ls-files', '--others', '--exclude-standard'],
                          capture_output=True, text=True, check=True).stdout.splitlines()
    return sorted({int(m.group(1)) for f in out if (m := SHARD_RE.search(f))})


def shard_path(i):
    return os.path.join(DATA, 'products.%d.json' % i)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--shards', type=int, nargs='*')
    ap.add_argument('--all', action='store_true')
    ap.add_argument('--full-render', action='store_true',
                    help='渲染闸门逐条全量(默认每片抽 1200 条)')
    a = ap.parse_args()

    fails = []
    if a.all:
        shards = [int(m.group(1)) for p in glob.glob(DATA + '/products.*.json')
                  if (m := SHARD_RE.search(p))]
    elif a.shards is not None:
        shards = a.shards
    else:
        shards = changed_shards()
    print('待审分片:', shards if shards else '(无改动，仅做全库索引/唯一性检查)')

    # 全库 code 唯一性（便宜，必做）
    codes = {}
    counts = []
    for p in sorted(glob.glob(DATA + '/products.*.json')):
        if not SHARD_RE.search(p):
            continue
        rs = json.load(open(p, encoding='utf-8'))
        counts.append(len(rs))
        for r in rs:
            c = r.get('code', '')
            if c in codes:
                fails.append('重复 code %s (分片 %s / %s)' % (c, codes[c], p))
            codes[c] = p
    print('全库记录 %d | 唯一 code %d' % (sum(counts), len(codes)))

    for i in shards:
        path = shard_path(i)
        rs = json.load(open(path, encoding='utf-8'))
        bad_shape = bad_code = bad_name = 0
        samples = []
        for r in rs:
            c = r.get('code', '')
            if not CODE_RE.match(c):
                bad_shape += 1
                if len(samples) < 3: samples.append(('CODE-SHAPE', c, r.get('n', '')))
            ik = inchikey_from_smiles(r.get('s', ''))
            if not ik:
                fails.append('分片 %d: InChIKey 计算失败 %s' % (i, r.get('s', '')[:60]))
                continue
            base = make_code(ik)
            if c != base and not c.startswith(base + '-'):
                bad_code += 1
                if len(samples) < 3: samples.append(('CODE-SEMANTIC', c, base))
            for key in ('n', 'e'):
                if BAD_NAME_RE.search(r.get(key, '')):
                    bad_name += 1
                    if len(samples) < 3: samples.append(('NAME-' + key, r.get(key, '')[:100], ''))
        smis, ens = {r['s'] for r in rs}, {r.get('e', '') for r in rs}
        dup_name = len(ens) != len(smis)
        print('分片 %d: %d 条 | 码形异常 %d | 码语义异常 %d | 名称夹带原料/CAS %d | 英文名与结构非一一 %s'
              % (i, len(rs), bad_shape, bad_code, bad_name, dup_name))
        for s in samples: print('   ', s)
        if bad_shape or bad_code or bad_name or dup_name:
            fails.append('分片 %d 未通过（见上）' % i)

    man = json.load(open(DATA + '/products.manifest.json', encoding='utf-8'))
    fpc = json.load(open(DATA + '/fp_chunks.json', encoding='utf-8'))
    st = json.load(open(DATA + '/stats.json', encoding='utf-8'))
    idx_ok = (man['count'] == fpc['numRecords'] == st['total'] == sum(counts)
              and len(man['shards']) == len(counts))
    print('索引对账 manifest/fp/stats/分片合计: %s (%d/%d/%d/%d)'
          % ('PASS' if idx_ok else 'FAIL', man['count'], fpc['numRecords'], st['total'], sum(counts)))
    if not idx_ok:
        fails.append('索引不一致')

    # ---- 渲染闸门（2026-11-26 起强制）----
    # 起因: 用户发现卡片上异丙基的键穿过苯并噻唑环, 而当时没有任何自动检查能拦住这种"数据对但画错"的产物。
    # 口径: 待审分片逐片抽样, 用线上同一个渲染器(search/app.js::molToSVG)实渲, 再用独立实现的几何度量判:
    #   键线交叉=0 / 键长比 0.70~1.65 / 非键最近距离≥0.40×中位 / 环内角平均偏差≤8°。
    # 任一不过 → FAIL → 禁止 push。**这是硬闸门, 没有跳过开关**: 画不出来的库不许上传。
    if shards:
        gate = os.path.join(REPO, 'render_gate.js')
        if not os.path.exists(gate):
            fails.append('缺少 render_gate.js，渲染闸门无法执行')
        else:
            cmd = ['node', gate, '--shards'] + [str(i) for i in shards]
            if a.full_render:
                cmd.append('--full')
            print('\n[渲染闸门] ' + ' '.join(cmd))
            try:
                rc = subprocess.run(cmd, cwd=REPO).returncode
            except FileNotFoundError:
                fails.append('找不到 node，渲染闸门无法执行（不许绕过）')
            else:
                if rc != 0:
                    fails.append('渲染闸门未通过（见上），本批不许上传')

    if fails:
        print('\nPREPUSH AUDIT FAIL:')
        for f in fails: print(' -', f)
        return 1
    print('\nPREPUSH AUDIT PASS')
    return 0


if __name__ == '__main__':
    sys.exit(main())
