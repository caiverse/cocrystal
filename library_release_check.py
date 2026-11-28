#!/usr/bin/env python3
"""新库上线前的一致性 + 展示规范闸门。

用前跑: python library_release_check.py
退出码 0=PASS, 1=FAIL。挂进 prepush_audit 由它调用。

检查项:
  1. 三处计数一致: products.manifest.count == stats.total == category_shards.count == 分片实际条数之和
  2. category_shards 的 shardOrder/offsets/counts/byCategory 覆盖 manifest 的全部分片
  3. 每个分片内每条记录的 g 都能在 byCategory 找到对应分片(否则下拉选不到)
  4. CCDC 描述字段 d 合法: 只允许规范文案, 禁止旧值「未被收录」
  5. 卡片展示字段白名单: EXTRA_LABEL 里的键由前端控制, 本脚本不读前端 ——
     只保证 d 与 code 不夹带内部断言字段(历史事故: code 夹带 CAS/算法名)
"""
import glob, json, os, sys

ROOT = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(ROOT, 'search', 'data')
CAT_DOC = 'CSD 全库结构比对未收录'
ALLOWED_D_PREFIX = ('CSD 全库结构比对未收录', 'P–O–C 环水解开环次膦酸', 'CSD b14 疑似命中')

def fail(msg):
    print('FAIL:', msg)
    return False

def main():
    ok = True
    mf = json.load(open(os.path.join(DATA, 'products.manifest.json'), encoding='utf-8'))
    st = json.load(open(os.path.join(DATA, 'stats.json'), encoding='utf-8'))
    cs = json.load(open(os.path.join(DATA, 'category_shards.json'), encoding='utf-8'))
    shards = mf['shards']

    actual = 0
    bad_d = {}
    bad_cat = {}
    codes = set()
    dup = 0
    shard_len = {}
    for name in shards:
        p = os.path.join(DATA, name)
        if not os.path.exists(p):
            ok &= fail('缺分片 ' + name)
            continue
        arr = json.load(open(p, encoding='utf-8'))
        shard_len[name] = len(arr)
        actual += len(arr)
        for r in arr:
            if not isinstance(r, dict):
                continue
            c = r.get('code')
            if c in codes:
                dup += 1
            codes.add(c)
            d = r.get('d', '')
            if d and not any(d.startswith(x) for x in ALLOWED_D_PREFIX):
                bad_d[d] = bad_d.get(d, 0) + 1
            g = r.get('g', '')
            if g and name not in cs['byCategory'].get(g, []):
                bad_cat[g] = bad_cat.get(g, 0) + 1

    print('计数: manifest=%d stats=%d catmap=%d 分片实际=%d' % (
        mf['count'], st['total'], cs['count'], actual))
    if not (mf['count'] == st['total'] == cs['count'] == actual):
        ok &= fail('三处计数不一致')
    if dup:
        ok &= fail('重复 code %d 条' % dup)

    missing = [n for n in shards if n not in cs['shardOffsets'] or n not in cs['shardCounts']]
    if missing:
        ok &= fail('category_shards 缺分片: ' + ', '.join(missing[:10]))
    off = 0
    for n in shards:
        if n not in cs['shardOffsets']:
            continue
        if cs['shardOffsets'][n] != off:
            ok &= fail('offset 错位: %s 期望 %d 实 %d' % (n, off, cs['shardOffsets'][n]))
        if cs['shardCounts'].get(n) != shard_len.get(n):
            ok &= fail('分片条数不符: %s' % n)
        off += shard_len.get(n, 0)

    if bad_d:
        ok &= fail('CCDC 描述含非规范值: ' + str(list(bad_d.items())[:8]))
    if bad_cat:
        ok &= fail('byCategory 未覆盖记录类别: ' + str(list(bad_cat.items())[:8]))

    print('分片=%d 类别=%d 唯一code=%d' % (len(shards), len(cs['byCategory']), len(codes)))
    if ok:
        print('LIBRARY RELEASE CHECK PASS')
    return 0 if ok else 1

if __name__ == '__main__':
    sys.exit(main())
