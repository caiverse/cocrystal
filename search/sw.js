/* 产物库检索站 — 分片缓存 Service Worker
 * 只缓存 data/products.*.json 与 data/fp.*.bin 这类内容不变的分片文件;
 * HTML/JS/manifest/stats/category_shards 一律走网络, 保证代码与元数据即时更新。
 * 缓存名带 CACHE_VER: 站点结构升级时改这个数字即可整体失效旧缓存。
 */
'use strict';
const CACHE_VER = 2;   // 2026-11-28: 33 个分片 d 字段统一, 旧缓存必须整体失效
const CACHE_NAME = 'cocrystal-shards-v' + CACHE_VER;
const DATA_RE = /\/data\/(products\.\d+\.json|fp\.\d+\.bin)$/;

self.addEventListener('install', (ev) => {
  self.skipWaiting();
});

self.addEventListener('activate', (ev) => {
  ev.waitUntil((async () => {
    const keys = await caches.keys();
    await Promise.all(keys.filter(k => k.startsWith('cocrystal-shards-') && k !== CACHE_NAME).map(k => caches.delete(k)));
    await self.clients.claim();
  })());
});

self.addEventListener('fetch', (ev) => {
  const req = ev.request;
  if (req.method !== 'GET') return;
  const url = new URL(req.url);
  if (url.origin !== self.location.origin) return;
  if (!DATA_RE.test(url.pathname)) return;   // 非分片: 不介入, 走浏览器默认缓存

  ev.respondWith((async () => {
    const cache = await caches.open(CACHE_NAME);
    const hit = await cache.match(req, { ignoreSearch: true });
    if (hit) return hit;                     // 分片内容不可变, cache-first 安全
    const res = await fetch(req);
    if (res.ok) cache.put(req, res.clone()); // clone 后再返回给页面
    return res;
  })());
});
