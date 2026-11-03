// 首页两侧滚动结构推送（产物库真实结构）。由 gen_side_rails.py 生成。
(function(){
  var RAILS = [{"cat": "席夫碱", "name": "N-(4-(甲基氨基)苯基)亚甲基-2,6-二氟苯胺", "code": "CC-vVnLAZjZu", "svg": "side-structures/struct_01.svg"}, {"cat": "N-取代地平（二甲酯）", "name": "N-五氟苯基-2,6-二甲基-4-(4-溴-1H-吡咯-2-基)-1,4-二氢吡啶-3,5-二甲酸二甲酯", "code": "CC-oj8Z9bP5j", "svg": "side-structures/struct_02.svg"}, {"cat": "胺取代吖啶二酮", "name": "N-(4-氟苯基)-9-(3-甲基-1-苯并噻吩-2-基)-1,2,3,4,5,6,7,8-八氢吖啶-1,8-二酮", "code": "CC-tEZcKNdpl", "svg": "side-structures/struct_03.svg"}, {"cat": "酰腙", "name": "4-氰基苯甲酰肼缩4-(2-甲氧基乙基)苯甲醛酰腙", "code": "CC-zvKJjGDWc", "svg": "side-structures/struct_04.svg"}, {"cat": "腙", "name": "N'-(5-氟-2-甲基苯基)-N-(3,5-二溴水杨亚甲基)肼", "code": "CC-LrKGHs3nf", "svg": "side-structures/struct_05.svg"}, {"cat": "酞嗪二酮", "name": "5,6,7,8-四溴-2-(1,3-苯并噁唑-2-基)-2,3-二氢酞嗪-1,4-二酮", "code": "CC-DftVcv83k", "svg": "side-structures/struct_06.svg"}, {"cat": "2-取代-异吲哚-1,3-二酮", "name": "N-(2,5-二甲基苯基)-1H-异吲哚-1,3(2H)-二酮", "code": "CC-dLgVqcKAW", "svg": "side-structures/struct_07.svg"}, {"cat": "席夫碱", "name": "N-(2-甲氧基萘-1-基亚甲基)-4-氟苯胺", "code": "CC-p288LoMDn", "svg": "side-structures/struct_08.svg"}, {"cat": "N-取代地平（二甲酯）", "name": "N-(萘-2-基)-2,6-二甲基-4-(1-苯并噻吩-2-基)-1,4-二氢吡啶-3,5-二甲酸二甲酯", "code": "CC-fEGseYO7z", "svg": "side-structures/struct_09.svg"}, {"cat": "胺取代吖啶二酮", "name": "N-(3-甲基苯基)-9-(2,3-二羟基-4-甲氧基苯基)-1,2,3,4,5,6,7,8-八氢吖啶-1,8-二酮", "code": "CC-wj1RSyAyw", "svg": "side-structures/struct_10.svg"}, {"cat": "酰腙", "name": "3-(三氟甲基)苯甲酰肼缩2-氰基苯甲醛酰腙", "code": "CC-geL90DXKE", "svg": "side-structures/struct_11.svg"}, {"cat": "腙", "name": "N'-(2-溴-5-氟苯基)-N-(2-溴-6-羟基苯基亚甲基)肼", "code": "CC-Mg01GHBlL", "svg": "side-structures/struct_12.svg"}, {"cat": "酞嗪二酮", "name": "2-(5-氟-2-甲基苯基)-2,3-二氢吡啶并[2,3-d]哒嗪-1,4-二酮", "code": "CC-4ZesQhRv6", "svg": "side-structures/struct_13.svg"}, {"cat": "2-取代-异吲哚-1,3-二酮", "name": "N-(苯甲酰肼)-1H-异吲哚-1,3(2H)-二酮", "code": "CC-rSu8HZhUx", "svg": "side-structures/struct_14.svg"}];

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
