// 首页两侧滚动结构推送（产物库真实结构）。由 gen_side_rails.py 生成。
(function(){
  var RAILS = [{"cat": "席夫碱", "name": "N-(5-氯-2,4-二羟基苯基)亚甲基-3-甲基苯胺", "code": "CC-2I8zxxJmA", "svg": "side-structures/struct_01.svg"}, {"cat": "N-取代地平（二甲酯）", "name": "N-(3-(三氟甲基)苯基)-2,6-二甲基-4-(5-乙氧基呋喃-2-基)-1,4-二氢吡啶-3,5-二甲酸二甲酯", "code": "CC-ELLk2CRoi", "svg": "side-structures/struct_02.svg"}, {"cat": "胺取代吖啶二酮", "name": "N-(3,5-二溴苯基)-9-(4-乙酰苯基)-1,2,3,4,5,6,7,8-八氢吖啶-1,8-二酮", "code": "CC-bEuBVoQFS", "svg": "side-structures/struct_03.svg"}, {"cat": "酰腙", "name": "水杨酰肼缩4-溴噻吩-2-甲醛酰腙", "code": "CC-mRx17fOXH", "svg": "side-structures/struct_04.svg"}, {"cat": "腙", "name": "N'-([trifluorometh]肼（待核实）)-N-(3-羧基苯甲亚甲基)肼", "code": "CC-HR8BrrnY5", "svg": "side-structures/struct_05.svg"}, {"cat": "酞嗪二酮", "name": "5,6,7,8-四氯-2-(4-乙基苯基)-2,3-二氢酞嗪-1,4-二酮", "code": "CC-GF1TLL42I", "svg": "side-structures/struct_06.svg"}, {"cat": "2-取代-异吲哚-1,3-二酮", "name": "2-(2-氯苯甲酰基)-2,3-二氢吡啶并[2,3-d]哒嗪-1,4-二酮", "code": "CC-HuU2g06cT", "svg": "side-structures/struct_07.svg"}, {"cat": "席夫碱", "name": "N-(3-甲基呋喃-2-基亚甲基)-2,4-二甲基苯胺", "code": "CC-FVEwrOq5u", "svg": "side-structures/struct_08.svg"}, {"cat": "N-取代地平（二甲酯）", "name": "N-(2-羟基苯基)-2,6-二甲基-4-(5-甲基噻吩-3-基)-1,4-二氢吡啶-3,5-二甲酸二甲酯", "code": "CC-pFKxJ39xQ", "svg": "side-structures/struct_09.svg"}, {"cat": "胺取代吖啶二酮", "name": "N-(2,3-二氟苯基)-9-(1,2,4-噻二唑-3-基)-1,2,3,4,5,6,7,8-八氢吖啶-1,8-二酮", "code": "CC-qnIslewtq", "svg": "side-structures/struct_10.svg"}, {"cat": "酰腙", "name": "3,5-二硝基苯甲酰肼缩3-(4-氯苯基)苯甲醛酰腙", "code": "CC-vR1Vdc5SI", "svg": "side-structures/struct_11.svg"}, {"cat": "腙", "name": "N'-(2-硝基苯基)-N-(3-溴柳亚甲基)肼", "code": "CC-bGmr8GZBT", "svg": "side-structures/struct_12.svg"}, {"cat": "酞嗪二酮", "name": "6-(苯乙炔基)-2-(4-溴-2-氟苯基)-2,3-二氢酞嗪-1,4-二酮", "code": "CC-NqwQHOIFV", "svg": "side-structures/struct_13.svg"}, {"cat": "2-取代-异吲哚-1,3-二酮", "name": "2-(2-甲氧基苯甲酰基)-2,3-二氢吡啶并[2,3-d]哒嗪-1,4-二酮", "code": "CC-5EzjVe4QH", "svg": "side-structures/struct_14.svg"}];

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
