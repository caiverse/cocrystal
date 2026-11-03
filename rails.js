// 首页两侧滚动结构推送（产物库真实结构）。由 gen_side_rails.py 生成。
(function(){
  var RAILS = [{"cat": "席夫碱", "name": "4-[4-羟基喹啉-3-基亚甲基氨基]苯甲酸", "code": "CC-uIOmuC1ns", "svg": "side-structures/struct_01.svg"}, {"cat": "N-取代地平（二甲酯）", "name": "N-甲基-2,6-二甲基-4-(1H-吡咯-2-基)-1,4-二氢吡啶-3,5-二甲酸二甲酯", "code": "CC-9OoxWYIUS", "svg": "side-structures/struct_02.svg"}, {"cat": "胺取代吖啶二酮", "name": "N-甲基-9-(2H-三唑-4-基)-1,2,3,4,5,6,7,8-八氢吖啶-1,8-二酮", "code": "CC-S9kJglEL2", "svg": "side-structures/struct_03.svg"}, {"cat": "酰腙", "name": "1-萘甲酰肼缩1-乙基-1H-吡咯-2-甲醛酰腙", "code": "CC-ReNXyBl42", "svg": "side-structures/struct_04.svg"}, {"cat": "腙", "name": "N'-(待核实)-N-(4-喹啉甲亚甲基)肼", "code": "CC-L2dlrScDw", "svg": "side-structures/struct_05.svg"}, {"cat": "酞嗪二酮", "name": "2-(喹啉-6-基)-2,3-二氢酞嗪-1,4-二酮", "code": "CC-t5hhuv8Gy", "svg": "side-structures/struct_06.svg"}, {"cat": "2-取代-异吲哚-1,3-二酮", "name": "N-(4-氯-3-(三氟甲基)苯基)-1H-异吲哚-1,3(2H)-二酮", "code": "CC-6nicDhKig", "svg": "side-structures/struct_07.svg"}, {"cat": "席夫碱", "name": "4-[2-甲基苯基亚甲基氨基]苯甲酸", "code": "CC-Tya24Pk4r", "svg": "side-structures/struct_08.svg"}, {"cat": "N-取代地平（二甲酯）", "name": "N-苯基-2,6-二甲基-4-苯基-1,4-二氢吡啶-3,5-二甲酸二甲酯", "code": "CC-Fz4yQydP5", "svg": "side-structures/struct_09.svg"}, {"cat": "胺取代吖啶二酮", "name": "N-(2-甲氧基乙基)-9-(3-氯-4-羟基-5-甲氧基苯基)-1,2,3,4,5,6,7,8-八氢吖啶-1,8-二酮", "code": "CC-grXwt6umF", "svg": "side-structures/struct_10.svg"}, {"cat": "酰腙", "name": "1-萘甲酰肼缩1,2,3-三唑-5-甲醛酰腙", "code": "CC-wFejWzfY2", "svg": "side-structures/struct_11.svg"}, {"cat": "腙", "name": "N'-(4-氯苯基)-N-(4-溴-2-羟基苯基亚甲基)肼", "code": "CC-SLlJFF1fH", "svg": "side-structures/struct_12.svg"}, {"cat": "酞嗪二酮", "name": "2-((4-氟苯基)甲基)-2,3-二氢酞嗪-1,4-二酮", "code": "CC-Oz39SQRV9", "svg": "side-structures/struct_13.svg"}, {"cat": "2-取代-异吲哚-1,3-二酮", "name": "N-(3-氟苯基)-1H-异吲哚-1,3(2H)-二酮", "code": "CC-y95SBbppd", "svg": "side-structures/struct_14.svg"}];

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
