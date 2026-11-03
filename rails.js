// 首页两侧滚动结构推送（产物库真实结构）。由 gen_side_rails.py 生成。
(function(){
  var RAILS = [{"cat": "席夫碱", "name": "N-(2-甲基噻吩-3-基亚甲基)-2,3-二氟苯胺", "code": "CC-c9jTlJfmJ", "svg": "side-structures/struct_01.svg"}, {"cat": "N-取代地平（二甲酯）", "name": "N-(3,4-二氟苯基)-2,6-二甲基-4-(呋喃-3-基)-1,4-二氢吡啶-3,5-二甲酸二甲酯", "code": "CC-ivPBhwhUK", "svg": "side-structures/struct_02.svg"}, {"cat": "胺取代吖啶二酮", "name": "N-(2-羟基苯基)-9-(4-(4-硝基苯基)苯基)-1,2,3,4,5,6,7,8-八氢吖啶-1,8-二酮", "code": "CC-Yt6ug297O", "svg": "side-structures/struct_03.svg"}, {"cat": "酰腙", "name": "2-甲氧基苯甲酰肼缩2,3,5,6-四氟苯甲醛酰腙", "code": "CC-VMflURel9", "svg": "side-structures/struct_04.svg"}, {"cat": "腙", "name": "N'-(—)-N-(对羟基苯甲亚甲基)肼", "code": "CC-8saQvQ9CS", "svg": "side-structures/struct_05.svg"}, {"cat": "酞嗪二酮", "name": "5-甲基-2-(2-(4-氟苯氧基)乙酰基肼基)-2,3-二氢酞嗪-1,4-二酮", "code": "CC-MMTnARBt7", "svg": "side-structures/struct_06.svg"}, {"cat": "2-取代-异吲哚-1,3-二酮", "name": "2-(2-羟基苯甲酰基)-2,3-二氢酞嗪-1,4-二酮", "code": "CC-wLyaiiJKa", "svg": "side-structures/struct_07.svg"}, {"cat": "席夫碱", "name": "4-[(4-(4-硝基苯基)苯基)亚甲基氨基]吡啶", "code": "CC-h1RnoqFsu", "svg": "side-structures/struct_08.svg"}, {"cat": "N-取代地平（二甲酯）", "name": "N-(3-氰基苯基)-2,6-二甲基-4-(3-(4-氯苯基)苯基)-1,4-二氢吡啶-3,5-二甲酸二甲酯", "code": "CC-vAd2baldd", "svg": "side-structures/struct_09.svg"}, {"cat": "胺取代吖啶二酮", "name": "N-正丁基-9-(1-苯并噻吩-4-基)-1,2,3,4,5,6,7,8-八氢吖啶-1,8-二酮", "code": "CC-ULtItIHqR", "svg": "side-structures/struct_10.svg"}, {"cat": "酰腙", "name": "异烟酰肼缩5-氮杂吲哚-3-甲醛酰腙", "code": "CC-cqOO4a4td", "svg": "side-structures/struct_11.svg"}, {"cat": "腙", "name": "N'-(3,4-二甲基苯基)-N-((2,4-二羟基苯基)亚甲基)肼", "code": "CC-ZcDSosTeF", "svg": "side-structures/struct_12.svg"}, {"cat": "酞嗪二酮", "name": "6-氯-2-(4-溴-2-氟苯基)-2,3-二氢酞嗪-1,4-二酮", "code": "CC-7kNhz74Sg", "svg": "side-structures/struct_13.svg"}, {"cat": "2-取代-异吲哚-1,3-二酮", "name": "N-(3-甲氧基苯基)-1H-异吲哚-1,3(2H)-二酮", "code": "CC-5bOl0AFum", "svg": "side-structures/struct_14.svg"}];

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
