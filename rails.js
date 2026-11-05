// 首页两侧滚动结构推送（产物库真实结构）。由 search/_gen_rails.js 生成（与检索站同一个渲染器）。
(function(){
  var RAILS = [{"cat":"席夫碱","name":"N-(5-氯-2,4-二羟基苯基)亚甲基-3-甲基苯胺","code":"CC-2I8zxxJmA","svg":"side-structures/struct_01.svg"},{"cat":"N-取代地平（二甲酯）","name":"N-(3-(三氟甲基)苯基)-2,6-二甲基-4-(5-乙氧基呋喃-2-基)-1,4-二氢吡啶-3,5-二甲酸二甲酯","code":"CC-ELLk2CRoi","svg":"side-structures/struct_02.svg"},{"cat":"胺取代吖啶二酮","name":"N-(3,5-二溴苯基)-9-(4-乙酰苯基)-1,2,3,4,5,6,7,8-八氢吖啶-1,8-二酮","code":"CC-bEuBVoQFS","svg":"side-structures/struct_03.svg"},{"cat":"酰腙","name":"水杨酰肼缩4-溴噻吩-2-甲醛酰腙","code":"CC-mRx17fOXH","svg":"side-structures/struct_04.svg"},{"cat":"腙","name":"N'-([trifluorometh]肼（待核实）)-N-(3-羧基苯甲亚甲基)肼","code":"CC-HR8BrrnY5","svg":"side-structures/struct_05.svg"},{"cat":"酞嗪二酮","name":"5,6,7,8-四氯-2-(4-乙基苯基)-2,3-二氢酞嗪-1,4-二酮","code":"CC-GF1TLL42I","svg":"side-structures/struct_06.svg"},{"cat":"2-取代-异吲哚-1,3-二酮","name":"2-(2-氯苯甲酰基)-2,3-二氢吡啶并[2,3-d]哒嗪-1,4-二酮","code":"CC-HuU2g06cT","svg":"side-structures/struct_07.svg"},{"cat":"DPPH醛肟酯","name":"2-氟苯甲醛（O-二苯基膦酰醛肟）","code":"CC-DPPH0009","svg":"side-structures/struct_08.svg"},{"cat":"DOPO-席夫碱","name":"DOPO-P-席夫碱（苯甲醛/甲胺）","code":"CC-DOPO000048","svg":"side-structures/struct_09.svg"},{"cat":"DOPO-酰腙","name":"DOPO-P-酰腙（苯甲醛/乙酰肼）","code":"CC-DOPO005032","svg":"side-structures/struct_10.svg"},{"cat":"DOPO-腙","name":"DOPO-P-腙（甲醛/苯基肼）","code":"CC-DOPO010031","svg":"side-structures/struct_11.svg"},{"cat":"DOPO-缩氨脲","name":"DOPO-P-缩氨脲（苯甲醛/N-甲基缩氨脲）","code":"CC-DOPO014447","svg":"side-structures/struct_12.svg"},{"cat":"氨基苯胂酸席夫碱","name":"4-((3,4-二甲氧基萘基-1-基)亚甲基氨基)苯胂酸","code":"CC-K5vges6do","svg":"side-structures/struct_13.svg"},{"cat":"氨基苯胂酸席夫碱","name":"4-(6-溴萘基-2-基亚甲基氨基)苯胂酸","code":"CC-eO0xzj59y","svg":"side-structures/struct_14.svg"}];

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
