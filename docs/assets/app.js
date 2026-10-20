import {DATA} from './structures.js';
import {mountViewer} from './renderer.js';
const $=s=>document.querySelector(s);
const select=$('#challenge'),cap=$('#caption'),pp=$('#playpause');
function onDesign(d,i,n){
  if(!d)return;
  if(select.value==='egfr')cap.textContent='EGFR · design '+(i+1)+' of '+n+' · '+d.name+' · Boltz-2 ipTM '+d.ipTM.toFixed(2);
  else if(d.name==='tnf_01')cap.textContent='TNF α · tnf_01 · submitted candidate · predicted crop complex (AF2-Multimer)';
  else cap.textContent='TNF α · '+d.name+' · '+d.tag+' · failed screening gate';
}
const api=mountViewer($('canvas'),DATA,()=>$('#error').hidden=false,onDesign);
function update(){
  api.select?.(select.value);
  if(select.value==='egfr'){
    $('#caveat').textContent='Predicted complexes of the 20 draft sequences, ranked by ipTM. Binders drawn as ribbon cartoons from the predicted backbone. Not experimentally validated.';
    $('#source').href='https://github.com/pepert03/protein-binder-design';
  }else{
    $('#caveat').textContent='Predicted AF2-Multimer crop complexes, rigidly aligned onto the 1TNF trimer for context (not intact-trimer predictions). tnf_01 is the submitted candidate; its intact-trimer re-prediction had a severe steric clash (minimum heavy-atom distance 0.27 A) and is not shown. The 24 session-3 redesigns did not pass the screening gate (best norm-iPAE 0.597). Submitted ranks 2-20 (65-aa family) showed no real interface in the intact audit and are not shown. Not experimentally validated.';
    $('#source').href='https://www.rcsb.org/structure/1TNF';
  }
}
function label(){const on=api.isPlaying();pp.textContent=on?'Pause':'Play';pp.setAttribute('aria-pressed',String(!on));}
select.addEventListener('change',update);
$('#reset').addEventListener('click',()=>api.reset?.());
pp.addEventListener('click',()=>{api.toggle();label();});
update();label();
