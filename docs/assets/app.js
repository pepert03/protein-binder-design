import {DATA} from './structures.js';
import {mountViewer} from './renderer.js';
const $=s=>document.querySelector(s);
const select=$('#challenge'),cap=$('#caption'),pp=$('#playpause');
let isEgfr=false;
function onDesign(d,i,n){
  if(!isEgfr||!d)return;
  cap.textContent='EGFR · design '+(i+1)+' of '+n+' · '+d.name+' · Boltz-2 ipTM '+d.ipTM.toFixed(2);
}
const api=mountViewer($('canvas'),DATA,()=>$('#error').hidden=false,onDesign);
if(new URLSearchParams(location.search).get('challenge')==='egfr')select.value='egfr';
function update(){
  isEgfr=select.value==='egfr';
  api.select?.(select.value);
  if(isEgfr){
    $('#caveat').textContent='Predicted complexes of the 20 draft sequences, ranked by ipTM. Binders drawn as a smooth backbone tube along the CA trace. Not experimentally validated.';
    $('#source').href='https://github.com/pepert03/protein-binder-design';
  }else{
    cap.textContent='TNF α · 1TNF · Target structure only';
    $('#caveat').textContent='Candidate complex coordinates for this challenge were not recovered. No binder shown.';
    $('#source').href='https://www.rcsb.org/structure/1TNF';
  }
}
function label(){const on=api.isPlaying();pp.textContent=on?'Pause':'Play';pp.setAttribute('aria-pressed',String(!on));}
select.addEventListener('change',update);
$('#reset').addEventListener('click',()=>api.reset?.());
pp.addEventListener('click',()=>{api.toggle();label();});
update();label();
