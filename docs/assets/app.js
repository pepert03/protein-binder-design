import {DATA} from './structures.js';
import {mountViewer} from './renderer.js';
const $=s=>document.querySelector(s);
const cap=$('#caption'),pp=$('#playpause');
const dd=$('#challenge-dd'),btn=$('#challenge-btn'),menu=$('#challenge-menu'),lab=$('#challenge-label');
let current='egfr';
function setChallenge(v){current=v;const opt=menu.querySelector(`[data-value="${v}"]`);lab.textContent=opt.textContent;menu.querySelectorAll('[role=option]').forEach(o=>o.setAttribute('aria-selected',String(o===opt)));update();}
function openMenu(){menu.hidden=false;btn.setAttribute('aria-expanded','true');(menu.querySelector('[aria-selected=true]')||menu.firstElementChild).focus?.();}
function closeMenu(){menu.hidden=true;btn.setAttribute('aria-expanded','false');}
btn.addEventListener('click',()=>{menu.hidden?openMenu():closeMenu();});
menu.addEventListener('click',e=>{const o=e.target.closest('[role=option]');if(o){setChallenge(o.dataset.value);closeMenu();btn.focus();}});
dd.addEventListener('keydown',e=>{const opts=[...menu.querySelectorAll('[role=option]')];const i=opts.indexOf(document.activeElement);if(e.key==='Escape'){closeMenu();btn.focus();}else if(e.key==='ArrowDown'){e.preventDefault();if(menu.hidden)openMenu();else (opts[i+1]||opts[0]).focus();}else if(e.key==='ArrowUp'){e.preventDefault();if(!menu.hidden)(opts[i-1]||opts[opts.length-1]).focus();}else if((e.key==='Enter'||e.key===' ')&&document.activeElement!==btn&&i>=0){e.preventDefault();setChallenge(opts[i].dataset.value);closeMenu();btn.focus();}});
document.addEventListener('pointerdown',e=>{if(!dd.contains(e.target))closeMenu();});
function onDesign(d,i,n){
  if(!d)return;
  if(current==='egfr')cap.textContent='EGFR · design '+(i+1)+' of '+n+' · '+d.name+' · Boltz-2 ipTM '+d.ipTM.toFixed(2);
  else if(d.name==='tnf_01')cap.textContent='TNF α · tnf_01 · submitted candidate · predicted crop complex (AF2-Multimer)';
  else cap.textContent='TNF α · '+d.name+' · '+d.tag+' · failed screening gate';
}
const api=mountViewer($('canvas'),DATA,()=>$('#error').hidden=false,onDesign);
function update(){
  api.select?.(current);
  if(current==='egfr'){
    $('#caveat').textContent='Predicted complexes of the 20 draft sequences, ranked by ipTM. Binders drawn as ribbon cartoons from the predicted backbone. Not experimentally validated.';
    $('#source').href='https://github.com/pepert03/protein-binder-design';
  }else{
    $('#caveat').textContent='Predicted AF2-Multimer crop complexes, rigidly aligned onto the 1TNF trimer for context (not intact-trimer predictions). tnf_01 is the submitted candidate; its intact-trimer re-prediction had a severe steric clash (minimum heavy-atom distance 0.27 A) and is not shown. The 24 session-3 redesigns did not pass the screening gate (best norm-iPAE 0.597). Submitted ranks 2-20 (65-aa family) showed no real interface in the intact audit and are not shown. Not experimentally validated.';
    $('#source').href='https://www.rcsb.org/structure/1TNF';
  }
}
function label(){const on=api.isPlaying();pp.textContent=on?'Pause':'Play';pp.setAttribute('aria-pressed',String(!on));}
$('#reset').addEventListener('click',()=>api.reset?.());
pp.addEventListener('click',()=>{api.toggle();label();});
update();label();
