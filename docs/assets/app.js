import {DATA} from './structures.js';
import {mountViewer} from './renderer.js';
const api=mountViewer(document.querySelector('canvas'),DATA,()=>document.querySelector('#error').hidden=false);
const select=document.querySelector('#challenge');if(new URLSearchParams(location.search).get('challenge')==='egfr')select.value='egfr';
function update(){api.select?.(select.value);const egfr=select.value==='egfr';document.querySelector('#caption').textContent=egfr?'EGFR + r02_l07 · In silico model':'TNF α · 1TNF · Target structure only';document.querySelector('#caveat').textContent=egfr?'Computational prediction; ribbon cartoon derived from the CA trace. Not experimentally validated.':'Candidate coordinates unavailable. No binder shown.';document.querySelector('#source').href=egfr?'https://github.com/pepert03/protein-binder-design':'https://www.rcsb.org/structure/1TNF';}
select.addEventListener('change',update);document.querySelector('#reset').addEventListener('click',()=>api.reset?.());update();
