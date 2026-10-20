
function decodeMesh(gl, p) {
    const raw = atob(p.b64);
    const bytes = new Uint8Array(raw.length);
    for (let i = 0; i < raw.length; i++) bytes[i] = raw.charCodeAt(i);
    let off = 0;
    const q = new Int16Array(bytes.buffer, off, p.nv * 3); off += p.nv * 6;
    const nn = new Int8Array(bytes.buffer, off, p.nv * 3); off += p.nv * 3;
    let cc = null;
    if (p.has_col) { cc = new Uint8Array(bytes.buffer, off, p.nv * 3); off += p.nv * 3; }
    off = Math.ceil(off / 4) * 4; const idx32 = new Uint32Array(bytes.buffer, off, p.nf * 3);
    const pos = new Float32Array(p.nv * 3);
    const scale = p.span / 32000;
    let maxR = 0; let cx = 0, cy = 0, cz = 0;
    for (let i = 0; i < p.nv; i++) {
        const x = q[i * 3] * scale + p.lo[0], y = q[i * 3 + 1] * scale + p.lo[1], z = q[i * 3 + 2] * scale + p.lo[2];
        pos[i * 3] = x; pos[i * 3 + 1] = y; pos[i * 3 + 2] = z;
        cx += x; cy += y; cz += z;
    }
    cx /= p.nv; cy /= p.nv; cz /= p.nv;
    for (let i = 0; i < p.nv; i++) {
        const dx = pos[i * 3], dy = pos[i * 3 + 1], dz = pos[i * 3 + 2];
        const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
        if (r > maxR) maxR = r;
    }
    const nrm = new Float32Array(p.nv * 3);
    for (let i = 0; i < p.nv * 3; i++) nrm[i] = nn[i] / 120;
    const colf = cc ? new Float32Array(p.nv * 3) : null;
    if (cc && colf) for (let i = 0; i < p.nv * 3; i++) colf[i] = cc[i] / 255;
    const idx16 = new Uint16Array(p.nf * 3);
    for (let i = 0; i < p.nf * 3; i++) idx16[i] = idx32[i];
    const mk = (arr, type) => {
        const b = gl.createBuffer(); if (!b) throw new Error('buffer');
        gl.bindBuffer(type, b); gl.bufferData(type, arr, gl.STATIC_DRAW); return b;
    };
    return { pos: mk(pos, gl.ARRAY_BUFFER), norm: mk(nrm, gl.ARRAY_BUFFER), col: colf ? mk(colf, gl.ARRAY_BUFFER) : null, idx: mk(idx16, gl.ELEMENT_ARRAY_BUFFER), count: p.nf * 3, radius: maxR, centroid: [cx, cy, cz] };
}

// --- minimal mat4 (column-major) ---
function mPerspective(fov, aspect, near, far) {
    const f = 1 / Math.tan(fov / 2), m = new Float32Array(16);
    m[0] = f / aspect; m[5] = f; m[10] = (far + near) / (near - far); m[11] = -1; m[14] = (2 * far * near) / (near - far);
    return m;
}
function mLookAt(eye, c, up) {
    let zx = eye[0] - c[0], zy = eye[1] - c[1], zz = eye[2] - c[2];
    let l = Math.hypot(zx, zy, zz); zx /= l; zy /= l; zz /= l;
    let xx = up[1] * zz - up[2] * zy, xy = up[2] * zx - up[0] * zz, xz = up[0] * zy - up[1] * zx;
    l = Math.hypot(xx, xy, xz); xx /= l; xy /= l; xz /= l;
    const yx = zy * xz - zz * xy, yy = zz * xx - zx * xz, yz = zx * xy - zy * xx;
    const m = new Float32Array(16);
    m[0] = xx; m[1] = yx; m[2] = zx; m[4] = xy; m[5] = yy; m[6] = zy; m[8] = xz; m[9] = yz; m[10] = zz; m[15] = 1;
    m[12] = -(xx * eye[0] + xy * eye[1] + xz * eye[2]);
    m[13] = -(yx * eye[0] + yy * eye[1] + yz * eye[2]);
    m[14] = -(zx * eye[0] + zy * eye[1] + zz * eye[2]);
    return m;
}
function mRotAxis(axis, ang) {
    const [ax, ay, az] = axis, c = Math.cos(ang), s = Math.sin(ang), t = 1 - c;
    return new Float32Array([
        t * ax * ax + c, t * ax * ay + s * az, t * ax * az - s * ay, 0,
        t * ax * ay - s * az, t * ay * ay + c, t * ay * az + s * ax, 0,
        t * ax * az + s * ay, t * ay * az - s * ax, t * az * az + c, 0,
        0, 0, 0, 1]);
}
function mMul(a, b) {
    const o = new Float32Array(16);
    for (let c = 0; c < 4; c++) for (let r = 0; r < 4; r++) {
        o[c * 4 + r] = a[r] * b[c * 4] + a[4 + r] * b[c * 4 + 1] + a[8 + r] * b[c * 4 + 2] + a[12 + r] * b[c * 4 + 3];
    }
    return o;
}
function mTranslate(v) {
    const m = new Float32Array(16); m[0] = m[5] = m[10] = m[15] = 1; m[12] = v[0]; m[13] = v[1]; m[14] = v[2]; return m;
}

const VS = `
attribute vec3 aPos; attribute vec3 aNorm; attribute vec3 aCol;
uniform mat4 uProj, uView, uModel;
varying vec3 vN; varying vec3 vW; varying vec3 vC;
void main(){
  vec4 w = uModel * vec4(aPos, 1.0);
  vW = w.xyz; vN = mat3(uModel[0].xyz, uModel[1].xyz, uModel[2].xyz) * aNorm; vC = aCol;
  gl_Position = uProj * uView * w;
}`;
const FS = `
precision mediump float;
varying vec3 vN; varying vec3 vW; varying vec3 vC;
uniform vec3 uCam, uFog, uBase, uLight;
uniform float uUseCol, uAlpha, uFogN, uFogF, uRim;
void main(){
  vec3 N = normalize(vN); vec3 V = normalize(uCam - vW);
  vec3 L = normalize(uCam + vec3(15.0, 30.0, 15.0));
  if (dot(N, V) < 0.0) N = -N;
  float dif = max(dot(N, L), 0.0);
  vec3 H = normalize(L + V);
  float spec = pow(max(dot(N, H), 0.0), 60.0) * 0.16;
  vec3 base = mix(uBase, vC, uUseCol);
  vec3 col = base * (0.58 + 0.42 * dif) + vec3(spec);
  float fres = pow(1.0 - max(dot(N, V), 0.0), 3.0);
  col += fres * uRim * vec3(0.75, 0.85, 1.0);
  vec3 L2 = normalize(vec3(-0.5, -0.35, -0.6));
  col += base * max(dot(N, L2), 0.0) * 0.18;
  float d = length(uCam - vW);
  float f = smoothstep(uFogN, uFogF, d);
  col = mix(col, uFog, f);
  gl_FragColor = vec4(col, uAlpha);
}`;

function makeProgram(gl) {
    const mk = (type, src) => {
        const s = gl.createShader(type); if (!s) throw new Error('shader');
        gl.shaderSource(s, src); gl.compileShader(s);
        if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s) || 'compile');
        return s;
    };
    const p = gl.createProgram(); if (!p) throw new Error('program');
    gl.attachShader(p, mk(gl.VERTEX_SHADER, VS));
    gl.attachShader(p, mk(gl.FRAGMENT_SHADER, FS));
    gl.bindAttribLocation(p, 0, 'aPos'); gl.bindAttribLocation(p, 1, 'aNorm'); gl.bindAttribLocation(p, 2, 'aCol');
    gl.linkProgram(p);
    if (!gl.getProgramParameter(p, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(p) || 'link');
    return p;
}

function hexToRgb(h) {
    const m = h.replace('#', '');
    if (m.length < 6) return [0.956, 0.949, 0.925];
    return [parseInt(m.slice(0, 2), 16) / 255, parseInt(m.slice(2, 4), 16) / 255, parseInt(m.slice(4, 6), 16) / 255];
}
function ease(t) { return t < 0.5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2; }


export function mountViewer(canvas, data, onError) {
 const gl = canvas.getContext('webgl',{antialias:true,alpha:false,preserveDrawingBuffer:true});
 if(!gl){onError();return ()=>{};}
 const prog=makeProgram(gl);gl.useProgram(prog);
 const u={};for(const n of ['Proj','View','Model','Cam','Fog','Base','Light','UseCol','Alpha','FogN','FogF','Rim'])u[n]=gl.getUniformLocation(prog,'u'+n);
 const loaded={};for(const key of Object.keys(data)){loaded[key]={rec:decodeMesh(gl,data[key].receptor),binder:data[key].binder?decodeMesh(gl,data[key].binder):null};}
 let selected='tnf',yaw=.62,pitch=.12,zoom=3.3,drag=false,lx=0,ly=0,raf,stopped=false;
 gl.enable(gl.DEPTH_TEST);gl.disable(gl.CULL_FACE);gl.clearColor(1,1,1,1);
 function draw(mesh,base){gl.uniformMatrix4fv(u.Model,false,mTranslate([0,0,0]));gl.uniform3fv(u.Base,base);gl.uniform1f(u.UseCol,0);gl.uniform1f(u.Alpha,1);gl.uniform1f(u.Rim,.1);gl.bindBuffer(gl.ARRAY_BUFFER,mesh.pos);gl.vertexAttribPointer(0,3,gl.FLOAT,false,0,0);gl.enableVertexAttribArray(0);gl.bindBuffer(gl.ARRAY_BUFFER,mesh.norm);gl.vertexAttribPointer(1,3,gl.FLOAT,false,0,0);gl.enableVertexAttribArray(1);gl.disableVertexAttribArray(2);gl.vertexAttrib3f(2,...base);gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER,mesh.idx);gl.drawElements(gl.TRIANGLES,mesh.count,gl.UNSIGNED_SHORT,0);}
 function frame(){if(stopped)return;const w=canvas.clientWidth,h=canvas.clientHeight,dpr=Math.min(devicePixelRatio||1,2);if(canvas.width!==Math.round(w*dpr)||canvas.height!==Math.round(h*dpr)){canvas.width=Math.round(w*dpr);canvas.height=Math.round(h*dpr);}gl.viewport(0,0,canvas.width,canvas.height);const {rec,binder}=loaded[selected];const radius=rec.radius*zoom*Math.max(1,.95/(w/h));const eye=[radius*Math.cos(pitch)*Math.sin(yaw),radius*Math.sin(pitch),radius*Math.cos(pitch)*Math.cos(yaw)];gl.uniformMatrix4fv(u.Proj,false,mPerspective(.62,w/h,1,radius*8));gl.uniformMatrix4fv(u.View,false,mLookAt(eye,[0,0,0],[0,1,0]));gl.uniform3fv(u.Cam,eye);gl.uniform3f(u.Fog,1,1,1);gl.uniform1f(u.FogN,radius*2);gl.uniform1f(u.FogF,radius*5);gl.uniform3f(u.Light,.5,.8,.7);gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER_BIT);draw(rec,[.72,.75,.76]);if(binder)draw(binder,[.27,.40,.95]);raf=requestAnimationFrame(frame);}
 const down=e=>{drag=true;lx=e.clientX;ly=e.clientY;try{canvas.setPointerCapture(e.pointerId)}catch(_){}};
 const move=e=>{if(!drag)return;yaw-=(e.clientX-lx)*.008;pitch=Math.max(-1.45,Math.min(1.45,pitch+(e.clientY-ly)*.008));lx=e.clientX;ly=e.clientY;};const up=()=>{drag=false;};const wheel=e=>{e.preventDefault();zoom=Math.max(2.2,Math.min(5.5,zoom*(1+e.deltaY*.001)));};
 canvas.addEventListener('pointerdown',down);canvas.addEventListener('pointermove',move);canvas.addEventListener('pointerup',up);canvas.addEventListener('pointercancel',up);canvas.addEventListener('wheel',wheel,{passive:false});
 const api={select(key){selected=key;zoom=key==='tnf'?3.3:3.1;yaw=key==='tnf'?.62:-2.11;pitch=key==='tnf'?.12:-.43;},reset(){api.select(selected);},dispose(){stopped=true;cancelAnimationFrame(raf);canvas.removeEventListener('pointerdown',down);canvas.removeEventListener('pointermove',move);canvas.removeEventListener('pointerup',up);canvas.removeEventListener('pointercancel',up);canvas.removeEventListener('wheel',wheel);}};
 frame();return api;
}
