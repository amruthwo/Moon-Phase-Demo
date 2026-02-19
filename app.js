(() => {
  const SYNODIC_DAYS = 29.53;
  const phaseNames = ["New Moon","Waxing Crescent","First Quarter","Waxing Gibbous","Full Moon","Waning Gibbous","Last Quarter","Waning Crescent"];
  const phases = phaseNames.map((name,i)=>({name, t:i/8, day:SYNODIC_DAYS*(i/8)}));

  const elMain = document.getElementById('main');
  const elInset = document.getElementById('inset');
  const elSlider = document.getElementById('slider');
  const elTicks = document.getElementById('ticks');
  const elDayTxt = document.getElementById('dayTxt');
  const elSynTxt = document.getElementById('synTxt');
  const elPhaseTxt = document.getElementById('phaseTxt');
  const elPhaseChip = document.getElementById('phaseChip');
  const elPlayBtn = document.getElementById('playBtn');
  const elCamBtn = document.getElementById('camBtn');
  const elCamTxt = document.getElementById('camTxt');
  const elBodyBtn = document.getElementById('bodyBtn');
  const elBodyTxt = document.getElementById('bodyTxt');
  const elPlayTxt = document.getElementById('playTxt');
  const elErr = document.getElementById('err');

  elSynTxt.textContent = SYNODIC_DAYS.toFixed(2);
  elSlider.max = String(SYNODIC_DAYS);

  function showErr(msg){ elErr.hidden = false; elErr.textContent = msg; }
  window.addEventListener('error', (e)=> showErr("JS error: " + (e.message || e.error || "Unknown")));
  window.addEventListener('unhandledrejection', (e)=> showErr("Promise rejection: " + String((e.reason && (e.reason.stack||e.reason.message)) || e.reason)));

  function buildTicks(){
    elTicks.innerHTML = "";
    for(const p of phases){
      const x = (p.day / SYNODIC_DAYS) * 100;
      const t = document.createElement('div');
      t.className = 'tick';
      t.style.left = x + '%';
      const lbl = document.createElement('div');
      lbl.className = 'tickLbl';
      lbl.style.left = x + '%';
      lbl.innerHTML = `<div class="nm">${p.name}</div><div class="dy">${p.day.toFixed(1)}d</div>`;
      elTicks.appendChild(t);
      elTicks.appendChild(lbl);
    }
  }
  buildTicks();

  // --- Minimal math (column-major matrices) ---
  const V3 = (x=0,y=0,z=0)=>({x,y,z});
  const v3_add=(a,b)=>V3(a.x+b.x,a.y+b.y,a.z+b.z);
  const v3_sub=(a,b)=>V3(a.x-b.x,a.y-b.y,a.z-b.z);
  const v3_mul=(a,s)=>V3(a.x*s,a.y*s,a.z*s);
  const v3_dot=(a,b)=>a.x*b.x+a.y*b.y+a.z*b.z;
  const v3_len=(a)=>Math.hypot(a.x,a.y,a.z);
  const v3_norm=(a)=>{ const l=v3_len(a)||1; return V3(a.x/l,a.y/l,a.z/l); };
  const v3_cross=(a,b)=>V3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);

  const M4 = ()=>new Float32Array(16);
  const m4_ident=(m)=>{ m.set([1,0,0,0,  0,1,0,0,  0,0,1,0,  0,0,0,1]); return m; };
  const m4_mul=(a,b,out)=>{
    const o=out||M4();
    const a00=a[0],a01=a[1],a02=a[2],a03=a[3];
    const a10=a[4],a11=a[5],a12=a[6],a13=a[7];
    const a20=a[8],a21=a[9],a22=a[10],a23=a[11];
    const a30=a[12],a31=a[13],a32=a[14],a33=a[15];

    const b00=b[0],b01=b[1],b02=b[2],b03=b[3];
    const b10=b[4],b11=b[5],b12=b[6],b13=b[7];
    const b20=b[8],b21=b[9],b22=b[10],b23=b[11];
    const b30=b[12],b31=b[13],b32=b[14],b33=b[15];

    o[0]=a00*b00 + a10*b01 + a20*b02 + a30*b03;
    o[1]=a01*b00 + a11*b01 + a21*b02 + a31*b03;
    o[2]=a02*b00 + a12*b01 + a22*b02 + a32*b03;
    o[3]=a03*b00 + a13*b01 + a23*b02 + a33*b03;

    o[4]=a00*b10 + a10*b11 + a20*b12 + a30*b13;
    o[5]=a01*b10 + a11*b11 + a21*b12 + a31*b13;
    o[6]=a02*b10 + a12*b11 + a22*b12 + a32*b13;
    o[7]=a03*b10 + a13*b11 + a23*b12 + a33*b13;

    o[8]=a00*b20 + a10*b21 + a20*b22 + a30*b23;
    o[9]=a01*b20 + a11*b21 + a21*b22 + a31*b23;
    o[10]=a02*b20 + a12*b21 + a22*b22 + a32*b23;
    o[11]=a03*b20 + a13*b21 + a23*b22 + a33*b23;

    o[12]=a00*b30 + a10*b31 + a20*b32 + a30*b33;
    o[13]=a01*b30 + a11*b31 + a21*b32 + a31*b33;
    o[14]=a02*b30 + a12*b31 + a22*b32 + a32*b33;
    o[15]=a03*b30 + a13*b31 + a23*b32 + a33*b33;
    return o;
  };
  const m4_translate=(m,v)=>{ const t=M4(); m4_ident(t); t[12]=v.x; t[13]=v.y; t[14]=v.z; return m4_mul(m,t,m); };
  const m4_scale=(m,v)=>{ const s=M4(); m4_ident(s); s[0]=v.x; s[5]=v.y; s[10]=v.z; return m4_mul(m,s,m); };
  const m4_rotY=(m,rad)=>{ const r=M4(); m4_ident(r); const c=Math.cos(rad), s=Math.sin(rad); r[0]=c; r[2]=-s; r[8]=s; r[10]=c; return m4_mul(m,r,m); };
  const m4_lookAt=(eye,at,up)=>{
    const f = v3_norm(v3_sub(at, eye));
    const s = v3_norm(v3_cross(f, up));
    const u = v3_cross(s, f);
    const m=M4(); m4_ident(m);
    m[0]=s.x; m[1]=u.x; m[2]=-f.x; m[3]=0;
    m[4]=s.y; m[5]=u.y; m[6]=-f.y; m[7]=0;
    m[8]=s.z; m[9]=u.z; m[10]=-f.z; m[11]=0;
    m[12]=-v3_dot(s,eye);
    m[13]=-v3_dot(u,eye);
    m[14]= v3_dot(f,eye);
    m[15]=1;
    return m;
  };
  const m4_persp=(fovy,aspect,near,far)=>{
    const f=1/Math.tan(fovy/2);
    const nf=1/(near-far);
    const m=M4();
    m[0]=f/aspect; m[1]=0; m[2]=0; m[3]=0;
    m[4]=0; m[5]=f; m[6]=0; m[7]=0;
    m[8]=0; m[9]=0; m[10]=(far+near)*nf; m[11]=-1;
    m[12]=0; m[13]=0; m[14]=(2*far*near)*nf; m[15]=0;
    return m;
  };

  // --- WebGL helpers ---
  function glInit(canvas){
    const gl = canvas.getContext('webgl', {antialias:true, alpha:true, premultipliedAlpha:false});
    if(!gl) throw new Error('WebGL not available');
    return gl;
  }
  function compile(gl, type, src){
    const sh=gl.createShader(type);
    gl.shaderSource(sh, src);
    gl.compileShader(sh);
    if(!gl.getShaderParameter(sh, gl.COMPILE_STATUS)){
      throw new Error(gl.getShaderInfoLog(sh) || 'Shader compile failed');
    }
    return sh;
  }
  function link(gl, vsSrc, fsSrc){
    const p=gl.createProgram();
    gl.attachShader(p, compile(gl, gl.VERTEX_SHADER, vsSrc));
    gl.attachShader(p, compile(gl, gl.FRAGMENT_SHADER, fsSrc));
    gl.linkProgram(p);
    if(!gl.getProgramParameter(p, gl.LINK_STATUS)){
      throw new Error(gl.getProgramInfoLog(p) || 'Program link failed');
    }
    return p;
  }
  function makeSphere(lat=48, lon=72){
    const pos=[], nrm=[], uv=[], idx=[];
    for(let y=0;y<=lat;y++){
      const v=y/lat;
      const th=v*Math.PI;
      const st=Math.sin(th), ct=Math.cos(th);
      for(let x=0;x<=lon;x++){
        const u=x/lon;
        const ph=u*2*Math.PI;
        const sp=Math.sin(ph), cp=Math.cos(ph);
        const nx=cp*st, ny=ct, nz=sp*st;
        pos.push(nx,ny,nz);
        nrm.push(nx,ny,nz);
        uv.push(u, 1-v);
      }
    }
    for(let y=0;y<lat;y++){
      for(let x=0;x<lon;x++){
        const i0=y*(lon+1)+x;
        const i1=i0+1;
        const i2=i0+(lon+1);
        const i3=i2+1;
        idx.push(i0,i2,i1, i1,i2,i3);
      }
    }
    return {pos:new Float32Array(pos), nrm:new Float32Array(nrm), uv:new Float32Array(uv), idx:new Uint16Array(idx)};
  }

  const VS = `
    attribute vec3 aPos;
    attribute vec3 aNrm;
    attribute vec2 aUv;
    uniform mat4 uMVP;
    uniform mat4 uM;
    varying vec3 vN;
    varying vec3 vW;
    varying vec2 vUv;
    void main(){
      vec4 w = uM * vec4(aPos,1.0);
      vW = w.xyz;
      vN = mat3(uM) * aNrm;
      vUv = aUv;
      gl_Position = uMVP * vec4(aPos,1.0);
    }
  `;

  const FS_EARTH = `
    precision highp float;
    varying vec3 vN;
    varying vec3 vW;
    varying vec2 vUv;
    uniform vec3 uSunPos;
    uniform vec3 uCamPos;
    uniform vec3 uMoonDir;
    uniform float uTexRot;

    float hash(vec2 p){ return fract(sin(dot(p, vec2(127.1,311.7))) * 43758.5453123); }
    float noise(vec2 p){
      vec2 i=floor(p), f=fract(p);
      float a=hash(i);
      float b=hash(i+vec2(1,0));
      float c=hash(i+vec2(0,1));
      float d=hash(i+vec2(1,1));
      vec2 u=f*f*(3.0-2.0*f);
      return mix(a,b,u.x) + (c-a)*u.y*(1.0-u.x) + (d-b)*u.x*u.y;
    }
    float fbm(vec2 p){
      float v=0.0, a=0.5;
      for(int i=0;i<5;i++){
        v += a*noise(p);
        p *= 2.0;
        a *= 0.5;
      }
      return v;
    }
    vec2 rot2(vec2 p, float a){
      float c=cos(a), s=sin(a);
      return vec2(c*p.x - s*p.y, s*p.x + c*p.y);
    }
    float circle(vec2 p, vec2 c, float r){
      float d = length(p-c);
      return smoothstep(r, r-0.015, d);
    }

    void main(){
      vec3 N = normalize(vN);
      vec3 L = normalize(uSunPos - vW);
      vec3 V = normalize(uCamPos - vW);
      float wrap = 0.35;
      float diff = clamp((dot(N,L)+wrap)/(1.0+wrap), 0.0, 1.0);

      vec2 uv = vUv;
      uv.x = fract(uv.x + uTexRot/(2.0*3.14159265));

      vec2 p = vec2(uv.x*6.0, uv.y*3.0);
      p = rot2(p, 0.6);
      float n = fbm(p);
      float land = smoothstep(0.52, 0.60, n);

      vec3 oceanDeep = vec3(0.03,0.24,0.50);
      vec3 ocean = vec3(0.06,0.35,0.65);
      vec3 landG = vec3(0.10,0.55,0.26);
      vec3 landB = vec3(0.18,0.38,0.20);

      vec3 base = mix(oceanDeep, ocean, smoothstep(0.0,1.0, uv.y));
      base = mix(base, mix(landB, landG, smoothstep(0.0,1.0, uv.y)), land);

      vec3 night = base*0.15;
      vec3 col = mix(night, base, diff);

      vec3 H = normalize(L+V);
      float spec = pow(max(dot(N,H),0.0), 64.0) * 0.35;
      float fres = pow(1.0 - max(dot(N,V),0.0), 3.0);
      col += spec + fres*0.20;

      // Eyes stable in normal-space basis facing Moon (not tied to texture rotation)
      vec3 f = normalize(uMoonDir);
      vec3 up = vec3(0.0,1.0,0.0);
      // Tangent basis around face direction. Use right = f x up for a consistent handedness.
      vec3 r = normalize(cross(f, up));
      vec3 u = cross(r, f);

      float face = dot(N, f);
      if(face > 0.10){
        vec2 q = vec2(dot(N, r), dot(N, u));
        float sep = 0.25;
        float yoff = 0.18;
        vec2 e1 = vec2(+sep, yoff);
        vec2 e2 = vec2(-sep, yoff);
        float w = smoothstep(0.10, 0.50, face);
        float eyeW = max(circle(q, e1, 0.20), circle(q, e2, 0.20));
        vec2 p1 = e1 + vec2(-0.06, -0.03);
        vec2 p2 = e2 + vec2(+0.06, -0.03);
        float pup = max(circle(q, p1, 0.08), circle(q, p2, 0.08));

        vec3 eyeWhite = vec3(0.98,0.98,0.98);
        vec3 pupCol = vec3(0.05,0.05,0.06);
        col = mix(col, eyeWhite, eyeW*w);
        col = mix(col, pupCol, pup*w);
      }

      gl_FragColor = vec4(col, 1.0);
    }
  `;

  const FS_MOON = `
    precision highp float;
    varying vec3 vN;
    varying vec3 vW;
    uniform vec3 uSunPos;
    uniform vec3 uCamPos;
    void main(){
      // Outward surface normal (do NOT negate). Negating flips day/night
      // and makes the sun-facing side go dark in the main scene.
      vec3 N = normalize(vN);
      vec3 L = normalize(uSunPos - vW);
      vec3 V = normalize(uCamPos - vW);
      float ndl = max(dot(N,L), 0.0);
      float lit = smoothstep(0.02, 0.35, ndl);
      float fres = pow(1.0 - max(dot(N,V),0.0), 2.0);
      vec3 base = vec3(0.70,0.70,0.72);
      vec3 night = vec3(0.08,0.08,0.10);
      vec3 col = mix(night, base, lit);
      col += fres*0.08;
      gl_FragColor = vec4(col,1.0);
    }
  `;

  const FS_SUN = `
    precision highp float;
    varying vec3 vN;
    varying vec3 vW;
    uniform vec3 uCamPos;
    void main(){
      vec3 N = normalize(vN);
      vec3 V = normalize(uCamPos - vW);
      float rim = pow(1.0 - max(dot(N,V),0.0), 2.5);
      vec3 core = vec3(1.00,0.92,0.55);
      vec3 edge = vec3(1.00,0.58,0.18);
      vec3 col = mix(core, edge, rim);
      gl_FragColor = vec4(col, 1.0);
    }
  `;


const FS_SOLID = `precision highp float;
varying vec3 vW;
varying vec3 vN;
uniform vec3 uSunPos;
uniform vec3 uCamPos;
uniform vec3 uBaseColor;
void main(){
  vec3 N = normalize(vN);
  vec3 L = normalize(uSunPos - vW);
  vec3 V = normalize(uCamPos - vW);
  float diff = max(dot(N,L), 0.0);
  vec3 R = reflect(-L, N);
  float spec = pow(max(dot(R,V), 0.0), 32.0);
  vec3 col = uBaseColor * (0.28 + 0.92*diff) + vec3(1.0)*0.18*spec;
  // mild ambient lift
  col = mix(col, vec3(0.06,0.08,0.12), 0.05);
  gl_FragColor = vec4(pow(col, vec3(0.4545)), 1.0);
}`; 

  class Renderer {
    constructor(canvas){
      this.canvas = canvas;
      this.gl = glInit(canvas);
      this.dpr = 1;
      this.w = 1; this.h = 1;
      const gl=this.gl;
      gl.enable(gl.DEPTH_TEST);
      // Disable face culling to avoid inside-out sphere issues across browsers.
      // For closed spheres, back faces will be depth-occluded anyway.
      gl.disable(gl.CULL_FACE);
      gl.clearColor(0,0,0,0);

      this.progEarth = link(gl, VS, FS_EARTH);
      this.progMoon  = link(gl, VS, FS_MOON);
      this.progSun   = link(gl, VS, FS_SUN);

      this.progSolid = link(gl, VS, FS_SOLID);

      this.mesh = makeSphere(48,72);
      this.buffers = { pos:gl.createBuffer(), nrm:gl.createBuffer(), uv:gl.createBuffer(), idx:gl.createBuffer() };

      gl.bindBuffer(gl.ARRAY_BUFFER, this.buffers.pos);
      gl.bufferData(gl.ARRAY_BUFFER, this.mesh.pos, gl.STATIC_DRAW);
      gl.bindBuffer(gl.ARRAY_BUFFER, this.buffers.nrm);
      gl.bufferData(gl.ARRAY_BUFFER, this.mesh.nrm, gl.STATIC_DRAW);
      gl.bindBuffer(gl.ARRAY_BUFFER, this.buffers.uv);
      gl.bufferData(gl.ARRAY_BUFFER, this.mesh.uv, gl.STATIC_DRAW);
      gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.buffers.idx);
      gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, this.mesh.idx, gl.STATIC_DRAW);
    }
    resize(){
      const r = this.canvas.getBoundingClientRect();
      const dpr = Math.max(1, Math.min(2.25, window.devicePixelRatio || 1));
      const w = Math.max(2, Math.floor(r.width * dpr));
      const h = Math.max(2, Math.floor(r.height * dpr));
      if (w !== this.w || h !== this.h || dpr !== this.dpr){
        this.dpr = dpr; this.w = w; this.h = h;
        this.canvas.width = w; this.canvas.height = h;
        this.gl.viewport(0,0,w,h);
      }
    }
    clear(){ this.gl.clear(this.gl.COLOR_BUFFER_BIT | this.gl.DEPTH_BUFFER_BIT); }
    drawSphere(prog, model, view, proj, uniforms){
      const gl=this.gl;
      // Guard against invalid/foreign programs (can happen after context loss or if
      // a program created on another canvas/context is accidentally passed in).
      if(!prog || (gl.isProgram && !gl.isProgram(prog))) return;
      gl.useProgram(prog);

      const aPos = gl.getAttribLocation(prog, 'aPos');
      const aNrm = gl.getAttribLocation(prog, 'aNrm');
      const aUv  = gl.getAttribLocation(prog, 'aUv');

      gl.bindBuffer(gl.ARRAY_BUFFER, this.buffers.pos);
      gl.enableVertexAttribArray(aPos);
      gl.vertexAttribPointer(aPos, 3, gl.FLOAT, false, 0, 0);

      gl.bindBuffer(gl.ARRAY_BUFFER, this.buffers.nrm);
      gl.enableVertexAttribArray(aNrm);
      gl.vertexAttribPointer(aNrm, 3, gl.FLOAT, false, 0, 0);

      gl.bindBuffer(gl.ARRAY_BUFFER, this.buffers.uv);
      gl.enableVertexAttribArray(aUv);
      gl.vertexAttribPointer(aUv, 2, gl.FLOAT, false, 0, 0);

      gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.buffers.idx);

      const mvp = m4_mul(proj, m4_mul(view, model, M4()), M4());

      const uMVP = gl.getUniformLocation(prog, 'uMVP');
      const uM = gl.getUniformLocation(prog, 'uM');
      gl.uniformMatrix4fv(uMVP, false, mvp);
      gl.uniformMatrix4fv(uM, false, model);

      for (const k in uniforms){
        const loc = gl.getUniformLocation(prog, k);
        if (loc === null) continue;
        const v = uniforms[k];
        if (typeof v === 'number') {
          gl.uniform1f(loc, v);
        } else if (v && v.x !== undefined) {
          gl.uniform3f(loc, v.x, v.y, v.z);
        } else if (Array.isArray(v) || (v instanceof Float32Array)) {
          if (v.length === 2) gl.uniform2f(loc, v[0], v[1]);
          else if (v.length === 3) gl.uniform3f(loc, v[0], v[1], v[2]);
          else if (v.length === 4) gl.uniform4f(loc, v[0], v[1], v[2], v[3]);
          else if (v.length === 16) gl.uniformMatrix4fv(loc, false, v);
        }
      }

      gl.drawElements(gl.TRIANGLES, this.mesh.idx.length, gl.UNSIGNED_SHORT, 0);
    }
  }

  const Rmain = new Renderer(elMain);
  const Rinset = new Renderer(elInset);

  const earthR = 1.25, moonR = 0.42, orbitR = 3, sunR = 1.95;

  const earthPos = V3(0,0,0);
  const sunPos   = V3(12.0, 0.0, 0.0);
  
// --- Debug hooks (expose live positions to DevTools) ---

const moonPos  = V3(orbitR, 0.0, 0.0);
  const moonDir  = V3(1,0,0);

  // Main camera orbit/zoom/pan controls
  const camUp = V3(0,1,0);
  const cam = {
    target: V3(0,0,0),
    radius: 12.0,
    theta: Math.PI/2,   // +Z for perfect side-on
    phi: Math.PI/2,
    minR: 3.2,
    maxR: 18.0
  };
  const clamp=(v,a,b)=>Math.max(a,Math.min(b,v));
  const camPosFromSpherical=()=>{
    const sinPhi = Math.sin(cam.phi);
    return V3(
      cam.target.x + cam.radius * sinPhi * Math.cos(cam.theta),
      cam.target.y + cam.radius * Math.cos(cam.phi),
      cam.target.z + cam.radius * sinPhi * Math.sin(cam.theta)
    );
  };


  // Camera presets: true side-on and true top-down-from-+Y
  const CAM_PRESET = { SIDE: "side", TOP: "top" };
  let camPreset = CAM_PRESET.SIDE;

  const setCamPreset = (preset) => {
    camPreset = preset;
    if (preset === CAM_PRESET.SIDE) {
      // Perfect side-on: viewer at +Z looking to origin
      cam.theta = Math.PI/2;
      cam.phi   = Math.PI/2;
      elCamTxt.textContent = "Side";
    } else {
      // True top-down: camera above (+Y) looking down (-Y)
      // Keep a tiny tilt to avoid singularity
      cam.theta = Math.PI/2;
      cam.phi   = 0.18;
      elCamTxt.textContent = "Top";
    }
  };

  const input = {dragging:false, mode:"orbit", px:0, py:0};
  elMain.addEventListener('contextmenu', (e)=>e.preventDefault());
  elMain.addEventListener('pointerdown', (e)=>{
    elMain.setPointerCapture(e.pointerId);
    input.dragging = true;
    input.px = e.clientX; input.py = e.clientY;
    const panMode = (e.button === 2) || e.shiftKey;
    input.mode = panMode ? "pan" : "orbit";
  });
  elMain.addEventListener('pointermove', (e)=>{
    if(!input.dragging) return;
    const dx = e.clientX - input.px;
    const dy = e.clientY - input.py;
    input.px = e.clientX; input.py = e.clientY;

    const r = elMain.getBoundingClientRect();
    const scale = 1 / Math.max(180, Math.min(900, r.width));
    if(input.mode === "orbit"){
      cam.theta += dx * (2.2*scale);
      cam.phi   += dy * (2.2*scale);
      cam.phi = clamp(cam.phi, 0.15, Math.PI-0.15);
    } else {
      const cp = camPosFromSpherical();
      const forward = v3_norm(v3_sub(cam.target, cp));
      const right = v3_norm(v3_cross(forward, camUp));
      const up = v3_norm(v3_cross(right, forward));
      const panScale = cam.radius * 0.0014;
      cam.target = v3_add(cam.target, v3_add(v3_mul(right, -dx*panScale), v3_mul(up, dy*panScale)));
    }
  });
  elMain.addEventListener('pointerup', ()=>{ input.dragging = false; });
  elMain.addEventListener('wheel', (e)=>{
    e.preventDefault();
    const delta = Math.sign(e.deltaY);
    const zoomFactor = 1.0 + delta*0.08;
    cam.radius = clamp(cam.radius * zoomFactor, cam.minR, cam.maxR);
  }, {passive:false});

  // Time, play, and snap
    let bodyMode = 'earth'; // 'earth' | 'person'
let day = 0.0;
  let playing = false;
  let lastT = performance.now();
  const autoplayDaysPerSec = 2.0;

  const nearestPhaseDay=(d)=>{
    let best = phases[0].day, bd=1e9;
    for(const p of phases){
      const dd = Math.abs(d - p.day);
      if(dd < bd){ bd=dd; best=p.day; }
    }
    return best;
  };
  const phaseNameForDay=(d)=>{
    const t = ((d % SYNODIC_DAYS) + SYNODIC_DAYS) % SYNODIC_DAYS / SYNODIC_DAYS;
    const i = (Math.round(t*8) % 8 + 8) % 8;
    return phases[i].name;
  };
  const updateUI=()=>{
    elDayTxt.textContent = "Day " + day.toFixed(1);
    const nm = phaseNameForDay(day);
    elPhaseTxt.textContent = nm;
    elPhaseChip.textContent = nm;
  };

  elSlider.addEventListener('input', ()=>{ day = parseFloat(elSlider.value); updateUI(); });
  elPlayBtn.addEventListener('click', ()=>{
    playing = !playing;
    elPlayTxt.textContent = playing ? "Pause" : "Play";
    elPlayBtn.querySelector('.tri').style.transform = playing ? "rotate(90deg)" : "none";
    lastT = performance.now();
  });
  elCamBtn.addEventListener('click', ()=>{
    setCamPreset(camPreset === CAM_PRESET.SIDE ? CAM_PRESET.TOP : CAM_PRESET.SIDE);
  });

  function syncBodyLabel(){
    elBodyTxt.textContent = (bodyMode === 'earth') ? 'Earth' : 'Person';
  }
  elBodyBtn.addEventListener('click', ()=>{
    bodyMode = (bodyMode === 'earth') ? 'person' : 'earth';
    syncBodyLabel();
  });
  syncBodyLabel();


  window.addEventListener('keydown', (e)=>{
    if(e.code === 'Space'){ e.preventDefault(); elPlayBtn.click(); }
    if(e.code === 'KeyT'){ e.preventDefault(); elCamBtn.click(); }
    if(e.code === 'KeyP'){ e.preventDefault(); elBodyBtn.click(); }
    if(e.code === 'ArrowLeft'){ day = Math.max(0, day-0.1); elSlider.value = day; updateUI(); }
    if(e.code === 'ArrowRight'){ day = Math.min(SYNODIC_DAYS, day+0.1); elSlider.value = day; updateUI(); }
  }, {passive:false});

  const resizeAll=()=>{ Rmain.resize(); Rinset.resize(); };
  window.addEventListener('resize', resizeAll);
  resizeAll();

  function frame(tNow){
    const dt = Math.min(0.05, (tNow - lastT)/1000);
    lastT = tNow;

    if(playing){
      day += dt * autoplayDaysPerSec;
      if(day > SYNODIC_DAYS) day -= SYNODIC_DAYS;
      const snap = nearestPhaseDay(day);
      const dist = Math.abs(day - snap);
      if(dist < 0.20) day = day + (snap - day) * 0.12;
      elSlider.value = day.toFixed(2);
      updateUI();
    }

    // Moon orbit: X–Z plane (in/out of screen), CCW from +Y
    // Coordinate convention (as agreed):
    // Earth at (0,0,0), Sun at +X, viewer at +Z, top of screen +Y.
    // Moon orbits in the X–Z plane. Day 0 (New Moon) is at +X.
    // As day increases, we want the Moon to go *behind* the Earth first
    // (negative Z), then to -X (Full Moon), then in front (+Z) at Last Quarter.
    // That corresponds to a negative angle (clockwise when viewed from +Y).
    const theta = -(day / SYNODIC_DAYS) * Math.PI * 2.0;
    moonPos.x = Math.cos(theta) * orbitR;
    moonPos.y = 0.0;
    moonPos.z = Math.sin(theta) * orbitR;
    const md = v3_norm(v3_sub(moonPos, earthPos));
    moonDir.x = md.x; moonDir.y = md.y; moonDir.z = md.z;

    // Earth texture rotation: 1 rev/day, reversed
    const texRot = day * (Math.PI*2.0);

    // Camera matrices
    const camPos = camPosFromSpherical();
const Vmain = m4_lookAt(camPos, cam.target, camUp);
    const Pmain = m4_persp(55*Math.PI/180, Rmain.w/Rmain.h, 0.1, 80.0);

    // Models
    let Mearth = M4(); m4_ident(Mearth);
    Mearth = m4_translate(Mearth, earthPos);
    Mearth = m4_scale(Mearth, V3(earthR, earthR, earthR));

    let Mmoon = M4(); m4_ident(Mmoon);
    Mmoon = m4_translate(Mmoon, moonPos);
    Mmoon = m4_scale(Mmoon, V3(moonR, moonR, moonR));

    let Msun = M4(); m4_ident(Msun);
    Msun = m4_translate(Msun, sunPos);
    Msun = m4_scale(Msun, V3(sunR, sunR, sunR));

    // Main render
    Rmain.clear();
    Rmain.drawSphere(Rmain.progSun, Msun, Vmain, Pmain, {uCamPos: camPos});
    if(bodyMode === 'earth'){
      Rmain.drawSphere(Rmain.progEarth, Mearth, Vmain, Pmain, {uSunPos:sunPos, uCamPos:camPos, uMoonDir:moonDir, uTexRot:texRot});
    }else{
      // Procedural cartoon "person" made from scaled spheres (ellipsoids).
      const yaw = (-Math.atan2(moonDir.z, moonDir.x) + Math.PI/2); // face +X toward moon direction projected in XZ (fixed dir + 90°)
      const pr = earthR;
      const skin = [1.0, 0.82, 0.70];
      const shirt = [0.98, 0.98, 0.98];
      const shorts = [0.45, 0.55, 0.72];
      const shoes = [0.96, 0.96, 0.98];
      const hair = [0.62, 0.18, 0.12];
      const dress = [0.72, 0.48, 0.88];
      const lips1 = [0.0, 0.0, 0.0];
      const lips2 = [0.2, 0.1, 0.1];


      const drawPart = (tx,ty,tz, sx,sy,sz, col)=>{
        let M = M4(); m4_ident(M);
        M = m4_translate(M, earthPos);
        M = m4_rotY(M, yaw);
        M = m4_translate(M, V3(tx,ty,tz));
        M = m4_scale(M, V3(sx,sy,sz));
        Rmain.drawSphere(Rmain.progSolid, M, Vmain, Pmain, {uSunPos:sunPos, uCamPos:camPos, uBaseColor:col});
      };

      // Head / hair
      drawPart(0, 1.20*pr, 0, 0.55*pr, 0.55*pr, 0.55*pr, skin);
      drawPart(0, 1.28*pr, -0.06*pr, 0.62*pr, 0.60*pr, 0.54*pr, hair);
      drawPart(0.65*pr, 1.5*pr, -0.06*pr, 0.3*pr, 0.3*pr, 0.3*pr, hair);
      drawPart(-0.65*pr, 1.5*pr, -0.06*pr, 0.3*pr, 0.3*pr, 0.3*pr, hair);

      // Face
      // NOTE: the procedural character's "forward" axis is +Z in its local space
      // (because of the +90° yaw offset used to aim the model at the Moon).
      // Place facial features on the +Z side of the head so they are visible in the
      // default SIDE view (camera is +Z looking toward origin).
      const eyeWhite = [1.00, 1.00, 1.00];
      const eyeBlack = [0.08, 0.08, 0.10];
      // Eyes (two whites + pupils)
      const faceZ  = 0.4*pr;           // sit slightly above the head surface
      const eyeX   = 0.16*pr;
      const eyeY   = 1.1*pr;
      drawPart(  eyeX, eyeY, faceZ + 0.07*pr, 0.08*pr, 0.08*pr, 0.08*pr, eyeBlack);
      drawPart( -eyeX, eyeY, faceZ + 0.07*pr, 0.08*pr, 0.08*pr, 0.08*pr, eyeBlack);


      // Torso / dress
drawPart(0, 0.25*pr, 0, 0.58*pr, 0.78*pr, 0.38*pr, dress);

      // Nose
drawPart(0, 1.0*pr, 0.2*pr, 0.2*pr, 0.1*pr, 0.38*pr, skin);

      // Mouth
drawPart(0, 0.91*pr, 0.43*pr, 0.15*pr, 0.06*pr, 0.01*pr, lips1);
// drawPart(0, 0.87*pr, 0.091*pr, 0.26*pr, 0.06*pr, 0.36*pr, lips2);


// Shirt emblem (a star-ish dot)
const emblem = [0.98, 0.55, 0.68];
drawPart(0.32*pr, 0.99*pr, 0.00*pr, 0.12*pr, 0.12*pr, 0.12*pr, emblem);

// Simple dress / skirt (three layered ellipsoids to feel wider at the bottom)

drawPart(0, -0.40*pr, 0, 0.58*pr, 0.32*pr, 0.45*pr, shirt);
drawPart(0, -0.65*pr, 0, 0.80*pr, 0.45*pr, 0.60*pr, dress);
drawPart(0, -0.9*pr, 0, 0.90*pr, 0.15*pr, 0.70*pr, dress);
drawPart(0, -0.94*pr, 0, 0.95*pr, 0.1*pr, 0.75*pr, dress);

      // Arms (simple) — moved inward by ~another 1/2 arm width
      drawPart(-0.5*pr, 0.275*pr, 0, 0.18*pr, 0.55*pr, 0.18*pr, skin);
      drawPart( 0.5*pr, 0.275*pr, 0, 0.18*pr, 0.55*pr, 0.18*pr, skin);
      drawPart(-0.45*pr, -0.275*pr, 0.1*pr, 0.18*pr, 0.18*pr, 0.18*pr, skin);
      drawPart( 0.45*pr, -0.275*pr, 0.1*pr, 0.18*pr, 0.18*pr, 0.18*pr, skin);

      // Legs
      drawPart(-0.20*pr, -0.85*pr, 0, 0.18*pr, 0.60*pr, 0.18*pr, skin);
      drawPart( 0.20*pr, -0.85*pr, 0, 0.18*pr, 0.60*pr, 0.18*pr, skin);

      // Shoes
      drawPart(-0.20*pr, -1.30*pr, 0.10*pr, 0.22*pr, 0.14*pr, 0.35*pr, shoes);
      drawPart( 0.20*pr, -1.30*pr, 0.10*pr, 0.22*pr, 0.14*pr, 0.35*pr, shoes);
    }
    Rmain.drawSphere(Rmain.progMoon,  Mmoon, Vmain, Pmain, {uSunPos:sunPos, uCamPos:camPos});

    // Inset: what Earth would see.
    // Render Moon centered at origin and express camera/sun positions relative to the Moon.
    Rinset.clear();
    const earthFromMoon = v3_sub(earthPos, moonPos);
    const sunFromMoon   = v3_sub(sunPos,  moonPos);
    const insetEyeDir   = v3_norm(earthFromMoon);
    const insetSunDir   = v3_norm(sunFromMoon);
    const insetEye      = v3_mul(insetEyeDir, 4.8);
    const insetSunPos   = v3_mul(insetSunDir, 8.0);
    const Vinset        = m4_lookAt(insetEye, V3(0,0,0), V3(0,1,0));
    const Pinset        = m4_persp(28*Math.PI/180, Rinset.w/Rinset.h, 0.1, 80.0);
    let MinsetMoon = M4(); m4_ident(MinsetMoon);
    MinsetMoon = m4_scale(MinsetMoon, V3(moonR*1.05, moonR*1.05, moonR*1.05));
    Rinset.drawSphere(Rinset.progMoon, MinsetMoon, Vinset, Pinset, {uSunPos:insetSunPos, uCamPos:insetEye});

    requestAnimationFrame(frame);
  }

  setCamPreset(CAM_PRESET.SIDE);
  updateUI();
  elSlider.value = day.toFixed(2);
  requestAnimationFrame(frame);

  // DEBUG: expose live positions to DevTools
  const __MP_DBG_clone = (v) => {
    if(!v) return null;
    if(v.x !== undefined) return ({x:v.x, y:v.y, z:v.z});
    if(v.length === 3) return ({x:+v[0], y:+v[1], z:+v[2]});
    return null;
  };
  window.__MP_DEBUG = {
    getSunPos:  () => __MP_DBG_clone(sunPos),
    getMoonPos: () => __MP_DBG_clone(moonPos),
    getCamPos:  () => __MP_DBG_clone(camPosFromSpherical()),
    getEyeTarget: () => {
      const md = v3_norm(v3_sub(moonPos, earthPos));
      return __MP_DBG_clone(md);
    },
    getEarthRotation: () => {
      const earthSpinRad = day * 2*Math.PI;
      const texRotRad = -earthSpinRad;
      return { day, earthSpinRad, texRotRad };
    },
  };

})();
