// src/render.js  (V3)
// ============================================================

import { NX, NY, N } from "./lbm.js";

// ---------------------------------------------------------
// COLORMAP — cubic Hermite, 512 steps (kept from V2)
// ---------------------------------------------------------
const COLORMAP_STEPS = 512;
const STOPS = [
  [0,   0,   180],
  [0,   210, 255],
  [0,   255, 80 ],
  [255, 220, 0  ],
  [255, 20,  0  ],
];
const COLORMAP = buildColormap(COLORMAP_STEPS);
const drawnOverlay = new Uint8Array(N);

function cubicHermite(a, b, c, d, t) {
  const t2 = t * t, t3 = t2 * t;
  return (
    (-0.5*a + 1.5*b - 1.5*c + 0.5*d) * t3 +
    (a - 2.5*b + 2*c - 0.5*d) * t2 +
    (-0.5*a + 0.5*c) * t +
    b
  );
}

function buildColormap(steps) {
  const table = new Uint8ClampedArray(steps * 3);
  const nStops = STOPS.length;
  for (let s = 0; s < steps; s++) {
    const t      = s / (steps - 1);
    const scaled = t * (nStops - 1);
    const i1     = Math.min(Math.floor(scaled), nStops - 2);
    const local  = scaled - i1;
    const i0     = Math.max(i1 - 1, 0);
    const i2     = Math.min(i1 + 1, nStops - 1);
    const i3     = Math.min(i1 + 2, nStops - 1);
    for (let c = 0; c < 3; c++) {
      table[s * 3 + c] = Math.round(
        cubicHermite(STOPS[i0][c], STOPS[i1][c], STOPS[i2][c], STOPS[i3][c], local)
      );
    }
  }
  return table;
}

const SOLID_R = 40, SOLID_G = 40, SOLID_B = 45;

// ---------------------------------------------------------
// SCRATCH BUFFERS — allocated once, reused every frame
// ---------------------------------------------------------
const scalarField = new Float32Array(N);
const scalarTmp   = new Float32Array(N);
const uxS         = new Float32Array(N);   
const uyS         = new Float32Array(N);
const uxT         = new Float32Array(N);  
const uyT         = new Float32Array(N);

// ---------------------------------------------------------
// Renderer(displayCanvas)
// ---------------------------------------------------------
export function Renderer(displayCanvas) {
  this.buf     = document.createElement('canvas');
  this.buf.width  = NX;
  this.buf.height = NY;
  this.bufCx   = this.buf.getContext('2d');
  this.img     = this.bufCx.createImageData(NX, NY);

  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  this.dpr   = dpr;
  displayCanvas.width  = displayCanvas.clientWidth  * dpr || NX * 4;
  displayCanvas.height = displayCanvas.clientHeight * dpr || NY * 4;
  this.ctx   = displayCanvas.getContext('2d');
  this.W     = displayCanvas.width;
  this.H     = displayCanvas.height;

  this.smoothedScale = 0.2;
}

// ---------------------------------------------------------
// Renderer.prototype.draw(fields, solid, rho, options)
// ---------------------------------------------------------
Renderer.prototype.draw = function (fields, solid, rho, options, paintedNodes) {
  const mode      = (options && options.mode)      || 'speed';
  const autoScale = (options && options.autoScale) || false;
  const showVectors  = (options && options.showVectors)  || false;  // ← add
  let   scale     = (options && options.scale)     || 0.2;
  const { ux, uy } = fields;
  const d = this.img.data;

  // ---- 1. Pre-smooth velocity for vorticity only ----
  if (mode === 'vort') smoothVelocity(ux, uy, solid, 2);

  // ---- 2. Extract scalar field ----
  for (let n = 0; n < N; n++) {
    if (solid[n]) { scalarField[n] = -1; continue; }

    if (mode === 'vort') {
      const x = n % NX, y = (n / NX) | 0;
      if (x > 0 && x < NX-1 && y > 0 && y < NY-1) {
        const curl = (uyS[n+1] - uyS[n-1]) - (uxS[n+NX] - uxS[n-NX]);
        scalarField[n] = curl;
      } else {
        scalarField[n] = 0;
      }
    } else if (mode === 'press') {
      const cp = (rho[n] - 1.0) / 3.0;
      scalarField[n] = cp;
    } else {
      scalarField[n] = Math.sqrt(ux[n]*ux[n] + uy[n]*uy[n]);
    }
  }

  // ---- 3. Fill solid cells from fluid neighbours ----
  fillSolids(scalarField, solid);

  // ---- 4. Two passes of [1,2,1]/4 scalar blur ----
  blurField(scalarField, scalarTmp, solid, 2);

  // ---- 5. Auto-scale ----
  if (autoScale && mode === 'speed') {
    let frameMax = 0;
    for (let n = 0; n < N; n++) {
      if (!solid[n] && scalarField[n] > frameMax) frameMax = scalarField[n];
    }
    this.smoothedScale = 0.95 * this.smoothedScale + 0.05 * frameMax;
    scale = this.smoothedScale * 0.9 || scale;
  }

  // ---- 6. Map scalar → pixels ----
  const vortGain  = 15;
  const pressGain = 0.55;
  const invScale  = 1 / scale;

  for (let n = 0; n < N; n++) {
    const p = n * 4;

    if (solid[n] && !(paintedNodes && paintedNodes.has(n))) {
      d[p] = SOLID_R; d[p+1] = SOLID_G; d[p+2] = SOLID_B; d[p+3] = 255;
      continue;
    }
    let t;
    if (mode === 'vort') {
      t = scalarField[n] * vortGain;
      t = (Math.max(-1, Math.min(1, t)) + 1) * 0.5;
    } else if (mode === 'press') {
      t = scalarField[n] * pressGain;
      t = (Math.max(-1, Math.min(1, t)) + 1) * 0.5;
    } else {
      t = scalarField[n] * invScale;
      if (t < 0) t = 0; else if (t > 1) t = 1;
    }

    const step = (t * (COLORMAP_STEPS - 1)) | 0;
    const c    = step * 3;
    d[p] = COLORMAP[c]; d[p+1] = COLORMAP[c+1]; d[p+2] = COLORMAP[c+2]; d[p+3] = 255;
  }

  // ---- 7. Write small buffer, GPU-upscale to display canvas ----
  this.bufCx.putImageData(this.img, 0, 0);
  this.ctx.imageSmoothingEnabled = true;
  this.ctx.imageSmoothingQuality = 'high';
  this.ctx.drawImage(this.buf, 0, 0, this.W, this.H);


  // ---- 8. Vector field overlay ----
  if (showVectors && mode === 'speed') {
    this.drawVectorField(fields.ux, fields.uy, solid, {
      spacing: 8,
      scale,
      maxLen: 1.2,
    });
  }

};
// ---------------------------------------------------------
// smoothVelocity(ux, uy, solid, passes)
// ---------------------------------------------------------
function smoothVelocity(ux, uy, solid, passes) {
  uxS.set(ux); uyS.set(uy);
  for (let p = 0; p < passes; p++) {
    for (let y = 1; y < NY-1; y++) {
      const row = y * NX;
      for (let x = 1; x < NX-1; x++) {
        const n = row + x;
        if (solid[n]) { uxT[n] = uxS[n]; uyT[n] = uyS[n]; continue; }
        uxT[n] = (uxS[n-1] + 2*uxS[n] + uxS[n+1]) * 0.25;
        uyT[n] = (uyS[n-1] + 2*uyS[n] + uyS[n+1]) * 0.25;
      }
    }
    for (let y = 1; y < NY-1; y++) {
      const row = y * NX;
      for (let x = 1; x < NX-1; x++) {
        const n = row + x;
        if (solid[n]) { uxS[n] = uxT[n]; uyS[n] = uyT[n]; continue; }
        uxS[n] = (uxT[n-NX] + 2*uxT[n] + uxT[n+NX]) * 0.25;
        uyS[n] = (uyT[n-NX] + 2*uyT[n] + uyT[n+NX]) * 0.25;
      }
    }
  }
}

// ---------------------------------------------------------
// fillSolids(field, solid)
// ---------------------------------------------------------
function fillSolids(field, solid) {
  for (let y = 1; y < NY-1; y++) {
    for (let x = 1; x < NX-1; x++) {
      const n = y * NX + x;
      if (!solid[n]) continue;
      let sum = 0, cnt = 0;
      if (!solid[n-1])   { sum += field[n-1];   cnt++; }
      if (!solid[n+1])   { sum += field[n+1];   cnt++; }
      if (!solid[n-NX])  { sum += field[n-NX];  cnt++; }
      if (!solid[n+NX])  { sum += field[n+NX];  cnt++; }
      if (!solid[n-NX-1]){ sum += field[n-NX-1]; cnt++; }
      if (!solid[n-NX+1]){ sum += field[n-NX+1]; cnt++; }
      if (!solid[n+NX-1]){ sum += field[n+NX-1]; cnt++; }
      if (!solid[n+NX+1]){ sum += field[n+NX+1]; cnt++; }
      field[n] = cnt > 0 ? sum / cnt : 0;
    }
  }
}

// ---------------------------------------------------------
// blurField(field, tmp, solid, passes)
// ---------------------------------------------------------
function blurField(field, tmp, solid, passes) {
  for (let p = 0; p < passes; p++) {
    // horizontal
    for (let y = 0; y < NY; y++) {
      const row = y * NX;
      tmp[row] = field[row];
      tmp[row + NX - 1] = field[row + NX - 1];
      for (let x = 1; x < NX-1; x++) {
        const n = row + x;
        tmp[n] = (field[n-1] + 2*field[n] + field[n+1]) * 0.25;
      }
    }
    // vertical
    for (let x = 0; x < NX; x++) {
      field[x] = tmp[x];
      field[(NY-1)*NX + x] = tmp[(NY-1)*NX + x];
      for (let y = 1; y < NY-1; y++) {
        const n = y * NX + x;
        field[n] = (tmp[n-NX] + 2*tmp[n] + tmp[n+NX]) * 0.25;
      }
    }
  }
}
let bodyOutlineCache = null;
let bodyCanvasCtx   = null;

export function setBodyCanvas(canvas) {
  bodyCanvasCtx = canvas.getContext('2d');
}
export function drawBodyVector(params) {
  if (!bodyCanvasCtx) return;
  const ctx = bodyCanvasCtx;
  const W   = ctx.canvas.width;
  const H   = ctx.canvas.height;
  const sx  = W / NX;
  const sy  = H / NY;

  ctx.clearRect(0, 0, W, H);
  ctx.save();
  ctx.shadowColor = 'rgba(130,190,255,0.5)';
  ctx.shadowBlur  = 14 * Math.min(sx, sy);

  if (params.type === 'cylinder') {
    const pad = 1.5;
    const rx = params.radius * sx;
    const ry = params.radius * sy;

    const grad = ctx.createRadialGradient(
      params.cx * sx, params.cy * sy, 0,
      params.cx * sx, params.cy * sy, rx
    );
    grad.addColorStop(0,   '#48506a');
    grad.addColorStop(0.5, '#1a1f2e');
    grad.addColorStop(1,   '#0a0d16');

    ctx.beginPath();
    ctx.ellipse(params.cx * sx, params.cy * sy, rx, ry, 0, 0, Math.PI * 2);
    ctx.fillStyle = grad;
    ctx.fill();
    ctx.strokeStyle = 'rgba(175,198,230,0.8)';
    ctx.lineWidth   = 1.5;
    ctx.stroke();

  } else if (params.type === 'airfoil') {
    const pts = nacaContourExpanded(params, 0.5);

    ctx.beginPath();
    ctx.moveTo(pts[0][0] * sx, pts[0][1] * sy);
    for (let k = 1; k < pts.length; k++) {
      ctx.lineTo(pts[k][0] * sx, pts[k][1] * sy);
    }
    ctx.closePath();
    const grad = ctx.createLinearGradient(
      0, (params.cy - params.chord * 0.25) * sy,
      0, (params.cy + params.chord * 0.25) * sy
    );
    grad.addColorStop(0,    '#48506a');
    grad.addColorStop(0.45, '#1a1f2e');
    grad.addColorStop(1,    '#0a0d16');

    ctx.fillStyle   = grad;
    ctx.strokeStyle = 'rgba(175,198,230,0.8)';
    ctx.lineWidth   = 1.5;
    ctx.fill();
    ctx.stroke();

  } else if (params.type === 'blocks') {
    const { blocks, blockSize } = params;
    const pad  = 1.0;  
    const half = blockSize / 2 + pad;

    for (const [bx, by] of blocks) {
      const pw = half * 2 * sx;
      const ph = half * 2 * sy;
      const px = bx * sx - half * sx;
      const py = by * sy - half * sy;
      const r  = Math.min(pw, ph) * 0.2;

      const grad = ctx.createRadialGradient(
        bx * sx, by * sy, 0,
        bx * sx, by * sy, Math.max(pw, ph) * 0.7
      );
      grad.addColorStop(0,   '#48506a');
      grad.addColorStop(0.5, '#1a1f2e');
      grad.addColorStop(1,   '#0a0d16');

      ctx.beginPath();
      ctx.roundRect(px, py, pw, ph, r);
      ctx.fillStyle = grad;
      ctx.fill();

      ctx.strokeStyle = 'rgba(175,198,230,0.8)';
      ctx.lineWidth   = 1.2;
      ctx.stroke();
    }
  }
  ctx.restore();
}

function nacaContourExpanded(params, padLattice) {
  const raw = nacaContour(params);        
  const n   = raw.length;
  const out = [];

  for (let i = 0; i < n; i++) {
    const prev = raw[(i - 1 + n) % n];
    const next = raw[(i + 1) % n];
    const tx = next[0] - prev[0];
    const ty = next[1] - prev[1];
    const len = Math.sqrt(tx * tx + ty * ty) || 1;
    const nx = -ty / len;
    const ny =  tx / len;

    out.push([
      raw[i][0] + nx * padLattice,
      raw[i][1] + ny * padLattice,
    ]);
  }

  return out;
}

function nacaContour(params) {
  const { cx, cy, chord, alpha } = params;
  const t    = 0.12;
  const cosA = Math.cos(alpha);
  const sinA = Math.sin(alpha);
  const nPts = 120;
  const upper = [], lower = [];

  for (let s = 0; s <= nPts; s++) {
    const xn = 0.5 * (1 - Math.cos(Math.PI * s / nPts));
    let yt = (t / 0.2) * chord * (
       0.2969 * Math.sqrt(xn)
      - 0.1260 * xn
      - 0.3516 * xn * xn
      + 0.2843 * xn * xn * xn
      - 0.1015 * xn * xn * xn * xn
    );
    if (xn >= 1.0) yt = 0;
    const lx  =  (xn - 0.5) * chord;
    const lyU =  yt;
    const lyL = -yt;
    upper.push([
      cx + lx * cosA - lyU * sinA,
      cy + lx * sinA + lyU * cosA,
    ]);
    lower.push([
      cx + lx * cosA - lyL * sinA,
      cy + lx * sinA + lyL * cosA,
    ]);
  }
  return [...upper, ...[...lower].reverse()];
}




// ---------------------------------------------------------
// Renderer.prototype.drawVectorField(ux, uy, solid, options)
// ---------------------------------------------------------
// Draws velocity arrows on a coarse subgrid over the display
// canvas. Arrows are scaled by local speed and colored white
// with opacity proportional to magnitude so slow/recirculating
// regions read as faint and fast regions read as bright.
// Drawn on the display canvas (not the physics buffer) so
// arrow geometry stays sharp at any zoom level.
// Params:
//   ux, uy  : Float32Array(N) velocity components
//   solid   : Uint8Array(N)
//   options : {
//     spacing  : arrow grid spacing in lattice units (default 8)
//     scale    : same colormap scale as draw() — used to
//                normalize arrow length so the longest arrow
//                fills ~1 grid cell at the max expected speed
//     maxLen   : max arrow length as fraction of grid spacing
//                in display pixels (default 0.85)
//   }
Renderer.prototype.drawVectorField = function (ux, uy, solid, options) {
  const spacing = (options && options.spacing) || 8;
  const scale   = (options && options.scale)   || 0.15;
  const maxLen  = (options && options.maxLen)  || 0.85;

  const ctx = this.ctx;
  const W   = this.W;
  const H   = this.H;

  // Pixel size of one lattice cell on the display canvas
  const cellW = W / NX;
  const cellH = H / NY;

  // Max arrow length in display pixels — one grid spacing × maxLen
  const maxPx = Math.min(cellW, cellH) * spacing * maxLen;

  ctx.save();
  ctx.lineCap = 'round';

  for (let gy = Math.floor(spacing / 2); gy < NY; gy += spacing) {
    for (let gx = Math.floor(spacing / 2); gx < NX; gx += spacing) {
      const n = gy * NX + gx;
      if (solid[n]) continue;

      const u = ux[n];
      const v = uy[n];
      const spd = Math.sqrt(u * u + v * v);
      if (spd < 1e-6) continue;

      // Normalized speed [0,1] relative to colormap scale
      const t = Math.min(spd / scale, 1.0);

      // Arrow length in display pixels — proportional to speed
      const len = t * maxPx;

      // Arrow direction in display space
      // Note: uy positive = downward in screen coords (y increases downward)
      const dx = (u / spd) * len;
      const dy = (v / spd) * len;  // no flip needed — grid y=0 is top

      // Arrow base position in display pixels
      const px = (gx + 0.5) * cellW;
      const py = (gy + 0.5) * cellH;

      // Tip position
      const tx = px + dx;
      const ty = py + dy;

      // Opacity scales with speed so near-zero vectors are invisible
      const alpha = 0.25 + 0.75 * t;
      ctx.strokeStyle = `rgba(255,255,255,${alpha.toFixed(2)})`;
      ctx.lineWidth   = Math.max(0.8, cellW * 0.18);

      // Shaft
      ctx.beginPath();
      ctx.moveTo(px, py);
      ctx.lineTo(tx, ty);
      ctx.stroke();

      // Arrowhead — two short lines branching back from the tip
      if (len > 3) {
        const headLen   = Math.min(len * 0.38, cellW * 1.2);
        const headAngle = Math.PI / 6; // 30 degrees
        const angle     = Math.atan2(dy, dx);

        ctx.beginPath();
        ctx.moveTo(tx, ty);
        ctx.lineTo(
          tx - headLen * Math.cos(angle - headAngle),
          ty - headLen * Math.sin(angle - headAngle)
        );
        ctx.moveTo(tx, ty);
        ctx.lineTo(
          tx - headLen * Math.cos(angle + headAngle),
          ty - headLen * Math.sin(angle + headAngle)
        );
        ctx.stroke();
      }
    }
  }

  ctx.restore();
};