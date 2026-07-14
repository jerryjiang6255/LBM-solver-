// src/render.js  (V3)
// ============================================================

import { NX, NY, N, f, Q } from "./lbm.js";

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

  const showGrid    = (options && options.showGrid)    || false;
  const gridFine    = (options && options.gridFine     !== undefined) ? options.gridFine : true;
  const showDist    = (options && options.showDist)    || false;

  // ---- 1. Pre-smooth velocity for Q-criterion ----
  if (mode === 'qcrit') smoothVelocity(ux, uy, solid, 2);

  // ---- 2. Extract scalar field ----
  for (let n = 0; n < N; n++) {
    if (solid[n]) { scalarField[n] = -1; continue; }

    if (mode === 'qcrit') {
      const x = n % NX, y = (n / NX) | 0;
      if (x > 0 && x < NX-1 && y > 0 && y < NY-1) {
        const dudx = (uxS[n+1]  - uxS[n-1])  * 0.5;
        const dudy = (uxS[n+NX] - uxS[n-NX]) * 0.5;
        const dvdx = (uyS[n+1]  - uyS[n-1])  * 0.5;
        const dvdy = (uyS[n+NX] - uyS[n-NX]) * 0.5;
        const omega   = dvdx - dudy;
        const strain2 = dudx*dudx + dvdy*dvdy
                      + 0.5*(dudy + dvdx)*(dudy + dvdx);
        scalarField[n] = 0.5 * (omega*omega - 2.0*strain2);
      } else {
        scalarField[n] = 0;
      }

    } else if (mode === 'press') {
      scalarField[n] = (rho[n] - 1.0) / 3.0;

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
  const qGain = 5000;  // tune this — increase if field looks flat
  const pressGain = 0.55;
  const invScale  = 1 / scale;

  for (let n = 0; n < N; n++) {
    const p = n * 4;

    if (solid[n] && !(paintedNodes && paintedNodes.has(n))) {
      d[p] = SOLID_R; d[p+1] = SOLID_G; d[p+2] = SOLID_B; d[p+3] = 255;
      continue;
    }
    let t;
    if (mode === 'qcrit') {
      t = scalarField[n] * qGain;
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


  // ---- 10. Grid + distribution overlay ----
  if (showGrid) {
    this.drawGrid(solid, {
      showFine:    gridFine,
      showCoarse:  true,
      coarseStep:  8,
      showDist:    showDist,
      distSpacing: 16,
      distScale:   48,
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
    } else if (params.type === 'nozzle') {
      const { cy, diameter, length } = params;
      const halfD         = diameter / 2;
      const wallThickness = 4;
      const sx_           = W / NX;
      const sy_           = H / NY;

      const topY1 = (cy - halfD - wallThickness) * sy_;
      const topY2 = (cy - halfD) * sy_;
      const botY1 = (cy + halfD) * sy_;
      const botY2 = (cy + halfD + wallThickness) * sy_;
const endX = (length + 0.9) * sx_; // add 3 lattice units of padding

      // Gradient matching cylinder/airfoil style
      const gradTop = ctx.createLinearGradient(0, topY1, 0, topY2);
      gradTop.addColorStop(0,   '#48506a');
      gradTop.addColorStop(0.5, '#1a1f2e');
      gradTop.addColorStop(1,   '#0a0d16');

      const gradBot = ctx.createLinearGradient(0, botY1, 0, botY2);
      gradBot.addColorStop(0,   '#0a0d16');
      gradBot.addColorStop(0.5, '#1a1f2e');
      gradBot.addColorStop(1,   '#48506a');

      const startX = -10; // bleed well past the left canvas edge
      // Top wall
      ctx.beginPath();
      ctx.rect(startX, topY1, endX - startX, topY2 - topY1);
      ctx.fillStyle = gradTop;
      ctx.fill();
      ctx.strokeStyle = 'rgba(175,198,230,0.8)';
      ctx.lineWidth   = 1.2;
      ctx.stroke();

      // Bottom wall
      ctx.beginPath();
      ctx.rect(startX, botY1, endX - startX, botY2 - botY1);
      ctx.fillStyle = gradBot;
      ctx.fill();
      ctx.strokeStyle = 'rgba(175,198,230,0.8)';
      ctx.lineWidth   = 1.2;
      ctx.stroke();

      // Exit lip — right end only, no left stroke
      ctx.beginPath();
      ctx.moveTo(endX, topY1);
      ctx.lineTo(endX, topY2);
      ctx.moveTo(endX, botY1);
      ctx.lineTo(endX, botY2);
      ctx.strokeStyle = 'rgba(175,198,230,0.8)';
      ctx.lineWidth   = 1.5;
      ctx.stroke();
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

    for (let gy = 0; gy < NY; gy += spacing) {
      for (let gx = 0; gx < NX; gx += spacing) {
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
      const px = gx * cellW;
      const py = gy * cellH;

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






// ---------------------------------------------------------
// Renderer.prototype.drawGrid(solid, options)
// ---------------------------------------------------------
// Draws the D2Q9 lattice grid as thin lines on the display
// canvas. Two levels of detail:
//   - Fine grid  : every lattice cell boundary (1px, faint)
//   - Coarse grid: every `coarseStep` cells (slightly brighter)
// Also optionally draws the f distribution arrows at sampled
// nodes showing the 9 D2Q9 populations as line segments whose
// length is proportional to f[n*Q+i], giving a direct visual
// of the distribution function state at that node.
// Params:
//   solid   : Uint8Array(N)
//   options : {
//     showFine      : bool   — draw every cell boundary (default true)
//     showCoarse    : bool   — draw coarser major gridlines (default true)
//     coarseStep    : int    — major grid every N cells (default 8)
//     showDist      : bool   — draw f distribution arrows (default false)
//     distSpacing   : int    — sample f every N cells (default 16)
//     distScale     : float  — scale factor for f arrow lengths
//   }
Renderer.prototype.drawGrid = function (solid, options) {
  const showFine    = (options && options.showFine    !== undefined) ? options.showFine    : true;
  const showCoarse  = (options && options.showCoarse  !== undefined) ? options.showCoarse  : true;
  const coarseStep  = (options && options.coarseStep) || 8;
  const showDist    = (options && options.showDist)   || false;
  const distSpacing = (options && options.distSpacing)|| 16;
  const distScale   = (options && options.distScale)  || 80;

  const ctx  = this.ctx;
  const W    = this.W;
  const H    = this.H;
  const cellW = W / NX;
  const cellH = H / NY;

  ctx.save();

  // ---- Fine grid: every lattice cell ----
  if (showFine) {
    ctx.strokeStyle = 'rgba(255,255,255,0.06)';
    ctx.lineWidth   = 0.5;
    ctx.beginPath();

    // Vertical lines
    for (let x = 0; x <= NX; x++) {
      const px = x * cellW;
      ctx.moveTo(px, 0);
      ctx.lineTo(px, H);
    }
    // Horizontal lines
    for (let y = 0; y <= NY; y++) {
      const py = y * cellH;
      ctx.moveTo(0, py);
      ctx.lineTo(W, py);
    }
    ctx.stroke();
  }

  // ---- Coarse grid: every coarseStep cells ----
  if (showCoarse) {
    ctx.strokeStyle = 'rgba(255,255,255,0.18)';
    ctx.lineWidth   = 0.75;
    ctx.beginPath();

    for (let x = 0; x <= NX; x += coarseStep) {
      const px = x * cellW;
      ctx.moveTo(px, 0);
      ctx.lineTo(px, H);
    }
    for (let y = 0; y <= NY; y += coarseStep) {
      const py = y * cellH;
      ctx.moveTo(0, py);
      ctx.lineTo(W, py);
    }
    ctx.stroke();
  }

  // ---- Distribution function arrows at sampled nodes ----
  // For each sampled fluid node, draws 9 line segments — one per
  // D2Q9 direction — with length proportional to f[n*Q+i].
  // The rest direction (i=0) is shown as a dot at the node center.
  // This gives a direct visual of the microscopic state: you can
  // see how the distribution is skewed toward the flow direction
  // (more mass in directions aligned with ux/uy) and how it
  // deviates from equilibrium near the obstacle.
  if (showDist) {
    // D2Q9 velocity vectors, same order as lbm.js CX/CY
    const CX_D = [ 0,  1,  0, -1,  0,  1, -1, -1,  1];
    const CY_D = [ 0,  0, -1,  0,  1, -1, -1,  1,  1]; 

    // Equilibrium weights for visual reference baseline
    // (direction thickness scales with how far f deviates from W[i])
    const W_D  = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];

    // Color per direction group:
    //   rest     : white
    //   axis     : cyan
    //   diagonal : yellow
    const DIR_COLORS = [
      'rgba(255,255,255,0.9)',  // 0 rest
      'rgba(100,220,255,0.85)', // 1 E
      'rgba(100,220,255,0.85)', // 2 N
      'rgba(100,220,255,0.85)', // 3 W
      'rgba(100,220,255,0.85)', // 4 S
      'rgba(255,220,80,0.8)',   // 5 NE
      'rgba(255,220,80,0.8)',   // 6 NW
      'rgba(255,220,80,0.8)',   // 7 SW
      'rgba(255,220,80,0.8)',   // 8 SE
    ];

    ctx.lineWidth = Math.max(0.6, cellW * 0.12);
    ctx.lineCap   = 'round';

    for (let gy = 0; gy < NY; gy += distSpacing) {
      for (let gx = 0; gx < NX; gx += distSpacing) {
        const n = gy * NX + gx;
        if (solid[n]) continue;

        // Node center in display pixels
        const px = gx * cellW;
        const py = gy * cellH;

        // Draw a small dot at the node center
        ctx.fillStyle = 'rgba(255,255,255,0.6)';
        ctx.beginPath();
        ctx.arc(px, py, Math.max(1, cellW * 0.15), 0, Math.PI * 2);
        ctx.fill();

        // Draw each of the 9 direction arms
        for (let i = 0; i < 9; i++) {
          const fi = f[n * Q + i];

          // Arrow length: fi * distScale in display pixels
          // Scaled by cell size so it looks reasonable at any zoom
          const armLen = fi * distScale * Math.min(cellW, cellH);

          if (i === 0) {
            // Rest direction — draw a circle whose radius ~ f0
            // so you can see how much mass is "resting"
            const r = Math.min(armLen * 0.5, cellW * 0.4);
            if (r > 0.5) {
              ctx.strokeStyle = DIR_COLORS[0];
              ctx.beginPath();
              ctx.arc(px, py, r, 0, Math.PI * 2);
              ctx.stroke();
            }
            continue;
          }

          // Moving direction — draw a line from node center
          // in the lattice direction, length proportional to fi
          const ex = px + CX_D[i] * armLen;
          const ey = py + CY_D[i] * armLen;

          ctx.strokeStyle = DIR_COLORS[0];
          ctx.beginPath();
          ctx.moveTo(px, py);
          ctx.lineTo(ex, ey);
          ctx.stroke();

          // Small arrowhead at tip if arm is long enough
          if (armLen > cellW * 0.4) {
            const angle   = Math.atan2(CY_D[i], CX_D[i]);
            const headLen = armLen * 0.3;
            const headAng = Math.PI / 5;
            ctx.beginPath();
            ctx.moveTo(ex, ey);
            ctx.lineTo(
              ex - headLen * Math.cos(angle - headAng),
              ey - headLen * Math.sin(angle - headAng)
            );
            ctx.moveTo(ex, ey);
            ctx.lineTo(
              ex - headLen * Math.cos(angle + headAng),
              ey - headLen * Math.sin(angle + headAng)
            );
            ctx.stroke();
          }
        }
      }
    }
  }

  ctx.restore();
};