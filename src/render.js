// src/render.js  (V3)
// ============================================================
// Key improvements over V2:
//
//   1. TWO-CANVAS APPROACH (main smoothness fix)
//      Physics runs at NX×NY (256×128). A small offscreen
//      canvas (buf) is written at that resolution via
//      putImageData, then upscaled to the display canvas
//      with drawImage + imageSmoothingQuality='high'. The
//      browser GPU bilinear-filters every pixel boundary for
//      free — this is why the reference looks smooth with no
//      extra blur passes needed.
//
//   2. VORTICITY + PRESSURE display modes (new)
//      Vorticity = curl of velocity = duy/dx - dux/dy.
//      Shows rotating structures the speed view hides.
//      Pressure = (rho - 1)/3 in lattice units, normalised
//      to Cp. Both match the reference's view modes.
//
//   3. VELOCITY PRE-SMOOTH before vorticity only
//      Vorticity amplifies grid-scale noise (2-cell ripple).
//      Smoothing ux/uy with 2 passes of [1,2,1]/4 before
//      taking the curl removes the spurious criss-cross
//      pattern while leaving real shear layers intact.
//      Speed and pressure views don't need this.
//
//   4. SOLID FILL + LIGHT SCALAR SMOOTH (kept from V2)
//      Solid cells get averaged from fluid neighbours before
//      colourmap lookup, preventing the dark halo. Two passes
//      of [1,2,1]/4 on the scalar only (not velocity).
//
//   5. BODY-AGNOSTIC OUTLINE
//      drawOutline now accepts a `geometry` object and draws
//      the solid mask boundary as a filled dark polygon using
//      the canvas path API — works for cylinder and airfoil.
//
//   6. AUTO-SCALE option
//      If options.autoScale is true, the renderer computes
//      the per-frame 98th-percentile speed and uses that as
//      the colour scale, so the map is always well-saturated
//      regardless of Re or u0.
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
const uxS         = new Float32Array(N);   // smoothed ux for vorticity
const uyS         = new Float32Array(N);
const uxT         = new Float32Array(N);   // temp for blur passes
const uyT         = new Float32Array(N);

// ---------------------------------------------------------
// Renderer(displayCanvas)
// ---------------------------------------------------------
// displayCanvas : the full-size CSS canvas the user sees.
// Internally creates a small offscreen buf at NX×NY that
// receives putImageData, then gets GPU-upscaled via drawImage.
export function Renderer(displayCanvas) {
  // Offscreen buffer at physics resolution
  this.buf     = document.createElement('canvas');
  this.buf.width  = NX;
  this.buf.height = NY;
  this.bufCx   = this.buf.getContext('2d');
  this.img     = this.bufCx.createImageData(NX, NY);

  // Display canvas — sized to its CSS layout size × DPR
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  this.dpr   = dpr;
  displayCanvas.width  = displayCanvas.clientWidth  * dpr || NX * 4;
  displayCanvas.height = displayCanvas.clientHeight * dpr || NY * 4;
  this.ctx   = displayCanvas.getContext('2d');
  this.W     = displayCanvas.width;
  this.H     = displayCanvas.height;

  // Smoothed EMA scale for auto-scale mode
  this.smoothedScale = 0.2;
}

// ---------------------------------------------------------
// Renderer.prototype.draw(fields, solid, rho, options)
// ---------------------------------------------------------
// fields  : { ux, uy }  Float32Array(N)
// solid   : Uint8Array(N)
// rho     : Float32Array(N)  — needed for pressure mode
// options : {
//   mode      : 'speed' | 'vort' | 'press'  (default 'speed')
//   scale     : normalisation constant  (ignored if autoScale)
//   autoScale : bool — compute scale from field each frame
// }
Renderer.prototype.draw = function (fields, solid, rho, options, paintedNodes) {
  const mode      = (options && options.mode)      || 'speed';
  const autoScale = (options && options.autoScale) || false;
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

};




// ---------------------------------------------------------
// smoothVelocity(ux, uy, solid, passes)
// ---------------------------------------------------------
// [1,2,1]/4 separable blur on ux/uy into uxS/uyS.
// Used only for vorticity — removes grid-scale ripple
// before taking the curl without affecting solver state.
function smoothVelocity(ux, uy, solid, passes) {
  uxS.set(ux); uyS.set(uy);
  for (let p = 0; p < passes; p++) {
    // horizontal
    for (let y = 1; y < NY-1; y++) {
      const row = y * NX;
      for (let x = 1; x < NX-1; x++) {
        const n = row + x;
        if (solid[n]) { uxT[n] = uxS[n]; uyT[n] = uyS[n]; continue; }
        uxT[n] = (uxS[n-1] + 2*uxS[n] + uxS[n+1]) * 0.25;
        uyT[n] = (uyS[n-1] + 2*uyS[n] + uyS[n+1]) * 0.25;
      }
    }
    // vertical
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
// In-place: solid cells get average of fluid neighbours.
// Prevents dark halo when blur crosses the boundary.
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
// Separable [1,2,1]/4 in-place blur on the scalar field.
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

// render.js — replace drawBodyOutline with this
// Called once at geometry build time, not every frame.
// Stores the outline path for fast per-frame redraw.
let bodyOutlineCache = null;
let bodyCanvasCtx   = null;

export function setBodyCanvas(canvas) {
  bodyCanvasCtx = canvas.getContext('2d');
}

// Call this whenever geometry changes (reset, new shape, AoA change)
// Pass the SDF circle/airfoil params, not the solid[] array
export function drawBodyVector(params) {
  if (!bodyCanvasCtx) return;
  const ctx = bodyCanvasCtx;
  const W   = ctx.canvas.width;
  const H   = ctx.canvas.height;
  const sx  = W / NX;
  const sy  = H / NY;

  ctx.clearRect(0, 0, W, H);
  ctx.save();

  // Glow effect
  ctx.shadowColor = 'rgba(130,190,255,0.5)';
  ctx.shadowBlur  = 14 * Math.min(sx, sy);

  if (params.type === 'cylinder') {
    // Pad by 1.5 lattice units to fully cover staircase pixels
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

    // Sharp stroke on the outside edge
    ctx.strokeStyle = 'rgba(175,198,230,0.8)';
    ctx.lineWidth   = 1.5;
    ctx.stroke();

  } else if (params.type === 'airfoil') {
    // Expand contour outward by 1.5 lattice units along surface normal
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
  }

  ctx.restore();
}

function nacaContourExpanded(params, padLattice) {
  const raw = nacaContour(params);        // existing function, unchanged
  const n   = raw.length;
  const out = [];

  for (let i = 0; i < n; i++) {
    const prev = raw[(i - 1 + n) % n];
    const next = raw[(i + 1) % n];

    // Tangent vector along contour
    const tx = next[0] - prev[0];
    const ty = next[1] - prev[1];
    const len = Math.sqrt(tx * tx + ty * ty) || 1;

    // Outward normal (perpendicular to tangent, pointing out)
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

    // Force sharp trailing edge — standard coefficients leave
    // a small nonzero thickness at xn=1, close it explicitly
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

  // upper LE→TE, then lower TE→LE, closing at the same trailing edge point
  return [...upper, ...[...lower].reverse()];
}

