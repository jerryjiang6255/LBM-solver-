// src/render.js
// ============================================================
// Field renderer only: scalar field, vector overlay, grid,
// and distribution overlay.
// Compare-mode ready: no global solver state imports.
// Updated for global fidelity presets: renderer owns NX/NY/N/Q
// and per-instance scratch buffers.
// ============================================================

import { Q } from "./lbm.js";

// ---------------------------------------------------------
// COLORMAP — cubic Hermite, 512 steps
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

function cubicHermite(a, b, c, d, t) {
  const t2 = t * t, t3 = t2 * t;
  return (
    (-0.5 * a + 1.5 * b - 1.5 * c + 0.5 * d) * t3 +
    (a - 2.5 * b + 2 * c - 0.5 * d) * t2 +
    (-0.5 * a + 0.5 * c) * t +
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
// Renderer(displayCanvas, NX, NY)
// ---------------------------------------------------------
export function Renderer(displayCanvas, NX, NY) {
  this.NX = NX;
  this.NY = NY;
  this.N  = NX * NY;
  this.Q  = Q;

  this.buf = document.createElement("canvas");
  this.buf.width  = NX;
  this.buf.height = NY;
  this.bufCx = this.buf.getContext("2d");
  this.img   = this.bufCx.createImageData(NX, NY);

  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  this.dpr = dpr;

  displayCanvas.width  = displayCanvas.clientWidth  * dpr || NX * 4;
  displayCanvas.height = displayCanvas.clientHeight * dpr || NY * 4;

  this.ctx = displayCanvas.getContext("2d");
  this.W   = displayCanvas.width;
  this.H   = displayCanvas.height;

  // -------------------------------------------------------
  // Per-instance scratch buffers
  // -------------------------------------------------------
  this.scalarField = new Float32Array(this.N);
  this.scalarTmp   = new Float32Array(this.N);
  this.uxS         = new Float32Array(this.N);
  this.uyS         = new Float32Array(this.N);
  this.uxT         = new Float32Array(this.N);
  this.uyT         = new Float32Array(this.N);

  this.smoothedScale = 0.2;
}

// ---------------------------------------------------------
// Renderer.prototype.draw(fields, solid, rho, options, paintedNodes)
// fields = { ux, uy, f? }
// ---------------------------------------------------------
Renderer.prototype.draw = function (fields, solid, rho, options, paintedNodes) {
  const NX = this.NX;
  const NY = this.NY;
  const N  = this.N;

  const mode        = (options && options.mode) || "speed";
  const autoScale   = (options && options.autoScale) || false;
  const showVectors = (options && options.showVectors) || false;
  let scale         = (options && options.scale) || 0.2;

  const { ux, uy, f: fField } = fields;
  const d = this.img.data;

  const showGrid = (options && options.showGrid) || false;
  const gridFine = (options && options.gridFine !== undefined) ? options.gridFine : true;
  const showDist = (options && options.showDist) || false;

  if (mode === "qcrit") this.smoothVelocity(ux, uy, solid, 2);

  for (let n = 0; n < N; n++) {
    if (solid[n]) {
      this.scalarField[n] = -1;
      continue;
    }

    if (mode === "qcrit") {
      const x = n % NX;
      const y = (n / NX) | 0;

      if (x > 0 && x < NX - 1 && y > 0 && y < NY - 1) {
        const dudx = (this.uxS[n + 1] - this.uxS[n - 1]) * 0.5;
        const dudy = (this.uxS[n + NX] - this.uxS[n - NX]) * 0.5;
        const dvdx = (this.uyS[n + 1] - this.uyS[n - 1]) * 0.5;
        const dvdy = (this.uyS[n + NX] - this.uyS[n - NX]) * 0.5;
        const omega   = dvdx - dudy;
        const strain2 = dudx * dudx + dvdy * dvdy + 0.5 * (dudy + dvdx) * (dudy + dvdx);
        this.scalarField[n] = 0.5 * (omega * omega - 2.0 * strain2);
      } else {
        this.scalarField[n] = 0;
      }

    } else if (mode === "press") {
      this.scalarField[n] = (rho[n] - 1.0) / 3.0;

    } else {
      this.scalarField[n] = Math.sqrt(ux[n] * ux[n] + uy[n] * uy[n]);
    }
  }

  this.fillSolids(this.scalarField, solid);
  this.blurField(this.scalarField, this.scalarTmp, 2);

  if (autoScale && mode === "speed") {
    let frameMax = 0;
    for (let n = 0; n < N; n++) {
      if (!solid[n] && this.scalarField[n] > frameMax) frameMax = this.scalarField[n];
    }
    this.smoothedScale = 0.95 * this.smoothedScale + 0.05 * frameMax;
    scale = this.smoothedScale * 0.9 || scale;
  }

  const qGain     = 5000;
  const pressGain = 0.55;
  const invScale  = 1 / scale;

  for (let n = 0; n < N; n++) {
    const p = n * 4;

    if (solid[n] && !(paintedNodes && paintedNodes.has(n))) {
      d[p]     = SOLID_R;
      d[p + 1] = SOLID_G;
      d[p + 2] = SOLID_B;
      d[p + 3] = 255;
      continue;
    }

    let t;
    if (mode === "qcrit") {
      t = this.scalarField[n] * qGain;
      t = (Math.max(-1, Math.min(1, t)) + 1) * 0.5;
    } else if (mode === "press") {
      t = this.scalarField[n] * pressGain;
      t = (Math.max(-1, Math.min(1, t)) + 1) * 0.5;
    } else {
      t = this.scalarField[n] * invScale;
      if (t < 0) t = 0;
      else if (t > 1) t = 1;
    }

    const step = (t * (COLORMAP_STEPS - 1)) | 0;
    const c = step * 3;
    d[p]     = COLORMAP[c];
    d[p + 1] = COLORMAP[c + 1];
    d[p + 2] = COLORMAP[c + 2];
    d[p + 3] = 255;
  }

  this.bufCx.putImageData(this.img, 0, 0);
  this.ctx.imageSmoothingEnabled = true;
  this.ctx.imageSmoothingQuality = "high";
  this.ctx.drawImage(this.buf, 0, 0, this.W, this.H);

  if (showVectors && mode === "speed") {
    this.drawVectorField(ux, uy, solid, {
      spacing: 8,
      scale,
      maxLen: 1.2,
    });
  }

  if (showGrid) {
    this.drawGrid(solid, fField, {
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
// Helpers
// ---------------------------------------------------------
Renderer.prototype.smoothVelocity = function (ux, uy, solid, passes) {
  const NX = this.NX;
  const NY = this.NY;

  this.uxS.set(ux);
  this.uyS.set(uy);

  for (let p = 0; p < passes; p++) {
    for (let y = 1; y < NY - 1; y++) {
      const row = y * NX;
      for (let x = 1; x < NX - 1; x++) {
        const n = row + x;
        if (solid[n]) {
          this.uxT[n] = this.uxS[n];
          this.uyT[n] = this.uyS[n];
          continue;
        }
        this.uxT[n] = (this.uxS[n - 1] + 2 * this.uxS[n] + this.uxS[n + 1]) * 0.25;
        this.uyT[n] = (this.uyS[n - 1] + 2 * this.uyS[n] + this.uyS[n + 1]) * 0.25;
      }
    }

    for (let y = 1; y < NY - 1; y++) {
      const row = y * NX;
      for (let x = 1; x < NX - 1; x++) {
        const n = row + x;
        if (solid[n]) {
          this.uxS[n] = this.uxT[n];
          this.uyS[n] = this.uyT[n];
          continue;
        }
        this.uxS[n] = (this.uxT[n - NX] + 2 * this.uxT[n] + this.uxT[n + NX]) * 0.25;
        this.uyS[n] = (this.uyT[n - NX] + 2 * this.uyT[n] + this.uyT[n + NX]) * 0.25;
      }
    }
  }
};

Renderer.prototype.fillSolids = function (field, solid) {
  const NX = this.NX;
  const NY = this.NY;

  for (let y = 1; y < NY - 1; y++) {
    for (let x = 1; x < NX - 1; x++) {
      const n = y * NX + x;
      if (!solid[n]) continue;

      let sum = 0, cnt = 0;
      if (!solid[n - 1])      { sum += field[n - 1];      cnt++; }
      if (!solid[n + 1])      { sum += field[n + 1];      cnt++; }
      if (!solid[n - NX])     { sum += field[n - NX];     cnt++; }
      if (!solid[n + NX])     { sum += field[n + NX];     cnt++; }
      if (!solid[n - NX - 1]) { sum += field[n - NX - 1]; cnt++; }
      if (!solid[n - NX + 1]) { sum += field[n - NX + 1]; cnt++; }
      if (!solid[n + NX - 1]) { sum += field[n + NX - 1]; cnt++; }
      if (!solid[n + NX + 1]) { sum += field[n + NX + 1]; cnt++; }

      field[n] = cnt > 0 ? sum / cnt : 0;
    }
  }
};

Renderer.prototype.blurField = function (field, tmp, passes) {
  const NX = this.NX;
  const NY = this.NY;

  for (let p = 0; p < passes; p++) {
    for (let y = 0; y < NY; y++) {
      const row = y * NX;
      tmp[row] = field[row];
      tmp[row + NX - 1] = field[row + NX - 1];
      for (let x = 1; x < NX - 1; x++) {
        const n = row + x;
        tmp[n] = (field[n - 1] + 2 * field[n] + field[n + 1]) * 0.25;
      }
    }

    for (let x = 0; x < NX; x++) {
      field[x] = tmp[x];
      field[(NY - 1) * NX + x] = tmp[(NY - 1) * NX + x];
      for (let y = 1; y < NY - 1; y++) {
        const n = y * NX + x;
        field[n] = (tmp[n - NX] + 2 * tmp[n] + tmp[n + NX]) * 0.25;
      }
    }
  }
}; 
// ---------------------------------------------------------
// Vector overlay
// ---------------------------------------------------------
Renderer.prototype.drawVectorField = function (ux, uy, solid, options) {
  const NX = this.NX;
  const NY = this.NY;

  const spacing = (options && options.spacing) || 8;
  const scale   = (options && options.scale) || 0.15;
  const maxLen  = (options && options.maxLen) || 0.85;

  const ctx = this.ctx;
  const Wc  = this.W;
  const Hc  = this.H;

  const cellW = Wc / NX;
  const cellH = Hc / NY;
  const maxPx = Math.min(cellW, cellH) * spacing * maxLen;

  ctx.save();
  ctx.lineCap = "round";

  for (let gy = 0; gy < NY; gy += spacing) {
    for (let gx = 0; gx < NX; gx += spacing) {
      const n = gy * NX + gx;
      if (solid[n]) continue;

      const u = ux[n];
      const v = uy[n];
      const spd = Math.sqrt(u * u + v * v);
      if (spd < 1e-6) continue;

      const t = Math.min(spd / scale, 1.0);
      const len = t * maxPx;

      const dx = (u / spd) * len;
      const dy = (v / spd) * len;

      const px = gx * cellW;
      const py = gy * cellH;
      const tx = px + dx;
      const ty = py + dy;

      const alpha = 0.25 + 0.75 * t;
      ctx.strokeStyle = `rgba(255,255,255,${alpha.toFixed(2)})`;
      ctx.lineWidth = Math.max(0.8, cellW * 0.18);

      ctx.beginPath();
      ctx.moveTo(px, py);
      ctx.lineTo(tx, ty);
      ctx.stroke();

      if (len > 3) {
        const headLen   = Math.min(len * 0.38, cellW * 1.2);
        const headAngle = Math.PI / 6;
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
// Grid / distribution overlay
// ---------------------------------------------------------
Renderer.prototype.drawGrid = function (solid, fField, options) {
  const NX = this.NX;
  const NY = this.NY;
  const Q  = this.Q;

  const showFine    = (options && options.showFine !== undefined) ? options.showFine : true;
  const showCoarse  = (options && options.showCoarse !== undefined) ? options.showCoarse : true;
  const coarseStep  = (options && options.coarseStep) || 8;
  const showDist    = (options && options.showDist) || false;
  const distSpacing = (options && options.distSpacing) || 16;
  const distScale   = (options && options.distScale) || 80;

  const ctx   = this.ctx;
  const Wc    = this.W;
  const Hc    = this.H;
  const cellW = Wc / NX;
  const cellH = Hc / NY;

  ctx.save();

  if (showFine) {
    ctx.strokeStyle = "rgba(255,255,255,0.06)";
    ctx.lineWidth = 0.5;
    ctx.beginPath();

    for (let x = 0; x <= NX; x++) {
      const px = x * cellW;
      ctx.moveTo(px, 0);
      ctx.lineTo(px, Hc);
    }
    for (let y = 0; y <= NY; y++) {
      const py = y * cellH;
      ctx.moveTo(0, py);
      ctx.lineTo(Wc, py);
    }
    ctx.stroke();
  }

  if (showCoarse) {
    ctx.strokeStyle = "rgba(255,255,255,0.18)";
    ctx.lineWidth = 0.75;
    ctx.beginPath();

    for (let x = 0; x <= NX; x += coarseStep) {
      const px = x * cellW;
      ctx.moveTo(px, 0);
      ctx.lineTo(px, Hc);
    }
    for (let y = 0; y <= NY; y += coarseStep) {
      const py = y * cellH;
      ctx.moveTo(0, py);
      ctx.lineTo(Wc, py);
    }
    ctx.stroke();
  }

  if (showDist && fField) {
    const CX_D = [0, 1, 0, -1, 0, 1, -1, -1, 1];
    const CY_D = [0, 0, -1, 0, 1, -1, -1, 1, 1];

    const DIR_COLORS = [
      "rgba(255,255,255,0.9)",
      "rgba(100,220,255,0.85)",
      "rgba(100,220,255,0.85)",
      "rgba(100,220,255,0.85)",
      "rgba(100,220,255,0.85)",
      "rgba(255,220,80,0.8)",
      "rgba(255,220,80,0.8)",
      "rgba(255,220,80,0.8)",
      "rgba(255,220,80,0.8)",
    ];

    ctx.lineWidth = Math.max(0.6, cellW * 0.12);
    ctx.lineCap = "round";

    for (let gy = 0; gy < NY; gy += distSpacing) {
      for (let gx = 0; gx < NX; gx += distSpacing) {
        const n = gy * NX + gx;
        if (solid[n]) continue;

        const px = gx * cellW;
        const py = gy * cellH;

        ctx.fillStyle = "rgba(255,255,255,0.6)";
        ctx.beginPath();
        ctx.arc(px, py, Math.max(1, cellW * 0.15), 0, Math.PI * 2);
        ctx.fill();

        for (let i = 0; i < 9; i++) {
          const fi = fField[n * Q + i];
          const armLen = fi * distScale * Math.min(cellW, cellH);

          if (i === 0) {
            const r = Math.min(armLen * 0.5, cellW * 0.4);
            if (r > 0.5) {
              ctx.strokeStyle = DIR_COLORS[0];
              ctx.beginPath();
              ctx.arc(px, py, r, 0, Math.PI * 2);
              ctx.stroke();
            }
            continue;
          }

          const ex = px + CX_D[i] * armLen;
          const ey = py + CY_D[i] * armLen;

          ctx.strokeStyle = DIR_COLORS[i];
          ctx.beginPath();
          ctx.moveTo(px, py);
          ctx.lineTo(ex, ey);
          ctx.stroke();

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
