// src/render.js  (V2 — enhanced display)
// ============================================================
// Same core as before: converts macroscopic fields to pixels.
// New additions per the display enhancement list:
//   1. Separable [1,2,1]/4 binomial blur on the scalar field
//      (display only — solver arrays untouched)
//   2. Solid cells filled from fluid neighbor average before blur
//      (removes dark halo at obstacle boundary)
//   3. Analytic circle outline drawn as a thin overlay
//   4. Smoother colormap using cubic Hermite interpolation
//      between hue stops instead of linear segments
// ============================================================

import { NX, NY, N } from "./lbm.js";

// ---------------------------------------------------------
// COLORMAP
// ---------------------------------------------------------
// 5 control-point hues (blue->cyan->green->yellow->red),
// sampled at COLORMAP_STEPS via cubic Hermite interpolation
// between adjacent stops. Produces smoother gradients than
// the original piecewise-linear version, especially visible
// at low colormap-scale values where the full range is
// compressed into a narrow speed band.
const COLORMAP_STEPS = 512; // higher resolution, same lookup cost

// Control points: [r, g, b] at t=0, 0.25, 0.5, 0.75, 1.0
const STOPS = [
  [0,   0,   180],  // deep blue
  [0,   210, 255],  // cyan
  [0,   255, 80 ],  // green
  [255, 220, 0  ],  // yellow
  [255, 20,  0  ],  // red
];

const COLORMAP = buildColormap(COLORMAP_STEPS);


function cubicHermite(a, b, c, d, t) {
  // Catmull-Rom tangent blending — smoother than linear
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
    const t = s / (steps - 1);
    // which segment?
    const scaled = t * (nStops - 1);
    const i1 = Math.min(Math.floor(scaled), nStops - 2);
    const local = scaled - i1;

    // 4 control points for Catmull-Rom (clamp at edges)
    const i0 = Math.max(i1 - 1, 0);
    const i2 = Math.min(i1 + 1, nStops - 1);
    const i3 = Math.min(i1 + 2, nStops - 1);

    for (let c = 0; c < 3; c++) {
      table[s * 3 + c] = Math.round(
        cubicHermite(STOPS[i0][c], STOPS[i1][c], STOPS[i2][c], STOPS[i3][c], local)
      );
    }
  }
  return table;
}

// Solid cell display color (dark gray)
const SOLID_R = 55, SOLID_G = 55, SOLID_B = 60;

// ---------------------------------------------------------
// SCRATCH BUFFERS
// ---------------------------------------------------------
// Allocated once, reused every frame.
// scalarField : extracted display values before blur, Float32Array(N)
// blurBuf     : intermediate buffer for separable blur, Float32Array(N)
// filledField : scalarField with solid cells filled, Float32Array(N)
const scalarField  = new Float32Array(N);
const blurBuf      = new Float32Array(N);
const filledField  = new Float32Array(N);

// ---------------------------------------------------------
// Renderer(canvas, cylCx, cylCy, cylRadius)
// ---------------------------------------------------------
// Params:
//   canvas    : HTMLCanvasElement
//   cylCx, cylCy, cylRadius : cylinder geometry in lattice coords,
//     used by drawOutline() to render a clean analytic circle
//     instead of relying on the rasterized solid mask pixels.
export function Renderer(canvas, cylCx, cylCy, cylRadius) {
  canvas.width  = NX;
  canvas.height = NY;
  this.ctx       = canvas.getContext("2d");
  this.imageData = this.ctx.createImageData(NX, NY);
  this.pixels    = this.imageData.data;
  this.cylCx     = cylCx;
  this.cylCy     = cylCy;
  this.cylRadius = cylRadius;
}

// ---------------------------------------------------------
// Renderer.prototype.draw(fields, solid, options)
// ---------------------------------------------------------
// Same interface as before, now with blur + fill + outline.
// Params:
//   fields  : { ux, uy } Float32Array(N) each
//   solid   : Uint8Array(N)
//   options : { mode, scale }
Renderer.prototype.draw = function (fields, solid, options) {
  const mode     = (options && options.mode)  || "speed";
  const scale    = (options && options.scale) || 0.15;
  const invScale = 1 / scale;
  const pixels   = this.pixels;
  const { ux, uy } = fields;

  // 1. Extract raw scalar field
  for (let n = 0; n < N; n++) {
    if (mode === "speed") {
      const u = ux[n], v = uy[n];
      scalarField[n] = Math.sqrt(u * u + v * v);
    } else {
      scalarField[n] = 0;
    }
  }

  // 2. Fill solid cells from fluid neighbor average
  //    (prevents dark halo when blur crosses solid/fluid boundary)
  fillSolidsFromNeighbors(scalarField, solid, filledField);

  // 3. Two passes of separable [1,2,1]/4 binomial blur
  //    applied to filledField -> blurBuf -> filledField
  blurH(filledField, blurBuf);
  blurV(blurBuf, filledField);
  blurH(filledField, blurBuf);
  blurV(blurBuf, filledField);

  // 4. Write pixels from smoothed field
  for (let n = 0; n < N; n++) {
    const p = n * 4;

    if (solid[n]) {
      pixels[p]     = SOLID_R;
      pixels[p + 1] = SOLID_G;
      pixels[p + 2] = SOLID_B;
      pixels[p + 3] = 255;
      continue;
    }

    let t = filledField[n] * invScale;
    if (t < 0) t = 0;
    else if (t > 1) t = 1;

    const step = ((t * (COLORMAP_STEPS - 1)) | 0);
    const c = step * 3;
    pixels[p]     = COLORMAP[c];
    pixels[p + 1] = COLORMAP[c + 1];
    pixels[p + 2] = COLORMAP[c + 2];
    pixels[p + 3] = 255;
  }

  this.ctx.putImageData(this.imageData, 0, 0);

  // 5. Draw analytic cylinder outline on top
  drawOutline(this.ctx, this.cylCx, this.cylCy, this.cylRadius);
};

// ---------------------------------------------------------
// fillSolidsFromNeighbors(src, solid, dst)
// ---------------------------------------------------------
// Copies src into dst, but for solid cells, replaces the value
// with the average of up to 4 orthogonal fluid neighbors
// (or the cell's own value if all neighbors are also solid).
// This is done before blurring so the kernel never smears a
// hard solid=0 boundary into adjacent fluid pixels.
function fillSolidsFromNeighbors(src, solid, dst) {
  for (let y = 0; y < NY; y++) {
    for (let x = 0; x < NX; x++) {
      const n = y * NX + x;
      if (!solid[n]) {
        dst[n] = src[n];
        continue;
      }
      let sum = 0, count = 0;
      if (x > 0       && !solid[n - 1 ]) { sum += src[n - 1];  count++; }
      if (x < NX - 1  && !solid[n + 1 ]) { sum += src[n + 1];  count++; }
      if (y > 0       && !solid[n - NX]) { sum += src[n - NX]; count++; }
      if (y < NY - 1  && !solid[n + NX]) { sum += src[n + NX]; count++; }
      dst[n] = count > 0 ? sum / count : src[n];
    }
  }
}

// ---------------------------------------------------------
// blurH(src, dst) / blurV(src, dst)
// ---------------------------------------------------------
// Separable [1,2,1]/4 binomial blur, horizontal then vertical.
// Applied to the display scalar only — solver arrays are never
// touched. Edge pixels just copy from the nearest valid neighbor
// (clamped boundary), which is unnoticeable given the domain
// edges are covered by the Dirichlet BC color anyway.
function blurH(src, dst) {
  for (let y = 0; y < NY; y++) {
    const row = y * NX;
    dst[row] = (src[row] * 3 + src[row + 1]) / 4;
    for (let x = 1; x < NX - 1; x++) {
      dst[row + x] = (src[row + x - 1] + 2 * src[row + x] + src[row + x + 1]) / 4;
    }
    dst[row + NX - 1] = (src[row + NX - 2] + src[row + NX - 1] * 3) / 4;
  }
}

function blurV(src, dst) {
  for (let x = 0; x < NX; x++) {
    dst[x] = (src[x] * 3 + src[x + NX]) / 4;
    for (let y = 1; y < NY - 1; y++) {
      const n = y * NX + x;
      dst[n] = (src[n - NX] + 2 * src[n] + src[n + NX]) / 4;
    }
    const last = (NY - 1) * NX + x;
    dst[last] = (src[last - NX] + src[last] * 3) / 4;
  }
}

// ---------------------------------------------------------
// drawOutline(ctx, cx, cy, radius)
// ---------------------------------------------------------
// Draws a clean anti-aliased analytic circle outline over the
// canvas after putImageData. Since putImageData replaces pixels
// including anti-aliasing, the outline is drawn as a 2D canvas
// path on top — always a mathematically perfect circle
// regardless of how coarsely the solid mask rasterizes it,
// and drawn on the CSS-scaled canvas so it stays sharp at any
// display zoom level.
function drawOutline(ctx, cx, cy, radius) {
  ctx.save();
  ctx.beginPath();
  ctx.arc(cx + 0.5, cy + 0.5, radius, 0, 2 * Math.PI);
  ctx.strokeStyle = "rgba(255, 255, 255, 0.75)";
  ctx.lineWidth   = 1.0;
  ctx.stroke();
  ctx.restore();
}