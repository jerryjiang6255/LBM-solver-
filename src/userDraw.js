// src/userDraw.js
// ============================================================

import { NX, NY, N, Q, CX, CY, W, f, rho, ux, uy } from "./lbm.js";
const BRUSH_R = 1; 
export const paintedNodes = new Set();
export const strokes = [];        
let currentStroke = null;        
let lastLX = null;
let lastLY = null;
let wallCanvasCtx = null;

export function setWallCanvas(canvas) {
  wallCanvasCtx = canvas;
}
export function beginStroke(lx, ly) {
  lastLX = Math.round(lx);
  lastLY = Math.round(ly);
  currentStroke = [{ x: lx, y: ly }];  
}

export function paintStroke(lx, ly, solid, mode) {
  const x1 = Math.round(lx);
  const y1 = Math.round(ly);

  if (lastLX === null) {
    paintBrush(x1, y1, solid, mode);
  } else {
    bresenham(lastLX, lastLY, x1, y1, solid, mode);
  }

  lastLX = x1;
  lastLY = y1;

  if (currentStroke && mode === 'draw') {
    currentStroke.push({ x: lx, y: ly });
  }
}

export function endStroke(mode) {
  lastLX = null;
  lastLY = null;
  if (currentStroke && currentStroke.length > 0 && mode === 'draw') {
    strokes.push(currentStroke);  // commit completed stroke
  }
  currentStroke = null;
}

function bresenham(x0, y0, x1, y1, solid, mode) {
  const dx =  Math.abs(x1 - x0);
  const dy = -Math.abs(y1 - y0);
  const sx = x0 < x1 ? 1 : -1;
  const sy = y0 < y1 ? 1 : -1;
  let err = dx + dy;

  let x = x0, y = y0;
  while (true) {
    paintBrush(x, y, solid, mode);
    if (x === x1 && y === y1) break;
    const e2 = 2 * err;
    if (e2 >= dy) { err += dy; x += sx; }
    if (e2 <= dx) { err += dx; y += sy; }
  }
}

// ---------------------------------------------------------
// paintBrush(lx, ly, solid, mode)
// ---------------------------------------------------------
export function paintBrush(lx, ly, solid, mode) {
  const x0 = Math.round(lx);
  const y0 = Math.round(ly);

  for (let dy = -BRUSH_R; dy <= BRUSH_R; dy++) {
    for (let dx = -BRUSH_R; dx <= BRUSH_R; dx++) {
      if (dx * dx + dy * dy > BRUSH_R * BRUSH_R) continue;
      const x = x0 + dx;
      const y = y0 + dy;
      if (x < 1 || x >= NX - 1 || y < 1 || y >= NY - 1) continue;
      const n = y * NX + x;
      if (mode === 'draw') {
        if (solid[n]) continue;         
        solid[n] = 1;
        paintedNodes.add(n);
        _reequilibrate(n);

      } else {                          // 'erase'
        if (!paintedNodes.has(n)) continue; // never erase SDF geometry
        solid[n] = 0;
        paintedNodes.delete(n);
        _seedFromNeighbors(n, x, y, solid);
      }
    }
  }
}

// ---------------------------------------------------------
// clearDrawing(solid, solidBase)
// ---------------------------------------------------------
export function clearDrawing(solid, solidBase) {
  for (const n of paintedNodes) {
    solid[n] = solidBase[n];
  }
  paintedNodes.clear();
  strokes.length = 0;       
  currentStroke  = null;
}

// ---------------------------------------------------------
// _reequilibrate(n)
// ---------------------------------------------------------
function _reequilibrate(n) {
  rho[n] = 1.0;
  ux[n]  = 0;
  uy[n]  = 0;
  const base = n * Q;
  for (let i = 0; i < Q; i++) {
    f[base + i] = W[i];           
  }
}

// ---------------------------------------------------------
// _seedFromNeighbors(n, x, y, solid)
// ---------------------------------------------------------
function _seedFromNeighbors(n, x, y, solid) {
  let rhoSum = 0, uxSum = 0, uySum = 0, count = 0;
  for (let i = 1; i < Q; i++) {
    const nx = x + CX[i];
    const ny = y + CY[i];
    if (nx < 0 || nx >= NX || ny < 0 || ny >= NY) continue;
    const nb = ny * NX + nx;
    if (solid[nb]) continue;
    rhoSum += rho[nb];
    uxSum  += ux[nb];
    uySum  += uy[nb];
    count++;
  }
  const r = count > 0 ? rhoSum / count : 1.0;
  const u = count > 0 ? uxSum  / count : 0;
  const v = count > 0 ? uySum  / count : 0;
  rho[n] = r; ux[n] = u; uy[n] = v;

  const v215 = 1.5 * (u * u + v * v);
  const base = n * Q;
  for (let i = 0; i < Q; i++) {
    const cu = CX[i] * u + CY[i] * v;
    f[base + i] = W[i] * r * (1 + 3*cu + 4.5*cu*cu - v215);  // was: f[i*N+n]
  }
}

export function redrawWalls() {
  if (!wallCanvasCtx) return;
  const cnv  = wallCanvasCtx;           
  const ctx  = cnv.getContext('2d');
  const cW   = cnv.width;              
  const cH   = cnv.height;
  const sx   = cW / NX;
  const sy   = cH / NY;

  ctx.clearRect(0, 0, cW, cH);

  const allStrokes = currentStroke
    ? [...strokes, currentStroke]
    : strokes;

  if (allStrokes.length === 0) return;

  const lineW = sx * (BRUSH_R + 0.5);

  for (const pts of allStrokes) {
    if (pts.length === 0) continue;

    // Pass 1 — glow: widest, semi-transparent cyan, drawn first (bottom)
    ctx.save();
    ctx.strokeStyle = 'rgba(175, 198, 230, 0.75)';
    ctx.lineWidth   = lineW * 1.25;
    ctx.lineCap     = 'round';
    ctx.lineJoin    = 'round';
    _drawCatmullRom(ctx, pts, sx, sy);
    ctx.stroke();
    ctx.restore();

    // Pass 2 — body: medium width, dark grey fill (middle)
    ctx.save();
    ctx.strokeStyle = '#292a2f';   // dark grey matching your current color
    ctx.lineWidth   = lineW;
    ctx.lineCap     = 'round';
    ctx.lineJoin    = 'round';
    _drawCatmullRom(ctx, pts, sx, sy);
    ctx.stroke();
    ctx.restore();
   }
}

// Draws a Catmull-Rom spline through pts array of {x,y} lattice coords.
function _drawCatmullRom(ctx, pts, sx, sy) {
  if (pts.length === 1) {
    // Single point — draw a circle
    const x = pts[0].x * sx;
    const y = pts[0].y * sy;
    ctx.beginPath();
    ctx.arc(x, y, ctx.lineWidth / 2, 0, Math.PI * 2);
    ctx.fill();
    return;
  }

  ctx.beginPath();
  ctx.moveTo(pts[0].x * sx, pts[0].y * sy);

  for (let i = 1; i < pts.length - 1; i++) {
    const mx = (pts[i].x + pts[i + 1].x) / 2 * sx;
    const my = (pts[i].y + pts[i + 1].y) / 2 * sy;
    ctx.quadraticCurveTo(
      pts[i].x * sx,
      pts[i].y * sy,
      mx, my
    );
  }

  const last = pts[pts.length - 1];
  ctx.lineTo(last.x * sx, last.y * sy);
}
