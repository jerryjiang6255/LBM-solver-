// src/boundary.js  (V2 — fused single-pass edges)
// ============================================================

import { NX, NY, N, Q, f, rho, ux, uy, idx, feqAt } from "./lbm.js";

// ---------------------------------------------------------
// Sponge weights — precomputed once at module load.
// ---------------------------------------------------------

const SPONGE_WIDTH    = 28;
const SPONGE_STRENGTH = 0.6;
const spongeW = new Float32Array(NX);

for (let x = NX - 1 - SPONGE_WIDTH; x < NX; x++) {
  const t = (x - (NX - 1 - SPONGE_WIDTH)) / SPONGE_WIDTH;
  spongeW[x] = SPONGE_STRENGTH * t * t;
}

// ---------------------------------------------------------
// Nozzle inlet config
// ---------------------------------------------------------
let _nozzleActive  = false;
let _nozzleCy      = 64;
let _nozzleHalfD   = 16;
let _uCoflow       = 0.05;  // slow ambient freestream (Ma≈0.005)

export function setNozzleInlet(cy, diameter, uCoflow) {
  _nozzleActive = true;
  _nozzleCy     = cy;
  _nozzleHalfD  = diameter / 2;
  _uCoflow      = uCoflow;
}

export function clearNozzleInlet() {
  _nozzleActive = false;
}
// ---------------------------------------------------------
// applyBoundaryConditions(u0, solid)
// ---------------------------------------------------------

export function applyBoundaryConditions(u0, solid) {

  // ---- left column: split inlet for nozzle, uniform for everything else ----
  for (let y = 0; y < NY; y++) {
    const n = idx(0, y);
    if (solid[n]) continue; // nozzle wall nodes — skip

    let inletU;
    if (_nozzleActive) {
      const distFromCenter = Math.abs(y - _nozzleCy);
      inletU = distFromCenter <= _nozzleHalfD ? u0 : _uCoflow;
    } else {
      inletU = u0;
    }

    for (let i = 0; i < Q; i++) {
      f[n * Q + i] = feqAt(i, 1.0, inletU, 0.0);
    }
  }

  // ---- right columns: sponge ----
  const xSpongeStart = NX - 1 - SPONGE_WIDTH;
  for (let x = xSpongeStart; x < NX - 1; x++) {
    const s = spongeW[x];
    for (let y = 0; y < NY; y++) {
      const n = idx(x, y);
      if (solid[n]) continue;
      const r = rho[n], u = ux[n], v = uy[n];
      for (let i = 0; i < Q; i++) {
        const off = n * Q + i;
        f[off] = (1 - s) * f[off] + s * feqAt(i, r, u, v);
      }
    }
  }

  // ---- outlet zero-gradient copy ----
  for (let y = 0; y < NY; y++) {
    const nSrc = idx(NX - 2, y);
    const nDst = idx(NX - 1, y);
    for (let i = 0; i < Q; i++) {
      f[nDst * Q + i] = f[nSrc * Q + i];
    }
  }

  // ---- top and bottom: coflow or freestream depending on mode ----
  for (let x = 0; x < NX; x++) {
    const nb  = idx(x, 0);
    const nt  = idx(x, NY - 1);
    const wallU = _nozzleActive ? _uCoflow : u0;
    for (let i = 0; i < Q; i++) {
      const eq = feqAt(i, 1.0, wallU, 0.0);
      f[nb * Q + i] = eq;
      f[nt * Q + i] = eq;
    }
  }
}