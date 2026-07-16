// src/boundary.js
// ============================================================
// Compare-mode ready boundary conditions.
// Updated for global fidelity presets: no fixed NX/NY imports.
// Change 1: supports pulsing inlet velocity via sim.flow.
// ============================================================

import { feqAt, getInletVelocity } from "./simInstance.js";

// ---------------------------------------------------------
// Sponge config
// ---------------------------------------------------------
const SPONGE_WIDTH    = 28;
const SPONGE_STRENGTH = 0.6;

// Cache sponge weights by NX so they are built once per grid width
const spongeCache = new Map();

// ---------------------------------------------------------
// getSpongeWeights(NX)
// ---------------------------------------------------------
function getSpongeWeights(NX) {
  let spongeW = spongeCache.get(NX);
  if (spongeW) return spongeW;

  spongeW = new Float32Array(NX);

  const width = Math.min(SPONGE_WIDTH, Math.max(1, NX - 2));
  const start = NX - 1 - width;

  for (let x = start; x < NX; x++) {
    const t = (x - start) / width;
    spongeW[x] = SPONGE_STRENGTH * t * t;
  }

  spongeCache.set(NX, spongeW);
  return spongeW;
}

// ---------------------------------------------------------
// setNozzleInlet(sim, cy, diameter, uCoflow)
// ---------------------------------------------------------
export function setNozzleInlet(sim, cy, diameter, uCoflow) {
  sim.nozzle.active  = true;
  sim.nozzle.cy      = cy;
  sim.nozzle.halfD   = diameter / 2;
  sim.nozzle.uCoflow = uCoflow;
}

// ---------------------------------------------------------
// clearNozzleInlet(sim)
// ---------------------------------------------------------
export function clearNozzleInlet(sim) {
  sim.nozzle.active = false;
}

// ---------------------------------------------------------
// applyBoundaryConditions(sim, u0)
// ---------------------------------------------------------
// u0 is kept for backward compatibility, but the actual inlet
// velocity comes from getInletVelocity(sim) when pulsation is on.
export function applyBoundaryConditions(sim, u0) {
  const { NX, NY, Q, f, rho, ux, uy, solid, nozzle } = sim;
  const idx = (x, y) => y * NX + x;
  const spongeW = getSpongeWeights(NX);

  const inletBaseU = getInletVelocity(sim);

  // ---- left column: split inlet for nozzle, uniform for everything else ----
  for (let y = 0; y < NY; y++) {
    const n = idx(0, y);
    if (solid[n]) continue;

    let inletU;
    if (nozzle.active) {
      const distFromCenter = Math.abs(y - nozzle.cy);
      inletU = distFromCenter <= nozzle.halfD ? inletBaseU : nozzle.uCoflow;
    } else {
      inletU = inletBaseU;
    }

    const base = n * Q;
    for (let i = 0; i < Q; i++) {
      f[base + i] = feqAt(i, 1.0, inletU, 0.0);
    }
  }

  // ---- right columns: sponge ----
  const width = Math.min(SPONGE_WIDTH, Math.max(1, NX - 2));
  const xSpongeStart = NX - 1 - width;

  for (let x = xSpongeStart; x < NX - 1; x++) {
    const s = spongeW[x];

    for (let y = 0; y < NY; y++) {
      const n = idx(x, y);
      if (solid[n]) continue;

      const r = rho[n];
      const u = ux[n];
      const v = uy[n];
      const base = n * Q;

      for (let i = 0; i < Q; i++) {
        const off = base + i;
        f[off] = (1 - s) * f[off] + s * feqAt(i, r, u, v);
      }
    }
  }

  // ---- outlet zero-gradient copy ----
  if (NX >= 2) {
    for (let y = 0; y < NY; y++) {
      const nSrc = idx(NX - 2, y);
      const nDst = idx(NX - 1, y);

      for (let i = 0; i < Q; i++) {
        f[nDst * Q + i] = f[nSrc * Q + i];
      }
    }
  }

  // ---- top and bottom: coflow or freestream depending on mode ----
  const wallU = nozzle.active ? nozzle.uCoflow : inletBaseU;

  for (let x = 0; x < NX; x++) {
    const nb = idx(x, 0);
    const nt = idx(x, NY - 1);

    const baseB = nb * Q;
    const baseT = nt * Q;

    for (let i = 0; i < Q; i++) {
      const eq = feqAt(i, 1.0, wallU, 0.0);
      f[baseB + i] = eq;
      f[baseT + i] = eq;
    }
  }
}


