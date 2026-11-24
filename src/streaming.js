// src/streaming.js
// ============================================================
// Compare-mode ready streaming.
// Updated for global fidelity presets: no fixed NX/NY imports.
//
// This version keeps streaming + wall reconstruction only.
// Force / moment / coefficient calculation is no longer done
// here; it is handled separately by forceSurface.js.
// ============================================================

import { CX, CY, OPP } from "./lbm.js";

// ---------------------------------------------------------
// stream(sim)
// ---------------------------------------------------------
export function stream(sim) {
  // Keep the post-collision field: BFL and the force both need f*,
  // and streamInPlace overwrites it (the GPU reads the same values from fB).
  if (sim.fSnapshot) sim.fSnapshot.set(sim.f);
  streamInPlace(sim);
  applyBFL(sim);
  bounceBackNodes(sim);
}

// ---------------------------------------------------------
// streamInPlace(sim)
// ---------------------------------------------------------
function streamInPlace(sim) {
  const { NX, NY, Q, f } = sim;

  // Pass 1: N(2) + NW(6), sweep top→bottom, left→right
  for (let y = NY - 2; y >= 1; y--) {
    const rowDst = y * NX;
    const rowSrc = (y - 1) * NX;
    for (let x = 1; x < NX - 1; x++) {
      f[(rowDst + x) * Q + 2] = f[(rowSrc + x    ) * Q + 2];
      f[(rowDst + x) * Q + 6] = f[(rowSrc + x + 1) * Q + 6];
    }
  }

  // Pass 2: E(1) + NE(5), sweep top→bottom, right→left
  for (let y = NY - 2; y >= 1; y--) {
    const rowDst = y * NX;
    const rowSrc = (y - 1) * NX;
    for (let x = NX - 2; x >= 1; x--) {
      f[(rowDst + x) * Q + 1] = f[(rowDst + x - 1) * Q + 1];
      f[(rowDst + x) * Q + 5] = f[(rowSrc + x - 1) * Q + 5];
    }
  }

  // Pass 3: S(4) + SE(8), sweep bottom→top, right→left
  for (let y = 1; y < NY - 1; y++) {
    const rowDst = y * NX;
    const rowSrc = (y + 1) * NX;
    for (let x = NX - 2; x >= 1; x--) {
      f[(rowDst + x) * Q + 4] = f[(rowSrc + x    ) * Q + 4];
      f[(rowDst + x) * Q + 8] = f[(rowSrc + x - 1) * Q + 8];
    }
  }

  // Pass 4: W(3) + SW(7), sweep bottom→top, left→right
  for (let y = 1; y < NY - 1; y++) {
    const rowDst = y * NX;
    const rowSrc = (y + 1) * NX;
    for (let x = 1; x < NX - 1; x++) {
      f[(rowDst + x) * Q + 3] = f[(rowDst + x + 1) * Q + 3];
      f[(rowDst + x) * Q + 7] = f[(rowSrc + x + 1) * Q + 7];
    }
  }
}

// ---------------------------------------------------------
// applyBFL(sim)
// ---------------------------------------------------------
// BFL (Bouzidi) interpolated bounce-back, same logic as the GPU
// stream_main. Link i at a fluid node x means x + c_i is solid and the
// wall sits a fraction q along c_i. After streaming, slot i at x is
// known (it came from the fluid side) and is left alone; the unknown is
// the population coming back out of the wall, o = OPP[i]. It is rebuilt
// from the post-collision field f* (fSnapshot, taken in stream()):
//   q <  1/2: f_o(x) = 2q f*_i(x) + (1 - 2q) f*_i(x - c_i)
//   q >= 1/2: f_o(x) = f*_i(x) / (2q) + (1 - 1/(2q)) f*_o(x)
// If x - c_i is outside the domain or solid, f*_i(x) is used instead.
// No force accumulation is performed here.
// ---------------------------------------------------------
function applyBFL(sim) {
  const { NX, NY, Q, f, fSnapshot, linkMask, linkQ, fluidSurfaceList, solid } = sim;

  if (!fluidSurfaceList || !linkMask || !linkQ || !fSnapshot) return;

  for (let k = 0; k < fluidSurfaceList.length; k++) {
    const n = fluidSurfaceList[k];
    const x = n % NX;
    const y = (n / NX) | 0;
    const baseN = n * Q;

    for (let i = 1; i < Q; i++) {
      const off = baseN + i;
      if (!linkMask[off]) continue;

      const q = linkQ[off];
      const o = OPP[i];
      const fiHere = fSnapshot[off];

      if (q < 0.5) {
        const xb = x - CX[i];
        const yb = y - CY[i];
        let fiBehind = fiHere;
        if (xb >= 0 && xb < NX && yb >= 0 && yb < NY) {
          const nb = yb * NX + xb;
          if (!solid[nb]) fiBehind = fSnapshot[nb * Q + i];
        }
        f[baseN + o] = 2.0 * q * fiHere + (1.0 - 2.0 * q) * fiBehind;
      } else {
        const inv2q = 1.0 / (2.0 * q);
        f[baseN + o] = inv2q * fiHere + (1.0 - inv2q) * fSnapshot[baseN + o];
      }
    }
  }
}

// ---------------------------------------------------------
// bounceBackNodes(sim)
// ---------------------------------------------------------
function bounceBackNodes(sim) {
  const { Q, f, paintedNodes } = sim;

  for (const n of paintedNodes) {
    const base = n * Q;
    let tmp;

    tmp = f[base + 1]; f[base + 1] = f[base + 3]; f[base + 3] = tmp;
    tmp = f[base + 2]; f[base + 2] = f[base + 4]; f[base + 4] = tmp;
    tmp = f[base + 5]; f[base + 5] = f[base + 7]; f[base + 7] = tmp;
    tmp = f[base + 6]; f[base + 6] = f[base + 8]; f[base + 8] = tmp;
  }
}

