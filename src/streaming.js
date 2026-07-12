// src/streaming.js  (V3 — SDF BFL reconstruction)
// ============================================================

import { NX, NY, N, Q, CX, CY, OPP, f } from "./lbm.js";
const fSnapshot = new Float32Array(Q * N);

// ---------------------------------------------------------
// stream(solid, solidList, linkMask, linkQ)
// ---------------------------------------------------------
export function stream(solid, solidList, linkMask, linkQ, fluidSurfaceList, paintedNodes) {
  streamInPlace();
  applyBFL(linkMask, linkQ, fluidSurfaceList);
  bounceBackNodes(paintedNodes);
}
// ---------------------------------------------------------
// streamInPlace()
// ---------------------------------------------------------
function streamInPlace() {
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
// applyBFL(solid, linkMask, linkQ)
// ---------------------------------------------------------
function applyBFL(linkMask, linkQ, fluidSurfaceList) {
  fSnapshot.set(f);

  for (let k = 0; k < fluidSurfaceList.length; k++) {
    const n = fluidSurfaceList[k];
    const x = n % NX;
    const y = (n / NX) | 0;
    const baseN = n * Q;

    for (let i = 1; i < Q; i++) {
      const off = baseN + i;
      if (!linkMask[off]) continue;

      const q   = linkQ[off];
      const opp = OPP[i];
      const fOppHere = fSnapshot[baseN + opp];

      if (Math.abs(q - 0.5) < 1e-6) {
        f[off] = fOppHere;

      } else if (q < 0.5) {
        const xb = x - CX[i];
        const yb = y - CY[i];
        const fOppBehind = (xb >= 0 && xb < NX && yb >= 0 && yb < NY)
          ? fSnapshot[(yb * NX + xb) * Q + opp]
          : fOppHere;
        const inv1p2q = 1.0 / (1.0 + 2.0 * q);
        f[off] = inv1p2q * fOppHere + (2.0 * q * inv1p2q) * fOppBehind;

      } else {
        const xa = x + CX[i];
        const ya = y + CY[i];
        const fOppAhead = (xa >= 0 && xa < NX && ya >= 0 && ya < NY)
          ? fSnapshot[(ya * NX + xa) * Q + opp]
          : fOppHere;
        const inv2q = 1.0 / (2.0 * q);
        f[off] = inv2q * fOppHere + (1.0 - inv2q) * fOppAhead;
      }
    }
  }
}

function bounceBackNodes(paintedNodes) {
  for (const n of paintedNodes) {
    const base = n * Q;
    let tmp;
    tmp = f[base+1]; f[base+1] = f[base+3]; f[base+3] = tmp;
    tmp = f[base+2]; f[base+2] = f[base+4]; f[base+4] = tmp;
    tmp = f[base+5]; f[base+5] = f[base+7]; f[base+7] = tmp;
    tmp = f[base+6]; f[base+6] = f[base+8]; f[base+8] = tmp;
  }
}