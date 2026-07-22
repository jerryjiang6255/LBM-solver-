// src/streaming.js
// ============================================================
// Compare-mode ready streaming.
// Updated for global fidelity presets: no fixed NX/NY imports.
// Added: standard Ladd-style momentum exchange accumulation
// for fixed SDF/BFL body links using the existing wall-link loop.
// Added: drag/lift/moment coefficient calculation.
// ============================================================

import { CX, CY, OPP } from "./lbm.js";

// ---------------------------------------------------------
// stream(sim)
// ---------------------------------------------------------
export function stream(sim) {
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
function applyBFL(sim) {
  const { NX, NY, Q, f, fSnapshot, linkMask, linkQ, fluidSurfaceList } = sim;

  ensureForceStorage(sim);

  if (!fluidSurfaceList || !linkMask || !linkQ) {
    sim.force.fx = 0.0;
    sim.force.fy = 0.0;
    sim.force.mz = 0.0;

    sim.coeff.cd = 0.0;
    sim.coeff.cl = 0.0;
    sim.coeff.cm = 0.0;
    return;
  }

  fSnapshot.set(f);

  // -------------------------------------------------------
  // Standard Ladd-style force accumulation
  // Uses existing active wall-link loop (no extra pass).
  // Reference point for moment comes from sim.forceRef if present.
  // -------------------------------------------------------
  let fx = 0.0;
  let fy = 0.0;
  let mz = 0.0;

  const xRef = sim.forceRef?.x ?? 0.0;
  const yRef = sim.forceRef?.y ?? 0.0;

  for (let k = 0; k < fluidSurfaceList.length; k++) {
    const n = fluidSurfaceList[k];
    const x = n % NX;
    const y = (n / NX) | 0;
    const baseN = n * Q;

    for (let i = 1; i < Q; i++) {
      const off = baseN + i;
      if (!linkMask[off]) continue;

      const q = linkQ[off];
      const opp = OPP[i];
      const fOppHere = fSnapshot[baseN + opp];

      // ---------------------------------------------------
      // Apply BFL boundary update
      // ---------------------------------------------------
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

      // ---------------------------------------------------
      // Standard Ladd-style momentum exchange approximation
      // ---------------------------------------------------
      const fRet = f[off];

      const dFx = (fRet + fOppHere) * CX[i];
      const dFy = (fRet + fOppHere) * CY[i];

      fx += dFx;
      fy += dFy;

      const xw = x + q * CX[i];
      const yw = y + q * CY[i];

      mz += (xw - xRef) * dFy - (yw - yRef) * dFx;
    }
  }

  sim.force.fx = fx;
  sim.force.fy = fy;
  sim.force.mz = mz;

  // -------------------------------------------------------
  // Coefficients
  //
  // cd = Fx / forceScale
  // cl = Fy / forceScale
  // cm = Mz / momentScale
  //
  // forceScale and momentScale should be set elsewhere
  // (e.g. in index.html when geometry is built/applied).
  // -------------------------------------------------------
  const forceScale = sim.forceScale ?? 1.0;
  const momentScale = sim.momentScale ?? forceScale * Math.max(sim.L || 1.0, 1.0);

  sim.coeff.cd = forceScale !== 0 ? fx / forceScale : 0.0;
  sim.coeff.cl = forceScale !== 0 ? fy / forceScale : 0.0;
  sim.coeff.cm = momentScale !== 0 ? mz / momentScale : 0.0;
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

// ---------------------------------------------------------
// ensureForceStorage(sim)
// ---------------------------------------------------------
function ensureForceStorage(sim) {
  if (!sim.force) {
    sim.force = {
      fx: 0.0,
      fy: 0.0,
      mz: 0.0,
    };
  }

  if (!sim.coeff) {
    sim.coeff = {
      cd: 0.0,
      cl: 0.0,
      cm: 0.0,
    };
  }
}


