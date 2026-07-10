// src/collisionFused.js
// ============================================================
// Single-pass fused collision, modeled directly on the aero-
// tunnel reference's collide(): macros, non-equilibrium stress,
// analytic Smagorinsky-consistent tauEff, and BGK relaxation are
// all computed in ONE loop over the grid, instead of the 6
// separate full-array passes the MRT/LES pipeline used
// (computeMacros, computeEquilibriumAndNonEquilibrium,
// computeStrainRate, computeSmagorinskyViscosity, applyWallModel,
// collideMRT). This is the main lever for matching the
// reference's FPS — most of the previous cost was from repeated
// full-grid sweeps and the 9x9 MRT moment transform, not from
// BGK vs MRT itself.
//
// This intentionally drops: the wall-model log-law correction
// (no longer computed — the reference doesn't have one either,
// per its own paper: "no wall model... the grid is far from LES
// resolution requirements"), and the standing fEq/fNeq/strainMag/
// nuT arrays (no longer needed since nothing outside this
// function consumes them anymore). If you still want to *view*
// strain rate or eddy viscosity in the render dropdown, see the
// note at the bottom for how to expose them cheaply.
// ============================================================

import { N, Q, CX, CY, W, f, rho, ux, uy, NU0 } from "./lbm.js";

// ---------------------------------------------------------
// Stability constants, taken directly from the reference impl.
// ---------------------------------------------------------
const CS2 = 0.03;          // Smagorinsky constant squared (reference's Cs2)
const U_CAP = 0.28;        // hard velocity-magnitude cap before feq
const RHO_MIN = 0.35;
const RHO_MAX = 3.0;

// Reused scratch buffer for equilibrium values at one node —
// allocated once, not per node/frame.
const feqScratch = new Float32Array(Q);

// ---------------------------------------------------------
// collideFused(solid, u0)
// ---------------------------------------------------------
// Params:
//   solid : Uint8Array(N), skipped.
//   u0    : inlet velocity, lattice units — used as the
//           re-equilibration target if a cell's rho/velocity
//           goes non-finite or out of the stable range (matches
//           reference's instant per-cell re-equilibration inside
//           collide(), rather than a separate post-hoc watchdog
//           pass).
// Side effects: overwrites rho, ux, uy, and f in place for every
// fluid node.
export function collideFused(solid, u0) {
  const tauBase = 3 * NU0 + 0.5;

  for (let n = 0; n < N; n++) {
    if (solid[n]) continue;

    // ---- macros ----
    let r = 0, jx = 0, jy = 0;
    for (let i = 0; i < Q; i++) {
      const fi = f[n * Q + i];
      r += fi;
      jx += fi * CX[i];
      jy += fi * CY[i];
    }

    if (!Number.isFinite(r) || r < RHO_MIN || r > RHO_MAX) {
      reequilibrate(n, u0, 0, 1.0);
      continue;
    }

    let u = jx / r;
    let v = jy / r;

    if (!Number.isFinite(u) || !Number.isFinite(v)) {
      reequilibrate(n, u0, 0, 1.0);
      continue;
    }

    const uu = u * u + v * v;
    if (uu > U_CAP * U_CAP) {
      const scale = U_CAP / Math.sqrt(uu);
      u *= scale;
      v *= scale;
    }

    rho[n] = r;
    ux[n] = u;
    uy[n] = v;

    // ---- equilibrium + non-equilibrium stress, same pass ----
    const vx3 = 3 * u, vy3 = 3 * v;
    const vx2 = u * u, vy2 = v * v;
    const vxvy2 = 2 * u * v;
    const v2 = vx2 + vy2;
    const v215 = 1.5 * v2;

    let Pxx = 0, Pyy = 0, Pxy = 0;

    for (let i = 0; i < Q; i++) {
      const cx = CX[i], cy = CY[i];
      const cu = cx * u + cy * v;
      const eq = W[i] * r * (1 + 3 * cu + 4.5 * cu * cu - v215);
      feqScratch[i] = eq;

      const fneq = f[n * Q + i] - eq;
      Pxx += cx * cx * fneq;
      Pyy += cy * cy * fneq;
      Pxy += cx * cy * fneq;
    }

    const Qbar = Math.sqrt(Pxx * Pxx + Pyy * Pyy + 2 * Pxy * Pxy);

    // ---- analytic Smagorinsky-consistent tau (reference Eq. 7) ----
    const tl = tauBase; // sponge is applied separately in boundary.js
    const tauEff = 0.5 * (tl + Math.sqrt(tl * tl + (18 * CS2 * Qbar) / r));
    const omega = 1 / tauEff;

    // ---- BGK relaxation ----
    for (let i = 0; i < Q; i++) {
      const off = n * Q + i;
      f[off] += omega * (feqScratch[i] - f[off]);
    }
  }
}

// ---------------------------------------------------------
// reequilibrate(n, u, v, r)
// ---------------------------------------------------------
// Instantly resets one node to equilibrium at a known-good
// state. Called inline the moment a bad rho/velocity is
// detected, rather than waiting for a separate watchdog pass —
// same reasoning as the reference's per-cell checks in its own
// collide().
function reequilibrate(n, u, v, r) {
  rho[n] = r;
  ux[n] = u;
  uy[n] = v;
  for (let i = 0; i < Q; i++) {
    const cx = CX[i], cy = CY[i];
    const cu = cx * u + cy * v;
    const uu = u * u + v * v;
    f[n * Q + i] = W[i] * r * (1 + 3 * cu + 4.5 * cu * cu - 1.5 * uu);
  }
}