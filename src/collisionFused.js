// src/collisionFused.js  (compare-mode ready, minimal refactor)
// ============================================================

import { CX, CY, W } from "./lbm.js";

// ---------------------------------------------------------
// LES constants
// ---------------------------------------------------------
const CS          = 0.17;
const DELTA       = 1.0;
const CS_DELTA_SQ = CS * CS * DELTA * DELTA;

// ---------------------------------------------------------
// Wall model constants
// ---------------------------------------------------------
const KAPPA      = 0.41;
const B_LOG      = 5.2;
const A_PLUS     = 25.0;
const WALL_ITERS = 6;
const WALL_TOL   = 1e-4;
const Y_MIN      = 0.5;

// ---------------------------------------------------------
// Stability caps
// ---------------------------------------------------------
const U_CAP   = 0.28;
const RHO_MIN = 0.35;
const RHO_MAX = 3.0;
const TAU_MIN = 0.55;
const TAU_MAX = 5.0;

// ---------------------------------------------------------
// TRT Magic Parameter
// ---------------------------------------------------------
const LAMBDA = 3.0 / 16.0;
const TRT_PAIRS = new Uint8Array([1, 3,  2, 4,  5, 7,  6, 8]);

// feq scratch remains module-local since both sims run sequentially,
// not concurrently. Q is fixed at 9 in current architecture.
const feqScratch = new Float32Array(9);

// ---------------------------------------------------------
// collideFused(sim, u0)
// ---------------------------------------------------------
export function collideFused(sim, u0) {
  const { N, Q, f, rho, ux, uy, solid, wallAdjacent, wallDist, NU0 } = sim;

  const STRAIN_PRE = -1.5; // -3/2 factor in Sij = -3/(2*rho)*sum(...)

  for (let n = 0; n < N; n++) {
    if (solid[n]) continue;

    const base = n * Q;

    // =====================================================
    // Direction loop 1: macroscopic moments
    // =====================================================
    let r = 0.0, jx = 0.0, jy = 0.0;
    for (let i = 0; i < Q; i++) {
      const fi = f[base + i];
      r  += fi;
      jx += fi * CX[i];
      jy += fi * CY[i];
    }

    if (!Number.isFinite(r) || r < RHO_MIN || r > RHO_MAX) {
      reequilibrate(sim, n, u0, 0.0, 1.0);
      continue;
    }

    let u = jx / r;
    let v = jy / r;

    if (!Number.isFinite(u) || !Number.isFinite(v)) {
      reequilibrate(sim, n, u0, 0.0, 1.0);
      continue;
    }

    const spd2 = u * u + v * v;
    if (spd2 > U_CAP * U_CAP) {
      const sc = U_CAP / Math.sqrt(spd2);
      u *= sc;
      v *= sc;
    }

    rho[n] = r;
    ux[n]  = u;
    uy[n]  = v;

    // =====================================================
    // Direction loop 2: feq + non-equilibrium stress
    // =====================================================
    const v215         = 1.5 * (u * u + v * v);
    const strainFactor = STRAIN_PRE / r; // -3/(2*rho)

    let Pxx = 0.0, Pyy = 0.0, Pxy = 0.0;

    for (let i = 0; i < Q; i++) {
      const cx = CX[i];
      const cy = CY[i];
      const cu = cx * u + cy * v;
      const eq = W[i] * r * (1.0 + 3.0 * cu + 4.5 * cu * cu - v215);
      feqScratch[i] = eq;

      const fneq = f[base + i] - eq;
      Pxx += cx * cx * fneq;
      Pyy += cy * cy * fneq;
      Pxy += cx * cy * fneq;
    }

    // =====================================================
    // Scalar: full LES strain rate + Smagorinsky nu_t
    // =====================================================
    const Sxx  = strainFactor * Pxx;
    const Syy  = strainFactor * Pyy;
    const Sxy  = strainFactor * Pxy;
    const Smag = Math.sqrt(2.0 * (Sxx * Sxx + Syy * Syy + 2.0 * Sxy * Sxy));
    let nuT    = CS_DELTA_SQ * Smag;

    // =====================================================
    // Scalar: wall model (log-law + Van Driest)
    // =====================================================
    if (wallAdjacent[n]) {
      const y    = Math.max(wallDist[n], Y_MIN);
      const uMag = Math.sqrt(u * u + v * v);

      let uStar = Math.max(0.1 * uMag, 1e-6);
      for (let k = 0; k < WALL_ITERS; k++) {
        const yuOvNu = Math.max((y * uStar) / NU0, 1e-3);
        const lhs    = (1.0 / KAPPA) * Math.log(yuOvNu) + B_LOG;
        const R      = uStar * lhs - uMag;
        const Rp     = lhs + 1.0 / KAPPA;
        if (Math.abs(Rp) < 1e-10) break;

        let next = uStar - R / Rp;
        if (!Number.isFinite(next) || next <= 0.0) next = 1e-6;
        if (Math.abs(next - uStar) < WALL_TOL) {
          uStar = next;
          break;
        }
        uStar = next;
      }

      const yPlus = (y * uStar) / NU0;
      const fDamp = 1.0 - Math.exp(-yPlus / A_PLUS);
      nuT = nuT * fDamp;
    }

    // =====================================================
    // Scalar: TRT relaxation rates
    // =====================================================
    const nuTot = NU0 + nuT;
    let tauPlus = 3.0 * nuTot + 0.5;
    if (tauPlus < TAU_MIN) tauPlus = TAU_MIN;
    if (tauPlus > TAU_MAX) tauPlus = TAU_MAX;

    const omegaPlus  = 1.0 / tauPlus;
    const tauMinus   = LAMBDA / (tauPlus - 0.5) + 0.5;
    const omegaMinus = 1.0 / tauMinus;

    // =====================================================
    // Direction loop 3: TRT relax + recombine
    // =====================================================
    f[base] -= omegaPlus * (f[base] - feqScratch[0]);

    for (let p = 0; p < 8; p += 2) {
      const i  = TRT_PAIRS[p];
      const ib = TRT_PAIRS[p + 1];

      const fi   = f[base + i];
      const fib  = f[base + ib];
      const eqi  = feqScratch[i];
      const eqib = feqScratch[ib];

      const fPlus   = 0.5 * (fi + fib);
      const fMinus  = 0.5 * (fi - fib);
      const eqPlus  = 0.5 * (eqi + eqib);
      const eqMinus = 0.5 * (eqi - eqib);

      const fPlusPost  = fPlus  - omegaPlus  * (fPlus  - eqPlus);
      const fMinusPost = fMinus - omegaMinus * (fMinus - eqMinus);

      f[base + i]  = fPlusPost + fMinusPost;
      f[base + ib] = fPlusPost - fMinusPost;
    }
  }
}

// ---------------------------------------------------------
// reequilibrate(sim, n, u, v, r)
// ---------------------------------------------------------
function reequilibrate(sim, n, u, v, r) {
  const { Q, f, rho, ux, uy } = sim;

  rho[n] = r;
  ux[n]  = u;
  uy[n]  = v;

  const base = n * Q;
  const v215 = 1.5 * (u * u + v * v);

  for (let i = 0; i < Q; i++) {
    const cu = CX[i] * u + CY[i] * v;
    f[base + i] = W[i] * r * (1.0 + 3.0 * cu + 4.5 * cu * cu - v215);
  }
}



