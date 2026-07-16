// src/collisionFused.js  (compare-mode ready, V6 physics-ready)
// ============================================================

import { CX, CY, W } from "./lbm.js";

// ---------------------------------------------------------
// LES constants
// ---------------------------------------------------------
const CS          = 0.17;   // Smagorinsky constant
const CW          = 0.325;  // WALE constant
const DELTA       = 1.0;
const CS_DELTA_SQ = CS * CS * DELTA * DELTA;
const CW_DELTA_SQ = CW * CW * DELTA * DELTA;
const WALE_EPS    = 1e-20;

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

const feqScratch = new Float32Array(9);

// ---------------------------------------------------------
// collideFused(sim, u0)
// ---------------------------------------------------------
// u0 should be the instantaneous inlet velocity if pulsation is active.
export function collideFused(sim, u0) {
  const { N, Q, NX, NY, f, rho, ux, uy, solid, wallAdjacent, wallDist, NU0 } = sim;
  const physics = sim.physics || {
    collisionModel: "trt",
    lesModel: "smagorinsky",
    wallModel: true,
  };

  const STRAIN_PRE = -1.5;

  for (let n = 0; n < N; n++) {
    if (solid[n]) continue;

    const base = n * Q;

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

    const v215         = 1.5 * (u * u + v * v);
    const strainFactor = -1.5 / r;

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

    let nuT = computeSubgridViscosity(sim, physics, n, u, v, strainFactor, Pxx, Pyy, Pxy, NX, NY);
    nuT = applyWallModelIfEnabled(physics, wallAdjacent, wallDist, n, u, v, NU0, nuT);

    const nuTot = NU0 + nuT;

    if (physics.collisionModel === "bgk") {
      relaxBGK(base, Q, f, feqScratch, nuTot);
    } else {
      relaxTRT(base, f, feqScratch, nuTot);
    }
  }
}

function computeSubgridViscosity(sim, physics, n, u, v, strainFactor, Pxx, Pyy, Pxy, NX, NY) {
  if (physics.lesModel === "off") return 0.0;

  if (physics.lesModel === "wale") {
    return computeWaleNuT(sim, n, NX, NY);
  }

  const Sxx  = strainFactor * Pxx;
  const Syy  = strainFactor * Pyy;
  const Sxy  = strainFactor * Pxy;
  const Smag = Math.sqrt(2.0 * (Sxx * Sxx + Syy * Syy + 2.0 * Sxy * Sxy));
  return CS_DELTA_SQ * Smag;
}

function computeWaleNuT(sim, n, NX, NY) {
  const { ux, uy, solid } = sim;

  const x = n % NX;
  const y = (n / NX) | 0;

  if (x <= 0 || x >= NX - 1 || y <= 0 || y >= NY - 1) return 0.0;

  const nL = n - 1;
  const nR = n + 1;
  const nD = n - NX;
  const nU = n + NX;

  if (solid[nL] || solid[nR] || solid[nD] || solid[nU]) return 0.0;

  const dudx = 0.5 * (ux[nR] - ux[nL]);
  const dudy = 0.5 * (ux[nU] - ux[nD]);
  const dvdx = 0.5 * (uy[nR] - uy[nL]);
  const dvdy = 0.5 * (uy[nU] - uy[nD]);

  const Sxx = dudx;
  const Syy = dvdy;
  const Sxy = 0.5 * (dudy + dvdx);

  const IS = Sxx * Sxx + Syy * Syy + 2.0 * Sxy * Sxy;

  const Gxx = dudx * dudx + dudy * dvdx;
  const Gxy = dudy * (dudx + dvdy);
  const Gyx = dvdx * (dudx + dvdy);
  const Gyy = dvdx * dudy + dvdy * dvdy;

  const trG = Gxx + Gyy;
  const Sd_xx = Gxx - (1.0 / 3.0) * trG;
  const Sd_yy = Gyy - (1.0 / 3.0) * trG;
  const Sd_xy = 0.5 * (Gxy + Gyx);

  const ID = Sd_xx * Sd_xx + Sd_yy * Sd_yy + 2.0 * Sd_xy * Sd_xy;

  if (ID <= WALE_EPS) return 0.0;

  const num = Math.pow(ID, 1.5);
  const den = Math.pow(IS, 2.5) + Math.pow(ID, 1.25) + WALE_EPS;

  return CW_DELTA_SQ * (num / den);
}

function applyWallModelIfEnabled(physics, wallAdjacent, wallDist, n, u, v, NU0, nuT) {
  if (!physics.wallModel) return nuT;
  if (physics.lesModel === "off") return nuT;
  if (!wallAdjacent[n]) return nuT;

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
  return nuT * fDamp;
}

function relaxTRT(base, f, feqScratch, nuTot) {
  let tauPlus = 3.0 * nuTot + 0.5;
  if (tauPlus < TAU_MIN) tauPlus = TAU_MIN;
  if (tauPlus > TAU_MAX) tauPlus = TAU_MAX;

  const omegaPlus  = 1.0 / tauPlus;
  const tauMinus   = LAMBDA / (tauPlus - 0.5) + 0.5;
  const omegaMinus = 1.0 / tauMinus;

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

function relaxBGK(base, Q, f, feqScratch, nuTot) {
  let tau = 3.0 * nuTot + 0.5;
  if (tau < TAU_MIN) tau = TAU_MIN;
  if (tau > TAU_MAX) tau = TAU_MAX;

  const omega = 1.0 / tau;

  for (let i = 0; i < Q; i++) {
    f[base + i] -= omega * (f[base + i] - feqScratch[i]);
  }
}

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



