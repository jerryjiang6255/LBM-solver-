// src/forceSurface.js
// ============================================================
// Force via momentum exchange (same formula as the GPU force_main).
//
// For each blocked link i at fluid node n (n + c_i solid):
//   dF = c_i * (f*_i(n) + f_o(n) - 2 W[i] rho0),   o = OPP[i]
// f*_i is the post-collision population going into the wall
// (sim.fSnapshot, taken in streaming.js stream() before streaming) and
// f_o is the population BFL sends back (sim.f after stream()).
// The rest state 2 W[i] rho0 sums to zero over a closed body; removing
// it avoids cancellation between the large front and rear terms.
//
// Requires:
//   sim.fSnapshot  — post-collision f, captured in streaming.js stream()
//   sim.f          — post-stream f after BFL
//   sim.RHO0       — background density (default 1.0)
// ============================================================

import { CX, CY, W, OPP } from "./lbm.js";

const ALL_DIRS = [1, 2, 3, 4, 5, 6, 7, 8];

// ---------------------------------------------------------
// Public API
// ---------------------------------------------------------
export function computeSurfaceForces(sim) {
  ensureForceStorage(sim);

  if (!sim.fluidSurfaceList || !sim.linkMask || !sim.fSnapshot) {
    zeroForces(sim);
    return;
  }

  const result = momentumExchange(sim);

  sim.force.fx = result.fx;
  sim.force.fy = result.fy;
  sim.force.mz = result.mz;

  const { forceScale, momentScale } = getScaling(sim);

  // Cl negated: y=0 is at top of domain so upward lift = negative Fy.
  // Negate to match aerodynamic convention (positive Cl = upward lift).
  sim.coeff.cd =  forceScale  !== 0 ?  result.fx / forceScale  : 0.0;
  sim.coeff.cl =  forceScale  !== 0 ? -result.fy / forceScale  : 0.0;
  sim.coeff.cm =  momentScale !== 0 ?  result.mz / momentScale : 0.0;
}

// ---------------------------------------------------------
// Momentum exchange
// ---------------------------------------------------------
function momentumExchange(sim) {
  const { Q, NX, f, fSnapshot, fluidSurfaceList, linkMask } = sim;
  const rho0 = sim.RHO0 ?? 1.0;
  const xRef = sim.forceRef?.x ?? 0.0;
  const yRef = sim.forceRef?.y ?? 0.0;

  let fx = 0.0, fy = 0.0, mz = 0.0;

  for (let k = 0; k < fluidSurfaceList.length; k++) {
    const n    = fluidSurfaceList[k];
    const xn   = n % NX;
    const yn   = (n / NX) | 0;
    const base = n * Q;

    for (let d = 0; d < ALL_DIRS.length; d++) {
      const i   = ALL_DIRS[d];
      const off = base + i;
      if (!linkMask[off]) continue;

      const cx = CX[i];
      const cy = CY[i];
      const ex = fSnapshot[off] + f[base + OPP[i]] - 2.0 * W[i] * rho0;
      const dFx = cx * ex;
      const dFy = cy * ex;

      const ax = xn + 0.5 * cx - xRef;
      const ay = yn + 0.5 * cy - yRef;

      fx += dFx;
      fy += dFy;
      mz += ax * dFy - ay * dFx;
    }
  }

  return { fx, fy, mz };
}

// ---------------------------------------------------------
// Scaling
// ---------------------------------------------------------
function getScaling(sim) {
  // Dynamic pressure uses the CURRENT reference inlet velocity so the
  // coefficients stay correct when the inlet-velocity slider moves.
  // (Pulsation is deliberately excluded: u0 is the mean reference speed.)
  const rho0        = sim.RHO0 ?? 1.0;
  const U           = sim.params?.u0 ?? sim.U0 ?? 0.1;
  const L           = sim.forceRefLength ?? sim.L ?? 1.0;
  const qDyn        = 0.5 * rho0 * U * U;
  const forceScale  = qDyn * L;
  const momentScale = qDyn * L * L;

  sim.forceScale  = forceScale;
  sim.momentScale = momentScale;

  return { rho0, U, L, forceScale, momentScale };
}

function ensureForceStorage(sim) {
  if (!sim.force) sim.force = { fx: 0.0, fy: 0.0, mz: 0.0 };
  if (!sim.coeff) sim.coeff = { cd: 0.0, cl: 0.0, cm: 0.0 };
}

function zeroForces(sim) {
  sim.force.fx = sim.force.fy = sim.force.mz = 0.0;
  sim.coeff.cd = sim.coeff.cl = sim.coeff.cm = 0.0;
}

// ---------------------------------------------------------
// Debug
// ---------------------------------------------------------
export function debugPressureDistribution(sim) {
  if (!sim.fluidSurfaceList || !sim.linkMask || !sim.fSnapshot) {
    console.warn("[MEM debug] missing required sim fields"); return;
  }

  const { Q, f, fSnapshot, fluidSurfaceList, linkMask } = sim;
  const { rho0, U, L, forceScale } = getScaling(sim);
  const result = momentumExchange(sim);

  const byDir = {};
  for (const i of ALL_DIRS) byDir[i] = { links: 0, fx: 0, fy: 0 };

  for (let k = 0; k < fluidSurfaceList.length; k++) {
    const base = fluidSurfaceList[k] * Q;
    for (const i of ALL_DIRS) {
      const off = base + i;
      if (!linkMask[off]) continue;
      const ex = fSnapshot[off] + f[base + OPP[i]] - 2.0 * W[i] * rho0;
      byDir[i].links++;
      byDir[i].fx += CX[i] * ex;
      byDir[i].fy += CY[i] * ex;
    }
  }

  console.group(`[MEM debug] ${sim.id ?? "sim"}`);
  console.log("Scaling:", { rho0, U, L, forceScale });
  console.log("Net forces:", result);
  console.log("Cd:", result.fx / forceScale, "  Cl:", -result.fy / forceScale);
  console.table(ALL_DIRS.map(i => ({
    dir: i, cx: CX[i], cy: CY[i],
    links: byDir[i].links,
    fx: +byDir[i].fx.toFixed(6),
    fy: +byDir[i].fy.toFixed(6),
  })));
  console.groupEnd();
}
