// src/forceSurface.js
// ============================================================
// Force via Normal-Filtered Fluctuation Momentum Exchange Method.
//
// Formula per blocked link i at fluid node n:
//   ΔF = 2 * c_i * (f_i*(n) - W[i]*rho0)
//
// Link i is included only if:
//   dot(c_i, n_body) > 0
//
// where n_body = local outward normal of the solid body,
// estimated from the SDF gradient at node n using finite differences
// on the sdf array (which geometry.js already computes and stores).
//
// WHY NORMAL FILTERING:
//   - All-directions MEM overcounts on flat surfaces (square at 0°):
//     diagonal links from a flat front face have CX=+1 and pile onto
//     drag from the perpendicular cardinal link, giving Cd ~2.5× too high.
//   - Cardinal-only MEM undercounts on curved/angled surfaces (airfoil):
//     a surface at 45° has NO cardinal blocked links, only diagonals,
//     so cardinal-only misses all the real momentum transfer there.
//   - Normal filtering solves both: it includes exactly the links that
//     point "into" the body at that surface location, regardless of
//     whether they are cardinal or diagonal. On a flat wall perpendicular
//     to x, only i=1 passes (correct). On a 45° surface, only i=5 or i=8
//     passes (also correct).
//
// WHY SUBTRACT W[i]*rho0 (fluctuation form):
//   f_i* ≈ W[i]*rho0 + δf_i where δf_i << W[i]*rho0.
//   Summing the full f_i* gives front contribution ~+7 and rear ~-7,
//   with the physical drag signal buried in their difference (~0.025).
//   Subtracting the common background W[i]*rho0 before summing removes
//   the ±7 offset, giving direct access to the ~0.025 signal without
//   catastrophic floating-point cancellation.
//
// Requires:
//   sim.fSnapshot  — post-streaming f, captured in streaming.js applyBFL()
//                    before the BFL wall reconstruction overwrites links
//   sim.sdf        — signed distance field from geometry builder
//   sim.RHO0       — background density (default 1.0)
// ============================================================

import { CX, CY, W } from "./lbm.js";

const ALL_DIRS = [1, 2, 3, 4, 5, 6, 7, 8];

// Cosine threshold for including a link.
// dot(c_i_normalized, n_body) > COS_THRESH means the link points
// meaningfully into the body (not tangential or away from it).
// 0.0 = include anything with any inward component (hemisphere).
// Using a small positive value avoids counting tangential links.
const COS_THRESH = 0.1;

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
// Normal-filtered fluctuation MEM
// ---------------------------------------------------------
function momentumExchange(sim) {
  const { NX, NY, Q, fSnapshot, fluidSurfaceList, linkMask, sdf } = sim;
  const rho0 = sim.RHO0 ?? 1.0;
  const xRef = sim.forceRef?.x ?? 0.0;
  const yRef = sim.forceRef?.y ?? 0.0;

  // Whether we have an SDF to compute local normals from.
  // Without SDF we fall back to accepting all inward-facing links
  // (equivalent to all-directions MEM — better than nothing).
  const hasSdf = sdf && sdf.length > 0;

  let fx = 0.0, fy = 0.0, mz = 0.0;

  for (let k = 0; k < fluidSurfaceList.length; k++) {
    const n    = fluidSurfaceList[k];
    const xn   = n % NX;
    const yn   = (n / NX) | 0;
    const base = n * Q;

    // Estimate local body outward normal from SDF gradient.
    // grad(SDF) points from solid toward fluid = outward body normal.
    let nx_body = 0.0, ny_body = 0.0;
    if (hasSdf) {
      const xl = xn > 0      ? xn - 1 : xn;
      const xr = xn < NX - 1 ? xn + 1 : xn;
      const yd = yn > 0      ? yn - 1 : yn;
      const yu = yn < NY - 1 ? yn + 1 : yn;

      const gx = (sdf[yn  * NX + xr] - sdf[yn  * NX + xl]) / (xr - xl);
      const gy = (sdf[yu  * NX + xn] - sdf[yd  * NX + xn]) / (yu - yd);
      const glen = Math.sqrt(gx * gx + gy * gy);
      if (glen > 1e-10) {
        nx_body = gx / glen;
        ny_body = gy / glen;
      }
    }

    for (let d = 0; d < ALL_DIRS.length; d++) {
      const i   = ALL_DIRS[d];
      const off = base + i;
      if (!linkMask[off]) continue;

      const cx = CX[i];
      const cy = CY[i];

      // Normal filter: only count this link if it points into the body.
      // The lattice velocity c_i = (cx, cy) points FROM fluid TOWARD solid.
      // The body outward normal n_body points FROM solid TOWARD fluid.
      // For a valid inward-momentum link: dot(c_i, -n_body) > 0
      //   i.e. the link direction has a component OPPOSING the body normal
      //   (pointing into the body surface, not away from it or tangentially).
      if (hasSdf) {
        // Magnitude of c_i: 1 for cardinal, sqrt(2) for diagonal
        const cMag  = Math.sqrt(cx * cx + cy * cy);
        // dot(c_i, -n_body) / |c_i| = -(cx*nx_body + cy*ny_body) / cMag
        const cosA  = -(cx * nx_body + cy * ny_body) / cMag;
        if (cosA <= COS_THRESH) continue;
      }

      // Fluctuation MEM: subtract rest-equilibrium background
      const fi  = fSnapshot[off] - W[i] * rho0;
      const dFx = 2.0 * cx * fi;
      const dFy = 2.0 * cy * fi;

      const ax  = xn + 0.5 * cx - xRef;
      const ay  = yn + 0.5 * cy - yRef;

      fx  += dFx;
      fy  += dFy;
      mz  += ax * dFy - ay * dFx;
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

  const { NX, NY, Q, fSnapshot, fluidSurfaceList, linkMask, sdf } = sim;
  const { rho0, U, L, forceScale } = getScaling(sim);
  const hasSdf = sdf && sdf.length > 0;
  const result = momentumExchange(sim);

  const byDir = {};
  for (const i of ALL_DIRS) byDir[i] = { included: 0, excluded: 0, fx: 0, fy: 0 };

  for (let k = 0; k < fluidSurfaceList.length; k++) {
    const n    = fluidSurfaceList[k];
    const xn   = n % NX;
    const yn   = (n / NX) | 0;
    const base = n * Q;

    let nx_body = 0, ny_body = 0;
    if (hasSdf) {
      const xl = xn > 0      ? xn - 1 : xn;
      const xr = xn < NX - 1 ? xn + 1 : xn;
      const yd = yn > 0      ? yn - 1 : yn;
      const yu = yn < NY - 1 ? yn + 1 : yn;
      const gx = (sdf[yn * NX + xr] - sdf[yn * NX + xl]) / (xr - xl);
      const gy = (sdf[yu * NX + xn] - sdf[yd * NX + xn]) / (yu - yd);
      const glen = Math.sqrt(gx * gx + gy * gy);
      if (glen > 1e-10) { nx_body = gx / glen; ny_body = gy / glen; }
    }

    for (const i of ALL_DIRS) {
      const off = base + i;
      if (!linkMask[off]) continue;
      const cx = CX[i], cy = CY[i];
      const cMag = Math.sqrt(cx * cx + cy * cy);
      const cosA = hasSdf ? -(cx * nx_body + cy * ny_body) / cMag : 1.0;

      if (cosA <= COS_THRESH) {
        byDir[i].excluded++;
        continue;
      }
      const fi = fSnapshot[off] - W[i] * rho0;
      byDir[i].included++;
      byDir[i].fx += 2 * cx * fi;
      byDir[i].fy += 2 * cy * fi;
    }
  }

  console.group(`[MEM debug] ${sim.id ?? "sim"}`);
  console.log("Scaling:", { rho0, U, L, forceScale });
  console.log("Net forces:", result);
  console.log("Cd:", result.fx / forceScale, "  Cl:", -result.fy / forceScale);
  console.log("SDF available:", hasSdf, "  COS_THRESH:", COS_THRESH);
  console.table(ALL_DIRS.map(i => ({
    dir: i, cx: CX[i], cy: CY[i],
    included: byDir[i].included,
    excluded: byDir[i].excluded,
    fx: +byDir[i].fx.toFixed(6),
    fy: +byDir[i].fy.toFixed(6),
  })));
  console.groupEnd();
}