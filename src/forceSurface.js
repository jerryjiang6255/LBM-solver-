// src/forceSurface.js
// ============================================================
// Surface-force integration.
//
// Hybrid approach:
// - Airfoils / NACA2412: integrate over ACTUAL wall-link points
//   from the simulated BFL geometry (fluidSurfaceList/linkMask/linkQ)
// - Cylinder / square / tandem cylinders: analytic contour integration
//
// Force model:
//   traction = -p * n + tau · n
//   p = cs^2 * rho = rho / 3
//
// Notes:
// - For airfoils, this avoids mismatch between an analytic contour
//   and the actual wall used by the solver.
// - This file is self-contained so future stress-model upgrades stay
//   here without requiring changes elsewhere.
// ============================================================

import { CX, CY } from "./lbm.js";

// ---------------------------------------------------------
// Configuration
// ---------------------------------------------------------
const CS2 = 1.0 / 3.0;
const FORCE_MODEL = "pressure+viscous";

// Sample slightly outside the body along outward normal
const SAMPLE_OFFSET = 0.35;

// Default contour resolution for analytic-shape integration
const DEFAULTS = {
  cylinderPoints: 180,
  squarePointsPerEdge: 60,
  tandemCylinderPoints: 160,
};

// ---------------------------------------------------------
// Public API
// ---------------------------------------------------------
export function computeSurfaceForces(sim) {
  ensureForceStorage(sim);

  const bodyType = sim?.params?.bodyType;

  let result;

  if (bodyType === "airfoil" || bodyType === "naca2412") {
    result = integrateForcesFromWallLinks(sim);
  } else {
    const contour = buildContourForSim(sim);

    if (
      !contour ||
      (
        !contour.multi &&
        (!Array.isArray(contour.points) || contour.points.length < 2)
      ) ||
      (
        contour.multi &&
        (!Array.isArray(contour.loops) ||
         contour.loops.length === 0 ||
         contour.loops.every(loop => !loop.points || loop.points.length < 2))
      )
    ) {
      result = { fx: 0.0, fy: 0.0, mz: 0.0 };
    } else {
      result = integrateSurfaceTraction(sim, contour);
    }
  }

  sim.force.fx = result.fx;
  sim.force.fy = result.fy;
  sim.force.mz = result.mz;

  const forceScale = sim.forceScale ?? 1.0;
  const momentScale =
    sim.momentScale ?? forceScale * Math.max(sim.L || 1.0, 1.0);

  sim.coeff.cd = forceScale !== 0 ? -result.fx / forceScale : 0.0;
  sim.coeff.cl = forceScale !== 0 ?  result.fy / forceScale : 0.0;
  sim.coeff.cm = momentScale !== 0 ? result.mz / momentScale : 0.0;
}

// ---------------------------------------------------------
// Airfoil / NACA force integration from actual wall links
// ---------------------------------------------------------
function integrateForcesFromWallLinks(sim) {
  const {
    NX, NY, Q,
    fluidSurfaceList,
    linkMask,
    linkQ,
    wallDist,
    forceRef,
  } = sim;

  if (!fluidSurfaceList || !linkMask || !linkQ || fluidSurfaceList.length === 0) {
    return { fx: 0.0, fy: 0.0, mz: 0.0 };
  }

  let fx = 0.0;
  let fy = 0.0;
  let mz = 0.0;

  const xRef = forceRef?.x ?? 0.0;
  const yRef = forceRef?.y ?? 0.0;

  // -------------------------------------------------------
  // Build traction contributions directly from each active
  // wall link:
  //   wall point: xw = x + q c_i
  //   normal:     n ~ outward from body into fluid = -c_i/|c_i|
  //
  // Weighting:
  //   Since each link is a directional wall intercept rather than
  //   an explicit segment, we assign a simple geometric weight ds:
  //   - cardinal links: ds = 1
  //   - diagonal links: ds = sqrt(2)
  //
  // This is a practical first-pass weighting consistent with
  // lattice spacing in D2Q9.
  // -------------------------------------------------------
  for (let k = 0; k < fluidSurfaceList.length; k++) {
    const n = fluidSurfaceList[k];
    const x = n % NX;
    const y = (n / NX) | 0;
    const base = n * Q;

    for (let i = 1; i < Q; i++) {
      const off = base + i;
      if (!linkMask[off]) continue;

      const q = linkQ[off];

      const cx = CX[i];
      const cy = CY[i];
      const cLen = Math.sqrt(cx * cx + cy * cy);

      // Wall point on actual solver-used boundary link
      const xw = x + q * cx;
      const yw = y + q * cy;

      // Outward normal from body into fluid
      // Direction i points fluid -> solid, so outward is -c_i
      const nx = -cx / cLen;
      const ny = -cy / cLen;

      // Approximate local segment weight
      const ds = cLen;

      const traction = computeTractionAtPoint(sim, xw, yw, nx, ny);

      const dFx = traction.tx * ds;
      const dFy = traction.ty * ds;

      fx += dFx;
      fy += dFy;
      mz += (xw - xRef) * dFy - (yw - yRef) * dFx;
    }
  }

  return { fx, fy, mz };
} 
// ---------------------------------------------------------
// Traction model
// ---------------------------------------------------------
function computeTractionAtPoint(sim, x, y, nx, ny) {
  const xs = x + SAMPLE_OFFSET * nx;
  const ys = y + SAMPLE_OFFSET * ny;

  const rho = sampleScalar(sim.rho, sim.NX, sim.NY, xs, ys);
  const p = CS2 * rho;

  let tx = -p * nx;
  let ty = -p * ny;

  if (FORCE_MODEL === "pressure+viscous") {
    const visc = computeViscousTraction(sim, xs, ys, nx, ny);
    tx += visc.tx;
    ty += visc.ty;
  }

  return { tx, ty, p };
}

// ---------------------------------------------------------
// Viscous traction
// tau = mu (grad u + grad u^T)
// traction_visc = tau · n
// ---------------------------------------------------------
function computeViscousTraction(sim, x, y, nx, ny) {
  const rho = sampleScalar(sim.rho, sim.NX, sim.NY, x, y);
  const nuEff = sim.NU0 ?? 0.0;
  const mu = rho * nuEff;

  const { dudx, dudy, dvdx, dvdy } = velocityGradients(sim, x, y);

  const tau_xx = 2.0 * mu * dudx;
  const tau_yy = 2.0 * mu * dvdy;
  const tau_xy = mu * (dudy + dvdx);

  return {
    tx: tau_xx * nx + tau_xy * ny,
    ty: tau_xy * nx + tau_yy * ny,
  };
}

// ---------------------------------------------------------
// Geometry dispatch for analytic-shape integration
// ---------------------------------------------------------
function buildContourForSim(sim) {
  const bodyType = sim?.params?.bodyType;

  if (bodyType === "cylinder") {
    return buildCylinderContour(sim);
  }

  if (bodyType === "square-cylinder") {
    return buildSquareContour(sim);
  }

  if (bodyType === "tandem-cylinders") {
    return buildTandemCylinderContour(sim);
  }

  return { points: [], closed: true };
}

// ---------------------------------------------------------
// Cylinder contour
// ---------------------------------------------------------
function buildCylinderContour(sim) {
  const NX = sim.NX;
  const NY = sim.NY;

  const cx = Math.floor(NX / 4);
  const cy = Math.floor(NY / 2);
  const radius = Math.max(4, Math.floor(NY / 8));

  return {
    points: buildCirclePoints(cx, cy, radius, DEFAULTS.cylinderPoints),
    closed: true,
  };
}

// ---------------------------------------------------------
// Tandem cylinders contour
// ---------------------------------------------------------
function buildTandemCylinderContour(sim) {
  const NX = sim.NX;
  const NY = sim.NY;

  const cx = Math.floor(NX * 0.30);
  const cy = Math.floor(NY / 2);
  const radius = Math.max(4, Math.floor(NY / 10));
  const spacing =
    sim.params?.tandemSpacing ?? Math.max(12, Math.floor(NX * 0.12));

  const cx1 = cx - spacing / 2;
  const cx2 = cx + spacing / 2;

  const c1 = buildCirclePoints(cx1, cy, radius, DEFAULTS.tandemCylinderPoints);
  const c2 = buildCirclePoints(cx2, cy, radius, DEFAULTS.tandemCylinderPoints);

  return {
    multi: true,
    loops: [
      { points: c1, closed: true },
      { points: c2, closed: true },
    ],
  };
}

// ---------------------------------------------------------
// Square contour
// ---------------------------------------------------------
function buildSquareContour(sim) {
  const NX = sim.NX;
  const NY = sim.NY;

  const cx = Math.floor(NX / 4);
  const cy = Math.floor(NY / 2);
  const size = Math.max(8, Math.floor(NY / 4));
  const alpha = (sim.params?.aoaDeg ?? 0) * Math.PI / 180;

  const half = size / 2;

  const cornersLocal = [
    { x: -half, y: -half },
    { x:  half, y: -half },
    { x:  half, y:  half },
    { x: -half, y:  half },
  ];

  const corners = cornersLocal.map(p =>
    rotateTranslate(p.x, p.y, cx, cy, alpha)
  );

  return {
    points: samplePolylineClosed(corners, DEFAULTS.squarePointsPerEdge),
    closed: true,
  };
}

// ---------------------------------------------------------
// Surface integration entry wrapper
// Handles single-loop and multi-loop analytic bodies
// ---------------------------------------------------------
function integrateSurfaceTraction(sim, contour) {
  if (contour.multi && Array.isArray(contour.loops)) {
    let fx = 0.0, fy = 0.0, mz = 0.0;
    for (const loop of contour.loops) {
      if (!loop.points || loop.points.length < 2) continue;
      const r = integrateSingleLoop(sim, loop);
      fx += r.fx;
      fy += r.fy;
      mz += r.mz;
    }
    return { fx, fy, mz };
  }

  return integrateSingleLoop(sim, contour);
}

function integrateSingleLoop(sim, contour) {
  const { points, closed = true } = contour;
  const { x: xRef = 0.0, y: yRef = 0.0 } = sim.forceRef || {};

  let fx = 0.0;
  let fy = 0.0;
  let mz = 0.0;

  const count = points.length;
  if (count < 2) return { fx, fy, mz };

  const segCount = closed ? count : count - 1;

  for (let k = 0; k < segCount; k++) {
    const a = points[k];
    const b = points[(k + 1) % count];

    const dx = b.x - a.x;
    const dy = b.y - a.y;
    const ds = Math.sqrt(dx * dx + dy * dy);
    if (ds <= 1e-12) continue;

    const xs = 0.5 * (a.x + b.x);
    const ys = 0.5 * (a.y + b.y);

    // For CCW contour in screen/lattice coordinates (y downward),
    // outward normal is (dy, -dx)/|t|
    const nx =  dy / ds;
    const ny = -dx / ds;

    const traction = computeTractionAtPoint(sim, xs, ys, nx, ny);

    const dFx = traction.tx * ds;
    const dFy = traction.ty * ds;

    fx += dFx;
    fy += dFy;
    mz += (xs - xRef) * dFy - (ys - yRef) * dFx;
  }

  return { fx, fy, mz };
}

// ---------------------------------------------------------
// Bilinear interpolation
// ---------------------------------------------------------
function sampleScalar(field, NX, NY, x, y) {
  const xc = clamp(x, 0, NX - 1.001);
  const yc = clamp(y, 0, NY - 1.001);

  const x0 = Math.floor(xc);
  const y0 = Math.floor(yc);
  const x1 = Math.min(x0 + 1, NX - 1);
  const y1 = Math.min(y0 + 1, NY - 1);

  const tx = xc - x0;
  const ty = yc - y0;

  const i00 = y0 * NX + x0;
  const i10 = y0 * NX + x1;
  const i01 = y1 * NX + x0;
  const i11 = y1 * NX + x1;

  const v00 = field[i00];
  const v10 = field[i10];
  const v01 = field[i01];
  const v11 = field[i11];

  const a = v00 * (1 - tx) + v10 * tx;
  const b = v01 * (1 - tx) + v11 * tx;
  return a * (1 - ty) + b * ty;
}

// ---------------------------------------------------------
// Velocity gradients by central difference on interpolated field
// ---------------------------------------------------------
function velocityGradients(sim, x, y) {
  const h = 1.0;

  const ux_xp = sampleScalar(sim.ux, sim.NX, sim.NY, x + h, y);
  const ux_xm = sampleScalar(sim.ux, sim.NX, sim.NY, x - h, y);
  const ux_yp = sampleScalar(sim.ux, sim.NX, sim.NY, x, y + h);
  const ux_ym = sampleScalar(sim.ux, sim.NX, sim.NY, x, y - h);

  const uy_xp = sampleScalar(sim.uy, sim.NX, sim.NY, x + h, y);
  const uy_xm = sampleScalar(sim.uy, sim.NX, sim.NY, x - h, y);
  const uy_yp = sampleScalar(sim.uy, sim.NX, sim.NY, x, y + h);
  const uy_ym = sampleScalar(sim.uy, sim.NX, sim.NY, x, y - h);

  return {
    dudx: 0.5 * (ux_xp - ux_xm) / h,
    dudy: 0.5 * (ux_yp - ux_ym) / h,
    dvdx: 0.5 * (uy_xp - uy_xm) / h,
    dvdy: 0.5 * (uy_yp - uy_ym) / h,
  };
} 
// ---------------------------------------------------------
// Helpers
// ---------------------------------------------------------
function ensureForceStorage(sim) {
  if (!sim.force) {
    sim.force = { fx: 0.0, fy: 0.0, mz: 0.0 };
  }
  if (!sim.coeff) {
    sim.coeff = { cd: 0.0, cl: 0.0, cm: 0.0 };
  }
}

function clamp(v, lo, hi) {
  return Math.max(lo, Math.min(hi, v));
}

function rotateTranslate(x, y, cx, cy, alpha) {
  const c = Math.cos(alpha);
  const s = Math.sin(alpha);
  return {
    x: cx + x * c - y * s,
    y: cy + x * s + y * c,
  };
}

function buildCirclePoints(cx, cy, r, Np) {
  const pts = [];
  for (let k = 0; k < Np; k++) {
    const th = (2.0 * Math.PI * k) / Np;
    pts.push({
      x: cx + r * Math.cos(th),
      y: cy + r * Math.sin(th),
    });
  }
  return pts;
}

function samplePolylineClosed(corners, ptsPerEdge) {
  const pts = [];
  const m = corners.length;

  for (let e = 0; e < m; e++) {
    const a = corners[e];
    const b = corners[(e + 1) % m];

    for (let j = 0; j < ptsPerEdge; j++) {
      const t = j / ptsPerEdge;
      pts.push({
        x: a.x * (1 - t) + b.x * t,
        y: a.y * (1 - t) + b.y * t,
      });
    }
  }

  return pts;
} 
