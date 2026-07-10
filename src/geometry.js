// src/geometry.js  (V2)
// ============================================================
// Cylinder geometry, analytic (closed-form) wall distance, and
// analytic BFL q-values. No SDF grid, no ray-marching — since
// V2's obstacle is a circle, every geometric query has a direct
// formula, so we skip the general SDF infrastructure the big
// V2 README used (that's for arbitrary future shapes; V2 here
// only needs a cylinder).
// ============================================================


import { NX, NY, N, Q, CX, CY, idx, setCharacteristicLength } from "./lbm.js";


// ---------------------------------------------------------
// buildCylinder(cx, cy, radius)
// ---------------------------------------------------------
export function buildCylinder(cx, cy, radius) {
  const solid = new Uint8Array(N);
  const dist = new Float32Array(N);


  for (let y = 0; y < NY; y++) {
    for (let x = 0; x < NX; x++) {
      const n = idx(x, y);
      const d = Math.sqrt((x - cx) ** 2 + (y - cy) ** 2) - radius; // signed dist to circle
      dist[n] = -d; // flip sign: positive = fluid, negative = solid
      solid[n] = d <= 0 ? 1 : 0;
    }
  }


  setCharacteristicLength(2 * radius);
  return { solid, dist };
}


// ---------------------------------------------------------
// defaultCylinder()
// ---------------------------------------------------------
export function defaultCylinder() {
  const cx = Math.floor(NX / 4);
  const cy = Math.floor(NY / 2);
  const radius = Math.max(4, Math.floor(NY / 8));
  return buildCylinder(cx, cy, radius);
}


// ---------------------------------------------------------
// buildWallDistanceField(dist, solid)
// ---------------------------------------------------------
export function buildWallDistanceField(dist, solid) {
  const wallDist = new Float32Array(N);
  for (let n = 0; n < N; n++) {
    wallDist[n] = solid[n] ? 0 : dist[n];
  }
  return wallDist;
}


// ---------------------------------------------------------
// findWallAdjacent(solid)
// ---------------------------------------------------------
export function findWallAdjacent(solid) {
  const adj = new Uint8Array(N);
  for (let y = 0; y < NY; y++) {
    for (let x = 0; x < NX; x++) {
      const n = idx(x, y);
      if (solid[n]) continue;

      for (let i = 1; i < Q; i++) {
        const nx = x + CX[i];
        const ny = y + CY[i];
        if (nx < 0 || nx >= NX || ny < 0 || ny >= NY) continue;
        if (solid[idx(nx, ny)]) {
          adj[n] = 1;
          break;
        }
      }
    }
  }
  return adj;
}

// buildSolidList(solid)

export function buildSolidList(solid) {
  const tmp = [];
  for (let y = 1; y < NY - 1; y++) {
    for (let x = 1; x < NX - 1; x++) {
      const n = idx(x, y);
      if (!solid[n]) continue;
      // only include if at least one axis-aligned neighbor is fluid
      if (
        !solid[n - 1]    ||   // W
        !solid[n + 1]    ||   // E
        !solid[n - NX]   ||   // S
        !solid[n + NX]        // N
      ) {
        tmp.push(n);
      }
    }
  }
  return new Uint32Array(tmp);
}

// ---------------------------------------------------------
// buildCylinder(cx, cy, radius)
// ---------------------------------------------------------
// Builds the solid mask + per-node signed distance to the
// circle (positive = fluid, negative = solid, matching the
// sign convention other V2 modules expect).
// Params: cx, cy (center, lattice units), radius (lattice units).
// Returns: { solid: Uint8Array(N), dist: Float32Array(N) }
// Also calls setCharacteristicLength(2*radius) so lbm.js's
// NU0/TAU0 derive from the actual obstacle diameter.

// ---------------------------------------------------------
// defaultCylinder()
// ---------------------------------------------------------
// Same default placement as V1: centered vertically, quarter
// of the way in from the left, radius ~NY/8.
// Returns: { solid, dist } (see buildCylinder).

// ---------------------------------------------------------
// buildWallDistanceField(dist, solid)
// ---------------------------------------------------------
// Wall distance for fluid cells = their signed distance to the
// circle (already positive by construction in buildCylinder).
// Solid cells get 0. Used by turbulence.js's wall model (y in
// the log-law formula — see workflow doc, "Distance to wall").
// Params: dist, solid — outputs of buildCylinder.
// Returns: Float32Array(N) wall distance, 0 inside solid.

// ---------------------------------------------------------
// findWallAdjacent(solid)
// ---------------------------------------------------------
// Flags fluid cells that have at least one solid orthogonal
// neighbor — these are the cells the wall model applies to
// (workflow doc: "applies only in fluid cells adjacent to solid
// walls... same cells as bounce back fluid cells").
// Params: solid : Uint8Array(N).
// Returns: Uint8Array(N), 1 = wall-adjacent fluid cell.


// Returns a compact Uint32Array of node indices where solid[n] === 1
// AND at least one face-connected neighbor is fluid.
// This is the V2 equivalent of the reference's barList —
// used by streaming.js's bounce-back so it only visits the
// ~500 surface cells instead of all N nodes.