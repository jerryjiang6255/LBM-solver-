// src/lbm.js  (V2)
// ============================================================
// Core LBM data + math layer for the V2 solver: D2Q9 lattice,
// MRT moment-space transform, equilibrium distributions AND
// equilibrium moments, macroscopic recovery, and all shared
// state arrays. This is the "foundation" module — geometry.js,
// turbulence.js, collision.js, streaming.js, and boundary.js
// all read from and write into the arrays defined here, but
// this file itself depends on nothing else in the project.
//
// V2 SCOPE vs V1:
//   - Adds the MRT moment matrix M and its (precomputed) inverse,
//     used by collision.js for moment-space relaxation instead of
//     V1's single-rate BGK relaxation.
//   - Adds fEq/fNeq scratch arrays, needed by turbulence.js to
//     compute the strain-rate tensor from non-equilibrium parts
//     of the distribution (fneq = f - feq).
//   - Adds moment scratch arrays (m, meq, mpost), needed by
//     collision.js to do relax-in-moment-space without allocating
//     new arrays every frame.
//   - Splits what V1 did inside one collide() function into
//     separate computeMacros() / computeEquilibrium() steps,
//     because V2's per-frame workflow needs macroscopic fields
//     and non-equilibrium data available to turbulence.js BEFORE
//     collision.js runs (V1 could compute+collide in one pass
//     since there was no LES/wall-model step in between).
// ============================================================


// ---------------------------------------------------------
// GRID CONSTANTS
// ---------------------------------------------------------
export const NX = 256;
export const NY = 128;
export const N = NX * NY;
export const Q = 9;


export const RHO0 = 1.0;
export const U0 = 0.06; // lattice-unit inlet velocity target


export let L = Math.max(4, Math.floor(NY / 8)) * 2; // placeholder: 2*default radius
export let RE = 200;
export let NU0 = (U0 * L) / RE;
export let TAU0 = 3 * NU0 + 0.5;
export let OMEGA0 = 1 / TAU0;


// ---------------------------------------------------------
// setCharacteristicLength(diameter) / setReynolds(re)
// ---------------------------------------------------------
export function setCharacteristicLength(diameter) {
  L = diameter;
  recomputeBaseViscosity();
}
export function setReynolds(re) {
  RE = re;
  recomputeBaseViscosity();
}
function recomputeBaseViscosity() {
  NU0 = (U0 * L) / RE;
  TAU0 = 3 * NU0 + 0.5;
  OMEGA0 = 1 / TAU0;
}

// ---------------------------------------------------------
// D2Q9 LATTICE VECTORS (unchanged from V1)
// ---------------------------------------------------------
export const CX = new Int8Array([0, 1, 0, -1, 0, 1, -1, -1, 1]);
export const CY = new Int8Array([0, 0, 1, 0, -1, 1, 1, -1, -1]);
export const W = new Float32Array([
  4 / 9,
  1 / 9, 1 / 9, 1 / 9, 1 / 9,
  1 / 36, 1 / 36, 1 / 36, 1 / 36,
]);
export const OPP = new Uint8Array([0, 3, 4, 1, 2, 7, 8, 5, 6]);

// ---------------------------------------------------------
// STATE ARRAYS
// ---------------------------------------------------------
export const f = new Float32Array(Q * N);
// export const fNew = new Float32Array(Q * N);
export const rho = new Float32Array(N);
export const ux = new Float32Array(N);
export const uy = new Float32Array(N);
// export const fEq = new Float32Array(Q * N);
// export const fNeq = new Float32Array(Q * N);
// export const m = new Float32Array(Q * N);
// export const meq = new Float32Array(Q * N);
// export const mpost = new Float32Array(Q * N);


// ---------------------------------------------------------
// idx(x, y)
// ---------------------------------------------------------
// Same as V1: flat node index from 2D lattice coordinates.
export function idx(x, y) {
  return y * NX + x;
}

// ---------------------------------------------------------
// feqAt(i, r, u, v)
// ---------------------------------------------------------
export function feqAt(i, r, u, v) {
  const cu = CX[i] * u + CY[i] * v;
  const uu = u * u + v * v;
  return W[i] * r * (1 + 3 * cu + 4.5 * cu * cu - 1.5 * uu);
}


// ---------------------------------------------------------
// initField(u0, solid)
// ---------------------------------------------------------
export function initField(u0, solid) {
  for (let n = 0; n < N; n++) {
    rho[n] = RHO0;
    if (solid[n]) {
      ux[n] = 0;
      uy[n] = 0;
    } else {
      ux[n] = u0;
      uy[n] = 0;
    }
  }
  for (let n = 0; n < N; n++) {
    for (let i = 0; i < Q; i++) {
      f[n * Q + i] = feqAt(i, rho[n], ux[n], uy[n]);
    }
  }
}

// ---------------------------------------------------------
// computeMacros(solid)
// ---------------------------------------------------------
export function computeMacros(solid) {
  for (let n = 0; n < N; n++) {
    if (solid[n]) {
      rho[n] = RHO0;
      ux[n] = 0;
      uy[n] = 0;
      continue;
    }
    let r = 0, jx = 0, jy = 0;
    for (let i = 0; i < Q; i++) {
      const fi = f[n * Q + i];
      r += fi;
      jx += fi * CX[i];
      jy += fi * CY[i];
    }
    rho[n] = r;
    ux[n] = jx / r;
    uy[n] = jy / r;
  }
}


// ---------------------------------------------------------
// computeEquilibriumAndNonEquilibrium()
// ---------------------------------------------------------
export function computeEquilibriumAndNonEquilibrium() {
  for (let n = 0; n < N; n++) {
    const r = rho[n];
    const u = ux[n];
    const v = uy[n];
    for (let i = 0; i < Q; i++) {
      const off = n * Q + i;
      const eq = feqAt(i, r, u, v);
      fEq[off] = eq;
      fNeq[off] = f[off] - eq;
    }
  }
}


// ---------------------------------------------------------
// checkStability(solid)
// ---------------------------------------------------------
export function checkStability(solid) {
  let corrected = 0;
  for (let n = 0; n < N; n++) {
    if (solid[n]) continue;
    const r = rho[n];
    if (!Number.isFinite(r) || r < 0.2 || r > 5) {
      rho[n] = 1.0;
      ux[n] = 0;
      uy[n] = 0;
      for (let i = 0; i < Q; i++) {
        f[n * Q + i] = feqAt(i, 1.0, 0, 0);
      }
      corrected++;
    }
  }
  return corrected;
}

// ---------------------------------------------------------
// GRID CONSTANTS
// ---------------------------------------------------------
// NX, NY : lattice dimensions (columns, rows). Same 256x128
//          domain as V1 unless you resize for a specific case.
// N      : total node count.
// Q      : D2Q9 velocity count.
// RHO0   : reference/rest density used in equilibrium-moment
//          formulas (jx^2/RHO0 etc). Kept as 1.0 unless you're
//          deliberately modeling compressible effects.
// U0     : characteristic inlet velocity, lattice units. Used
//          only to derive NU0 from the target Reynolds number —
//          the actual per-node inlet velocity is still whatever
//          boundary.js's Zou-He inlet prescribes each frame.
// L      : characteristic length scale, lattice units (here,
//          the cylinder diameter is the natural choice since
//          Reynolds number for flow past a cylinder is
//          conventionally defined off diameter, not domain
//          height). Set to match geometry.js's cylinder radius
//          at startup (see NOTE below).
// RE     : target Reynolds number. This is now a first-class
//          physical parameter in V2 (V1 only had TAU directly);
//          exposing Re instead of tau is more physically
//          meaningful and matches how the workflow doc defines
//          nu0 = U*L/Re.
// NU0    : molecular (laminar) kinematic viscosity, derived from
//          Re, U0, L. This feeds turbulence.js's total-viscosity
//          calculation (nu_tot = nu0 + nu_t) and is the baseline
//          relaxation rate before any LES/wall-model correction.
// TAU0, OMEGA0 : baseline relaxation time/rate corresponding to
//          NU0 alone (nu = (tau-0.5)/3 => tau = 3*nu + 0.5). These
//          are what collision.js falls back to for the CONSERVED
//          moments' relaxation rates (which are always 0, i.e.
//          never relaxed) and give turbulence.js a "no LES" tau
//          to start its correction from.
//
// NOTE: L defaults to a placeholder here; index.html or main.js
// should overwrite it (via setCharacteristicLength()) right after
// geometry.js builds the cylinder, so NU0/TAU0 reflect the actual
// obstacle size rather than a guess made before geometry exists.

// ---------------------------------------------------------
// setCharacteristicLength(diameter) / setReynolds(re)
// ---------------------------------------------------------
// Recompute NU0/TAU0/OMEGA0 after geometry or target Re changes.
// Kept as explicit setters (rather than recomputing NU0 inline
// everywhere) so there is exactly one place this derivation
// logic lives — call these once at startup and again whenever a
// V2 UI slider for Reynolds number changes.
// Params:
//   diameter : cylinder diameter in lattice units (2*radius)
//   re       : target Reynolds number

// ---------------------------------------------------------
// D2Q9 LATTICE VECTORS (unchanged from V1)
// ---------------------------------------------------------
// Direction layout: 0 rest, 1 E, 2 N, 3 W, 4 S, 5 NE, 6 NW,
// 7 SW, 8 SE. This exact ordering is NOT arbitrary in V2 — the
// MRT moment matrix M below is only valid for this specific
// CX/CY/W ordering (its rows encode Hermite-polynomial moments
// evaluated at these exact velocity vectors). If you ever
// reorder directions, M and MRT_MINV must be rebuilt to match.

// ---------------------------------------------------------
// STATE ARRAYS
// ---------------------------------------------------------
// f, fNew : distribution functions, flat Float32Array(Q*N),
//           f[i*N+n]. Same layout/rationale as V1 (double
//           buffering needed for streaming; see streaming.js).
// rho, ux, uy : macroscopic fields, Float32Array(N) each.
// fEq     : equilibrium distribution at the CURRENT macroscopic
//           state, Float32Array(Q*N). Needed (unlike V1) as a
//           standing array rather than a scratch value, because
//           turbulence.js needs fEq subtracted from f to get
//           fneq for the strain-rate tensor, computed in a
//           dedicated pass before collision.js runs.
// fNeq    : non-equilibrium part, f - fEq, Float32Array(Q*N).
//           Read by turbulence.js's strain-rate calculation.
// m, meq, mpost : moment-space scratch arrays, Float32Array(Q*N)
//           each, m[k*N+n] = k-th moment at node n. Used by
//           collision.js to avoid allocating temporary 9-element
//           arrays per node per frame (a correctness/performance
//           risk explicitly flagged in the V2 README under
//           "collision/mrt.js allocates JS arrays in hot loops" —
//           preallocating these here avoids that mistake from
//           the start).

// ---------------------------------------------------------
// feqAt(i, r, u, v)
// ---------------------------------------------------------
// Same second-order equilibrium expansion as V1 — this formula
// is unchanged between V1 and V2; MRT changes HOW the collision
// relaxes toward equilibrium (per-moment rates instead of one
// global rate), not what the equilibrium itself is.
// Params: i (direction), r (rho), u/v (velocity components).
// Returns: equilibrium distribution value for that direction.

// ---------------------------------------------------------
// initField(u0, solid)
// ---------------------------------------------------------
// Initializes the lattice per the workflow doc's t=0 spec:
//   Fluid cells: rho=1, u=uinlet
//   Solid cells: rho=1, u=0
// then sets f_i = f_i^eq everywhere (hydrodynamic equilibrium
// start state — no shock on frame 1).
// Params:
//   u0    : target inlet horizontal velocity, lattice units.
//   solid : Uint8Array(N) solid mask from geometry.js. Solid
//           nodes get zero velocity rather than u0, matching the
//           doc's explicit "solid cells ... u=0" initialization
//           (distinct from V1, which didn't special-case solid
//           nodes here since collide() skipped them anyway — in
//           V2, fEq/fNeq are computed over ALL nodes including
//           solid ones for array-uniformity in the strain-rate
//           pass, so solid nodes need a sane, non-garbage initial
//           state too).
// Side effects: fills rho, ux, uy, f in place for every node.

// ---------------------------------------------------------
// computeMacros(solid)
// ---------------------------------------------------------
// Recovers density and velocity from the current distribution
// field. Split out as its own function (rather than fused into
// collision, as V1 did) because V2's per-frame order needs
// rho/ux/uy finalized BEFORE turbulence.js's strain-rate and
// wall-model steps run, and those in turn must finish before
// collision.js does its moment-space relaxation.
// Params:
//   solid : Uint8Array(N). Solid nodes are set to the same rest
//           state used at init (rho=1, u=0) rather than left
//           stale or computed from meaningless streamed-in data —
//           this keeps fEq/fNeq well-defined everywhere, which
//           turbulence.js relies on when it scans wall-adjacent
//           cells generically.
// Side effects: overwrites rho, ux, uy in place.


// ---------------------------------------------------------
// computeEquilibriumAndNonEquilibrium()
// ---------------------------------------------------------
// Fills fEq (equilibrium distribution at the CURRENT rho/ux/uy)
// and fNeq (= f - fEq) for every node/direction. This is the
// direct analog of the workflow doc's "compute equilibrium
// distributions... fineq = fi - fieq" step, run right after
// computeMacros() and right before turbulence.js's strain-rate
// calculation (which consumes fNeq).
// Params: none — reads current f/rho/ux/uy module state.
// Side effects: overwrites fEq and fNeq in place, every node
// and direction.
// Note: this is NOT redundant with collision.js's own
// equilibrium evaluation — collision.js needs f_i^eq expressed
// in MOMENT space (via equilibriumMoments below) for the MRT
// relax step, while this function produces the DISTRIBUTION-
// space equilibrium needed for the strain-rate tensor. Both are
// derived from the same rho/ux/uy but serve different consumers.

// ---------------------------------------------------------
// computeMoments(solid)
// ---------------------------------------------------------
// Transforms the current distribution field into moment space:
// m[k*N+n] = sum_i MRT_M[k*Q+i] * f[i*N+n], for every node and
// every one of the 9 moments. This is the m = M*f step from the
// workflow doc, computed once per frame for ALL fluid nodes so
// collision.js can just read m[][] rather than recomputing the
// transform itself.
// Params:
//   solid : Uint8Array(N); solid nodes are skipped (left at
//           whatever m held previously — collision.js also
//           skips solid nodes, matching V1's convention of never
//           touching solid-node populations there, since
//           boundary.js/geometry.js's BFL step owns them).
// Side effects: overwrites m in place for fluid nodes.

// ---------------------------------------------------------
// computeEquilibriumMoments(solid)
// ---------------------------------------------------------
// Fills meq using the CLOSED-FORM analytic equilibrium-moment
// formulas from the workflow doc, rather than transforming
// fEq through M (mathematically equivalent, but the closed form
// is direct arithmetic — no 9x9 matrix-vector product needed —
// and matches exactly what the doc specifies):
//   m0eq = rho
//   m1eq = -2*rho + 3*(jx^2+jy^2)/rho
//   m2eq =  rho  - 3*(jx^2+jy^2)/rho
//   m3eq =  jx
//   m4eq = -jx
//   m5eq =  jy
//   m6eq = -jy
//   m7eq = (jx^2 - jy^2)/rho
//   m8eq =  jx*jy/rho
// where jx = rho*ux, jy = rho*uy (momentum, not velocity).
// Params:
//   solid : Uint8Array(N); solid nodes skipped, same convention
//           as computeMoments.
// Side effects: overwrites meq in place for fluid nodes.

// ---------------------------------------------------------
// momentsToF(solid)
// ---------------------------------------------------------
// Transforms relaxed moments back into distribution space:
// fNew_i[n] = sum_k MRT_MINV[i*Q+k] * mpost[k*N+n]. This is the
// f_post = M^-1 * m_post step from the workflow doc. Writes into
// `f` directly (post-collision, pre-streaming state) — this is
// the LAST step collision.js needs from lbm.js each frame; after
// this call, streaming.js takes over.
// Params:
//   solid : Uint8Array(N); solid nodes left untouched (their f
//           values are governed entirely by geometry.js's BFL
//           reconstruction after streaming, not by collision).
// Side effects: overwrites f in place for fluid nodes.

// ---------------------------------------------------------
// checkStability(solid)
// ---------------------------------------------------------
// Scans rho for NaN/Infinity or values outside a physically
// sane range, and re-equilibrates any offending cell back to
// rest state. This stops a local instability from being able to
// grow across frames, regardless of where it starts — a direct,
// general safety net rather than patching specific known trouble
// spots one at a time.
// Params: solid : Uint8Array(N)
// Returns: number of cells that were corrected (0 = all healthy)