// src/collisionFused.js  (V2 — full LES + wall model + TRT collision)
// ============================================================
// TRT replaces BGK as the collision operator. Physics motivation:
//
//   BGK: one relaxation rate omega for ALL deviations from eq.
//        Ties viscosity (physical) and numerical damping (ghost)
//        to the same tau — no independent control over either.
//
//   MRT: 9 independent rates in moment space. Most accurate but
//        requires 2× full 9×9 matrix transforms per node per
//        frame (~162 multiply-adds just for the transforms).
//
//   TRT: two rates — omega+ (even/symmetric, sets viscosity) and
//        omega- (odd/asymmetric, controls ghost/numerical modes).
//        No matrices — decomposition is just arithmetic on
//        opposite-direction pairs. Cost over BGK: ~9 extra
//        adds and 2 multiplies per node (negligible).
//        Stability over BGK: the Magic Parameter Lambda = 3/16
//        eliminates parasitic checkerboard and boundary-layer
//        artifacts that BGK develops at low viscosity.
//
// TRT even/odd decomposition (per opposite pair i / ī):
//   fi+  = (fi  + fī ) / 2   (symmetric part)
//   fi-  = (fi  - fī ) / 2   (antisymmetric part)
//   feqi+ = (feqi + feqī) / 2
//   feqi- = (feqi - feqī) / 2
//
// TRT relaxation:
//   fi+,post = fi+  - omega+ * (fi+  - feqi+)
//   fi-,post = fi-  - omega- * (fi-  - feqi-)
//
// Recombine:
//   fipost  = fi+,post + fi-,post
//   fīpost  = fi+,post - fi-,post
//   f0post  = f0 - omega+ * (f0 - f0eq)   (rest dir, symmetric only)
//
// omega+ is set by viscosity (same as BGK):
//   omega+ = 1 / (3*nu_tot + 0.5)
//
// omega- is set by the Magic Parameter Lambda = 3/16:
//   Lambda = (1/omega+ - 0.5) * (1/omega- - 0.5) = 3/16
//   => omega- = 1 / (Lambda/(1/omega+ - 0.5) + 0.5)
//
// Lambda = 3/16 is the unique value that makes bounce-back walls
// produce exact no-slip at the halfway point between nodes
// (eliminating the leading-order position error that BGK has).
//
// Loop structure: still 1N × 3Q direction ops total.
//   Loop 1 (Q): macros
//   Loop 2 (Q): feq + fneq stress (LES strain rate)
//   Scalar: LES nu_t, wall model, TRT omega+/omega-
//   Loop 3 (Q/2 pairs + rest): TRT relax + recombine
//     Loop 3 iterates 4 pairs (i=1..4) + rest (i=0) = 5 iters
//     but does 2 reads+writes per pair = equivalent to Q=9 ops.
// ============================================================

import { N, Q, CX, CY, W, OPP, f, rho, ux, uy, NU0 } from "./lbm.js";

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
// Lambda = (1/omega+ - 0.5) * (1/omega- - 0.5)
// Lambda = 3/16 is the "magic" value: makes the effective
// wall position of bounce-back independent of viscosity,
// eliminating the leading-order boundary error that BGK has.
// Other published choices: 1/4 (minimizes third-order error
// in bulk), 1/12 (exact advection). 3/16 is the standard
// choice for channel flow / bluff bodies — use it here.
const LAMBDA = 3.0 / 16.0;

// ---------------------------------------------------------
// D2Q9 opposite-pair table for TRT decomposition.
// Only the 4 pairs (i, OPP[i]) for i = 1..4 are needed;
// direction 0 (rest) is its own opposite and handled separately.
// Pairs: (1,3) E-W, (2,4) N-S, (5,7) NE-SW, (6,8) NW-SE
// Precomputed as a plain array of [i, opp_i] pairs — avoids
// any runtime branching or OPP[] lookup inside the hot loop.
// ---------------------------------------------------------
const TRT_PAIRS = new Uint8Array([1, 3,  2, 4,  5, 7,  6, 8]);
// 4 pairs × 2 indices = 8 entries

// Scratch buffer for feq — reused per node, no per-frame alloc.
const feqScratch = new Float32Array(Q);

// ---------------------------------------------------------
// collideFused(solid, wallAdjacent, wallDist, u0)
// ---------------------------------------------------------
export function collideFused(solid, wallAdjacent, wallDist, u0) {
  const tauBase    = 3.0 * NU0 + 0.5;
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
      reequilibrate(n, u0, 0.0, 1.0);
      continue;
    }

    let u = jx / r;
    let v = jy / r;

    if (!Number.isFinite(u) || !Number.isFinite(v)) {
      reequilibrate(n, u0, 0.0, 1.0);
      continue;
    }

    const spd2 = u * u + v * v;
    if (spd2 > U_CAP * U_CAP) {
      const sc = U_CAP / Math.sqrt(spd2);
      u *= sc; v *= sc;
    }

    rho[n] = r;
    ux[n]  = u;
    uy[n]  = v;

    // =====================================================
    // Direction loop 2: feq + non-equilibrium stress
    // Sij = -3/(2*rho) * sum_q cq_i*cq_j*fneq_q
    // =====================================================
    const v215         = 1.5 * (u * u + v * v);
    const strainFactor = STRAIN_PRE / r; // -3/(2*rho)

    let Pxx = 0.0, Pyy = 0.0, Pxy = 0.0;

    for (let i = 0; i < Q; i++) {
      const cx  = CX[i], cy = CY[i];
      const cu  = cx * u + cy * v;
      const eq  = W[i] * r * (1.0 + 3.0 * cu + 4.5 * cu * cu - v215);
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
    const Smag = Math.sqrt(2.0 * (Sxx*Sxx + Syy*Syy + 2.0*Sxy*Sxy));
    let   nuT  = CS_DELTA_SQ * Smag;

    // =====================================================
    // Scalar: wall model (log-law + Van Driest)
    // Only on wall-adjacent cells (~500 nodes for cylinder).
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
        if (Math.abs(next - uStar) < WALL_TOL) { uStar = next; break; }
        uStar = next;
      }

      const yPlus = (y * uStar) / NU0;
      const fDamp = 1.0 - Math.exp(-yPlus / A_PLUS);
      nuT = nuT * fDamp;
    }

    // =====================================================
    // Scalar: TRT relaxation rates
    //
    // omega+ (even, symmetric) — controls physical viscosity:
    //   omega+ = 1 / (3*nu_tot + 0.5)
    //
    // omega- (odd, asymmetric) — controls ghost/numerical modes:
    //   From Lambda = (tau+ - 0.5)*(tau- - 0.5) = 3/16
    //   => tau- - 0.5 = Lambda / (tau+ - 0.5)
    //   => tau-       = Lambda / (tau+ - 0.5) + 0.5
    //   => omega-     = 1 / tau-
    //
    // tau+ - 0.5 = 3*nu_tot = the "excess" above the stability
    // boundary. Flooring at a small epsilon prevents division
    // by zero if nu_tot is somehow driven to 0 by clamping.
    // =====================================================
    const nuTot   = NU0 + nuT;
    let   tauPlus = 3.0 * nuTot + 0.5;
    if (tauPlus < TAU_MIN) tauPlus = TAU_MIN;
    if (tauPlus > TAU_MAX) tauPlus = TAU_MAX;

    const omegaPlus = 1.0 / tauPlus;

    // tau- from Magic Parameter: tau- = Lambda/(tau+ - 0.5) + 0.5
    const tauMinus  = LAMBDA / (tauPlus - 0.5) + 0.5;
    const omegaMinus = 1.0 / tauMinus;

    // =====================================================
    // Direction loop 3: TRT relax + recombine
    //
    // For each opposite pair (i, ī):
    //   fi+     = (fi + fī) / 2          even part of f
    //   fi-     = (fi - fī) / 2          odd part of f
    //   feqi+   = (feqi + feqī) / 2      even part of feq
    //   feqi-   = (feqi - feqī) / 2      odd part of feq
    //
    //   fi+,post = fi+  - omega+ * (fi+  - feqi+)
    //   fi-,post = fi-  - omega- * (fi-  - feqi-)
    //
    //   fipost  = fi+,post + fi-,post
    //   fīpost  = fi+,post - fi-,post
    //
    // Rest direction (i=0, ī=0) is symmetric only:
    //   f0post  = f0 - omega+ * (f0 - feq0)
    //
    // Loop over 4 pairs (8 entries in TRT_PAIRS), then handle
    // rest direction separately. Total: 8+1 = 9 populations
    // updated, same as BGK, with 2 extra multiplies per pair
    // (omegaMinus) and a handful of extra adds for the split.
    // =====================================================

    // Rest direction (index 0) — symmetric, no opposite
    f[base] -= omegaPlus * (f[base] - feqScratch[0]);

    // 4 opposite pairs
    for (let p = 0; p < 8; p += 2) {
      const i  = TRT_PAIRS[p];
      const ib = TRT_PAIRS[p + 1]; // opposite of i

      const fi   = f[base + i ];
      const fib  = f[base + ib];
      const eqi  = feqScratch[i ];
      const eqib = feqScratch[ib];

      // even/odd split
      const fPlus  = 0.5 * (fi  + fib);  // fi+
      const fMinus = 0.5 * (fi  - fib);  // fi-
      const eqPlus = 0.5 * (eqi + eqib); // feqi+
      const eqMinus= 0.5 * (eqi - eqib); // feqi-

      // TRT relax
      const fPlusPost  = fPlus  - omegaPlus  * (fPlus  - eqPlus );
      const fMinusPost = fMinus - omegaMinus * (fMinus - eqMinus);

      // Recombine into i and ī
      f[base + i ] = fPlusPost + fMinusPost;
      f[base + ib] = fPlusPost - fMinusPost;
    }
  }
}

// ---------------------------------------------------------
// reequilibrate(n, u, v, r)
// ---------------------------------------------------------
function reequilibrate(n, u, v, r) {
  rho[n]      = r;
  ux[n]       = u;
  uy[n]       = v;
  const base  = n * Q;
  const v215  = 1.5 * (u * u + v * v);
  for (let i = 0; i < Q; i++) {
    const cu    = CX[i] * u + CY[i] * v;
    f[base + i] = W[i] * r * (1.0 + 3.0 * cu + 4.5 * cu * cu - v215);
  }
}