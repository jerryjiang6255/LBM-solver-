// src/boundary.js  (V2 — fused single-pass edges)

// ============================================================

// All boundary work is done in one loop per edge, never over

// the full grid. Sponge, outlet copy, inlet reset, and wall

// reset are all fused into the traversals that already visit

// those cells — no separate full-grid pass anywhere.

// ============================================================



import { NX, NY, N, Q, f, rho, ux, uy, idx, feqAt } from "./lbm.js";





// ---------------------------------------------------------

// Sponge weights — precomputed once at module load.

// Only the rightmost SPONGE_WIDTH columns need them.

// ---------------------------------------------------------

const SPONGE_WIDTH    = 28;

const SPONGE_STRENGTH = 0.6;

const spongeW = new Float32Array(NX);



for (let x = NX - 1 - SPONGE_WIDTH; x < NX; x++) {

  const t = (x - (NX - 1 - SPONGE_WIDTH)) / SPONGE_WIDTH;

  spongeW[x] = SPONGE_STRENGTH * t * t;

}



// ---------------------------------------------------------

// applyBoundaryConditions(u0, solid)

// ---------------------------------------------------------

// Four edge passes, each visiting only the cells that live on

// that boundary. No pass over the interior. Order matches the

// original: sponge + outlet are fused into the right-column

// pass; inlet is the left-column pass; top/bottom are their

// own row passes. Corner overlap is harmless — the last writer

// wins, which is fine since all corners get a hard eq reset.

// ---------------------------------------------------------

export function applyBoundaryConditions(u0, solid) {



  // ---- left column: inlet hard Dirichlet reset ----

  for (let y = 0; y < NY; y++) {

    const n = idx(0, y);

    for (let i = 0; i < 9; i++) {

      f[n * Q + i] = feqAt(i, 1.0, u0, 0.0);

    }

  }



  // ---- right columns: sponge then outlet copy, fused ----

  // Sponge runs over the last SPONGE_WIDTH columns.

  // Outlet copy (zero-gradient) runs only on the last column.

  const xSpongeStart = NX - 1 - SPONGE_WIDTH;



  for (let x = xSpongeStart; x < NX; x++) {

    const s = spongeW[x];

    for (let y = 0; y < NY; y++) {

      const n = idx(x, y);

      if (solid[n]) continue;

      const r = rho[n], u = ux[n], v = uy[n];

      for (let i = 0; i < 9; i++) {

        const off = n * Q + i;

        f[off] = (1 - s) * f[off] + s * feqAt(i, r, u, v);

      }

    }

  }



  // outlet zero-gradient copy: last column gets second-to-last

  for (let y = 0; y < NY; y++) {

    const nSrc = idx(NX - 2, y);

    const nDst = idx(NX - 1, y);

    for (let i = 0; i < 9; i++) {

      f[nDst * Q + i] = f[nSrc * Q + i];

    }

  }



  // ---- top and bottom rows: far-field Dirichlet reset ----

  for (let x = 0; x < NX; x++) {

    const nb = idx(x, 0);

    const nt = idx(x, NY - 1);

    for (let i = 0; i < 9; i++) {

      const eq = feqAt(i, 1.0, u0, 0.0);

      f[nb * Q + i] = eq;

      f[nt * Q + i] = eq;

    }

  }

}