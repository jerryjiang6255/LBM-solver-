// src/geometry.js  (V3 — 2-pass SDF geometry)
// ============================================================
// Pass 1: N nodes, scalar per node — no neighbor reads
//   sdf[n], solid[n], wallDist[n]
//
// Pass 2: N nodes × Q directions — neighbor reads
//   wallAdjacent[n], linkMask[i*N+n], linkQ[i*N+n], solidList
//
// Total: N + N*Q ops at startup. Zero per-frame geometry cost.
// ============================================================

import { NX, NY, N, Q, CX, CY, idx, setCharacteristicLength } from "./lbm.js";

// BFL q clamping
const Q_MIN = 0.01;
const Q_MAX = 1.0;

// ---------------------------------------------------------
// SDF FUNCTIONS
// Only this needs to change when adding new shapes in V4+.
// ---------------------------------------------------------
function sdfCircle(cx, cy, radius) {
  return (x, y) => Math.sqrt((x - cx) ** 2 + (y - cy) ** 2) - radius;
}
// Add this alongside sdfCircle at the top of geometry.js

// ---------------------------------------------------------
// NACA 4-digit thickness distribution
// Returns signed distance approximation for a NACA 00xx airfoil
// at angle of attack `alpha` (radians), chord length `chord`,
// centered at (cx, cy).
// ---------------------------------------------------------
function sdfNACA0012(cx, cy, chord, alpha) {
  const cos_a = Math.cos(-alpha);
  const sin_a = Math.sin(-alpha);

  return (x, y) => {
    // Translate and rotate into airfoil-local frame
    const dx =  (x - cx) * cos_a + (y - cy) * sin_a;
    const dy = -(x - cx) * sin_a + (y - cy) * cos_a;

    // Normalise to [0,1] chord
    const xn = dx / chord + 0.5;   // 0 = LE, 1 = TE
    const yn = dy / chord;

    // Outside chord extent — treat as far fluid
    if (xn < 0 || xn > 1) return Math.abs(dx) + Math.abs(dy);

    // NACA 00xx half-thickness formula (xx = 12)
    const t = 0.12;
    const yt = (t / 0.2) * chord * (
       0.2969 * Math.sqrt(xn)
      - 0.1260 * xn
      - 0.3516 * xn * xn
      + 0.2843 * xn * xn * xn
      - 0.1015 * xn * xn * xn * xn
    );

    // Signed distance: positive outside, negative inside
    return Math.abs(yn * chord) - yt;
  };
}

// ---------------------------------------------------------
// buildGeometryFromSDF(sdfFunc)
// ---------------------------------------------------------
// Core builder — takes any SDF function, returns all geometry
// arrays needed by solver + renderer. Called by buildCylinder
// and rebuildGeometry (V4 live-redraw hook).
//
// Pass 1 — N nodes, 1 sdfFunc call + scalar math per node:
//   sdf[n]      : signed distance (+ fluid, - solid)
//   solid[n]    : 1 if sdf <= 0
//   wallDist[n] : sdf value for fluid nodes, 0 for solid
//
// Pass 2 — N nodes × Q directions, neighbor reads:
//   wallAdjacent[n]  : 1 if any Q neighbor is solid
//   linkMask[i*N+n]  : 1 if direction i from n hits solid
//   linkQ[i*N+n]     : BFL q from SDF linear interpolation
//   solidList        : compact array of surface solid node indices
//
// The two passes can't be merged because Pass 2 reads
// solid[neighbor] — which requires solid[] to be fully
// written for all nodes before any neighbor check runs.
// That's the fundamental data dependency that makes 2 passes
// the minimum rather than 1.
function buildGeometryFromSDF(sdfFunc) {
  const sdf          = new Float32Array(N);
  const solid        = new Uint8Array(N);
  const wallDist     = new Float32Array(N);
  const wallAdjacent = new Uint8Array(N);
  const linkMask     = new Uint8Array(Q * N);
  const linkQ        = new Float32Array(Q * N);
  const solidTmp     = [];
  const fluidTmp     = [];   // ← moved inside the function

  // Pass 1
  for (let y = 0; y < NY; y++) {
    for (let x = 0; x < NX; x++) {
      const n     = y * NX + x;
      const d     = sdfFunc(x, y);
      sdf[n]      = d;
      solid[n]    = d <= 0 ? 1 : 0;
      wallDist[n] = d > 0 ? d : 0.0;
    }
  }

  // Pass 2
  for (let y = 0; y < NY; y++) {
    for (let x = 0; x < NX; x++) {
      const n       = y * NX + x;
      const isSolid = solid[n];

      if (isSolid) {
        if (
          x > 0 && x < NX - 1 &&
          y > 0 && y < NY - 1 &&
          (!solid[n - 1] || !solid[n + 1] || !solid[n - NX] || !solid[n + NX])
        ) {
          solidTmp.push(n);
        }
        continue;
      }

      // Fluid node
      const sdfF = sdf[n];
      let isAdj  = 0;

      for (let i = 1; i < Q; i++) {
        const nx = x + CX[i];
        const ny = y + CY[i];
        if (nx < 0 || nx >= NX || ny < 0 || ny >= NY) continue;
        const ns = ny * NX + nx;
        if (!solid[ns]) continue;

        isAdj = 1;

        const sdfS  = sdf[ns];
        const denom = sdfF - sdfS;
        let q = denom > 1e-10 ? sdfF / denom : 0.5;
        if (q < Q_MIN) q = Q_MIN;
        if (q > Q_MAX) q = Q_MAX;

        const off     = n * Q + i;
        linkMask[off] = 1;
        linkQ[off]    = q;
      }

      wallAdjacent[n] = isAdj;
      if (isAdj) fluidTmp.push(n);   // ← inside the loop, after isAdj is set
    }
  }

  return {
    solid,
    sdf,
    wallDist,
    wallAdjacent,
    linkMask,
    linkQ,
    solidList:        new Uint32Array(solidTmp),
    fluidSurfaceList: new Uint32Array(fluidTmp),
  };
}

// ---------------------------------------------------------
// buildCylinder(cx, cy, radius)
// ---------------------------------------------------------
export function buildCylinder(cx, cy, radius) {
  setCharacteristicLength(2 * radius);
  return buildGeometryFromSDF(sdfCircle(cx, cy, radius));
}
// ---------------------------------------------------------
// buildAirfoil(cx, cy, chord, alpha)
// ---------------------------------------------------------
// alpha : angle of attack in RADIANS
// chord : chord length in lattice units
export function buildAirfoil(cx, cy, chord, alpha) {
  setCharacteristicLength(chord);
  return buildGeometryFromSDF(sdfNACA0012(cx, cy, chord, alpha));
}

// ---------------------------------------------------------
// defaultCylinder()
// ---------------------------------------------------------
export function defaultCylinder() {
  const cx     = Math.floor(NX / 4);
  const cy     = Math.floor(NY / 2);
  const radius = Math.max(4, Math.floor(NY / 8));
  return buildCylinder(cx, cy, radius);
}

// ---------------------------------------------------------
// rebuildGeometry(sdfFunc)
// ---------------------------------------------------------
// V4 live-redraw hook: rebuild from any SDF function.
// Call from index.html when user draws a new shape.
// Returns the same object as buildCylinder — swap it in
// and the solver picks up the new geometry on the next step.
export function rebuildGeometry(sdfFunc) {
  return buildGeometryFromSDF(sdfFunc);
}