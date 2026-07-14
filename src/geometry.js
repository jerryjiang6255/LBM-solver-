// src/geometry.js  (V4)

import { NX, NY, N, Q, CX, CY, idx, setCharacteristicLength } from "./lbm.js";

// BFL q clamping
const Q_MIN = 0.01;
const Q_MAX = 1.0;

// ---------------------------------------------------------
// SDF FUNCTIONS
// ---------------------------------------------------------
function sdfCircle(cx, cy, radius) {
  return (x, y) => Math.sqrt((x - cx) ** 2 + (y - cy) ** 2) - radius;
}

// ---------------------------------------------------------
// Empty domain — no obstacle, pure free-stream
// ---------------------------------------------------------
function sdfEmpty() {
  return (x, y) => 1.0;
}

// ---------------------------------------------------------
// NACA 4-digit thickness distribution
// ---------------------------------------------------------
function sdfNACA0012(cx, cy, chord, alpha) {
  const cos_a = Math.cos(alpha);
  const sin_a = Math.sin(alpha);

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
// Slot nozzle extruding from the inlet wall
// ---------------------------------------------------------
function sdfNozzle(cy, diameter, length) {
  const halfD = diameter / 2;
  const wallThickness = 3; // lattice units — thickness of nozzle walls

  return (x, y) => {
    // Nozzle walls only exist from x=0 to x=length
    if (x > length) return 1.0;

    // Distance from centerline
    const distFromCenter = Math.abs(y - cy);

    // Inside the nozzle opening — fluid
    if (distFromCenter <= halfD) return 1.0;

    // Inside the wall thickness — solid
    if (distFromCenter <= halfD + wallThickness) return -1.0;

    // Outside the walls — fluid (coflow region)
    return 1.0;
  };
}

// ---------------------------------------------------------
// Matrix of square blocks
// ---------------------------------------------------------
function sdfExplicitBlocks(blocks, blockSize) {
  const half = blockSize / 2;

  return (x, y) => {
    let minDist = 1.0;
    for (let k = 0; k < blocks.length; k++) {
      const lx = Math.abs(x - blocks[k][0]) - half;
      const ly = Math.abs(y - blocks[k][1]) - half;
      // Box SDF: negative inside, positive outside
      const dist = Math.max(lx, ly);
      if (dist < minDist) minDist = dist;
    }
    return minDist;
  };
}

// ---------------------------------------------------------
// buildGeometryFromSDF(sdfFunc)
// ---------------------------------------------------------
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
export function buildAirfoil(cx, cy, chord, alpha) {
  setCharacteristicLength(chord);
  return buildGeometryFromSDF(sdfNACA0012(cx, cy, chord, alpha));
}

export function buildEmpty() {
  setCharacteristicLength(NY); // use domain height as reference
  return buildGeometryFromSDF(sdfEmpty());
}

// ---------------------------------------------------------
// defaultBlockMatrix()
// ---------------------------------------------------------
export function defaultBlockMatrix() {
  const blockSize = 4;

  const colXs = [
    Math.floor(NX * 0.30),
    Math.floor(NX * 0.43),
    Math.floor(NX * 0.56),
    Math.floor(NX * 0.69),
  ];

  const rows4 = [
    Math.floor(NY * 0.18),
    Math.floor(NY * 0.38),
    Math.floor(NY * 0.62),
    Math.floor(NY * 0.82),
  ];

  const rows3 = [
    Math.floor(NY * 0.28),
    Math.floor(NY * 0.50),
    Math.floor(NY * 0.72),
  ];

  const blocks = [];
  for (const x of [colXs[0], colXs[2]]) {
    for (const y of rows4) blocks.push([x, y]);
  }
  for (const x of [colXs[1], colXs[3]]) {
    for (const y of rows3) blocks.push([x, y]);
  }

  setCharacteristicLength(blockSize * 4);
  const geometry = buildGeometryFromSDF(sdfExplicitBlocks(blocks, blockSize));

  // Return blocks and blockSize alongside geometry so renderer can draw them
  return { geometry, blocks, blockSize };
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

export function buildNozzle(cy, diameter) {
  const length = Math.floor(NX * 0.15); // nozzle extends 15% into domain
  setCharacteristicLength(diameter);
  const geometry = buildGeometryFromSDF(sdfNozzle(cy, diameter, length));
  return { geometry, cy, diameter, length };
}

export function defaultNozzle() {
  const cy       = Math.floor(NY / 2);
  const diameter = Math.floor(NY / 4); // 25% of domain height
  return buildNozzle(cy, diameter);
}