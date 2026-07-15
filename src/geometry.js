// src/geometry.js
// ============================================================
// Compare-mode ready geometry builders.
// Updated for global fidelity presets: geometry no longer imports
// fixed NX/NY/N from lbm.js. Dimensions are passed in explicitly.
// ============================================================

import { Q, CX, CY } from "./lbm.js";

// BFL q clamping
const Q_MIN = 0.01;
const Q_MAX = 1.0;

// ---------------------------------------------------------
// SDF FUNCTIONS
// ---------------------------------------------------------
function sdfCircle(cx, cy, radius) {
  return (x, y) => Math.sqrt((x - cx) ** 2 + (y - cy) ** 2) - radius;
}

function sdfEmpty() {
  return () => 1.0;
}

function sdfNACA0012(cx, cy, chord, alpha) {
  const cos_a = Math.cos(alpha);
  const sin_a = Math.sin(alpha);

  return (x, y) => {
    const dx =  (x - cx) * cos_a + (y - cy) * sin_a;
    const dy = -(x - cx) * sin_a + (y - cy) * cos_a;

    const xn = dx / chord + 0.5;
    const yn = dy / chord;

    if (xn < 0 || xn > 1) return Math.abs(dx) + Math.abs(dy);

    const t = 0.12;
    const yt = (t / 0.2) * chord * (
       0.2969 * Math.sqrt(xn)
      - 0.1260 * xn
      - 0.3516 * xn * xn
      + 0.2843 * xn * xn * xn
      - 0.1015 * xn * xn * xn * xn
    );

    return Math.abs(yn * chord) - yt;
  };
}

function sdfNozzle(cy, diameter, length) {
  const halfD = diameter / 2;
  const wallThickness = 3;

  return (x, y) => {
    if (x > length) return 1.0;

    const distFromCenter = Math.abs(y - cy);

    if (distFromCenter <= halfD) return 1.0;
    if (distFromCenter <= halfD + wallThickness) return -1.0;

    return 1.0;
  };
}

function sdfExplicitBlocks(blocks, blockSize) {
  const half = blockSize / 2;

  return (x, y) => {
    let minDist = 1.0;
    for (let k = 0; k < blocks.length; k++) {
      const lx = Math.abs(x - blocks[k][0]) - half;
      const ly = Math.abs(y - blocks[k][1]) - half;
      const dist = Math.max(lx, ly);
      if (dist < minDist) minDist = dist;
    }
    return minDist;
  };
}

// ---------------------------------------------------------
// buildGeometryFromSDF(dims, sdfFunc)
// ---------------------------------------------------------
function buildGeometryFromSDF(dims, sdfFunc) {
  const { NX, NY } = dims;
  const N = NX * NY;

  const sdf          = new Float32Array(N);
  const solid        = new Uint8Array(N);
  const wallDist     = new Float32Array(N);
  const wallAdjacent = new Uint8Array(N);
  const linkMask     = new Uint8Array(Q * N);
  const linkQ        = new Float32Array(Q * N);
  const solidTmp     = [];
  const fluidTmp     = [];

  for (let y = 0; y < NY; y++) {
    for (let x = 0; x < NX; x++) {
      const n     = y * NX + x;
      const d     = sdfFunc(x, y);
      sdf[n]      = d;
      solid[n]    = d <= 0 ? 1 : 0;
      wallDist[n] = d > 0 ? d : 0.0;
    }
  }

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
      if (isAdj) fluidTmp.push(n);
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
// Public builders
// ---------------------------------------------------------
export function buildCylinder(dims, cx, cy, radius) {
  return {
    geometry: buildGeometryFromSDF(dims, sdfCircle(cx, cy, radius)),
    charLength: 2 * radius,
  };
}

export function buildAirfoil(dims, cx, cy, chord, alpha) {
  return {
    geometry: buildGeometryFromSDF(dims, sdfNACA0012(cx, cy, chord, alpha)),
    charLength: chord,
  };
}

export function buildEmpty(dims) {
  return {
    geometry: buildGeometryFromSDF(dims, sdfEmpty()),
    charLength: dims.NY,
  };
}

export function defaultBlockMatrix(dims) {
  const { NX, NY } = dims;
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

  return {
    geometry: buildGeometryFromSDF(dims, sdfExplicitBlocks(blocks, blockSize)),
    blocks,
    blockSize,
    charLength: blockSize * 4,
  };
}

export function defaultCylinder(dims) {
  const { NX, NY } = dims;
  const cx     = Math.floor(NX / 4);
  const cy     = Math.floor(NY / 2);
  const radius = Math.max(4, Math.floor(NY / 8));
  return buildCylinder(dims, cx, cy, radius);
}

export function buildNozzle(dims, cy, diameter) {
  const { NX } = dims;
  const length = Math.floor(NX * 0.15);

  return {
    geometry: buildGeometryFromSDF(dims, sdfNozzle(cy, diameter, length)),
    cy,
    diameter,
    length,
    charLength: diameter,
  };
}

export function defaultNozzle(dims) {
  const { NY } = dims;
  const cy       = Math.floor(NY / 2);
  const diameter = Math.floor(NY / 4);
  return buildNozzle(dims, cy, diameter);
}


