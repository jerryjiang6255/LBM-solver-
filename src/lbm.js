// src/lbm.js  (V4)
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

