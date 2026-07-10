// src/streaming.js  (V2 — 4-pass in-place streaming)

// ============================================================

// Replaces the 9×N scratch-buffer approach with 4 ordered

// directional passes, identical in structure to the aero-tunnel

// reference's stream(). No fNew buffer needed — each pass sweeps

// the grid in the direction that guarantees the source cell is

// read before it could be overwritten by the same pass.

//

// D2Q9 direction index convention (must match lbm.js):

//   0=rest  1=E(+1,0)  2=N(0,+1)  3=W(-1,0)  4=S(0,-1)

//   5=NE(+1,+1)  6=NW(-1,+1)  7=SW(-1,-1)  8=SE(+1,-1)

//

// Pass groupings — paired by shared safe sweep direction:

//   Pass 1  N(2) + NW(6)  : sweep top→bottom, left→right

//   Pass 2  E(1) + NE(5)  : sweep top→bottom, right→left

//   Pass 3  S(4) + SE(8)  : sweep bottom→top, right→left

//   Pass 4  W(3) + SW(7)  : sweep bottom→top, left→right

//   Rest(0) never moves   : no pass needed

// ============================================================



import { NX, NY, N, f, Q } from "./lbm.js";



// Direction-index offsets into the f array.

// Each slice f[i*N ... i*N+N-1] holds all nodes for direction i.

const I_E  = 1, I_N  = 2, I_W  = 3, I_S  = 4;

const I_NE = 5, I_NW = 6, I_SW = 7, I_SE = 8;



// Precomputed base offsets — avoids i*N multiplication inside loops.

const OFF_E  = I_E  * N, OFF_N  = I_N  * N, OFF_W  = I_W  * N, OFF_S  = I_S  * N;

const OFF_NE = I_NE * N, OFF_NW = I_NW * N, OFF_SW = I_SW * N, OFF_SE = I_SE * N;



export function stream(solid, solidList) {

    // Pass 1 — N(i=2) and NW(i=6), sweep top→bottom, left→right

    for (let y = NY - 2; y >= 1; y--) {

    const rowDst = y * NX;

    const rowSrc = (y - 1) * NX;

    for (let x = 1; x < NX - 1; x++) {

        f[(rowDst + x) * Q + 2] = f[(rowSrc + x    ) * Q + 2];  // N

        f[(rowDst + x) * Q + 6] = f[(rowSrc + x + 1) * Q + 6];  // NW

    }

    }



    // Pass 2 — E(i=1) and NE(i=5), sweep top→bottom, right→left

    for (let y = NY - 2; y >= 1; y--) {

    const rowDst = y * NX;

    const rowSrc = (y - 1) * NX;

    for (let x = NX - 2; x >= 1; x--) {

        f[(rowDst + x) * Q + 1] = f[(rowDst + x - 1) * Q + 1];  // E

        f[(rowDst + x) * Q + 5] = f[(rowSrc + x - 1) * Q + 5];  // NE

    }

    }



    // Pass 3 — S(i=4) and SE(i=8), sweep bottom→top, right→left

    for (let y = 1; y < NY - 1; y++) {

    const rowDst = y * NX;

    const rowSrc = (y + 1) * NX;

    for (let x = NX - 2; x >= 1; x--) {

        f[(rowDst + x) * Q + 4] = f[(rowSrc + x    ) * Q + 4];  // S

        f[(rowDst + x) * Q + 8] = f[(rowSrc + x - 1) * Q + 8];  // SE

    }

    }



    // Pass 4 — W(i=3) and SW(i=7), sweep bottom→top, left→right

    for (let y = 1; y < NY - 1; y++) {

    const rowDst = y * NX;

    const rowSrc = (y + 1) * NX;

    for (let x = 1; x < NX - 1; x++) {

        f[(rowDst + x) * Q + 3] = f[(rowDst + x + 1) * Q + 3];  // W

        f[(rowDst + x) * Q + 7] = f[(rowSrc + x + 1) * Q + 7];  // SW

    }

    }



  // Rest direction (index 0) never moves — no pass needed.



  bounceBack(solidList);

}



function bounceBack(solidList) {

  for (let k = 0; k < solidList.length; k++) {

    const n = solidList[k];

    swapPair(n, 1, 3);   // E  <-> W

    swapPair(n, 2, 4);   // N  <-> S

    swapPair(n, 5, 7);   // NE <-> SW

    swapPair(n, 6, 8);   // NW <-> SE

  }

}



// swapPair — n*Q+i instead of i*N+n

function swapPair(n, a, b) {

  const ia = n * Q + a, ib = n * Q + b;  // was a*N+n, b*N+n

  const tmp = f[ia]; f[ia] = f[ib]; f[ib] = tmp;

}