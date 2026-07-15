// src/simInstance.js
// ============================================================
// Minimal per-simulation state container for compare mode.
//
// Purpose:
// - Replace the single global mutable sim state with one object per sim
// - Keep your current architecture mostly intact
// - Let solver modules operate on `sim` instead of shared globals
//
// Notes:
// - We still import immutable lattice constants/helpers from lbm.js
// - Mutable field arrays (f, rho, ux, uy) now live on each sim instance
// - Geometry/state/display/runtime params are also stored per sim
// - For global fidelity presets, NX/NY are now passed in at creation time
// ============================================================

import { Q, RHO0, U0, CX, CY, W } from "./lbm.js";

// ---------------------------------------------------------
// feqAt(i, r, u, v)
// ---------------------------------------------------------
// Local equilibrium population for D2Q9.
// Kept here so sim instances don't depend on global mutable arrays.
export function feqAt(i, r, u, v) {
  const cu = CX[i] * u + CY[i] * v;
  const uu = u * u + v * v;
  return W[i] * r * (1 + 3 * cu + 4.5 * cu * cu - 1.5 * uu);
}

// ---------------------------------------------------------
// recomputeBaseViscosity(sim)
// ---------------------------------------------------------
export function recomputeBaseViscosity(sim) {
  sim.NU0    = (sim.params.u0 * sim.L) / sim.RE;
  sim.TAU0   = 3 * sim.NU0 + 0.5;
  sim.OMEGA0 = 1 / sim.TAU0;
}

// ---------------------------------------------------------
// setCharacteristicLength(sim, diameter)
// ---------------------------------------------------------
export function setCharacteristicLength(sim, diameter) {
  sim.L = diameter;
  recomputeBaseViscosity(sim);
}

// ---------------------------------------------------------
// setReynolds(sim, re)
// ---------------------------------------------------------
export function setReynolds(sim, re) {
  sim.RE = re;
  recomputeBaseViscosity(sim);
}

// ---------------------------------------------------------
// initField(sim, u0, solid)
// ---------------------------------------------------------
// Initializes macroscopic fields and populations for one sim.
export function initField(sim, u0 = sim.params.u0, solid = sim.solid) {
  for (let n = 0; n < sim.N; n++) {
    sim.rho[n] = RHO0;

    if (solid && solid[n]) {
      sim.ux[n] = 0.0;
      sim.uy[n] = 0.0;
    } else {
      sim.ux[n] = u0;
      sim.uy[n] = 0.0;
    }
  }

  for (let n = 0; n < sim.N; n++) {
    const base = n * sim.Q;
    const r = sim.rho[n];
    const u = sim.ux[n];
    const v = sim.uy[n];

    for (let i = 0; i < sim.Q; i++) {
      sim.f[base + i] = feqAt(i, r, u, v);
    }
  }
}

// ---------------------------------------------------------
// createSimInstance(options)
// ---------------------------------------------------------
// Creates one independent simulation object.
// This is the core object you will pass into collision/stream/boundary.
//
// For global fidelity presets, pass NX and NY from the selected preset.
// Example:
//   createSimInstance({ id: "sim1", NX: 160, NY: 80, ... })
export function createSimInstance(options = {}) {
  const {
    // -----------------------------------------------------
    // Grid dimensions (required for fidelity switching)
    // -----------------------------------------------------
    NX = 256,
    NY = 128,

    id = "sim",
    bodyType = "airfoil",
    aoaDeg = 10,
    nozzleDiameter = 16,
    u0 = 0.1732,
    reynolds = 5000,

    displayMode = "speed",
    colorScale = 0.5,
    showVectors = false,
    showGrid = false,
    gridFine = true,
    showDist = false,

    stepsPerFrame = 10,
    running = true,

    // nozzle defaults
    nozzleCy = Math.floor(NY / 2),
    uCoflow = 0.005,
  } = options;

  const N = NX * NY;

  const sim = {
    // -----------------------------------------------------
    // Identity
    // -----------------------------------------------------
    id,

    // -----------------------------------------------------
    // Per-sim dimensions
    // -----------------------------------------------------
    NX,
    NY,
    N,
    Q,

    // -----------------------------------------------------
    // Base physical references
    // -----------------------------------------------------
    RHO0,
    U0, // reference velocity currently used for base viscosity scaling

    // -----------------------------------------------------
    // Per-sim Reynolds/viscosity state
    // -----------------------------------------------------
    L: Math.max(4, Math.floor(NY / 8)) * 2,
    RE: reynolds,
    NU0: 0,
    TAU0: 0,
    OMEGA0: 0,

    // -----------------------------------------------------
    // Core field arrays (independent per sim)
    // -----------------------------------------------------
    f: new Float32Array(Q * N),
    rho: new Float32Array(N),
    ux: new Float32Array(N),
    uy: new Float32Array(N),

    // -----------------------------------------------------
    // Geometry/state arrays (filled later by geometry build)
    // -----------------------------------------------------
    solid: null,
    solidBase: null,
    wallDist: null,
    wallAdjacent: null,
    solidList: null,
    fluidSurfaceList: null,
    linkMask: null,
    linkQ: null,

    // -----------------------------------------------------
    // Scratch buffers that should not be shared between sims
    // -----------------------------------------------------
    fSnapshot: new Float32Array(Q * N),

    // -----------------------------------------------------
    // Interaction state
    // -----------------------------------------------------
    paintedNodes: new Set(),

    // Optional future drawing state
    strokes: [],
    currentStroke: null,

    // -----------------------------------------------------
    // Parameter state
    // -----------------------------------------------------
    params: {
      u0,
      reynolds,
      bodyType,
      aoaDeg,
      nozzleDiameter,
    },

    // -----------------------------------------------------
    // Display state
    // -----------------------------------------------------
    display: {
      mode: displayMode,
      colorScale,
      showVectors,
      showGrid,
      gridFine,
      showDist,
    },

    // -----------------------------------------------------
    // Runtime state
    // -----------------------------------------------------
    runtime: {
      running,
      frameCount: 0,
      stepsPerFrame,
    },

    // -----------------------------------------------------
    // Boundary/nozzle config
    // -----------------------------------------------------
    nozzle: {
      active: false,
      cy: nozzleCy,
      halfD: nozzleDiameter / 2,
      uCoflow,
    },

    // -----------------------------------------------------
    // Rendering-related handles (optional, assigned externally)
    // -----------------------------------------------------
    renderer: null,
    flowCanvas: null,
    bodyCanvas: null,
    wallCanvas: null,
  };

  recomputeBaseViscosity(sim);
  return sim;
}


