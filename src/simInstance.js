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
// - Version 6: physics config added per sim for solver toggles
// - Change 1: pulsing inlet flow state added per sim
// - Added: per-sim force / moment storage for Ladd momentum exchange
// - Added: per-sim coefficient storage and scaling references
// ============================================================

import { Q, RHO0, U0, CX, CY, W } from "./lbm.js";

// ---------------------------------------------------------
// feqAt(i, r, u, v)
// ---------------------------------------------------------
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

  // Reset force outputs on initialization
  if (sim.force) {
    sim.force.fx = 0.0;
    sim.force.fy = 0.0;
    sim.force.mz = 0.0;
  }

  // Reset coefficient outputs on initialization
  if (sim.coeff) {
    sim.coeff.cd = 0.0;
    sim.coeff.cl = 0.0;
    sim.coeff.cm = 0.0;
  }
}

// ---------------------------------------------------------
// getInletVelocity(sim)
// ---------------------------------------------------------
// Returns the instantaneous inlet velocity.
// Pulsation uses multiplicative sinusoidal forcing:
//   u(t) = u0 * (1 + amp * sin(2π f step))
// and is clamped to remain non-negative.
export function getInletVelocity(sim) {
  const base = sim.params.u0;
  const flow = sim.flow;

  if (!flow || !flow.pulsationOn) return base;

  const step = sim.runtime?.stepCount || 0;
  const amp  = flow.pulsationAmp;
  const freq = flow.pulsationFreq;

  const u = base * (1 + amp * Math.sin(2 * Math.PI * freq * step));
  return Math.max(0.0, u);
}

// ---------------------------------------------------------
// createSimInstance(options)
// ---------------------------------------------------------
export function createSimInstance(options = {}) {
  const {
    // -----------------------------------------------------
    // Grid dimensions
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

    // -----------------------------------------------------
    // Version 6 physics defaults
    // -----------------------------------------------------
    collisionModel = "trt",
    lesModel = "smagorinsky",
    wallModel = true,

    // -----------------------------------------------------
    // Change 1: inlet pulsation defaults
    // -----------------------------------------------------
    pulsationOn = false,
    pulsationAmp = 0.10,
    pulsationFreq = 0.01,
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
    U0,

    // -----------------------------------------------------
    // Per-sim Reynolds/viscosity state
    // -----------------------------------------------------
    L: Math.max(4, Math.floor(NY / 8)) * 2,
    RE: reynolds,
    NU0: 0,
    TAU0: 0,
    OMEGA0: 0,

    // -----------------------------------------------------
    // Core field arrays
    // -----------------------------------------------------
    f: new Float32Array(Q * N),
    rho: new Float32Array(N),
    ux: new Float32Array(N),
    uy: new Float32Array(N),

    // -----------------------------------------------------
    // Geometry/state arrays
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
    // Scratch buffers
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
    // Physics state (Version 6)
    // -----------------------------------------------------
    physics: {
      collisionModel,
      lesModel,
      wallModel,
    },

    // -----------------------------------------------------
    // Flow state
    // Change 1: pulsing inlet controls
    // -----------------------------------------------------
    flow: {
      pulsationOn,
      pulsationAmp,
      pulsationFreq,
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
      stepCount: 0,
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
    // Force / moment storage
    // Standard Ladd-style momentum exchange output:
    // - fx: drag-like force (streamwise)
    // - fy: lift-like force (cross-stream)
    // - mz: moment about forceRef
    //
    // forceRef is the chosen reference point for moment:
    // - cylinder: usually center
    // - airfoil: usually quarter-chord or center
    // -----------------------------------------------------
    force: {
      fx: 0.0,
      fy: 0.0,
      mz: 0.0,
    },

    forceRef: {
      x: 0.0,
      y: 0.0,
    },

    // -----------------------------------------------------
    // Coefficient storage
    // - cd: drag coefficient
    // - cl: lift coefficient
    // - cm: moment coefficient
    //
    // forceScale and momentScale are set elsewhere
    // (typically when geometry/reference quantities are applied).
    // -----------------------------------------------------
    coeff: {
      cd: 0.0,
      cl: 0.0,
      cm: 0.0,
    },

    forceScale: 1.0,
    momentScale: 1.0,

    // -----------------------------------------------------
    // Rendering-related handles
    // -----------------------------------------------------
    renderer: null,
    flowCanvas: null,
    bodyCanvas: null,
    wallCanvas: null,
  };

  recomputeBaseViscosity(sim);
  return sim;
}


