// src/convergenceMonitor.js
// ============================================================
// Lightweight, scalable monitor for per-sim tracking
// with optional popup plotting on an overlay canvas.
//
// Plot modes:
// - convergence : meanSpeed, rmsDeltaU, meanRhoDev (normalized)
// - forces      : Cd, Cl, Cm (absolute, signed)
//
// Notes:
// - Histories are stored together so switching modes preserves data
// - Only convergence metrics are normalized
// - Force mode uses coefficient histories only
// - Moment/Cm can be hidden for selected geometry types
// - Force histories are sourced directly from sim.coeff so the
//   monitor matches the solver's actual force/coefficient outputs
// ============================================================

const DEFAULT_SAMPLE_EVERY_FRAMES = 1;
const DEFAULT_MAX_HISTORY = 600;

// ---------------------------------------------------------
// Metric registry
// ---------------------------------------------------------
const METRICS = {
  meanSpeed: {
    label: "Mean Speed",
    color: "#38c8ff",
    compute(sim, ctx) {
      return ctx.fluidCount > 0 ? ctx.sumSpeed / ctx.fluidCount : 0.0;
    },
  },

  rmsDeltaU: {
    label: "RMS ΔU",
    color: "#ffb366",
    compute(sim, ctx, state) {
      if (!state.hasPrevSample || ctx.fluidCount === 0) return 0.0;
      return Math.sqrt(ctx.sumDeltaUSq / ctx.fluidCount);
    },
  },

  meanRhoDev: {
    label: "Mean |ρ-1|",
    color: "#c78cff",
    compute(sim, ctx) {
      return ctx.fluidCount > 0 ? ctx.sumRhoDev / ctx.fluidCount : 0.0;
    },
  },

  cd: {
    label: "Cd",
    color: "#62d8ff",
    compute(sim) {
      return sim.coeff?.cd ?? 0.0;
    },
  },

  cl: {
    label: "Cl",
    color: "#ffc27a",
    compute(sim) {
      return sim.coeff?.cl ?? 0.0;
    },
  },

  cm: {
    label: "Cm",
    color: "#d8a2ff",
    compute(sim) {
      return sim.coeff?.cm ?? 0.0;
    },
  },
};

const MODE_CONFIG = {
  convergence: ["meanSpeed", "rmsDeltaU", "meanRhoDev"],
  forces: ["cd", "cl", "cm"],
};

// Geometries for which moment / Cm are meaningful enough to show
function shouldShowMoment(sim) {
  const type = sim?.params?.bodyType;
  return (
    type === "airfoil" ||
    type === "naca2412" ||
    type === "cylinder" ||
    type === "square-cylinder"
  );
}

// ---------------------------------------------------------
// createConvergenceMonitor(sim, options?)
// ---------------------------------------------------------
export function createConvergenceMonitor(sim, options = {}) {
  const sampleEveryFrames =
    Number.isInteger(options.sampleEveryFrames) && options.sampleEveryFrames > 0
      ? options.sampleEveryFrames
      : DEFAULT_SAMPLE_EVERY_FRAMES;

  const maxHistory =
    Number.isInteger(options.maxHistory) && options.maxHistory > 0
      ? options.maxHistory
      : DEFAULT_MAX_HISTORY;

  const state = {
    hasPrevSample: false,
    prevUx: new Float32Array(sim.N),
    prevUy: new Float32Array(sim.N),
    lastSampledFrame: -1,
    visible: false,
    canvas: null,
    ctx: null,
    dpr: 1,
    mode: "convergence",
  };

  const history = {
    frame: [],
    step: [],

    meanSpeed: [],
    rmsDeltaU: [],
    meanRhoDev: [],

    cd: [],
    cl: [],
    cm: [],
  }; 
 function activeMetricKeys() {
    if (state.mode === "forces") {
      return shouldShowMoment(sim)
        ? MODE_CONFIG.forces
        : ["cd", "cl"];
    }
    return MODE_CONFIG.convergence;
  }

  function setMode(mode) {
    if (mode === "forces" || mode === "convergence") {
      state.mode = mode;
    }
  }

  function getMode() {
    return state.mode;
  }

  function maybeSample() {
    const frame = sim.runtime?.frameCount ?? 0;
    if (frame === state.lastSampledFrame) return false;
    if (frame % sampleEveryFrames !== 0) return false;

    sampleNow();
    state.lastSampledFrame = frame;
    return true;
  }

  function sampleNow() {
    const ctx = buildSharedContext(sim, state);

    history.frame.push(sim.runtime?.frameCount ?? 0);
    history.step.push(sim.runtime?.stepCount ?? 0);

    history.meanSpeed.push(METRICS.meanSpeed.compute(sim, ctx, state));
    history.rmsDeltaU.push(METRICS.rmsDeltaU.compute(sim, ctx, state));
    history.meanRhoDev.push(METRICS.meanRhoDev.compute(sim, ctx, state));

    // Use solver-computed coefficients directly
    history.cd.push(METRICS.cd.compute(sim));
    history.cl.push(METRICS.cl.compute(sim));
    history.cm.push(METRICS.cm.compute(sim));

    updatePreviousSample(sim, state);
    trimHistory(history, maxHistory);
  }

  function reset() {
    state.hasPrevSample = false;
    state.lastSampledFrame = -1;
    state.prevUx.fill(0);
    state.prevUy.fill(0);

    for (const key in history) {
      history[key].length = 0;
    }

    clearCanvas();
  }

  function attachCanvas(canvas) {
    state.canvas = canvas;
    state.ctx = canvas ? canvas.getContext("2d") : null;
    resize();
  }

  function setVisible(flag) {
    state.visible = !!flag;
    if (!state.visible) clearCanvas();
  }

  function isVisible() {
    return state.visible;
  }

  function resize() {
    if (!state.canvas) return;

    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    state.dpr = dpr;

    const cssW = state.canvas.clientWidth || 260;
    const cssH = state.canvas.clientHeight || 150;

    state.canvas.width = Math.round(cssW * dpr);
    state.canvas.height = Math.round(cssH * dpr);
  } 
 function clearCanvas() {
    if (!state.ctx || !state.canvas) return;
    state.ctx.clearRect(0, 0, state.canvas.width, state.canvas.height);
  }

  function draw() {
    if (!state.visible || !state.ctx || !state.canvas) return;

    const ctx = state.ctx;
    const W = state.canvas.width;
    const H = state.canvas.height;
    const dpr = state.dpr || 1;

    ctx.clearRect(0, 0, W, H);
    ctx.save();
    ctx.scale(dpr, dpr);

    const w = W / dpr;
    const h = H / dpr;

    // Panel background
    roundRect(ctx, 0, 0, w, h, 10);
    ctx.fillStyle = "rgba(10,16,26,0.72)";
    ctx.fill();
    ctx.strokeStyle = "rgba(80,160,255,0.18)";
    ctx.lineWidth = 1;
    ctx.stroke();

    // Plot region
    const left = 34;
    const right = w - 10;
    const top = 18;
    const bottom = h - 30;
    const plotW = right - left;
    const plotH = bottom - top;

    // Grid
    ctx.strokeStyle = "rgba(255,255,255,0.06)";
    ctx.lineWidth = 1;
    ctx.beginPath();
    for (let i = 0; i <= 4; i++) {
      const y = top + (plotH * i) / 4;
      ctx.moveTo(left, y);
      ctx.lineTo(right, y);
    }
    for (let i = 0; i <= 4; i++) {
      const x = left + (plotW * i) / 4;
      ctx.moveTo(x, top);
      ctx.lineTo(x, bottom);
    }
    ctx.stroke();

    // Axes
    ctx.strokeStyle = "rgba(255,255,255,0.16)";
    ctx.beginPath();
    ctx.moveTo(left, top);
    ctx.lineTo(left, bottom);
    ctx.lineTo(right, bottom);
    ctx.stroke();

    const keys = activeMetricKeys();

    // -----------------------------------------------------
    // Vertical scale
    // convergence -> normalized 0..1
    // forces      -> signed symmetric scale around zero
    // -----------------------------------------------------
    let globalYMax = 1.0;

    if (state.mode === "forces") {
      globalYMax = 1e-12;
      for (let i = 0; i < keys.length; i++) {
        const arr = history[keys[i]];
        if (!arr || arr.length === 0) continue;
        for (let j = 0; j < arr.length; j++) {
          const a = Math.abs(arr[j]);
          if (a > globalYMax) globalYMax = a;
        }
      }
      if (globalYMax <= 0) globalYMax = 1.0;
    }

    // Y-axis labels
    ctx.fillStyle = "rgba(255,255,255,0.45)";
    ctx.font = '500 8px "JetBrains Mono", monospace';
    for (let i = 0; i <= 4; i++) {
      const y = top + (plotH * i) / 4;
      let label;
      if (state.mode === "convergence") {
        label = (1.0 - i / 4).toFixed(2);
      } else {
        const frac = 1.0 - i / 2.0; // 1, 0.5, 0, -0.5, -1
        const value = frac * globalYMax;

        if (Math.abs(globalYMax) >= 1) {
          label = value.toFixed(1);
        } else if (Math.abs(globalYMax) >= 0.1) {
          label = value.toFixed(2);
        } else if (Math.abs(globalYMax) >= 0.01) {
          label = value.toFixed(3);
        } else {
          label = value.toFixed(4);
        }
      }
      ctx.fillText(label, 6, y + 3);
    }

    // Zero line for force coefficients
    if (state.mode === "forces") {
      const yZero = top + 0.5 * plotH;
      ctx.strokeStyle = "rgba(255,255,255,0.14)";
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(left, yZero);
      ctx.lineTo(right, yZero);
      ctx.stroke();
    }

    // Draw lines
    for (let i = 0; i < keys.length; i++) {
      const key = keys[i];
      const metric = METRICS[key];
      const arr = history[key];
      if (!arr || arr.length < 2) continue;

      let localMax = 1e-12;

      if (state.mode === "convergence") {
        for (let j = 0; j < arr.length; j++) {
          const a = Math.abs(arr[j]);
          if (a > localMax) localMax = a;
        }
        if (localMax <= 0) localMax = 1.0;
      }

      ctx.strokeStyle = metric.color;
      ctx.lineWidth = 1.6;
      ctx.beginPath();

      for (let j = 0; j < arr.length; j++) {
        const x = left + (plotW * j) / Math.max(arr.length - 1, 1);

        let y;
        if (state.mode === "convergence") {
          const yn = localMax > 0 ? Math.abs(arr[j]) / localMax : 0.0;
          y = bottom - plotH * yn;
        } else {
          const v = globalYMax > 0 ? arr[j] / globalYMax : 0.0;
          const yn = 0.5 - 0.5 * v; // +1 top, 0 middle, -1 bottom
          y = top + plotH * yn;
        }

        if (j === 0) ctx.moveTo(x, y);
        else ctx.lineTo(x, y);
      }

      ctx.stroke();
    }

    // Legend
    ctx.font = '500 10px "JetBrains Mono", monospace';
    let ly = h - 10;
    const legendWidths = keys.map(key => ctx.measureText(METRICS[key].label).width + 28);
    const totalLegendWidth = legendWidths.reduce((a, b) => a + b, 0);
    let lx = (w - totalLegendWidth) * 0.5;

    for (let i = 0; i < keys.length; i++) {
      const key = keys[i];
      const metric = METRICS[key];

      ctx.fillStyle = metric.color;
      ctx.fillRect(lx, ly - 7, 8, 2);

      ctx.fillStyle = "#9fb4d3";
      ctx.fillText(metric.label, lx + 12, ly);

      lx += ctx.measureText(metric.label).width + 28;
    } 

  // ---- coefficient + raw force readout — top right corner ----
    if (state.mode === "forces" && history.cd.length > 0 && history.cl.length > 0) {
      const lastCd = history.cd[history.cd.length - 1];
      const lastCl = history.cl[history.cl.length - 1];
      const lastCm = history.cm.length > 0 ? history.cm[history.cm.length - 1] : null;
      const ratio  = Math.abs(lastCd) > 1e-6 ? lastCl / lastCd : 0;

      const fx = sim.force?.fx ?? 0.0;
      const fy = sim.force?.fy ?? 0.0;
      const mz = sim.force?.mz ?? 0.0;

      ctx.font = '500 9px "JetBrains Mono", monospace';
      ctx.textAlign = 'right';

      const lines = [
        { label: "Cd", val: lastCd, color: METRICS.cd.color },
        { label: "Cl", val: lastCl, color: METRICS.cl.color },
      ];

      if (lastCm !== null && shouldShowMoment(sim)) {
        lines.push({ label: "Cm", val: lastCm, color: METRICS.cm.color });
      }

      lines.push(
        { label: "Cl/Cd", val: ratio, color: "#ffffff" },
        { label: "Fx", val: fx, color: "#7fdfff" },
        { label: "Fy", val: fy, color: "#ffd08a" },
        { label: "Mz", val: mz, color: "#dcb0ff" },
      );

      let ry = top + 10;
      for (const line of lines) {
        ctx.fillStyle = line.color;
        ctx.fillText(`${line.label}: ${line.val.toFixed(4)}`, right - 2, ry);
        ry += 12;
      }

      ctx.textAlign = "left"; // reset
    }



    ctx.restore();
  }

  return {
    history,
    sampleEveryFrames,
    maxHistory,
    maybeSample,
    sampleNow,
    reset,
    attachCanvas,
    setVisible,
    isVisible,
    resize,
    draw,
    setMode,
    getMode,
    activeMetricKeys,
  };
}

// ---------------------------------------------------------
// buildSharedContext(sim, state)
// ---------------------------------------------------------
function buildSharedContext(sim, state) {
  const { N, ux, uy, rho, solid } = sim;

  let fluidCount  = 0;
  let sumSpeed    = 0.0;
  let sumRhoDev   = 0.0;
  let sumDeltaUSq = 0.0;

  for (let n = 0; n < N; n++) {
    if (solid[n]) continue;

    const u = ux[n];
    const v = uy[n];
    const r = rho[n];
    const spd = Math.sqrt(u * u + v * v);

    fluidCount++;
    sumSpeed  += spd;
    sumRhoDev += Math.abs(r - 1.0);

    if (state.hasPrevSample) {
      const du = u - state.prevUx[n];
      const dv = v - state.prevUy[n];
      sumDeltaUSq += du * du + dv * dv;
    }
  }

  return {
    fluidCount,
    sumSpeed,
    sumRhoDev,
    sumDeltaUSq,
  };
}

// ---------------------------------------------------------
// updatePreviousSample(sim, state)
// ---------------------------------------------------------
function updatePreviousSample(sim, state) {
  state.prevUx.set(sim.ux);
  state.prevUy.set(sim.uy);
  state.hasPrevSample = true;
}

// ---------------------------------------------------------
// trimHistory(history, maxHistory)
// ---------------------------------------------------------
function trimHistory(history, maxHistory) {
  const over = history.frame.length - maxHistory;
  if (over <= 0) return;

  for (const key in history) {
    history[key].splice(0, over);
  }
}

// ---------------------------------------------------------
// roundRect(ctx, x, y, w, h, r)
// ---------------------------------------------------------
function roundRect(ctx, x, y, w, h, r) {
  const rr = Math.min(r, 0.5 * w, 0.5 * h);
  ctx.beginPath();
  ctx.moveTo(x + rr, y);
  ctx.arcTo(x + w, y, x + w, y + h, rr);
  ctx.arcTo(x + w, y + h, x, y + h, rr);
  ctx.arcTo(x, y + h, x, y, rr);
  ctx.arcTo(x, y, x + w, y, rr);
  ctx.closePath();
} 
