// src/convergenceMonitor.js
// ============================================================
// Lightweight, scalable convergence monitor for per-sim tracking
// with optional popup plotting on an overlay canvas.
//
// Current metrics:
// - meanSpeed     : mean |u| over fluid nodes
// - rmsDeltaU     : RMS change in velocity since previous sampled frame
// - meanRhoDev    : mean |rho - 1| over fluid nodes
//
// Design goals:
// - Minimal solver intrusion
// - Read-only over sim fields
// - Easy to add/remove metrics later
// - Internal sampling cadence (default every 3 frames)
// - Built-in lightweight themed plotting for popup overlay
// - Normalized plotting for mixed-magnitude metrics
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
};

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

  const enabledMetrics = options.enabledMetrics || [
    "meanSpeed",
    "rmsDeltaU",
    "meanRhoDev",
  ];

  const state = {
    hasPrevSample: false,
    prevUx: new Float32Array(sim.N),
    prevUy: new Float32Array(sim.N),
    lastSampledFrame: -1,
    visible: false,
    canvas: null,
    ctx: null,
    dpr: 1,
  };

  const history = {
    frame: [],
    step: [],
  };

  for (let i = 0; i < enabledMetrics.length; i++) {
    history[enabledMetrics[i]] = [];
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

    for (let i = 0; i < enabledMetrics.length; i++) {
      const key = enabledMetrics[i];
      const metric = METRICS[key];
      if (!metric) continue;
      history[key].push(metric.compute(sim, ctx, state));
    }

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

    // Title
    ctx.fillStyle = "#c8d8f0";
    ctx.font = '500 11px "JetBrains Mono", monospace';
    ctx.fillText("Convergence", 12, 16);

    // Plot region
    const left = 34;
    const right = w - 10;
    const top = 26;
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

    // Axis labels
    // ctx.fillStyle = "#9fb4d3";
    // ctx.font = '500 9px "JetBrains Mono", monospace';
    // ctx.fillText("norm", 4, top + 8);

    // ctx.save();
    // ctx.translate(w * 0.5, h - 8);
    // ctx.textAlign = "center";
    // ctx.fillText("iteration history", 0, 0);
    // ctx.restore();

    // Y tick labels
    ctx.fillStyle = "rgba(255,255,255,0.45)";
    ctx.font = '500 8px "JetBrains Mono", monospace';
    for (let i = 0; i <= 4; i++) {
      const y = bottom - (plotH * i) / 4;
      const label = (i / 4).toFixed(2);
      ctx.fillText(label, 8, y + 3);
    }

    // Draw normalized lines
    for (let i = 0; i < enabledMetrics.length; i++) {
      const key = enabledMetrics[i];
      const metric = METRICS[key];
      const arr = history[key];
      if (!arr || arr.length < 2) continue;

      let localMax = 1e-12;
      for (let j = 0; j < arr.length; j++) {
        if (arr[j] > localMax) localMax = arr[j];
      }

      ctx.strokeStyle = metric.color;
      ctx.lineWidth = 1.6;
      ctx.beginPath();

      for (let j = 0; j < arr.length; j++) {
        const x = left + (plotW * j) / Math.max(arr.length - 1, 1);
        const yn = localMax > 0 ? arr[j] / localMax : 0;
        const y = bottom - plotH * yn;

        if (j === 0) ctx.moveTo(x, y);
        else ctx.lineTo(x, y);
      }

      ctx.stroke();
    }

    // Legend
   let ly = h - 10;
    const legendWidths = enabledMetrics.map(key => ctx.measureText(METRICS[key].label).width + 28);
    const totalLegendWidth = legendWidths.reduce((a, b) => a + b, 0);
    let lx = (w - totalLegendWidth) * 0.5;

    ctx.font = '500 10px "JetBrains Mono", monospace';

    for (let i = 0; i < enabledMetrics.length; i++) {
      const key = enabledMetrics[i];
      const metric = METRICS[key];

      ctx.fillStyle = metric.color;
      ctx.fillRect(lx, ly - 7, 8, 2);

      ctx.fillStyle = "#9fb4d3";
      ctx.fillText(metric.label, lx + 12, ly);

      lx += ctx.measureText(metric.label).width + 28;
    }

    ctx.restore();
  }

  return {
    history,
    enabledMetrics,
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
  };
}

// ---------------------------------------------------------
// buildSharedContext(sim, state)
// ---------------------------------------------------------
function buildSharedContext(sim, state) {
  const { N, ux, uy, rho, solid } = sim;

  let fluidCount = 0;
  let sumSpeed = 0.0;
  let sumRhoDev = 0.0;
  let sumDeltaUSq = 0.0;

  for (let n = 0; n < N; n++) {
    if (solid[n]) continue;

    const u = ux[n];
    const v = uy[n];
    const r = rho[n];

    const spd = Math.sqrt(u * u + v * v);

    fluidCount++;
    sumSpeed += spd;
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



