// src/renderBody.js
// ============================================================
// Body overlay renderer only.
// Compare-mode ready: pass a canvas per sim, no module-global ctx.
// Updated for global fidelity presets: renderer owns NX/NY.
// ============================================================

export class BodyRenderer {
    
  constructor(canvas, NX, NY) {
    this.canvas = canvas;
    this.ctx = canvas.getContext("2d");
    this.NX = NX;
    this.NY = NY;
  }

  draw(params) {
    const ctx = this.ctx;
    if (!ctx) return;

    const NX = this.NX;
    const NY = this.NY;

    const Wc = ctx.canvas.width;
    const Hc = ctx.canvas.height;  // ← read, don't set

    const sx = Wc / this.NX;
    const sy = Hc / this.NY;

    ctx.clearRect(0, 0, Wc, Hc);
    ctx.save();
    ctx.shadowColor = "rgba(130,190,255,0.5)";
    ctx.shadowBlur  = 14 * Math.min(sx, sy);

    if (params.type === "cylinder") {
      const rx = params.radius * sx;
      const ry = params.radius * sy;

      const grad = ctx.createRadialGradient(
        params.cx * sx, params.cy * sy, 0,
        params.cx * sx, params.cy * sy, rx
      );
      grad.addColorStop(0,   "#48506a");
      grad.addColorStop(0.5, "#1a1f2e");
      grad.addColorStop(1,   "#0a0d16");

      ctx.beginPath();
      ctx.ellipse(params.cx * sx, params.cy * sy, rx, ry, 0, 0, Math.PI * 2);
      ctx.fillStyle = grad;
      ctx.fill();
      ctx.strokeStyle = "rgba(175,198,230,0.8)";
      ctx.lineWidth = 1.5;
      ctx.stroke();

    } else if (params.type === "airfoil") {
      const pts = nacaContourExpanded(params, 0.5);

      ctx.beginPath();
      ctx.moveTo(pts[0][0] * sx, pts[0][1] * sy);
      for (let k = 1; k < pts.length; k++) {
        ctx.lineTo(pts[k][0] * sx, pts[k][1] * sy);
      }
      ctx.closePath();

      const grad = ctx.createLinearGradient(
        0, (params.cy - params.chord * 0.25) * sy,
        0, (params.cy + params.chord * 0.25) * sy
      );
      grad.addColorStop(0,    "#48506a");
      grad.addColorStop(0.45, "#1a1f2e");
      grad.addColorStop(1,    "#0a0d16");

      ctx.fillStyle   = grad;
      ctx.strokeStyle = "rgba(175,198,230,0.8)";
      ctx.lineWidth   = 1.5;
      ctx.fill();
      ctx.stroke();

    } else if (params.type === "blocks") {
      const { blocks, blockSize } = params;
      const pad = 1.0;
      const half = blockSize / 2 + pad;

      for (const [bx, by] of blocks) {
        const pw = half * 2 * sx;
        const ph = half * 2 * sy;
        const px = bx * sx - half * sx;
        const py = by * sy - half * sy;
        const r  = Math.min(pw, ph) * 0.2;

        const grad = ctx.createRadialGradient(
          bx * sx, by * sy, 0,
          bx * sx, by * sy, Math.max(pw, ph) * 0.7
        );
        grad.addColorStop(0,   "#48506a");
        grad.addColorStop(0.5, "#1a1f2e");
        grad.addColorStop(1,   "#0a0d16");

        ctx.beginPath();
        ctx.roundRect(px, py, pw, ph, r);
        ctx.fillStyle = grad;
        ctx.fill();

        ctx.strokeStyle = "rgba(175,198,230,0.8)";
        ctx.lineWidth = 1.2;
        ctx.stroke();
      }

    } else if (params.type === "nozzle") {
      const { cy, diameter, length } = params;
      const halfD = diameter / 2;
      const wallThickness = 4;

      const topY1 = (cy - halfD - wallThickness) * sy;
      const topY2 = (cy - halfD) * sy;
      const botY1 = (cy + halfD) * sy;
      const botY2 = (cy + halfD + wallThickness) * sy;
      const endX  = (length + 0.9) * sx;

      const gradTop = ctx.createLinearGradient(0, topY1, 0, topY2);
      gradTop.addColorStop(0,   "#48506a");
      gradTop.addColorStop(0.5, "#1a1f2e");
      gradTop.addColorStop(1,   "#0a0d16");

      const gradBot = ctx.createLinearGradient(0, botY1, 0, botY2);
      gradBot.addColorStop(0,   "#0a0d16");
      gradBot.addColorStop(0.5, "#1a1f2e");
      gradBot.addColorStop(1,   "#48506a");

      const startX = -10;

      ctx.beginPath();
      ctx.rect(startX, topY1, endX - startX, topY2 - topY1);
      ctx.fillStyle = gradTop;
      ctx.fill();
      ctx.strokeStyle = "rgba(175,198,230,0.8)";
      ctx.lineWidth = 1.2;
      ctx.stroke();

      ctx.beginPath();
      ctx.rect(startX, botY1, endX - startX, botY2 - botY1);
      ctx.fillStyle = gradBot;
      ctx.fill();
      ctx.strokeStyle = "rgba(175,198,230,0.8)";
      ctx.lineWidth = 1.2;
      ctx.stroke();

      ctx.beginPath();
      ctx.moveTo(endX, topY1);
      ctx.lineTo(endX, topY2);
      ctx.moveTo(endX, botY1);
      ctx.lineTo(endX, botY2);
      ctx.strokeStyle = "rgba(175,198,230,0.8)";
      ctx.lineWidth = 1.5;
      ctx.stroke();
    }

    ctx.restore();
  }

  clear() {
    const ctx = this.ctx;
    if (!ctx) return;
    ctx.clearRect(0, 0, ctx.canvas.width, ctx.canvas.height);
  }
}

function nacaContourExpanded(params, padLattice) {
  const raw = nacaContour(params);
  const n   = raw.length;
  const out = [];

  for (let i = 0; i < n; i++) {
    const prev = raw[(i - 1 + n) % n];
    const next = raw[(i + 1) % n];
    const tx = next[0] - prev[0];
    const ty = next[1] - prev[1];
    const len = Math.sqrt(tx * tx + ty * ty) || 1;
    const nx = -ty / len;
    const ny =  tx / len;

    out.push([
      raw[i][0] + nx * padLattice,
      raw[i][1] + ny * padLattice,
    ]);
  }

  return out;
}

function nacaContour(params) {
  const { cx, cy, chord, alpha } = params;
  const t    = 0.12;
  const cosA = Math.cos(alpha);
  const sinA = Math.sin(alpha);
  const nPts = 120;
  const upper = [], lower = [];

  for (let s = 0; s <= nPts; s++) {
    const xn = 0.5 * (1 - Math.cos(Math.PI * s / nPts));
    let yt = (t / 0.2) * chord * (
       0.2969 * Math.sqrt(xn)
      - 0.1260 * xn
      - 0.3516 * xn * xn
      + 0.2843 * xn * xn * xn
      - 0.1015 * xn * xn * xn * xn
    );
    if (xn >= 1.0) yt = 0;

    const lx  = (xn - 0.5) * chord;
    const lyU = yt;
    const lyL = -yt;

    upper.push([
      cx + lx * cosA - lyU * sinA,
      cy + lx * sinA + lyU * cosA,
    ]);
    lower.push([
      cx + lx * cosA - lyL * sinA,
      cy + lx * sinA + lyL * cosA,
    ]);
  }

  return [...upper, ...[...lower].reverse()];
}


