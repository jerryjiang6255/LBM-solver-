# Lattice Boltzmann Flow Visualizer · V5

A real-time 2D fluid dynamics simulator running entirely in the browser, built on the Lattice Boltzmann Method (LBM). No server, no installs — just open the link and the physics runs instantly.

**[Live Demo →](https://jerryjiang6255.github.io/LBM-solver-/)**

---

## What is LBM?

Traditional CFD solvers (Fluent, OpenFOAM) discretise the Navier-Stokes equations directly — solving for pressure and velocity fields on a mesh. This requires complex mesh generation, and remeshing every time geometry changes.

LBM takes a fundamentally different approach: instead of solving macroscopic fluid equations, it simulates the statistical behaviour of particle distribution functions on a discrete lattice. At each timestep, particles **collide** (redistribute momentum locally) and **stream** (propagate to neighbouring nodes). The macroscopic fluid behaviour — velocity, pressure, vorticity — emerges from these two simple steps.

This makes LBM:
- **Naturally parallel** — every node updates independently
- **Mesh-free for geometry changes** — new obstacles just flip nodes from fluid to solid, no remeshing needed
- **Efficient in the browser** — the entire solver runs in typed JavaScript arrays with no GPU required

---

## Solver Details

| Component | Implementation |
|---|---|
| Lattice | D2Q9 (2D, 9 velocity directions) |
| Grid | 256 × 128 nodes (32,768 total) |
| Collision | TRT (Two-Relaxation-Time) |
| Turbulence | Smagorinsky LES eddy viscosity |
| Boundary (geometry) | BFL interpolated bounce-back with SDF q values |
| Boundary (inlet) | Hard Dirichlet reset |
| Boundary (outlet) | Zero-gradient copy |
| Boundary (top/bottom) | Far-field Dirichlet |
| Geometry representation | Signed Distance Functions (SDF) |

### TRT Collision
Unlike BGK (single relaxation time), TRT uses two rates: **ω+** controls physical viscosity and **ω−** controls numerical stability. The Magic Parameter Λ = 3/16 ties them together, eliminating the viscosity-dependent boundary error that BGK suffers from near solid walls.

### Smagorinsky LES
Adds an eddy viscosity term νt = (Cs·Δ)²·|S| to the base fluid viscosity. This damps sub-grid turbulent fluctuations without resolving them explicitly, allowing stable simulation at higher Reynolds numbers than a laminar solver would permit.

### BFL Boundary
The Bouzidi-Firdaouss-Lallemand boundary condition uses the fractional wall distance **q** (computed from the SDF at startup) to reconstruct distribution functions with sub-lattice accuracy. This places the effective no-slip wall at the true geometric boundary, not at a half-node offset.

---

## Features

**Geometry options**
- NACA 0012 airfoil with adjustable angle of attack (−20° to +20°)
- Cylinder
- Staggered block matrix (porous media / turbulence generation)
- Empty domain

**Display modes**
- Speed magnitude
- Vorticity (curl of velocity field)

**Sandbox mode**
Toggle sandbox to draw arbitrary walls directly on the canvas with click-and-drag. The flow field responds in real time — no reset, no remesh. Walls can also be erased. This is the clearest demonstration of LBM's geometric flexibility over traditional CFD.

**Live controls**
- Inlet velocity (lattice units)
- Reynolds number (live update, no restart needed)
- Colormap scale
- Sim speed (steps per rendered frame)
- Pause / Resume / Reset

---

## Running Locally

No build step needed. Serve the repo root with any static file server:

```bash
# Python
python -m http.server 8080

# Node (npx)
npx serve .

# VS Code
# Use the Live Server extension
```

Then open `http://localhost:8080` in Chrome or Firefox.

> The solver uses ES modules (`import`/`export`) so it must be served over HTTP — opening `index.html` directly as a `file://` URL will not work.

---

## Project Structure

```
LBM-solver-/
├── index.html          # UI, wiring, frame loop
├── style.css           # All styling
├── preview.png         # Open Graph preview image
└── src/
    ├── lbm.js          # Lattice constants, shared arrays, field init
    ├── geometry.js     # SDF geometry builders (cylinder, airfoil, blocks)
    ├── collisionFused.js  # TRT collision + Smagorinsky LES + wall model
    ├── streaming.js    # BFL interpolated streaming
    ├── boundary.js     # Inlet/outlet/wall boundary conditions
    ├── render.js       # Colormap, offscreen canvas upscale, body vector draw
    └── userDraw.js     # Sandbox drawing layer, Catmull-Rom spline render
```

---

## Physics Limitations

This is an educational simulator, not a production CFD tool. Known limitations:

- **Low Mach number only** — LBM is weakly compressible; the solver becomes unstable above u₀ ≈ 0.3–0.4 lattice units
- **2D only** — real aerodynamic flows are 3D; lift/drag values here are qualitative, not quantitative
- **No wall-resolved turbulence** — the Smagorinsky model stabilises the solver but does not replace proper LES resolution near walls
- **Coarse grid** — 256 × 128 is chosen for real-time browser performance; a production simulation would use 10–100× more nodes

---

## Built By

**Jerry Jiang** · [LinkedIn](https://www.linkedin.com/in/jerry-jiang-a38127285/)

Any feedback or advice is greatly appreciated.
