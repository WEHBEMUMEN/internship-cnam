# Digital Twin & Reduced Order Modelling for Nonlinear Dynamics

[![GitHub Pages](https://img.shields.io/badge/Live_Demo-GitHub_Pages-blue?logo=github)](https://wehbemumen.github.io/internship-cnam/)
[![Institution](https://img.shields.io/badge/Cnam-EPN04_Paris-red)](https://www.cnam.fr/)
[![Methodology](https://img.shields.io/badge/Methodology-IGA_%2B_ROM-green)]()

**Author:** Mumen Wehbe  
**Supervisor:** Dr. Christophe Hoareau  
**Institution:** Conservatoire National des Arts et Métiers (Cnam), Équipement pour la Performance et le Numérique (EPN04), Paris  
**Duration:** March – September 2026  

---

## 📌 Executive Summary

This project presents the design and implementation of a real-time, web-based **Digital Twin** platform for geometrically nonlinear structural mechanics. By bridging **Isogeometric Analysis (IGA)** with projection-based **Reduced Order Modelling (ROM)** techniques, the platform executes high-fidelity physics natively inside modern web browsers at sub-10 ms solve times.

The platform serves as both a research validation engine and an interactive pedagogical tool for advanced computational mechanics courses (e.g. Non-linear Dynamics, Fluid-Structure Interaction).

---

## 🔬 Key Scientific Pillars

1. **Isogeometric Analysis (IGA):**
   - Exact CAD geometry discretization using **Non-Uniform Rational B-Splines (NURBS)**.
   - Elimination of mesh generation errors.
   - Built-in support for $h$-refinement (knot insertion), $p$-refinement (degree elevation), and $k$-refinement (maximum inter-element smoothness $C^{p-1}$).

2. **Geometrically Nonlinear Continuum Mechanics:**
   - Total Lagrangian formulation on a fixed bi-unit master domain $\hat{\Omega} = [0,1]^2$.
   - **Green–Lagrange strain tensor** $\mathbf{E}$ and **Second Piola–Kirchhoff stress** $\mathbf{S}$.
   - Newton–Raphson nonlinear solver with tangent stiffness matrix updates.

3. **Projection-Based Reduced Order Modelling (ROM):**
   - Snapshot generation across parametric configurations.
   - **Proper Orthogonal Decomposition (POD)** via Singular Value Decomposition (SVD).
   - Galerkin projection reducing systems with hundreds of DOFs to a compact reduced basis ($\sim 5\text{--}10$ modes).
   - Hyper-reduction concepts (DEIM / ECSW) for nonlinear internal force evaluation.

---

## 🎮 Interactive Simulation Suite

- **[Essential IGA Playground](app/playground.html):** A clean, unified sandbox focusing on core milestones:
  - **1D IGA Rod Baseline & Refinements:** Phases 1.1, 1.2 ($h$), 1.3 ($p$), 1.5 ($k$)
  - **2D Elastic Plate Benchmark:** Phase 2.1
  - **Infinite Plate with Circular Hole:** Full Order Model (Phase 5.3) & Reduced Order Model (Phase 5.3a)
- **[Full Development Versions Catalog](app/index.html):** All iterative research prototypes and flowcharts across Phases 1 through 5.

---

## 🏗️ Repository Architecture

```
├── app/
│   ├── playground.html       # Streamlined 1-page interactive playground
│   ├── index.html            # Complete catalog of all development versions
│   ├── simulations/          # Modular standalone IGA & ROM simulation engines
│   └── shared/               # Reusable math, NURBS, and solver utilities
├── report/                   # Complete LaTeX master thesis report & chapters
│   ├── main.tex              # Master LaTeX document
│   └── chapters/             # Chapters 1-4, intro, conclusion, titlepage
├── theory/                   # Interactive theory documentation & glossary
│   ├── overview.html         # Research abstract & project executive summary
│   ├── timeline.html         # Internship Gantt chart & roadmap
│   ├── tasks.html            # Work package task tracker
│   └── dictionary.html       # Mathematical glossary & definition cards
├── index.html                # Main project landing page
├── index.css                 # Unified styling and design tokens
└── .nojekyll                 # GitHub Pages direct deployment bypass
```

---

## 🚀 Running Locally

You can run the web platform with any static local web server:

```bash
# Python 3
python -m http.server 8000

# Or Node.js
npx serve .
```

Open `http://localhost:8000` in your web browser.

---

## 📄 License & Attribution
© 2026 Mumen Wehbe & EPN04, Cnam (Paris). All rights reserved.