# Project Architecture: Phase 2 (Computational Core)

This document outlines the technical architecture for Phase 2 of the internship, focusing on the transition from 1D prototypes to 2D/3D Reduced Order Models (ROM).

## Technical Stack

### 1. High-Fidelity Solvers (Offline Phase)
- **Language:** Python 3.10+
- **IGA Library:** [Nutils](https://nutils.org/) or [GeoPDEs](http://rafavzqz.github.io/geopdes/) (final choice pending supervisor meeting).
- **Matrix Operations:** NumPy / SciPy.
- **Snapshot Generation:** Parallel execution of Full Order Model (FOM) across the sampled parameter space.

### 2. Reduced Order Modeling (Offline-to-Online)
- **Library:** [pyMOR](https://pymor.org/) for Model Order Reduction (MOR) routines.
- **Methods:** Proper Orthogonal Decomposition (POD) for basis extraction.
- **Hyper-reduction:** Empirical Cubature Method (ECM) or Energy Conserving Sampling and Weighting (ECSW) for nonlinear terms.

### 3. Digital Twin Frontend (Online Phase)
- **Language:** JavaScript (ES6+).
- **Rendering:** HTML5 Canvas (2D) / Three.js (3D).
- **Communication:** Static JSON Exports (for browser-only demos) or Flask/FastAPI (for real-time server-side querying).

## System Workflow

1.  **Sampling**: Define parameter ranges and generate samples using Latin Hypercube Sampling (LHS).
2.  **FOM Simulation**: Run high-fidelity IGA simulations for each sample.
3.  **SVD Basis**: Perform Singular Value Decomposition on the snapshot matrix to find the POD basis.
4.  **Galerkin Projection**: Project the governing equations onto the reduced basis.
5.  **Online Query**: Use the pre-computed basis and reduced operators to solve for new parameters in < 100ms.
