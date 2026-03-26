# Parameter Space Definition

This document defines the parameters ($\mu$) that will be varied during the Phase 2 Computational Core simulations to generate the snapshot database for the Reduced Order Model.

## Parameter Vector ($\mu$)

The behavior of the structural system is governed by a vector of parameters $\mu \in \mathcal{P}$, where $\mathcal{P}$ is the parameter space.

| Parameter | Symbol | Description | Initial Range |
| :--- | :--- | :--- | :--- |
| **Young's Modulus** | $E$ | Material stiffness (Linear/Nonlinear) | $[1 \times 10^6, 2 \times 10^8]$ Pa |
| **Poisson's Ratio** | $\nu$ | Transverse strain ratio | $[0.25, 0.49]$ |
| **Load Intensity** | $F$ | Magnitude of the applied point/surface load | $[0, 1000]$ N |
| **Geometry (L/H)** | $\alpha$ | Aspect ratio of the beam/plate | $[5, 20]$ |
| **Damping Ratio** | $\zeta$ | Rayleigh damping coefficients | $[0.01, 0.1]$ |

## Sampling Strategy

- **Method**: Latin Hypercube Sampling (LHS).
- **Constraint**: Uniform distribution across all ranges initially.
- **Sample Size**: 25 snapshots for initial POD validation; 100+ for final hyper-reduction (ECSW).

## 2D Test Case Geometry

- **Type**: Cantilever Plate (Physical model counterpart ready for 3D printing).
- **Parametrization**: Single-patch NURBS surface.
- **Boundary Conditions**: Fixed (Clamped) on the left edge ($x=0$), Free elsewhere.
