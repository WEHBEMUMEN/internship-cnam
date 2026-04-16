---
layout: page
title: "Weekly Log & Planner"
---
# 📝 Internship Weekly Log & Planner

Use this document to track weekly goals, daily tasks, and notes from meetings with the supervisor (Christophe Hoareau).

---

## ❓ Open Questions for Supervisor — Priority List
*To be answered at the **first meeting**. These block the start of the Computational Core phase.*

| # | Question | My Proposed Default | Blocks |
|---|---|---|---|
| Q1 | What parameters should be controllable in real-time? | E, ν, F, L/H, ζ — as defined in `docs/parameter_space.md`. Should I add/remove any? | Parameter space definition |
| Q2 | Is the 2D test case a cantilever beam, a plate, or a 2D fluid tank cross-section? | Cantilever plate (single-patch NURBS). I have a working 2D surface viewer in `phase-2-core/`. | FOM geometry & boundary conditions |
| Q3 | What is the acceptable ROM accuracy tolerance vs. FOM? | < 1% relative L2 error, 99.9% energy capture via SVD (as in `docs/architecture.md`). | Number of POD modes & snapshots needed |
| Q4 | What is the final Digital Twin deliverable? | Web app (the entire hub is already web-based with Three.js for 3D). | Entire frontend architecture |
| Q5 | Are there existing IGA codes in the lab I can reuse or extend? | I've explored `slosh_ml`. If none, I'll proceed with Nutils (Python IGA). | Choice of IGA solver & language |

**Answers (fill in during meeting):**
- **A1:**
- **A2:**
- **A3:**
- **A4:**
- **A5:**

---

## 📅 Week 1: March 23 – March 29, 2026
**Primary Goal:** Repository setup, initial literature review, and understanding the digital twin architecture.

### 📋 Task Tracker
- [x] Set up GitHub repository and file structure.
- [x] Create project Gantt chart (timeline.html) and deploy GitHub Pages site.
- [x] Configure Sveltia CMS with tasks and downloads collections.
- [x] Restructure literature reference catalog into 5 annotated categories.
- [x] Read Hoareau et al. 2025 (supervisor's paper) — **first priority**.
- [x] Read Cottrell et al. 2009 chapters 1–4 (IGA foundations).
- [x] Check "Shared resources Mumen" link and catalog contents.
- [x] Explore `slosh_ml` GitHub repository.
- [x] **Milestone:** Completed Phase 1.6 IGA Simulator & Phase 1.7 ROM Benchmark.

### 🤝 Supervisor Meeting Notes
**Date:** *(Fill in when meeting is scheduled)*

**Pre-meeting prep:**
- Read at least the abstract and Sec. 1–2 of Hoareau 2025 before meeting.
- Prepare answers/thoughts on the Q1–Q5 above.

**Notes/Answers from meeting:**
*(Write notes here during the meeting)*

**Action items for Week 2:**
- [ ] *(Fill in after meeting)*

---

## 📅 Week 2: March 30 – April 5, 2026
**Primary Goal:** Complete literature review of primary papers + lock in software stack decision.

### 📋 Task Tracker
  - [x] Finish reading Hoareau et al. 2025 (nonlinear foundations).
  - [x] Write `docs/architecture.md` with the chosen software stack.
  - [x] Write `docs/parameter_space.md`.
  - [x] Update `_data/tasks.yml` with all Phase 2 sub-tasks.
  - [x] **Milestone:** Updated Project Branding (WEHBE Mumen & EPN04).
  - [x] **Milestone:** Phase 1.7 Engine-Theory Alignment (Gauss Quadrature, Hyper-reduction, Newmark-beta Dynamics).
  - [x] **Milestone:** Completed Phase 1.8 (Circle) & Phase 1.NL (Nonlinear IGA/ROM 1D).
  - [x] **Milestone:** Phase 2.0 2D NURBS Surface Mapping & Verification Lab.

### 🤝 Supervisor Meeting Notes
**Date:** *(Fill in)*

**Notes:**
*(Fill in during meeting)*

**Action items for Week 3:**
- [ ] *(Fill in after meeting)*

---

## 📅 Week 3: April 6 – April 12, 2026
**Primary Goal:** Transition to 2D Physics (Stiffness Assembly) and Snapshot Database Generation.

### 📋 Task Tracker
- [x] **Shift to Pure JavaScript Core**: Removed Python/Nutils backend dependency; all IGA assembly and solving now runs natively in the browser.
- [x] **Kirsch Benchmark Implementation**: Assembled 2D Stiffness Matrix (K) for the 'Plate with a Hole' validation case.
- [x] **Penalty Method Resolution**: Solved the coincident corner points issue using $10^{12}$ penalty constraints to prevent numerical locking.
- [x] **Stress Field Recovery**: Implemented post-processing for $\sigma_{xx}$, $\sigma_{yy}$, and $\sigma_{xy}$ with analytical comparison.
- [x] **Snapshot Generation (2C)**: Developed native JS Latin Hypercube Sampler (LHS) and generated initial snapshot database.
- [x] **Internal Fixes**: Resolved path issues in the Kanban board and Repository links for correct Netlify deployment.

### 🤝 Supervisor Meeting Notes
**Date:** *(Fill in)*

**Notes:**
*(Fill in during meeting)*

**Action items for Week 4:**
- [ ] *(Fill in after meeting)*

---

## 📅 Week 4: April 13 – April 19, 2026
**Primary Goal:** Finalize 2D Structural Nonlinearity (Geometric) and expand Theory Docs for 2D ROM.

### 📋 Task Tracker
- [x] **Geometric Nonlinearity (3A)**: Implemented a 2D Isogeometric Nonlinear engine using Green-Lagrange strain and Newton-Raphson iteration.
- [x] **Penalty Stabilization**: Deployed robust penalty constraints to handle geometric singularities in corner-degenerate NURBS patches.
- [x] **Interactive Analytics**: Integrated `chartjs-plugin-zoom` for professional-grade Stress-Strain and Force-Displacement telemetry.
- [x] **Theory Suite Expansion**: Authored high-fidelity documentation for **2D POD (theory-2d1.html)** and **2D Projection (theory-2d2.html)**.
- [x] **Hub Synchronization**: Fully updated the simulation hub with Phase 3A integration and logic-path flowcharts.
- [x] **GitHub Deployment**: Synchronized the local codebase with the remote repository, ensuring full architectural integrity.

### 🤝 Supervisor Meeting Notes
**Date:** *(Scheduled for end of week)*

**Pre-meeting prep:**
- Demonstrate the 60 FPS nonlinear structural divergence on the Plate with Hole benchmark.
- Show the Energy Capture (Scree Plot) in the new SVD Lab documentation.

**Notes/Answers from meeting:**
*(Fill in during meeting)*

**Action items for Week 5:**
- [ ] Begin Phase 3B: Implementation of Hyperelastic (Neo-Hookean) material models.

---

*(Copy the week block above to create new weeks as the internship progresses.)*
