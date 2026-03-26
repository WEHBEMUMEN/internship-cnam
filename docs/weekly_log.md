---
layout: page
title: "Weekly Log & Planner"
---
# 📝 Internship Weekly Log & Planner

Use this document to track weekly goals, daily tasks, and notes from meetings with the supervisor (Christophe Hoareau).

---

## ❓ Open Questions for Supervisor — Priority List
*To be answered at the **first meeting**. These block the start of the Computational Core phase.*

| # | Question | Blocks |
|---|---|---|
| Q1 | What parameters should be controllable in real-time? (stiffness, damping, geometry dimensions, load?) | Parameter space definition (`docs/parameter_space.md`) |
| Q2 | Is the 2D test case a cantilever beam, a plate, or a 2D fluid tank cross-section? | FOM geometry & boundary conditions (`src/fom/`) |
| Q3 | What is the acceptable ROM accuracy tolerance vs. FOM? (e.g., < 1% relative L2 error?) | Number of POD modes & snapshots needed |
| Q4 | What is the final Digital Twin deliverable? (Web app? Jupyter notebook? Dash/Streamlit app?) | Entire frontend architecture |
| Q5 | Are there existing IGA codes in the lab (MATLAB, Python, Fortran) that I can reuse or extend? | Choice of IGA solver and implementation language |

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
- [ ] Finish reading Hoareau et al. 2025 with full annotations in `reference_catalog.md`.
- [ ] Finish reading Cottrell et al. 2009 (focus chapters on IGA for structural mechanics).
- [ ] Start reading Farhat et al. 2015 (ECSW — hyper-reduction foundations).
- [ ] Decide on IGA solver: test Nutils installation locally.
- [ ] Install pyMOR and run a simple POD example.
- [ ] Write `docs/architecture.md` with the chosen software stack.
- [ ] Write `docs/parameter_space.md` (based on Q1–Q2 answers from supervisor).
- [ ] Update `_data/tasks.yml` with all Phase 2 sub-tasks.

### 🤝 Supervisor Meeting Notes
**Date:** *(Fill in)*

**Notes:**
*(Fill in during meeting)*

**Action items for Week 3:**
- [ ] *(Fill in after meeting)*

---

*(Copy the week block above to create new weeks as the internship progresses.)*
