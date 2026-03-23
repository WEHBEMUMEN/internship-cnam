---
layout: page
title: "Reference Catalog"
---
# 📚 Literature & Resource Catalog

Annotated knowledge base for the 2026 CNAM internship on Digital Twin & Reduced Order Modelling.

> **Status key:** `🔴 To Read` | `🟡 In Progress` | `🟢 Done` | `⚡ Skim Only`

---

## 🔗 Quick Links
* **Shared resources Mumen:** [Link](https://lightcoral-goat-369310.hostingersite.com/)
* **Rixen textbook (Nonlinear ROM geometric parameters):** [Google Books](https://books.google.fr/books?hl=fr&lr=&id=LjXYEAAAQBAJ&oi=fnd&pg=PA189&dq=Rixen+nonlinear+reduced+order+model+geometric+parameter+circle&ots=Kxp-5jB0ex&sig=Kok9nJTqY1afuQUeFXLCsP5Qqug&redir_esc=y#v=onepage&q&f=false)

---

## ⭐ Category 1: IGA Foundations
*Isogeometric Analysis — the discretization method used for the Full Order Model.*

### Cottrell, J. A., Hughes, T. J. R., & Bazilevs, Y. (2009)
*Isogeometric Analysis: Toward Integration of CAD and FEA* | **🟡 In Progress**

- **Why it matters:** The foundational textbook — defines the NURBS basis functions used in the IGA FOM.
- **Key Topics:** B-splines, NURBS, k-refinement, T-splines, isoparametric concept, patch coupling.
- **Key Takeaways:** *(Fill in while reading — chapter by chapter)*
  - Ch. 1–2: *(mesh-free concept, motivation over classical FEM)*
  - Ch. 3–4: *(NURBS basis functions, geometric exactness)*
  - Ch. 5+: *(structural mechanics applications)*
- **Connection to project:** Direct foundation for FOM assembly in `src/fom/beam_iga.py`.
- **Cross-links:** Hoareau 2025 uses IGA as the FOM backbone.

---

## ⚡ Category 2: ROM & Hyper-reduction
*Reduced Order Modelling methods — POD, DEIM, ECS, ECSW.*

### Benner, P., Gugercin, S., & Willcox, K. (2015)
*A Survey of Projection-Based Model Reduction Methods for Parametric Dynamical Systems* | **🔴 To Read**

- **Why it matters:** Comprehensive overview of the entire ROM landscape — gives vocabulary and context.
- **Key Topics:** POD, balanced truncation, reduced basis methods, parametric dependencies.
- **Key Takeaways:** *(Fill in after reading)*
- **Connection to project:** Background reading for Phases 2C–2E. Read for overview, not for implementation details.
- **Cross-links:** All methods surveyed are candidates for Phase 2D.

---

### Farhat, C., Avery, P., Chapman, T., & Cortial, J. (2015)
*Structure-Preserving, Stability, and Accuracy Properties of the Energy-Conserving Sampling and Weighting (ECSW) Method for Hyper-Reduction of Nonlinear Finite Element Dynamic Models* | **🔴 To Read**

- **Why it matters:** Defines the **ECSW** method — one of the three hyper-reduction approaches to be compared.
- **Key Topics:** Energy-conserving sampling, non-negative weights, stability in dynamics.
- **Key Takeaways:** *(Fill in after reading)*
- **Connection to project:** Phase 2E — implementation and benchmarking of ECSW.
- **Cross-links:** Hoareau 2025 references this approach. Compare against DEIM (Chaturantabut 2010).

---

### Fonn, E., et al. (2025)
*Least-Squares Projected Models for Non-Intrusive Affinization* | **🔴 To Read**

- **Why it matters:** Covers non-intrusive affine decomposition — relevant if we cannot access FOM assembly code directly.
- **Key Topics:** Reduced basis, least-squares projection, non-intrusive ROM.
- **Key Takeaways:** *(Fill in after reading)*
- **Connection to project:** Phase 2F (online phase) — may simplify the parameter interpolation step.
- **Cross-links:** Complementary to pyMOR's non-intrusive pipeline.

---

## 🔥 Category 3: FSI & Sloshing (Supervisor's Work)
*The direct research lineage of this internship. Read these first.*

### Hoareau, C., et al. (2025)
*Projection-Based ROM and Hyper-reduction of Linear Sloshing...* | **🔴 To Read** ⚠️ **READ FIRST**

- **Why it matters:** **This is the supervisor's paper.** It defines the exact methodology this internship extends.
- **Key Topics:** IGA FOM for sloshing, POD-based ROM, hyper-reduction applied to FSI, error analysis.
- **Key Takeaways:** *(Fill in page by page)*
  - Method used: *(DEIM? ECSW? Note here after reading)*
  - Test case geometry: *(Cylindrical tank? Rectangular? Note here)*
  - Parameter space: *(What parameters vary? Note here)*
  - Achieved speedup: *(Note the reported computational gain)*
- **Connection to project:** This is the starting point — the internship extends this to nonlinear dynamics and geometrical parameters.
- **Cross-links:** Uses IGA (Cottrell 2009) + ROM (Benner 2015 survey) + hyper-reduction (Farhat 2015).

---

### LGST-LAB / slosh_ml (GitHub Repository)
*Machine Learning application for sloshing dynamics* | **🔴 To Explore**

- **Link:** [https://github.com/LGST-LAB/slosh_ml](https://github.com/LGST-LAB/slosh_ml)
- **Why it matters:** May contain reference implementations or datasets relevant to the sloshing problem.
- **Key Topics:** ML surrogate for sloshing, data-driven ROM.
- **Key Takeaways:** *(Fill in after exploring codebase)*
  - Code language: *(Python/MATLAB?)*
  - Existing datasets: *(Snapshot data available? Note here)*
- **Connection to project:** Potential data source or methodology inspiration.

---

## 📐 Category 4: Nonlinear Structural Dynamics
*FEM theory for large-displacement nonlinear problems — the physics foundation.*

### Belytschko, T., Liu, W. K., & Moran, B. (2014)
*Nonlinear Finite Elements for Continua and Structures* | **🔴 To Read** *(selective chapters)*

- **Why it matters:** Reference for the nonlinear FEM formulation used in the FOM.
- **Key Topics:** Updated/Total Lagrangian formulation, geometric nonlinearity, Newton-Raphson iteration, Newmark-β time integration.
- **Key Takeaways:** *(Fill in after reading selected chapters)*
  - Ch. on geometric nonlinearity: *(Note the stiffness tangent formulation)*
  - Newmark-β: *(Note numerical stability conditions)*
- **Connection to project:** Phase 2B — FOM assembly with large-displacement kinematics.
- **Reading strategy:** Do NOT read cover-to-cover. Only the chapters on geometric nonlinearity and time integration are needed.

---

### Kim, Y., et al. (2022)
*Improved Nonlinear Analysis of a Propeller Blade Using ROM* | **⚡ Skim Only**

- **Why it matters:** Example of ROM applied to a nonlinear structural problem (propeller blade).
- **Key Topics:** Geometric nonlinearity, ROM application in aeroelastics.
- **Key Takeaways:** *(Skim for the ROM accuracy and speedup results — compare to our targets)*
- **Connection to project:** Validation benchmark — our results should be comparable in accuracy.

---

## 🌐 Category 5: Digital Twin Architecture
*Methods and frameworks for real-time virtual/physical coupling.*

*(No papers added yet — search for references on: "Digital Twin real-time simulation", "Digital Twin FSI", "Industry 4.0 Digital Twin structural mechanics")*

**Possible references to add:**
- Grieves, M. & Vickers, J. (2017) — *Digital Twin: Mitigating Unpredictable, Undesirable Emergent Behavior in Complex Systems*
- Tao, F., et al. (2019) — *Digital Twins and Cyber-Physical Systems toward Smart Manufacturing*

---

## 🛠️ How to Add a New Reference

Copy this template and paste it into the correct category:

```markdown
### Author(s) (Year)
*Title* | **🔴 To Read**

- **Why it matters:**
- **Key Topics:**
- **Key Takeaways:** *(Fill in after reading)*
- **Connection to project:**
- **Cross-links:**
```
