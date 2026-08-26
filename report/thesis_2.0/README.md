# Thesis 2.0: Isogeometric Analysis & Reduced Order Modeling for Structural Dynamics

> **Version**: 2.0 (Visual-First, Scannable & Human-Friendly Edition)  
> **Author**: Mumen Wehbe  
> **Supervisor**: Christophe Hoareau  

---

## 📌 Project Overview & Purpose

`thesis_2.0` is the modernized, visual-first edition of the Master 2 thesis manuscript:  
*_Isogeometric Analysis and Reduced Order Modeling for Parameterized Structural Dynamics_*.

The primary objective of Version 2.0 is to **eliminate dense wall-of-text paragraphs**, elevate **visual clarity**, and make the thesis **highly readable and scannable** for evaluators without sacrificing mathematical rigor.

---

## 📂 Directory Architecture

```
thesis_2.0/
├── README.md                      # [THIS FILE] Context & guidance for future AI agents
├── main.tex                       # Root LaTeX compiler file
├── CNAM_Logo.pdf                  # University / Lab Logo
├── references.bib                 # Bibliography dataset
└── chapters/                      # Chapter source files
    ├── 0_frontmatter/             # Title page, abstract, acknowledgements, nomenclature
    ├── 1_introduction/            # [COMPLETED] Visual-first Intro & State of the Art
    ├── 2_theory/                  # [PLACEHOLDER] Problem formulation & continuum mechanics
    ├── 3_discretization/          # [PLACEHOLDER] NURBS & Isogeometric Analysis (IGA)
    ├── 4_pullback/                # [PLACEHOLDER] Reference configuration pullback mapping
    ├── 5_simulation/              # [PLACEHOLDER] Full Order Model (FOM) & time integration
    ├── 6_rom/                     # [PLACEHOLDER] POD, ECSW hyper-reduction & Digital Twin
    ├── 7_conclusion/              # [PLACEHOLDER] Synthesis, key results & future work
    └── 8_appendix/                # [PLACEHOLDER] Detailed proofs & numerical tables
```

---

## 🤖 Instructions & Rules for Future AI Agents

When editing or updating chapters in `thesis_2.0`, **strictly follow these 6 design rules**:

### 1. **Visual-First Layout**
* **Tables over Text**: Whenever comparing 2 or more methods (e.g. FEM vs IGA, POD vs DEIM vs ECSW), use a `booktabs` table (`\toprule`, `\midrule`, `\bottomrule`) instead of long paragraphs.
* **TikZ Diagrams**: Use TikZ flowcharts for all multi-step algorithms, domain transformations, and system architectures.
* **Callout Highlights**: Wrap executive takeaways, key findings, and performance speedups in `tcolorbox` blocks.

### 2. **Human-Centric & Active Tone**
* Avoid dry, overly passive academic prose.
* Use active, direct phrasing (e.g., *"We extract a low-dimensional subspace..."* instead of *"It has been shown that a low-dimensional subspace was extracted by the algorithm..."*).
* Use bold lead-in bullet points (`\item \textbf{Keyword}: ...`) for lists.

### 3. **Inline Math Annotations**
* Annotate complex mathematical terms directly using `\underbrace{...}_{\text{explanation}}` rather than writing paragraphs explaining terms below.

### 4. **Algorithm Pseudocode**
* Use `\begin{algorithm}` + `\begin{algorithmic}` for numerical loops (time-stepping, SVD extraction, NNLS weight training).

### 5. **Hyperref & PDF String Compatibility**
* Always wrap math or Greek symbols in section titles with `\texorpdfstring{$\beta$}{beta}` to prevent hyperref PDF bookmark warnings.

### 6. **Compilation Safety**
* Always run `pdflatex -interaction=nonstopmode main.tex` to verify zero errors after making edits.

---

## 📊 Current Implementation Status

| Chapter | Status | Visual Elements Added |
| :--- | :--- | :--- |
| **0. Frontmatter** | 🟡 Placeholder | Structured Abstract & Acknowledgements |
| **1. Introduction** | 🟢 **Fully Refactored** | Roadmap Box, State-of-the-Art Comparison Table, TikZ Chapter Flowchart |
| **2. Theory** | 🟡 Placeholder | Structural template ready for text insertion |
| **3. Discretization** | 🟡 Placeholder | Structural template ready for text insertion |
| **4. Pullback** | 🟡 Placeholder | Structural template ready for text insertion |
| **5. Simulation** | 🟡 Placeholder | Structural template ready for text insertion |
| **6. ROM & Digital Twin** | 🟡 Placeholder | Structural template ready for text insertion |
| **7. Conclusion** | 🟡 Placeholder | Structural template ready for text insertion |
| **8. Appendix** | 🟡 Placeholder | Structural template ready for text insertion |
