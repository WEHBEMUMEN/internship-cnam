# Thesis Improvement Suggestions

A collection of ideas to further strengthen the report before final submission.

---

## Writing & Structure

- **Abstract**: Make sure the abstract mentions all three pillars (IGA + Pullback + ECSW) and includes at least one concrete number (e.g., `64.5×` speedup, `<1.5 ms`).
- **Chapter transitions**: Each chapter summary (Key Takeaway box) already previews the next chapter — consider also adding a single sentence at the very start of each chapter that references the *previous* chapter, creating a smoother reading flow.
- **Avoid repetition**: The phrase *"reference domain"* appears very frequently. In some places, a pronoun or shorthand like *"$\hat{\Omega}$"* alone would be sufficient.
- **Consistent tense**: Most of the report uses present tense, but some sentences in Chapter 7 switch to past tense. Unify throughout.

---

## Content & Technical Depth

- **Convergence study**: Add a table or plot showing how the $L_2$ error changes as the number of training snapshots increases — this strengthens the offline training section.
- **NNLS tolerance sensitivity**: Briefly discuss how the choice of `ε_tol = 1×10⁻⁴` affects the final number of selected elements `E*`. Even one sentence acknowledging this would add robustness.
- **Condition number**: Mention that the penalty stiffness `kp = 10³⁰` can affect the condition number of `K_eff`, and that this is an accepted trade-off for simplicity.
- **Parameter space coverage**: In Chapter 6, note explicitly which `(x_c, r)` parameter combinations were used for training snapshots and whether a grid or Latin Hypercube sampling was used.
- **WebGL limitations**: In Chapter 7 (Future Perspectives), add a bullet point about the single-precision floating-point limitation of WebGL and how it could affect accuracy for large-scale problems.

---

## Figures & Tables

- **Figure captions**: Some captions are very brief (e.g., *"Displacement field showing deflection…"*). Consider adding one sentence explaining what the reader should take away from the figure.
- **Color consistency**: Make sure the same color scheme (blue = FOM, orange = ECSW ROM) is used consistently across all plots (Figures 6.3 and 6.5).
- **Table column widths**: Table 6.2 has uneven column widths when printed — consider adjusting `p{...}` values so all columns are visually balanced.
- **Algorithm numbering**: Double-check that Algorithm 1 (FOM Solver) and Algorithm 2 (ECSW Training) are cross-referenced in the text where they are first introduced.

---

## Bibliography & References

- **Hughes et al. (2005)**: This is cited in Chapter 3 but not formally in the bibliography — verify the full reference is listed in `main.bib`.
- **ECSW original paper**: Make sure Farhat et al. (the original ECSW paper) is cited in Section 6.4 when ECSW is first introduced formally.
- **IGA textbook**: Consider adding the Cottrell, Hughes & Bazilevs *"Isogeometric Analysis"* book (2009) as a general reference for Chapter 3.

---

## Formatting

- **Nomenclature page**: Check that all new symbols introduced in the latest edits (e.g., `J_{μ}`, `d(μ)`) are listed in the nomenclature table.
- **Overfull hboxes**: A few section titles overflow the text width (e.g., Chapter 3 Section 3.3 title). These can be fixed with `\texorpdfstring` or manual line break hints.
- **Appendix references**: In the main body, add at least one `\ref{...}` pointer to Appendix A and Appendix B so readers know where to find the detailed derivations.
