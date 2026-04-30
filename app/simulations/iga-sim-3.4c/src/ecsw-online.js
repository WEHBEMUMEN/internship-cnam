/**
 * ECSW Online Engine
 * Runtime assembly and solving.
 */
class ECSWOnlineEngine {
    constructor() {
        this.sampleElements = [];
        this.Kp_red = null;
        this.Phi = null;
        this.k = 0;
        this.nDofs = 0;
    }

    /**
     * Set the trained state.
     */
    init(trainedData, Phi) {
        this.sampleElements = trainedData.sampleElements;
        this.Kp_red = trainedData.Kp_red;
        this.Phi = Phi;
        this.k = Phi.columns;
        this.nDofs = Phi.rows;
    }

    /**
     * Weighted Assembly Online.
     */
    assembleReducedSystem(fomSolver, patch, ur, u_full) {
        const F_red_assembly = new Float64Array(this.k).fill(0);
        const Kt_red_assembly = Array.from({ length: this.k }, () => new Float64Array(this.k).fill(0));
        
        const { p, q } = patch;
        const gRule = window.GaussQuadrature2D.getPoints(Math.max(p, q) + 1);

        this.sampleElements.forEach(el => {
            const { f_e, k_e, activeDofs } = window.ECSWCore._assembleElementPhysics(fomSolver, patch, el, u_full, gRule);
            
            const nLocal = activeDofs.length;
            const weight = el.weight;

            // 1. Reduced Force: F_red += w * Phi_e^T * f_e
            for (let i = 0; i < this.k; i++) {
                let dot = 0;
                for (let a = 0; a < nLocal; a++) {
                    dot += this.Phi.get(activeDofs[a], i) * f_e[a];
                }
                F_red_assembly[i] += weight * dot;
            }

            // 2. Reduced Tangent: Kt_red += w * Phi_e^T * k_e * Phi_e
            const Phi_e = Array.from({ length: nLocal }, () => new Float64Array(this.k));
            for (let a = 0; a < nLocal; a++) {
                for (let i = 0; i < this.k; i++) {
                    Phi_e[a][i] = this.Phi.get(activeDofs[a], i);
                }
            }

            // Matrix multiplication: (Phi_e^T * k_e) * Phi_e
            const temp = Array.from({ length: nLocal }, () => new Float64Array(this.k));
            for (let a = 0; a < nLocal; a++) {
                for (let i = 0; i < this.k; i++) {
                    let val = 0;
                    for (let b = 0; b < nLocal; b++) {
                        val += k_e[a][b] * Phi_e[b][i];
                    }
                    temp[a][i] = val;
                }
            }

            for (let i = 0; i < this.k; i++) {
                for (let j = 0; j < this.k; j++) {
                    let val = 0;
                    for (let a = 0; a < nLocal; a++) {
                        val += Phi_e[a][i] * temp[a][j];
                    }
                    Kt_red_assembly[i][j] += weight * val;
                }
            }
        });

        // 3. Add Precomputed Reduced Penalty (linear)
        for (let i = 0; i < this.k; i++) {
            for (let j = 0; j < this.k; j++) {
                F_red_assembly[i] += this.Kp_red[i][j] * ur[j];
                Kt_red_assembly[i][j] += this.Kp_red[i][j];
            }
        }

        return { F_red: F_red_assembly, Kt_red: Kt_red_assembly };
    }

    /**
     * Main Online Solver (Newton-Raphson)
     */
    solveReduced(fomSolver, patch, loads, options = {}) {
        const { steps = 1, iterations = 15, tolerance = 1e-6 } = options;
        const k = this.k;
        const nDofs = this.nDofs;
        const nV = patch.controlPoints[0].length;

        let ur = new Float64Array(k);
        const residualHistory = [];

        // Precompute projected external force (Linear)
        const F_ext_total = new Float64Array(nDofs);
        loads.forEach(l => {
            const idx = (l.i * nV + l.j) * 2;
            F_ext_total[idx] += l.fx;
            F_ext_total[idx+1] += l.fy;
        });

        const F_ext_red = new Float64Array(k);
        for (let i = 0; i < k; i++) {
            for (let d = 0; d < nDofs; d++) {
                F_ext_red[i] += this.Phi.get(d, i) * F_ext_total[d];
            }
        }

        for (let s = 1; s <= steps; s++) {
            const frac = s / steps;
            
            for (let iter = 0; iter < iterations; iter++) {
                // Pre-expand ur to full u for assembly and diagnostics
                const u_full = new Float64Array(this.nDofs);
                for (let j = 0; j < this.k; j++) {
                    const val = ur[j];
                    if (val === 0) continue;
                    for (let d = 0; d < this.nDofs; d++) u_full[d] += this.Phi.get(d, j) * val;
                }

                // 1. Online Assembly
                const { F_red, Kt_red } = this.assembleReducedSystem(fomSolver, patch, ur, u_full);
                
                // --- DIAGNOSTIC: Compare with full projection for audit ---
                if (window.ECSWAudit) {
                    const F_int_full = fomSolver.calculateInternalForce(patch, u_full);
                    fomSolver.applyPenaltyConstraints(null, F_int_full, u_full, patch);
                    const F_red_exact = new Float64Array(this.k);
                    for (let i = 0; i < this.k; i++) {
                        for (let d = 0; d < this.nDofs; d++) F_red_exact[i] += this.Phi.get(d, i) * F_int_full[d];
                    }
                    
                    let err2 = 0, norm2 = 0;
                    for (let i = 0; i < this.k; i++) {
                        err2 += (F_red[i] - F_red_exact[i])**2;
                        norm2 += F_red_exact[i]**2;
                    }
                    const relError = Math.sqrt(err2) / Math.max(Math.sqrt(norm2), 1e-12);
                    
                    const residualNorm = Math.sqrt(F_red.reduce((a, b) => a + b*b, 0));
                    window.ECSWAudit.logIteration(iter, residualNorm, relError);
                }
                // ---------------------------------------------------------

                // 2. Reduced Residual
                const R_red = F_ext_red.map((f, i) => f * frac - F_red[i]);
                const norm = Math.sqrt(R_red.reduce((a, b) => a + b*b, 0));

                if (norm < tolerance && iter > 0) break;
                if (isNaN(norm)) {
                    console.error("ECSW: Newton exploded (NaN residual)");
                    break;
                }

                // 3. Solve Kt_red * dur = R_red
                const dur = fomSolver.gaussianElimination(Kt_red, Array.from(R_red));
                
                // 4. Update
                for (let i = 0; i < k; i++) ur[i] += dur[i];

                residualHistory.push({ step: s, iter, norm });
            }
        }

        // Reconstruct full displacement field
        const u = new Float64Array(nDofs);
        for (let j = 0; j < k; j++) {
            const val = ur[j];
            for (let d = 0; d < nDofs; d++) {
                u[d] += this.Phi.get(d, j) * val;
            }
        }

        return { u, ur, residualHistory };
    }
}

window.ECSWOnlineEngine = ECSWOnlineEngine;
