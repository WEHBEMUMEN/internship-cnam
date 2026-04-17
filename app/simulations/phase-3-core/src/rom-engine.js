/**
 * 2D Isogeometric Reduced Order Model Engine (Phase 3.2)
 * Handles Snapshot collection, POD basis extraction, and Reduced NR solvers.
 */

class ROMEngine {
    constructor(solverFOM) {
        this.fom = solverFOM; // IGANonlinearSolver instance
        this.Phi = null;      // Reduced Basis [nDofs x k]
        this.snapshots = [];  // List of Float64Arrays
        this.singularValues = [];
    }

    /**
     * Store a displacement field from the FOM
     */
    addSnapshot(u) {
        this.snapshots.push(new Float64Array(u));
    }

    clearSnapshots() {
        this.snapshots = [];
    }

    /**
     * Compute POD Basis via SVD
     * @param {number} k Number of modes to retain
     */
    computePOD(k = 5) {
        if (this.snapshots.length === 0) throw new Error("No snapshots available");
        
        const { Matrix, SVD } = window.mlMatrix;
        
        // Build Snapshot Matrix S [nSamples x nDofs]
        const S = new Matrix(this.snapshots);
        
        // Transpose S to [nDofs x nSamples] for standard POD
        const ST = S.transpose();
        
        // Perform SVD: ST = U * Σ * V^T
        const svd = new SVD(ST, { computeLeftSingularVectors: true, computeRightSingularVectors: false });
        
        this.singularValues = svd.diagonal;
        
        // Left singular vectors are the POD modes
        const fullModes = svd.leftSingularVectors; // [nDofs x nSamples]
        
        // Truncate to k modes
        this.Phi = fullModes.subMatrix(0, fullModes.rows - 1, 0, Math.min(k, fullModes.columns) - 1);
        
        const totalEnergy = this.singularValues.reduce((s, v) => s + v*v, 0);
        const retainedEnergy = this.singularValues.slice(0, k).reduce((s, v) => s + v*v, 0);
        
        return {
            energy: retainedEnergy / totalEnergy,
            singularValues: this.singularValues.slice(0, Math.min(k * 2, this.singularValues.length))
        };
    }

    /**
     * Project Full Order Residual and Tangent into Reduced Space
     * @param {Float64Array} R Full Order Residual [nDofs]
     * @param {Array<Float64Array>} Kt Full Order Tangent [nDofs x nDofs]
     */
    projectSystem(R, Kt) {
        const { Matrix } = window.mlMatrix;
        const PhiT = this.Phi.transpose();
        const R_mat = new Matrix([Array.from(R)]).transpose();
        const Kt_mat = new Matrix(Kt);
        
        const R_red = PhiT.mmul(R_mat).to2DArray().map(r => r[0]);
        const Kt_red = PhiT.mmul(Kt_mat).mmul(this.Phi).to2DArray();
        
        return { R_red, Kt_red };
    }

    /**
     * Nonlinear ROM Solver (Reduced Newton-Raphson)
     */
    solveReduced(patch, bcs, loads, options = {}) {
        const { iterations = 10, tolerance = 1e-6, steps = 1, onProgress } = options;
        if (!this.Phi) throw new Error("POD Basis not computed");

        const { Matrix } = window.mlMatrix;
        const k = this.Phi.columns;
        const nDofs = this.Phi.rows;
        
        let ur = new Float64Array(k).fill(0); // Reduced coordinates
        const residualHistory = [];

        // Total External Load (Full)
        const F_ext_total = new Float64Array(nDofs).fill(0);
        const nV = patch.controlPoints[0].length;
        loads.forEach(load => {
            const idx = (load.i * nV + load.j) * 2;
            F_ext_total[idx] += load.fx;
            F_ext_total[idx + 1] += load.fy;
        });

        // Fixed DOFs handling in ROM is usually done by projecting the BC-constrained FOM
        // or by using Lagrange multipliers. For simplicity here, we assume the basis Phi
        // already satisfies homogeneous BCs (which it does if snapshots did).
        
        for (let s = 1; s <= steps; s++) {
            const F_ext = F_ext_total.map(f => f * (s / steps));

            for (let iter = 0; iter < iterations; iter++) {
                // 1. Expand to Full Space: u = Phi * ur
                const ur_mat = new Matrix([Array.from(ur)]).transpose();
                const u_full_mat = this.Phi.mmul(ur_mat);
                const u_full = new Float64Array(u_full_mat.to2DArray().map(r => r[0]));

                // 2. Evaluate FOM Residual and Tangent
                const F_int = this.fom.calculateInternalForce(patch, u_full);
                const R_full = F_ext.map((f, i) => f - F_int[i]);
                
                const Kt_full = this.fom.calculateTangentStiffness(patch, u_full);
                this.fom.applyPenaltyConstraints(Kt_full, patch);

                // 3. Project to Reduced Space
                const { R_red, Kt_red } = this.projectSystem(R_full, Kt_full);

                // 4. Convergence Check (Reduced Norm)
                let resNorm = 0;
                R_red.forEach(ri => resNorm += ri * ri);
                const norm = Math.sqrt(resNorm);
                residualHistory.push({step: s, iter, norm});

                if (onProgress) onProgress({ step: s, iter, norm, isReduced: true });
                if (norm < tolerance && iter > 0) break;

                // 5. Solve Reduced System: Kt_red * dur = R_red
                const dur = this.fom.gaussianElimination(Kt_red, R_red);
                for (let i = 0; i < k; i++) ur[i] += dur[i];
            }
        }

        // Return reconstructed full displacement
        const final_ur_mat = new Matrix([Array.from(ur)]).transpose();
        const final_u_full = this.Phi.mmul(final_ur_mat).to2DArray().map(r => r[0]);
        
        return { 
            u: new Float64Array(final_u_full), 
            ur: ur,
            residualHistory 
        };
    }
}

// Expose globally for script-tag loading compatibility
window.ROMEngine = ROMEngine;
