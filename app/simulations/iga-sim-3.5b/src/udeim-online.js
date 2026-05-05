/**
 * Phase 3.4b: U-DEIM Online Solver
 */

UDEIMEngine.prototype.solveReduced = function(fomSolver, romEngine, patch, bcs, loads, options = {}) {
    const { iterations = 15, tolerance = 1e-6, steps = 1, quiet = false } = options;
    const Phi = romEngine.Phi;
    const k = Phi.columns;
    const nDofs = Phi.rows;
    const nV = patch.controlPoints[0].length;
    const { Matrix } = window.mlMatrix;

    let ur = new Float64Array(k);
    const residualHistory = [];

    // Build full external force vector
    const F_ext_total = new Float64Array(nDofs);
    loads.forEach(l => {
        const idx = (l.i * nV + l.j) * 2;
        F_ext_total[idx] += l.fx;
        F_ext_total[idx + 1] += l.fy;
    });

    const PhiT = Phi.transpose();
    const F_ext_red_total = new Float64Array(k);
    for (let i = 0; i < k; i++) {
        let dot = 0;
        for (let d = 0; d < nDofs; d++) dot += PhiT.get(i, d) * F_ext_total[d];
        F_ext_red_total[i] = dot;
    }

    const gRule = window.GaussQuadrature2D.getPoints(Math.max(patch.p, patch.q) + 1);


    for (let s = 1; s <= steps; s++) {
        const loadFraction = s / steps;
        
        // Tangent matrix evaluation moved inside the iter loop for Full Newton-Raphson stability

        for (let iter = 0; iter < iterations; iter++) {
            const u_full = new Float64Array(nDofs);
            for (let d = 0; d < nDofs; d++) {
                for (let j = 0; j < k; j++) u_full[d] += Phi.get(d, j) * ur[j];
            }

            // Recompute full tangent for accuracy (Galerkin-style) 
            // since hyper-reduction for tangent is not yet implemented for U-DEIM
            const Kt_full = fomSolver.calculateTangentStiffness(patch, u_full);
            fomSolver.applyPenaltyConstraints(Kt_full, null, u_full, patch, bcs);
            const Kt_red = PhiT.mmul(new Matrix(Kt_full.map(r => Array.from(r)))).mmul(Phi).to2DArray();

            // 1. Evaluate only active elements
            const f_sampled = this.calculateSampledInternalForce(fomSolver, patch, u_full);

            // 2. Project directly to reduced space: F_red_int = M * f_sampled
            const F_red_int = new Float64Array(k);
            for (let i = 0; i < k; i++) {
                let sum = 0;
                for (let j = 0; j < this.m; j++) sum += this.M[i][j] * f_sampled[j];
                F_red_int[i] = sum;
            }

            // 3. Project Penalty Constraints (Stabilization)
            const F_penalty_full = new Float64Array(nDofs);
            fomSolver.applyPenaltyConstraints(null, F_penalty_full, u_full, patch, bcs);
            for (let i = 0; i < k; i++) {
                for (let d = 0; d < nDofs; d++) F_red_int[i] += PhiT.get(i, d) * F_penalty_full[d];
            }

            // 4. Compute residual
            const R_red = new Float64Array(k);
            let norm = 0;
            for (let i = 0; i < k; i++) {
                R_red[i] = F_ext_red_total[i] * loadFraction - F_red_int[i];
                norm += R_red[i] * R_red[i];
            }
            norm = Math.sqrt(norm);
            
            // DIAGNOSTIC: Compare F_red_int vs True Galerkin Projection
            const F_int_true = fomSolver.calculateInternalForce(patch, u_full);
            const F_red_int_true = new Float64Array(k);
            let diffNorm = 0, trueNorm = 0;
            for (let i = 0; i < k; i++) {
                for (let d = 0; d < nDofs; d++) F_red_int_true[i] += PhiT.get(i, d) * F_int_true[d];
                const diff = F_red_int_true[i] - F_red_int[i];
                diffNorm += diff * diff;
                trueNorm += F_red_int_true[i] * F_red_int_true[i];
            }
            const relErrForce = (Math.sqrt(trueNorm) > 1e-12) ? (Math.sqrt(diffNorm) / Math.sqrt(trueNorm) * 100) : 0;
            if (!quiet && (iter < 3 || iter === iterations - 1)) {
                console.log(`   [Iter ${iter}] |R_red|: ${norm.toExponential(3)} | Force Projection Error: ${relErrForce.toFixed(4)}%`);
                if (relErrForce > 50) console.warn(`   ⚠️ WARNING: Extreme force reconstruction error detected! Check interpolation conditioning.`);
            }

            residualHistory.push({ step: s, iter, norm });

            if (norm < tolerance && iter > 0) break;
            if (isNaN(norm) || norm > 1e25) break;

            const Kt_red_copy = Kt_red.map(row => Array.from(row));
            const dur = fomSolver.gaussianElimination(Kt_red_copy, Array.from(R_red));
            for (let i = 0; i < k; i++) ur[i] += dur[i];
        }
    }

    const u = new Float64Array(nDofs);
    for (let d = 0; d < nDofs; d++) {
        let sum = 0;
        for (let j = 0; j < k; j++) sum += Phi.get(d, j) * ur[j];
        u[d] = sum;
    }

    return { u, residualHistory, sampledCount: this.m, totalDofs: nDofs };
};

UDEIMEngine.prototype.updatePoints = function(romEngine, m, kf) {
    const { Matrix, SVD } = window.mlMatrix;
    this.m = m;
    this.kf = kf;
    
    // Subset U_f and indices
    const kf_clamped = Math.min(kf, this.U_f.columns);
    const U_f_subset = this.U_f.subMatrix(0, this.U_f.rows - 1, 0, kf_clamped - 1);
    const indices_subset = this.indices.slice(0, m);
    
    // Recompute Interpolation Matrix P^T U_f
    const PtU_arr = Array.from({ length: m }, () => new Float64Array(kf_clamped));
    for (let i = 0; i < m; i++) {
        for (let j = 0; j < kf_clamped; j++) PtU_arr[i][j] = U_f_subset.get(indices_subset[i], j);
    }
    const PtU_mat = new Matrix(PtU_arr.map(r => Array.from(r)));
    
    // Recompute Pseudo-Inverse
    const svd = new SVD(PtU_mat);
    const U_i = svd.leftSingularVectors;
    const S_i = svd.diagonal;
    const V_i = svd.rightSingularVectors;
    const S_inv = Matrix.zeros(V_i.columns, U_i.columns);
    const threshold = Math.max(...S_i) * 1e-12;
    for (let i = 0; i < S_i.length; i++) if (S_i[i] > threshold) S_inv.set(i, i, 1.0 / S_i[i]);
    const PtU_pinv = V_i.mmul(S_inv).mmul(U_i.transpose());

    // Recompute M = Phi^T * A * U_f_subset * PtU_pinv
    const Phi = romEngine.Phi;
    const nGlobalDofs = Phi.rows;
    const A_Uf = new Matrix(Array.from({length: nGlobalDofs}, () => new Array(kf_clamped).fill(0)));
    for (let j = 0; j < kf_clamped; j++) {
        for (let e = 0; e < this.numElements; e++) {
            const map = this.elementDofMap[e];
            for (let ld = 0; ld < this.numLocalDofs; ld++) {
                const globalDof = map[ld];
                A_Uf.set(globalDof, j, A_Uf.get(globalDof, j) + U_f_subset.get(e * this.numLocalDofs + ld, j));
            }
        }
    }
    
    const M_matrix = Phi.transpose().mmul(A_Uf).mmul(PtU_pinv);
    this.M = M_matrix.to2DArray();
    this.PtU_pinv = PtU_pinv.to2DArray();
};

UDEIMEngine.prototype.calculateSampledInternalForce = function(fomSolver, patch, u_full) {
    const gRule = window.GaussQuadrature2D.getPoints(Math.max(patch.p, patch.q) + 1);
    const activeElementForces = {};
    for (const el of this.activeElements) {
        activeElementForces[el.index] = this._computeLocalForce(fomSolver, patch, el, u_full, gRule);
    }
    const f_sampled = new Float64Array(this.m);
    for (let i = 0; i < this.m; i++) {
        const { e, ld } = this.sampledDofs[i];
        f_sampled[i] = activeElementForces[e][ld];
    }
    return f_sampled;
};

UDEIMEngine.prototype.reconstructFullForce = function(f_sampled) {
    const { Matrix } = window.mlMatrix;
    const f_sampled_mat = new Matrix([Array.from(f_sampled)]).transpose();
    
    if (!this.PtU_pinv) return new Float64Array(0);
    const PtU_pinv_mat = new Matrix(this.PtU_pinv);
    const c = PtU_pinv_mat.mmul(f_sampled_mat).to1DArray();
    
    const Phi = window.app.romEngine.Phi;
    const N = Phi.rows;
    const f_full = new Float64Array(N);
    
    const N_u = this.U_f.rows;
    const f_u = new Float64Array(N_u);
    for (let i = 0; i < N_u; i++) {
        for (let j = 0; j < this.kf; j++) {
            f_u[i] += this.U_f.get(i, j) * c[j];
        }
    }
    
    for (let e = 0; e < this.numElements; e++) {
        const map = this.elementDofMap[e];
        for (let ld = 0; ld < this.numLocalDofs; ld++) {
            const globalDof = map[ld];
            f_full[globalDof] += f_u[e * this.numLocalDofs + ld];
        }
    }
    
    return f_full;
};
