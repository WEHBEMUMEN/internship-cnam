/**
 * Phase 3.4b: U-DEIM Online Solver
 */

UDEIMEngine.prototype.solveReduced = function(fomSolver, romEngine, patch, bcs, loads, options = {}) {
    const { iterations = 15, tolerance = 1e-6, steps = 1 } = options;
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

            // Recompute full tangent per iteration for nonlinear stability (Full Newton-Raphson)
            const Kt_full = fomSolver.calculateTangentStiffness(patch, u_full);
            const Kt_red = PhiT.mmul(new Matrix(Kt_full.map(r => Array.from(r)))).mmul(Phi).to2DArray();


            // 1. Evaluate only active elements
            const activeElementForces = {};
            for (const el of this.activeElements) {
                activeElementForces[el.index] = this._computeLocalForce(fomSolver, patch, el, u_full, gRule);
            }

            // 2. Extract sampled DOFs
            const f_sampled = new Float64Array(this.m);
            for (let i = 0; i < this.m; i++) {
                const { e, ld } = this.sampledDofs[i];
                f_sampled[i] = activeElementForces[e][ld];
            }

            // 3. Project directly to reduced space: F_red_int = M * f_sampled
            const F_red_int = new Float64Array(k);
            for (let i = 0; i < k; i++) {
                let sum = 0;
                for (let j = 0; j < this.m; j++) sum += this.M[i][j] * f_sampled[j];
                F_red_int[i] = sum;
            }

            // Penalty constraints removed from internal force to prevent Projection Explosion
            // const F_penalty_full = new Float64Array(nDofs);
            // fomSolver.applyPenaltyConstraints(null, F_penalty_full, u_full, patch, bcs);
            // for (let i = 0; i < k; i++) {
            //     for (let d = 0; d < nDofs; d++) F_red_int[i] += PhiT.get(i, d) * F_penalty_full[d];
            // }



            // 4. Compute residual
            const R_red = new Float64Array(k);
            let norm = 0;
            for (let i = 0; i < k; i++) {
                R_red[i] = F_ext_red_total[i] * loadFraction - F_red_int[i];
                norm += R_red[i] * R_red[i];
            }
            norm = Math.sqrt(norm);
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

    return { u, ur, residualHistory, sampledCount: this.activeElements.length };
};
