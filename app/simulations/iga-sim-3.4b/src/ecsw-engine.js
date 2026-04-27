/**
 * Phase 3.4b: ECSW (Energy-Conserving Sampling and Weighting)
 * 
 * Logic:
 * 1. Collect element-wise internal force snapshots: f_{e,s} = Phi^T * f_e(u_s)
 * 2. Solve NNLS: argmin || sum(w_e * f_{e,s}) - F_{int,s}^{reduced} ||_2  subject to w_e >= 0
 * 3. Result: A sparse set of elements with weights w_e.
 */
class ECSWEngine {
    constructor() {
        this.weights = null;      // Map of elementIndex -> weight
        this.sampleElements = []; // List of active element objects
        this.phi = null;          // Displacement basis
        this.Kp_red = null;       // Reduced Penalty Stiffness [k x k]
        this.nDofs = 0;
        this.k = 0;
    }

    /**
     * Train ECSW weights using element-wise force snapshots.
     * @param {IGANonlinearSolver} fomSolver 
     * @param {ROMEngine} romEngine 
     * @param {Object} patch 
     * @param {Array} snapshotDisplacements Array of Float64Arrays
     */
    async train(fomSolver, romEngine, patch, snapshotDisplacements) {
        this.phi = romEngine.Phi; // Matrix [nDofs x k]
        this.k = this.phi.columns;
        this.nDofs = this.phi.rows;
        
        const nSnaps = snapshotDisplacements.length;
        const { uniqueU, uniqueV, elements: allElements } = this._getAllElements(patch);
        const nElements = allElements.length;

        console.log(`ECSW: Training on ${nElements} elements across ${nSnaps} snapshots...`);

        // 0. Precompute all element projected forces (Big memory, but fast)
        // O(nElements * nSnaps * k)
        console.log("ECSW: Precomputing element force projections...");
        const elementProjections = allElements.map(el => 
            this._computeElementProjectedForce(fomSolver, patch, el, snapshotDisplacements)
        );

        // 1. Construct the Target Matrix B [k*nSnaps x 1]
        // B contains the projected assembly internal forces (EXCLUDING penalty)
        const B = new Float64Array(this.k * nSnaps);
        for (let s = 0; s < nSnaps; s++) {
            const F_int_assembly = fomSolver.calculateInternalForce(patch, snapshotDisplacements[s]);
            
            // Project to reduced space: F_red = Phi^T * F_int
            for (let i = 0; i < this.k; i++) {
                let dot = 0;
                for (let d = 0; d < this.nDofs; d++) dot += this.phi.get(d, i) * F_int_assembly[d];
                B[s * this.k + i] = dot;
            }
        }

        const weights = new Float64Array(nElements).fill(0);
        const { Matrix } = window.mlMatrix;

        // Build System Matrix A [k*nSnaps x nElements]
        const A_data = [];
        for (let row = 0; row < this.k * nSnaps; row++) {
            const a_row = new Float64Array(nElements);
            for (let e = 0; e < nElements; e++) a_row[e] = elementProjections[e][row];
            A_data.push(Array.from(a_row));
        }
        
        const A_mat = new Matrix(A_data);
        const B_mat = new Matrix([Array.from(B)]).transpose();

        const AT_A = A_mat.transpose().mmul(A_mat);
        const AT_B = A_mat.transpose().mmul(B_mat);
        
        const w = new Float64Array(nElements).fill(1.0); 
        for (let iter = 0; iter < 200; iter++) {
            let maxChange = 0;
            for (let i = 0; i < nElements; i++) {
                let grad = AT_B.get(i, 0);
                for (let j = 0; j < nElements; j++) grad -= AT_A.get(i, j) * w[j];
                const oldW = w[i];
                const diag = AT_A.get(i, i);
                if (diag > 1e-12) w[i] = Math.max(0, w[i] + grad / diag);
                maxChange = Math.max(maxChange, Math.abs(w[i] - oldW));
            }
            if (maxChange < 1e-8) break;
        }

        const activeIndices = [];
        for (let i = 0; i < nElements; i++) {
            if (w[i] > 1e-6) {
                weights[i] = w[i];
                activeIndices.push(i);
            }
        }

        const final_B = A_mat.mmul(new Matrix([Array.from(weights)]).transpose());
        const bNorm = Math.sqrt(B.reduce((a, b) => a + b*b, 0));
        let resNormFinal = 0;
        for(let i=0; i<B.length; i++) resNormFinal += (B[i] - final_B.get(i,0))**2;
        resNormFinal = Math.sqrt(resNormFinal);

        console.group("ECSW Training Debug");
        console.log(`Relative Training Error: ${(resNormFinal/bNorm*100).toFixed(6)}%`);
        console.log(`Weights: Count=${activeIndices.length}/${nElements}`);
        console.groupEnd();

        this.weights = weights;
        this.sampleElements = activeIndices.map(idx => ({
            ...allElements[idx],
            weight: weights[idx]
        }));

        // 3. Precompute Reduced Penalty Stiffness
        // Penalty is linear: F_p = K_p * u. We project K_p into reduced space.
        console.log("ECSW: Precomputing Reduced Penalty Matrix...");
        const Kp_full = Array.from({ length: this.nDofs }, () => new Float64Array(this.nDofs));
        // We use a dummy solver call to extract the penalty matrix
        fomSolver.applyPenaltyConstraints(Kp_full, null, new Float64Array(this.nDofs), patch);
        
        const Kp_mat = new Matrix(Kp_full);
        this.Kp_red = this.phi.transpose().mmul(Kp_mat).mmul(this.phi).to2DArray();

        console.log(`ECSW Trained: ${this.sampleElements.length} elements selected.`);
        return {
            elementCount: this.sampleElements.length,
            totalElements: nElements
        };
    }

    /**
     * Weighted Assembly Online.
     */
    assembleReducedForce(fomSolver, patch, ur) {
        const F_red_assembly = new Float64Array(this.k).fill(0);
        const F_red_penalty = new Float64Array(this.k).fill(0);
        const { p, q } = patch;
        const gRule = window.GaussQuadrature2D.getPoints(Math.max(p, q) + 1);

        const u_full = new Float64Array(this.nDofs);
        for (let d = 0; d < this.nDofs; d++) 
            for (let j = 0; j < this.k; j++) u_full[d] += this.phi.get(d, j) * ur[j];

        this.sampleElements.forEach(el => {
            const f_e_full = new Float64Array(this.nDofs).fill(0);
            this._assembleSingleElement(fomSolver, patch, el, u_full, f_e_full, gRule);
            
            for (let i = 0; i < this.k; i++) {
                let dot = 0;
                for (let d = 0; d < this.nDofs; d++) {
                    if (f_e_full[d] !== 0) dot += this.phi.get(d, i) * f_e_full[d];
                }
                F_red_assembly[i] += el.weight * dot;
            }
        });

        for (let i = 0; i < this.k; i++) {
            for (let j = 0; j < this.k; j++) {
                F_red_penalty[i] += this.Kp_red[i][j] * ur[j];
            }
        }

        const normA = Math.sqrt(F_red_assembly.reduce((a, b) => a + b*b, 0));
        const normP = Math.sqrt(F_red_penalty.reduce((a, b) => a + b*b, 0));
        if (normA > 1e12 || normP > 1e12) {
            console.warn(`ECSW: Large forces detected! Assembly=${normA.toExponential(2)}, Penalty=${normP.toExponential(2)}`);
        }

        return { 
            F_red: F_red_assembly.map((f, i) => f + F_red_penalty[i]),
            normAssembly: normA,
            normPenalty: normP
        };
    }

    solveReduced(fomSolver, patch, loads, options = {}) {
        const { steps = 10, iterations = 15, tolerance = 1e-6 } = options;
        const k = this.k;
        const nDofs = this.nDofs;
        const nV = patch.controlPoints[0].length;

        let ur = new Float64Array(k);
        const residualHistory = [];

        // Build projected external force
        const F_ext_total = new Float64Array(nDofs);
        loads.forEach(l => {
            const idx = (l.i * nV + l.j) * 2;
            F_ext_total[idx] += l.fx;
            F_ext_total[idx+1] += l.fy;
        });
        const F_ext_red = new Float64Array(k);
        for (let i = 0; i < k; i++) {
            for (let d = 0; d < nDofs; d++) F_ext_red[i] += this.phi.get(d, i) * F_ext_total[d];
        }

        const PhiT = this.phi.transpose();

        for (let s = 1; s <= steps; s++) {
            const frac = s / steps;
            
            for (let iter = 0; iter < iterations; iter++) {
                // --- Tangent Update (Full Newton-Raphson) ---
                const u_full_curr = new Float64Array(nDofs);
                for (let d = 0; d < nDofs; d++) 
                    for (let j = 0; j < k; j++) u_full_curr[d] += this.phi.get(d, j) * ur[j];
                
                const Kt_full = fomSolver.calculateTangentStiffness(patch, u_full_curr);
                fomSolver.applyPenaltyConstraints(Kt_full, null, u_full_curr, patch);
                
                const nV = patch.controlPoints[0].length;
                const penalty = 1e12;
                for (let j = 0; j < nV; j++) {
                    const idx = (0 * nV + j) * 2;
                    Kt_full[idx][idx] += penalty;
                    Kt_full[idx+1][idx+1] += penalty;
                }

                const Kt_red_mat = PhiT.mmul(new window.mlMatrix.Matrix(Kt_full)).mmul(this.phi);
                const Kt_red = Kt_red_mat.to2DArray();

                // 1. Current state evaluation
                const state = this.assembleReducedForce(fomSolver, patch, ur);
                let R_red = F_ext_red.map((f, i) => f * frac - state.F_red[i]);
                let norm = Math.sqrt(R_red.reduce((a, b) => a + b*b, 0));
                
                if (iter === 0 || iter % 5 === 0) {
                    console.log(`Step ${s} Iter ${iter}: Res=${norm.toExponential(2)}, F_asm=${state.normAssembly.toExponential(2)}`);
                }

                if (norm < tolerance && iter > 0) break;
                if (isNaN(norm) || norm > 1e25) break;

                // 2. Solve for increment
                const Kt_copy = Kt_red.map(row => new Float64Array(row));
                let dur = fomSolver.gaussianElimination(Kt_copy, Array.from(R_red));
                
                // 3. BACKTRACKING LINE SEARCH
                let damping = 0.8;
                let stepAccepted = false;
                
                for (let backtrack = 0; backtrack < 5; backtrack++) {
                    const next_ur = new Float64Array(k);
                    for (let i = 0; i < k; i++) next_ur[i] = ur[i] + dur[i] * damping;
                    
                    const next_state = this.assembleReducedForce(fomSolver, patch, next_ur);
                    const next_R = F_ext_red.map((f, i) => f * frac - next_state.F_red[i]);
                    const next_norm = Math.sqrt(next_R.reduce((a, b) => a + b*b, 0));
                    
                    if (next_norm < norm || iter === 0) {
                        // Step improves residual (or it's the very first step)
                        ur = next_ur;
                        stepAccepted = true;
                        break;
                    } else {
                        // Step made it worse! Backtrack.
                        damping *= 0.5;
                        if (damping < 0.05) break;
                    }
                }

                if (!stepAccepted) {
                    console.warn(`ECSW: Line search failed to find descent direction at Step ${s} Iter ${iter}`);
                    // Still take a tiny step to try and escape
                    for (let i = 0; i < k; i++) ur[i] += dur[i] * 0.01;
                }

                residualHistory.push({ step: s, iter, norm });
            }
        }

        const u = new Float64Array(nDofs);
        for (let d = 0; d < nDofs; d++) 
            for (let j = 0; j < k; j++) u[d] += this.phi.get(d, j) * ur[j];

        return { u, ur, residualHistory };
    }

    _getAllElements(patch) {
        const { U, V } = patch;
        const uniqueU = [...new Set(U)], uniqueV = [...new Set(V)];
        const elements = [];
        for (let i = 0; i < uniqueU.length - 1; i++) {
            for (let j = 0; j < uniqueV.length - 1; j++) {
                elements.push({ 
                    i, j, 
                    uMin: uniqueU[i], uMax: uniqueU[i+1],
                    vMin: uniqueV[j], vMax: uniqueV[j+1]
                });
            }
        }
        return { uniqueU, uniqueV, elements };
    }

    _computeElementProjectedForce(fomSolver, patch, el, snapshotDisplacements) {
        const nSnaps = snapshotDisplacements.length;
        const Ge = new Float64Array(this.k * nSnaps);
        const gRule = window.GaussQuadrature2D.getPoints(Math.max(patch.p, patch.q) + 1);

        for (let s = 0; s < nSnaps; s++) {
            const f_e_full = new Float64Array(this.nDofs);
            this._assembleSingleElement(fomSolver, patch, el, snapshotDisplacements[s], f_e_full, gRule);
            
            for (let i = 0; i < this.k; i++) {
                let dot = 0;
                for (let d = 0; d < this.nDofs; d++) {
                    if (f_e_full[d] !== 0) dot += this.phi.get(d, i) * f_e_full[d];
                }
                Ge[s * this.k + i] = dot;
            }
        }
        return Ge;
    }

    _assembleSingleElement(fomSolver, patch, el, u_disp, f_out, gRule) {
        const { uMin, uMax, vMin, vMax } = el;
        for (let gu = 0; gu < gRule.points.length; gu++) {
            const u = ((uMax - uMin) * gRule.points[gu] + (uMax + uMin)) / 2;
            const wu = gRule.weights[gu] * (uMax - uMin) / 2;
            for (let gv = 0; gv < gRule.points.length; gv++) {
                const v = ((vMax - vMin) * gRule.points[gv] + (vMax + vMin)) / 2;
                const wv = gRule.weights[gv] * (vMax - vMin) / 2;

                const deriv = fomSolver.engine.getSurfaceDerivatives(patch, u, v);
                const { grads: B_param, detJ, activeIndices } = fomSolver.getBParametric(patch, u, v, deriv);
                
                let dudx = 0, dudy = 0, dvdx = 0, dvdy = 0;
                for (let a = 0; a < activeIndices.length; a++) {
                    const k = activeIndices[a];
                    dudx += B_param[k][0] * u_disp[k * 2];
                    dudy += B_param[k][1] * u_disp[k * 2];
                    dvdx += B_param[k][0] * u_disp[k * 2 + 1];
                    dvdy += B_param[k][1] * u_disp[k * 2 + 1];
                }

                const Exx = dudx + 0.5 * (dudx*dudx + dvdx*dvdx);
                const Eyy = dvdy + 0.5 * (dudy*dudy + dvdy*dvdy);
                const Exy2 = (dudy + dvdx) + (dudx*dudy + dvdx*dvdy);
                
                const D = fomSolver.getPlaneStressD();
                const Sxx = D[0][0]*Exx + D[0][1]*Eyy;
                const Syy = D[1][0]*Exx + D[1][1]*Eyy;
                const Sxy = D[2][2]*Exy2;

                const factor = detJ * wu * wv * fomSolver.thickness;

                for (let a = 0; a < activeIndices.length; a++) {
                    const k = activeIndices[a];
                    const dRdx = B_param[k][0], dRdy = B_param[k][1];
                    
                    const bexx_u = (1 + dudx) * dRdx, bexx_v = (dvdx) * dRdx;
                    const beyy_u = (dudy) * dRdy, beyy_v = (1 + dvdy) * dRdy;
                    const bexy_u = (1 + dudx)*dRdy + dudy*dRdx, bexy_v = (1 + dvdy)*dRdx + dvdx*dRdy;

                    f_out[k * 2] += (bexx_u * Sxx + beyy_u * Syy + bexy_u * Sxy) * factor;
                    f_out[k * 2 + 1] += (bexx_v * Sxx + beyy_v * Syy + bexy_v * Sxy) * factor;
                }
            }
        }
    }
}

window.ECSWEngine = ECSWEngine;
