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
    async train(fomSolver, romEngine, patch, snapshotDisplacements, tolerance = 1e-4) {
        this.phi = romEngine.Phi; // Matrix [nDofs x k]
        this.k = this.phi.columns;
        this.nDofs = this.phi.rows;
        
        const nSnaps = snapshotDisplacements.length;
        const { uniqueU, uniqueV, elements: allElements } = this._getAllElements(patch);
        const nElements = allElements.length;

        console.log(`ECSW: Training on ${nElements} elements across ${nSnaps} snapshots with tol=${tolerance}...`);

        // 0. Precompute all element projected forces (Big memory, but fast)
        // O(nElements * nSnaps * k)
        console.log("ECSW: Precomputing element force projections...");
        const elementProjections = allElements.map(el => 
            this._computeElementProjectedForce(fomSolver, patch, el, snapshotDisplacements)
        );

        // 1. Construct the Target Matrix B [k*nSnaps x 1]
        // B contains the projected assembly internal forces (EXCLUDING penalty)
        const B = new Float64Array(this.k * nSnaps);
        const constrainedDofs = fomSolver._getConstrainedDofs ? fomSolver._getConstrainedDofs() : [];
        for (let s = 0; s < nSnaps; s++) {
            const F_int_assembly = fomSolver.calculateInternalForce(patch, snapshotDisplacements[s]);
            constrainedDofs.forEach(d => F_int_assembly[d] = 0);
            
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

        // --- NEW AUDIT LOGS ---
        let normB = 0;
        for (let r=0; r<B_mat.rows; r++) normB += B_mat.get(r,0)**2;
        normB = Math.sqrt(normB);
        console.log(`   [ECSW Audit] Target ||B||: ${normB.toExponential(4)}`);
        console.log(`   [ECSW Audit] Matrix A    : ${this.k * nSnaps}x${nElements}`);
        // ----------------------

        // Greedy ECSW (Non-Negative Orthogonal Matching Pursuit)
        const w = new Float64Array(nElements).fill(0.0);
        const activeSet = new Set();
        const residual = new Float64Array(B); // Initially B
        
        for (let step = 0; step < nElements; step++) {
            // 1. Find element that maximizes projection of residual
            let best_e = -1;
            let max_proj = -1e-12;
            for (let e = 0; e < nElements; e++) {
                if (activeSet.has(e)) continue;
                let proj = 0;
                for (let r = 0; r < this.k * nSnaps; r++) proj += A_mat.get(r, e) * residual[r];
                if (proj > max_proj) {
                    max_proj = proj;
                    best_e = e;
                }
            }
            
            if (best_e === -1) {
                console.log(`[Greedy ECSW] Stopping: No more elements improve the residual.`);
                break;
            }
            
            activeSet.add(best_e);
            const activeArr = Array.from(activeSet);
            
            // 2. Solve unconstrained NNLS on active set using tight coordinate descent
            for (let inner = 0; inner < 500; inner++) {
                let maxChange = 0;
                for (let i = 0; i < activeArr.length; i++) {
                    const e = activeArr[i];
                    let pred = 0;
                    for (let j = 0; j < activeArr.length; j++) pred += AT_A.get(e, activeArr[j]) * w[activeArr[j]];
                    const grad = AT_B.get(e, 0) - pred;
                    const diag = AT_A.get(e, e);
                    const oldW = w[e];
                    if (diag > 1e-12) w[e] = Math.max(0, oldW + grad / diag);
                    maxChange = Math.max(maxChange, Math.abs(w[e] - oldW));
                }
                if (maxChange < 1e-8) break; // Converged inner loop
            }
            
            // 3. Update global residual
            let resNorm = 0;
            for (let r = 0; r < this.k * nSnaps; r++) {
                let p = 0;
                for (let i = 0; i < activeArr.length; i++) p += A_mat.get(r, activeArr[i]) * w[activeArr[i]];
                residual[r] = B[r] - p;
                resNorm += residual[r]*residual[r];
            }
            resNorm = Math.sqrt(resNorm);
            const relError = resNorm / normB;
            
            console.log(`   [Greedy Step ${step+1}] Added El: ${best_e}, Active: ${activeSet.size}, RelError: ${(relError*100).toFixed(4)}%`);
            
            if (relError < tolerance) {
                console.log(`[Greedy ECSW] Converged to tolerance ${tolerance} with ${activeSet.size} elements.`);
                break;
            }
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
        this.activeElements = activeIndices.map(idx => ({
            ...allElements[idx],
            weight: weights[idx]
        }));
        this.sampleElements = this.activeElements;

        // Generate dummy history for Explorer (mapping element index to a representative DOF)
        const nV = patch.controlPoints[0].length;
        const p = patch.p, q = patch.q;
        this.history = this.activeElements.map((el, step) => {
            // Find a DOF in this span
            const uMid = (el.uMin + el.uMax) / 2;
            const vMid = (el.vMin + el.vMax) / 2;
            const spanU = fomSolver.engine.findSpan(patch.controlPoints.length - 1, p, uMid, patch.U);
            const spanV = fomSolver.engine.findSpan(nV - 1, q, vMid, patch.V);
            const cpIdx = spanU * nV + spanV;
            return {
                step: step + 1,
                point: cpIdx * 2, // Representative X-DOF
                residual: new Array(this.nDofs).fill(0),
                maxVal: el.weight
            };
        });

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
     * Assembles BOTH reduced force and reduced tangent stiffness directly.
     */
    assembleReducedSystem(fomSolver, patch, ur) {
        const F_red_assembly = new Float64Array(this.k).fill(0);
        const Kt_red_assembly = Array.from({ length: this.k }, () => new Float64Array(this.k).fill(0));
        
        const { p, q } = patch;
        const gRule = window.GaussQuadrature2D.getPoints(Math.max(p, q) + 1);

        // Pre-expand ur to full displacement ONLY ONCE per assembly
        const u_full = new Float64Array(this.nDofs);
        for (let j = 0; j < this.k; j++) {
            const val = ur[j];
            if (val === 0) continue;
            for (let d = 0; d < this.nDofs; d++) u_full[d] += this.phi.get(d, j) * val;
        }

        this.sampleElements.forEach(el => {
            const { f_e, k_e, activeDofs } = this._assembleElementPhysics(fomSolver, patch, el, u_full, gRule);
            
            const nLocal = activeDofs.length;
            const weight = el.weight;

            // F_red += w * Phi_e^T * f_e
            for (let i = 0; i < this.k; i++) {
                let dot = 0;
                for (let a = 0; a < nLocal; a++) {
                    dot += this.phi.get(activeDofs[a], i) * f_e[a];
                }
                F_red_assembly[i] += weight * dot;
            }

            // Kt_red += w * Phi_e^T * k_e * Phi_e
            const Phi_e = Array.from({ length: nLocal }, () => new Float64Array(this.k));
            for (let a = 0; a < nLocal; a++) {
                for (let i = 0; i < this.k; i++) Phi_e[a][i] = this.phi.get(activeDofs[a], i);
            }

            const temp = Array.from({ length: nLocal }, () => new Float64Array(this.k));
            for (let a = 0; a < nLocal; a++) {
                for (let i = 0; i < this.k; i++) {
                    let val = 0;
                    for (let b = 0; b < nLocal; b++) val += k_e[a][b] * Phi_e[b][i];
                    temp[a][i] = val;
                }
            }

            for (let i = 0; i < this.k; i++) {
                for (let j = 0; j < this.k; j++) {
                    let val = 0;
                    for (let a = 0; a < nLocal; a++) val += Phi_e[a][i] * temp[a][j];
                    Kt_red_assembly[i][j] += weight * val;
                }
            }
        });

        // Add precomputed Reduced Penalty
        for (let i = 0; i < this.k; i++) {
            for (let j = 0; j < this.k; j++) {
                F_red_assembly[i] += this.Kp_red[i][j] * ur[j];
                Kt_red_assembly[i][j] += this.Kp_red[i][j];
            }
        }

        return { F_red: F_red_assembly, Kt_red: Kt_red_assembly };
    }

    solveReduced(fomSolver, patch, loads, options = {}) {
        const { steps = 10, iterations = 15, tolerance = 1e-6 } = options;
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
            for (let d = 0; d < nDofs; d++) F_ext_red[i] += this.phi.get(d, i) * F_ext_total[d];
        }

        // Timing diagnostics
        let totalAssemblyTime = 0, totalSolveTime = 0, totalIters = 0;

        for (let s = 1; s <= steps; s++) {
            const frac = s / steps;
            
            for (let iter = 0; iter < iterations; iter++) {
                // 1. Online Assembly (Reduced only!)
                const t0 = performance.now();
                const { F_red, Kt_red } = this.assembleReducedSystem(fomSolver, patch, ur);
                totalAssemblyTime += performance.now() - t0;
                
                // 2. Residual
                const R_red = F_ext_red.map((f, i) => f * frac - F_red[i]);
                const norm = Math.sqrt(R_red.reduce((a, b) => a + b*b, 0));

                if (norm < tolerance && iter > 0) break;
                if (isNaN(norm) || norm > 1e25) break;

                // 3. Solve for increment dur = Kt_red^-1 * R_red [k x k]
                const t1 = performance.now();
                const dur = fomSolver.gaussianElimination(Kt_red, Array.from(R_red));
                totalSolveTime += performance.now() - t1;
                
                // 4. Update
                for (let i = 0; i < k; i++) ur[i] += dur[i];

                residualHistory.push({ step: s, iter, norm });
                totalIters++;
            }
        }

        // Performance Report
        console.log(`%cECSW Online Performance: ${totalIters} iters | Assembly: ${totalAssemblyTime.toFixed(1)}ms | Solve(${k}x${k}): ${totalSolveTime.toFixed(1)}ms | Elements: ${this.sampleElements.length}`, 'color:#10b981;font-weight:bold');

        const u = new Float64Array(nDofs);
        for (let j = 0; j < k; j++) {
            const val = ur[j];
            for (let d = 0; d < nDofs; d++) u[d] += this.phi.get(d, j) * val;
        }

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
            const { f_e, activeDofs } = this._assembleElementPhysics(fomSolver, patch, el, snapshotDisplacements[s], gRule);
            
            for (let i = 0; i < this.k; i++) {
                let dot = 0;
                for (let a = 0; a < activeDofs.length; a++) {
                    dot += this.phi.get(activeDofs[a], i) * f_e[a];
                }
                Ge[s * this.k + i] = dot;
            }
        }
        return Ge;
    }

    _assembleElementPhysics(fomSolver, patch, el, u_disp, gRule) {
        const { uMin, uMax, vMin, vMax } = el;
        
        // Find local active DOFs for this element using a center-point probe
        const derivCenter = fomSolver.engine.getSurfaceDerivatives(patch, (uMin+uMax)/2, (vMin+vMax)/2);
        const { activeIndices } = fomSolver.getBParametric(patch, (uMin+uMax)/2, (vMin+vMax)/2, derivCenter);
        
        const nLocalBasis = activeIndices.length;
        const nLocalDof = nLocalBasis * 2;
        const f_e = new Float64Array(nLocalDof);
        const k_e = Array.from({ length: nLocalDof }, () => new Float64Array(nLocalDof));
        const globalDofMap = new Uint32Array(nLocalDof);
        for(let i=0; i<nLocalBasis; i++) {
            globalDofMap[i*2] = activeIndices[i] * 2;
            globalDofMap[i*2+1] = activeIndices[i] * 2 + 1;
        }

        const D = fomSolver.getPlaneStressD();

        for (let gu = 0; gu < gRule.points.length; gu++) {
            const u = ((uMax - uMin) * gRule.points[gu] + (uMax + uMin)) / 2;
            const wu = gRule.weights[gu] * (uMax - uMin) / 2;
            for (let gv = 0; gv < gRule.points.length; gv++) {
                const v = ((vMax - vMin) * gRule.points[gv] + (vMax + vMin)) / 2;
                const wv = gRule.weights[gv] * (vMax - vMin) / 2;

                const deriv = fomSolver.engine.getSurfaceDerivatives(patch, u, v);
                const { grads: B_param, detJ } = fomSolver.getBParametric(patch, u, v, deriv);
                
                let dudx = 0, dudy = 0, dvdx = 0, dvdy = 0;
                for (let a = 0; a < nLocalBasis; a++) {
                    const k = activeIndices[a];
                    const dRdx = B_param[k][0], dRdy = B_param[k][1];
                    dudx += dRdx * u_disp[k * 2];
                    dudy += dRdy * u_disp[k * 2];
                    dvdx += dRdx * u_disp[k * 2 + 1];
                    dvdy += dRdy * u_disp[k * 2 + 1];
                }

                const Exx = dudx + 0.5 * (dudx*dudx + dvdx*dvdx);
                const Eyy = dvdy + 0.5 * (dudy*dudy + dvdy*dvdy);
                const Exy2 = (dudy + dvdx) + (dudx*dudy + dvdx*dvdy);
                
                const Sxx = D[0][0]*Exx + D[0][1]*Eyy;
                const Syy = D[1][0]*Exx + D[1][1]*Eyy;
                const Sxy = D[2][2]*Exy2;

                const factor = detJ * wu * wv * fomSolver.thickness;

                const B_NL = Array.from({ length: nLocalBasis }, () => [new Float64Array(2), new Float64Array(2), new Float64Array(2)]);
                for (let a = 0; a < nLocalBasis; a++) {
                    const k = activeIndices[a];
                    const dRdx = B_param[k][0], dRdy = B_param[k][1];
                    const bexx_u = (1 + dudx) * dRdx, bexx_v = (dvdx) * dRdx;
                    const beyy_u = (dudy) * dRdy, beyy_v = (1 + dvdy) * dRdy;
                    const bexy_u = (1 + dudx)*dRdy + dudy*dRdx, bexy_v = (1 + dvdy)*dRdx + dvdx*dRdy;

                    f_e[a * 2] += (bexx_u * Sxx + beyy_u * Syy + bexy_u * Sxy) * factor;
                    f_e[a * 2 + 1] += (bexx_v * Sxx + beyy_v * Syy + bexy_v * Sxy) * factor;

                    B_NL[a][0][0] = bexx_u; B_NL[a][0][1] = bexx_v;
                    B_NL[a][1][0] = beyy_u; B_NL[a][1][1] = beyy_v;
                    B_NL[a][2][0] = bexy_u; B_NL[a][2][1] = bexy_v;
                }

                for (let a = 0; a < nLocalBasis; a++) {
                    for (let b = 0; b < nLocalBasis; b++) {
                        for (let i = 0; i < 2; i++) {
                            for (let j = 0; j < 2; j++) {
                                let kab = 0;
                                for (let r = 0; r < 3; r++) {
                                    for (let c = 0; c < 3; c++) kab += B_NL[a][r][i] * D[r][c] * B_NL[b][c][j];
                                }
                                k_e[a * 2 + i][b * 2 + j] += kab * factor;
                            }
                        }
                    }
                }

                for (let a = 0; a < nLocalBasis; a++) {
                    for (let b = 0; b < nLocalBasis; b++) {
                        const dRdx_a = B_param[activeIndices[a]][0], dRdy_a = B_param[activeIndices[a]][1];
                        const dRdx_b = B_param[activeIndices[b]][0], dRdy_b = B_param[activeIndices[b]][1];
                        const k_geo = (dRdx_a * Sxx * dRdx_b + dRdy_a * Syy * dRdy_b + dRdx_a * Sxy * dRdy_b + dRdy_a * Sxy * dRdx_b) * factor;
                        k_e[a * 2][b * 2] += k_geo;
                        k_e[a * 2 + 1][b * 2 + 1] += k_geo;
                    }
                }
            }
        }
        return { f_e, k_e, activeDofs: globalDofMap };
    }
}

window.ECSWEngine = ECSWEngine;
