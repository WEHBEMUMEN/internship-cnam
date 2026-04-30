/**
 * ECSW Core Math Engine
 * Implements the Energy-Conserving Sampling and Weighting (ECSW) logic.
 */
window.ECSWCore = {
    /**
     * Solve Non-Negative Least Squares (NNLS)
     * argmin || A*w - B ||_2  subject to w >= 0
     * uses an iterative active-set / gradient descent hybrid for browser performance.
     */
    solveNNLS: function(A_mat, B_mat, options = {}) {
        const { tolerance = 1e-4, maxIters = 500 } = options;
        const nElements = A_mat.columns;
        const { Matrix } = window.mlMatrix;

        console.log(`ECSW Core: Solving NNLS for ${nElements} elements...`);

        // Compute AT_A and AT_B for the normal equations
        const AT = A_mat.transpose();
        const AT_A = AT.mmul(A_mat);
        const AT_B = AT.mmul(B_mat);

        const w = new Float64Array(nElements).fill(1.0); // Start with uniform weights
        
        for (let iter = 0; iter < maxIters; iter++) {
            let maxChange = 0;
            for (let i = 0; i < nElements; i++) {
                // Gradient of 0.5 * ||Aw - B||^2 is AT*(Aw - B)
                let grad = AT_B.get(i, 0);
                for (let j = 0; j < nElements; j++) {
                    grad -= AT_A.get(i, j) * w[j];
                }

                const oldW = w[i];
                const diag = AT_A.get(i, i);
                
                // Projected Gradient step
                if (diag > 1e-12) {
                    w[i] = Math.max(0, w[i] + grad / diag);
                }
                
                maxChange = Math.max(maxChange, Math.abs(w[i] - oldW));
            }
            
            if (maxChange < tolerance) {
                console.log(`ECSW Core: NNLS converged in ${iter} iterations.`);
                break;
            }
        }

        return w;
    },

    /**
     * Project element-wise forces into the reduced basis Phi.
     * Ge,s = Phi_e^T * f_e(u_s)
     */
    computeElementProjectedForce: function(fomSolver, patch, el, snapshotDisplacements, Phi) {
        const nSnaps = snapshotDisplacements.length;
        const k = Phi.columns;
        const Ge = new Float64Array(k * nSnaps);
        const gRule = window.GaussQuadrature2D.getPoints(Math.max(patch.p, patch.q) + 1);

        for (let s = 0; s < nSnaps; s++) {
            const { f_e, activeDofs } = this._assembleElementPhysics(fomSolver, patch, el, snapshotDisplacements[s], gRule);
            
            for (let i = 0; i < k; i++) {
                let dot = 0;
                for (let a = 0; a < activeDofs.length; a++) {
                    dot += Phi.get(activeDofs[a], i) * f_e[a];
                }
                Ge[s * k + i] = dot;
            }
        }
        return Ge;
    },

    /**
     * Helper: Assemble local element physics (force and tangent)
     */
    _assembleElementPhysics: function(fomSolver, patch, el, u_disp, gRule) {
        const { uMin, uMax, vMin, vMax } = el;
        
        // Probe center for active indices
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
                    const k_idx = activeIndices[a];
                    const dRdx = B_param[k_idx][0], dRdy = B_param[k_idx][1];
                    dudx += dRdx * u_disp[k_idx * 2];
                    dudy += dRdy * u_disp[k_idx * 2];
                    dvdx += dRdx * u_disp[k_idx * 2 + 1];
                    dvdy += dRdy * u_disp[k_idx * 2 + 1];
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
                    const k_idx = activeIndices[a];
                    const dRdx = B_param[k_idx][0], dRdy = B_param[k_idx][1];
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
};
