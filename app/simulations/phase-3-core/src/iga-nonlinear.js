/**
 * 2D Isogeometric Analysis Nonlinear Solver (Phase 3 Core)
 * Handles Large Deflections (Geometric Nonlinearity) using Newton-Raphson.
 */

// Conditional guard: only define if not already loaded by phase-2-core
if (typeof GaussQuadrature2D === 'undefined') {
    class GaussQuadrature2D {
        static getPoints(n = 3) {
            if (n === 2) {
                const p = 1.0 / Math.sqrt(3);
                return { points: [-p, p], weights: [1, 1] };
            }
            if (n === 3) {
                return { points: [-Math.sqrt(0.6), 0, Math.sqrt(0.6)], weights: [5/9, 8/9, 5/9] };
            }
            const p1 = Math.sqrt((3 - 2 * Math.sqrt(1.2)) / 7);
            const p2 = Math.sqrt((3 + 2 * Math.sqrt(1.2)) / 7);
            const w1 = (18 + Math.sqrt(30)) / 36;
            const w2 = (18 - Math.sqrt(30)) / 36;
            return { points: [-p2, -p1, p1, p2], weights: [w2, w1, w1, w2] };
        }
    }
    window.GaussQuadrature2D = GaussQuadrature2D;
}

class IGANonlinearSolver {
    constructor(nurbsEngine) {
        this.engine = nurbsEngine;
        this.E = 100000; // Benchmark (MPa)
        this.nu = 0.3;
        this.thickness = 1.0;
    }

    getPlaneStressD() {
        const factor = this.E / (1 - this.nu * this.nu);
        return [
            [factor, factor * this.nu, 0],
            [factor * this.nu, factor, 0],
            [0, 0, factor * (1 - this.nu) / 2]
        ];
    }

    getBParametric(patch, u, v, deriv) {
        const { p, q, U, V, weights, controlPoints } = patch;
        const nU = controlPoints.length;
        const nV = controlPoints[0].length;
        const J = [[deriv.dU.x, deriv.dV.x], [deriv.dU.y, deriv.dV.y]];
        let detJ_2D = J[0][0] * J[1][1] - J[0][1] * J[1][0];
        if (isNaN(detJ_2D) || Math.abs(detJ_2D) < 1e-12) detJ_2D = (detJ_2D >= 0) ? 1e-12 : -1e-12;
        const J_inv = [[J[1][1]/detJ_2D, -J[0][1]/detJ_2D], [-J[1][0]/detJ_2D, J[0][0]/detJ_2D]];
        
        const W = deriv.W, Wu = deriv.Wu, Wv = deriv.Wv;
        const grads = Array(nU * nV).fill(0).map(() => [0, 0]);
        const activeIndices = [];

        const spanU = this.engine.findSpan(nU - 1, p, u, U);
        const spanV = this.engine.findSpan(nV - 1, q, v, V);
        const dersU = this.engine.basisFunsDerivs(spanU, u, p, U, 1);
        const dersV = this.engine.basisFunsDerivs(spanV, v, q, V, 1);

        for (let i = 0; i <= p; i++) {
            const Ni = dersU[0][i], dNi = dersU[1][i];
            const cpI = spanU - p + i;
            for (let j = 0; j <= q; j++) {
                const Mj = dersV[0][j], dMj = dersV[1][j];
                const cpJ = spanV - q + j;
                const w = weights[cpI][cpJ];
                
                const dRdu = ((dNi * Mj * w) * W - (Ni * Mj * w) * Wu) / (W * W);
                const dRdv = ((Ni * dMj * w) * W - (Ni * Mj * w) * Wv) / (W * W);
                const dRdx = J_inv[0][0] * dRdu + J_inv[1][0] * dRdv;
                const dRdy = J_inv[0][1] * dRdu + J_inv[1][1] * dRdv;
                
                const idx = cpI * nV + cpJ;
                grads[idx] = [dRdx, dRdy];
                activeIndices.push(idx);
            }
        }
        return { grads, detJ: detJ_2D, activeIndices };
    }

    calculateInternalForce(patch, u_disp) {
        const { p, q, U, V, controlPoints } = patch;
        const nBasisU = controlPoints.length;
        const nBasisV = controlPoints[0].length;
        const nDofs = nBasisU * nBasisV * 2;
        const F_int = new Float64Array(nDofs).fill(0);

        const uniqueU = [...new Set(U)], uniqueV = [...new Set(V)];
        const gRule = GaussQuadrature2D.getPoints(Math.max(p, q) + 1);

        for (let i = 0; i < uniqueU.length - 1; i++) {
            const uMin = uniqueU[i], uMax = uniqueU[i+1];
            if (uMax - uMin < 1e-10) continue;
            for (let j = 0; j < uniqueV.length - 1; j++) {
                const vMin = uniqueV[j], vMax = uniqueV[j+1];
                if (vMax - vMin < 1e-10) continue;

                for (let gu = 0; gu < gRule.points.length; gu++) {
                    const u = ((uMax - uMin) * gRule.points[gu] + (uMax + uMin)) / 2;
                    const wu = gRule.weights[gu] * (uMax - uMin) / 2;
                    for (let gv = 0; gv < gRule.points.length; gv++) {
                        const v = ((vMax - vMin) * gRule.points[gv] + (vMax + vMin)) / 2;
                        const wv = gRule.weights[gv] * (vMax - vMin) / 2;

                        // Consistency: use same J calculations
                        const deriv = this.engine.getSurfaceDerivatives(patch, u, v);
                        const { grads: B_param, detJ, activeIndices } = this.getBParametric(patch, u, v, deriv);
                        
                        let dudx = 0, dudy = 0, dvdx = 0, dvdy = 0;
                        for (let a = 0; a < activeIndices.length; a++) {
                            const k = activeIndices[a];
                            const dRdx = B_param[k][0], dRdy = B_param[k][1];
                            dudx += dRdx * u_disp[k * 2];
                            dudy += dRdy * u_disp[k * 2];
                            dvdx += dRdx * u_disp[k * 2 + 1];
                            dvdy += dRdy * u_disp[k * 2 + 1];
                        }

                        const Exx = dudx + 0.5 * (dudx*dudx + dvdx*dvdx);
                        const Eyy = dvdy + 0.5 * (dudy*dudy + dvdy*dvdy);
                        const Exy2 = (dudy + dvdx) + (dudx*dudy + dvdx*dvdy); // 2Exy
                        
                        const D = this.getPlaneStressD();
                        const Sxx = D[0][0]*Exx + D[0][1]*Eyy;
                        const Syy = D[1][0]*Exx + D[1][1]*Eyy;
                        const Sxy = D[2][2]*Exy2;

                        const factor = detJ * wu * wv * this.thickness;

                        for (let a = 0; a < activeIndices.length; a++) {
                            const k = activeIndices[a];
                            const dRdx = B_param[k][0], dRdy = B_param[k][1];
                            const bexx_u = (1 + dudx) * dRdx;
                            const bexx_v = (dvdx) * dRdx;
                            const beyy_u = (dudy) * dRdy;
                            const beyy_v = (1 + dvdy) * dRdy;
                            const bexy_u = (1 + dudx)*dRdy + dudy*dRdx;
                            const bexy_v = (1 + dvdy)*dRdx + dvdx*dRdy;

                            const val_u = (bexx_u * Sxx + beyy_u * Syy + bexy_u * Sxy) * factor;
                            const val_v = (bexx_v * Sxx + beyy_v * Syy + bexy_v * Sxy) * factor;

                            if (!isNaN(val_u)) F_int[k * 2] += val_u;
                            if (!isNaN(val_v)) F_int[k * 2 + 1] += val_v;
                        }
                    }
                }
            }
        }
        return F_int;
    }

    calculateTangentStiffness(patch, u_disp) {
        const { p, q, U, V, controlPoints } = patch;
        const nBasisU = controlPoints.length;
        const nBasisV = controlPoints[0].length;
        const nDofs = nBasisU * nBasisV * 2;
        const Kt = Array.from({ length: nDofs }, () => new Float64Array(nDofs).fill(0));

        const uniqueU = [...new Set(U)], uniqueV = [...new Set(V)];
        const gRule = GaussQuadrature2D.getPoints(Math.max(p, q) + 1);

        for (let i = 0; i < uniqueU.length - 1; i++) {
            const uMin = uniqueU[i], uMax = uniqueU[i+1];
            if (uMax - uMin < 1e-10) continue;
            for (let j = 0; j < uniqueV.length - 1; j++) {
                const vMin = uniqueV[j], vMax = uniqueV[j+1];
                if (vMax - vMin < 1e-10) continue;

                for (let gu = 0; gu < gRule.points.length; gu++) {
                    const u = ((uMax - uMin) * gRule.points[gu] + (uMax + uMin)) / 2;
                    const wu = gRule.weights[gu] * (uMax - uMin) / 2;
                    for (let gv = 0; gv < gRule.points.length; gv++) {
                        const v = ((vMax - vMin) * gRule.points[gv] + (vMax + vMin)) / 2;
                        const wv = gRule.weights[gv] * (vMax - vMin) / 2;

                        const deriv = this.engine.getSurfaceDerivatives(patch, u, v);
                        const { grads: B_param, detJ, activeIndices } = this.getBParametric(patch, u, v, deriv);

                        let dudx = 0, dudy = 0, dvdx = 0, dvdy = 0;
                        for (let a = 0; a < activeIndices.length; a++) {
                            const k = activeIndices[a];
                            dudx += B_param[k][0] * u_disp[k * 2];
                            dudy += B_param[k][1] * u_disp[k * 2];
                            dvdx += B_param[k][0] * u_disp[k * 2 + 1];
                            dvdy += B_param[k][1] * u_disp[k * 2 + 1];
                        }

                        const D = this.getPlaneStressD();
                        const factor = detJ * wu * wv * this.thickness;
                        
                        // Assemble Linear & Material Nonlinear
                        for (let aIdx = 0; aIdx < activeIndices.length; aIdx++) {
                            const a = activeIndices[aIdx];
                            const dRdx_a = B_param[a][0], dRdy_a = B_param[a][1];
                            const B_a = [
                                [(1 + dudx) * dRdx_a, (dvdx) * dRdx_a],
                                [(dudy) * dRdy_a, (1 + dvdy) * dRdy_a],
                                [(1 + dudx) * dRdy_a + dudy * dRdx_a, (1 + dvdy) * dRdx_a + dvdx * dRdy_a]
                            ];
                            for (let bIdx = 0; bIdx < activeIndices.length; bIdx++) {
                                const b = activeIndices[bIdx];
                                const dRdx_b = B_param[b][0], dRdy_b = B_param[b][1];
                                const B_b = [
                                    [(1 + dudx) * dRdx_b, (dvdx) * dRdx_b],
                                    [(dudy) * dRdy_b, (1 + dvdy) * dRdy_b],
                                    [(1 + dudx) * dRdy_b + dudy * dRdx_b, (1 + dvdy) * dRdx_b + dvdx * dRdy_b]
                                ];

                                for (let ri = 0; ri < 2; ri++) {
                                    for (let rj = 0; rj < 2; rj++) {
                                        let kab = 0;
                                        for (let r = 0; r < 3; r++) {
                                            for (let c = 0; c < 3; c++) {
                                                kab += B_a[r][ri] * D[r][c] * B_b[c][rj];
                                            }
                                        }
                                        Kt[a * 2 + ri][b * 2 + rj] += kab * factor;
                                    }
                                }
                            }
                        }

                        // Assemble Geometric Stiffness
                        const Exx = dudx + 0.5 * (dudx*dudx + dvdx*dvdx);
                        const Eyy = dvdy + 0.5 * (dudy*dudy + dvdy*dvdy);
                        const Exy2 = (dudy + dvdx) + (dudx*dudy + dvdx*dvdy);
                        const Sxx = D[0][0]*Exx + D[0][1]*Eyy;
                        const Syy = D[1][0]*Exx + D[1][1]*Eyy;
                        const Sxy = D[2][2]*Exy2;

                        for (let aIdx = 0; aIdx < activeIndices.length; aIdx++) {
                            const a = activeIndices[aIdx];
                            for (let bIdx = 0; bIdx < activeIndices.length; bIdx++) {
                                const b = activeIndices[bIdx];
                                const dRdx_a = B_param[a][0], dRdy_a = B_param[a][1];
                                const dRdx_b = B_param[b][0], dRdy_b = B_param[b][1];
                                
                                const k_geo = (dRdx_a * Sxx * dRdx_b + dRdy_a * Syy * dRdy_b + 
                                               dRdx_a * Sxy * dRdy_b + dRdy_a * Sxy * dRdx_b) * factor;
                                
                                Kt[a * 2][b * 2] += isNaN(k_geo) ? 0 : k_geo;
                                Kt[a * 2 + 1][b * 2 + 1] += isNaN(k_geo) ? 0 : k_geo;
                            }
                        }
                    }
                }
            }
        }
        return Kt;
    }

    applyPenaltyConstraints(Kt, F_int, u, patch) {
        const { controlPoints } = patch;
        const nU = controlPoints.length;
        const nV = controlPoints[0].length;
        const pts = [];
        
        for (let i = 0; i < nU; i++) {
            for (let j = 0; j < nV; j++) {
                pts.push({ id: (i * nV + j) * 2, cp: controlPoints[i][j] });
            }
        }
        
        const penalty = 1e9; // Fixed penalty for C0 coupling
        for (let i = 0; i < pts.length; i++) {
            for (let j = i + 1; j < pts.length; j++) {
                const p1 = pts[i].cp, p2 = pts[j].cp;
                const d2 = (p1.x - p2.x)**2 + (p1.y - p2.y)**2;
                
                if (d2 < 1e-18) { // Overlapping periodic or degenerate nodes
                    const id1 = pts[i].id, id2 = pts[j].id;
                    if (Kt) {
                        Kt[id1][id1] += penalty; Kt[id2][id2] += penalty;
                        Kt[id1][id2] -= penalty; Kt[id2][id1] -= penalty;
                        Kt[id1+1][id1+1] += penalty; Kt[id2+1][id2+1] += penalty;
                        Kt[id1+1][id2+1] -= penalty; Kt[id2+1][id1+1] -= penalty;
                    }
                    if (F_int && u) {
                        const dux = u[id1] - u[id2], duy = u[id1+1] - u[id2+1];
                        F_int[id1] += penalty * dux; F_int[id2] -= penalty * dux;
                        F_int[id1+1] += penalty * duy; F_int[id2+1] -= penalty * duy;
                    }
                }
            }
        }
    }

    /**
     * Calculate Follower Forces (Traction that turns with the body)
     * For 3.3 Annulus: Applies tangential traction on the outer boundary.
     */
    calculateFollowerForces(patch, totalTorque, uDisp) {
        const { controlPoints } = patch;
        const nU = controlPoints.length;
        const nV = controlPoints[0].length;
        const F_follow = new Float64Array(nU * nV * 2).fill(0);
        
        for (let i = 0; i < nU; i++) {
            const cp = controlPoints[i][nV - 1]; // Outer boundary
            let px = cp.x + uDisp[(i * nV + nV - 1) * 2];
            let py = cp.y + uDisp[(i * nV + nV - 1) * 2 + 1];
            
            const r = Math.sqrt(px*px + py*py);
            if (r < 1e-12) continue;
            
            const tx = -py / r; // Tangential unit vector
            const ty = px / r;
            
            let nodeForce = totalTorque / (nU - 1);
            if (i === 0 || i === nU - 1) nodeForce /= 2;
            
            F_follow[(i * nV + nV - 1) * 2] = nodeForce * tx;
            F_follow[(i * nV + nV - 1) * 2 + 1] = nodeForce * ty;
        }
        return F_follow;
    }

    getNumericalStress(patch, u_disp, u, v, E, nu, isLinear = false) {
        this.E = E;
        this.nu = nu;
        
        const { controlPoints } = patch;
        const nBasisU = controlPoints.length;
        const nBasisV = controlPoints[0].length;
        
        const deriv = this.engine.getSurfaceDerivatives(patch, u, v);
        const { grads: B_param } = this.getBParametric(patch, u, v, deriv);

        let dudx = 0, dudy = 0, dvdx = 0, dvdy = 0;
        for (let k = 0; k < nBasisU * nBasisV; k++) {
            dudx += B_param[k][0] * u_disp[k * 2];
            dudy += B_param[k][1] * u_disp[k * 2];
            dvdx += B_param[k][0] * u_disp[k * 2 + 1];
            dvdy += B_param[k][1] * u_disp[k * 2 + 1];
        }

        // Strain calculation
        let Exx, Eyy, Exy2;
        if (isLinear) {
            // Infinitesimal Strain (Linear)
            Exx = dudx;
            Eyy = dvdy;
            Exy2 = (dudy + dvdx);
        } else {
            // Green-Lagrange Strain (Nonlinear)
            Exx = dudx + 0.5 * (dudx*dudx + dvdx*dvdx);
            Eyy = dvdy + 0.5 * (dudy*dudy + dvdy*dvdy);
            Exy2 = (dudy + dvdx) + (dudx*dudy + dvdx*dvdy); // 2Exy
        }

        // Stress Calculation
        const D = this.getPlaneStressD();
        const Sxx = D[0][0]*Exx + D[0][1]*Eyy;
        const Syy = D[1][0]*Exx + D[1][1]*Eyy;
        const Sxy = D[2][2]*Exy2;

        let vonMises = Math.sqrt(Sxx*Sxx - Sxx*Syy + Syy*Syy + 3*Sxy*Sxy);
        if (isNaN(vonMises)) vonMises = 0;
        
        // Equivalent (von Mises) Strain
        const Exy = Exy2 / 2;
        let equivalentStrain = Math.sqrt(Exx*Exx - Exx*Eyy + Eyy*Eyy + 3*Exy*Exy);
        if (isNaN(equivalentStrain)) equivalentStrain = 0;

        return {
            sxx: Sxx,
            syy: Syy,
            sxy: Sxy,
            vonMises: vonMises,
            vonMisesStrain: equivalentStrain
        };
    }

    solveNonlinear(patch, bcs, loads, options = {}) {
        const { iterations = 10, tolerance = 1e-6, steps = 1, isLinear = false, initialU = null, onProgress } = options;
        const nU = patch.controlPoints.length;
        const nV = patch.controlPoints[0].length;
        const nDofs = nU * nV * 2;
        
        let u = initialU ? new Float64Array(initialU) : new Float64Array(nDofs).fill(0);
        const residualHistory = [];

        // Apply specific loads
        const F_ext_total = new Float64Array(nDofs).fill(0);
        loads.forEach(load => {
            if (load.type === 'nodal') {
                const idx = (load.i * nV + load.j) * 2;
                F_ext_total[idx] += load.fx;
                F_ext_total[idx + 1] += load.fy;
            }
        });

        const fixedDofs = new Map();
        bcs.forEach(bc => {
            const baseIdx = (bc.i * nV + bc.j) * 2;
            const val = bc.value || 0;
            if (bc.axis === 'x' || bc.axis === 'both') fixedDofs.set(baseIdx, val);
            if (bc.axis === 'y' || bc.axis === 'both') fixedDofs.set(baseIdx + 1, val);
        });

        // Preallocate for solver to avoid garbage collection
        const nFree = nDofs - fixedDofs.size;
        const freeIndices = [];
        for (let i = 0; i < nDofs; i++) if (!fixedDofs.has(i)) freeIndices.push(i);

        const Kt_red = Array.from({ length: nFree }, () => new Float64Array(nFree));
        const R_red = new Float64Array(nFree);
        const Kt_lin = isLinear ? this.calculateTangentStiffness(patch, new Float64Array(nDofs)) : null;
        if (isLinear && Kt_lin) this.applyPenaltyConstraints(Kt_lin, null, u, patch);

        for (let s = 1; s <= steps; s++) {
            const fraction = s / steps;
            
            // Incrementally enforce prescribed non-zero displacements
            fixedDofs.forEach((val, idx) => { u[idx] = val * fraction; });

            for (let iter = 0; iter < iterations; iter++) {
                // Update Follower Loads (dynamic) if requested
                let F_ext_current = F_ext_total.map(f => f * fraction);
                if (options.updateLoad) {
                    const F_follow = options.updateLoad(u, fraction);
                    for (let i = 0; i < nDofs; i++) F_ext_current[i] += F_follow[i];
                }

                let F_int = isLinear ? new Float64Array(nDofs) : this.calculateInternalForce(patch, u);
                if (isLinear && Kt_lin) {
                    for (let i = 0; i < nDofs; i++) {
                        for (let j = 0; j < nDofs; j++) F_int[i] += Kt_lin[i][j] * u[j];
                    }
                }
                
                this.applyPenaltyConstraints(null, F_int, u, patch);
                const R = F_ext_current.map((f, i) => f - F_int[i]);

                let resNorm = 0;
                fixedDofs.forEach((_, idx) => R[idx] = 0);
                R.forEach(ri => resNorm += ri * ri);
                const norm = Math.sqrt(resNorm);
                
                if (isNaN(norm) || norm > 1e15) {
                    console.error("Nonlinear Solver: Divergence detected.");
                    return { u, residualHistory, diverged: true };
                }
                
                residualHistory.push({step: s, iter, norm});
                if (onProgress) onProgress({ step: s, iter, norm });
                if (norm < tolerance && iter > 0) break;
                if (isLinear && iter > 0) break;

                const Kt = isLinear ? Kt_lin : this.calculateTangentStiffness(patch, u);
                this.applyPenaltyConstraints(Kt, null, u, patch);

                for (let i = 0; i < nFree; i++) {
                    const row = freeIndices[i];
                    R_red[i] = R[row];
                    for (let j = 0; j < nFree; j++) Kt_red[i][j] = Kt[row][freeIndices[j]];
                }

                const du_red = this.gaussianElimination(Kt_red, R_red);
                for (let i = 0; i < nFree; i++) u[freeIndices[i]] += du_red[i];
            }
        }
        return { u, residualHistory };
    }

    gaussianElimination(A, b) {
        const n = b.length;
        for (let i = 0; i < n; i++) {
            let max = i;
            for (let j = i + 1; j < n; j++) if (Math.abs(A[j][i]) > Math.abs(A[max][i])) max = j;
            [A[i], A[max]] = [A[max], A[i]]; [b[i], b[max]] = [b[max], b[i]];
            if (Math.abs(A[i][i]) < 1e-20) A[i][i] = 1e-20;
            for (let j = i + 1; j < n; j++) {
                const f = A[j][i] / A[i][i];
                b[j] -= f * b[i];
                for (let k = i; k < n; k++) A[j][k] -= f * A[i][k];
            }
        }
        const x = new Float64Array(n);
        for (let i = n - 1; i >= 0; i--) {
            let s = 0;
            for (let j = i + 1; j < n; j++) s += A[i][j] * x[j];
            const div = A[i][i];
            x[i] = (Math.abs(div) < 1e-25) ? 0 : (b[i] - s) / div;
            if (isNaN(x[i])) x[i] = 0;
        }
        return x;
    }
}

// Expose globally for script-tag loading compatibility
window.IGANonlinearSolver = IGANonlinearSolver;
