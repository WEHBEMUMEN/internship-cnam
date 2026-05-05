/**
 * Dynamics Core — Nonlinear Transient Solver (Phase 4.0)
 * Implements Newmark-beta time integration for IGA elements.
 */

class DynamicsSolver {
    constructor(patch, solverFOM) {
        this.patch = patch; // Geometry data
        this.fom = solverFOM; // IGANonlinearSolver instance
        this.engine = solverFOM.engine; // NURBS2D engine
        
        const nU = patch.controlPoints.length;
        const nV = patch.controlPoints[0].length;
        this.dofs = nU * nV * 2;
        
        // Dynamic state vectors
        this.u = new Float64Array(this.dofs);
        this.v = new Float64Array(this.dofs);
        this.a = new Float64Array(this.dofs);
        
        // Matrices
        this.M = null; // Global Mass
        
        // Parameters
        this.rho = 0.0001; // Density (kg/mm2)
        this.alpha = 0.1;  // Rayleigh Mass damping
        this.betaR = 0.001; // Rayleigh Stiffness damping
        this.beta = 0.25;  // Newmark beta (Average Acceleration)
        this.gamma = 0.5;  // Newmark gamma
    }

    /**
     * Assemble the Global Consistent Mass Matrix M
     */
    assembleMass() {
        const { p, q, U, V, weights, controlPoints } = this.patch;
        const nBasisU = controlPoints.length;
        const nBasisV = controlPoints[0].length;
        const nDofs = nBasisU * nBasisV * 2;
        
        this.M = Array.from({ length: nDofs }, () => new Float64Array(nDofs).fill(0));

        const uniqueU = [...new Set(U)], uniqueV = [...new Set(V)];
        const gRule = window.GaussQuadrature2D ? window.GaussQuadrature2D.getPoints(Math.max(p, q) + 1) : { points: [-0.7745, 0, 0.7745], weights: [0.555, 0.888, 0.555] };

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

                        const detJ = this.engine.getJacobianDeterminant(this.patch, u, v);
                        const factor = detJ * wu * wv * this.rho * (this.fom.thickness || 1.0);

                        const spanU = this.engine.findSpan(nBasisU - 1, p, u, U);
                        const spanV = this.engine.findSpan(nBasisV - 1, q, v, V);
                        const NU = this.engine.basisFuns(spanU, u, p, U);
                        const MV = this.engine.basisFuns(spanV, v, q, V);
                        const W = this.engine.computeDenominator(this.patch, u, v);

                        const activeIndices = [];
                        const R = [];
                        for (let ki = 0; ki <= p; ki++) {
                            for (let kj = 0; kj <= q; kj++) {
                                const idx = (spanU - p + ki) * nBasisV + (spanV - q + kj);
                                activeIndices.push(idx);
                                const w = weights[spanU - p + ki][spanV - q + kj];
                                R.push(NU[ki] * MV[kj] * w / W);
                            }
                        }

                        for (let a = 0; a < activeIndices.length; a++) {
                            for (let b = 0; b < activeIndices.length; b++) {
                                const val = R[a] * R[b] * factor;
                                const ia = activeIndices[a], ib = activeIndices[b];
                                this.M[ia * 2][ib * 2] += val;
                                this.M[ia * 2 + 1][ib * 2 + 1] += val;
                            }
                        }
                    }
                }
            }
        }
    }

    /**
     * Solve one Newmark-beta time step
     */
    async solveStep(dt, Fext, options = {}) {
        const tol = options.tol || 1e-6;
        const maxIters = options.maxIters || 15;
        
        const u_n = new Float64Array(this.u);
        const v_n = new Float64Array(this.v);
        const a_n = new Float64Array(this.a);

        // Prediction (Newmark constants)
        const a0 = 1 / (this.beta * dt * dt);
        const a1 = this.gamma / (this.beta * dt);
        const a2 = 1 / (this.beta * dt);
        const a3 = 1 / (2 * this.beta) - 1;
        const a4 = this.gamma / this.beta - 1;
        const a5 = (dt / 2) * (this.gamma / this.beta - 2);

        let u_next = new Float64Array(u_n);
        
        for (let iter = 0; iter < maxIters; iter++) {
            const a_next = new Float64Array(this.dofs);
            const v_next = new Float64Array(this.dofs);
            for (let i = 0; i < this.dofs; i++) {
                a_next[i] = a0 * (u_next[i] - u_n[i]) - a2 * v_n[i] - a3 * a_n[i];
                v_next[i] = v_n[i] + (1 - this.gamma) * dt * a_n[i] + this.gamma * dt * a_next[i];
            }

            const Fint = this.fom.calculateInternalForce(this.patch, u_next);
            const Kt = this.fom.calculateTangentStiffness(this.patch, u_next);
            
            // Effective Stiffness Keff = Kt + a0*M + a1*C
            // where C = alpha*M + betaR*Kt
            const R_res = new Float64Array(this.dofs);
            const Keff = Array.from({ length: this.dofs }, () => new Float64Array(this.dofs).fill(0));

            for (let i = 0; i < this.dofs; i++) {
                let Ma = 0, Cv = 0;
                for (let j = 0; j < this.dofs; j++) {
                    const mij = this.M[i][j];
                    const kij = Kt[i][j];
                    const cij = this.alpha * mij + this.betaR * kij;
                    
                    Ma += mij * a_next[j];
                    Cv += cij * v_next[j];
                    Keff[i][j] = kij + a0 * mij + a1 * cij;
                }
                R_res[i] = Ma + Cv + Fint[i] - Fext[i];
            }

            // Boundary Conditions (Fixed Left End)
            const nV = this.patch.controlPoints[0].length;
            const bcs = [];
            for (let j = 0; j < nV; j++) bcs.push({ i: 0, j: j, axis: 'both', value: 0 });
            
            this.fom.applyPenaltyConstraints(Keff, R_res, u_next, this.patch, bcs);

            // Solve step: Keff * du = -R_res
            const du = this.fom.gaussianElimination(Keff, R_res.map(x => -x));
            
            let norm = 0;
            for (let i = 0; i < this.dofs; i++) {
                u_next[i] += du[i];
                norm += du[i] * du[i];
            }

            if (Math.sqrt(norm) < tol) {
                this.u = u_next;
                this.v = v_next;
                this.a = a_next;
                return { iters: iter + 1, converged: true };
            }
        }

        this.u = u_next;
        return { iters: maxIters, converged: false };
    }
}

window.DynamicsSolver = DynamicsSolver;
