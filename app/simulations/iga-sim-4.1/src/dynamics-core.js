/**
 * Dynamics Core — Nonlinear Transient Solver (Phase 4.0)
 * Implements Newmark-beta time integration for IGA elements.
 */

class DynamicsSolver {
    constructor(patch, solverFOM) {
        this.patch = patch;
        this.fom = solverFOM; // Reuse the nonlinear FOM solver for K and Fint
        this.dofs = patch.controlPoints.length * patch.controlPoints[0].length * 2;
        
        // Dynamic state vectors
        this.u = new Float64Array(this.dofs);
        this.v = new Float64Array(this.dofs);
        this.a = new Float64Array(this.dofs);
        
        // Matrices
        this.M = null; // Global Mass
        this.C = null; // Global Damping
        
        // Parameters
        this.rho = 0.001; // Density (scaled for visualization)
        this.alpha = 0.5; // Rayleigh alpha (Mass)
        this.betaR = 0.001; // Rayleigh beta (Stiffness)
        this.beta = 0.25; // Newmark beta
        this.gamma = 0.5; // Newmark gamma
    }

    /**
     * Assemble the Global Consistent Mass Matrix M
     */
    assembleMass() {
        const { p, q, U, V, controlPoints } = this.patch;
        const nU = controlPoints.length;
        const nV = controlPoints[0].length;
        this.M = new Array(this.dofs).fill(0).map(() => new Float64Array(this.dofs));

        const uniqueU = [...new Set(U)];
        const uniqueV = [...new Set(V)];

        for (let i = 0; i < uniqueU.length - 1; i++) {
            for (let j = 0; j < uniqueV.length - 1; j++) {
                if (uniqueU[i+1] <= uniqueU[i] || uniqueV[j+1] <= uniqueV[j]) continue;

                // Gauss Integration (p+1 points)
                const gPoints = this.fom.gaussPoints || 3;
                const { weights, coords } = this.getGaussRule(gPoints);

                for (let wu = 0; wu < gPoints; wu++) {
                    for (let wv = 0; wv < gPoints; wv++) {
                        const u = ((uniqueU[i+1] - uniqueU[i]) * coords[wu] + (uniqueU[i+1] + uniqueU[i])) / 2;
                        const v = ((uniqueV[j+1] - uniqueV[j]) * coords[wv] + (uniqueV[j+1] + uniqueV[j])) / 2;
                        const weight = weights[wu] * weights[wv] * (uniqueU[i+1] - uniqueU[i]) * (uniqueV[j+1] - uniqueV[j]) / 4;

                        const spanU = window.nurbsUtils.findSpan(nU - 1, p, u, U);
                        const spanV = window.nurbsUtils.findSpan(nV - 1, q, v, V);
                        const ders = window.nurbsUtils.dersBasisFuns(spanU, p, u, U, 0);
                        const dersV = window.nurbsUtils.dersBasisFuns(spanV, q, v, V, 0);

                        // Basis functions N
                        const N = [];
                        for (let ki = 0; ki <= p; ki++) {
                            for (let kj = 0; kj <= q; kj++) {
                                N.push(ders[0][ki] * dersV[0][kj]);
                            }
                        }

                        const Jac = this.fom.calculateJacobian(spanU, spanV, u, v, ders, dersV);
                        const detJ = Jac.determinant();
                        const dV = detJ * weight * this.rho;

                        // Local Mass Matrix: m_ij = rho * N_i * N_j
                        const activeIndices = this.fom.getActiveIndices(spanU, spanV, p, q, nV);
                        for (let r = 0; r < N.length; r++) {
                            for (let c = 0; c < N.length; c++) {
                                const val = N[r] * N[c] * dV;
                                // Add to x and y DOFs
                                const idxR = activeIndices[r] * 2;
                                const idxC = activeIndices[c] * 2;
                                this.M[idxR][idxC] += val;
                                this.M[idxR+1][idxC+1] += val;
                            }
                        }
                    }
                }
            }
        }
    }

    getGaussRule(n) {
        if (n === 3) return { weights: [5/9, 8/9, 5/9], coords: [-Math.sqrt(0.6), 0, Math.sqrt(0.6)] };
        return { weights: [1, 1], coords: [-1/Math.sqrt(3), 1/Math.sqrt(3)] };
    }

    /**
     * Perform one transient step using Newmark-beta
     */
    async solveStep(dt, Fext, options = {}) {
        const tol = options.tol || 1e-6;
        const maxIters = options.maxIters || 10;
        
        // Predictors
        const u_n = new Float64Array(this.u);
        const v_n = new Float64Array(this.v);
        const a_n = new Float64Array(this.a);

        // Constants for Newmark
        const a0 = 1 / (this.beta * dt * dt);
        const a1 = this.gamma / (this.beta * dt);
        const a2 = 1 / (this.beta * dt);
        const a3 = 1 / (2 * this.beta) - 1;
        const a4 = this.gamma / this.beta - 1;
        const a5 = (dt / 2) * (this.gamma / this.beta - 2);

        // Initial guess for u_{n+1}
        let u_next = new Float64Array(u_n);
        
        for (let iter = 0; iter < maxIters; iter++) {
            // Compute accelerations and velocities from current u_next
            // a_{n+1} = a0*(u_{n+1} - u_n) - a2*v_n - a3*a_n
            const a_next = new Float64Array(this.dofs);
            const v_next = new Float64Array(this.dofs);
            for (let i = 0; i < this.dofs; i++) {
                a_next[i] = a0 * (u_next[i] - u_n[i]) - a2 * v_n[i] - a3 * a_n[i];
                v_next[i] = v_n[i] + (1 - this.gamma) * dt * a_n[i] + this.gamma * dt * a_next[i];
            }

            // Internal Forces and Tangent Stiffness at u_next
            const { Fint, Kt } = this.fom.computeInternalForces(u_next);
            
            // Damping Matrix C = alpha*M + betaR*Kt
            // Residual R = M*a_next + C*v_next + Fint - Fext
            const R = new Float64Array(this.dofs);
            const Keff = new Array(this.dofs).fill(0).map(() => new Float64Array(this.dofs));

            for (let i = 0; i < this.dofs; i++) {
                // M*a_next component
                let Ma = 0;
                for (let j = 0; j < this.dofs; j++) Ma += this.M[i][j] * a_next[j];

                // C*v_next component
                let Cv = 0;
                for (let j = 0; j < this.dofs; j++) {
                    const Cij = this.alpha * this.M[i][j] + this.betaR * Kt[i][j];
                    Cv += Cij * v_next[j];
                    
                    // Effective Tangent: Keff = Kt + a0*M + a1*C
                    Keff[i][j] = Kt[i][j] + a0 * this.M[i][j] + a1 * Cij;
                }

                R[i] = Ma + Cv + Fint[i] - Fext[i];
            }

            // Apply Boundary Conditions to R and Keff (Penalty or Direct)
            this.fom.applyBCs(Keff, R);

            // Solve Keff * du = -R
            const du = this.fom.linearSolver(Keff, R.map(x => -x));

            // Update
            let norm = 0;
            for (let i = 0; i < this.dofs; i++) {
                u_next[i] += du[i];
                norm += du[i] * du[i];
            }

            if (Math.sqrt(norm) < tol) {
                // Converged! Update state
                this.u = u_next;
                this.a = a_next;
                this.v = v_next;
                return { iters: iter + 1, converged: true };
            }
        }

        return { iters: maxIters, converged: false };
    }
}
