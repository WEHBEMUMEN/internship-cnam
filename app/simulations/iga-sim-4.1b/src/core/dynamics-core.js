/**
 * Dynamics Core — Nonlinear Transient Solver (Phase 4.0)
 * Implements Newmark-beta time integration for IGA elements.
 */

class DynamicsSolver {
    constructor(patch, solverFOM) {
        this.patch = patch;
        this.fom = solverFOM; 
        this.load = patch.load;
        this.fixedDOFs = new Set();
        this.dofs = patch.controlPoints.length * patch.controlPoints[0].length * 2;
        
        // Dynamic state vectors
        this.u = new Float64Array(this.dofs);
        this.v = new Float64Array(this.dofs);
        this.a = new Float64Array(this.dofs);
        
        // Matrices
        this.M = null; // Global Mass
        this.C = null; // Global Damping
        
        // Parameters (SI Units: m, kg, N)
        this.rho = 1000.0; // Density (kg/m3)
        this.alpha = 0.5;  // Rayleigh Mass damping
        this.betaR = 0.0;  // Rayleigh Stiffness damping
        this.beta = 0.25;  // Newmark beta
        this.gamma = 0.5;  // Newmark gamma
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

                        const spanU = this.fom.engine.findSpan(nU - 1, p, u, U);
                        const spanV = this.fom.engine.findSpan(nV - 1, q, v, V);
                        const ders = this.fom.engine.basisFunsDerivs(spanU, u, p, U, 0);
                        const dersV = this.fom.engine.basisFunsDerivs(spanV, v, q, V, 0);

                        // Basis functions N
                        const N = [];
                        for (let ki = 0; ki <= p; ki++) {
                            for (let kj = 0; kj <= q; kj++) {
                                N.push(ders[0][ki] * dersV[0][kj]);
                            }
                        }

                        // Jacobian Determinant
                        const deriv = this.fom.engine.getSurfaceDerivatives(this.patch, u, v);
                        const detJ = Math.abs(deriv.dU.x * deriv.dV.y - deriv.dU.y * deriv.dV.x);
                        const dV = detJ * weight * this.rho;

                        // Local indices and Mass assembly
                        const activeIndices = [];
                        for (let ki = 0; ki <= p; ki++) {
                            for (let kj = 0; kj <= q; kj++) {
                                activeIndices.push((spanU - p + ki) * nV + (spanV - q + kj));
                            }
                        }

                        for (let r = 0; r < N.length; r++) {
                            for (let c = 0; c < N.length; c++) {
                                const val = N[r] * N[c] * dV;
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
            const Fint = this.fom.calculateInternalForce(this.patch, u_next);
            const Kt = this.fom.calculateTangentStiffness(this.patch, u_next);
            
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

            // Apply Boundary Conditions to R and Keff
            this.fom.applyPenaltyConstraints(Keff, R, u_next, this.patch, []); 

            // Solve Keff * du = -R
            const du = this.fom.gaussianElimination(Keff, R.map(x => -x));

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
