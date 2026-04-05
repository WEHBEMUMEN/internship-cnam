/**
 * 2D IGA Solver Module
 * Handles Stiffness Assembly, Numerical Integration, and Boundary Conditions.
 * Phase 2B.1 | Computational Core
 */

class GaussQuadrature2D {
    /**
     * Get Gauss points and weights for an n x n rule on [-1, 1]
     */
    static getPoints(n = 3) {
        if (n === 2) {
            const p = 1.0 / Math.sqrt(3);
            return {
                points: [-p, p],
                weights: [1, 1]
            };
        }
        if (n === 3) {
            return {
                points: [-Math.sqrt(0.6), 0, Math.sqrt(0.6)],
                weights: [5/9, 8/9, 5/9]
            };
        }
        // $n=4$
        const p1 = Math.sqrt((3 - 2 * Math.sqrt(1.2)) / 7);
        const p2 = Math.sqrt((3 + 2 * Math.sqrt(1.2)) / 7);
        const w1 = (18 + Math.sqrt(30)) / 36;
        const w2 = (18 - Math.sqrt(30)) / 36;
        return {
            points: [-p2, -p1, p1, p2],
            weights: [w2, w1, w1, w2]
        };
    }
}

class IGA2DSolver {
    constructor(engine) {
        this.engine = engine;
        this.E = 200000; // Young's Modulus (MPa) - Default: Steel
        this.nu = 0.3;   // Poisson's Ratio
        this.thickness = 1.0;
    }

    /**
     * Assemble Linear Stiffness Matrix (K)
     * K = sum( Ke )
     */
    assembleStiffness(patch) {
        const { p, q, U, V, weights, controlPoints } = patch;
        const nBasisU = controlPoints.length;
        const nBasisV = controlPoints[0].length;
        const nDofs = nBasisU * nBasisV * 2; // 2 DOF per control point (u, v)
        
        const K = Array.from({ length: nDofs }, () => new Float64Array(nDofs).fill(0));
        
        // Element strategy: knot spans
        const uniqueU = [...new Set(U)];
        const uniqueV = [...new Set(V)];
        
        const gRule = GaussQuadrature2D.getPoints(Math.max(p, q) + 1);
        
        for (let i = 0; i < uniqueU.length - 1; i++) {
            const uMin = uniqueU[i];
            const uMax = uniqueU[i+1];
            if (uMax - uMin < 1e-10) continue;

            for (let j = 0; j < uniqueV.length - 1; j++) {
                const vMin = uniqueV[j];
                const vMax = uniqueV[j+1];
                if (vMax - vMin < 1e-10) continue;

                // --- Element Integration ---
                for (let gu = 0; gu < gRule.points.length; gu++) {
                    // Map Gauss point to [uMin, uMax]
                    const u = ((uMax - uMin) * gRule.points[gu] + (uMax + uMin)) / 2;
                    const wu = gRule.weights[gu] * (uMax - uMin) / 2;

                    for (let gv = 0; gv < gRule.points.length; gv++) {
                        const v = ((vMax - vMin) * gRule.points[gv] + (vMax + vMin)) / 2;
                        const wv = gRule.weights[gv] * (vMax - vMin) / 2;

                        const detJ = this.engine.getJacobianDeterminant(patch, u, v);
                        const deriv = this.engine.getSurfaceDerivatives(patch, u, v);
                        
                        // We need Inverse Jacobian mapping for physical derivatives
                        // For a 2D surface in 3D, the "inverse" is the Moore-Penrose pseudo-inverse
                        // or we solve for the physical plane components.
                        // For Phase 2B.1, we assume a planar patch (X-Y plane) for initial verification.
                        
                        // Local Tangent basis
                        const ex = deriv.dU;
                        const ey = deriv.dV;
                        
                        // Compute Local D (Plane Stress)
                        const D = this.getPlaneStressD();
                        
                        // Compute B Matrix (Strain-Displacement)
                        const B = this.getBMatrix(patch, u, v, deriv);
                        
                        // Contribution: B' * D * B * detJ * weight
                        const factor = detJ * wu * wv * this.thickness;
                        
                        this.accumulateContribution(K, B, D, factor, nBasisU, nBasisV);
                    }
                }
            }
        }
        return K;
    }

    getPlaneStressD() {
        const factor = this.E / (1 - this.nu * this.nu);
        return [
            [factor, factor * this.nu, 0],
            [factor * this.nu, factor, 0],
            [0, 0, factor * (1 - this.nu) / 2]
        ];
    }

    /**
     * B Matrix for 2D IGA
     * Maps control point [ux, uy] to [eps_x, eps_y, gamma_xy]
     */
    getBMatrix(patch, u, v, deriv) {
        const { p, q, U, V, weights, controlPoints } = patch;
        const nU = controlPoints.length;
        const nV = controlPoints[0].length;
        
        // 1. Calculate physical derivatives dRi/dx, dRi/dy
        // Using J_inv * [dRi/du, dRi/dv]'
        const J = [
            [deriv.dU.x, deriv.dV.x],
            [deriv.dU.y, deriv.dV.y]
        ];
        const detJ_2D = J[0][0] * J[1][1] - J[0][1] * J[1][0];
        const J_inv = [
            [ J[1][1]/detJ_2D, -J[0][1]/detJ_2D],
            [-J[1][0]/detJ_2D,  J[0][0]/detJ_2D]
        ];

        const B = []; // Array of matrices for each CP [B1, B2, ...]
        
        // Denominator derivatives
        const W = deriv.W;
        const Wu = deriv.Wu;
        const Wv = deriv.Wv;

        for (let i = 0; i < nU; i++) {
            const Ni = this.engine.basis1D(i, p, U, u);
            const dNi = this.engine.basis1DDeriv(i, p, U, u);
            
            for (let j = 0; j < nV; j++) {
                const Mj = this.engine.basis1D(j, q, V, v);
                const dMj = this.engine.basis1DDeriv(j, q, V, v);
                
                const w = weights[i][j];
                
                // Basis derivatives w.r.t. parametric u, v
                const dRdu = ( (dNi * Mj * w) * W - (Ni * Mj * w) * Wu ) / (W * W);
                const dRdv = ( (Ni * dMj * w) * W - (Ni * Mj * w) * Wv ) / (W * W);
                
                // Basis derivatives w.r.t. physical x, y
                const dRdx = J_inv[0][0] * dRdu + J_inv[0][1] * dRdv;
                const dRdy = J_inv[1][0] * dRdu + J_inv[1][1] * dRdv;

                // [B_ij] = [ dRdx  0   ]
                //          [ 0     dRdy]
                //          [ dRdy  dRdx]
                B.push([
                    [dRdx, 0],
                    [0, dRdy],
                    [dRdy, dRdx]
                ]);
            }
        }
        return B;
    }

    accumulateContribution(K, B, D, factor, nU, nV) {
        const nPoints = nU * nV;
        for (let a = 0; a < nPoints; a++) {
            for (let b = 0; b < nPoints; b++) {
                // Ke_ab = Ba' * D * Bb
                const Ba = B[a];
                const Bb = B[b];

                for (let i = 0; i < 2; i++) { // DOF a
                    for (let j = 0; j < 2; j++) { // DOF b
                        
                        let kab = 0;
                        // Manual matrix multiplication: row_i(Ba') * D * col_j(Bb)
                        // Equivalent to: Column_i(Ba) * D * Column_j(Bb)
                        
                        // D is 3x3, Ba is 3x2, Bb is 3x2
                        for (let r = 0; r < 3; r++) {
                            for (let c = 0; c < 3; c++) {
                                kab += Ba[r][i] * D[r][c] * Bb[c][j];
                            }
                        }
                        
                        K[a * 2 + i][b * 2 + j] += kab * factor;
                    }
                }
            }
        }
    }

    /**
     * Solve Linear System: K * d = f
     * Applies Dirichlet BCs (constrained DOFs) and Neumann BCs (load vector)
     */
    solve(patch, bcs, loads) {
        const nU = patch.controlPoints.length;
        const nV = patch.controlPoints[0].length;
        const nDofs = nU * nV * 2;
        
        const K_full = this.assembleStiffness(patch);
        const F = new Float64Array(nDofs).fill(0);
        
        // 1. Apply Point Loads (Neumann)
        loads.forEach(load => {
            // load: { i, j, fx, fy }
            const idx = (load.i * nV + load.j) * 2;
            F[idx] += load.fx;
            F[idx + 1] += load.fy;
        });

        // 2. Apply Boundary Conditions (Dirichlet)
        // Fixed DoFs approach
        const activeDofs = [];
        const fixedValues = new Map();
        
        bcs.forEach(bc => {
            // bc: { i, j, axis: 'x'|'y'|'both', value }
            const baseIdx = (bc.i * nV + bc.j) * 2;
            if (bc.axis === 'x' || bc.axis === 'both') fixedValues.set(baseIdx, bc.value);
            if (bc.axis === 'y' || bc.axis === 'both') fixedValues.set(baseIdx + 1, bc.value);
        });

        // Reduction
        const nFree = nDofs - fixedValues.size;
        const freeIndices = [];
        for (let i = 0; i < nDofs; i++) if (!fixedValues.has(i)) freeIndices.push(i);
        
        const K_red = Array.from({ length: nFree }, () => new Float64Array(nFree).fill(0));
        const F_red = new Float64Array(nFree);

        for (let i = 0; i < nFree; i++) {
            const row = freeIndices[i];
            let fi = F[row];
            // Adjust force for non-zero Dirichlet
            fixedValues.forEach((val, col) => {
                fi -= K_full[row][col] * val;
            });
            F_red[i] = fi;

            for (let j = 0; j < nFree; j++) {
                const col = freeIndices[j];
                K_red[i][j] = K_full[row][col];
            }
        }

        // 3. Dense Solver (Gaussian Elimination)
        const sol_red = this.gaussianElimination(K_red, F_red);
        
        const displacement = new Float64Array(nDofs);
        let freeIdx = 0;
        for (let i = 0; i < nDofs; i++) {
            if (fixedValues.has(i)) {
                displacement[i] = fixedValues.get(i);
            } else {
                displacement[i] = sol_red[freeIdx++];
            }
        }

        return displacement;
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
            x[i] = (b[i] - s) / A[i][i];
        }
        return x;
    }
}

// Export for browser
if (typeof module !== 'undefined' && module.exports) {
    module.exports = { IGA2DSolver, GaussQuadrature2D };
} else {
    window.IGA2DSolver = IGA2DSolver;
    window.GaussQuadrature2D = GaussQuadrature2D;
}
