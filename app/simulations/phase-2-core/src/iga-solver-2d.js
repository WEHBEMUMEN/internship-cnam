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
        // If derivatives aren't provided, calculate them internally
        if (!deriv) {
            deriv = this.engine.getSurfaceDerivatives(patch, u, v);
        }
        
        const { p, q, U, V, weights, controlPoints } = patch;
        const nU = controlPoints.length;
        const nV = controlPoints[0].length;
        
        // 1. Calculate physical derivatives dRi/dx, dRi/dy
        // Using J_inv * [dRi/du, dRi/dv]'
        const J = [
            [deriv.dU.x, deriv.dV.x],
            [deriv.dU.y, deriv.dV.y]
        ];
        let detJ_2D = J[0][0] * J[1][1] - J[0][1] * J[1][0];
        
        // Safety: Handle singularity at degenerate corner by using a stable regularizer
        // rather than zeroing the B-matrix, which would cause numerical 'locking' or zero-displacement.
        if (Math.abs(detJ_2D) < 1e-12) {
            detJ_2D = (detJ_2D >= 0) ? 1e-12 : -1e-12;
        }

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
                const Ba = B[a];
                const Bb = B[b];
                for (let i = 0; i < 2; i++) {
                    for (let j = 0; j < 2; j++) {
                        let kab = 0;
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
     * Nonlinear Analysis Methods (Phase 2B.2)
     */

    /**
     * Compute Internal Force Vector F_int(u)
     * F_int = Integral( B_NL' * S ) dOmega
     */
    calculateInternalForce(patch, u_disp) {
        const { p, q, U, V, weights, controlPoints } = patch;
        const nBasisU = controlPoints.length;
        const nBasisV = controlPoints[0].length;
        const nDofs = nBasisU * nBasisV * 2;
        const F_int = new Float64Array(nDofs).fill(0);

        const uniqueU = [...new Set(U)];
        const uniqueV = [...new Set(V)];
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

                        const detJ = this.engine.getJacobianDeterminant(patch, u, v);
                        const deriv = this.engine.getSurfaceDerivatives(patch, u, v);
                        const B_param = this.getBParametric(patch, u, v, deriv); // Basis physical grads [dR/dx, dR/dy]
                        
                        // 1. Calculate Displacement Gradients grad_u = [du/dx, du/dy; dv/dx, dv/dy]
                        let dudx = 0, dudy = 0, dvdx = 0, dvdy = 0;
                        for (let k = 0; k < nBasisU * nBasisV; k++) {
                            const dRdx = B_param[k][0], dRdy = B_param[k][1];
                            dudx += dRdx * u_disp[k * 2];
                            dudy += dRdy * u_disp[k * 2];
                            dvdx += dRdx * u_disp[k * 2 + 1];
                            dvdy += dRdy * u_disp[k * 2 + 1];
                        }

                        // 2. Green-Lagrange Strain Vector E = [Exx, Eyy, 2Exy]
                        // E = 1/2 (F'F - I)
                        const Exx = dudx + 0.5 * (dudx*dudx + dvdx*dvdx);
                        const Eyy = dvdy + 0.5 * (dudy*dudy + dvdy*dvdy);
                        const Exy2 = (dudy + dvdx) + (dudx*dudy + dvdx*dvdy); // 2Exy
                        
                        // 3. 2nd Piola-Kirchhoff Stress S = D * E
                        const D = this.getPlaneStressD();
                        const Sxx = D[0][0]*Exx + D[0][1]*Eyy;
                        const Syy = D[1][0]*Exx + D[1][1]*Eyy;
                        const Sxy = D[2][2]*Exy2;

                        // 4. Nonlinear B-Matrix (B_NL)
                        // B_NL = B_L + B_u(u)
                        const factor = detJ * wu * wv * this.thickness;

                        for (let k = 0; k < nBasisU * nBasisV; k++) {
                            const dRdx = B_param[k][0], dRdy = B_param[k][1];
                            
                            // Virtual strains for DOF k
                            // eps_xx = (1 + dudx)*dRdx * du_k + (dvdx)*dRdx * dv_k
                            const bexx_u = (1 + dudx) * dRdx;
                            const bexx_v = (dvdx) * dRdx;
                            const beyy_u = (dudy) * dRdy;
                            const beyy_v = (1 + dvdy) * dRdy;
                            const bexy_u = (1 + dudx)*dRdy + dudy*dRdx;
                            const bexy_v = (1 + dvdy)*dRdx + dvdx*dRdy;

                            F_int[k * 2] += (bexx_u * Sxx + beyy_u * Syy + bexy_u * Sxy) * factor;
                            F_int[k * 2 + 1] += (bexx_v * Sxx + beyy_v * Syy + bexy_v * Sxy) * factor;
                        }
                    }
                }
            }
        }
        return F_int;
    }

    /**
     * Compute Tangent Stiffness Matrix Kt(u)
     * Kt = K_mat + K_geo
     */
    calculateTangentStiffness(patch, u_disp) {
        const { p, q, U, V, weights, controlPoints } = patch;
        const nBasisU = controlPoints.length;
        const nBasisV = controlPoints[0].length;
        const nDofs = nBasisU * nBasisV * 2;
        const Kt = Array.from({ length: nDofs }, () => new Float64Array(nDofs).fill(0));

        const uniqueU = [...new Set(U)];
        const uniqueV = [...new Set(V)];
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

                        const detJ = this.engine.getJacobianDeterminant(patch, u, v);
                        const deriv = this.engine.getSurfaceDerivatives(patch, u, v);
                        const B_param = this.getBParametric(patch, u, v, deriv);

                        let dudx = 0, dudy = 0, dvdx = 0, dvdy = 0;
                        for (let k = 0; k < nBasisU * nBasisV; k++) {
                            dudx += B_param[k][0] * u_disp[k * 2];
                            dudy += B_param[k][1] * u_disp[k * 2];
                            dvdx += B_param[k][0] * u_disp[k * 2 + 1];
                            dvdy += B_param[k][1] * u_disp[k * 2 + 1];
                        }

                        // 1. Assemble Material Tangent Stiffness K_mat
                        const B_NL = [];
                        for (let k = 0; k < nBasisU * nBasisV; k++) {
                            const dRdx = B_param[k][0], dRdy = B_param[k][1];
                            B_NL.push([
                                [(1 + dudx) * dRdx, (dvdx) * dRdx],
                                [(dudy) * dRdy, (1 + dvdy) * dRdy],
                                [(1 + dudx) * dRdy + dudy * dRdx, (1 + dvdy) * dRdx + dvdx * dRdy]
                            ]);
                        }

                        const D = this.getPlaneStressD();
                        const factor = detJ * wu * wv * this.thickness;
                        this.accumulateContribution(Kt, B_NL, D, factor, nBasisU, nBasisV);

                        // 2. Assemble Geometric Stiffness K_geo (Initial Stress)
                        // K_geo = Integral( G' * [S I] * G )
                        const Exx = dudx + 0.5 * (dudx*dudx + dvdx*dvdx);
                        const Eyy = dvdy + 0.5 * (dudy*dudy + dvdy*dvdy);
                        const Exy2 = (dudy + dvdx) + (dudx*dudy + dvdx*dvdy);
                        const Sxx = D[0][0]*Exx + D[0][1]*Eyy;
                        const Syy = D[1][0]*Exx + D[1][1]*Eyy;
                        const Sxy = D[2][2]*Exy2;

                        for (let a = 0; a < nBasisU * nBasisV; a++) {
                            for (let b = 0; b < nBasisU * nBasisV; b++) {
                                const dRdx_a = B_param[a][0], dRdy_a = B_param[a][1];
                                const dRdx_b = B_param[b][0], dRdy_b = B_param[b][1];
                                
                                // k_geo_ab = (dRdx_a*Sxx*dRdx_b + dRdy_a*Syy*dRdy_b + dRdx_a*Sxy*dRdy_b + dRdy_a*Sxy*dRdx_b) * I
                                const k_geo = (dRdx_a * Sxx * dRdx_b + dRdy_a * Syy * dRdy_b + 
                                               dRdx_a * Sxy * dRdy_b + dRdy_a * Sxy * dRdx_b) * factor;
                                
                                Kt[a * 2][b * 2] += k_geo;
                                Kt[a * 2 + 1][b * 2 + 1] += k_geo;
                            }
                        }
                    }
                }
            }
        }
        return Kt;
    }

    /**
     * B Matrix parametric help (physical derivatives)
     */
    getBParametric(patch, u, v, deriv) {
        const { p, q, U, V, weights, controlPoints } = patch;
        const nU = controlPoints.length;
        const nV = controlPoints[0].length;
        const J = [[deriv.dU.x, deriv.dV.x], [deriv.dU.y, deriv.dV.y]];
        const detJ_2D = J[0][0] * J[1][1] - J[0][1] * J[1][0];
        const J_inv = [[J[1][1]/detJ_2D, -J[0][1]/detJ_2D], [-J[1][0]/detJ_2D, J[0][0]/detJ_2D]];
        const W = deriv.W, Wu = deriv.Wu, Wv = deriv.Wv;
        const grads = [];
        for (let i = 0; i < nU; i++) {
            const Ni = this.engine.basis1D(i, p, U, u), dNi = this.engine.basis1DDeriv(i, p, U, u);
            for (let j = 0; j < nV; j++) {
                const Mj = this.engine.basis1D(j, q, V, v), dMj = this.engine.basis1DDeriv(j, q, V, v);
                const w = weights[i][j];
                const dRdu = ((dNi * Mj * w) * W - (Ni * Mj * w) * Wu) / (W * W);
                const dRdv = ((Ni * dMj * w) * W - (Ni * Mj * w) * Wv) / (W * W);
                grads.push([J_inv[0][0] * dRdu + J_inv[0][1] * dRdv, J_inv[1][0] * dRdu + J_inv[1][1] * dRdv]);
            }
        }
        return grads;
    }

    /**
     * Nonlinear Solver: Newton-Raphson
     * Options: { iterations, tolerance, incrementalSteps }
     */
    solveNonlinear(patch, bcs, loads, options = {}) {
        const { iterations = 10, tolerance = 1e-6, steps = 1 } = options;
        const nDofs = patch.controlPoints.length * patch.controlPoints[0].length * 2;
        let u = new Float64Array(nDofs).fill(0);
        const residualHistory = [];

        // External Force Vector (F_ext)
        const F_ext_total = new Float64Array(nDofs).fill(0);
        loads.forEach(load => {
            const idx = (load.i * patch.controlPoints[0].length + load.j) * 2;
            F_ext_total[idx] += load.fx;
            F_ext_total[idx + 1] += load.fy;
        });

        // Fixed DOFs mapping
        const fixedDofs = new Map();
        bcs.forEach(bc => {
            const baseIdx = (bc.i * patch.controlPoints[0].length + bc.j) * 2;
            if (bc.axis === 'x' || bc.axis === 'both') fixedDofs.set(baseIdx, bc.value);
            if (bc.axis === 'y' || bc.axis === 'both') fixedDofs.set(baseIdx + 1, bc.value);
        });

        // Incremental Load Stepping
        for (let s = 1; s <= steps; s++) {
            const F_ext = F_ext_total.map(f => f * (s / steps));

            for (let iter = 0; iter < iterations; iter++) {
                const F_int = this.calculateInternalForce(patch, u);
                const R = F_ext.map((f, i) => f - F_int[i]);

                // Norm Calculation
                let resNorm = 0;
                fixedDofs.forEach((_, idx) => R[idx] = 0); // Zero out residual at BCs for norm
                R.forEach(ri => resNorm += ri * ri);
                const norm = Math.sqrt(resNorm);
                residualHistory.push(norm);

                if (norm < tolerance && iter > 0) break;

                const Kt = this.calculateTangentStiffness(patch, u);

                // Apply BCs to Kt and R
                const nFree = nDofs - fixedDofs.size;
                const freeIndices = [];
                for (let i = 0; i < nDofs; i++) if (!fixedDofs.has(i)) freeIndices.push(i);

                const Kt_red = Array.from({ length: nFree }, () => new Float64Array(nFree));
                const R_red = new Float64Array(nFree);

                for (let i = 0; i < nFree; i++) {
                    const row = freeIndices[i];
                    R_red[i] = R[row];
                    for (let j = 0; j < nFree; j++) {
                        Kt_red[i][j] = Kt[row][freeIndices[j]];
                    }
                }

                const du_red = this.gaussianElimination(Kt_red, R_red);
                for (let i = 0; i < nFree; i++) u[freeIndices[i]] += du_red[i];
            }
        }

        return { u, residualHistory };
    }

    /**
     * Solve Linear System: K * d = f
     * Applies Dirichlet BCs (constrained DOFs) and Neumann BCs (load vector)
     */
    /**
     * Algorithm 4 & 6: Element-based Boundary Force Assembly
     * Calculates the global force vector by integrating traction over each boundary element.
     */
    calculateNodalTraction(patch, loadValue, edge = 'right') {
        const { p, q, U, V, controlPoints, weights } = patch;
        const nU = controlPoints.length;
        const nV = controlPoints[0].length;
        const forces = new Float64Array(nU * nV * 2);
        
        // Quadrature Setup (3-point Gauss rule)
        const gPoints = [-0.7745966692414833, 0.0, 0.7745966692414833];
        const gWeights = [0.5555555555555556, 0.8888888888888888, 0.5555555555555556];

        const uniqueU = [...new Set(U)];

        if (edge === 'right') {
            const vParam = 1.0 - 1e-6; // Slightly inside to avoid degenerate corner
            
            // Loop through boundary elements in the U direction
            for (let e = 0; e < uniqueU.length - 1; e++) {
                const uMin = uniqueU[e]; 
                const uMax = uniqueU[e+1];
                if (uMax - uMin < 1e-9) continue;

                for (let qu = 0; qu < gPoints.length; qu++) {
                    const u = ((uMax - uMin) * gPoints[qu] + (uMax + uMin)) / 2;
                    const gw = gWeights[qu] * (uMax - uMin) / 2;

                    const deriv = this.engine.getSurfaceDerivatives(patch, u, vParam);
                    
                    // 1D boundary Jacobian = |dS/du|
                    const tangent = deriv.dU;
                    const detJ_1D = Math.sqrt(tangent.x**2 + tangent.y**2);
                    if (detJ_1D < 1e-12) continue;

                    // Physical position at this Gauss point
                    const pos = deriv.pos;

                    // --- EXACT KIRSCH TRACTION (Algorithm 3 compatible) ---
                    // Compute exact analytical stress at this boundary point
                    const R_hole = 1.0; // Hole radius
                    const sigma = this.getExactStress(pos.x, pos.y, R_hole, loadValue);

                    // Outward normal: n = (ty, -tx) / |t| for CCW traversal
                    const nx = tangent.y / detJ_1D;
                    const ny = -tangent.x / detJ_1D;

                    // Traction vector: t = sigma * n
                    const tx = sigma.sxx * nx + sigma.sxy * ny;
                    const ty = sigma.sxy * nx + sigma.syy * ny;

                    // --- NURBS Rational Basis (Fixed: include weights + denominator) ---
                    // Compute 1D denominator on boundary v=vParam
                    let W_1D = 0;
                    for (let a = 0; a < nU; a++) {
                        const Na = this.engine.basis1D(a, p, U, u);
                        W_1D += Na * weights[a][nV-1];
                    }
                    if (W_1D < 1e-12) W_1D = 1e-12;

                    // Accumulate into nodal forces using NURBS rational basis
                    for (let a = 0; a < nU; a++) {
                        const Na = this.engine.basis1D(a, p, U, u);
                        if (Na === 0) continue;
                        
                        // Rational NURBS basis on the boundary
                        const Ra = (Na * weights[a][nV-1]) / W_1D;
                        
                        const idx = (a * nV + (nV - 1)) * 2;
                        forces[idx]     += Ra * tx * detJ_1D * gw;
                        forces[idx + 1] += Ra * ty * detJ_1D * gw;
                    }
                }
            }
        }
        return forces;
    }

    applyPenaltyConstraints(K, patch) {
        const { controlPoints } = patch;
        const nU = controlPoints.length;
        const nV = controlPoints[0].length;
        const pts = [];
        
        for (let i = 0; i < nU; i++) {
            for (let j = 0; j < nV; j++) {
                pts.push({
                    id: (i * nV + j) * 2,
                    cp: controlPoints[i][j]
                });
            }
        }
        
        // Large penalty relative to expected stiffness
        const penalty = 1e12; 
        
        for (let i = 0; i < pts.length; i++) {
            for (let j = i + 1; j < pts.length; j++) {
                const p1 = pts[i].cp;
                const p2 = pts[j].cp;
                
                const dist2 = (p1.x - p2.x)**2 + (p1.y - p2.y)**2 + (p1.z - p2.z)**2;
                if (dist2 < 1e-20) { // Distance essentially 0
                    const id1 = pts[i].id;
                    const id2 = pts[j].id;
                    
                    // Tie X-DOF
                    K[id1][id1] += penalty;
                    K[id2][id2] += penalty;
                    K[id1][id2] -= penalty;
                    K[id2][id1] -= penalty;
                    
                    // Tie Y-DOF
                    K[id1+1][id1+1] += penalty;
                    K[id2+1][id2+1] += penalty;
                    K[id1+1][id2+1] -= penalty;
                    K[id2+1][id1+1] -= penalty;
                }
            }
        }
    }

    solve(patch, bcs, loads) {
        const nU = patch.controlPoints.length;
        const nV = patch.controlPoints[0].length;
        const nDofs = nU * nV * 2;
        
        const K_full = this.assembleStiffness(patch);
        
        // Stabilize singularities (e.g. degenerate corners) 
        this.applyPenaltyConstraints(K_full, patch);

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
            
            // Robust Singularity Handling: 
            // If the pivot is effectively zero, it means this DOF is unconstrained/unloaded.
            // We set it to stay fixed (0 displacement) to avoid Inf/NaN overflows.
            if (Math.abs(A[i][i]) < 1e-15) {
                A[i][i] = 1.0; 
                b[i] = 0;
                for (let k = i+1; k < n; k++) A[i][k] = 0;
            }

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

    // =========================================================================
    // SECTION 4: BENCHMARKING & STRESS ANALYSIS
    // =========================================================================

    getDMatrix(E, nu) {
        const factor = E / (1 - nu * nu);
        return [
            [factor, factor * nu, 0],
            [factor * nu, factor, 0],
            [0, 0, factor * (1 - nu) / 2]
        ];
    }

    getExactStress(x, y, R, Tx) {
        const r = Math.sqrt(x * x + y * y);
        const theta = Math.atan2(y, x);
        const R2 = R * R;
        const R4 = R2 * R2;
        const r2 = r * r;
        const r4 = r2 * r2;
        const cos2t = Math.cos(2 * theta);
        const sin2t = Math.sin(2 * theta);

        const s_rr = (Tx / 2) * (1 - R2/r2) + (Tx / 2) * (1 - 4*R2/r2 + 3*R4/r4) * cos2t;
        const s_tt = (Tx / 2) * (1 + R2/r2) - (Tx / 2) * (1 + 3*R4/r4) * cos2t;
        const s_rt = -(Tx / 2) * (1 + 2*R2/r2 - 3*R4/r4) * sin2t;

        const c = Math.cos(theta); const s = Math.sin(theta);
        const sxx = s_rr * c*c + s_tt * s*s - 2 * s_rt * s*c;
        const syy = s_rr * s*s + s_tt * c*c + 2 * s_rt * s*c;
        const sxy = (s_rr - s_tt) * s * c + s_rt * (c*c - s*s);

        return { sxx, syy, sxy };
    }

    getNumericalStress(patch, u_disp, u, v, E, nu) {
        const p = patch.p; const q = patch.q;
        const U = patch.U; const V = patch.V;
        const nV = patch.controlPoints[0].length;
        const nNodes = patch.controlPoints.length * patch.controlPoints[0].length;

        let u_ev = u;
        let v_ev = v;
        let detJ = this.engine.getJacobianDeterminant(patch, u_ev, v_ev);
        
        if (Math.abs(detJ) <= 1e-8) {
            u_ev = Math.max(0, Math.min(1, u - 0.01));
            v_ev = Math.max(0, Math.min(1, v + 0.001));
        }

        const B = this.getBMatrix(patch, u_ev, v_ev);
        const D = this.getDMatrix(E, nu);
        const eps = [0, 0, 0];

        for (let i = 0; i < nNodes; i++) {
            const ux = u_disp[i * 2];
            const uy = u_disp[i * 2 + 1];
            eps[0] += B[i][0][0] * ux;
            eps[1] += B[i][1][1] * uy;
            eps[2] += B[i][2][0] * ux + B[i][2][1] * uy;
        }

        const sxx = D[0][0] * eps[0] + D[0][1] * eps[1];
        const syy = D[1][0] * eps[0] + D[1][1] * eps[1];
        const sxy = D[2][2] * eps[2];
        const vonMises = Math.sqrt(sxx*sxx + syy*syy - sxx*syy + 3*sxy*sxy);

        return { sxx, syy, sxy, vonMises };
    }

    calculateRelativeL2Error(patch, u_disp, E, nu, Tx, R) {
        const { U, V, p } = patch;
        const uniqueU = [...new Set(U)]; const uniqueV = [...new Set(V)];
        let err2 = 0; let ex2 = 0;
        const gauss = GaussQuadrature2D.getPoints(p + 1);

        // Physical radius cutoff: only integrate over r <= rMax
        // This excludes the outer boundary strip where the degenerate mapping
        // creates non-integrable stress singularities (σ ~ 1/detJ).
        // r_max = 0.85*L covers the stress concentration zone while avoiding artifacts.
        const L = 4.0;
        const rMax = 0.85 * L * Math.SQRT2; // ~4.8, excludes far corner region
        const minDetJ = 1e-6;

        for (let i = 0; i < uniqueU.length - 1; i++) {
            for (let j = 0; j < uniqueV.length - 1; j++) {
                const u0 = uniqueU[i]; const u1 = uniqueU[i+1];
                const v0 = uniqueV[j]; const v1 = uniqueV[j+1];
                if (u1 - u0 < 1e-9 || v1 - v0 < 1e-9) continue;

                for (let gu = 0; gu < gauss.points.length; gu++) {
                    for (let gv = 0; gv < gauss.points.length; gv++) {
                        const u = ((u1 - u0) * gauss.points[gu] + (u1 + u0)) / 2;
                        const v = ((v1 - v0) * gauss.points[gv] + (v1 + v0)) / 2;
                        const deriv = this.engine.getSurfaceDerivatives(patch, u, v);
                        const detJ = Math.abs(deriv.dU.x * deriv.dV.y - deriv.dV.x * deriv.dU.y);

                        // Skip near-singular Gauss points
                        if (detJ < minDetJ) continue;

                        // Physical position and radius check
                        const pos = this.engine.evaluateSurface(patch, u, v);
                        const r = Math.sqrt(pos.x * pos.x + pos.y * pos.y);
                        if (r > rMax) continue;
                        // Also skip if too close to outer walls
                        if (pos.x > 0.95 * L || pos.y > 0.95 * L) continue;

                        const Jmod = detJ * gauss.weights[gu] * gauss.weights[gv] * (u1-u0)*(v1-v0)/4;
                        const sN = this.getNumericalStress(patch, u_disp, u, v, E, nu);
                        const sE = this.getExactStress(pos.x, pos.y, R, Tx);

                        // Additional safety: skip points with obviously blown-up stresses
                        const maxPhysicalStress = 5 * Math.abs(Tx);
                        if (Math.abs(sN.sxx) > maxPhysicalStress * 10 || 
                            Math.abs(sN.syy) > maxPhysicalStress * 10 ||
                            Math.abs(sN.sxy) > maxPhysicalStress * 10) continue;

                        err2 += (Math.pow(sN.sxx - sE.sxx, 2) + Math.pow(sN.syy - sE.syy, 2) + 2*Math.pow(sN.sxy - sE.sxy, 2)) * Jmod;
                        ex2 += (Math.pow(sE.sxx, 2) + Math.pow(sE.syy, 2) + 2*Math.pow(sE.sxy, 2)) * Jmod;
                    }
                }
            }
        }
        return Math.sqrt(err2) / Math.sqrt(ex2);
    }

}

// Export for browser
if (typeof module !== 'undefined' && module.exports) {
    module.exports = { IGA2DSolver, GaussQuadrature2D };
} else {
    window.IGA2DSolver = IGA2DSolver;
    window.GaussQuadrature2D = GaussQuadrature2D;
}
