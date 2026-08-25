/**
 * Phase 5.3 - Reference Configuration Mapping Solver
 * Parallel physical (Naive) and pullback (Mapped) linear elastic structural dynamics.
 */

class MappingSolver {
    constructor(engine) {
        this.engine = engine;
        this.E = 200000;      // Young's Modulus (MPa)
        this.nu = 0.3;        // Poisson's Ratio
        this.thickness = 1.0;
        this.rho = 1.0;       // Damping/Mass density scaler

        // Cantilever Beam dimensions
        this.L = 10.0;
        this.H = 2.0;

        // Shape parameter limits
        this.mu = 5.0;        // Notch center
        this.r = 0.3;         // Notch depth
        this.sigma = 0.6;     // Notch width factor

        // State trackers
        this.referencePatch = null;
        this.precomputedGauss = null;

        // Dynamic Solver Parameters
        this.dt = 0.05;
        this.time = 0.0;
        this.alphaM = 0.08;   // Rayleigh mass damping
        this.betaK = 0.003;   // Rayleigh stiffness damping
        
        // Dynamic States
        this.uNaive = null;
        this.vNaive = null;
        this.aNaive = null;
        this.uMapped = null;
        this.vMapped = null;
        this.aMapped = null;

        // External Force
        this.forceAmp = 150.0;
        this.forceFreq = 1.2;
    }

    /**
     * Initialize reference patch and precomputed Gauss points
     */
    initializeReference(hLevel = 0, pLevel = 0) {
        console.log(`[Solver] Initializing Reference Configuration (h=${hLevel}, p=${pLevel})...`);
        // Reference patch is flat (r = 0)
        const flatPatch = window.GeometryFactory.generateNotchedBeam(this.L, this.H, this.mu, 0.0);
        this.referencePatch = window.GeometryFactory.refine(flatPatch, hLevel, pLevel);

        const patch = this.referencePatch;
        const { p, q, U, V, weights, controlPoints } = patch;
        const nU = controlPoints.length;
        const nV = controlPoints[0].length;
        const nDofs = nU * nV * 2;

        // Precompute Gauss quadrature evaluation data
        const uniqueU = [...new Set(U)];
        const uniqueV = [...new Set(V)];
        const gRule = GaussQuadrature2D.getPoints(Math.max(p, q) + 1);

        const list = [];

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
                        const v_val = ((vMax - vMin) * gRule.points[gv] + (vMax + vMin)) / 2;
                        const wv = gRule.weights[gv] * (vMax - vMin) / 2;

                        const detJ_ref = this.engine.getJacobianDeterminant(patch, u, v_val);
                        const deriv_ref = this.engine.getSurfaceDerivatives(patch, u, v_val);

                        // Reference physical coordinate (flat space)
                        const hx = deriv_ref.pos.x;
                        const hy = deriv_ref.pos.y;

                        // Basis gradients on reference patch
                        const grads = [];
                        const W = deriv_ref.W, Wu = deriv_ref.Wu, Wv = deriv_ref.Wv;

                        const J = [[deriv_ref.dU.x, deriv_ref.dV.x], [deriv_ref.dU.y, deriv_ref.dV.y]];
                        let detJ_2D = J[0][0] * J[1][1] - J[0][1] * J[1][0];
                        if (Math.abs(detJ_2D) < 1e-12) detJ_2D = (detJ_2D >= 0) ? 1e-12 : -1e-12;
                        const J_inv = [
                            [ J[1][1]/detJ_2D, -J[0][1]/detJ_2D],
                            [-J[1][0]/detJ_2D,  J[0][0]/detJ_2D]
                        ];

                        const spanU = this.engine.findSpan(nU - 1, p, u, U);
                        const spanV = this.engine.findSpan(nV - 1, q, v_val, V);

                        const dersU = this.engine.basisFunsDerivs(spanU, u, p, U, 1);
                        const dersV = this.engine.basisFunsDerivs(spanV, v_val, q, V, 1);

                        // Extract only active shape functions inside this knot span
                        const activeNodes = [];
                        for (let iu = 0; iu <= p; iu++) {
                            const cpI = spanU - p + iu;
                            const Ni = dersU[0][iu];
                            const dNi = dersU[1][iu];

                            for (let jv = 0; jv <= q; jv++) {
                                const cpJ = spanV - q + jv;
                                const Mj = dersV[0][jv];
                                const dMj = dersV[1][jv];

                                const w = weights[cpI][cpJ];
                                const dRdu = ((dNi * Mj * w) * W - (Ni * Mj * w) * Wu) / (W * W);
                                const dRdv = ((Ni * dMj * w) * W - (Ni * Mj * w) * Wv) / (W * W);

                                const dRdx = J_inv[0][0] * dRdu + J_inv[1][0] * dRdv;
                                const dRdy = J_inv[0][1] * dRdu + J_inv[1][1] * dRdv;

                                activeNodes.push({
                                    dofIdx: (cpI * nV + cpJ) * 2,
                                    nodeIdx: cpI * nV + cpJ,
                                    dRdx_ref: dRdx,
                                    dRdy_ref: dRdy
                                });
                            }
                        }

                        list.push({
                            hx, hy,
                            detJ_ref,
                            weightFactor: wu * wv,
                            activeNodes
                        });
                    }
                }
            }
        }

        this.precomputedGauss = list;
        this.resetDynamicStates();
    }

    resetDynamicStates() {
        const nDofs = this.referencePatch.controlPoints.length * this.referencePatch.controlPoints[0].length * 2;
        this.uNaive = new Float64Array(nDofs).fill(0);
        this.vNaive = new Float64Array(nDofs).fill(0);
        this.aNaive = new Float64Array(nDofs).fill(0);

        this.uMapped = new Float64Array(nDofs).fill(0);
        this.vMapped = new Float64Array(nDofs).fill(0);
        this.aMapped = new Float64Array(nDofs).fill(0);

        this.time = 0.0;
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
     * Physical stiffness assembly on morphed mesh Ω(μ, r)
     * Direct integration on physical coordinates
     */
    assembleNaiveStiffness(mu, r) {
        // Step 1: morph physical grid
        const morphedPatchBase = window.GeometryFactory.generateNotchedBeam(this.L, this.H, mu, r);
        // Find active refinement levels of the reference patch
        // Refine this morphed patch exactly to match
        const nU_ref = this.referencePatch.controlPoints.length;
        const nV_ref = this.referencePatch.controlPoints[0].length;
        
        let refinedPatch = morphedPatchBase;
        // Simple division logic based on control points
        const hSteps = Math.round(Math.log2((nU_ref - 3) / 18)); // Base nU = 21, base internal = 18.
        const pSteps = this.referencePatch.p - 2;

        refinedPatch = window.GeometryFactory.refine(morphedPatchBase, Math.max(0, hSteps), Math.max(0, pSteps));

        const solver = new window.IGA2DSolver(this.engine);
        solver.E = this.E;
        solver.nu = this.nu;
        solver.thickness = this.thickness;

        const K = solver.assembleStiffness(refinedPatch);
        solver.applyPenaltyConstraints(K, refinedPatch);
        return K;
    }

    /**
     * Mapped stiffness assembly via coordinates pullback integration on reference configuration \hat{Ω}
     */
    assembleMappedStiffness(mu, r) {
        const nU = this.referencePatch.controlPoints.length;
        const nV = this.referencePatch.controlPoints[0].length;
        const nDofs = nU * nV * 2;

        const K = Array.from({ length: nDofs }, () => new Float64Array(nDofs).fill(0));
        const D = this.getPlaneStressD();

        const H = this.H;
        const sigma = Math.max(r, 0.4);

        // Loop over precomputed reference Gauss points
        for (let gIdx = 0; gIdx < this.precomputedGauss.length; gIdx++) {
            const gp = this.precomputedGauss[gIdx];
            const { hx, hy, detJ_ref, weightFactor, activeNodes } = gp;

            // 1. Analytical coordinate mapping pullback formulas
            const dist = Math.abs(hx - mu);
            const disp = r * Math.exp(-0.5 * Math.pow(dist / sigma, 2));
            const dispDeriv = -((hx - mu) / (sigma * sigma)) * disp;

            const detJphi = 1.0 - (2.0 * disp / H);
            
            // Singular protection
            const safeDetJphi = Math.abs(detJphi) < 1e-12 ? 1e-12 : detJphi;

            const shearTerm = (2.0 / H) * hy * dispDeriv / safeDetJphi;
            const yTerm = 1.0 / safeDetJphi;

            // 2. Physical gradients calculation using coordinate pullback transpose Jacobian
            const B_phys = Array(nU * nV).fill(0).map(() => [[0, 0], [0, 0], [0, 0]]);

            for (let a = 0; a < activeNodes.length; a++) {
                const node = activeNodes[a];
                const dRdx_phys = node.dRdx_ref + shearTerm * node.dRdy_ref;
                const dRdy_phys = yTerm * node.dRdy_ref;

                const bidx = node.nodeIdx;
                B_phys[bidx][0][0] = dRdx_phys;
                B_phys[bidx][1][1] = dRdy_phys;
                B_phys[bidx][2][0] = dRdy_phys;
                B_phys[bidx][2][1] = dRdx_phys;
            }

            // 3. Integrate with analytical determinant pullback scale
            const factor = safeDetJphi * detJ_ref * weightFactor * this.thickness;

            // 4. Accumulate element contribution
            for (let a = 0; a < activeNodes.length; a++) {
                const na = activeNodes[a];
                const Ba = B_phys[na.nodeIdx];

                for (let b = 0; b < activeNodes.length; b++) {
                    const nb = activeNodes[b];
                    const Bb = B_phys[nb.nodeIdx];

                    for (let i = 0; i < 2; i++) {
                        for (let j = 0; j < 2; j++) {
                            let kab = 0;
                            for (let row = 0; row < 3; row++) {
                                for (let col = 0; col < 3; col++) {
                                    kab += Ba[row][i] * D[row][col] * Bb[col][j];
                                }
                            }
                            K[na.dofIdx + i][nb.dofIdx + j] += kab * factor;
                        }
                    }
                }
            }
        }

        // Apply boundary penalty stabilization
        const solver = new window.IGA2DSolver(this.engine);
        solver.applyPenaltyConstraints(K, this.referencePatch);

        return K;
    }

    /**
     * Compute a simple diagonal lumped mass matrix
     */
    assembleMassMatrix() {
        const nU = this.referencePatch.controlPoints.length;
        const nV = this.referencePatch.controlPoints[0].length;
        const nDofs = nU * nV * 2;

        const M = new Float64Array(nDofs).fill(0);
        
        // Total physical beam mass (approximate)
        const volume = this.L * this.H * this.thickness;
        const totalMass = this.rho * volume;
        const nodalMass = totalMass / (nU * nV);

        for (let i = 0; i < nDofs; i++) {
            M[i] = nodalMass;
        }

        return M;
    }

    /**
     * Apply Dirichlet Clamping boundary conditions (left edge x=0)
     */
    applyClamping(K_eff, R) {
        const nV = this.referencePatch.controlPoints[0].length;
        // Fix all DOFs for control points along i=0
        for (let j = 0; j < nV; j++) {
            const baseDof = j * 2; // (i * nV + j) * 2 with i=0
            
            // X DOF
            K_eff[baseDof][baseDof] = 1e30;
            R[baseDof] = 0.0;

            // Y DOF
            K_eff[baseDof + 1][baseDof + 1] = 1e30;
            R[baseDof + 1] = 0.0;
        }
    }

    /**
     * Assemble dynamic forces at right edge (axial and bending tension)
     */
    calculateExternalForceVector(t) {
        const nU = this.referencePatch.controlPoints.length;
        const nV = this.referencePatch.controlPoints[0].length;
        const nDofs = nU * nV * 2;
        const F = new Float64Array(nDofs).fill(0);

        // Apply a dynamic shear harmonic load on the right edge
        const dynamicLoad = this.forceAmp * Math.sin(this.forceFreq * t);

        for (let j = 0; j < nV; j++) {
            const idx = ((nU - 1) * nV + j) * 2; // Right edge i = nU - 1
            F[idx + 1] = dynamicLoad / nV; // Divide force equally among V nodes
        }
        return F;
    }

    /**
     * Implicit Newmark-beta Transient Dynamics Step
     */
    stepDynamics(mu, r, runNaive = true, runMapped = true) {
        this.time += this.dt;

        // Step 1: Assemble Stiffness matrices
        let tNaiveAssembly = 0;
        let KNaive = null;
        if (runNaive) {
            const tNaiveStart = performance.now();
            KNaive = this.assembleNaiveStiffness(mu, r);
            tNaiveAssembly = performance.now() - tNaiveStart;
        }

        let tMappedAssembly = 0;
        let KMapped = null;
        if (runMapped) {
            const tMappedStart = performance.now();
            KMapped = this.assembleMappedStiffness(mu, r);
            tMappedAssembly = performance.now() - tMappedStart;
        }

        // Step 2: Mass & Damping matrices
        const M = this.assembleMassMatrix();
        const nDofs = M.length;

        // Implicit parameters
        const beta = 0.25;
        const gamma = 0.5;

        // Solve dynamics for both engines
        const solveEngine = (K, u, v, a) => {
            // Predictor phase
            const uPred = new Float64Array(nDofs);
            const vPred = new Float64Array(nDofs);
            for (let i = 0; i < nDofs; i++) {
                uPred[i] = u[i] + this.dt * v[i] + (this.dt * this.dt) * (0.5 - beta) * a[i];
                vPred[i] = v[i] + this.dt * (1.0 - gamma) * a[i];
            }

            // Assemble Damping matrix C = alpha * M + beta * K
            const C_uPred = new Float64Array(nDofs);
            const K_uPred = new Float64Array(nDofs);

            for (let i = 0; i < nDofs; i++) {
                let k_up = 0;
                for (let j = 0; j < nDofs; j++) k_up += K[i][j] * uPred[j];
                K_uPred[i] = k_up;
                C_uPred[i] = this.alphaM * M[i] * vPred[i] + this.betaK * k_up; // Rayleigh damping damping
            }

            // Effective Force
            const F_ext = this.calculateExternalForceVector(this.time);
            const R = new Float64Array(nDofs);
            for (let i = 0; i < nDofs; i++) {
                R[i] = F_ext[i] - C_uPred[i] - K_uPred[i] - M[i] * a[i] * (1.0 - 2.0 * beta) / (2.0 * beta); // standard prediction
            }

            // Effective Tangent Matrix K_eff = K * (1 + gamma * betaK / (beta * dt)) + M * (1 / (beta * dt^2) + gamma * alphaM / (beta * dt))
            const K_eff = Array.from({ length: nDofs }, () => new Float64Array(nDofs).fill(0));
            const factorK = 1.0 + (gamma * this.betaK) / (beta * this.dt);
            const factorM = 1.0 / (beta * this.dt * this.dt) + (gamma * this.alphaM) / (beta * this.dt);

            for (let i = 0; i < nDofs; i++) {
                for (let j = 0; j < nDofs; j++) {
                    K_eff[i][j] = K[i][j] * factorK;
                }
                K_eff[i][i] += M[i] * factorM;
            }

            // Apply clamps BCs
            this.applyClamping(K_eff, R);

            // Solve system
            const solver = new window.IGA2DSolver(this.engine);
            const du = solver.gaussianElimination(K_eff, R);

            // Corrector phase
            const uNext = new Float64Array(nDofs);
            const aNext = new Float64Array(nDofs);
            const vNext = new Float64Array(nDofs);

            for (let i = 0; i < nDofs; i++) {
                uNext[i] = uPred[i] + du[i];
                aNext[i] = (uNext[i] - u[i] - this.dt * v[i]) / (beta * this.dt * this.dt) - ((1.0 - 2.0 * beta) / (2.0 * beta)) * a[i];
                vNext[i] = vPred[i] + gamma * this.dt * aNext[i];
            }

            return { u: uNext, v: vNext, a: aNext };
        };

        if (runNaive && KNaive) {
            const resNaive = solveEngine(KNaive, this.uNaive, this.vNaive, this.aNaive);
            this.uNaive = resNaive.u;
            this.vNaive = resNaive.v;
            this.aNaive = resNaive.a;
        }

        if (runMapped && KMapped) {
            const resMapped = solveEngine(KMapped, this.uMapped, this.vMapped, this.aMapped);
            this.uMapped = resMapped.u;
            this.vMapped = resMapped.v;
            this.aMapped = resMapped.a;
        }

        return {
            time: this.time,
            tNaiveAssembly,
            tMappedAssembly,
            uNaive: this.uNaive,
            uMapped: this.uMapped
        };
    }
}

window.MappingSolver = MappingSolver;
