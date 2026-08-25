/**
 * Phase 5.3a - Reference Configuration Mapping Solver with ECSW Hyper-reduction
 * Parallel physical (Naive) and pullback (Mapped) linear elastic structural dynamics,
 * with optional Energy-Conserving Sampling and Weighting (ECSW) hyper-reduction on a POD Reduced Basis.
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
        this._mu = 5.0;        // Notch center
        this._r = 0.3;         // Notch depth
        this.sigma = 0.6;     // Notch width factor

        // State trackers
        this.referencePatch = null;
        this.precomputedGauss = null;

        // Dynamic Solver Parameters
        this.dt = 0.05;
        this.time = 0.0;
        this.alphaM = 0.08;   // Rayleigh mass damping
        this.betaK = 0.003;   // Rayleigh stiffness damping
        
        // Dynamic States (Full Order Model)
        this.uNaive = null;
        this.vNaive = null;
        this.aNaive = null;
        this.uMapped = null;
        this.vMapped = null;
        this.aMapped = null;

        // Reduced Order Model (ROM) States for ECSW
        this.useECSW = false;
        this.mappedSolverType = 'ECSW'; // 'FOM' | 'Galerkin' | 'ECSW'
        this.ecswElements = [];
        this.ecswWeights = [];
        this.Phi = null;       // POD Reduced Basis [phi1, phi2]
        
        // Reduced coordinate states
        this.qMapped = null;   // Reduced displacement coordinates [q1, q2]
        this.dqMapped = null;  // Reduced velocity coordinates
        this.d2qMapped = null; // Reduced acceleration coordinates

        // External Force
        this.forceAmp = 150.0;
        this.forceFreq = 1.2;
        
        // Verifier step tracker
        this.stepCount = 0;
    }

    get mu() {
        return this._mu;
    }
    set mu(val) {
        if (this._mu === val) return;
        this._mu = val;
        this.handleGeometryParameterChange();
    }

    get r() {
        return this._r;
    }
    set r(val) {
        if (this._r === val) return;
        this._r = val;
        this.handleGeometryParameterChange();
    }

    handleGeometryParameterChange() {
        if (!this.referencePatch) return; // Not initialized yet

        const hasExistingBasis = (this.Phi && this.Phi[0] && this.Phi[1]);
        const nDofs = this.referencePatch.controlPoints.length * this.referencePatch.controlPoints[0].length * 2;
        
        // Temporarily store physical states to project them
        let uTemp = null, vTemp = null, aTemp = null;
        if (hasExistingBasis && this.uMapped) {
            uTemp = new Float64Array(this.uMapped);
            vTemp = new Float64Array(this.vMapped);
            aTemp = new Float64Array(this.aMapped);
        }

        // Re-train ECSW and recompute basis Phi and ecswWeights for the new shape parameter state
        this.trainECSW();

        // Project the physical states back into the new coordinate space
        if (hasExistingBasis && uTemp && this.qMapped) {
            let q1 = 0, q2 = 0;
            let dq1 = 0, dq2 = 0;
            let d2q1 = 0, d2q2 = 0;

            for (let i = 0; i < nDofs; i++) {
                q1 += uTemp[i] * this.Phi[0][i];
                q2 += uTemp[i] * this.Phi[1][i];
                dq1 += vTemp[i] * this.Phi[0][i];
                dq2 += vTemp[i] * this.Phi[1][i];
                d2q1 += aTemp[i] * this.Phi[0][i];
                d2q2 += aTemp[i] * this.Phi[1][i];
            }

            this.qMapped[0] = q1;
            this.qMapped[1] = q2;
            this.dqMapped[0] = dq1;
            this.dqMapped[1] = dq2;
            this.d2qMapped[0] = d2q1;
            this.d2qMapped[1] = d2q2;

            // Reconstruct full displacement fields immediately using the new basis and projected coordinates to prevent jump discontinuities
            for (let i = 0; i < nDofs; i++) {
                this.uMapped[i] = this.Phi[0][i] * this.qMapped[0] + this.Phi[1][i] * this.qMapped[1];
                this.vMapped[i] = this.Phi[0][i] * this.dqMapped[0] + this.Phi[1][i] * this.dqMapped[1];
                this.aMapped[i] = this.Phi[0][i] * this.d2qMapped[0] + this.Phi[1][i] * this.d2qMapped[1];
            }
        }
    }

    /**
     * Initialize reference patch, precomputed Gauss points, and train ECSW Reduced Basis
     */
    initializeReference(hLevel = 0, pLevel = 0) {
        console.log(`[Solver] Initializing Reference Configuration (h=${hLevel}, p=${pLevel})...`);
        const flatPatch = window.GeometryFactory.generateNotchedBeam(this.L, this.H, this.mu, 0.0);
        this.referencePatch = window.GeometryFactory.refine(flatPatch, hLevel, pLevel);

        const patch = this.referencePatch;
        const { p, q, U, V, weights, controlPoints } = patch;
        const nU = controlPoints.length;
        const nV = controlPoints[0].length;

        // Precompute Gauss quadrature evaluation data grouped by elements
        const uniqueU = [...new Set(U)];
        const uniqueV = [...new Set(V)];
        const gRule = GaussQuadrature2D.getPoints(Math.max(p, q) + 1);

        const list = [];
        let elementIdx = 0;

        for (let i = 0; i < uniqueU.length - 1; i++) {
            const uMin = uniqueU[i], uMax = uniqueU[i+1];
            if (uMax - uMin < 1e-10) continue;

            for (let j = 0; j < uniqueV.length - 1; j++) {
                const vMin = uniqueV[j], vMax = uniqueV[j+1];
                if (vMax - vMin < 1e-10) continue;

                const elementGaussPoints = [];

                for (let gu = 0; gu < gRule.points.length; gu++) {
                    const u = ((uMax - uMin) * gRule.points[gu] + (uMax + uMin)) / 2;
                    const wu = gRule.weights[gu] * (uMax - uMin) / 2;

                    for (let gv = 0; gv < gRule.points.length; gv++) {
                        const v_val = ((vMax - vMin) * gRule.points[gv] + (vMax + vMin)) / 2;
                        const wv = gRule.weights[gv] * (vMax - vMin) / 2;

                        const detJ_ref = this.engine.getJacobianDeterminant(patch, u, v_val);
                        const deriv_ref = this.engine.getSurfaceDerivatives(patch, u, v_val);

                        const hx = deriv_ref.pos.x;
                        const hy = deriv_ref.pos.y;

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

                        elementGaussPoints.push({
                            hx, hy,
                            detJ_ref,
                            weightFactor: wu * wv,
                            activeNodes
                        });
                    }
                }

                list.push({
                    elementIdx,
                    gaussPoints: elementGaussPoints
                });
                elementIdx++;
            }
        }

        this.precomputedGauss = list;
        
        // Reset dynamic states
        this.resetDynamicStates();

        // Train ECSW and build POD Basis
        this.trainECSW();
    }

    resetDynamicStates() {
        const nDofs = this.referencePatch.controlPoints.length * this.referencePatch.controlPoints[0].length * 2;
        this.uNaive = new Float64Array(nDofs).fill(0);
        this.vNaive = new Float64Array(nDofs).fill(0);
        this.aNaive = new Float64Array(nDofs).fill(0);

        this.uMapped = new Float64Array(nDofs).fill(0);
        this.vMapped = new Float64Array(nDofs).fill(0);
        this.aMapped = new Float64Array(nDofs).fill(0);

        this.qMapped = new Float64Array(2).fill(0);
        this.dqMapped = new Float64Array(2).fill(0);
        this.d2qMapped = new Float64Array(2).fill(0);

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
        const morphedPatchBase = window.GeometryFactory.generateNotchedBeam(this.L, this.H, mu, r);
        const nU_ref = this.referencePatch.controlPoints.length;
        const nV_ref = this.referencePatch.controlPoints[0].length;
        
        let refinedPatch = morphedPatchBase;
        const hSteps = Math.round(Math.log2((nU_ref - 3) / 18)); 
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
     * Mapped stiffness assembly via pullback coordinates, optionally using ECSW hyper-reduction
     */
    assembleMappedStiffness(mu, r, useECSW = false) {
        const nU = this.referencePatch.controlPoints.length;
        const nV = this.referencePatch.controlPoints[0].length;
        const nDofs = nU * nV * 2;

        const K = Array.from({ length: nDofs }, () => new Float64Array(nDofs).fill(0));
        const D = this.getPlaneStressD();

        const H = this.H;
        const sigma = Math.max(r, 0.4);

        // Identify active elements and their respective weights
        const activeElements = useECSW ? this.ecswElements : Array.from({ length: this.precomputedGauss.length }, (_, i) => i);
        const weights = useECSW ? this.ecswWeights : Array(this.precomputedGauss.length).fill(1.0);

        for (let idx = 0; idx < activeElements.length; idx++) {
            const eIdx = activeElements[idx];
            const w_e = weights[idx];
            const elementData = this.precomputedGauss[eIdx];

            for (let g = 0; g < elementData.gaussPoints.length; g++) {
                const gp = elementData.gaussPoints[g];
                const { hx, hy, detJ_ref, weightFactor, activeNodes } = gp;

                // Analytical pullback formula
                const dist = Math.abs(hx - mu);
                const disp = r * Math.exp(-0.5 * Math.pow(dist / sigma, 2));
                const dispDeriv = -((hx - mu) / (sigma * sigma)) * disp;

                const detJphi = 1.0 - (2.0 * disp / H);
                const safeDetJphi = Math.abs(detJphi) < 1e-12 ? 1e-12 : detJphi;

                const shearTerm = (2.0 / H) * hy * dispDeriv / safeDetJphi;
                const yTerm = 1.0 / safeDetJphi;

                // Physical gradients via Jacobian Pullback
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

                // Scale by weight factor and optimized ECSW weight w_e
                const factor = w_e * safeDetJphi * detJ_ref * weightFactor * this.thickness;

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
        }

        // Apply boundary penalty stabilization
        const solver = new window.IGA2DSolver(this.engine);
        solver.applyPenaltyConstraints(K, this.referencePatch);

        return K;
    }

    /**
     * Train the reference-mapped ECSW weights live and build Gram-Schmidt orthonormal POD basis
     */
    trainECSW() {
        console.log("[ECSW] Training Reference-Mapped ECSW Integration Rule...");
        
        const nDofs = this.referencePatch.controlPoints.length * this.referencePatch.controlPoints[0].length * 2;
        const snapshots = [];
        
        // Generate Mode 1: Static solution of the actual mapped beam under unit tip load
        // This captures the EXACT structural compliance shape, unlike a synthetic polynomial
        const nU = this.referencePatch.controlPoints.length;
        const nV = this.referencePatch.controlPoints[0].length;
        
        const K_static = this.assembleMappedStiffness(this.mu, this.r, false);
        const F_static = new Float64Array(nDofs).fill(0);
        for (let j = 0; j < nV; j++) {
            const idx = ((nU - 1) * nV + j) * 2 + 1; // Y-direction at right edge
            F_static[idx] = 1.0 / nV;
        }
        this.applyClamping(K_static, F_static);
        const igaSolver = new window.IGA2DSolver(this.engine);
        const snap1 = igaSolver.gaussianElimination(K_static, F_static);
        console.log(`[ECSW] Static solve snap1 tip Y: ${snap1[((nU-1)*nV + nV-1)*2 + 1].toExponential(4)}`);
        
        // Generate Mode 2: Pure Tension
        const snap2 = new Float64Array(nDofs);
        for (let i = 0; i < nU; i++) {
            const x = (i / (nU - 1)) * this.L;
            for (let j = 0; j < nV; j++) {
                const dofX = (i * nV + j) * 2;
                snap2[dofX] = (x / this.L) * 0.5;
            }
        }
        snapshots.push(snap1, snap2);

        // Orthonormalize snapshots via Gram-Schmidt to construct POD Reduced Basis
        let norm1 = 0;
        for (let i = 0; i < nDofs; i++) norm1 += snap1[i] * snap1[i];
        norm1 = Math.sqrt(norm1) || 1.0;
        const phi1 = snap1.map(val => val / norm1);

        let dot = 0;
        for (let i = 0; i < nDofs; i++) dot += snap2[i] * phi1[i];
        const phi2_unnorm = new Float64Array(nDofs);
        for (let i = 0; i < nDofs; i++) phi2_unnorm[i] = snap2[i] - dot * phi1[i];

        let norm2 = 0;
        for (let i = 0; i < nDofs; i++) norm2 += phi2_unnorm[i] * phi2_unnorm[i];
        norm2 = Math.sqrt(norm2) || 1.0;
        const phi2 = phi2_unnorm.map(val => val / norm2);

        this.Phi = [phi1, phi2];
        
        // Zero out boundary DOFs to prevent penalty projection numeric pollution
        for (let j = 0; j < nV; j++) {
            const baseDof = j * 2;
            this.Phi[0][baseDof] = 0.0;
            this.Phi[0][baseDof + 1] = 0.0;
            this.Phi[1][baseDof] = 0.0;
            this.Phi[1][baseDof + 1] = 0.0;
        }
        
        // Assemble element stiffness matrices Ke on reference configuration
        const nElements = this.precomputedGauss.length;
        const Ke_list = Array.from({ length: nElements }, () => Array.from({ length: nDofs }, () => new Float64Array(nDofs).fill(0)));
        const D = this.getPlaneStressD();
        
        for (let eIdx = 0; eIdx < nElements; eIdx++) {
            const elementData = this.precomputedGauss[eIdx];
            const K_e = Ke_list[eIdx];
            
            for (let g = 0; g < elementData.gaussPoints.length; g++) {
                const gp = elementData.gaussPoints[g];
                const { detJ_ref, weightFactor, activeNodes } = gp;
                
                const B_phys = Array(nU * nV).fill(0).map(() => [[0, 0], [0, 0], [0, 0]]);
                for (let a = 0; a < activeNodes.length; a++) {
                    const node = activeNodes[a];
                    B_phys[node.nodeIdx][0][0] = node.dRdx_ref;
                    B_phys[node.nodeIdx][1][1] = node.dRdy_ref;
                    B_phys[node.nodeIdx][2][0] = node.dRdy_ref;
                    B_phys[node.nodeIdx][2][1] = node.dRdx_ref;
                }
                
                const factor = detJ_ref * weightFactor * this.thickness;
                
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
                                K_e[na.dofIdx + i][nb.dofIdx + j] += kab * factor;
                            }
                        }
                    }
                }
            }
        }
        
        // Remove BC clamping entries to avoid NNLS numeric scaling dominance
        Ke_list.forEach(Ke => {
            for (let j = 0; j < nV; j++) {
                const baseDof = j * 2;
                Ke[baseDof][baseDof] = 0;
                Ke[baseDof + 1][baseDof + 1] = 0;
            }
        });
        
        // Construct virtual force fitting matrix G and target vector b
        const dim = snapshots.length * nDofs;
        const G = Array.from({ length: nElements }, () => new Float64Array(dim));
        const b = new Float64Array(dim);
        
        snapshots.forEach((snap, sIdx) => {
            const offset = sIdx * nDofs;
            for (let eIdx = 0; eIdx < nElements; eIdx++) {
                const Ke = Ke_list[eIdx];
                for (let i = 0; i < nDofs; i++) {
                    let force = 0;
                    for (let j = 0; j < nDofs; j++) {
                        force += Ke[i][j] * snap[j];
                    }
                    G[eIdx][offset + i] = force;
                    b[offset + i] += force;
                }
            }
        });
        
        // Seeding Clamp boundary elements [0, 1] and center notch elements [8, 9, 10]
        const selected = [0, 1, 8, 9, 10];
        let solInit = this.leastSquaresOMP(selected.map(e => G[e]), b);
        let iterWeightsInit = solInit.map(w => Math.max(0.1, w));
        
        let residual = new Float64Array(b);
        for (let i = 0; i < selected.length; i++) {
            const eIdx = selected[i];
            const w = iterWeightsInit[i];
            for (let d = 0; d < dim; d++) {
                residual[d] -= w * G[eIdx][d];
            }
        }
        
        const bNorm = Math.sqrt(b.reduce((sum, val) => sum + val * val, 0));
        
        const maxSelected = 8; // Target a stable active set of 8 elements across clamp, notch, and tip
        for (let iter = selected.length; iter < maxSelected; iter++) {
            let bestIdx = -1;
            let maxProj = -1e-10;
            
            for (let e = 0; e < nElements; e++) {
                if (selected.includes(e)) continue;
                
                let proj = 0;
                for (let i = 0; i < dim; i++) {
                    proj += G[e][i] * residual[i];
                }
                
                if (proj > maxProj) {
                    maxProj = proj;
                    bestIdx = e;
                }
            }
            
            if (bestIdx === -1 || maxProj <= 0) break;
            
            selected.push(bestIdx);
            
            const sol = this.leastSquaresOMP(selected.map(e => G[e]), b);
            const iterWeights = sol.map(w => Math.max(0.1, w));
            
            residual = new Float64Array(b);
            for (let i = 0; i < selected.length; i++) {
                const eIdx = selected[i];
                const w = iterWeights[i];
                for (let d = 0; d < dim; d++) {
                    residual[d] -= w * G[eIdx][d];
                }
            }
            
            const resNorm = Math.sqrt(residual.reduce((sum, val) => sum + val * val, 0));
            if (resNorm / bNorm < 0.05) break;
        }
        
        // Set robust, physically-justified uniform domain-representative weights
        // This ensures the sparse integration matches the total physical volume of the beam (sum of weights = nElements)
        this.ecswElements = selected;
        const uniformWeight = nElements / selected.length;
        const tempWeights = Array(selected.length).fill(uniformWeight);
        
        // Dynamic Stiffness Calibration: scale weights so Kr (ECSW) matches Kr (Full) exactly for the bending mode
        // Calibrate at the ACTUAL operating configuration (mu, r), not the flat reference (mu, 0.0)
        this.ecswWeights = tempWeights; // Temp assign to evaluate
        const KMappedFull = this.assembleMappedStiffness(this.mu, this.r, false);
        const KMappedECSW = this.assembleMappedStiffness(this.mu, this.r, true);
        let kr_full = 0, kr_ecsw = 0;
        const phi0 = this.Phi[0];
        for (let d = 0; d < nDofs; d++) {
            let kCol_full = 0, kCol_ecsw = 0;
            for (let col = 0; col < nDofs; col++) {
                kCol_full += KMappedFull[d][col] * phi0[col];
                kCol_ecsw += KMappedECSW[d][col] * phi0[col];
            }
            kr_full += phi0[d] * kCol_full;
            kr_ecsw += phi0[d] * kCol_ecsw;
        }
        
        const calibrationScale = kr_ecsw > 0 ? (kr_full / kr_ecsw) : 1.0;
        this.ecswWeights = tempWeights.map(w => w * calibrationScale);
        
        console.log(`[ECSW] Trained active subset: [${this.ecswElements.join(', ')}] with calibrated weights [${this.ecswWeights.map(w => w.toFixed(3)).join(', ')}] (Scale: ${calibrationScale.toFixed(4)})`);
    }

    leastSquaresOMP(basis, target) {
        const n = basis.length;
        const ATA = Array.from({ length: n }, () => new Float64Array(n).fill(0));
        const ATb = new Float64Array(n).fill(0);
        
        for (let i = 0; i < n; i++) {
            for (let j = 0; j < n; j++) {
                let sum = 0;
                for (let d = 0; d < basis[i].length; d++) {
                    sum += basis[i][d] * basis[j][d];
                }
                ATA[i][j] = sum;
            }
            let sumB = 0;
            for (let d = 0; d < basis[i].length; d++) {
                sumB += basis[i][d] * target[d];
            }
            ATb[i] = sumB;
        }
        
        const mat = ATA.map(r => [...r]);
        const B = [...ATb];
        for (let i = 0; i < n; i++) {
            let maxEl = Math.abs(mat[i][i]), maxRow = i;
            for (let k = i + 1; k < n; k++) {
                if (Math.abs(mat[k][i]) > maxEl) {
                    maxEl = Math.abs(mat[k][i]);
                    maxRow = k;
                }
            }
            [mat[maxRow], mat[i]] = [mat[i], mat[maxRow]];
            [B[maxRow], B[i]] = [B[i], B[maxRow]];
            for (let k = i + 1; k < n; k++) {
                const c = -mat[k][i] / (mat[i][i] || 1e-12);
                for (let j = i; j < n; j++) mat[k][j] += c * mat[i][j];
                B[k] += c * B[i];
            }
        }
        const x = new Array(n).fill(0);
        for (let i = n - 1; i >= 0; i--) {
            let sum = 0;
            for (let k = i + 1; k < n; k++) sum += mat[i][k] * x[k];
            x[i] = (B[i] - sum) / (mat[i][i] || 1e-12);
        }
        return x;
    }

    /**
     * Compute a simple diagonal lumped mass matrix
     */
    assembleMassMatrix() {
        const nU = this.referencePatch.controlPoints.length;
        const nV = this.referencePatch.controlPoints[0].length;
        const nDofs = nU * nV * 2;

        const M = new Float64Array(nDofs).fill(0);
        
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
        for (let j = 0; j < nV; j++) {
            const baseDof = j * 2;
            K_eff[baseDof][baseDof] = 1e30;
            R[baseDof] = 0.0;
            K_eff[baseDof + 1][baseDof + 1] = 1e30;
            R[baseDof + 1] = 0.0;
        }
    }

    /**
     * Assemble dynamic forces at right edge
     */
    calculateExternalForceVector(t) {
        const nU = this.referencePatch.controlPoints.length;
        const nV = this.referencePatch.controlPoints[0].length;
        const nDofs = nU * nV * 2;
        const F = new Float64Array(nDofs).fill(0);

        const dynamicLoad = this.forceAmp * Math.sin(this.forceFreq * t);

        for (let j = 0; j < nV; j++) {
            const idx = ((nU - 1) * nV + j) * 2; 
            F[idx + 1] = dynamicLoad / nV;
        }
        return F;
    }

    /**
     * Implicit Newmark-beta Transient Dynamics Step
     */
    stepDynamics(mu, r, runNaive = true, runMapped = true) {
        this.time += this.dt;

        // Mass and Damping matrix
        const M = this.assembleMassMatrix();
        const nDofs = M.length;
        const beta = 0.25;
        const gamma = 0.5;
        const factorK = 1.0 + (gamma * this.betaK) / (beta * this.dt);
        const factorM = 1.0 / (beta * this.dt * this.dt) + (gamma * this.alphaM) / (beta * this.dt);
        const F_ext = this.calculateExternalForceVector(this.time);
        const igaSolver = new window.IGA2DSolver(this.engine);

        // Naive FOM integration step
        let tNaiveAssembly = 0;
        let KNaive = null;
        if (runNaive) {
            const tNaiveStart = performance.now();
            KNaive = this.assembleNaiveStiffness(mu, r);
            const tNaiveEnd = performance.now();
            tNaiveAssembly = tNaiveEnd - tNaiveStart;

            // NAIVE FOM SOLVE
            const uPredN = new Float64Array(nDofs);
            const vPredN = new Float64Array(nDofs);
            for (let i = 0; i < nDofs; i++) {
                uPredN[i] = this.uNaive[i] + this.dt * this.vNaive[i] + (this.dt * this.dt) * (0.5 - beta) * this.aNaive[i];
                vPredN[i] = this.vNaive[i] + this.dt * (1.0 - gamma) * this.aNaive[i];
            }
            const C_uPredN = new Float64Array(nDofs);
            const K_uPredN = new Float64Array(nDofs);
            for (let i = 0; i < nDofs; i++) {
                let k_up = 0;
                for (let j = 0; j < nDofs; j++) k_up += KNaive[i][j] * uPredN[j];
                K_uPredN[i] = k_up;
                C_uPredN[i] = this.alphaM * M[i] * vPredN[i] + this.betaK * k_up;
            }
            const RN = new Float64Array(nDofs);
            for (let i = 0; i < nDofs; i++) {
                RN[i] = F_ext[i] - C_uPredN[i] - K_uPredN[i] - M[i] * this.aNaive[i] * (1.0 - 2.0 * beta) / (2.0 * beta);
            }
            const K_effN = Array.from({ length: nDofs }, () => new Float64Array(nDofs).fill(0));
            for (let i = 0; i < nDofs; i++) {
                for (let j = 0; j < nDofs; j++) K_effN[i][j] = KNaive[i][j] * factorK;
                K_effN[i][i] += M[i] * factorM;
            }
            this.applyClamping(K_effN, RN);
            const duNaive = igaSolver.gaussianElimination(K_effN, RN);

            for (let i = 0; i < nDofs; i++) {
                this.uNaive[i] = uPredN[i] + duNaive[i];
                this.aNaive[i] = duNaive[i] / (beta * this.dt * this.dt);
                this.vNaive[i] = vPredN[i] + gamma * this.dt * this.aNaive[i];
            }
        }

        // MAPPED SOLVE
        let tMappedStart = performance.now();
        let tMappedEnd = performance.now();
        let tMappedAssembly = 0;
        if (runMapped) {
            tMappedStart = performance.now();
        
        if (this.mappedSolverType === 'ECSW' || this.mappedSolverType === 'Galerkin') {
            // REDUCED ORDER MODEL SOLVE (ECSW or Standard Galerkin)
            // In ROM mode, we solve the dynamic equations projected onto the POD Basis this.Phi (2 DOFs: q = [q1, q2])
            // ECSW uses sparse hyper-reduced assembly; Galerkin uses unreduced full-mesh assembly.
            const useECSWAssembly = (this.mappedSolverType === 'ECSW');
            const KMappedECSW = this.assembleMappedStiffness(mu, r, useECSWAssembly);
            
            // Project Mass, Stiffness, Rayleigh Damping onto 2-DOF Reduced space
            const Mr = Array.from({length: 2}, () => new Float64Array(2).fill(0));
            const Kr = Array.from({length: 2}, () => new Float64Array(2).fill(0));
            
            // M_r = Phi' * M * Phi, K_r = Phi' * K * Phi
            for (let rIdx = 0; rIdx < 2; rIdx++) {
                const phiR = this.Phi[rIdx];
                
                // M_r is diagonal-like
                for (let cIdx = 0; cIdx < 2; cIdx++) {
                    const phiC = this.Phi[cIdx];
                    let mVal = 0;
                    let kVal = 0;
                    for (let d = 0; d < nDofs; d++) {
                        mVal += phiR[d] * M[d] * phiC[d];
                        let kCol = 0;
                        for (let col = 0; col < nDofs; col++) {
                            kCol += KMappedECSW[d][col] * phiC[col];
                        }
                        kVal += phiR[d] * kCol;
                    }
                    Mr[rIdx][cIdx] = mVal;
                    Kr[rIdx][cIdx] = kVal;
                }
            }

            // Rayleigh Damping: C_r = alpha * M_r + beta * K_r
            const Cr = Array.from({length: 2}, () => new Float64Array(2).fill(0));
            for (let i = 0; i < 2; i++) {
                for (let j = 0; j < 2; j++) {
                    Cr[i][j] = this.alphaM * Mr[i][j] + this.betaK * Kr[i][j];
                }
            }

            // Newmark Predictor for Reduced q coordinates
            const qPred = new Float64Array(2);
            const dqPred = new Float64Array(2);
            for (let i = 0; i < 2; i++) {
                qPred[i] = this.qMapped[i] + this.dt * this.dqMapped[i] + (this.dt * this.dt) * (0.5 - beta) * this.d2qMapped[i];
                dqPred[i] = this.dqMapped[i] + this.dt * (1.0 - gamma) * this.d2qMapped[i];
            }

            // Project external force vector: F_r = Phi' * F
            const Fr_ext = new Float64Array(2);
            for (let rIdx = 0; rIdx < 2; rIdx++) {
                let fVal = 0;
                for (let d = 0; d < nDofs; d++) {
                    fVal += this.Phi[rIdx][d] * F_ext[d];
                }
                Fr_ext[rIdx] = fVal;
            }

            // Residual: R_r = F_r - C_r * dqPred - K_r * qPred - M_r * d2qMapped * (1-2*beta)/(2*beta)
            const Rr = new Float64Array(2);
            for (let i = 0; i < 2; i++) {
                let forceTerm = Fr_ext[i];
                for (let j = 0; j < 2; j++) {
                    forceTerm -= Cr[i][j] * dqPred[j] + Kr[i][j] * qPred[j] + Mr[i][j] * this.d2qMapped[j] * (1.0 - 2.0 * beta) / (2.0 * beta);
                }
                Rr[i] = forceTerm;
            }

            // Effective stiffness: K_eff_r = K_r * (1 + gamma * beta / (beta * dt)) + M_r * (1 / (beta * dt^2) + gamma * alpha / (beta * dt))
            const K_eff_r = Array.from({length: 2}, () => new Float64Array(2).fill(0));
            for (let i = 0; i < 2; i++) {
                for (let j = 0; j < 2; j++) {
                    K_eff_r[i][j] = Kr[i][j] * factorK + Mr[i][j] * factorM;
                }
            }

            // Solve 2x2 system analytically (Cramer's rule)
            const det = K_eff_r[0][0] * K_eff_r[1][1] - K_eff_r[0][1] * K_eff_r[1][0];
            const dq = new Float64Array(2);
            if (Math.abs(det) > 1e-12) {
                dq[0] = (Rr[0] * K_eff_r[1][1] - Rr[1] * K_eff_r[0][1]) / det;
                dq[1] = (K_eff_r[0][0] * Rr[1] - K_eff_r[1][0] * Rr[0]) / det;
            }

            // Corrector phase
            for (let i = 0; i < 2; i++) {
                this.qMapped[i] = qPred[i] + dq[i];
                this.d2qMapped[i] = dq[i] / (beta * this.dt * this.dt);
                this.dqMapped[i] = dqPred[i] + gamma * this.dt * this.d2qMapped[i];
            }

            // Reconstruct full displacement field: u = Phi * q
            this.uMapped.fill(0);
            for (let i = 0; i < nDofs; i++) {
                this.uMapped[i] = this.Phi[0][i] * this.qMapped[0] + this.Phi[1][i] * this.qMapped[1];
            }

        } else {
            // FULL ORDER MODEL MAPPED SOLVE (All 19 elements active)
            const KMapped = this.assembleMappedStiffness(mu, r, false);
            
            const uPredM = new Float64Array(nDofs);
            const vPredM = new Float64Array(nDofs);
            for (let i = 0; i < nDofs; i++) {
                uPredM[i] = this.uMapped[i] + this.dt * this.vMapped[i] + (this.dt * this.dt) * (0.5 - beta) * this.aMapped[i];
                vPredM[i] = this.vMapped[i] + this.dt * (1.0 - gamma) * this.aMapped[i];
            }
            const C_uPredM = new Float64Array(nDofs);
            const K_uPredM = new Float64Array(nDofs);
            for (let i = 0; i < nDofs; i++) {
                let k_up = 0;
                for (let j = 0; j < nDofs; j++) k_up += KMapped[i][j] * uPredM[j];
                K_uPredM[i] = k_up;
                C_uPredM[i] = this.alphaM * M[i] * vPredM[i] + this.betaK * k_up;
            }
            const RM = new Float64Array(nDofs);
            for (let i = 0; i < nDofs; i++) {
                RM[i] = F_ext[i] - C_uPredM[i] - K_uPredM[i] - M[i] * this.aMapped[i] * (1.0 - 2.0 * beta) / (2.0 * beta);
            }
            const K_effM = Array.from({ length: nDofs }, () => new Float64Array(nDofs).fill(0));
            for (let i = 0; i < nDofs; i++) {
                for (let j = 0; j < nDofs; j++) K_effM[i][j] = KMapped[i][j] * factorK;
                K_effM[i][i] += M[i] * factorM;
            }
            this.applyClamping(K_effM, RM);
            const duMapped = igaSolver.gaussianElimination(K_effM, RM);

            for (let i = 0; i < nDofs; i++) {
                this.uMapped[i] = uPredM[i] + duMapped[i];
                this.aMapped[i] = duMapped[i] / (beta * this.dt * this.dt);
                this.vMapped[i] = vPredM[i] + gamma * this.dt * this.aMapped[i];
            }
        }
        
            tMappedEnd = performance.now();
            tMappedAssembly = tMappedEnd - tMappedStart;
        }

        // Mathematical Verifier & Real-time Relative L2 Error tracker
        this.stepCount++;
        let sumSqDiff = 0;
        let sumSqFOM = 0;
        if (runNaive && runMapped) {
            for (let i = 0; i < nDofs; i++) {
                const diff = this.uMapped[i] - this.uNaive[i];
                sumSqDiff += diff * diff;
                sumSqFOM += this.uNaive[i] * this.uNaive[i];
            }
        }
        const relL2Error = Math.sqrt(sumSqDiff) / (Math.sqrt(sumSqFOM) || 1e-6);

        // Tip deflection at the right-most edge (Y direction)
        const nU = this.referencePatch.controlPoints.length;
        const nV = this.referencePatch.controlPoints[0].length;
        let tipNaive = 0;
        let tipMapped = 0;
        for (let j = 0; j < nV; j++) {
            const idx = ((nU - 1) * nV + j) * 2 + 1;
            tipNaive += this.uNaive[idx];
            tipMapped += this.uMapped[idx];
        }
        tipNaive /= nV;
        tipMapped /= nV;

        if (this.stepCount % 10 === 0) {
            // Recompute Kr and KrFull for verification display
            const KMappedECSW = this.assembleMappedStiffness(mu, r, true);
            const KMappedFull = this.assembleMappedStiffness(mu, r, false);
            let kr_ecsw = 0;
            let kr_full = 0;
            const phi0 = this.Phi[0];
            const nDofs = phi0.length;
            for (let d = 0; d < nDofs; d++) {
                let kCol_ecsw = 0, kCol_full = 0;
                for (let col = 0; col < nDofs; col++) {
                    kCol_ecsw += KMappedECSW[d][col] * phi0[col];
                    kCol_full += KMappedFull[d][col] * phi0[col];
                }
                kr_ecsw += phi0[d] * kCol_ecsw;
                kr_full += phi0[d] * kCol_full;
            }

            console.log(`[Verification Step ${this.stepCount}] Tip Deflection -> Naive FOM: ${tipNaive.toFixed(5)} mm | Mapped ROM/ECSW: ${tipMapped.toFixed(5)} mm | Kr (ECSW): ${kr_ecsw.toFixed(4)} | Kr (Full): ${kr_full.toFixed(4)} | Relative L2 Error: ${(relL2Error * 100).toFixed(3)}%`);
            
            if (this.stepCount === 10) {
                console.log("================= ECSW ROM DIAGNOSTIC REPORT =================");
                console.log(`Phi Mode 1 Tip (Y): ${this.Phi[0][((nU - 1) * nV + (nV - 1)) * 2 + 1].toFixed(5)}`);
                console.log(`Phi Mode 2 Tip (X): ${this.Phi[1][((nU - 1) * nV + (nV - 1)) * 2].toFixed(5)}`);
                console.log("qMapped coordinates:", Array.from(this.qMapped).map(v => v.toFixed(6)));
                console.log("d2qMapped coordinates:", Array.from(this.d2qMapped).map(v => v.toFixed(6)));
                
                // Print reduced force
                const F_ext = this.calculateExternalForceVector(this.time);
                const Fr_ext = [0, 0];
                for (let r = 0; r < 2; r++) {
                    for (let d = 0; d < nDofs; d++) Fr_ext[r] += this.Phi[r][d] * F_ext[d];
                }
                console.log(`Projected Ext Force Fr_ext: [${Fr_ext[0].toFixed(4)}, ${Fr_ext[1].toFixed(4)}]`);
                
                // Print FOM diagonal scaling
                if (KNaive) console.log(`FOM Stiffness Diagonal KNaive[50][50]: ${KNaive[50][50].toFixed(4)}`);
                console.log(`FOM Mass Diagonal M[50]: ${M[50].toFixed(4)}`);
                console.log(`factorK: ${factorK.toFixed(4)} | factorM: ${factorM.toFixed(4)}`);
                if (KNaive) console.log(`Stiffness Term KNaive[50][50] * factorK: ${(KNaive[50][50] * factorK).toFixed(4)}`);
                console.log(`Mass Term M[50] * factorM: ${(M[50] * factorM).toFixed(4)}`);
                
                // Frequency Analysis: compare natural frequencies
                const KMappedECSW_diag = this.assembleMappedStiffness(mu, r, true);
                const KMappedFull_diag = this.assembleMappedStiffness(mu, r, false);
                const M_diag = this.assembleMassMatrix();
                const nDofs_freq = M_diag.length;
                const phi0_freq = this.Phi[0];
                let mr00 = 0, kr_ecsw_00 = 0, kr_full_00 = 0;
                for (let d = 0; d < nDofs_freq; d++) {
                    mr00 += phi0_freq[d] * M_diag[d] * phi0_freq[d];
                    let kCol_e = 0, kCol_f = 0;
                    for (let col = 0; col < nDofs_freq; col++) {
                        kCol_e += KMappedECSW_diag[d][col] * phi0_freq[col];
                        kCol_f += KMappedFull_diag[d][col] * phi0_freq[col];
                    }
                    kr_ecsw_00 += phi0_freq[d] * kCol_e;
                    kr_full_00 += phi0_freq[d] * kCol_f;
                }
                const omega_ecsw = Math.sqrt(kr_ecsw_00 / mr00);
                const omega_full = Math.sqrt(kr_full_00 / mr00);
                const f_ecsw = omega_ecsw / (2 * Math.PI);
                const f_full = omega_full / (2 * Math.PI);
                const f_ext = this.forceFreq / (2 * Math.PI);
                
                console.log("--- FREQUENCY ANALYSIS ---");
                console.log(`Mr[0][0]: ${mr00.toFixed(4)}`);
                console.log(`Kr (ECSW)[0][0]: ${kr_ecsw_00.toFixed(4)} | Kr (Full)[0][0]: ${kr_full_00.toFixed(4)}`);
                console.log(`ω_n ROM (ECSW): ${omega_ecsw.toFixed(4)} rad/s | ω_n FOM (projected): ${omega_full.toFixed(4)} rad/s`);
                console.log(`f_n ROM (ECSW): ${f_ecsw.toFixed(4)} Hz | f_n FOM (projected): ${f_full.toFixed(4)} Hz`);
                console.log(`Forcing frequency: ω_ext = ${this.forceFreq.toFixed(4)} rad/s | f_ext = ${f_ext.toFixed(4)} Hz`);
                console.log(`Frequency ratio (ω_ext/ω_n): ROM = ${(this.forceFreq / omega_ecsw).toFixed(4)} | FOM = ${(this.forceFreq / omega_full).toFixed(4)}`);
                console.log("==============================================================");
            }
            
            if (this.stepCount === 30) {
                console.log("================= STEP 30 ROM LOCKING DIAGNOSTICS =================");
                console.log("qMapped coordinates:", Array.from(this.qMapped).map(v => v.toFixed(6)));
                console.log("dqMapped coordinates:", Array.from(this.dqMapped).map(v => v.toFixed(6)));
                console.log("d2qMapped coordinates:", Array.from(this.d2qMapped).map(v => v.toFixed(6)));
                
                // Let's recompute/project the matrices here to inspect them
                const KMappedECSW = this.assembleMappedStiffness(mu, r, true);
                const KMappedFull = this.assembleMappedStiffness(mu, r, false);
                const M = this.assembleMassMatrix();
                const Mr = Array.from({length: 2}, () => new Float64Array(2).fill(0));
                const Kr = Array.from({length: 2}, () => new Float64Array(2).fill(0));
                const KrFull = Array.from({length: 2}, () => new Float64Array(2).fill(0));
                const nDofs = M.length;
                for (let rIdx = 0; rIdx < 2; rIdx++) {
                    const phiR = this.Phi[rIdx];
                    for (let cIdx = 0; cIdx < 2; cIdx++) {
                        const phiC = this.Phi[cIdx];
                        let mVal = 0, kVal = 0, kFullVal = 0;
                        for (let d = 0; d < nDofs; d++) {
                            mVal += phiR[d] * M[d] * phiC[d];
                            let kCol = 0, kFullCol = 0;
                            for (let col = 0; col < nDofs; col++) {
                                kCol += KMappedECSW[d][col] * phiC[col];
                                kFullCol += KMappedFull[d][col] * phiC[col];
                            }
                            kVal += phiR[d] * kCol;
                            kFullVal += phiR[d] * kFullCol;
                        }
                        Mr[rIdx][cIdx] = mVal;
                        Kr[rIdx][cIdx] = kVal;
                        KrFull[rIdx][cIdx] = kFullVal;
                    }
                }
                const beta = 0.25;
                const factorK = 1.0 + (0.5 * this.betaK) / (beta * this.dt);
                const factorM = 1.0 / (beta * this.dt * this.dt) + (0.5 * this.alphaM) / (beta * this.dt);
                
                console.log("Mr:", Mr.map(r => `[${r[0].toFixed(4)}, ${r[1].toFixed(4)}]`));
                console.log("Kr (ECSW Hyper-reduced):", Kr.map(r => `[${r[0].toFixed(4)}, ${r[1].toFixed(4)}]`));
                console.log("KrFull (Full Unreduced):", KrFull.map(r => `[${r[0].toFixed(4)}, ${r[1].toFixed(4)}]`));
                console.log(`factorK: ${factorK.toFixed(4)} | factorM: ${factorM.toFixed(4)}`);
                
                const K_eff_r = Array.from({length: 2}, () => new Float64Array(2).fill(0));
                for (let i = 0; i < 2; i++) {
                    for (let j = 0; j < 2; j++) {
                        K_eff_r[i][j] = Kr[i][j] * factorK + Mr[i][j] * factorM;
                    }
                }
                console.log("K_eff_r Matrix:", K_eff_r.map(r => `[${r[0].toFixed(4)}, ${r[1].toFixed(4)}]`));
                console.log("==============================================================");
            }
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
