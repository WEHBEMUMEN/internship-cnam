/**
 * NL App Physics Engine (Clean & Robust)
 * Implements 1D IGA & FEM for Bending and Axial physics.
 * No external dependencies (mathjs removed for compatibility).
 */

export class PhysicsEngine {
    constructor(nurbs) {
        this.nurbs = nurbs;
        this.E = 100.0; this.A = 1.0; this.I = 1.0; this.L = 1.0;
        this.physicsMode = 'bending';
    }

    get flexuralRigidity() { return this.E * this.I; }
    get axialRigidity() { return this.E * this.A; }

    solveTorqueToCircle(moment) {
        const kappa = moment / this.flexuralRigidity;
        const L = 1.0, numNodes = 100, result = [];
        for (let i = 0; i < numNodes; i++) {
            const s = (i / (numNodes - 1)) * L;
            if (Math.abs(kappa) < 1e-8) {
                result.push({ x: s, y: 0.5 });
            } else {
                const r = 1 / kappa, theta = s * kappa;
                result.push({ x: r * Math.sin(theta), y: 0.5 + r * (1 - Math.cos(theta)) });
            }
        }
        return result;
    }

    assembleStiffness() {
        const numCP = this.nurbs.controlPoints.length;
        const K = Array.from({ length: numCP }, () => new Float64Array(numCP));
        const samples = 100, dXi = 1.0 / samples;

        for (let s = 0; s <= samples; s++) {
            const xi = s * dXi;
            const w = (s === 0 || s === samples) ? 0.5 * dXi : dXi;
            if (this.physicsMode === 'bending') {
                const d2 = this.nurbs.evaluateAllBasisDerivatives(xi, 2);
                for (let i = 0; i < numCP; i++) for (let j = 0; j < numCP; j++) K[i][j] += d2[i] * d2[j] * this.flexuralRigidity * w;
            } else {
                const d1 = this.nurbs.evaluateAllBasisDerivatives(xi, 1);
                for (let i = 0; i < numCP; i++) for (let j = 0; j < numCP; j++) K[i][j] += d1[i] * d1[j] * this.axialRigidity * w;
            }
        }
        return K;
    }

    assembleIGALoad(loadPos, loadMag) {
        const numCP = this.nurbs.controlPoints.length;
        const F = new Float64Array(numCP).fill(0);
        const actualPos = this.physicsMode === 'axial' ? 1.0 : loadPos;
        const basis = this.nurbs.evaluateAllBasis(actualPos);
        for (let i = 0; i < numCP; i++) F[i] = basis[i] * loadMag;
        return F;
    }

    getModalBasis(modes, bcs, type = 'igen') {
        const numCP = this.nurbs.controlPoints.length;
        if (type === 'igen') return this.getTrueIGAModes(modes, bcs);

        // Improved Analytic Basis (Zero-Slope aware)
        const phi = Array.from({ length: numCP }, () => new Float64Array(modes));
        const isLeftFixed = bcs.some(bc => bc.index === 0);
        const isRightFixed = bcs.some(bc => bc.index >= numCP - 2);

        for (let m = 0; m < modes; m++) {
            for (let i = 0; i < numCP; i++) {
                const xi = i / (numCP - 1);
                if (isLeftFixed && !isRightFixed) {
                    // Cantilever: Clamped at 0, Free at 1. Use (1 - cos((m+0.5)PI*xi))
                    phi[i][m] = 1 - Math.cos((m + 0.5) * Math.PI * xi);
                } else if (isLeftFixed && isRightFixed) {
                    // Clamped-Clamped: Use (1 - cos(2*PI*(m+1)*xi))
                    phi[i][m] = 1 - Math.cos(2 * Math.PI * (m + 1) * xi);
                } else {
                    // Pinned-Pinned: Pure Sine
                    phi[i][m] = Math.sin((m + 1) * Math.PI * xi);
                }
            }
        }
        return phi;
    }

    getTrueIGAModes(modes, bcs) {
        const K_full = this.assembleStiffness();
        const numCP = K_full.length;
        
        // 1. Identify Free DOF house (Not in BCs)
        const fixedIndices = new Set(bcs.map(bc => bc.index));
        const freeIndices = [];
        for (let i = 0; i < numCP; i++) if (!fixedIndices.has(i)) freeIndices.push(i);
        const nFree = freeIndices.length;
        if (nFree === 0) return Array.from({ length: numCP }, () => new Float64Array(modes));

        // 2. Extract Constrained Submatrix A
        const A = Array.from({ length: nFree }, () => new Float64Array(nFree));
        for (let i = 0; i < nFree; i++) {
            for (let j = 0; j < nFree; j++) A[i][j] = K_full[freeIndices[i]][freeIndices[j]];
        }

        const phiFull = Array.from({ length: numCP }, () => new Float64Array(modes));
        
        // 3. Power Iteration + Hotelling Deflation for first M modes
        // Note: For stiffness K, the smallest eigenvalues are the lowest frequencies.
        // To use power iteration, we need to invert A or use the shifting/inverse iteration trick.
        // For simplicity and 1D performance, we'll invert A once (Gaussian).
        const I = Array.from({ length: nFree }, (_, i) => {
            const row = new Float64Array(nFree); row[i] = 1.0; return row;
        });
        const A_inv = I.map(row => this.gaussianElimination(Array.from(A, r => new Float64Array(r)), row));

        let Current_A = A_inv;
        for (let m = 0; m < modes; m++) {
            let v = new Float64Array(nFree).fill(1).map(() => Math.random());
            let lambda = 0;

            // Power iteration
            for (let iter = 0; iter < 50; iter++) {
                const next_v = new Float64Array(nFree);
                for (let i = 0; i < nFree; i++) {
                    for (let j = 0; j < nFree; j++) next_v[i] += Current_A[i][j] * v[j];
                }
                const norm = Math.sqrt(next_v.reduce((acc, val) => acc + val * val, 0));
                for (let i = 0; i < nFree; i++) v[i] = next_v[i] / norm;
                lambda = norm;
            }

            // Map back to global coordinates
            for (let i = 0; i < nFree; i++) phiFull[freeIndices[i]][m] = v[i];

            // Hotelling Deflation: A = A - lambda * v * vT
            for (let i = 0; i < nFree; i++) {
                for (let j = 0; j < nFree; j++) {
                    Current_A[i][j] -= lambda * v[i] * v[j];
                }
            }
        }
        
        return phiFull;
    }

    solveROM(loadF, bcs, modes = 3) {
        const numCP = this.nurbs.controlPoints.length;
        const phi = this.getModalBasis(modes, bcs);
        const K = this.assembleStiffness();
        
        // Reduced Stiffness: Kr = phi' * K * phi
        const Kr = Array.from({ length: modes }, () => new Float64Array(modes));
        for (let i = 0; i < modes; i++) {
            for (let j = 0; j < modes; j++) {
                let val = 0;
                for (let m = 0; m < numCP; m++) {
                    let temp = 0;
                    for (let l = 0; l < numCP; l++) temp += K[m][l] * phi[l][j];
                    val += phi[m][i] * temp;
                }
                Kr[i][j] = val;
            }
        }

        // Reduced Force: Fr = phi' * F
        const Fr = new Float64Array(modes);
        for (let i = 0; i < modes; i++) {
            let val = 0;
            for (let j = 0; j < numCP; j++) val += phi[j][i] * loadF[j];
            Fr[i] = val;
        }

        const qr = this.gaussianElimination(Kr, Fr);
        const u = new Float64Array(numCP).fill(0);
        for (let i = 0; i < numCP; i++) {
            for (let j = 0; j < modes; j++) u[i] += phi[i][j] * qr[j];
        }
        return u;
    }

    solveStatics(loadF, bcs) {
        const numCP = this.nurbs.controlPoints.length;
        const K_orig = this.assembleStiffness();
        const K = Array.from({ length: numCP }, (_, i) => new Float64Array(K_orig[i]));
        const f = new Float64Array(loadF);

        bcs.forEach(bc => {
            const idx = bc.index;
            if (idx < 0 || idx >= numCP) return;
            for (let j = 0; j < numCP; j++) {
                if (j !== idx) {
                    f[j] -= K[j][idx] * bc.value;
                    K[idx][j] = 0; K[j][idx] = 0;
                }
            }
            K[idx][idx] = 1.0; f[idx] = bc.value;
        });
        return this.gaussianElimination(K, f);
    }

    assembleGeometricStiffness(u) {
        const numCP = this.nurbs.controlPoints.length;
        const Kg = Array.from({ length: numCP }, () => new Float64Array(numCP));
        const samples = 100, dXi = 1.0 / samples;
        const EA = this.axialRigidity;

        // 1. Calculate average axial strain from large deflections: 1/2 * integral( (v')^2 ) dx
        let stretchEnergy = 0;
        for (let s = 0; s <= samples; s++) {
            const xi = s * dXi;
            const w = (s === 0 || s === samples) ? 0.5 * dXi : dXi;
            const d1 = this.nurbs.evaluateAllBasisDerivatives(xi, 1);
            let slope = 0;
            for(let i=0; i<numCP; i++) slope += d1[i] * u[i];
            stretchEnergy += slope * slope * w;
        }

        // Induced axial tension (stress stiffening)
        const N_ax = (EA / (2 * this.L)) * stretchEnergy;

        // 2. Assemble Geometric Stiffness Matrix (Kg)
        for (let s = 0; s <= samples; s++) {
            const xi = s * dXi;
            const w = (s === 0 || s === samples) ? 0.5 * dXi : dXi;
            const d1 = this.nurbs.evaluateAllBasisDerivatives(xi, 1);
            for (let i = 0; i < numCP; i++) {
                for (let j = 0; j < numCP; j++) {
                    Kg[i][j] += d1[i] * d1[j] * N_ax * w;
                }
            }
        }
        return Kg;
    }

    solveNonLinearStatics(loadF, bcs, iterations = 20, method = 'newton', initialU = null) {
        const numCP = this.nurbs.controlPoints.length;
        let u = initialU ? new Float64Array(initialU) : new Float64Array(numCP).fill(0);
        let iterCount = 0;
        let residualHistory = [];
        
        const K_lin = this.assembleStiffness();
        let K_tangent = null;
        if (method === 'modified_newton' || method === 'incremental') {
            K_tangent = Array.from({ length: numCP }, (_, i) => new Float64Array(numCP));
            const K_geom_init = this.assembleGeometricStiffness(u);
            for (let i = 0; i < numCP; i++) {
                for (let j = 0; j < numCP; j++) K_tangent[i][j] = K_lin[i][j] + K_geom_init[i][j];
            }
        }

        for (let iter = 0; iter < iterations; iter++) {
            iterCount++;
            
            let K_curr;
            if (method === 'newton' || (method === 'picard' && iter % 2 === 0)) {
                const K_geom = this.assembleGeometricStiffness(u);
                K_curr = Array.from({ length: numCP }, (_, i) => new Float64Array(numCP));
                for (let i = 0; i < numCP; i++) {
                    for (let j = 0; j < numCP; j++) K_curr[i][j] = K_lin[i][j] + K_geom[i][j];
                }
            } else if (method === 'modified_newton' || method === 'incremental') {
                K_curr = Array.from({ length: numCP }, (_, i) => new Float64Array(K_tangent[i]));
            } else {
                K_curr = Array.from({ length: numCP }, (_, i) => new Float64Array(K_lin[i]));
            }

            const K_geom_curr = this.assembleGeometricStiffness(u);
            const F_int = new Float64Array(numCP).fill(0);
            for (let i = 0; i < numCP; i++) {
                for (let j = 0; j < numCP; j++) {
                    F_int[i] += (K_lin[i][j] + 0.5 * K_geom_curr[i][j]) * u[j];
                }
            }

            const R = new Float64Array(numCP);
            for (let i = 0; i < numCP; i++) R[i] = loadF[i] - F_int[i];

            let resNorm = 0;
            for (let i = 0; i < numCP; i++) resNorm += R[i] * R[i];
            const currentRes = Math.sqrt(resNorm);
            residualHistory.push(currentRes);

            if (currentRes < 1e-7 && iter > 0) break;
            if (method === 'incremental' && iter > 0) break;

            bcs.forEach(bc => {
                const idx = bc.index;
                if (idx < 0 || idx >= numCP) return;
                for (let j = 0; j < numCP; j++) {
                    if (j !== idx) {
                        R[j] -= K_curr[j][idx] * 0;
                        K_curr[idx][j] = 0; K_curr[j][idx] = 0;
                    }
                }
                K_curr[idx][idx] = 1.0; R[idx] = 0;
            });

            const du = this.gaussianElimination(K_curr, R);
            let updateNorm = 0;
            for (let i = 0; i < numCP; i++) {
                u[i] += du[i];
                updateNorm += du[i] * du[i];
            }
            if (Math.sqrt(updateNorm) < 1e-9) break;
        }

        return { u, iterations: iterCount, residualHistory };
    }

    solveNonLinearROM(loadF, bcs, modes = 3, iterations = 20, method = 'newton', initialQ = null, phi = null) {
        const numCP = this.nurbs.controlPoints.length;
        phi = phi || this.getModalBasis(modes, bcs);
        let q = initialQ ? new Float64Array(initialQ) : new Float64Array(modes).fill(0);
        let iterCount = 0;
        let residualHistory = [];

        const K_lin_full = this.assembleStiffness();
        let Kr_tangent = null;
        if (method === 'modified_newton' || method === 'incremental') {
            const K_geom_init = this.assembleGeometricStiffness(new Float64Array(numCP).fill(0));
            const K_T_init = Array.from({ length: numCP }, (_, i) => new Float64Array(numCP));
            for (let i = 0; i < numCP; i++) {
                for (let j = 0; j < numCP; j++) K_T_init[i][j] = K_lin_full[i][j] + K_geom_init[i][j];
            }
            Kr_tangent = Array.from({ length: modes }, () => new Float64Array(modes));
            for (let i = 0; i < modes; i++) {
                for (let j = 0; j < modes; j++) {
                    let val = 0;
                    for (let m = 0; m < numCP; m++) {
                        let temp = 0;
                        for (let l = 0; l < numCP; l++) temp += K_T_init[m][l] * phi[l][j];
                        val += phi[m][i] * temp;
                    }
                    Kr_tangent[i][j] = val;
                }
            }
        }

        for (let iter = 0; iter < iterations; iter++) {
            iterCount++;
            
            const u = new Float64Array(numCP).fill(0);
            for (let i = 0; i < numCP; i++) {
                for (let j = 0; j < modes; j++) u[i] += phi[i][j] * q[j];
            }

            const K_geom = this.assembleGeometricStiffness(u);
            const F_int = new Float64Array(numCP).fill(0);
            const R_full = new Float64Array(numCP);

            for (let i = 0; i < numCP; i++) {
                for (let j = 0; j < numCP; j++) {
                    F_int[i] += (K_lin_full[i][j] + 0.5 * K_geom[i][j]) * u[j];
                }
                R_full[i] = loadF[i] - F_int[i];
            }

            let Kr_curr;
            if (method === 'newton' || (method === 'picard' && iter % 2 === 0)) {
                const K_T = Array.from({ length: numCP }, () => new Float64Array(numCP));
                for (let i = 0; i < numCP; i++) {
                    for (let j = 0; j < numCP; j++) K_T[i][j] = K_lin_full[i][j] + K_geom[i][j];
                }
                Kr_curr = Array.from({ length: modes }, () => new Float64Array(modes));
                for (let i = 0; i < modes; i++) {
                    for (let j = 0; j < modes; j++) {
                        let val = 0;
                        for (let m = 0; m < numCP; m++) {
                            let temp = 0;
                            for (let l = 0; l < numCP; l++) temp += K_T[m][l] * phi[l][j];
                            val += phi[m][i] * temp;
                        }
                        Kr_curr[i][j] = val;
                    }
                }
            } else if (method === 'modified_newton' || method === 'incremental') {
                Kr_curr = Array.from({ length: modes }, (_, i) => new Float64Array(Kr_tangent[i]));
            } else {
                Kr_curr = Array.from({ length: modes }, () => new Float64Array(modes));
                for (let i = 0; i < modes; i++) {
                    for (let j = 0; j < modes; j++) {
                        let val = 0;
                        for (let m = 0; m < numCP; m++) {
                            let temp = 0;
                            for (let l = 0; l < numCP; l++) temp += K_lin_full[m][l] * phi[l][j];
                            val += phi[m][i] * temp;
                        }
                        Kr_curr[i][j] = val;
                    }
                }
            }

            const R_red = new Float64Array(modes).fill(0);
            for (let i = 0; i < modes; i++) {
                for (let j = 0; j < numCP; j++) R_red[i] += phi[j][i] * R_full[j];
            }

            let resNorm = 0;
            for (let i = 0; i < modes; i++) resNorm += R_red[i] * R_red[i];
            const currentRes = Math.sqrt(resNorm);
            residualHistory.push(currentRes);

            if (currentRes < 1e-7 && iter > 0) break;
            if (method === 'incremental' && iter > 0) break;

            const dq = this.gaussianElimination(Kr_curr, R_red);
            let updateNorm = 0;
            for (let i = 0; i < modes; i++) {
                q[i] += dq[i];
                updateNorm += dq[i] * dq[i];
            }
            if (Math.sqrt(updateNorm) < 1e-9) break;
        }

        const u_final = new Float64Array(numCP).fill(0);
        for (let i = 0; i < numCP; i++) {
            for (let j = 0; j < modes; j++) u_final[i] += phi[i][j] * q[j];
        }
        return { u: u_final, iterations: iterCount, residualHistory, q };
    }


    calculatePhysicsState(u) {
        const samples = 100, results = { stresses: [], maxStress: 0 };
        const E = this.E, I = this.I, A = this.A;
        const numCP = this.nurbs.controlPoints.length;

        for (let i = 0; i <= samples; i++) {
            const xi = i / samples;
            let stress = 0;
            
            if (this.physicsMode === 'bending') {
                const d2 = this.nurbs.evaluateAllBasisDerivatives(xi, 2);
                let curvature = 0;
                for (let j = 0 ; j < numCP; j++) curvature += d2[j] * u[j];
                // Fiber stress at top: sigma = E * y * kappa (y = half beam depth, assuming 0.1 for visual)
                stress = E * 0.1 * curvature;
            } else {
                const d1 = this.nurbs.evaluateAllBasisDerivatives(xi, 1);
                let strain = 0;
                for (let j = 0; j < numCP; j++) strain += d1[j] * u[j];
                stress = E * strain;
            }
            results.stresses.push(stress);
            if (Math.abs(stress) > results.maxStress) results.maxStress = Math.abs(stress);
        }
        return results;
    }

    solveIncremental(loadPos, totalMag, bcs, steps = 10, method = 'newton') {
        const results = { path: [] };
        const numCP = this.nurbs.controlPoints.length;
        let currentU = new Float64Array(numCP).fill(0);
        
        for (let s = 1; s <= steps; s++) {
            const stepMag = (s / steps) * totalMag;
            const stepLoadF = this.assembleIGALoad(loadPos, stepMag);
            const solveRes = this.solveNonLinearStatics(stepLoadF, bcs, 20, method, currentU);
            currentU = solveRes.u;
            
            // Track tip displacement (assuming right end is DOFs near numCP-1)
            const tipDisp = currentU[numCP - 1]; 
            results.path.push({ loadFactor: s / steps, displacement: tipDisp });
        }
        return results;
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

    solveTorqueToCircle(moment) {
        const kappa = moment / this.flexuralRigidity;
        const L = 1.0, numNodes = 100, result = [];
        for (let i = 0; i < numNodes; i++) {
            const s = (i / (numNodes - 1)) * L;
            if (Math.abs(kappa) < 1e-8) {
                result.push({ x: s, y: 0.5 });
            } else {
                const r = 1 / kappa, theta = s * kappa;
                result.push({ x: r * Math.sin(theta), y: 0.5 + r * (1 - Math.cos(theta)) });
            }
        }
        return result;
    }

    runBenchmark(loadMag, numRuns = 20) {
        const numCP = this.nurbs.controlPoints.length;
        const bcs = [{ index: 0, value: 0 }, { index: 1, value: 0 }];
        const loadF = this.assembleIGALoad(0.5, loadMag);
        const t0 = performance.now();
        for (let i = 0; i < numRuns; i++) this.solveStatics(loadF, bcs);
        const t1 = performance.now();
        return { iga: ((t1 - t0) / numRuns).toFixed(3), dofs: numCP, mode: this.physicsMode };
    }

    getDegreesOfFreedom() {
        return this.nurbs.controlPoints.length;
    }
}

export class ReferenceFEM {
    constructor(numElements = 100) {
        this.numElements = numElements;
        this.E = 100.0; this.A = 1.0; this.I = 1.0;
    }

    solve(loadPos, loadMag, igaBCs, mode = 'bending', numCP = 0, p = 1) {
        const n = this.numElements, L = 1.0, h = L / n;
        
        // Force Bending to use Hermite (C1) to prevent singular stiffness matrices.
        // We will visualize p=1,2 by linearizing/quadratizing the interpolation later.
        const isHermite = (mode === 'bending');
        const numNodes = isHermite ? (n + 1) : (n * p + 1);
        const dofsPerNode = isHermite ? 2 : 1;
        const numDofs = numNodes * dofsPerNode;

        const K_orig = this.assembleFEMStiffness(h, mode, p);
        const F_orig = this.assembleFEMLoad(loadPos, loadMag, mode, p);
        
        if (!K_orig || K_orig.length === 0) return null;

        const K = Array.from({ length: numDofs }, (_, i) => new Float64Array(K_orig[i]));
        const F = new Float64Array(F_orig);

        const femBCs = new Set();
        if (igaBCs.some(bc => bc.index === 0)) femBCs.add(0); // Left displacement
        if (igaBCs.some(bc => bc.index === 1) && isHermite) femBCs.add(1); // Left rotation
        
        if (numCP > 0) {
            if (igaBCs.some(bc => bc.index === numCP - 1)) femBCs.add((numNodes - 1) * dofsPerNode); // Right displacement
            if (igaBCs.some(bc => bc.index === numCP - 2) && isHermite) femBCs.add((numNodes - 1) * dofsPerNode + 1); // Right rotation
        }

        [...femBCs].forEach(idx => {
            if (idx < 0 || idx >= numDofs) return;
            for (let j = 0; j < numDofs; j++) {
                if (j !== idx) {
                    F[j] -= K[j][idx] * 0;
                    K[idx][j] = 0; K[j][idx] = 0;
                }
            }
            K[idx][idx] = 1.0; F[idx] = 0;
        });

        const u = this.gaussianElimination(K, F);
        
        const result = [];
        const plotSamples = 100;
        for (let i = 0; i <= plotSamples; i++) {
            const xi = i / plotSamples;
            const x = xi * L;
            const element = Math.min(Math.floor(xi * n), n - 1);
            const localXi = (xi * n) - element;
            
            let y = 0;
            if (isHermite) {
                const baseDof = element * 2;
                const nextDof = (element + 1) * 2;
                if (p === 1) {
                    y = u[baseDof] * (1 - localXi) + u[nextDof] * localXi; // Linear visual
                } else if (p === 2) {
                    // Quadratic visual using midpoint
                    const midY = (u[baseDof] + u[nextDof]) / 2 + (u[baseDof+1] - u[nextDof+1]) * h / 8;
                    y = u[baseDof] * (1 - localXi) * (1 - 2 * localXi) +
                        midY * 4 * localXi * (1 - localXi) +
                        u[nextDof] * localXi * (2 * localXi - 1);
                } else {
                    const h1 = 1 - 3*localXi**2 + 2*localXi**3;
                    const h2 = h * (localXi - 2*localXi**2 + localXi**3);
                    const h3 = 3*localXi**2 - 2*localXi**3;
                    const h4 = h * (localXi**3 - localXi**2);
                    y = u[baseDof]*h1 + u[baseDof+1]*h2 + u[nextDof]*h3 + u[nextDof+1]*h4;
                }
            } else {
                const baseNode = element * p;
                if (p === 1) {
                    y = u[baseNode] * (1 - localXi) + u[baseNode + 1] * localXi;
                } else if (p === 2) {
                    y = u[baseNode] * (1 - localXi) * (1 - 2 * localXi) +
                        u[baseNode + 1] * 4 * localXi * (1 - localXi) +
                        u[baseNode + 2] * localXi * (2 * localXi - 1);
                } else if (p === 3) {
                    // Cubic Lagrangian
                    const n1 = -(9/2)*(localXi - 1/3)*(localXi - 2/3)*(localXi - 1);
                    const n2 = (27/2)*localXi*(localXi - 2/3)*(localXi - 1);
                    const n3 = -(27/2)*localXi*(localXi - 1/3)*(localXi - 1);
                    const n4 = (9/2)*localXi*(localXi - 1/3)*(localXi - 2/3);
                    y = u[baseNode]*n1 + u[baseNode+1]*n2 + u[baseNode+2]*n3 + u[baseNode+3]*n4;
                }
            }
            result.push({ x, y });
        }
        return result;
    }

    assembleFEMStiffness(h, mode, p) {
        const n = this.numElements;
        const isHermite = (mode === 'bending');
        const numNodes = isHermite ? (n + 1) : (n * p + 1);
        const dofsPerNode = isHermite ? 2 : 1;
        const numDofs = numNodes * dofsPerNode;
        const K = Array.from({ length: numDofs }, () => new Float64Array(numDofs));

        if (mode === 'axial') {
            const EA = this.E * this.A;
            let ke;
            if (p === 1) {
                const k = EA / h;
                ke = [[k, -k], [-k, k]];
            } else if (p === 2) {
                const k = EA / (3 * h);
                ke = [[7*k, -8*k, k], [-8*k, 16*k, -8*k], [k, -8*k, 7*k]];
            } else {
                const k = EA / (40 * h);
                ke = [[148*k, -189*k, 54*k, -13*k], [-189*k, 432*k, -297*k, 54*k], [54*k, -297*k, 432*k, -189*k], [-13*k, 54*k, -189*k, 148*k]];
            }
            for (let e = 0; e < n; e++) {
                const start = e * p;
                for (let i = 0; i <= p; i++) for (let j = 0; j <= p; j++) K[start + i][start + j] += ke[i][j];
            }
            return K;
        } else {
            const EI = this.E * this.I;
            const k_rel = EI / (h**3);
            const ke = [[12,6*h,-12,6*h], [6*h,4*h*h,-6*h,2*h*h], [-12,-6*h,12,-6*h], [6*h,2*h*h,-6*h,4*h*h]];
            for (let e = 0; e < n; e++) {
                const d = [e*2, e*2+1, (e+1)*2, (e+1)*2+1];
                for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) K[d[i]][d[j]] += ke[i][j] * k_rel;
            }
            return K;
        }
    }

    assembleFEMLoad(loadPos, loadMag, mode, p) {
        const n = this.numElements;
        const isHermite = (mode === 'bending');
        const numNodes = isHermite ? (n + 1) : (n * p + 1);
        const dofsPerNode = isHermite ? 2 : 1;
        const numDofs = numNodes * dofsPerNode;
        const F = new Float64Array(numDofs).fill(0);

        if (isHermite) {
            const loadNode = Math.round(loadPos * n);
            const loadDof = loadNode * 2;
            if (loadDof < F.length) F[loadDof] = loadMag;
        } else {
            const loadNode = Math.round(loadPos * n * p);
            if (loadNode < F.length) F[loadNode] = loadMag;
        }
        return F;
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

    getDegreesOfFreedom() { 
        return this.numElements * 2 + 2; // Approximation for beam
    }
}
