/**
 * Phase 1.7+ Physics Engine (Refined for Geometric Accuracy)
 * Implements 1D IGA & FEM with variable Jacobian support.
 * Aligns with Theory 1.7: Stiffness ~ 1/J^3 and Load ~ R*J.
 * No external dependencies.
 */

export class PhysicsEngine {
    constructor(nurbs) {
        this.nurbs = nurbs; // Expects object with controlPoints, knots, degree
        this.E = 100.0;
        this.A = 1.0;
        this.I = 1.0;
        this.L = 1.0;
        this.physicsMode = 'bending';
    }

    get flexuralRigidity() { return this.E * this.I; }
    get axialRigidity() { return this.E * this.A; }

    /**
     * Calculates the Jacobian J = dx/dxi at a given parametric point.
     * For a 2D curve, J is the magnitude of the tangent vector (stretching factor).
     */
    calculateJacobian(xi) {
        // Evaluate first derivatives of all basis functions
        const d1 = this.nurbs.evaluateAllBasisDerivatives(xi, 1);
        const pts = this.nurbs.controlPoints;
        let dx = 0, dy = 0;

        // Compute tangent vector components
        for (let i = 0; i < pts.length; i++) {
            dx += d1[i] * pts[i].x;
            dy += d1[i] * pts[i].y;
        }

        // J is the Euclidean norm of the tangent
        const J = Math.sqrt(dx * dx + dy * dy);
        return Math.max(J, 1e-6); // Guard against division by zero
    }

    /**
     * Assembles the Global Stiffness Matrix using Gauss Quadrature.
     * Implements the 1/J^3 transformation for Euler-Bernoulli bending.
     */
    assembleStiffness() {
        const numCP = this.nurbs.controlPoints.length;
        const K = Array.from({ length: numCP }, () => new Float64Array(numCP));

        // Use high-density sampling for integration to capture Jacobian variations
        const samples = 150, dXi = 1.0 / samples;

        for (let s = 0; s <= samples; s++) {
            const xi = Math.min(s * dXi, 0.9999);
            const weight = (s === 0 || s === samples) ? 0.5 * dXi : dXi;
            const J = this.calculateJacobian(xi);

            if (this.physicsMode === 'bending') {
                // Curvature derivatives (Bi)
                const B = this.nurbs.evaluateAllBasisDerivatives(xi, 2);

                // Theory 1.7: Curvature transforms as (1/J^2). 
                // Differential dx transforms as (J * dXi).
                // Result: K = integral( EI * (B/J^2) * (B/J^2) * J ) dXi = integral( EI * B^2 / J^3 ) dXi
                const stiffnessFactor = this.flexuralRigidity / (J * J * J);

                for (let i = 0; i < numCP; i++) {
                    for (let j = 0; j < numCP; j++) {
                        K[i][j] += B[i] * B[j] * stiffnessFactor * weight;
                    }
                }
            } else {
                // Axial strain derivatives (Bi)
                const B = this.nurbs.evaluateAllBasisDerivatives(xi, 1);

                // Theory 1.7: Axial strain transforms as (1/J).
                // Result: K = integral( EA * (B/J) * (B/J) * J ) dXi = integral( EA * B^2 / J ) dXi
                const stiffnessFactor = this.axialRigidity / J;

                for (let i = 0; i < numCP; i++) {
                    for (let j = 0; j < numCP; j++) {
                        K[i][j] += B[i] * B[j] * stiffnessFactor * weight;
                    }
                }
            }
        }
        return K;
    }

    /**
     * Assembles the Load Vector (F) for a point load.
     * Implements the J-scaling required for consistent physical magnitude.
     */
    assembleIGALoad(loadPos, loadMag) {
        const numCP = this.nurbs.controlPoints.length;
        const F = new Float64Array(numCP).fill(0);

        // For axial mode, we usually apply at the end (xi=1.0)
        const actualPos = this.physicsMode === 'axial' ? 0.9999 : loadPos;

        const basis = this.nurbs.evaluateAllBasis(actualPos);
        const J = this.calculateJacobian(actualPos);

        // Theory 1.7: Point load sifting property F = P * R(xi) * J(xi)
        // Scaling by J ensures the force work is invariant under parametric refinement
        for (let i = 0; i < numCP; i++) {
            F[i] = basis[i] * loadMag * J;
        }
        return F;
    }

    /**
     * Constructs a modal transformation matrix (Phi).
     * Acts as the basis for the Reduced Order Model.
     */
    getModalBasis(modes, bcs) {
        const numCP = this.nurbs.controlPoints.length;
        const phi = Array.from({ length: numCP }, () => new Float64Array(modes));

        // Simple heuristic for boundary conditions in the reduced space
        const isRightFixed = bcs.some(bc => bc.index >= numCP - 2);

        for (let m = 0; m < modes; m++) {
            for (let i = 0; i < numCP; i++) {
                const xi = i / (numCP - 1);
                // Sinusoidal shapes are an excellent general basis for beam modal analysis
                phi[i][m] = isRightFixed ? Math.sin((m + 1) * Math.PI * xi) : Math.sin((m + 0.5) * Math.PI * xi);
            }
        }
        return phi;
    }

    /**
     * Reduced Order Model (ROM) Solver
     * Projects the large stiffness matrix into a small modal subspace.
     */
    solveROM(loadF, bcs, modes = 3) {
        const numCP = this.nurbs.controlPoints.length;
        const phi = this.getModalBasis(modes, bcs);
        const K = this.assembleStiffness();

        // 1. Reduced Stiffness: Kr = Phi^T * K * Phi (Galerkin Projection)
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

        // 2. Reduced Force: Fr = Phi^T * F
        const Fr = new Float64Array(modes);
        for (let i = 0; i < modes; i++) {
            let val = 0;
            for (let j = 0; j < numCP; j++) val += phi[j][i] * loadF[j];
            Fr[i] = val;
        }

        // 3. Solve for modal amplitudes (qr)
        const qr = this.gaussianElimination(Kr, Fr);

        // 4. Reconstruct physical displacement: u = Phi * qr
        const u = new Float64Array(numCP).fill(0);
        for (let i = 0; i < numCP; i++) {
            for (let j = 0; j < modes; j++) u[i] += phi[i][j] * qr[j];
        }
        return u;
    }

    /**
     * Full Static Solver
     * Applies boundary conditions using penalty-style enforcement.
     */
    solveStatics(loadF, bcs) {
        const numCP = this.nurbs.controlPoints.length;
        const K_orig = this.assembleStiffness();
        const K = Array.from({ length: numCP }, (_, i) => new Float64Array(K_orig[i]));
        const f = new Float64Array(loadF);

        // Apply boundary conditions
        bcs.forEach(bc => {
            const idx = bc.index;
            if (idx < 0 || idx >= numCP) return;
            for (let j = 0; j < numCP; j++) {
                if (j !== idx) {
                    f[j] -= K[j][idx] * bc.value;
                    K[idx][j] = 0; K[j][idx] = 0;
                }
            }
            K[idx][idx] = 1.0;
            f[idx] = bc.value;
        });

        return this.gaussianElimination(K, f);
    }

    /**
     * Standard Gaussian Elimination with pivoting for matrix solving.
     */
    gaussianElimination(A, b) {
        const n = b.length;
        for (let i = 0; i < n; i++) {
            let max = i;
            for (let j = i + 1; j < n; j++) if (Math.abs(A[j][i]) > Math.abs(A[max][i])) max = j;
            [A[i], A[max]] = [A[max], A[i]];[b[i], b[max]] = [b[max], b[i]];

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