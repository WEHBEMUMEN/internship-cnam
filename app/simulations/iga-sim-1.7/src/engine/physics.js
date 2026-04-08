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
     * Assembles the Global Stiffness Matrix using Gauss-Legendre Quadrature.
     * Evaluates integrals exactly knot-span by knot-span.
     * Implements the 1/J^3 transformation for Euler-Bernoulli bending.
     */
    assembleStiffness() {
        const numCP = this.nurbs.controlPoints.length;
        const K = Array.from({ length: numCP }, () => new Float64Array(numCP));

        // Gauss-Legendre points and weights for [-1, 1]
        const gaussData = {
            1: { pts: [0], w: [2] },
            2: { pts: [-0.577350269, 0.577350269], w: [1, 1] },
            3: { pts: [-0.774596669, 0, 0.774596669], w: [0.555555556, 0.888888889, 0.555555556] },
            4: { pts: [-0.861136312, -0.339981044, 0.339981044, 0.861136312], w: [0.347854845, 0.652145155, 0.652145155, 0.347854845] },
            5: { pts: [-0.906179846, -0.538469310, 0, 0.538469310, 0.906179846], w: [0.236926885, 0.478628670, 0.568888889, 0.478628670, 0.236926885] },
            6: { pts: [-0.932469514,-0.661209386,-0.238619186,0.238619186,0.661209386,0.932469514], w: [0.171324492,0.360761573,0.467913935,0.467913935,0.360761573,0.171324492]}
        };

        const p = this.nurbs.degree;
        const knots = this.nurbs.knots;
        // Exact integration requires p+1 points
        const ng = Math.min(p + 1, 6); 
        const gPts = gaussData[ng].pts;
        const gW = gaussData[ng].w;

        // Iterate over non-empty knot spans
        for (let span = 0; span < knots.length - 1; span++) {
            const xiK = knots[span];
            const xiK1 = knots[span + 1];
            if (xiK1 - xiK < 1e-10) continue; // Skip zero-length knot spans

            // Map standard [-1, 1] to arbitrary [xiK, xiK1]
            const mapScale = (xiK1 - xiK) / 2.0;
            const mapOffset = (xiK1 + xiK) / 2.0;

            for (let g = 0; g < ng; g++) {
                const xi = gPts[g] * mapScale + mapOffset;
                const weight = gW[g] * mapScale; // scaled by domain mapping derivative
                const J = this.calculateJacobian(xi);

                if (this.physicsMode === 'bending') {
                    const B = this.nurbs.evaluateAllBasisDerivatives(xi, 2);
                    const stiffnessFactor = this.flexuralRigidity / (J * J * J);

                    for (let i = 0; i < numCP; i++) {
                        for (let j = 0; j < numCP; j++) {
                            K[i][j] += B[i] * B[j] * stiffnessFactor * weight;
                        }
                    }
                } else {
                    // Axial mode
                    const B = this.nurbs.evaluateAllBasisDerivatives(xi, 1);
                    const stiffnessFactor = this.axialRigidity / J;

                    for (let i = 0; i < numCP; i++) {
                        for (let j = 0; j < numCP; j++) {
                            K[i][j] += B[i] * B[j] * stiffnessFactor * weight;
                        }
                    }
                }
            }
        }
        return K;
    }

    /**
     * Assembles the Consistent Mass Matrix using Gauss-Legendre Quadrature.
     * M_ij = integral( rho * A * N_i * N_j * J ) dXi
     */
    assembleMass() {
        const numCP = this.nurbs.controlPoints.length;
        const M = Array.from({ length: numCP }, () => new Float64Array(numCP));

        // Using same Gauss data as stiffness
        const gaussData = {
            1: { pts: [0], w: [2] },
            2: { pts: [-0.577350269, 0.577350269], w: [1, 1] },
            3: { pts: [-0.774596669, 0, 0.774596669], w: [0.555555556, 0.888888889, 0.555555556] },
            4: { pts: [-0.861136312, -0.339981044, 0.339981044, 0.861136312], w: [0.347854845, 0.652145155, 0.652145155, 0.347854845] },
            5: { pts: [-0.906179846, -0.538469310, 0, 0.538469310, 0.906179846], w: [0.236926885, 0.478628670, 0.568888889, 0.478628670, 0.236926885] },
            6: { pts: [-0.932469514,-0.661209386,-0.238619186,0.238619186,0.661209386,0.932469514], w: [0.171324492,0.360761573,0.467913935,0.467913935,0.360761573,0.171324492]}
        };

        const p = this.nurbs.degree;
        const knots = this.nurbs.knots;
        const ng = Math.min(p + 1, 6);
        const gPts = gaussData[ng].pts;
        const gW = gaussData[ng].w;
        const rhoA = 50.0; // Increased artificial mass for slower, visual browser oscillations

        for (let span = 0; span < knots.length - 1; span++) {
            const xiK = knots[span];
            const xiK1 = knots[span + 1];
            if (xiK1 - xiK < 1e-10) continue;

            const mapScale = (xiK1 - xiK) / 2.0;
            const mapOffset = (xiK1 + xiK) / 2.0;

            for (let g = 0; g < ng; g++) {
                const xi = gPts[g] * mapScale + mapOffset;
                const weight = gW[g] * mapScale;
                const J = this.calculateJacobian(xi);
                const R = this.nurbs.evaluateAllBasis(xi);

                for (let i = 0; i < numCP; i++) {
                    for (let j = 0; j < numCP; j++) {
                        M[i][j] += R[i] * R[j] * rhoA * J * weight;
                    }
                }
            }
        }
        return M;
    }

    resetDynamics() {
        const numCP = this.nurbs.controlPoints.length;
        this.u_t = new Float64Array(numCP).fill(0);
        this.v_t = new Float64Array(numCP).fill(0);
        this.a_t = new Float64Array(numCP).fill(0);
        this.K_dyn = null;
        this.M_dyn = null;
        this.dynBCsLength = -1;
    }

    /**
     * Solves structural dynamics via Newmark-beta method.
     * Integrates generalized equation of motion: M*a + C*v + K*u = F
     */
    solveDynamicsNewmark(loadF, bcs, dt, useDamping) {
        const numCP = this.nurbs.controlPoints.length;

        // Init states if missing or changed topology
        if (!this.u_t || this.u_t.length !== numCP || this.dynBCsLength !== bcs.length || !this.K_dyn) {
            this.resetDynamics();
            this.K_dyn = this.assembleStiffness();
            this.M_dyn = this.assembleMass();
            this.dynBCsLength = bcs.length;

            // Apply penalty BCs to K and M
            bcs.forEach(bc => {
                const idx = bc.index;
                if (idx < 0 || idx >= numCP) return;
                for (let j = 0; j < numCP; j++) {
                    this.K_dyn[idx][j] = 0; this.K_dyn[j][idx] = 0;
                    this.M_dyn[idx][j] = 0; this.M_dyn[j][idx] = 0;
                }
                this.K_dyn[idx][idx] = 1.0;
                this.M_dyn[idx][idx] = 1.0; 
            });
        }

        // Apply BCs to the current load vector
        const F_t = new Float64Array(loadF);
        bcs.forEach(bc => {
            const idx = bc.index;
            if (idx >= 0 && idx < numCP) F_t[idx] = bc.value; 
        });

        // Newmark params (Average Acceleration Method)
        const beta = 0.25;
        const gamma = 0.5;

        // Rayleigh Damping matrix C = alpha*M + betaR*K
        const alpha = useDamping ? 0.2 : 0;
        const betaR = useDamping ? 0.002 : 0;

        const a0 = 1.0 / (beta * dt * dt);
        const a1 = gamma / (beta * dt);
        const a2 = 1.0 / (beta * dt);
        const a3 = 1.0 / (2 * beta) - 1.0;
        const a4 = gamma / beta - 1.0;
        const a5 = (dt / 2.0) * (gamma / beta - 2.0);

        const K_hat = Array.from({ length: numCP }, () => new Float64Array(numCP));
        const F_hat = new Float64Array(numCP);

        for (let i = 0; i < numCP; i++) {
            let sumM_u = 0, sumM_v = 0, sumM_a = 0;
            let sumC_u = 0, sumC_v = 0, sumC_a = 0;

            for (let j = 0; j < numCP; j++) {
                const M_ij = this.M_dyn[i][j];
                const K_ij = this.K_dyn[i][j];
                const C_ij = alpha * M_ij + betaR * K_ij;

                K_hat[i][j] = K_ij + a0 * M_ij + a1 * C_ij;

                sumM_u += M_ij * this.u_t[j];
                sumM_v += M_ij * this.v_t[j];
                sumM_a += M_ij * this.a_t[j];

                sumC_u += C_ij * this.u_t[j];
                sumC_v += C_ij * this.v_t[j];
                sumC_a += C_ij * this.a_t[j];
            }

            F_hat[i] = F_t[i] 
                + M_ij_term(a0, a2, a3, sumM_u, sumM_v, sumM_a) 
                + C_ij_term(a1, a4, a5, sumC_u, sumC_v, sumC_a);
        }

        // Helper closures for cleaner addition
        function M_ij_term(a_0, a_2, a_3, s_u, s_v, s_a) { return a_0 * s_u + a_2 * s_v + a_3 * s_a; }
        function C_ij_term(a_1, a_4, a_5, s_u, s_v, s_a) { return a_1 * s_u + a_4 * s_v + a_5 * s_a; }

        bcs.forEach(bc => {
            const idx = bc.index;
            if (idx < 0 || idx >= numCP) return;
            for (let j = 0; j < numCP; j++) {
                if (j !== idx) {
                    F_hat[j] -= K_hat[j][idx] * bc.value;
                    K_hat[idx][j] = 0; K_hat[j][idx] = 0;
                }
            }
            K_hat[idx][idx] = 1.0;
            F_hat[idx] = bc.value;
        });

        const u_next = this.gaussianElimination(K_hat, F_hat);

        for (let i = 0; i < numCP; i++) {
            const a_next = a0 * (u_next[i] - this.u_t[i]) - a2 * this.v_t[i] - a3 * this.a_t[i];
            const v_next = this.v_t[i] + dt * ((1 - gamma) * this.a_t[i] + gamma * a_next);

            this.a_t[i] = a_next;
            this.v_t[i] = v_next;
            this.u_t[i] = u_next[i];
        }

        return new Float64Array(u_next);
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

        // 2. Reduced Force: Fr = Phi^T * F (Hyper-reduction Selective Extraction)
        // Extract only the rows of Phi corresponding to the local support of the load
        const Fr = new Float64Array(modes);
        
        // Find active degrees of freedom (non-zero load entries)
        const activeIndices = [];
        for (let j = 0; j < numCP; j++) {
            if (Math.abs(loadF[j]) > 1e-10) {
                activeIndices.push(j);
            }
        }

        // Fast extraction inner product 
        for (let i = 0; i < modes; i++) {
            let val = 0;
            for (let k = 0; k < activeIndices.length; k++) {
                const j = activeIndices[k];
                val += phi[j][i] * loadF[j];
            }
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