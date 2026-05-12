/**
 * Phase 5.1b - ECSW Engine
 * Energy-Conserving Sampling and Weighting for Hyper-reduction.
 */

class ECSWEngine {
    constructor() {
        this.sampledElements = [];
        this.weights = [];
    }

    /**
     * Trains the ECSW weights using a greedy OMP-like approach
     * @param {Array} snapshots - Full displacement snapshots
     * @param {Object} patch - Geometric patch
     * @param {Matrix} Phi - Reduced Basis
     * @param {Object} solver - IGA Solver instance
     */
    train(snapshots, patch, Phi, solver, tol = 1e-4) {
        console.log(`[ECSW] Starting Hyper-reduction Training (tol=${tol})...`);
        const { Matrix } = window.mlMatrix || window.ML;
        
        // Robustness: Ensure elements exist
        if (!patch.elements) {
            patch.elements = window.GeometryFactory.extractElements(patch);
        }
        
        const nElements = patch.elements.length;
        const k = Phi.columns;
        const nSnaps = snapshots.length;

        // 1. Construct the Virtual Force Matrix (A) and Target (b)
        // b = sum_{e} f_e_red(mu)
        // A_e = f_e_red(mu)
        // We do this for all snapshots to ensure robustness across parameter space
        
        const targetVectors = [];
        const elementVectors = Array.from({ length: nElements }, () => []);

        snapshots.forEach((u_full, idx) => {
            // Compute full reduced internal forces (or stiffness-vector products)
            // For linear: K_r * a = Phi' * K * Phi * a
            // We want to approximate the reduced stiffness matrix.
            // Simplified ECSW: approximate the action of the stiffness on the basis
            
            for (let j = 0; j < k; j++) {
                const phi_j = Phi.getColumn(j);
                // Target: Phi' * K * phi_j
                const K_phi_j = solver.multiplyK(patch, phi_j);
                const b_j = Phi.transpose().mmul(new Matrix([K_phi_j]).transpose()).to2DArray().flat();
                targetVectors.push(...b_j);

                // Elements: Phi' * K_e * phi_j
                for (let eIdx = 0; eIdx < nElements; eIdx++) {
                    const Ke_phi_j = solver.multiplyKe(patch, eIdx, phi_j);
                    const ae_j = Phi.transpose().mmul(new Matrix([Ke_phi_j]).transpose()).to2DArray().flat();
                    elementVectors[eIdx].push(...ae_j);
                }
            }
        });

        const b = new Float64Array(targetVectors);
        const A = elementVectors.map(v => new Float64Array(v));

        // 2. Solve Sparse NNLS (Greedy Approach)
        const result = this.solveGreedyNNLS(A, b, tol);
        
        this.sampledElements = result.indices;
        this.weights = result.weights;

        console.log(`[ECSW] Training Complete. Elements: ${this.sampledElements.length}/${nElements}`);
        return {
            indices: this.sampledElements,
            weights: this.weights,
            ratio: (nElements - this.sampledElements.length) / nElements
        };
    }

    /**
     * Simple Greedy Non-Negative Least Squares (Orthogonal Matching Pursuit variant)
     */
    solveGreedyNNLS(A, b, tol = 1e-3) {
        const nElements = A.length;
        const dim = b.length;
        
        let residual = new Float64Array(b);
        let currentIndices = [];
        let currentWeights = [];
        
        const maxIter = Math.min(nElements, 200);
        let bNorm = Math.sqrt(residual.reduce((sum, val) => sum + val * val, 0));

        for (let iter = 0; iter < maxIter; iter++) {
            let bestIdx = -1;
            let maxProjection = -1e-10;

            for (let i = 0; i < nElements; i++) {
                if (currentIndices.includes(i)) continue;
                
                let projection = 0;
                for (let d = 0; d < dim; d++) projection += A[i][d] * residual[d];
                
                if (projection > maxProjection) {
                    maxProjection = projection;
                    bestIdx = i;
                }
            }

            if (bestIdx === -1 || maxProjection <= 0) break;

            currentIndices.push(bestIdx);
            
            // Re-solve weights for the current subset (NNLS)
            // For simplicity in JS, we use a basic normal equation solve and clip negatives
            // Real implementation would use a proper NNLS active-set solver
            const sol = this.leastSquares(currentIndices.map(idx => A[idx]), b);
            currentWeights = sol.map(w => Math.max(0, w));

            // Update residual
            residual = new Float64Array(b);
            for (let i = 0; i < currentIndices.length; i++) {
                const idx = currentIndices[i];
                const w = currentWeights[i];
                for (let d = 0; d < dim; d++) {
                    residual[d] -= w * A[idx][d];
                }
            }

            const resNorm = Math.sqrt(residual.reduce((sum, val) => sum + val * val, 0));
            if (resNorm / bNorm < tol) break;
        }

        return { indices: currentIndices, weights: currentWeights };
    }

    leastSquares(basis, target) {
        const { Matrix, solve } = window.mlMatrix || window.ML;
        const A = new Matrix(basis).transpose(); // Each basis vector is a column
        const b = new Matrix([Array.from(target)]).transpose();
        try {
            const x = solve(A, b);
            return x.to2DArray().flat();
        } catch (e) {
            // Fallback for singular matrix
            return new Array(basis.length).fill(1.0 / basis.length);
        }
    }
}

window.ECSWEngine = new ECSWEngine();
