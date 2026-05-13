/**
 * Phase 5.1b - Solver Bridge
 * Connects the App State to the optimized High-Performance Solvers.
 */

class SolverBridge {
    constructor(solver) {
        this.solver = solver;
    }

    /**
     * Solves the full system using standard Linear IGA
     */
    solveLinearStatic(patch, Tx = 100.0) {
        console.log(`[SolverBridge] Linear Solve started. Tx=${Tx}`);
        const nU = patch.controlPoints.length;
        const nV = patch.controlPoints[0].length;
        
        try {
            // Cantilever BCs: Clamped on the left (i=0)
            const clampedBcs = [];
            for (let j = 0; j < nV; j++) {
                clampedBcs.push({ i: 0, j: j, axis: 'both', value: 0 });
            }
            
            const u = this.solver.solve(patch, clampedBcs, []);
            
            // Manual super-position for traction if needed, 
            // but the solver should handle it if passed correctly.
            // Re-solving with nodal traction:
            const f_ext = this.solver.calculateNodalTraction(patch, Tx, 'right');
            const K_full = this.solver.assembleStiffness(patch);
            this.solver.applyPenaltyConstraints(K_full, patch);
            
            const fixedValues = new Map();
            clampedBcs.forEach(bc => {
                const baseIdx = (bc.i * nV + bc.j) * 2;
                if (bc.axis === 'x' || bc.axis === 'both') fixedValues.set(baseIdx, bc.value);
                if (bc.axis === 'y' || bc.axis === 'both') fixedValues.set(baseIdx + 1, bc.value);
            });

            const nDofs = nU * nV * 2;
            const freeIndices = [];
            for (let i = 0; i < nDofs; i++) if (!fixedValues.has(i)) freeIndices.push(i);
            const nFree = freeIndices.length;

            const K_red = Array.from({ length: nFree }, () => new Float64Array(nFree));
            const F_red = new Float64Array(nFree);

            for (let i = 0; i < nFree; i++) {
                const row = freeIndices[i];
                F_red[i] = f_ext[row];
                for (let j = 0; j < nFree; j++) {
                    K_red[i][j] = K_full[row][freeIndices[j]];
                }
            }

            const sol_red = this.solver.gaussianElimination(K_red, F_red);
            const u_full = new Float64Array(nDofs);
            let freeCounter = 0;
            for (let i = 0; i < nDofs; i++) {
                if (fixedValues.has(i)) u_full[i] = fixedValues.get(i);
                else                    u_full[i] = sol_red[freeCounter++];
            }

            return { u: u_full };
        } catch (err) {
            console.error("[SolverBridge] Linear Solve Failed:", err);
            throw err;
        }
    }

    /**
     * Solves the system using standard Galerkin ROM
     */
    solveReduced(patch, basis, Tx = 100.0) {
        console.log(`[SolverBridge] Reduced Solve started. Tx=${Tx}`);
        const nU = patch.controlPoints.length;
        const nV = patch.controlPoints[0].length;
        const nDofs = nU * nV * 2;
        const k = basis[0].length;
        const { Matrix } = window.mlMatrix || window.ML;
        const Phi = new Matrix(basis);

        try {
            // 1. Assemble Reduced Stiffness (Optimized element-wise)
            const Kr = this.solver.assembleReducedStiffness(patch, Phi);
            
            // 2. Reduce Force Vector
            const f_ext = this.solver.calculateNodalTraction(patch, Tx, 'right');
            const fr = Phi.transpose().mmul(new Matrix([Array.from(f_ext)]).transpose());

            // 3. Tiny k x k solve
            const ar = window.PODEngine.solveLinear(Kr.to2DArray(), fr.to2DArray().flat());
            
            // 4. Reconstruct
            const u_full = Phi.mmul(new Matrix([ar]).transpose()).to2DArray().flat();

            return { u: u_full };
        } catch (err) {
            console.error("[SolverBridge] Reduced Solve Failed:", err);
            throw err;
        }
    }

    /**
     * Solves the system using ECSW Hyper-reduction
     */
    solveECSW(patch, basis, ecswData, Tx = 100.0) {
        console.log(`[SolverBridge] ECSW Solve started. Tx=${Tx}`);
        const nU = patch.controlPoints.length;
        const nV = patch.controlPoints[0].length;
        const nDofs = nU * nV * 2;
        const k = basis[0].length;
        const { Matrix } = window.mlMatrix || window.ML;
        const Phi = new Matrix(basis);

        try {
            // 1. Direct Reduced Assembly (Hyper-reduction)
            const Kr = this.solver.assembleECSWReduced(patch, ecswData.indices, ecswData.weights, Phi);
            
            // 2. Reduce Force Vector
            const f_ext = this.solver.calculateNodalTraction(patch, Tx, 'right');
            const fr = Phi.transpose().mmul(new Matrix([Array.from(f_ext)]).transpose());

            // 3. Tiny k x k solve
            const ar = window.PODEngine.solveLinear(Kr.to2DArray(), fr.to2DArray().flat());
            
            // 4. Reconstruct
            const u_full = Phi.mmul(new Matrix([ar]).transpose()).to2DArray().flat();

            return { u: u_full };
        } catch (err) {
            console.error("[SolverBridge] ECSW Solve Failed:", err);
            throw err;
        }
    }
}

window.SolverBridge = SolverBridge;
