/**
 * Phase 5.1a - Solver Bridge
 * Maps Physics Boundary Conditions and Loads to the IGA Solver.
 */

class SolverBridge {
    constructor(solver) {
        this.solver = solver;
    }

    /**
     * Solves the linear static system for a given patch
     */
    solveLinearStatic(patch, Tx = 100.0) {
        console.log(`[SolverBridge] Linear Solve started. Tx=${Tx}`);
        const nU = patch.controlPoints.length;
        const nV = patch.controlPoints[0].length;

        try {
            const bcs = [];
            for (let j = 0; j < nV; j++) bcs.push({ i: 0, j: j, axis: 'y', value: 0 });
            for (let j = 0; j < nU; j++) bcs.push({ i: j, j: 0, axis: 'x', value: 0 }); // Symmetry check
            
            // Re-apply 2B.2 logic exactly
            const symmetryBcs = [];
            for (let j = 0; j < nV; j++) symmetryBcs.push({ i: 0, j: j, axis: 'y', value: 0 });
            for (let j = 0; j < nV; j++) symmetryBcs.push({ i: nU - 1, j: j, axis: 'x', value: 0 });

            const integratedForces = this.solver.calculateNodalTraction(patch, Tx, 'right');
            const loads = [];
            for (let i = 0; i < integratedForces.length; i++) {
                if (Math.abs(integratedForces[i]) > 1e-12) {
                    const nodeIdx = Math.floor(i / 2);
                    const a = Math.floor(nodeIdx / nV);
                    const b = nodeIdx % nV;
                    if (i % 2 === 0) loads.push({ i: a, j: b, fx: integratedForces[i], fy: 0 });
                    else            loads.push({ i: a, j: b, fx: 0, fy: integratedForces[i] });
                }
            }
            
            const u = this.solver.solve(patch, symmetryBcs, loads);
            let maxU = 0;
            for (let i = 0; i < u.length; i++) if (Math.abs(u[i]) > maxU) maxU = Math.abs(u[i]);
            console.log(`[SolverBridge] Solve successful. u_max=${maxU.toFixed(8)}`);
            
            return { u };
        } catch (err) {
            console.error("[SolverBridge] Linear Solve Failed:", err);
            throw err;
        }
    }

    /**
     * Solves the system using a Reduced Basis (ROM)
     */
    solveReduced(patch, basis, Tx = 100.0) {
        console.log(`[SolverBridge] Reduced Solve started. Tx=${Tx}`);
        if (!basis || basis.length === 0) throw new Error("No basis provided for reduced solve.");

        const nU = patch.controlPoints.length;
        const nV = patch.controlPoints[0].length;
        const nDofs = nU * nV * 2;
        const k = basis[0].length;
        
        const { Matrix } = window.mlMatrix || window.ML;

        try {
            // 1. Prep BCs and Loads (Must match Linear exactly)
            const symmetryBcs = [];
            for (let j = 0; j < nV; j++) symmetryBcs.push({ i: 0, j: j, axis: 'y', value: 0 });
            for (let j = 0; j < nV; j++) symmetryBcs.push({ i: nU - 1, j: j, axis: 'x', value: 0 });
            
            const f_ext = this.solver.calculateNodalTraction(patch, Tx, 'right');
            const K_full = this.solver.assembleStiffness(patch);
            this.solver.applyPenaltyConstraints(K_full, patch);

            // 2. Identify Free Degrees of Freedom
            const fixedValues = new Map();
            symmetryBcs.forEach(bc => {
                const baseIdx = (bc.i * nV + bc.j) * 2;
                if (bc.axis === 'x' || bc.axis === 'both') fixedValues.set(baseIdx, bc.value);
                if (bc.axis === 'y' || bc.axis === 'both') fixedValues.set(baseIdx + 1, bc.value);
            });
            const freeIndices = [];
            for (let i = 0; i < nDofs; i++) if (!fixedValues.has(i)) freeIndices.push(i);
            const nFree = freeIndices.length;

            console.log(`[SolverBridge] ROM Dimensions: nDofs=${nDofs}, nFree=${nFree}, k=${k}`);

            // 3. Subset and Validate Basis
            if (basis.length !== nDofs) {
                throw new Error(`Basis row count (${basis.length}) does not match current system DOFs (${nDofs}). Please re-train.`);
            }

            // Create Phi_red as a plain 2D array for stability with ml-matrix
            const Phi_red_data = Array.from({ length: nFree }, () => new Float64Array(k));
            for (let i = 0; i < nFree; i++) {
                const globalIdx = freeIndices[i];
                for (let j = 0; j < k; j++) {
                    Phi_red_data[i][j] = basis[globalIdx][j];
                }
            }
            const Phi_red = new Matrix(Phi_red_data);

            // Create K_red as a plain 2D array
            const K_red_data = Array.from({ length: nFree }, () => new Float64Array(nFree));
            const f_red_data = new Float64Array(nFree);
            for (let i = 0; i < nFree; i++) {
                const rowIdx = freeIndices[i];
                f_red_data[i] = f_ext[rowIdx];
                for (let j = 0; j < nFree; j++) {
                    K_red_data[i][j] = K_full[rowIdx][freeIndices[j]];
                }
            }
            const K_red_mat = new Matrix(K_red_data);
            const f_red_mat = new Matrix([Array.from(f_red_data)]).transpose();

            // 4. Projection
            console.log("[SolverBridge] Projecting to reduced space...");
            const Kr = Phi_red.transpose().mmul(K_red_mat).mmul(Phi_red);
            const fr = Phi_red.transpose().mmul(f_red_mat);

            // 5. Solve Reduced System
            console.log("[SolverBridge] Solving Kr * ar = fr...");
            const ar = window.PODEngine.solveLinear(Kr.to2DArray(), fr.to2DArray().flat());

            // 6. Reconstruct Full Solution
            const u_free = Phi_red.mmul(new Matrix([ar]).transpose()).to2DArray().flat();
            const u_full = new Float64Array(nDofs);
            let freeCounter = 0;
            for (let i = 0; i < nDofs; i++) {
                if (fixedValues.has(i)) {
                    u_full[i] = fixedValues.get(i);
                } else {
                    u_full[i] = u_free[freeCounter++];
                }
            }

            console.log(`[SolverBridge] Reduced Solve successful. u_max=${Math.max(...u_full.map(Math.abs)).toFixed(8)}`);
            return { u: u_full };
        } catch (err) {
            console.error("[SolverBridge] Reduced Solve Failed:", err);
            throw err;
        }
    }
}

window.SolverBridge = SolverBridge;
