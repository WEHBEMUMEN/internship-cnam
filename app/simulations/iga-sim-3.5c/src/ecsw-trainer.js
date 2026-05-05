/**
 * ECSW Offline Trainer
 * Handles snapshot processing and weight optimization.
 */
class ECSWTrainer {
    constructor() {
        this.allElements = [];
    }

    /**
     * Train ECSW weights.
     * @param {Object} context { fomSolver, romEngine, patch, snapshots, tolerance }
     */
    async train(context) {
        const { fomSolver, romEngine, patch, snapshots, tolerance = 1e-4 } = context;
        const Phi = romEngine.Phi;
        const k = Phi.columns;
        const nDofs = Phi.rows;
        const nSnaps = snapshots.length;

        this.allElements = this._getAllElements(patch);
        const nElements = this.allElements.length;

        console.log(`ECSW Trainer: Starting training on ${nElements} elements...`);

        // 1. Precompute all element projected forces Ge,s
        const elementProjections = this.allElements.map(el => 
            window.ECSWCore.computeElementProjectedForce(fomSolver, patch, el, snapshots, Phi)
        );

        // 2. Construct Target Matrix B [k*nSnaps x 1]
        // B contains projected assembly forces F_red,s = Phi^T * F_int,s
        const B = new Float64Array(k * nSnaps);
        const constrainedDofs = fomSolver._getConstrainedDofs ? fomSolver._getConstrainedDofs() : [];

        for (let s = 0; s < nSnaps; s++) {
            const F_int_assembly = fomSolver.calculateInternalForce(patch, snapshots[s]);
            
            // Zero out constrained DOFs to match snapshot training
            constrainedDofs.forEach(d => F_int_assembly[d] = 0);

            // Project to reduced space
            for (let i = 0; i < k; i++) {
                let dot = 0;
                for (let d = 0; d < nDofs; d++) {
                    dot += Phi.get(d, i) * F_int_assembly[d];
                }
                B[s * k + i] = dot;
            }
        }

        // 3. Build System Matrix A [k*nSnaps x nElements]
        const { Matrix } = window.mlMatrix;
        const A_data = [];
        for (let row = 0; row < k * nSnaps; row++) {
            const a_row = new Float64Array(nElements);
            for (let e = 0; e < nElements; e++) {
                a_row[e] = elementProjections[e][row];
            }
            A_data.push(Array.from(a_row));
        }
        
        const A_mat = new Matrix(A_data);
        const B_mat = new Matrix([Array.from(B)]).transpose();

        // 4. Solve NNLS via Core
        const w = window.ECSWCore.solveNNLS(A_mat, B_mat, { tolerance });

        // 5. Filter active elements
        const activeIndices = [];
        const weights = new Float64Array(nElements);
        for (let i = 0; i < nElements; i++) {
            if (w[i] > 1e-8) {
                weights[i] = w[i];
                activeIndices.push(i);
            }
        }

        const sampleElements = activeIndices.map(idx => ({
            ...this.allElements[idx],
            weight: weights[idx]
        }));

        // 6. Precompute Reduced Penalty (linear projection)
        console.log("ECSW Trainer: Projecting penalty constraints...");
        const Kp_full = Array.from({ length: nDofs }, () => new Float64Array(nDofs));
        fomSolver.applyPenaltyConstraints(Kp_full, null, new Float64Array(nDofs), patch);
        const Kp_mat = new Matrix(Kp_full);
        const Kp_red = Phi.transpose().mmul(Kp_mat).mmul(Phi).to2DArray();

        console.log(`ECSW Trainer: Training complete. ${sampleElements.length} elements selected.`);

        return {
            sampleElements,
            Kp_red,
            nDofs,
            k,
            weights
        };
    }

    _getAllElements(patch) {
        const { U, V } = patch;
        const uniqueU = [...new Set(U)], uniqueV = [...new Set(V)];
        const elements = [];
        for (let i = 0; i < uniqueU.length - 1; i++) {
            for (let j = 0; j < uniqueV.length - 1; j++) {
                elements.push({ 
                    i, j, 
                    uMin: uniqueU[i], uMax: uniqueU[i+1],
                    vMin: uniqueV[j], vMax: uniqueV[j+1]
                });
            }
        }
        return elements;
    }
}

window.ECSWTrainer = ECSWTrainer;
