/**
 * DEIM Engine — Setup & Pre-computation Component
 */

/**
 * Map DEIM indices to the minimal set of elements (knot spans) that must be assembled.
 * This is the secret to 100x speedup.
 */
DEIMEngine.prototype.computeActiveElements = function(patch) {
    const { p, q, U, V, controlPoints } = patch;
    const nV = controlPoints[0].length;
    const uniqueU = [...new Set(U)], uniqueV = [...new Set(V)];
    const elements = [];

    // For each DEIM index, find which elements contribute to it
    this.indices.forEach(idx => {
        const cpIdx = Math.floor(idx / 2);
        const cpI = Math.floor(cpIdx / nV);
        const cpJ = cpIdx % nV;

        // Control point (cpI, cpJ) is influenced by elements in spans:
        // i in [cpI-p, cpI] (where i is index of uniqueU span)
        // But we can just find which spans have this basis function non-zero.
        for (let i = 0; i < uniqueU.length - 1; i++) {
            const uMid = (uniqueU[i] + uniqueU[i + 1]) / 2;
            const spanU = window.nurbsUtils ? window.nurbsUtils.findSpan(controlPoints.length - 1, p, uMid, U) : i + p; // Fallback estimate

            // Basis function N_{cpI,p} is non-zero if spanU is in [cpI, cpI+p]
            // Correct logic: findSpan returns the index k such that u \in [u_k, u_{k+1}]
            // The basis functions non-zero on [u_k, u_{k+1}] are N_{k-p}, ..., N_k.

            for (let j = 0; j < uniqueV.length - 1; j++) {
                const vMid = (uniqueV[j] + uniqueV[j + 1]) / 2;

                // Optimization: Check if cpI is in the support of this element
                // For element (i, j) defined by [uniqueU[i], uniqueU[i+1]] x [uniqueV[j], uniqueV[j+1]]
                // we need to find the knot indices
                const kU = U.lastIndexOf(uniqueU[i]);
                const kV = V.lastIndexOf(uniqueV[j]);

                // Basis N_{A,p} is non-zero on [u_k, u_{k+1}] if k-p <= A <= k
                if (cpI >= kU - p && cpI <= kU && cpJ >= kV - q && cpJ <= kV) {
                    const elKey = `${i}-${j}`;
                    if (!elements.some(e => e.key === elKey)) {
                        elements.push({ i, j, key: elKey, uMin: uniqueU[i], uMax: uniqueU[i + 1], vMin: uniqueV[j], vMax: uniqueV[j + 1] });
                    }
                }
            }
        }
    });
    this.activeElements = elements;

    // Identify all DOFs touched by these elements
    const dofSet = new Set();
    elements.forEach(el => {
        // Knot span i, j corresponds to a set of basis functions
        // N_{k-p, p} ... N_{k, p}
        const kU = U.lastIndexOf(uniqueU[el.i]);
        const kV = V.lastIndexOf(uniqueV[el.j]);
        for (let ii = kU - p; ii <= kU; ii++) {
            for (let jj = kV - q; jj <= kV; jj++) {
                const cpIdx = ii * nV + jj;
                dofSet.add(cpIdx * 2);
                dofSet.add(cpIdx * 2 + 1);
            }
        }
    });
    this.activeDofs = Array.from(dofSet);
    const nDofs = patch.controlPoints.length * nV * 2;
    this._fBuf = new Float64Array(nDofs);

    console.log(`DEIM: Hyper-Reduction ready. Elements: ${elements.length}, Active DOFs: ${this.activeDofs.length}`);
};

/**
 * Pre-compute the reduced tangent stiffness from training snapshots.
 * Uses the average tangent over all snapshots for the online phase.
 * This is the key to avoiding O(N²) tangent assembly online.
 */
DEIMEngine.prototype.precomputeReducedTangent = function(fomSolver, romEngine, patch, snapU) {
    const { Matrix } = window.mlMatrix;
    const Phi = romEngine.Phi;
    const k = Phi.columns;
    const nDofs = Phi.rows;
    const PhiT = Phi.transpose();

    // Average the reduced tangent over multiple snapshots
    const Kt_red_avg = Array.from({ length: k }, () => new Float64Array(k));
    const nSnaps = snapU.length;

    for (let s = 0; s < nSnaps; s++) {
        const Kt_full = fomSolver.calculateTangentStiffness(patch, snapU[s]);
        fomSolver.applyPenaltyConstraints(Kt_full, null, snapU[s], patch);
        const Kt_mat = new Matrix(Kt_full);
        const Kt_red = PhiT.mmul(Kt_mat).mmul(Phi).to2DArray();
        for (let i = 0; i < k; i++)
            for (let j = 0; j < k; j++)
                Kt_red_avg[i][j] += Kt_red[i][j] / nSnaps;
    }

    this.Kt_red_ref = Kt_red_avg;
};

/**
 * Pre-compute the reduced penalty stiffness matrix.
 * Penalty is a linear operator — no DEIM approximation needed.
 * Kp_red = Φ^T * K_penalty * Φ  (exact, precomputed once)
 */
DEIMEngine.prototype.precomputeReducedPenalty = function(fomSolver, romEngine, patch) {
    const { Matrix } = window.mlMatrix;
    const Phi = romEngine.Phi;
    const k = Phi.columns;
    const nDofs = Phi.rows;
    const PhiT = Phi.transpose();

    // Build penalty-only stiffness (pass zero displacement → pure penalty structure)
    const Kp = Array.from({ length: nDofs }, () => new Float64Array(nDofs));
    fomSolver.applyPenaltyConstraints(Kp, null, new Float64Array(nDofs), patch);

    // Project to reduced space
    const Kp_mat = new Matrix(Kp.map(r => Array.from(r)));
    this.Kp_red = PhiT.mmul(Kp_mat).mmul(Phi).to2DArray();
    console.log(`DEIM: Reduced penalty matrix precomputed (${k}x${k})`);
};
