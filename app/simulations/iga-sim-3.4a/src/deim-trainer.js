/**
 * DEIM Engine — Offline Training Component
 */

DEIMEngine.prototype.train = function(forceSnapshots, m, kf = m, excludeDofs = []) {
    // Use full snapshots as requested to capture reaction forces
    const pod = DEIMEngine.podVectors(forceSnapshots, m);
    const sigmas = pod.sigmas;

    // Auto-truncate based on energy drop
    let effective_modes = 0;
    for (let i = 0; i < sigmas.length; i++) {
        if (sigmas[i] / sigmas[0] > 1e-12) effective_modes++;
    }
    if (effective_modes === 0) effective_modes = 1;

    this.m = Math.min(m, effective_modes);
    this.kf = Math.min(kf, this.m); 

    const U_greedy = pod.basis.subMatrix(0, pod.basis.rows - 1, 0, this.m - 1);
    const N = U_greedy.rows;

    // ═══════════════════════════════════════════════════════════════
    //  Q-DEIM: QR Factorization with Column Pivoting on U_f^T
    // ═══════════════════════════════════════════════════════════════
    let indices = [];
    this.history = [];

    // We want to pivot columns of U_greedy^T, which corresponds to rows of U_greedy.
    // Therefore, we perform pivoted Gram-Schmidt directly on the rows of U_greedy.
    const A = Array.from({ length: N }, (_, i) => {
        const row = new Float64Array(this.m);
        for (let j = 0; j < this.m; j++) row[j] = U_greedy.get(i, j);
        return row;
    });

    // Keep track of the squared 2-norm of each row
    const norms = new Float64Array(N);
    for (let i = 0; i < N; i++) {
        let sum = 0;
        for (let j = 0; j < this.m; j++) sum += A[i][j] * A[i][j];
        norms[i] = sum;
    }

    for (let l = 0; l < this.m; l++) {
        // Find row with max norm (Pivoting)
        let maxNorm = -1, pivot = -1;
        for (let i = 0; i < N; i++) {
            if (!indices.includes(i) && norms[i] > maxNorm && !excludeDofs.includes(i)) {
                maxNorm = norms[i];
                pivot = i;
            }
        }

        if (pivot === -1 || maxNorm < 1e-22) {
            console.warn(`DEIM: Residual dropped to noise floor (${maxNorm.toExponential(2)}) at step ${l + 1}. Truncating points.`);
            this.m = l;
            this.kf = Math.min(this.kf, this.m);
            break;
        }

        indices.push(pivot);
        this.history.push({ 
            step: l + 1, 
            point: pivot, 
            maxVal: Math.sqrt(maxNorm),
            residual: Array.from(norms).map(n => Math.sqrt(Math.max(0, n)))
        });

        // Orthogonalize remaining rows against the pivot row
        const pivotRow = A[pivot];
        const pNorm = Math.sqrt(maxNorm) || 1e-15;
        for (let j = 0; j < this.m; j++) pivotRow[j] /= pNorm; // Normalize pivot vector

        for (let i = 0; i < N; i++) {
            if (indices.includes(i)) continue;
            // dot product
            let dot = 0;
            for (let j = 0; j < this.m; j++) dot += A[i][j] * pivotRow[j];
            // subtract projection & recompute norm
            let newNorm = 0;
            for (let j = 0; j < this.m; j++) {
                A[i][j] -= dot * pivotRow[j];
                newNorm += A[i][j] * A[i][j];
            }
            norms[i] = newNorm;
        }
    }

    this.indices = indices;
    // Final reconstruction basis
    this.U_f = U_greedy.subMatrix(0, N - 1, 0, this.kf - 1);

    this._precomputeInterpolation();

    return {
        m: this.m,
        indices: [...this.indices],
        forceSigmas: pod.sigmas
    };
};

DEIMEngine.prototype._precomputeInterpolation = function() {
    // In Oversampled DEIM, PtU is m x kf. 
    // We precompute the pseudo-inverse: PtU_inv = (PtU^T PtU)^{-1} PtU^T   [kf x m]
    const m = this.m;
    const kf = this.kf;
    const { Matrix } = window.mlMatrix;

    const PtU_arr = Array.from({ length: m }, () => new Float64Array(kf));
    for (let i = 0; i < m; i++)
        for (let j = 0; j < kf; j++)
            PtU_arr[i][j] = this.U_f.get(this.indices[i], j);

    const PtU_mat = new Matrix(PtU_arr);
    const PtU_T = PtU_mat.transpose();

    // (PtU^T PtU)
    let M = PtU_T.mmul(PtU_mat);
    // Dynamic Tikhonov regularization to strictly cap the condition number at ~10^8
    let max_diag = 0;
    for (let i = 0; i < kf; i++) max_diag = Math.max(max_diag, M.get(i, i));
    const reg = Math.max(max_diag * 1e-8, 1e-14);
    for (let i = 0; i < kf; i++) M.set(i, i, M.get(i, i) + reg);

    // M_inv = (PtU^T PtU)^{-1}
    const M_inv = window.mlMatrix.inverse(M);

    // PtU_pinv = M_inv * PtU^T
    const PtU_pinv_mat = M_inv.mmul(PtU_T);

    // Store as regular array for fast access
    this.PtU_pinv = PtU_pinv_mat.to2DArray();
};
