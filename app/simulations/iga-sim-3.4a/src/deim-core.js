/**
 * DEIM Engine — Phase 3.4a
 * Discrete Empirical Interpolation Method for Hyper-Reduction
 * 
 * Core Component: State management and shared utilities.
 */

class DEIMEngine {
    constructor() {
        this.U_f = null;
        this.indices = null;
        this.m = 0;
        this.kf = 0;
        this.PtU_pinv = null;
        this.activeElements = null;
        this.activeDofs = null;      // Subset of DOFs touched by active elements
        this._fBuf = null;           // Reusable buffer for assembly
        this.history = [];           // Stores step-by-step greedy selection history
        this.Kp_red = null;          // Reduced penalty stiffness
        this.Kt_red_ref = null;      // Reduced tangent stiffness reference
    }

    // ═══════════════════════════════════════════════════════════════
    //  STATIC UTILITIES
    // ═══════════════════════════════════════════════════════════════

    /**
     * Proper Orthogonal Decomposition (SVD) on snapshots.
     */
    static podVectors(snapshots, m) {
        const { Matrix, SVD } = window.mlMatrix;
        const S = new Matrix(snapshots.map(s => Array.from(s))).transpose();
        const svd = new SVD(S, {
            computeLeftSingularVectors: true,
            computeRightSingularVectors: false
        });
        const U = svd.leftSingularVectors;
        const trunc = Math.min(m, U.columns);
        return {
            basis: U.subMatrix(0, U.rows - 1, 0, trunc - 1),
            sigmas: svd.diagonal.slice(0, trunc)
        };
    }

    /**
     * Fast Gaussian elimination for small reduced systems.
     */
    static _solveLinear(A, b) {
        const n = b.length;
        const M = A.map(r => Array.from(r));
        const B = Array.from(b);

        for (let i = 0; i < n; i++) {
            let mx = i;
            for (let j = i + 1; j < n; j++)
                if (Math.abs(M[j][i]) > Math.abs(M[mx][i])) mx = j;
            [M[i], M[mx]] = [M[mx], M[i]];
            [B[i], B[mx]] = [B[mx], B[i]];
            if (Math.abs(M[i][i]) < 1e-20) M[i][i] = 1e-20;
            for (let j = i + 1; j < n; j++) {
                const f = M[j][i] / M[i][i];
                B[j] -= f * B[i];
                for (let k = i; k < n; k++) M[j][k] -= f * M[i][k];
            }
        }

        const x = new Float64Array(n);
        for (let i = n - 1; i >= 0; i--) {
            let s = 0;
            for (let j = i + 1; j < n; j++) s += M[i][j] * x[j];
            x[i] = (B[i] - s) / (M[i][i] || 1e-20);
            if (isNaN(x[i])) x[i] = 0;
        }
        return x;
    }
}

window.DEIMEngine = DEIMEngine;
