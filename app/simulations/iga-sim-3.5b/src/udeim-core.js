/**
 * Phase 3.4b: U-DEIM Core Engine
 */
class UDEIMEngine {
    constructor() {
        this.U_f = null;       
        this.indices = null;   
        this.m = 0;            
        this.activeElements = null; 
        this.history = [];     
        this.M = null;         // Precomputed matrix M for F_red = M * f_sampled
        this.elementDofMap = [];
        this.numElements = 0;
        this.numLocalDofs = 0;
        this.sampledDofs = []; // List of { e: elementIndex, ld: localDofIndex }
        this.kf = 0;
        this.Kt_red_ref = null;
    }

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
        }
        return x;
    }
}
window.UDEIMEngine = UDEIMEngine;
