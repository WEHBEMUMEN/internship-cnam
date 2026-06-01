/**
 * Linear Algebra & Mathematical Utilities Library
 * Centralizes common mathematical functions used across solvers and simulation modules.
 */

class MathUtils {
    /**
     * Solve a linear system Ax = b using Gaussian elimination with partial pivoting.
     * @param {Array<Array<number>>} A Matrix
     * @param {Array<number>} b Vector
     * @returns {Array<number>} Solutions vector x
     */
    static solveGaussian(A, b) {
        const n = A.length;
        const mat = A.map(r => [...r]);
        const B = [...b];
        
        for (let i = 0; i < n; i++) {
            // Find pivot row
            let maxEl = Math.abs(mat[i][i]);
            let maxRow = i;
            for (let k = i + 1; k < n; k++) {
                if (Math.abs(mat[k][i]) > maxEl) {
                    maxEl = Math.abs(mat[k][i]);
                    maxRow = k;
                }
            }
            
            // Swap rows
            [mat[maxRow], mat[i]] = [mat[i], mat[maxRow]];
            [B[maxRow], B[i]] = [B[i], B[maxRow]];
            
            // Eliminate columns below pivot
            for (let k = i + 1; k < n; k++) {
                const c = -mat[k][i] / (mat[i][i] || 1e-12);
                for (let j = i; j < n; j++) {
                    mat[k][j] += c * mat[i][j];
                }
                B[k] += c * B[i];
            }
        }
        
        // Back substitution
        const x = new Array(n).fill(0);
        for (let i = n - 1; i >= 0; i--) {
            let sum = 0;
            for (let k = i + 1; k < n; k++) {
                sum += mat[i][k] * x[k];
            }
            x[i] = (B[i] - sum) / (mat[i][i] || 1e-12);
        }
        return x;
    }

    /**
     * Compute dot product of two vectors
     * @param {Array<number>} v1 
     * @param {Array<number>} v2 
     * @returns {number}
     */
    static dotProduct(v1, v2) {
        let sum = 0;
        const len = Math.min(v1.length, v2.length);
        for (let i = 0; i < len; i++) {
            sum += v1[i] * v2[i];
        }
        return sum;
    }

    /**
     * Matrix-vector multiplication (y = A * x)
     * @param {Array<Array<number>>} A 
     * @param {Array<number>} x 
     * @returns {Array<number>}
     */
    static matVecMultiply(A, x) {
        const n = A.length;
        const y = new Array(n).fill(0);
        for (let i = 0; i < n; i++) {
            let sum = 0;
            for (let j = 0; j < A[i].length; j++) {
                sum += A[i][j] * x[j];
            }
            y[i] = sum;
        }
        return y;
    }

    /**
     * Matrix transpose
     * @param {Array<Array<number>>} A 
     * @returns {Array<Array<number>>}
     */
    static transpose(A) {
        const rows = A.length;
        const cols = A[0].length;
        const AT = Array(cols).fill(0).map(() => new Array(rows));
        for (let i = 0; i < rows; i++) {
            for (let j = 0; j < cols; j++) {
                AT[j][i] = A[i][j];
            }
        }
        return AT;
    }
}

// Export for environment compatibility
if (typeof module !== 'undefined' && module.exports) {
    module.exports = MathUtils;
} else {
    window.MathUtils = MathUtils;
}
