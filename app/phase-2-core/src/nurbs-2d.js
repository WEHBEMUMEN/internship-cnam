/**
 * 2D NURBS Core Engine
 * Handles bivariate basis functions and surface mapping.
 * Phase 2.0 | Computational Core
 */

class NURBS2D {
    constructor() {
        console.log("NURBS2D Engine Initialized");
    }

    /**
     * Compute 1D basis functions (Cox-de Boor)
     * @param {number} i Basis index
     * @param {number} p Degree
     * @param {Array} U Knot vector
     * @param {number} xi Parametric coordinate
     */
    basis1D(i, p, U, xi) {
        if (p === 0) {
            return (U[i] <= xi && xi < U[i + 1]) ? 1.0 : 0.0;
        }

        let denom1 = U[i + p] - U[i];
        let term1 = 0;
        if (denom1 > 0) {
            term1 = ((xi - U[i]) / denom1) * this.basis1D(i, p - 1, U, xi);
        }

        let denom2 = U[i + p + 1] - U[i + 1];
        let term2 = 0;
        if (denom2 > 0) {
            term2 = ((U[i + p + 1] - xi) / denom2) * this.basis1D(i + 1, p - 1, U, xi);
        }

        return term1 + term2;
    }

    /**
     * Compute Bivariate Basis Function R(i,j)
     * @param {number} i Index in U direction
     * @param {number} j Index in V direction
     * @param {object} patch NURBS surface data
     * @param {number} xi Parametric coordinate u
     * @param {number} eta Parametric coordinate v
     */
    bivariateBasis(i, j, patch, xi, eta) {
        const { p, q, U, V, weights } = patch;
        
        // Compute numerator: N_i(xi) * M_j(eta) * w_ij
        const N = this.basis1D(i, p, U, xi);
        const M = this.basis1D(j, q, V, eta);
        const w = weights[i][j];
        
        const numerator = N * M * w;

        // Compute denominator: sum of all weighted products
        let denominator = 0;
        for (let k = 0; k < weights.length; k++) {
            for (let l = 0; l < weights[0].length; l++) {
                const Nk = this.basis1D(k, p, U, xi);
                const Ml = this.basis1D(l, q, V, eta);
                denominator += Nk * Ml * weights[k][l];
            }
        }

        return denominator > 0 ? numerator / denominator : 0;
    }

    /**
     * Map (xi, eta) to (x, y, z)
     * @param {object} patch NURBS surface data
     * @param {number} xi Parametric coordinate u
     * @param {number} eta Parametric coordinate v
     */
    evaluateSurface(patch, xi, eta) {
        const { controlPoints, weights } = patch;
        let point = { x: 0, y: 0, z: 0 };

        // Handle edge case xi=1 or eta=1 (Cox-de Boor range)
        const eps = 1e-10;
        if (xi >= 1) xi = 1 - eps;
        if (eta >= 1) eta = 1 - eps;

        for (let i = 0; i < weights.length; i++) {
            for (let j = 0; j < weights[0].length; j++) {
                const R = this.bivariateBasis(i, j, patch, xi, eta);
                const cp = controlPoints[i][j];
                
                point.x += R * cp.x;
                point.y += R * cp.y;
                point.z += R * cp.z;
            }
        }

        return point;
    }

    /**
     * Verification: Partition of Unity
     */
    verifyPartitionOfUnity(patch, xi, eta) {
        let sum = 0;
        for (let i = 0; i < patch.weights.length; i++) {
            for (let j = 0; j < patch.weights[0].length; j++) {
                sum += this.bivariateBasis(i, j, patch, xi, eta);
            }
        }
        return sum;
    }
}

// Export for use in browser
if (typeof module !== 'undefined' && module.exports) {
    module.exports = NURBS2D;
} else {
    window.NURBS2D = NURBS2D;
}
