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
     * @param {number} cachedDenominator Optional pre-computed denominator
     */
    bivariateBasis(i, j, patch, xi, eta, cachedDenominator = null) {
        const { p, q, U, V, weights } = patch;
        
        const N = this.basis1D(i, p, U, xi);
        const M = this.basis1D(j, q, V, eta);
        const w = weights[i][j];
        
        const numerator = N * M * w;

        if (cachedDenominator !== null) {
            return cachedDenominator > 0 ? numerator / cachedDenominator : 0;
        }

        let denominator = 0;
        for (let k = 0; k < weights.length; k++) {
            const Nk = this.basis1D(k, p, U, xi);
            if (Nk === 0) continue;
            for (let l = 0; l < weights[0].length; l++) {
                const Ml = this.basis1D(l, q, V, eta);
                denominator += Nk * Ml * weights[k][l];
            }
        }

        return denominator > 0 ? numerator / denominator : 0;
    }

    /**
     * Compute Surface Denominator
     */
    computeDenominator(patch, xi, eta) {
        const { p, q, U, V, weights } = patch;
        let denominator = 0;
        for (let k = 0; k < weights.length; k++) {
            const Nk = this.basis1D(k, p, U, xi);
            if (Nk === 0) continue;
            for (let l = 0; l < weights[0].length; l++) {
                const Ml = this.basis1D(l, q, V, eta);
                denominator += Nk * Ml * weights[k][l];
            }
        }
        return denominator;
    }

    /**
     * Map (xi, eta) to (x, y, z) with full state data
     */
    getSurfaceState(patch, xi, eta) {
        const { controlPoints, weights, p, q, U, V } = patch;
        
        // Handle edge case xi=1 or eta=1
        const eps = 1e-10;
        if (xi >= 1) xi = 1 - eps;
        if (eta >= 1) eta = 1 - eps;

        const denominator = this.computeDenominator(patch, xi, eta);
        
        let position = { x: 0, y: 0, z: 0 };
        // In a real implementation, we would also compute derivatives here.
        // For Phase 2.1, we'll implement position first, then add analytical derivatives in 2.2.
        
        for (let i = 0; i < weights.length; i++) {
            const Ni = this.basis1D(i, p, U, xi);
            if (Ni === 0) continue;
            for (let j = 0; j < weights[0].length; j++) {
                const Mj = this.basis1D(j, q, V, eta);
                if (Mj === 0) continue;
                
                const R = (Ni * Mj * weights[i][j]) / denominator;
                const cp = controlPoints[i][j];
                
                position.x += R * cp.x;
                position.y += R * cp.y;
                position.z += R * cp.z;
            }
        }

        return { position, denominator };
    }

    /**
     * Simple evaluate call for legacy support
     */
    evaluateSurface(patch, xi, eta) {
        return this.getSurfaceState(patch, xi, eta).position;
    }

    /**
     * Jacobian Determinant (Numerical Approximation for Phase 2.1)
     */
    getJacobianDeterminant(patch, xi, eta) {
        const h = 0.001;
        const p0 = this.evaluateSurface(patch, xi, eta);
        const pu = this.evaluateSurface(patch, Math.min(xi + h, 1), eta);
        const pv = this.evaluateSurface(patch, xi, Math.min(eta + h, 1));

        const tu = { x: (pu.x - p0.x)/h, y: (pu.y - p0.y)/h, z: (pu.z - p0.z)/h };
        const tv = { x: (pv.x - p0.x)/h, y: (pv.y - p0.y)/h, z: (pv.z - p0.z)/h };

        // Cross product for Area
        const cp = {
            x: tu.y * tv.z - tu.z * tv.y,
            y: tu.z * tv.x - tu.x * tv.z,
            z: tu.x * tv.y - tu.y * tv.x
        };

        return Math.sqrt(cp.x*cp.x + cp.y*cp.y + cp.z*cp.z);
    }

    /**
     * Verification: Partition of Unity
     */
    verifyPartitionOfUnity(patch, xi, eta) {
        const denom = this.computeDenominator(patch, xi, eta);
        let sum = 0;
        for (let i = 0; i < patch.weights.length; i++) {
            for (let j = 0; j < patch.weights[0].length; j++) {
                sum += this.bivariateBasis(i, j, patch, xi, eta, denom);
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
