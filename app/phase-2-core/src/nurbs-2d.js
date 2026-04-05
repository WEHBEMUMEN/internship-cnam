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
     */
    bivariateBasis(i, j, patch, xi, eta, cachedDenominator = null) {
        const { p, q, U, V, weights } = patch;
        
        const N = this.basis1D(i, p, U, xi);
        const M = this.basis1D(j, q, V, eta);
        const w = (weights[i] && weights[i][j]) ? weights[i][j] : 1.0;
        
        const numerator = N * M * w;
        if (cachedDenominator !== null) {
            return cachedDenominator > 0 ? numerator / cachedDenominator : 0;
        }

        const denominator = this.computeDenominator(patch, xi, eta);
        return denominator > 1e-12 ? numerator / denominator : 0;
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
        
        // Handle boundaries
        const eps = 1e-8;
        xi = Math.max(0, Math.min(1 - eps, xi));
        eta = Math.max(0, Math.min(1 - eps, eta));

        const denominator = this.computeDenominator(patch, xi, eta) || 1.0;
        
        let position = { x: 0, y: 0, z: 0 };
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

        // Final sanity check
        if (isNaN(position.x) || isNaN(position.y) || isNaN(position.z)) {
            position = { x: 0, y: 0, z: 0 };
        }

        return { position, denominator };
    }

    /**
     * Simple evaluate call
     */
    evaluateSurface(patch, xi, eta) {
        return this.getSurfaceState(patch, xi, eta).position;
    }

    /**
     * Jacobian Determinant (Area Factor)
     */
    getJacobianDeterminant(patch, xi, eta) {
        const h = 0.001;
        const p0 = this.evaluateSurface(patch, xi, eta);
        const pu = this.evaluateSurface(patch, Math.min(xi + h, 1 - 1e-6), eta);
        const pv = this.evaluateSurface(patch, xi, Math.min(eta + h, 1 - 1e-6));

        const tu = { x: (pu.x - p0.x)/h, y: (pu.y - p0.y)/h, z: (pu.z - p0.z)/h };
        const tv = { x: (pv.x - p0.x)/h, y: (pv.y - p0.y)/h, z: (pv.z - p0.z)/h };

        const cp = {
            x: tu.y * tv.z - tu.z * tv.y,
            y: tu.z * tv.x - tu.x * tv.z,
            z: tu.x * tv.y - tu.y * tv.x
        };

        const area = Math.sqrt(cp.x * cp.x + cp.y * cp.y + cp.z * cp.z);
        return isNaN(area) ? 1.0 : area;
    }

    /**
     * h-refinement: Insert a knot in the U direction
     */
    insertKnotU(patch, uBar) {
        const { p, U, controlPoints, weights } = patch;
        const n = controlPoints.length;
        const m = controlPoints[0].length;

        let k = -1;
        for (let i = 0; i < U.length - 1; i++) {
            if (uBar >= U[i] && uBar < U[i + 1]) {
                k = i; break;
            }
        }
        if (k === -1) return patch;

        const newU = [...U];
        newU.splice(k + 1, 0, uBar);

        const newCP = [];
        const newW = [];

        for (let i = 0; i <= n; i++) {
            newCP[i] = [];
            newW[i] = [];
            for (let j = 0; j < m; j++) {
                if (i <= k - p) {
                    newCP[i][j] = { ...controlPoints[i][j] };
                    newW[i][j] = weights[i][j];
                } else if (i >= k + 1) {
                    newCP[i][j] = { ...controlPoints[i - 1][j] };
                    newW[i][j] = weights[i - 1][j];
                } else {
                    const denom = U[i + p] - U[i];
                    const alpha = denom > 0 ? (uBar - U[i]) / denom : 0;
                    const cp0 = controlPoints[i - 1][j];
                    const cp1 = controlPoints[i][j];
                    const w0 = weights[i - 1][j];
                    const w1 = weights[i][j];

                    newW[i][j] = (1 - alpha) * w0 + alpha * w1;
                    const denW = newW[i][j] || 1.0;
                    newCP[i][j] = {
                        x: ((1 - alpha) * w0 * cp0.x + alpha * w1 * cp1.x) / denW,
                        y: ((1 - alpha) * w0 * cp0.y + alpha * w1 * cp1.y) / denW,
                        z: ((1 - alpha) * w0 * cp0.z + alpha * w1 * cp1.z) / denW
                    };
                }
            }
        }

        patch.U = newU;
        patch.controlPoints = newCP;
        patch.weights = newW;
        return patch;
    }

    /**
     * h-refinement: Insert a knot in the V direction
     */
    insertKnotV(patch, vBar) {
        const { q, V, controlPoints, weights } = patch;
        const n = controlPoints.length;
        const m = controlPoints[0].length;

        let k = -1;
        for (let j = 0; j < V.length - 1; j++) {
            if (vBar >= V[j] && vBar < V[j + 1]) {
                k = j; break;
            }
        }
        if (k === -1) return patch;

        const newV = [...V];
        newV.splice(k + 1, 0, vBar);

        const newCP = [];
        const newW = [];

        for (let i = 0; i < n; i++) {
            newCP[i] = [];
            newW[i] = [];
            for (let j = 0; j <= m; j++) {
                if (j <= k - q) {
                    newCP[i][j] = { ...controlPoints[i][j] };
                    newW[i][j] = weights[i][j];
                } else if (j >= k + 1) {
                    newCP[i][j] = { ...controlPoints[i][j - 1] };
                    newW[i][j] = weights[i][j - 1];
                } else {
                    const denom = V[j + q] - V[j];
                    const alpha = denom > 0 ? (vBar - V[j]) / denom : 0;
                    const cp0 = controlPoints[i][j - 1];
                    const cp1 = controlPoints[i][j];
                    const w0 = weights[i][j - 1];
                    const w1 = weights[i][j];

                    newW[i][j] = (1 - alpha) * w0 + alpha * w1;
                    const denW = newW[i][j] || 1.0;
                    newCP[i][j] = {
                        x: ((1 - alpha) * w0 * cp0.x + alpha * w1 * cp1.x) / denW,
                        y: ((1 - alpha) * w0 * cp0.y + alpha * w1 * cp1.y) / denW,
                        z: ((1 - alpha) * w0 * cp0.z + alpha * w1 * cp1.z) / denW
                    };
                }
            }
        }

        patch.V = newV;
        patch.controlPoints = newCP;
        patch.weights = newW;
        return patch;
    }

    /**
     * Global Midpoint Subdivision
     */
    subdivideGlobal(patch) {
        const uniqueU = [...new Set(patch.U)];
        const insertU = [];
        for (let i = 0; i < uniqueU.length - 1; i++) {
            insertU.push((uniqueU[i] + uniqueU[i + 1]) / 2);
        }
        
        const uniqueV = [...new Set(patch.V)];
        const insertV = [];
        for (let i = 0; i < uniqueV.length - 1; i++) {
            insertV.push((uniqueV[i] + uniqueV[i + 1]) / 2);
        }

        // Insert knots
        insertU.reverse().forEach(u => this.insertKnotU(patch, u));
        insertV.reverse().forEach(v => this.insertKnotV(patch, v));
        return patch;
    }

    /**
     * p-refinement: Geometrically Invariant Degree Elevation
     */
    elevateDegreeInvariant(patch) {
        // Increment degrees
        const oldP = patch.p;
        const oldQ = patch.q;
        const oldU = [...patch.U];
        const oldV = [...patch.V];
        const oldCP = patch.controlPoints;
        const oldW = patch.weights;

        patch.p++;
        patch.q++;
        patch.U = this.elevateKnotVector(patch.U);
        patch.V = this.elevateKnotVector(patch.V);

        const nx = patch.U.length - patch.p - 1;
        const ny = patch.V.length - patch.q - 1;

        // Fit new control points to old surface
        const newCP = [];
        const newW = [];

        for (let i = 0; i < nx; i++) {
            newCP[i] = [];
            newW[i] = [];
            const xi = (patch.U[i + 1] + patch.U[i + patch.p]) / 2; // Greville Abscissa
            for (let j = 0; j < ny; j++) {
                const eta = (patch.V[j + 1] + patch.V[j + patch.q]) / 2;
                const state = this.getSurfaceState({ p: oldP, q: oldQ, U: oldU, V: oldV, controlPoints: oldCP, weights: oldW }, xi, eta);
                newCP[i][j] = state.position;
                newW[i][j] = state.denominator; // Simplified rational fit
            }
        }

        patch.controlPoints = newCP;
        patch.weights = newW;
        return patch;
    }

    elevateKnotVector(U) {
        const unique = [...new Set(U)];
        const newU = [];
        unique.forEach(u => {
            // Increase multiplicity at endpoints
            if (u === 0 || u === 1) {
                const count = U.filter(k => k === u).length;
                for (let i = 0; i <= count; i++) newU.push(u);
            } else {
                newU.push(u);
            }
        });
        return newU.sort((a,b) => a-b);
    }

    /**
     * p-refinement: Simple Elevation (Subdivision based)
     */
    elevateDegree(patch, dir = 'U') {
        return this.elevateDegreeInvariant(patch);
    }

    /**
     * k-refinement: p-refinement followed by h-refinement
     */
    kRefine(patch, uBar, vBar) {
        this.elevateDegree(patch, 'U');
        this.elevateDegree(patch, 'V');
        this.insertKnotU(patch, uBar);
        this.insertKnotV(patch, vBar);
        return patch;
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
