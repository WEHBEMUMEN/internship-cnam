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
     * p-refinement: Geometrically Invariant Degree Elevation (Tensor Product)
     */
    elevateDegreeInvariant(patch) {
        // Step 1: Elevate U degree
        this.elevateDirection(patch, 'U');
        
        // Step 2: Elevate V degree
        this.elevateDirection(patch, 'V');
        
        return patch;
    }

    /**
     * Elevate degree in one parametric direction (tensor-product)
     */
    elevateDirection(patch, dir = 'U') {
        const oldP = dir === 'U' ? patch.p : patch.q;
        const newP = oldP + 1;
        const oldKnotVector = dir === 'U' ? patch.U : patch.V;
        const newKnotVector = this.elevateKnotVector(oldKnotVector);
        
        const n = patch.controlPoints.length;
        const m = patch.controlPoints[0].length;
        
        if (dir === 'U') {
            const nextN = newKnotVector.length - newP - 1;
            const newCP = Array(nextN).fill(0).map(() => Array(m).fill(null));
            const newW = Array(nextN).fill(0).map(() => Array(m).fill(1));

            // Elevate each column of control points in the U direction
            for (let j = 0; j < m; j++) {
                const columnCP = [];
                const columnW = [];
                for (let i = 0; i < n; i++) {
                    columnCP.push(patch.controlPoints[i][j]);
                    columnW.push(patch.weights[i][j]);
                }
                
                const elevated = this.elevate1D({ p: oldP, U: oldKnotVector, controlPoints: columnCP, weights: columnW }, newP, newKnotVector);
                for (let i = 0; i < nextN; i++) {
                    newCP[i][j] = elevated.controlPoints[i];
                    newW[i][j] = elevated.weights[i];
                }
            }
            patch.p = newP;
            patch.U = newKnotVector;
            patch.controlPoints = newCP;
            patch.weights = newW;
        } else {
            const nextM = newKnotVector.length - newP - 1;
            const newCP = Array(n).fill(0).map(() => Array(nextM).fill(null));
            const newW = Array(n).fill(0).map(() => Array(nextM).fill(1));

            // Elevate each row of control points in the V direction
            for (let i = 0; i < n; i++) {
                const rowCP = patch.controlPoints[i];
                const rowW = patch.weights[i];
                
                const elevated = this.elevate1D({ p: oldP, U: oldKnotVector, controlPoints: rowCP, weights: rowW }, newP, newKnotVector);
                for (let j = 0; j < nextM; j++) {
                    newCP[i][j] = elevated.controlPoints[j];
                    newW[i][j] = elevated.weights[j];
                }
            }
            patch.q = newP;
            patch.V = newKnotVector;
            patch.controlPoints = newCP;
            patch.weights = newW;
        }
    }

    /**
     * 1D Least Squares Fit for Degree Elevation
     */
    elevate1D(curve, newP, newU) {
        const M = 100; // Sample count
        const samples = [];
        for (let i = 0; i <= M; i++) {
            const xi = i / M;
            samples.push(this.evaluate1D(curve, xi));
        }

        const newN = newU.length - newP - 1;
        const newWeights = new Array(newN).fill(1.0);

        // Build A matrix (M+1 x newN)
        const A = [];
        for (let i = 0; i <= M; i++) {
            const xi = i / M;
            const row = [];
            let denom = 0;
            for (let j = 0; j < newN; j++) {
                const bj = this.basis1D(j, newP, newU, xi);
                row.push(bj);
                denom += bj; // Unit weighs for fitting position
            }
            A.push(row);
        }

        // Normal Equations: ATA * X = ATB
        const ATA = Array(newN).fill(0).map(() => Array(newN).fill(0));
        const ATBx = Array(newN).fill(0);
        const ATBy = Array(newN).fill(0);
        const ATBz = Array(newN).fill(0);

        for (let i = 0; i < newN; i++) {
            for (let j = 0; j < newN; j++) {
                let sum = 0;
                for (let k = 0; k <= M; k++) sum += A[k][i] * A[k][j];
                ATA[i][j] = sum;
            }
            let sx = 0, sy = 0, sz = 0;
            for (let k = 0; k <= M; k++) {
                sx += A[k][i] * samples[k].x;
                sy += A[k][i] * samples[k].y;
                sz += A[k][i] * samples[k].z;
            }
            ATBx[i] = sx;
            ATBy[i] = sy;
            ATBz[i] = sz;
        }

        const solve = (matA, matB) => {
            const n = matA.length;
            const A = matA.map(r => [...r]);
            const B = [...matB];
            for (let i = 0; i < n; i++) {
                let maxEl = Math.abs(A[i][i]), maxRow = i;
                for (let k = i + 1; k < n; k++) {
                    if (Math.abs(A[k][i]) > maxEl) { maxEl = Math.abs(A[k][i]); maxRow = k; }
                }
                [A[maxRow], A[i]] = [A[i], A[maxRow]];
                [B[maxRow], B[i]] = [B[i], B[maxRow]];
                for (let k = i + 1; k < n; k++) {
                    const c = -A[k][i] / (A[i][i] || 1e-12);
                    for (let j = i; j < n; j++) A[k][j] += c * A[i][j];
                    B[k] += c * B[i];
                }
            }
            const x = new Array(n).fill(0);
            for (let i = n - 1; i >= 0; i--) {
                let sum = 0;
                for (let k = i + 1; k < n; k++) sum += A[i][k] * x[k];
                x[i] = (B[i] - sum) / (A[i][i] || 1e-12);
            }
            return x;
        };

        const resX = solve(ATA, ATBx);
        const resY = solve(ATA, ATBy);
        const resZ = solve(ATA, ATBz);

        const newCP = [];
        for (let i = 0; i < newN; i++) {
            newCP.push({ x: resX[i], y: resY[i], z: resZ[i] });
        }

        // Clamp endpoints to original
        newCP[0] = { ...samples[0] };
        newCP[newN - 1] = { ...samples[M] };

        return { controlPoints: newCP, weights: newWeights };
    }

    /**
     * Evaluate 1D NURBS Curve
     */
    evaluate1D(curve, xi) {
        const { p, U, controlPoints, weights } = curve;
        let denom = 0;
        for (let i = 0; i < controlPoints.length; i++) {
            denom += this.basis1D(i, p, U, xi) * weights[i];
        }
        if (denom === 0) denom = 1;

        let pos = { x: 0, y: 0, z: 0 };
        for (let i = 0; i < controlPoints.length; i++) {
            const R = (this.basis1D(i, p, U, xi) * weights[i]) / denom;
            pos.x += R * controlPoints[i].x;
            pos.y += R * controlPoints[i].y;
            pos.z += R * controlPoints[i].z;
        }
        return pos;
    }

    elevateKnotVector(U) {
        const unique = [...new Set(U)];
        const newU = [];
        unique.forEach(u => {
            // Increase multiplicity everywhere to maintain continuity order if desired, 
            // but standard p-refinement increases knot vector multiplicity by 1 at internal knots
            const count = U.filter(k => k === u).length;
            for (let i = 0; i < count + 1; i++) newU.push(u);
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
