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
     * Piegl-Tiller Degree Elevation Workflow (Exact)
     */
    elevate1D(curve, newP, newU) {
        // Fast-Path: Direct Degree Elevation (DDE) for linear clamped NURBS (p=1 -> p=2)
        if (curve.p === 1 && newP === 2) {
            return this.applyDDE(curve);
        }

        // General Case: Piegl-Tiller Algorithm
        // 1. Decompose to Bezier segments (requires knot insertion)
        const segments = this.decomposeToBezier(curve);
        
        // 2. Elevate each Bezier segment
        const elevatedSegments = segments.map(seg => this.elevateBezierSegment(seg));
        
        // 3. Recompose into a single NURBS
        return this.recomposeFromBezier(elevatedSegments, curve.p, newP, newU);
    }

    /**
     * Direct Degree Elevation (DDE) optimization for linear curves
     * As requested: Q_{2i} = P_i, Q_{2i+1} = (P_i + P_{i+1})/2
     */
    applyDDE(curve) {
        const { controlPoints, weights } = curve;
        const n = controlPoints.length;
        const newCP = [];
        const newW = [];

        for (let i = 0; i < n; i++) {
            // Q_{2i} = P_i
            newCP.push({ ...controlPoints[i] });
            newW.push(weights[i]);

            // Q_{2i+1} = (P_i + P_{i+1})/2
            if (i < n - 1) {
                const p0 = controlPoints[i];
                const p1 = controlPoints[i+1];
                const w0 = weights[i];
                const w1 = weights[i+1];

                const midW = (w0 + w1) / 2;
                newCP.push({
                    x: (w0 * p0.x + w1 * p1.x) / (2 * midW),
                    y: (w0 * p0.y + w1 * p1.y) / (2 * midW),
                    z: (w0 * p0.z + w1 * p1.z) / (2 * midW)
                });
                newW.push(midW);
            }
        }
        return { controlPoints: newCP, weights: newW };
    }

    /**
     * Piegl-Tiller Step 1: Bézier Decomposition
     */
    decomposeToBezier(curve) {
        const { p, U, controlPoints, weights } = curve;
        let currentCurve = { ...curve };
        
        // Find internal knots with multiplicity < p
        const uniqueKnots = [...new Set(U.slice(p + 1, -p - 1))];
        for (const uBar of uniqueKnots) {
            const mult = U.filter(k => k === uBar).length;
            const insertCount = p - mult;
            for (let i = 0; i < insertCount; i++) {
                currentCurve = this.insertKnot1D(currentCurve, uBar);
            }
        }

        // Split into segments (each segment has p+1 points)
        const segments = [];
        const nSegments = uniqueKnots.length + 1;
        for (let s = 0; s < nSegments; s++) {
            const startIdx = s * p;
            segments.push({
                p,
                controlPoints: currentCurve.controlPoints.slice(startIdx, startIdx + p + 1),
                weights: currentCurve.weights.slice(startIdx, startIdx + p + 1)
            });
        }
        return segments;
    }

    /**
     * Piegl-Tiller Step 2: Segment Elevation
     * Q_i = (i/(p+1))*P_{i-1} + (1 - i/(p+1))*P_i
     */
    elevateBezierSegment(segment) {
        const { p, controlPoints, weights } = segment;
        const newP = p + 1;
        const newCP = [];
        const newW = [];

        // Convert to homogeneous coordinates
        const Pw = controlPoints.map((pt, i) => ({
            x: pt.x * weights[i],
            y: pt.y * weights[i],
            z: pt.z * weights[i],
            w: weights[i]
        }));

        for (let i = 0; i <= newP; i++) {
            let qw;
            if (i === 0) {
                qw = { ...Pw[0] };
            } else if (i === newP) {
                qw = { ...Pw[p] };
            } else {
                const alpha = i / newP;
                const p0 = Pw[i - 1];
                const p1 = Pw[i];
                qw = {
                    x: alpha * p0.x + (1 - alpha) * p1.x,
                    y: alpha * p0.y + (1 - alpha) * p1.y,
                    z: alpha * p0.z + (1 - alpha) * p1.z,
                    w: alpha * p0.w + (1 - alpha) * p1.w
                };
            }
            newCP.push({ x: qw.x / qw.w, y: qw.y / qw.w, z: qw.z / qw.w });
            newW.push(qw.w);
        }

        return { p: newP, controlPoints: newCP, weights: newW };
    }

    /**
     * Piegl-Tiller Step 3: Re-Composition
     * Note: For global elevation, we use the pre-computed elevated knot vector.
     */
    recomposeFromBezier(segments, oldP, newP, newU) {
        const newCP = [];
        const newW = [];

        for (let i = 0; i < segments.length; i++) {
            const seg = segments[i];
            const start = (i === 0) ? 0 : 1; // Avoid duplicate points at segment boundaries
            for (let j = start; j < seg.controlPoints.length; j++) {
                newCP.push(seg.controlPoints[j]);
                newW.push(seg.weights[j]);
            }
        }
        
        return { controlPoints: newCP, weights: newW };
    }

    /**
     * Internal 1D Knot Insertion for decomposition
     */
    insertKnot1D(curve, uBar) {
        const { p, U, controlPoints, weights } = curve;
        let k = -1;
        for (let i = 0; i < U.length - 1; i++) {
            if (uBar >= U[i] && uBar < U[i+1]) { k = i; break; }
        }
        if (k === -1) return curve;

        const newU = [...U];
        newU.splice(k + 1, 0, uBar);
        const newCP = [];
        const newW = [];

        for (let i = 0; i <= controlPoints.length; i++) {
            if (i <= k - p) {
                newCP.push({ ...controlPoints[i] });
                newW.push(weights[i]);
            } else if (i >= k + 1) {
                newCP.push({ ...controlPoints[i - 1] });
                newW.push(weights[i - 1]);
            } else {
                const alpha = (uBar - U[i]) / (U[i + p] - U[i]);
                const p0 = controlPoints[i - 1], p1 = controlPoints[i];
                const w0 = weights[i - 1], w1 = weights[i];
                const nw = (1 - alpha) * w0 + alpha * w1;
                newW.push(nw);
                newCP.push({
                    x: ((1 - alpha) * w0 * p0.x + alpha * w1 * p1.x) / (nw || 1),
                    y: ((1 - alpha) * w0 * p0.y + alpha * w1 * p1.y) / (nw || 1),
                    z: ((1 - alpha) * w0 * p0.z + alpha * w1 * p1.z) / (nw || 1)
                });
            }
        }
        return { p, U: newU, controlPoints: newCP, weights: newW };
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
