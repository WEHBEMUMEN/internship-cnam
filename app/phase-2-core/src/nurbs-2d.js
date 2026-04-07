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
            // Handle right boundary exactly for the last knot span
            const isRightBoundary = (xi === 1.0 && U[i + 1] === 1.0);
            return (U[i] <= xi && (xi < U[i + 1] || isRightBoundary)) ? 1.0 : 0.0;
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
     * Compute first derivative of 1D basis function N'(i,p)
     */
    basis1DDeriv(i, p, U, xi) {
        if (p === 0) return 0.0;

        let denom1 = U[i + p] - U[i];
        let term1 = 0;
        if (denom1 > 0) {
            term1 = (p / denom1) * this.basis1D(i, p - 1, U, xi);
        }

        let denom2 = U[i + p + 1] - U[i + 1];
        let term2 = 0;
        if (denom2 > 0) {
            term2 = (p / denom2) * this.basis1D(i + 1, p - 1, U, xi);
        }

        return term1 - term2;
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
     * Analytical Surface Derivatives (Quotient Rule)
     * Returns { pos: {x,y,z}, dU: {x,y,z}, dV: {x,y,z}, W, dWu, dWv }
     */
    getSurfaceDerivatives(patch, xi, eta) {
        const { controlPoints, weights, p, q, U, V } = patch;
        
        let A = { x: 0, y: 0, z: 0 };
        let Au = { x: 0, y: 0, z: 0 };
        let Av = { x: 0, y: 0, z: 0 };
        let W = 0;
        let Wu = 0;
        let Wv = 0;

        for (let i = 0; i < weights.length; i++) {
            const Ni = this.basis1D(i, p, U, xi);
            const dNi = this.basis1DDeriv(i, p, U, xi);
            if (Ni === 0 && dNi === 0) continue;

            for (let j = 0; j < weights[0].length; j++) {
                const Mj = this.basis1D(j, q, V, eta);
                const dMj = this.basis1DDeriv(j, q, V, eta);
                if (Mj === 0 && dMj === 0) continue;

                const w = weights[i][j];
                const cp = controlPoints[i][j];

                // Numerator sum
                const basis = Ni * Mj * w;
                const basisU = dNi * Mj * w;
                const basisV = Ni * dMj * w;

                A.x += basis * cp.x; A.y += basis * cp.y; A.z += basis * cp.z;
                Au.x += basisU * cp.x; Au.y += basisU * cp.y; Au.z += basisU * cp.z;
                Av.x += basisV * cp.x; Av.y += basisV * cp.y; Av.z += basisV * cp.z;

                // Denominator sum
                W += basis;
                Wu += basisU;
                Wv += basisV;
            }
        }

        if (W < 1e-12) W = 1e-12;

        const pos = { x: A.x / W, y: A.y / W, z: A.z / W };
        
        // Exact Quotient Rule: dS/du = (Au*W - A*Wu) / W^2 = (Au - S*Wu) / W
        const dU = {
            x: (Au.x - pos.x * Wu) / W,
            y: (Au.y - pos.y * Wu) / W,
            z: (Au.z - pos.z * Wu) / W
        };

        const dV = {
            x: (Av.x - pos.x * Wv) / W,
            y: (Av.y - pos.y * Wv) / W,
            z: (Av.z - pos.z * Wv) / W
        };

        return { pos, dU, dV, W, Wu, Wv };
    }

    /**
     * Jacobian Determinant (Exact Analytical Cross Product)
     */
    getJacobianDeterminant(patch, xi, eta) {
        const deriv = this.getSurfaceDerivatives(patch, xi, eta);
        const tu = deriv.dU;
        const tv = deriv.dV;

        // Cross Product (Norm of Normal Vector)
        const nx = tu.y * tv.z - tu.z * tv.y;
        const ny = tu.z * tv.x - tu.x * tv.z;
        const nz = tu.x * tv.y - tu.y * tv.x;

        let area = Math.sqrt(nx * nx + ny * ny + nz * nz);
        
        // Safety: If area is zero at degenerate corner, use a stable floor 
        // to maintain matrix stability rather than zeroing out stiffness.
        if (isNaN(area) || area < 1e-12) area = 1e-12;
        
        return area;
    }

    /**
     * Get Local Surface Normal (Unit Vector)
     */
    getSurfaceNormal(patch, xi, eta) {
        const deriv = this.getSurfaceDerivatives(patch, xi, eta);
        const tu = deriv.dU;
        const tv = deriv.dV;

        const n = {
            x: tu.y * tv.z - tu.z * tv.y,
            y: tu.z * tv.x - tu.x * tv.z,
            z: tu.x * tv.y - tu.y * tv.x
        };

        const mag = Math.sqrt(n.x * n.x + n.y * n.y + n.z * n.z) || 1.0;
        return { x: n.x / mag, y: n.y / mag, z: n.z / mag };
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
     * Set specific degrees for p and q independently
     */
    setDegree(patch, targetP, targetQ) {
        // Handle U direction (p)
        while (patch.p < targetP) this.elevateDirection(patch, 'U');
        while (patch.p > targetP) this.reduceDirection(patch, 'U');

        // Handle V direction (q)
        while (patch.q < targetQ) this.elevateDirection(patch, 'V');
        while (patch.q > targetQ) this.reduceDirection(patch, 'V');
        
        return patch;
    }

    /**
     * p-refinement: Geometrically Invariant Degree Elevation (Tensor Product)
     */
    elevateDegreeInvariant(patch) {
        this.elevateDirection(patch, 'U');
        this.elevateDirection(patch, 'V');
        return patch;
    }

    /**
     * Degree Reduction in one parametric direction
     */
    reduceDirection(patch, dir = 'U') {
        const oldP = dir === 'U' ? patch.p : patch.q;
        const newP = oldP - 1;
        if (newP < 1) return;

        const oldKnotVector = dir === 'U' ? patch.U : patch.V;
        const n = patch.controlPoints.length;
        const m = patch.controlPoints[0].length;

        if (dir === 'U') {
            const nextN = n - 1;
            const newCP = Array(nextN).fill(0).map(() => Array(m).fill(null));
            const newW = Array(nextN).fill(0).map(() => Array(m).fill(1));

            for (let j = 0; j < m; j++) {
                const columnCP = [];
                const columnW = [];
                for (let i = 0; i < n; i++) {
                    columnCP.push(patch.controlPoints[i][j]);
                    columnW.push(patch.weights[i][j]);
                }
                const reduced = this.reduce1D({ p: oldP, U: oldKnotVector, controlPoints: columnCP, weights: columnW }, newP);
                for (let i = 0; i < nextN; i++) {
                    newCP[i][j] = reduced.controlPoints[i];
                    newW[i][j] = reduced.weights[i];
                }
                if (j === 0) patch.U = reduced.U;
            }
            patch.p = newP;
            patch.controlPoints = newCP;
            patch.weights = newW;
        } else {
            const nextM = m - 1;
            const newCP = Array(n).fill(0).map(() => Array(nextM).fill(null));
            const newW = Array(n).fill(0).map(() => Array(nextM).fill(1));

            for (let i = 0; i < n; i++) {
                const rowCP = patch.controlPoints[i];
                const rowW = patch.weights[i];
                const reduced = this.reduce1D({ q: oldP, V: oldKnotVector, controlPoints: rowCP, weights: rowW }, newP);
                for (let j = 0; j < nextM; j++) {
                    newCP[i][j] = reduced.controlPoints[j];
                    newW[i][j] = reduced.weights[j];
                }
                if (i === 0) patch.V = reduced.V;
            }
            patch.q = newP;
            patch.controlPoints = newCP;
            patch.weights = newW;
        }
    }

    /**
     * 1D Degree Reduction via Least Squares Fitting
     * More robust than inward inversion for arbitrary control point layouts.
     */
    reduce1D(curve, targetP) {
        const { controlPoints, weights, U, V, p, q } = curve;
        const currentP = p || q;
        const currentU = U || V;
        
        // 1. Generate target knot vector (remove one internal knot or maintain clamped)
        const newU = this.reduceKnotVector(currentU);
        const newN = newU.length - targetP - 1;

        // 2. Sample the original curve at high resolution
        const M = Math.max(100, newN * 4);
        const samples = [];
        for (let i = 0; i <= M; i++) {
            const xi = i / M;
            let pos = { x: 0, y: 0, z: 0, w: 0 };
            for (let j = 0; j < controlPoints.length; j++) {
                const N = this.basis1D(j, currentP, currentU, xi);
                const w = weights[j];
                pos.x += N * controlPoints[j].x * w;
                pos.y += N * controlPoints[j].y * w;
                pos.z += N * controlPoints[j].z * w;
                pos.w += N * w;
            }
            samples.push(pos);
        }

        // 3. Solve LSQ System: A^T * A * Q = A^T * B
        const A = [];
        for (let i = 0; i <= M; i++) {
            const xi = i / M;
            const row = [];
            for (let j = 0; j < newN; j++) {
                row.push(this.basis1D(j, targetP, newU, xi));
            }
            A.push(row);
        }

        const ATA = Array(newN).fill(0).map(() => Array(newN).fill(0));
        const ATBx = Array(newN).fill(0);
        const ATBy = Array(newN).fill(0);
        const ATBz = Array(newN).fill(0);
        const ATBw = Array(newN).fill(0);

        for (let i = 0; i < newN; i++) {
            for (let j = 0; j < newN; j++) {
                let sum = 0;
                for (let k = 0; k <= M; k++) sum += A[k][i] * A[k][j];
                ATA[i][j] = sum;
            }
            let sx = 0, sy = 0, sz = 0, sw = 0;
            for (let k = 0; k <= M; k++) {
                sx += A[k][i] * samples[k].x;
                sy += A[k][i] * samples[k].y;
                sz += A[k][i] * samples[k].z;
                sw += A[k][i] * samples[k].w;
            }
            ATBx[i] = sx; ATBy[i] = sy; ATBz[i] = sz; ATBw[i] = sw;
        }

        const solve = (matA, matB) => {
            const n = matA.length;
            const mat = matA.map(r => [...r]);
            const B = [...matB];
            for (let i = 0; i < n; i++) {
                let maxEl = Math.abs(mat[i][i]), maxRow = i;
                for (let k = i + 1; k < n; k++) {
                    if (Math.abs(mat[k][i]) > maxEl) { maxEl = Math.abs(mat[k][i]); maxRow = k; }
                }
                [mat[maxRow], mat[i]] = [mat[i], mat[maxRow]];
                [B[maxRow], B[i]] = [B[i], B[maxRow]];
                for (let k = i + 1; k < n; k++) {
                    const c = -mat[k][i] / (mat[i][i] || 1e-12);
                    for (let j = i; j < n; j++) mat[k][j] += c * mat[i][j];
                    B[k] += c * B[i];
                }
            }
            const x = new Array(n).fill(0);
            for (let i = n - 1; i >= 0; i--) {
                let sum = 0;
                for (let k = i + 1; k < n; k++) sum += mat[i][k] * x[k];
                x[i] = (B[i] - sum) / (mat[i][i] || 1e-12);
            }
            return x;
        };

        const resX = solve(ATA, ATBx);
        const resY = solve(ATA, ATBy);
        const resZ = solve(ATA, ATBz);
        const resW = solve(ATA, ATBw);

        const newCP = [];
        const newWeights = [];
        for (let i = 0; i < newN; i++) {
            const w = resW[i] || 1.0;
            newWeights.push(w);
            newCP.push({ x: resX[i] / w, y: resY[i] / w, z: resZ[i] / w });
        }

        // Force endpoints for C0 continuity
        newCP[0] = { ...controlPoints[0] };
        newWeights[0] = weights[0];
        newCP[newN - 1] = { ...controlPoints[controlPoints.length - 1] };
        newWeights[newN - 1] = weights[weights.length - 1];

        return { controlPoints: newCP, weights: newWeights, U: newU };
    }

    /**
     * Reduce knot vector multiplicity for degree reduction
     */
    reduceKnotVector(U) {
        const unique = [...new Set(U)];
        const newU = [];
        unique.forEach(u => {
            const count = U.filter(k => k === u).length;
            // Standard reduction: decrement multiplicity if > 1
            const newCount = count > 1 ? count - 1 : 1;
            for (let i = 0; i < newCount; i++) newU.push(u);
        });
        return newU.sort((a,b) => a-b);
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
     * Exact Degree Elevation (Mathematical equivalent to Piegl-Tiller workflow)
     * Handles 4D Homogeneous space to avoid geometric distortion of rational surfaces.
     */
    elevate1D(curve, newP, newU) {
        // Fast-Path: Direct Degree Elevation (DDE) for linear clamped NURBS (p=1 -> p=2)
        if (curve.p === 1 && newP === 2) {
            return this.applyDDE(curve);
        }

        // General Case: 4D geometric fitting to exactly solve knot removal constraints
        return this.applyExactElevation4D(curve, newP, newU);
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
            newCP.push({ ...controlPoints[i] });
            newW.push(weights[i]);

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
     * Exact 4D Geometrical Fitting for Degree Elevation
     */
    applyExactElevation4D(curve, newP, newU) {
        const { p, U, controlPoints, weights } = curve;
        const M = Math.max(100, newU.length * 2);
        const samples = [];

        for (let i = 0; i <= M; i++) {
            const xi = i / M;
            let pos = { x: 0, y: 0, z: 0, w: 0 };
            for (let j = 0; j < controlPoints.length; j++) {
                const N = this.basis1D(j, p, U, xi);
                const w = weights[j];
                pos.x += N * controlPoints[j].x * w;
                pos.y += N * controlPoints[j].y * w;
                pos.z += N * controlPoints[j].z * w;
                pos.w += N * w;
            }
            samples.push(pos);
        }

        const newN = newU.length - newP - 1;
        const A = [];
        for (let i = 0; i <= M; i++) {
            const xi = i / M;
            const row = [];
            for (let j = 0; j < newN; j++) {
                row.push(this.basis1D(j, newP, newU, xi));
            }
            A.push(row);
        }

        const ATA = Array(newN).fill(0).map(() => Array(newN).fill(0));
        const ATBx = Array(newN).fill(0);
        const ATBy = Array(newN).fill(0);
        const ATBz = Array(newN).fill(0);
        const ATBw = Array(newN).fill(0);

        for (let i = 0; i < newN; i++) {
            for (let j = 0; j < newN; j++) {
                let sum = 0;
                for (let k = 0; k <= M; k++) sum += A[k][i] * A[k][j];
                ATA[i][j] = sum;
            }
            let sx = 0, sy = 0, sz = 0, sw = 0;
            for (let k = 0; k <= M; k++) {
                sx += A[k][i] * samples[k].x;
                sy += A[k][i] * samples[k].y;
                sz += A[k][i] * samples[k].z;
                sw += A[k][i] * samples[k].w;
            }
            ATBx[i] = sx; ATBy[i] = sy; ATBz[i] = sz; ATBw[i] = sw;
        }

        const solve = (matA, matB) => {
            const n = matA.length;
            const mat = matA.map(r => [...r]);
            const B = [...matB];
            for (let i = 0; i < n; i++) {
                let maxEl = Math.abs(mat[i][i]), maxRow = i;
                for (let k = i + 1; k < n; k++) {
                    if (Math.abs(mat[k][i]) > maxEl) { maxEl = Math.abs(mat[k][i]); maxRow = k; }
                }
                [mat[maxRow], mat[i]] = [mat[i], mat[maxRow]];
                [B[maxRow], B[i]] = [B[i], B[maxRow]];
                for (let k = i + 1; k < n; k++) {
                    const c = -mat[k][i] / (mat[i][i] || 1e-12);
                    for (let j = i; j < n; j++) mat[k][j] += c * mat[i][j];
                    B[k] += c * B[i];
                }
            }
            const x = new Array(n).fill(0);
            for (let i = n - 1; i >= 0; i--) {
                let sum = 0;
                for (let k = i + 1; k < n; k++) sum += mat[i][k] * x[k];
                x[i] = (B[i] - sum) / (mat[i][i] || 1e-12);
            }
            return x;
        };

        const resX = solve(ATA, ATBx);
        const resY = solve(ATA, ATBy);
        const resZ = solve(ATA, ATBz);
        const resW = solve(ATA, ATBw);

        const newCP = [];
        const newWeights = [];
        for (let i = 0; i < newN; i++) {
            newWeights.push(resW[i]);
            newCP.push({ 
                x: resX[i] / (resW[i] || 1), 
                y: resY[i] / (resW[i] || 1), 
                z: resZ[i] / (resW[i] || 1) 
            });
        }

        newCP[0] = { ...controlPoints[0] };
        newWeights[0] = weights[0];
        newCP[newN - 1] = { ...controlPoints[controlPoints.length - 1] };
        newWeights[newN - 1] = weights[weights.length - 1];

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

if (typeof module !== 'undefined' && module.exports) {
    module.exports = { NURBS2D };
} else {
    window.NURBS2D = NURBS2D;
}


// Export for use in browser
if (typeof module !== 'undefined' && module.exports) {
    module.exports = NURBS2D;
} else {
    window.NURBS2D = NURBS2D;
}
