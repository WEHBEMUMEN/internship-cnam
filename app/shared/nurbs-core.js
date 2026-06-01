/**
 * B-Spline and NURBS Mathematical Core Library
 * Contains analytical NURBS algorithms: Cox-de Boor evaluation, span finding,
 * surface mapping, and analytical Jacobian calculations.
 */

class NURBSCore {
    /**
     * Determine the active knot span index for a given parametric coordinate.
     * @param {number} n Number of basis functions - 1
     * @param {number} p Degree of the curve
     * @param {number} u Parametric coordinate
     * @param {Array<number>} U Knot vector
     * @returns {number} Knot span index
     */
    static findSpan(n, p, u, U) {
        if (u >= U[n + 1]) return n;
        if (u <= U[p]) return p;
        let low = p;
        let high = n + 1;
        let mid = Math.floor((low + high) / 2);
        while (u < U[mid] || u >= U[mid + 1]) {
            if (u < U[mid]) {
                high = mid;
            } else {
                low = mid;
            }
            mid = Math.floor((low + high) / 2);
        }
        return mid;
    }

    /**
     * Compute non-vanishing B-spline basis functions (Cox-de Boor recursive optimization).
     * @param {number} i Span index
     * @param {number} u Parametric coordinate
     * @param {number} p Degree
     * @param {Array<number>} U Knot vector
     * @returns {Float64Array} Basis function evaluations
     */
    static basisFuns(i, u, p, U) {
        const N = new Float64Array(p + 1);
        N[0] = 1.0;
        const left = new Float64Array(p + 1);
        const right = new Float64Array(p + 1);
        for (let j = 1; j <= p; j++) {
            left[j] = u - U[i + 1 - j];
            right[j] = U[i + j] - u;
            let saved = 0.0;
            for (let r = 0; r < j; r++) {
                const temp = N[r] / (right[r + 1] + left[j - r] || 1e-12);
                N[r] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
            N[j] = saved;
        }
        return N;
    }

    /**
     * Compute non-vanishing B-spline basis functions and their derivatives up to order n.
     * @param {number} i Span index
     * @param {number} u Parametric coordinate
     * @param {number} p Degree
     * @param {Array<number>} U Knot vector
     * @param {number} n Derivative order (usually 1)
     * @returns {Array<Float64Array>} Derivatives of basis functions
     */
    static basisFunsDerivs(i, u, p, U, n = 1) {
        const ndu = Array(p + 1).fill(0).map(() => new Float64Array(p + 1));
        ndu[0][0] = 1.0;
        const left = new Float64Array(p + 1);
        const right = new Float64Array(p + 1);
        for (let j = 1; j <= p; j++) {
            left[j] = u - U[i + 1 - j];
            right[j] = U[i + j] - u;
            let saved = 0.0;
            for (let r = 0; r < j; r++) {
                ndu[j][r] = right[r + 1] + left[j - r];
                const temp = ndu[r][j - 1] / (ndu[j][r] || 1e-12);
                ndu[r][j] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
            ndu[j][j] = saved;
        }

        const ders = Array(n + 1).fill(0).map(() => new Float64Array(p + 1));
        for (let j = 0; j <= p; j++) {
            ders[0][j] = ndu[j][p];
        }

        const a = Array(2).fill(0).map(() => new Float64Array(p + 1));
        for (let r = 0; r <= p; r++) {
            let s1 = 0, s2 = 1;
            a[0][0] = 1.0;
            for (let k = 1; k <= n; k++) {
                let d = 0.0;
                const rk = r - k;
                const pk = p - k;
                if (r >= k) {
                    a[s2][0] = a[s1][0] / (ndu[pk + 1][rk] || 1e-12);
                    d = a[s2][0] * ndu[rk][pk];
                }
                const j1 = (rk >= -1) ? 1 : -rk;
                const j2 = (r - 1 <= pk) ? k - 1 : p - r;
                for (let j = j1; j <= j2; j++) {
                    a[s2][j] = (a[s1][j] - a[s1][j - 1]) / (ndu[pk + 1][rk + j] || 1e-12);
                    d += a[s2][j] * ndu[rk + j][pk];
                }
                if (r <= pk) {
                    a[s2][k] = -a[s1][k - 1] / (ndu[pk + 1][r] || 1e-12);
                    d += a[s2][k] * ndu[r][pk];
                }
                ders[k][r] = d;
                const temp = s1; s1 = s2; s2 = temp;
            }
        }

        let rD = p;
        for (let k = 1; k <= n; k++) {
            for (let j = 0; j <= p; j++) {
                ders[k][j] *= rD;
            }
            rD *= (p - k);
        }
        return ders;
    }

    /**
     * Compute 1D basis function recursively (simple evaluation).
     */
    static basis1D(i, p, U, xi) {
        if (p === 0) {
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
     * Map parametric coordinates (xi, eta) to physical space (x, y, z) and return state.
     * @param {Object} patch NURBS patch definition { p, q, U, V, controlPoints, weights }
     * @param {number} xi Parametric coordinate u
     * @param {number} eta Parametric coordinate v
     * @returns {Object} { position: {x,y,z}, denominator }
     */
    static getSurfaceState(patch, xi, eta) {
        const { controlPoints, weights, p, q, U, V } = patch;
        const eps = 1e-12;
        xi = Math.max(0, Math.min(1 - eps, xi));
        eta = Math.max(0, Math.min(1 - eps, eta));

        const nU = controlPoints.length;
        const nV = controlPoints[0].length;

        const spanU = this.findSpan(nU - 1, p, xi, U);
        const spanV = this.findSpan(nV - 1, q, eta, V);

        const N = this.basisFuns(spanU, xi, p, U);
        const M = this.basisFuns(spanV, eta, q, V);

        const position = { x: 0, y: 0, z: 0 };
        let denominator = 0;

        for (let i = 0; i <= p; i++) {
            for (let j = 0; j <= q; j++) {
                const w = weights[spanU - p + i][spanV - q + j];
                const basis = N[i] * M[j] * w;
                denominator += basis;
                
                const cp = controlPoints[spanU - p + i][spanV - q + j];
                position.x += basis * cp.x;
                position.y += basis * cp.y;
                position.z += basis * cp.z;
            }
        }

        if (denominator > 0) {
            position.x /= denominator;
            position.y /= denominator;
            position.z /= denominator;
        }

        return { position, denominator };
    }

    /**
     * Analytical Surface Derivatives (Quotient Rule).
     */
    static getSurfaceDerivatives(patch, xi, eta) {
        const { controlPoints, weights, p, q, U, V } = patch;
        const eps = 1e-12;
        xi = Math.max(0, Math.min(1 - eps, xi));
        eta = Math.max(0, Math.min(1 - eps, eta));

        const nU = controlPoints.length;
        const nV = controlPoints[0].length;

        const spanU = this.findSpan(nU - 1, p, xi, U);
        const spanV = this.findSpan(nV - 1, q, eta, V);

        const dersU = this.basisFunsDerivs(spanU, xi, p, U, 1);
        const dersV = this.basisFunsDerivs(spanV, eta, q, V, 1);

        const A = { x: 0, y: 0, z: 0 };
        const Au = { x: 0, y: 0, z: 0 };
        const Av = { x: 0, y: 0, z: 0 };
        let W = 0;
        let Wu = 0;
        let Wv = 0;

        for (let i = 0; i <= p; i++) {
            const Ni = dersU[0][i];
            const dNi = dersU[1][i];
            for (let j = 0; j <= q; j++) {
                const Mj = dersV[0][j];
                const dMj = dersV[1][j];

                const w = weights[spanU - p + i][spanV - q + j];
                const cp = controlPoints[spanU - p + i][spanV - q + j];

                const basis = Ni * Mj * w;
                const basisU = dNi * Mj * w;
                const basisV = Ni * dMj * w;

                A.x += basis * cp.x; A.y += basis * cp.y; A.z += basis * cp.z;
                Au.x += basisU * cp.x; Au.y += basisU * cp.y; Au.z += basisU * cp.z;
                Av.x += basisV * cp.x; Av.y += basisV * cp.y; Av.z += basisV * cp.z;

                W += basis;
                Wu += basisU;
                Wv += basisV;
            }
        }

        if (W < 1e-12) W = 1e-12;

        const pos = { x: A.x / W, y: A.y / W, z: A.z / W };
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
     * Analytical cross product of derivatives to determine Jacobian area determinant.
     */
    static getJacobianDeterminant(patch, xi, eta) {
        const deriv = this.getSurfaceDerivatives(patch, xi, eta);
        const tu = deriv.dU;
        const tv = deriv.dV;

        const nx = tu.y * tv.z - tu.z * tv.y;
        const ny = tu.z * tv.x - tu.x * tv.z;
        const nz = tu.x * tv.y - tu.y * tv.x;

        let area = Math.sqrt(nx * nx + ny * ny + nz * nz);
        if (isNaN(area) || area < 1e-12) area = 1e-12;
        return area;
    }
}

// Export for environment compatibility
if (typeof module !== 'undefined' && module.exports) {
    module.exports = NURBSCore;
} else {
    window.NURBSCore = NURBSCore;
}
