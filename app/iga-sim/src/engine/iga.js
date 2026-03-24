/**
 * IGA Core Engine - 1D B-Spline Basis Functions
 * Implements Cox-de Boor recursion formula
 */

export class KnotVector {
    constructor(n, p) {
        this.n = n; // Number of control points
        this.p = p; // Degree
        this.values = this.generateUniformOpen();
    }

    /**
     * Generates an open uniform knot vector
     * [0, 0, ..., 0, 1/(n-p), 2/(n-p), ..., 1, 1, ..., 1]
     */
    generateUniformOpen() {
        const m = this.n + this.p + 1;
        const knots = new Array(m);
        
        // p + 1 zeros at the beginning
        for (let i = 0; i <= this.p; i++) {
            knots[i] = 0;
        }
        
        // Uniform internal knots
        const internalCount = this.n - this.p - 1;
        for (let i = 1; i <= internalCount; i++) {
            knots[this.p + i] = i / (internalCount + 1);
        }
        
        // p + 1 ones at the end
        for (let i = m - this.p - 1; i < m; i++) {
            knots[i] = 1;
        }
        
        return knots;
    }

    update(n, p) {
        this.n = n;
        this.p = p;
        this.values = this.generateUniformOpen();
    }
}

export class BasisFunctions {
    /**
     * Evaluates the i-th basis function of degree p at parameter xi
     * @param {number} i - Index of basis function (0 to n-1)
     * @param {number} p - Degree
     * @param {number} xi - Parameter value
     * @param {Array} U - Knot vector
     */
    static evaluate(i, p, xi, U) {
        // Degree 0 (Step function)
        if (p === 0) {
            if (xi >= U[i] && xi < U[i+1]) return 1.0;
            // Special case for the end of the interval
            if (xi === 1.0 && U[i+1] === 1.0 && U[i] < 1.0) return 1.0;
            return 0.0;
        }

        // Degree p > 0 (Recursion)
        let val = 0;
        
        const denom1 = U[i+p] - U[i];
        if (denom1 > 0) {
            val += ((xi - U[i]) / denom1) * this.evaluate(i, p - 1, xi, U);
        }
        
        const denom2 = U[i+p+1] - U[i+1];
        if (denom2 > 0) {
            val += ((U[i+p+1] - xi) / denom2) * this.evaluate(i + 1, p - 1, xi, U);
        }
        
        return val;
    }

    /**
     * Compute all non-zero basis functions at parameter xi
     * More efficient than evaluating individually if we need them all
     */
    static evaluateAll(n, p, xi, U) {
        const results = new Array(n).fill(0);
        for (let i = 0; i < n; i++) {
            results[i] = this.evaluate(i, p, xi, U);
        }
        return results;
    }
}
