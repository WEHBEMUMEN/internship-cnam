/**
 * IGA Core Engine — Phase 1.8: Circle Geometry + Refinement
 * Implements Cox-de Boor basis functions, NURBS curve evaluation,
 * and p-refinement, h-refinement, and k-refinement algorithms.
 */

export class BasisFunctions {
    /**
     * Cox-de Boor recursion: evaluate N_{i,p}(xi)
     */
    static evaluate(i, p, xi, U, memo = null) {
        if (p === 0) {
            if (xi >= U[i] && xi < U[i + 1]) return 1.0;
            if (xi === 1.0 && U[i + 1] === 1.0 && U[i] < 1.0) return 1.0;
            return 0.0;
        }

        const isRootCall = memo === null;
        if (isRootCall) memo = new Map();

        const key = `${i}_${p}`;
        if (memo.has(key)) return memo.get(key);

        let val = 0;
        const d1 = U[i + p] - U[i];
        if (d1 > 0) val += ((xi - U[i]) / d1) * this.evaluate(i, p - 1, xi, U, memo);

        const d2 = U[i + p + 1] - U[i + 1];
        if (d2 > 0) val += ((U[i + p + 1] - xi) / d2) * this.evaluate(i + 1, p - 1, xi, U, memo);

        memo.set(key, val);
        return val;
    }

    static evaluateAll(n, p, xi, U) {
        const results = new Array(n).fill(0);
        const memo = new Map();
        for (let i = 0; i < n; i++) results[i] = this.evaluate(i, p, xi, U, memo);
        return results;
    }
}

export class NURBS {
    /**
     * Evaluate all rational basis functions R_i(xi)
     */
    static evaluateAll(n, p, xi, U, weights) {
        const N = BasisFunctions.evaluateAll(n, p, xi, U);
        const R = new Array(n).fill(0);

        let sum = 0;
        for (let i = 0; i < n; i++) sum += N[i] * weights[i];

        if (sum > 0) {
            for (let i = 0; i < n; i++) R[i] = (N[i] * weights[i]) / sum;
        }
        return R;
    }
}

export class Curve {
    /**
     * Evaluate the curve point C(xi)
     */
    static evaluate(xi, p, U, points, weights) {
        const n = points.length;
        const R = NURBS.evaluateAll(n, p, xi, U, weights);

        let x = 0, y = 0;
        for (let i = 0; i < n; i++) {
            x += R[i] * points[i].x;
            y += R[i] * points[i].y;
        }
        return { x, y };
    }
}

/**
 * Standard NURBS circle definition (9 CPs, degree 2)
 * Uses the classic rational quadratic construction with w = cos(pi/4)
 */
export function createCircle(cx = 0, cy = 0, r = 1) {
    const w = Math.cos(Math.PI / 4); // = sqrt(2)/2 ≈ 0.7071

    const points = [
        { x: cx + r, y: cy },       // 0: right
        { x: cx + r, y: cy + r },   // 1: top-right (weighted)
        { x: cx,     y: cy + r },   // 2: top
        { x: cx - r, y: cy + r },   // 3: top-left (weighted)
        { x: cx - r, y: cy },       // 4: left
        { x: cx - r, y: cy - r },   // 5: bottom-left (weighted)
        { x: cx,     y: cy - r },   // 6: bottom
        { x: cx + r, y: cy - r },   // 7: bottom-right (weighted)
        { x: cx + r, y: cy },       // 8: right (closes the circle)
    ];

    const weights = [1, w, 1, w, 1, w, 1, w, 1];

    // Open knot vector for degree 2, 9 CPs: m = n + p + 1 = 9 + 2 + 1 = 12
    const knots = [0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1];

    return { degree: 2, points, weights, knots };
}

/**
 * Knot Insertion (h-refinement)
 * Inserts a single knot 'ubar' into the NURBS curve (Boehm's algorithm).
 * Works in homogeneous (projective) coordinates to handle rational curves correctly.
 */
export function insertKnot(degree, knots, points, weights, ubar) {
    const p = degree;
    const n = points.length;

    // Find span: k such that U[k] <= ubar < U[k+1]
    let k = -1;
    for (let i = 0; i < knots.length - 1; i++) {
        if (ubar >= knots[i] && ubar < knots[i + 1]) { k = i; break; }
    }
    if (k === -1) k = knots.length - p - 2; // clamp if at end

    // New knot vector
    const newKnots = [...knots.slice(0, k + 1), ubar, ...knots.slice(k + 1)];

    // Convert to homogeneous coords
    const Pw = points.map((pt, i) => ({ x: pt.x * weights[i], y: pt.y * weights[i], w: weights[i] }));

    // New control points in homogeneous space
    const newPw = [];
    for (let i = 0; i <= n; i++) {
        if (i <= k - p) {
            newPw.push(Pw[i]);
        } else if (i >= k + 1) {
            newPw.push(Pw[i - 1]);
        } else {
            // Blend
            const alpha = (ubar - knots[i]) / (knots[i + p] - knots[i]);
            newPw.push({
                x: (1 - alpha) * Pw[i - 1].x + alpha * Pw[i].x,
                y: (1 - alpha) * Pw[i - 1].y + alpha * Pw[i].y,
                w: (1 - alpha) * Pw[i - 1].w + alpha * Pw[i].w
            });
        }
    }

    // Back to Cartesian
    const newPoints = newPw.map(pw => ({ x: pw.x / pw.w, y: pw.y / pw.w }));
    const newWeights = newPw.map(pw => pw.w);

    return { degree: p, knots: newKnots, points: newPoints, weights: newWeights };
}

/**
 * h-refinement: insert multiple knots at midpoints of each unique knot span.
 */
export function hRefine(degree, knots, points, weights) {
    // Find unique knot spans
    const midpoints = [];
    for (let i = 0; i < knots.length - 1; i++) {
        if (knots[i + 1] - knots[i] > 1e-10) {
            const mid = (knots[i] + knots[i + 1]) / 2;
            // Avoid inserting duplicates
            if (!knots.includes(mid) && !midpoints.includes(mid)) {
                midpoints.push(mid);
            }
        }
    }
    midpoints.sort((a, b) => a - b);

    let result = { degree, knots: [...knots], points: [...points], weights: [...weights] };
    for (const u of midpoints) {
        result = insertKnot(result.degree, result.knots, result.points, result.weights, u);
    }
    return result;
}

/**
 * p-refinement (degree elevation)
 * Uses least-squares fitting in homogeneous coordinates to exactly preserve rational geometry.
 */
export function pRefine(degree, knots, points, weights) {
    const newP = degree + 1;
    const M = 200; // sampling density

    // Sample current curve in homogeneous coordinates
    const samplesW = [];
    for (let i = 0; i <= M; i++) {
        const xi = i / M;
        const N = BasisFunctions.evaluateAll(points.length, degree, xi, knots);
        let wx = 0, wy = 0, w = 0;
        for (let j = 0; j < points.length; j++) {
            const Pwj = points[j];
            const wj = weights[j];
            wx += N[j] * Pwj.x * wj;
            wy += N[j] * Pwj.y * wj;
            w  += N[j] * wj;
        }
        samplesW.push({ wx, wy, w });
    }

    // Build new knot vector for elevated degree
    const uniqueKnots = [];
    const multiplicities = [];
    for (let i = 0; i < knots.length; i++) {
        if (i === 0 || Math.abs(knots[i] - knots[i - 1]) > 1e-10) {
            uniqueKnots.push(knots[i]);
            multiplicities.push(1);
        } else {
            multiplicities[multiplicities.length - 1]++;
        }
    }

    // Elevate: each multiplicity increases by 1, capped at newP+1 for endpoints
    const newKnotsArr = [];
    for (let i = 0; i < uniqueKnots.length; i++) {
        const newMult = Math.min(multiplicities[i] + 1, newP + 1);
        for (let j = 0; j < newMult; j++) newKnotsArr.push(uniqueKnots[i]);
    }

    const newN = newKnotsArr.length - newP - 1;

    // Least-squares fit: A * Pw_new = samplesW
    const A = [];
    for (let i = 0; i <= M; i++) {
        const xi = i / M;
        A.push(BasisFunctions.evaluateAll(newN, newP, xi, newKnotsArr));
    }

    // ATA and ATB
    const ATA = Array.from({ length: newN }, () => new Array(newN).fill(0));
    const ATBx = new Array(newN).fill(0);
    const ATBy = new Array(newN).fill(0);
    const ATBw = new Array(newN).fill(0);

    for (let i = 0; i < newN; i++) {
        for (let j = 0; j < newN; j++) {
            let s = 0;
            for (let k = 0; k <= M; k++) s += A[k][i] * A[k][j];
            ATA[i][j] = s;
        }
        let sx = 0, sy = 0, sw = 0;
        for (let k = 0; k <= M; k++) {
            sx += A[k][i] * samplesW[k].wx;
            sy += A[k][i] * samplesW[k].wy;
            sw += A[k][i] * samplesW[k].w;
        }
        ATBx[i] = sx;
        ATBy[i] = sy;
        ATBw[i] = sw;
    }

    const solveGauss = (matA, matB) => {
        const n = matA.length;
        const A = matA.map(r => [...r]);
        const B = [...matB];
        for (let i = 0; i < n; i++) {
            let maxEl = Math.abs(A[i][i]), maxRow = i;
            for (let k = i + 1; k < n; k++) {
                if (Math.abs(A[k][i]) > maxEl) { maxEl = Math.abs(A[k][i]); maxRow = k; }
            }
            [A[i], A[maxRow]] = [A[maxRow], A[i]];
            [B[i], B[maxRow]] = [B[maxRow], B[i]];
            for (let k = i + 1; k < n; k++) {
                if (A[i][i] === 0) continue;
                const c = -A[k][i] / A[i][i];
                for (let j = i; j < n; j++) A[k][j] += c * A[i][j];
                B[k] += c * B[i];
            }
        }
        const x = new Array(n).fill(0);
        for (let i = n - 1; i >= 0; i--) {
            if (Math.abs(A[i][i]) < 1e-14) continue;
            let s = 0;
            for (let k = i + 1; k < n; k++) s += A[i][k] * x[k];
            x[i] = (B[i] - s) / A[i][i];
        }
        return x;
    };

    const newWx = solveGauss(ATA, ATBx);
    const newWy = solveGauss(ATA, ATBy);
    const newW  = solveGauss(ATA, ATBw);

    const newPoints = [];
    const newWeights = [];
    for (let i = 0; i < newN; i++) {
        const wi = newW[i];
        newWeights.push(wi);
        newPoints.push({ x: newWx[i] / wi, y: newWy[i] / wi });
    }

    // Pin endpoints
    if (newPoints.length > 0 && points.length > 0) {
        newPoints[0] = { x: points[0].x, y: points[0].y };
        newWeights[0] = weights[0];
        newPoints[newN - 1] = { x: points[points.length - 1].x, y: points[points.length - 1].y };
        newWeights[newN - 1] = weights[weights.length - 1];
    }

    return { degree: newP, knots: newKnotsArr, points: newPoints, weights: newWeights };
}

/**
 * k-refinement: first elevate degree, then insert knots.
 * This yields higher continuity than h-then-p because knot insertion
 * after degree elevation preserves C^{p_new - 1} inter-element continuity.
 */
export function kRefine(degree, knots, points, weights, elevations = 1) {
    let result = { degree, knots, points, weights };
    for (let i = 0; i < elevations; i++) {
        result = pRefine(result.degree, result.knots, result.points, result.weights);
    }
    return hRefine(result.degree, result.knots, result.points, result.weights);
}
