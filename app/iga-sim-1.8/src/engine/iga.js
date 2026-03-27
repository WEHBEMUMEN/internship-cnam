/**
 * IGA Core Engine — Phase 1.8: Circle Geometry + Refinement
 * Implements Cox-de Boor basis functions, NURBS curve evaluation,
 * and p-refinement, h-refinement, and k-refinement algorithms.
 */

export class BasisFunctions {
    /**
     * Cox-de Boor recursion: evaluate N_{i,p}(xi)
     */
    static evaluate(i, p, xi, U) {
        if (p === 0) {
            if (xi >= U[i] && xi < U[i + 1]) return 1.0;
            if (xi === 1.0 && U[i + 1] === 1.0 && U[i] < 1.0) return 1.0;
            return 0.0;
        }

        let val = 0;
        const d1 = U[i + p] - U[i];
        if (d1 > 0) val += ((xi - U[i]) / d1) * this.evaluate(i, p - 1, xi, U);

        const d2 = U[i + p + 1] - U[i + 1];
        if (d2 > 0) val += ((U[i + p + 1] - xi) / d2) * this.evaluate(i + 1, p - 1, xi, U);

        return val;
    }

    static evaluateAll(n, p, xi, U) {
        const results = new Array(n).fill(0);
        for (let i = 0; i < n; i++) results[i] = this.evaluate(i, p, xi, U);
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
 * Uses least-squares fitting to refit the curve at the new degree.
 * This preserves the geometric shape by sampling and refitting.
 */
export function pRefine(degree, knots, points, weights) {
    const newP = degree + 1;
    const M = 200; // sampling density

    // Sample current curve
    const samples = [];
    for (let i = 0; i <= M; i++) {
        const xi = i / M;
        samples.push(Curve.evaluate(xi, degree, knots, points, weights));
    }

    // Build new knot vector for elevated degree
    // Increase multiplicity of each unique knot by 1
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
    const newWeights = new Array(newN).fill(1.0);

    // Least-squares fit: A * P_new = samples
    const A = [];
    for (let i = 0; i <= M; i++) {
        const xi = i / M;
        A.push(NURBS.evaluateAll(newN, newP, xi, newKnotsArr, newWeights));
    }

    // ATA and ATB
    const ATA = Array.from({ length: newN }, () => new Array(newN).fill(0));
    const ATBx = new Array(newN).fill(0);
    const ATBy = new Array(newN).fill(0);

    for (let i = 0; i < newN; i++) {
        for (let j = 0; j < newN; j++) {
            let s = 0;
            for (let k = 0; k <= M; k++) s += A[k][i] * A[k][j];
            ATA[i][j] = s;
        }
        let sx = 0, sy = 0;
        for (let k = 0; k <= M; k++) {
            sx += A[k][i] * samples[k].x;
            sy += A[k][i] * samples[k].y;
        }
        ATBx[i] = sx;
        ATBy[i] = sy;
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

    const newX = solveGauss(ATA, ATBx);
    const newY = solveGauss(ATA, ATBy);

    const newPoints = [];
    for (let i = 0; i < newN; i++) newPoints.push({ x: newX[i], y: newY[i] });

    // Pin endpoints
    if (newPoints.length > 0 && samples.length > 0) {
        newPoints[0] = { x: samples[0].x, y: samples[0].y };
        newPoints[newN - 1] = { x: samples[M].x, y: samples[M].y };
    }

    return { degree: newP, knots: newKnotsArr, points: newPoints, weights: newWeights };
}

/**
 * k-refinement: first elevate degree, then insert knots.
 * This yields higher continuity than h-then-p because knot insertion
 * after degree elevation preserves C^{p_new - 1} inter-element continuity.
 */
export function kRefine(degree, knots, points, weights) {
    // Step 1: p-refine
    const elevated = pRefine(degree, knots, points, weights);
    // Step 2: h-refine the elevated curve
    return hRefine(elevated.degree, elevated.knots, elevated.points, elevated.weights);
}
