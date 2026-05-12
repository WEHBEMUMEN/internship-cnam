/**
 * Phase 5.1a - POD Engine
 */

class PODEngine {
    static computeBasis(snapshots, k = 10) {
        console.log("[PODEngine] computeBasis called with", snapshots.length, "snapshots, k =", k);
        if (!snapshots || snapshots.length === 0) {
            console.warn("[PODEngine] No snapshots provided");
            return null;
        }

        const nDofs = snapshots[0].length;
        const nSnaps = snapshots.length;
        console.log("[PODEngine] Matrix dimensions:", nDofs, "x", nSnaps);

        const lib = window.mlMatrix || window.ML;
        const { Matrix, SingularValueDecomposition } = lib;

        // 1. Build Snapshot Matrix S
        const S = new Matrix(nDofs, nSnaps);
        for (let j = 0; j < nSnaps; j++) {
            for (let i = 0; i < nDofs; i++) {
                S.set(i, j, snapshots[j][i]);
            }
        }

        // 2. SVD
        const svd = new SingularValueDecomposition(S);
        const U = svd.leftSingularVectors;
        const s = svd.diagonal;

        // 3. Truncate
        const actualK = Math.min(k, nSnaps, nDofs);
        const Phi = U.subMatrix(0, nDofs - 1, 0, actualK - 1);

        // 4. Energy Spectrum
        const totalVar = s.reduce((a, b) => a + b*b, 0);
        let cumVar = 0;
        const energy = s.slice(0, actualK).map(si => {
            cumVar += si*si;
            return cumVar / totalVar;
        });

        return {
            Phi: Phi.to2DArray(),
            energy: energy,
            singularValues: s.slice(0, actualK)
        };
    }

    /**
     * Solves a small linear system for ROM amplitudes using Math.js
     */
    static solveLinear(K_arr, f_arr) {
        try {
            // Using Math.js for robust linear solve (already loaded in index.html)
            // It handles singular or near-singular systems better in the browser.
            const sol = window.math.lusolve(K_arr, f_arr);
            // math.js returns a column vector (matrix)
            return Array.isArray(sol) ? sol.flat() : sol.toArray().flat();
        } catch (e) {
            console.warn("[PODEngine] math.js lusolve failed:", e);
            // Extreme fallback: Simple diagonal inverse if everything else fails
            const n = K_arr.length;
            const ar = new Float64Array(n);
            for(let i=0; i<n; i++) {
                ar[i] = f_arr[i] / (K_arr[i][i] || 1.0);
            }
            return Array.from(ar);
        }
    }
}

window.PODEngine = PODEngine;
