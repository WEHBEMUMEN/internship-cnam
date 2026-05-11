/**
 * Phase 5.1a - POD Engine
 */

class PODEngine {
    static computeBasis(snapshots, k = 10) {
        if (!snapshots || snapshots.length === 0) return null;

        const nDofs = snapshots[0].length;
        const nSnaps = snapshots.length;

        const { Matrix, SingularValueDecomposition } = window.mlMatrix;

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
}

window.PODEngine = PODEngine;
