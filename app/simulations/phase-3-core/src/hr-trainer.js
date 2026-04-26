/**
 * Hyper-Reduction Training Algorithms (Phase 3.4)
 * DEIM, UDEIM, Gappy POD, ECSW, ECM, Collocation (LSPG)
 */
class HRTrainer {
    /**
     * NNLS Solver (Lawson-Hanson active set method)
     * Solves: min ||A*x - b||^2  s.t. x >= 0
     */
    static nnls(A, b, maxIter = 300) {
        const m = A.length, n = A[0].length;
        const x = new Float64Array(n);
        const passive = new Set();
        const AtA = Array.from({length:n}, ()=>new Float64Array(n));
        const Atb = new Float64Array(n);
        for (let i = 0; i < n; i++) {
            for (let j = 0; j < n; j++) {
                let s = 0; for (let k = 0; k < m; k++) s += A[k][i]*A[k][j];
                AtA[i][j] = s;
            }
            let s = 0; for (let k = 0; k < m; k++) s += A[k][i]*b[k];
            Atb[i] = s;
        }
        for (let it = 0; it < maxIter; it++) {
            // Compute gradient w = At*(b - A*x)
            const w = new Float64Array(n);
            for (let j = 0; j < n; j++) {
                w[j] = Atb[j];
                for (let k = 0; k < n; k++) w[j] -= AtA[j][k]*x[k];
            }
            // Find max gradient in active set
            let maxW = -Infinity, maxIdx = -1;
            for (let j = 0; j < n; j++) {
                if (!passive.has(j) && w[j] > maxW) { maxW = w[j]; maxIdx = j; }
            }
            if (maxW <= 1e-10 || maxIdx < 0) break;
            passive.add(maxIdx);
            // Inner loop
            for (let inner = 0; inner < maxIter; inner++) {
                const pIdx = [...passive];
                const pn = pIdx.length;
                const subA = Array.from({length:pn}, ()=>new Float64Array(pn));
                const subB = new Float64Array(pn);
                for (let i = 0; i < pn; i++) {
                    for (let j = 0; j < pn; j++) subA[i][j] = AtA[pIdx[i]][pIdx[j]];
                    subB[i] = Atb[pIdx[i]];
                }
                const z = HRTrainer._solve(subA, subB);
                let allPos = true;
                for (let i = 0; i < pn; i++) if (z[i] <= 0) { allPos = false; break; }
                if (allPos) {
                    for (let i = 0; i < pn; i++) x[pIdx[i]] = z[i];
                    break;
                }
                // Find alpha
                let alpha = Infinity;
                for (let i = 0; i < pn; i++) {
                    if (z[i] <= 0) {
                        const a = x[pIdx[i]] / (x[pIdx[i]] - z[i]);
                        if (a < alpha) alpha = a;
                    }
                }
                for (let i = 0; i < pn; i++) x[pIdx[i]] += alpha * (z[i] - x[pIdx[i]]);
                // Remove zero entries from passive
                for (let i = pn - 1; i >= 0; i--) {
                    if (Math.abs(x[pIdx[i]]) < 1e-12) { x[pIdx[i]] = 0; passive.delete(pIdx[i]); }
                }
            }
        }
        return x;
    }

    static _solve(A, b) {
        const n = b.length;
        const M = A.map(r => Array.from(r)), B = Array.from(b);
        for (let i = 0; i < n; i++) {
            let mx = i;
            for (let j = i+1; j < n; j++) if (Math.abs(M[j][i]) > Math.abs(M[mx][i])) mx = j;
            [M[i],M[mx]] = [M[mx],M[i]]; [B[i],B[mx]] = [B[mx],B[i]];
            if (Math.abs(M[i][i]) < 1e-20) M[i][i] = 1e-20;
            for (let j = i+1; j < n; j++) {
                const f = M[j][i]/M[i][i]; B[j] -= f*B[i];
                for (let k = i; k < n; k++) M[j][k] -= f*M[i][k];
            }
        }
        const x = new Float64Array(n);
        for (let i = n-1; i >= 0; i--) {
            let s = 0; for (let j = i+1; j < n; j++) s += M[i][j]*x[j];
            x[i] = (B[i]-s) / (M[i][i] || 1e-20);
        }
        return x;
    }

    /** POD on a set of vector snapshots via SVD */
    static podVectors(snapshots, k) {
        const { Matrix, SVD } = window.mlMatrix;
        const S = new Matrix(snapshots.map(s => Array.from(s))).transpose();
        const svd = new SVD(S, { computeLeftSingularVectors: true, computeRightSingularVectors: false });
        const U = svd.leftSingularVectors;
        const trunc = Math.min(k, U.columns);
        return {
            basis: U.subMatrix(0, U.rows-1, 0, trunc-1),
            sigmas: svd.diagonal.slice(0, trunc)
        };
    }

    // ═══════════════════════════════════════════════════════════
    //  1. DEIM
    // ═══════════════════════════════════════════════════════════
    static trainDEIM(forceSnapshots, m) {
        const pod = HRTrainer.podVectors(forceSnapshots, m);
        const U_f = pod.basis; // [N x m]
        const N = U_f.rows;
        const indices = [];
        let maxVal = -1, maxIdx = -1;
        for (let i = 0; i < N; i++) {
            const val = Math.abs(U_f.get(i, 0));
            if (val > maxVal) { maxVal = val; maxIdx = i; }
        }
        indices.push(maxIdx);
        // Greedy selection
        for (let l = 1; l < m; l++) {
            // Solve P^T U_f[:,:l] c = P^T u_l
            const ul = []; for (let i = 0; i < N; i++) ul.push(U_f.get(i, l));
            const Psub = Array.from({length: l}, () => new Float64Array(l));
            const rhs = new Float64Array(l);
            for (let i = 0; i < l; i++) {
                for (let j = 0; j < l; j++) Psub[i][j] = U_f.get(indices[i], j);
                rhs[i] = ul[indices[i]];
            }
            const c = HRTrainer._solve(Psub, rhs);
            // Residual r = u_l - U_f[:,:l]*c
            const r = new Float64Array(N);
            for (let i = 0; i < N; i++) {
                r[i] = ul[i];
                for (let j = 0; j < l; j++) r[i] -= U_f.get(i, j) * c[j];
            }
            // Pick max |r| not already selected
            let maxR = -1, maxIdx = 0;
            for (let i = 0; i < N; i++) {
                if (!indices.includes(i) && Math.abs(r[i]) > maxR) { maxR = Math.abs(r[i]); maxIdx = i; }
            }
            indices.push(maxIdx);
        }
        return { U_f, indices, m, type: 'deim' };
    }

    /** Online DEIM reconstruction */
    static applyDEIM(F_partial, op) {
        const { U_f, indices, m } = op;
        const { Matrix } = window.mlMatrix;
        const N = U_f.rows;
        // Build P^T U_f
        const PtU = Array.from({length: m}, () => new Float64Array(m));
        for (let i = 0; i < m; i++)
            for (let j = 0; j < m; j++)
                PtU[i][j] = U_f.get(indices[i], j);
        const c = HRTrainer._solve(PtU, Array.from(F_partial));
        const result = new Float64Array(N);
        for (let i = 0; i < N; i++)
            for (let j = 0; j < m; j++)
                result[i] += U_f.get(i, j) * c[j];
        return result;
    }

    // ═══════════════════════════════════════════════════════════
    //  2. UDEIM (Unassembled DEIM)
    // ═══════════════════════════════════════════════════════════
    static trainUDEIM(elemForceSnapshots, nElems, nDofs, m) {
        // Stack unassembled: each snapshot -> [f_e1; f_e2; ...; f_eNel] of size nElems*nDofs
        const stacked = elemForceSnapshots.map(snap => {
            const v = new Float64Array(nElems * nDofs);
            for (let e = 0; e < nElems; e++)
                for (let d = 0; d < nDofs; d++)
                    v[e*nDofs + d] = snap[e][d];
            return v;
        });
        const pod = HRTrainer.podVectors(stacked, m);
        const U_e = pod.basis;
        const N = U_e.rows;
        const indices = [];
        let maxVal = -1, maxIdx = -1;
        for (let i = 0; i < N; i++) {
            const val = Math.abs(U_e.get(i, 0));
            if (val > maxVal) { maxVal = val; maxIdx = i; }
        }
        indices.push(maxIdx);
        for (let l = 1; l < m; l++) {
            const ul = []; for (let i = 0; i < N; i++) ul.push(U_e.get(i, l));
            const Psub = Array.from({length: l}, () => new Float64Array(l));
            const rhs = new Float64Array(l);
            for (let i = 0; i < l; i++) {
                for (let j = 0; j < l; j++) Psub[i][j] = U_e.get(indices[i], j);
                rhs[i] = ul[indices[i]];
            }
            const c = HRTrainer._solve(Psub, rhs);
            const r = new Float64Array(N);
            for (let i = 0; i < N; i++) {
                r[i] = ul[i]; for (let j = 0; j < l; j++) r[i] -= U_e.get(i, j)*c[j];
            }
            let maxR = -1, maxIdx = 0;
            for (let i = 0; i < N; i++)
                if (!indices.includes(i) && Math.abs(r[i]) > maxR) { maxR = Math.abs(r[i]); maxIdx = i; }
            indices.push(maxIdx);
        }
        // Map indices -> element IDs
        const elemIndices = [...new Set(indices.map(i => Math.floor(i / nDofs)))];
        return { U_e, indices, elemIndices, nElems, nDofs, m, type: 'udeim' };
    }

    static applyUDEIM(elemForces, op) {
        const { U_e, indices, nElems, nDofs, m } = op;
        const N = U_e.rows;
        // Stack partial
        const partial = new Float64Array(m);
        for (let i = 0; i < m; i++) {
            const eIdx = Math.floor(indices[i] / nDofs);
            const dIdx = indices[i] % nDofs;
            partial[i] = elemForces[eIdx][dIdx];
        }
        const PtU = Array.from({length:m}, ()=>new Float64Array(m));
        for (let i = 0; i < m; i++)
            for (let j = 0; j < m; j++) PtU[i][j] = U_e.get(indices[i], j);
        const c = HRTrainer._solve(PtU, Array.from(partial));
        // Reconstruct and assemble
        const assembled = new Float64Array(nDofs);
        for (let e = 0; e < nElems; e++)
            for (let d = 0; d < nDofs; d++) {
                let v = 0;
                for (let j = 0; j < m; j++) v += U_e.get(e*nDofs+d, j)*c[j];
                assembled[d] += v;
            }
        return assembled;
    }

    // ═══════════════════════════════════════════════════════════
    //  3. Gappy POD
    // ═══════════════════════════════════════════════════════════
    static trainGappyPOD(forceSnapshots, m) {
        const mOver = Math.min(Math.floor(m * 1.5), forceSnapshots[0].length);
        const pod = HRTrainer.podVectors(forceSnapshots, m);
        const U_f = pod.basis;
        const N = U_f.rows;
        // Select mOver indices via DEIM-like greedy
        const deimOp = HRTrainer.trainDEIM(forceSnapshots, mOver);
        return { U_f, indices: deimOp.indices.slice(0, mOver), m, mOver, type: 'gappy' };
    }

    static applyGappyPOD(F_partial, op) {
        const { U_f, indices, m, mOver } = op;
        const { Matrix, SVD } = window.mlMatrix;
        const N = U_f.rows;
        // P^T U_f [mOver x m]
        const PtU = new Matrix(mOver, m);
        for (let i = 0; i < mOver; i++)
            for (let j = 0; j < m; j++)
                PtU.set(i, j, U_f.get(indices[i], j));
        // Pseudo-inverse via SVD
        const svd = new SVD(PtU);
        const Fp = new Matrix([[...F_partial]]).transpose();
        // c = pinv(PtU) * F_partial
        const S = svd.diagonal;
        const Ut = svd.leftSingularVectors.transpose();
        const V = svd.rightSingularVectors;
        const Sinv = Array.from({length: Math.min(m, mOver)}, (_, i) =>
            Math.abs(S[i]) > 1e-12 ? 1/S[i] : 0);
        const temp = new Float64Array(m);
        for (let i = 0; i < Math.min(m, mOver); i++) {
            let dot = 0;
            for (let j = 0; j < mOver; j++) dot += Ut.get(i, j) * F_partial[j];
            for (let j = 0; j < m; j++) temp[j] += V.get(j, i) * Sinv[i] * dot;
        }
        const result = new Float64Array(N);
        for (let i = 0; i < N; i++)
            for (let j = 0; j < m; j++) result[i] += U_f.get(i, j) * temp[j];
        return result;
    }

    // ═══════════════════════════════════════════════════════════
    //  4. ECSW
    // ═══════════════════════════════════════════════════════════
    static trainECSW(Phi, elemForceSnapshots, nElems) {
        const Ns = elemForceSnapshots.length;
        const r = Phi.columns, nDofs = Phi.rows;
        // Build G [Ns*r x nElems]: projected element contributions
        const nRows = Ns * r;
        const G = Array.from({length: nRows}, () => new Float64Array(nElems));
        const bVec = new Float64Array(nRows);
        for (let s = 0; s < Ns; s++) {
            for (let e = 0; e < nElems; e++) {
                for (let i = 0; i < r; i++) {
                    let dot = 0;
                    for (let d = 0; d < nDofs; d++) dot += Phi.get(d, i) * elemForceSnapshots[s][e][d];
                    G[s*r + i][e] = dot;
                }
            }
            // b = Phi^T * sum_e f_e = Phi^T * F_int
            const Ftot = new Float64Array(nDofs);
            for (let e = 0; e < nElems; e++)
                for (let d = 0; d < nDofs; d++) Ftot[d] += elemForceSnapshots[s][e][d];
            for (let i = 0; i < r; i++) {
                let dot = 0;
                for (let d = 0; d < nDofs; d++) dot += Phi.get(d, i) * Ftot[d];
                bVec[s*r + i] = dot;
            }
        }
        const w = HRTrainer.nnls(G, bVec);
        const elements = [], weights = [];
        for (let e = 0; e < nElems; e++) {
            if (w[e] > 1e-10) { elements.push(e); weights.push(w[e]); }
        }
        return { elements, weights, type: 'ecsw' };
    }

    // ═══════════════════════════════════════════════════════════
    //  5. ECM (Empirical Cubature Method)
    // ═══════════════════════════════════════════════════════════
    static trainECM(Phi, elemForceSnapshots, nElems, tol = 1e-6) {
        const Ns = elemForceSnapshots.length;
        const r = Phi.columns, nDofs = Phi.rows;
        // Build contribution matrix C [r*Ns x nElems]
        const nRows = Ns * r;
        const C = Array.from({length: nRows}, () => new Float64Array(nElems));
        const target = new Float64Array(nRows);
        for (let s = 0; s < Ns; s++) {
            const Ftot = new Float64Array(nDofs);
            for (let e = 0; e < nElems; e++) {
                for (let d = 0; d < nDofs; d++) Ftot[d] += elemForceSnapshots[s][e][d];
                for (let i = 0; i < r; i++) {
                    let dot = 0;
                    for (let d = 0; d < nDofs; d++) dot += Phi.get(d, i) * elemForceSnapshots[s][e][d];
                    C[s*r + i][e] = dot;
                }
            }
            for (let i = 0; i < r; i++) {
                let dot = 0;
                for (let d = 0; d < nDofs; d++) dot += Phi.get(d, i) * Ftot[d];
                target[s*r + i] = dot;
            }
        }
        // Greedy column selection
        const selected = [];
        const residual = new Float64Array(target);
        let resNorm = Math.sqrt(residual.reduce((s,v)=>s+v*v, 0));
        const initNorm = resNorm;
        const maxPts = Math.min(nElems, r * Ns);
        const weights = [];
        for (let iter = 0; iter < maxPts && resNorm/initNorm > tol; iter++) {
            let bestCol = -1, bestCorr = -1;
            for (let e = 0; e < nElems; e++) {
                if (selected.includes(e)) continue;
                let corr = 0;
                for (let i = 0; i < nRows; i++) corr += C[i][e] * residual[i];
                if (Math.abs(corr) > bestCorr) { bestCorr = Math.abs(corr); bestCol = e; }
            }
            if (bestCol < 0) break;
            selected.push(bestCol);
            // Solve LSQ for weights
            const subG = Array.from({length: nRows}, (_, i) =>
                selected.map(e => C[i][e]));
            const w = HRTrainer.nnls(subG, Array.from(target));
            weights.length = 0;
            for (let i = 0; i < selected.length; i++) weights.push(w[i]);
            // Update residual
            for (let i = 0; i < nRows; i++) {
                residual[i] = target[i];
                for (let j = 0; j < selected.length; j++)
                    residual[i] -= C[i][selected[j]] * weights[j];
            }
            resNorm = Math.sqrt(residual.reduce((s,v)=>s+v*v, 0));
        }
        return { elements: selected, weights, type: 'ecm' };
    }

    // ═══════════════════════════════════════════════════════════
    //  6. Collocation (LSPG)
    // ═══════════════════════════════════════════════════════════
    static trainCollocation(forceSnapshots, Phi, m) {
        // Select collocation DOFs via DEIM on residual snapshots
        const deimOp = HRTrainer.trainDEIM(forceSnapshots, m);
        return { indices: deimOp.indices, m, type: 'collocation' };
    }
}

window.HRTrainer = HRTrainer;
