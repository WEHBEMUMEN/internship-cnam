/**
 * Hyper-Reduced ROM Solver (Phase 3.4)
 * Newton-Raphson in reduced space with hyper-reduced assembly.
 */
class HRSolver {
    constructor(solver, romEngine, hrElements) {
        this.fom = solver;
        this.rom = romEngine;
        this.hre = hrElements;
    }

    /**
     * Solve using element-sampling HR (ECSW or ECM)
     */
    solveElementSampling(patch, bcs, loads, hrOp, options = {}) {
        const { iterations = 15, tolerance = 1e-6, steps = 1 } = options;
        const Phi = this.rom.Phi;
        if (!Phi) throw new Error("POD basis not computed");
        const { Matrix } = window.mlMatrix;
        const k = Phi.columns, nDofs = Phi.rows;
        const nV = patch.controlPoints[0].length;
        let ur = new Float64Array(k);
        const residualHistory = [];
        const allElems = this.hre.getElements(patch);

        const F_ext_total = new Float64Array(nDofs);
        loads.forEach(l => {
            const idx = (l.i * nV + l.j) * 2;
            F_ext_total[idx] += l.fx;
            F_ext_total[idx + 1] += l.fy;
        });

        for (let s = 1; s <= steps; s++) {
            const F_ext = F_ext_total.map(f => f * (s/steps));
            for (let iter = 0; iter < iterations; iter++) {
                const ur_mat = new Matrix([Array.from(ur)]).transpose();
                const u_full = new Float64Array(Phi.mmul(ur_mat).to2DArray().map(r => r[0]));

                // HR assembly: only selected elements
                const F_int_proj = new Float64Array(k);
                const Kt_proj = Array.from({length:k}, () => new Float64Array(k));

                for (let idx = 0; idx < hrOp.elements.length; idx++) {
                    const eId = hrOp.elements[idx];
                    const w = hrOp.weights[idx];
                    const elem = allElems[eId];
                    const fe = this.hre.elementForce(patch, u_full, elem);
                    const ke = this.hre.elementTangent(patch, u_full, elem);

                    // Projected element force: Phi^T * fe * w
                    for (let i = 0; i < k; i++) {
                        let dot = 0;
                        for (let d = 0; d < nDofs; d++) dot += Phi.get(d, i) * fe[d];
                        F_int_proj[i] += dot * w;
                    }
                    // Projected element tangent: Phi^T * ke * Phi * w
                    for (let i = 0; i < k; i++) {
                        for (let j = 0; j < k; j++) {
                            let dot = 0;
                            for (let a = 0; a < nDofs; a++)
                                for (let b = 0; b < nDofs; b++)
                                    dot += Phi.get(a, i) * ke[a][b] * Phi.get(b, j);
                            Kt_proj[i][j] += dot * w;
                        }
                    }
                }

                // Projected external force
                const F_ext_proj = new Float64Array(k);
                for (let i = 0; i < k; i++) {
                    let dot = 0;
                    for (let d = 0; d < nDofs; d++) dot += Phi.get(d, i) * F_ext[d];
                    F_ext_proj[i] = dot;
                }

                const R_red = F_ext_proj.map((f, i) => f - F_int_proj[i]);
                let norm = Math.sqrt(R_red.reduce((s, v) => s + v*v, 0));
                residualHistory.push({step: s, iter, norm});
                if (norm < tolerance && iter > 0) break;

                const dur = this.fom.gaussianElimination(Kt_proj, Array.from(R_red));
                for (let i = 0; i < k; i++) ur[i] += dur[i];
            }
        }

        const final_ur = new Matrix([Array.from(ur)]).transpose();
        const u = new Float64Array(Phi.mmul(final_ur).to2DArray().map(r => r[0]));
        return { u, ur, residualHistory, sampledElements: hrOp.elements.length, totalElements: allElems.length };
    }

    /**
     * Solve using DEIM/Gappy/UDEIM
     */
    solveInterpolation(patch, bcs, loads, hrOp, method, options = {}) {
        const { iterations = 15, tolerance = 1e-6, steps = 1 } = options;
        const Phi = this.rom.Phi;
        if (!Phi) throw new Error("POD basis not computed");
        const { Matrix } = window.mlMatrix;
        const k = Phi.columns, nDofs = Phi.rows;
        const nV = patch.controlPoints[0].length;
        let ur = new Float64Array(k);
        const residualHistory = [];
        const allElems = this.hre.getElements(patch);

        const F_ext_total = new Float64Array(nDofs);
        loads.forEach(l => {
            const idx = (l.i * nV + l.j) * 2;
            F_ext_total[idx] += l.fx;
            F_ext_total[idx + 1] += l.fy;
        });

        for (let s = 1; s <= steps; s++) {
            const F_ext = F_ext_total.map(f => f * (s/steps));
            for (let iter = 0; iter < iterations; iter++) {
                const ur_mat = new Matrix([Array.from(ur)]).transpose();
                const u_full = new Float64Array(Phi.mmul(ur_mat).to2DArray().map(r => r[0]));

                let F_int_approx;
                if (method === 'deim') {
                    const F_int_full = this.fom.calculateInternalForce(patch, u_full);
                    this.fom.applyPenaltyConstraints(null, F_int_full, u_full, patch);
                    const partial = new Float64Array(hrOp.indices.length);
                    for (let i = 0; i < hrOp.indices.length; i++) partial[i] = F_int_full[hrOp.indices[i]];
                    F_int_approx = HRTrainer.applyDEIM(partial, hrOp);
                } else if (method === 'gappy') {
                    const F_int_full = this.fom.calculateInternalForce(patch, u_full);
                    this.fom.applyPenaltyConstraints(null, F_int_full, u_full, patch);
                    const partial = new Float64Array(hrOp.indices.length);
                    for (let i = 0; i < hrOp.indices.length; i++) partial[i] = F_int_full[hrOp.indices[i]];
                    F_int_approx = HRTrainer.applyGappyPOD(partial, hrOp);
                } else if (method === 'udeim') {
                    const elemForces = [];
                    for (const eId of hrOp.elemIndices) {
                        elemForces[eId] = this.hre.elementForce(patch, u_full, allElems[eId]);
                    }
                    // Fill missing with zeros
                    for (let e = 0; e < allElems.length; e++)
                        if (!elemForces[e]) elemForces[e] = new Float64Array(nDofs);
                    F_int_approx = HRTrainer.applyUDEIM(elemForces, hrOp);
                }

                // For tangent: use element sampling from associated elements
                const Kt_full = this.fom.calculateTangentStiffness(patch, u_full);
                this.fom.applyPenaltyConstraints(Kt_full, null, u_full, patch);
                const PhiT = Phi.transpose();
                const Kt_mat = new Matrix(Kt_full);
                const Kt_red = PhiT.mmul(Kt_mat).mmul(Phi).to2DArray();

                const R_full = F_ext.map((f, i) => f - F_int_approx[i]);
                const R_red_mat = PhiT.mmul(new Matrix([Array.from(R_full)]).transpose());
                const R_red = R_red_mat.to2DArray().map(r => r[0]);

                let norm = Math.sqrt(R_red.reduce((s, v) => s + v*v, 0));
                residualHistory.push({step: s, iter, norm});
                if (norm < tolerance && iter > 0) break;

                const dur = this.fom.gaussianElimination(Kt_red, R_red);
                for (let i = 0; i < k; i++) ur[i] += dur[i];
            }
        }

        const final_ur = new Matrix([Array.from(ur)]).transpose();
        const u = new Float64Array(Phi.mmul(final_ur).to2DArray().map(r => r[0]));
        return { u, ur, residualHistory, sampledIndices: hrOp.indices ? hrOp.indices.length : 0 };
    }

    /**
     * Solve using Collocation (LSPG)
     */
    solveCollocation(patch, bcs, loads, hrOp, options = {}) {
        const { iterations = 15, tolerance = 1e-6, steps = 1 } = options;
        const Phi = this.rom.Phi;
        if (!Phi) throw new Error("POD basis not computed");
        const { Matrix } = window.mlMatrix;
        const k = Phi.columns, nDofs = Phi.rows;
        const nV = patch.controlPoints[0].length;
        let ur = new Float64Array(k);
        const residualHistory = [];
        const m = hrOp.indices.length;

        const F_ext_total = new Float64Array(nDofs);
        loads.forEach(l => {
            const idx = (l.i * nV + l.j) * 2;
            F_ext_total[idx] += l.fx;
            F_ext_total[idx + 1] += l.fy;
        });

        for (let s = 1; s <= steps; s++) {
            const F_ext = F_ext_total.map(f => f * (s/steps));
            for (let iter = 0; iter < iterations; iter++) {
                const ur_mat = new Matrix([Array.from(ur)]).transpose();
                const u_full = new Float64Array(Phi.mmul(ur_mat).to2DArray().map(r => r[0]));

                const F_int = this.fom.calculateInternalForce(patch, u_full);
                this.fom.applyPenaltyConstraints(null, F_int, u_full, patch);
                const Kt = this.fom.calculateTangentStiffness(patch, u_full);
                this.fom.applyPenaltyConstraints(Kt, null, u_full, patch);

                // Sampled residual
                const R_sampled = new Float64Array(m);
                for (let i = 0; i < m; i++) R_sampled[i] = F_ext[hrOp.indices[i]] - F_int[hrOp.indices[i]];

                // Sampled Jacobian * Phi: [m x k]
                const JPhi = Array.from({length: m}, () => new Float64Array(k));
                for (let i = 0; i < m; i++)
                    for (let j = 0; j < k; j++) {
                        let dot = 0;
                        for (let d = 0; d < nDofs; d++) dot += Kt[hrOp.indices[i]][d] * Phi.get(d, j);
                        JPhi[i][j] = dot;
                    }

                // LSPG: dur = (JPhi^T JPhi)^{-1} JPhi^T R_sampled
                const JtJ = Array.from({length: k}, () => new Float64Array(k));
                const JtR = new Float64Array(k);
                for (let i = 0; i < k; i++) {
                    for (let j = 0; j < k; j++) {
                        let s2 = 0; for (let l = 0; l < m; l++) s2 += JPhi[l][i]*JPhi[l][j];
                        JtJ[i][j] = s2;
                    }
                    let s2 = 0; for (let l = 0; l < m; l++) s2 += JPhi[l][i]*R_sampled[l];
                    JtR[i] = s2;
                }

                let norm = Math.sqrt(JtR.reduce((s2, v) => s2 + v*v, 0));
                residualHistory.push({step: s, iter, norm});
                if (norm < tolerance && iter > 0) break;

                const dur = this.fom.gaussianElimination(JtJ, Array.from(JtR));
                for (let i = 0; i < k; i++) ur[i] += dur[i];
            }
        }

        const final_ur = new Matrix([Array.from(ur)]).transpose();
        const u = new Float64Array(Phi.mmul(final_ur).to2DArray().map(r => r[0]));
        return { u, ur, residualHistory, collocationPoints: m };
    }
}

window.HRSolver = HRSolver;
