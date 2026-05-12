/**
 * Phase 5.1b - ECSW Solver Extension
 * Adds element-wise assembly and multiplication for Hyper-reduction.
 */

class ECSWSolver extends window.IGA2DSolver {
    constructor(engine) {
        super(engine);
    }

    /**
     * Multiplies the global stiffness matrix by a vector WITHOUT assembling it.
     * K * v = sum( Ke * v )
     */
    multiplyK(patch, v) {
        const nDofs = v.length;
        const result = new Float64Array(nDofs).fill(0);
        if (!patch.elements) {
            patch.elements = window.GeometryFactory.extractElements(patch);
        }
        const nElements = patch.elements.length;

        for (let eIdx = 0; eIdx < nElements; eIdx++) {
            const Ke_v = this.multiplyKe(patch, eIdx, v);
            for (let i = 0; i < nDofs; i++) {
                result[i] += Ke_v[i];
            }
        }
        return result;
    }

    /**
     * Multiplies a single element stiffness matrix by a vector.
     * Ke * v
     */
    multiplyKe(patch, eIdx, v) {
        const { p, q, U, V, weights, controlPoints } = patch;
        const nU = controlPoints.length;
        const nV = controlPoints[0].length;
        const nDofs = nU * nV * 2;
        const result = new Float64Array(nDofs).fill(0);

        const element = patch.elements[eIdx];
        const uMin = element.uRange[0], uMax = element.uRange[1];
        const vMin = element.vRange[0], vMax = element.vRange[1];

        const gRule = GaussQuadrature2D.getPoints(Math.max(p, q) + 1);
        const D = this.getPlaneStressD();

        for (let gu = 0; gu < gRule.points.length; gu++) {
            const u = ((uMax - uMin) * gRule.points[gu] + (uMax + uMin)) / 2;
            const wu = gRule.weights[gu] * (uMax - uMin) / 2;
            for (let gv = 0; gv < gRule.points.length; gv++) {
                const v_val = ((vMax - vMin) * gRule.points[gv] + (vMax + vMin)) / 2;
                const wv = gRule.weights[gv] * (vMax - vMin) / 2;

                const detJ = this.engine.getJacobianDeterminant(patch, u, v_val);
                const deriv = this.engine.getSurfaceDerivatives(patch, u, v_val);
                const B = this.getBMatrix(patch, u, v_val, deriv);
                const factor = detJ * wu * wv * this.thickness;

                // Ke_v = B' * D * B * v * factor
                // 1. Strain eps = B * v
                const eps = [0, 0, 0];
                for (let a = 0; a < nU * nV; a++) {
                    eps[0] += B[a][0][0] * v[a * 2]     + B[a][0][1] * v[a * 2 + 1];
                    eps[1] += B[a][1][0] * v[a * 2]     + B[a][1][1] * v[a * 2 + 1];
                    eps[2] += B[a][2][0] * v[a * 2]     + B[a][2][1] * v[a * 2 + 1];
                }

                // 2. Stress sigma = D * eps
                const sig = [
                    D[0][0] * eps[0] + D[0][1] * eps[1],
                    D[1][0] * eps[0] + D[1][1] * eps[1],
                    D[2][2] * eps[2]
                ];

                // 3. Force f = B' * sigma * factor
                for (let a = 0; a < nU * nV; a++) {
                    result[a * 2]     += (B[a][0][0] * sig[0] + B[a][1][0] * sig[1] + B[a][2][0] * sig[2]) * factor;
                    result[a * 2 + 1] += (B[a][0][1] * sig[0] + B[a][1][1] * sig[1] + B[a][2][1] * sig[2]) * factor;
                }
            }
        }
        return result;
    }

    /**
     * Assembles the REDUCED stiffness matrix directly (Online ECSW)
     * Kr = sum_e w_e * Phi_e' * Ke * Phi_e
     */
    assembleECSWReduced(patch, sampledIndices, weights, Phi) {
        const k = Phi.columns;
        const { Matrix } = window.mlMatrix || window.ML;
        const Kr = new Matrix(k, k);

        const nU = patch.controlPoints.length;
        const nV = patch.controlPoints[0].length;

        for (let i = 0; i < sampledIndices.length; i++) {
            const eIdx = sampledIndices[i];
            const w = weights[i];
            
            // 1. Get element stiffness
            const Ke = this.assembleElementStiffness(patch, eIdx);
            
            // 2. Extract basis for this element's DOFs
            // For simplicity, we can project the whole Ke using Phi
            // But Ke is large (nDofs x nDofs) with most rows zero.
            // Optimized: Only multiply active rows/cols
            const Ke_Phi = new Matrix(Ke).mmul(Phi);
            const PhiT_Ke_Phi = Phi.transpose().mmul(Ke_Phi);
            
            // 3. Accumulate weighted contribution
            for (let r = 0; r < k; r++) {
                for (let c = 0; c < k; c++) {
                    const val = PhiT_Ke_Phi.get(r, c);
                    Kr.set(r, c, Kr.get(r, c) + val * w);
                }
            }
        }
        return Kr;
    }

    /**
     * Assembles the FULL reduced stiffness matrix (Standard Galerkin ROM)
     * Kr = Phi' * K * Phi = sum_e Phi_e' * Ke * Phi_e
     */
    assembleReducedStiffness(patch, Phi) {
        if (!patch.elements) {
            patch.elements = window.GeometryFactory.extractElements(patch);
        }
        const nElements = patch.elements.length;
        const weights = new Float64Array(nElements).fill(1.0);
        const indices = Array.from({ length: nElements }, (_, i) => i);
        return this.assembleECSWReduced(patch, indices, weights, Phi);
    }

    assembleElementStiffness(patch, eIdx) {
        const { p, q, U, V, weights, controlPoints } = patch;
        const nU = controlPoints.length;
        const nV = controlPoints[0].length;
        const nDofs = nU * nV * 2;
        const Ke = Array.from({ length: nDofs }, () => new Float64Array(nDofs).fill(0));

        const element = patch.elements[eIdx];
        const uMin = element.uRange[0], uMax = element.uRange[1];
        const vMin = element.vRange[0], vMax = element.vRange[1];

        const gRule = GaussQuadrature2D.getPoints(Math.max(p, q) + 1);
        const D = this.getPlaneStressD();

        for (let gu = 0; gu < gRule.points.length; gu++) {
            const u = ((uMax - uMin) * gRule.points[gu] + (uMax + uMin)) / 2;
            const wu = gRule.weights[gu] * (uMax - uMin) / 2;
            for (let gv = 0; gv < gRule.points.length; gv++) {
                const v_val = ((vMax - vMin) * gRule.points[gv] + (vMax + vMin)) / 2;
                const wv = gRule.weights[gv] * (vMax - vMin) / 2;

                const detJ = this.engine.getJacobianDeterminant(patch, u, v_val);
                const deriv = this.engine.getSurfaceDerivatives(patch, u, v_val);
                const B = this.getBMatrix(patch, u, v_val, deriv);
                const factor = detJ * wu * wv * this.thickness;

                this.accumulateContribution(Ke, B, D, factor, nU, nV);
            }
        }
        return Ke;
    }
}

window.ECSWSolver = ECSWSolver;
