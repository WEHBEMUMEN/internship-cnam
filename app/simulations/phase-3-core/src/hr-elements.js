/**
 * Element-Level Utilities for Hyper-Reduction (Phase 3.4)
 * Provides per-element F_int and K_T without modifying the FOM solver.
 */
class HRElements {
    constructor(solver, nurbsEngine) {
        this.solver = solver;
        this.engine = nurbsEngine;
    }

    /** Enumerate all non-degenerate knot-span elements */
    getElements(patch) {
        const uU = [...new Set(patch.U)], uV = [...new Set(patch.V)];
        const elems = [];
        for (let i = 0; i < uU.length - 1; i++) {
            if (uU[i+1] - uU[i] < 1e-10) continue;
            for (let j = 0; j < uV.length - 1; j++) {
                if (uV[j+1] - uV[j] < 1e-10) continue;
                elems.push({ uMin: uU[i], uMax: uU[i+1], vMin: uV[j], vMax: uV[j+1], id: elems.length });
            }
        }
        return elems;
    }

    /** Compute element-level internal force for ONE element */
    elementForce(patch, u_disp, elem) {
        const { p, q, U, V, controlPoints } = patch;
        const nU = controlPoints.length, nV = controlPoints[0].length;
        const nDofs = nU * nV * 2;
        const fe = new Float64Array(nDofs);
        const gRule = GaussQuadrature2D.getPoints(Math.max(p, q) + 1);

        for (let gu = 0; gu < gRule.points.length; gu++) {
            const u = ((elem.uMax - elem.uMin) * gRule.points[gu] + (elem.uMax + elem.uMin)) / 2;
            const wu = gRule.weights[gu] * (elem.uMax - elem.uMin) / 2;
            for (let gv = 0; gv < gRule.points.length; gv++) {
                const v = ((elem.vMax - elem.vMin) * gRule.points[gv] + (elem.vMax + elem.vMin)) / 2;
                const wv = gRule.weights[gv] * (elem.vMax - elem.vMin) / 2;

                const deriv = this.engine.getSurfaceDerivatives(patch, u, v);
                const { grads, detJ, activeIndices } = this.solver.getBParametric(patch, u, v, deriv);

                let dudx=0, dudy=0, dvdx=0, dvdy=0;
                for (let a = 0; a < activeIndices.length; a++) {
                    const k = activeIndices[a];
                    dudx += grads[k][0] * u_disp[k*2];
                    dudy += grads[k][1] * u_disp[k*2];
                    dvdx += grads[k][0] * u_disp[k*2+1];
                    dvdy += grads[k][1] * u_disp[k*2+1];
                }

                const Exx = dudx + 0.5*(dudx*dudx + dvdx*dvdx);
                const Eyy = dvdy + 0.5*(dudy*dudy + dvdy*dvdy);
                const Exy2 = (dudy + dvdx) + (dudx*dudy + dvdx*dvdy);

                const D = this.solver.getPlaneStressD();
                const Sxx = D[0][0]*Exx + D[0][1]*Eyy;
                const Syy = D[1][0]*Exx + D[1][1]*Eyy;
                const Sxy = D[2][2]*Exy2;
                const factor = detJ * wu * wv * this.solver.thickness;

                for (let a = 0; a < activeIndices.length; a++) {
                    const k = activeIndices[a];
                    const dRdx = grads[k][0], dRdy = grads[k][1];
                    const bu = (1+dudx)*dRdx, bv = dvdx*dRdx;
                    const cu = dudy*dRdy,     cv = (1+dvdy)*dRdy;
                    const eu = (1+dudx)*dRdy + dudy*dRdx;
                    const ev = (1+dvdy)*dRdx + dvdx*dRdy;

                    const vu = (bu*Sxx + cu*Syy + eu*Sxy) * factor;
                    const vv = (bv*Sxx + cv*Syy + ev*Sxy) * factor;
                    if (!isNaN(vu)) fe[k*2] += vu;
                    if (!isNaN(vv)) fe[k*2+1] += vv;
                }
            }
        }
        return fe;
    }

    /** Compute element-level tangent stiffness for ONE element */
    elementTangent(patch, u_disp, elem) {
        const { p, q, U, V, controlPoints } = patch;
        const nU = controlPoints.length, nV = controlPoints[0].length;
        const nDofs = nU * nV * 2;
        const ke = Array.from({length: nDofs}, () => new Float64Array(nDofs));
        const gRule = GaussQuadrature2D.getPoints(Math.max(p, q) + 1);
        const nPts = nU * nV;

        for (let gu = 0; gu < gRule.points.length; gu++) {
            const u = ((elem.uMax - elem.uMin) * gRule.points[gu] + (elem.uMax + elem.uMin)) / 2;
            const wu = gRule.weights[gu] * (elem.uMax - elem.uMin) / 2;
            for (let gv = 0; gv < gRule.points.length; gv++) {
                const v = ((elem.vMax - elem.vMin) * gRule.points[gv] + (elem.vMax + elem.vMin)) / 2;
                const wv = gRule.weights[gv] * (elem.vMax - elem.vMin) / 2;

                const deriv = this.engine.getSurfaceDerivatives(patch, u, v);
                const { grads, detJ, activeIndices } = this.solver.getBParametric(patch, u, v, deriv);

                let dudx=0, dudy=0, dvdx=0, dvdy=0;
                for (let aIdx = 0; aIdx < activeIndices.length; aIdx++) {
                    const k = activeIndices[aIdx];
                    dudx += grads[k][0] * u_disp[k*2];
                    dudy += grads[k][1] * u_disp[k*2];
                    dvdx += grads[k][0] * u_disp[k*2+1];
                    dvdy += grads[k][1] * u_disp[k*2+1];
                }

                const D = this.solver.getPlaneStressD();
                const factor = detJ * wu * wv * this.solver.thickness;

                // Material tangent
                for (let aIdx = 0; aIdx < activeIndices.length; aIdx++) {
                    const a = activeIndices[aIdx];
                    const ax = grads[a][0], ay = grads[a][1];
                    const Ba = [[(1+dudx)*ax, dvdx*ax], [dudy*ay, (1+dvdy)*ay],
                                [(1+dudx)*ay+dudy*ax, (1+dvdy)*ax+dvdx*ay]];
                    for (let bIdx = 0; bIdx < activeIndices.length; bIdx++) {
                        const b = activeIndices[bIdx];
                        const bx = grads[b][0], by = grads[b][1];
                        const Bb = [[(1+dudx)*bx, dvdx*bx], [dudy*by, (1+dvdy)*by],
                                    [(1+dudx)*by+dudy*bx, (1+dvdy)*bx+dvdx*by]];
                        for (let ri = 0; ri < 2; ri++) {
                            for (let rj = 0; rj < 2; rj++) {
                                let kab = 0;
                                for (let r = 0; r < 3; r++)
                                    for (let c = 0; c < 3; c++)
                                        kab += Ba[r][ri] * D[r][c] * Bb[c][rj];
                                ke[a*2+ri][b*2+rj] += kab * factor;
                            }
                        }
                    }
                }

                // Geometric stiffness
                const Exx = dudx + 0.5*(dudx*dudx + dvdx*dvdx);
                const Eyy = dvdy + 0.5*(dudy*dudy + dvdy*dvdy);
                const Exy2 = (dudy+dvdx) + (dudx*dudy + dvdx*dvdy);
                const Sxx = D[0][0]*Exx + D[0][1]*Eyy;
                const Syy = D[1][0]*Exx + D[1][1]*Eyy;
                const Sxy = D[2][2]*Exy2;

                for (let aIdx = 0; aIdx < activeIndices.length; aIdx++) {
                    const a = activeIndices[aIdx];
                    for (let bIdx = 0; bIdx < activeIndices.length; bIdx++) {
                        const b = activeIndices[bIdx];
                        const kg = (grads[a][0]*Sxx*grads[b][0] + grads[a][1]*Syy*grads[b][1] +
                                    grads[a][0]*Sxy*grads[b][1] + grads[a][1]*Sxy*grads[b][0]) * factor;
                        const v = isNaN(kg) ? 0 : kg;
                        ke[a*2][b*2] += v;
                        ke[a*2+1][b*2+1] += v;
                    }
                }
            }
        }
        return ke;
    }

    /** Collect all element forces as an array of Float64Arrays */
    allElementForces(patch, u_disp) {
        const elems = this.getElements(patch);
        return elems.map(e => this.elementForce(patch, u_disp, e));
    }
}

window.HRElements = HRElements;
