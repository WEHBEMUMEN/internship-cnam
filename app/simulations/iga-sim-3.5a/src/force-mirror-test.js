/**
 * FORCE MIRROR DIAGNOSTIC
 * Compares FOM and DEIM force assembly for a single DOF to find scaling/integration errors.
 */
DEIMBenchmarkApp.prototype.runDiagnostic = function() {
    console.log("%c--- DEIM FORCE MIRROR DIAGNOSTIC ---", "font-weight:bold; color:#8b5cf6; font-size:1.2em;");
    
    if (!this.isTrained || !this.snapDisp.length) {
        console.error("Please train DEIM first!");
        return;
    }

    const patch = this.patch;
    const u = this.snapDisp[0]; // Test on first snapshot
    const sensorIdx = this.deimEngine.indices[0]; // Test on first sensor
    
    console.log(`Testing Sensor: DOF ${sensorIdx}`);

    // 1. FOM Calculation for this DOF (Full Global Assembly loop for 1 DOF)
    let f_fom = 0;
    const nU = patch.controlPoints.length, nV = patch.controlPoints[0].length;
    const uniqueU = [...new Set(patch.U)], uniqueV = [...new Set(patch.V)];
    const gRule = GaussQuadrature2D.getPoints(Math.max(patch.p, patch.q) + 1);

    for (let i = 0; i < uniqueU.length - 1; i++) {
        for (let j = 0; j < uniqueV.length - 1; j++) {
            const uMin = uniqueU[i], uMax = uniqueU[i+1];
            const vMin = uniqueV[j], vMax = uniqueV[j+1];
            if (uMax - uMin < 1e-10 || vMax - vMin < 1e-10) continue;

            for (let gu = 0; gu < gRule.points.length; gu++) {
                const uVal = ((uMax - uMin) * gRule.points[gu] + (uMax + uMin)) / 2;
                const wu = gRule.weights[gu] * (uMax - uMin) / 2;
                for (let gv = 0; gv < gRule.points.length; gv++) {
                    const vVal = ((vMax - vMin) * gRule.points[gv] + (vMax + vMin)) / 2;
                    const wv = gRule.weights[gv] * (vMax - vMin) / 2;

                    const deriv = this.solverFOM.engine.getSurfaceDerivatives(patch, uVal, vVal);
                    const { grads: B_param, detJ, activeIndices } = this.solverFOM.getBParametric(patch, uVal, vVal, deriv);
                    
                    // Check if our sensor DOF is in this element
                    const cpIdx = Math.floor(sensorIdx / 2);
                    const localIdx = activeIndices.indexOf(cpIdx);
                    if (localIdx === -1) continue;

                    // Calculate stress at this point
                    let dudx=0, dudy=0, dvdx=0, dvdy=0;
                    activeIndices.forEach(k => {
                        dudx += B_param[k][0] * u[k*2];
                        dudy += B_param[k][1] * u[k*2];
                        dvdx += B_param[k][0] * u[k*2+1];
                        dvdy += B_param[k][1] * u[k*2+1];
                    });

                    const Exx = dudx + 0.5 * (dudx*dudx + dvdx*dvdx);
                    const Eyy = dvdy + 0.5 * (dudy*dudy + dvdy*dvdy);
                    const Exy2 = (dudy + dvdx) + (dudx*dudy + dvdx*dvdy);
                    
                    const D = this.solverFOM.getPlaneStressD();
                    const Sxx = D[0][0]*Exx + D[0][1]*Eyy;
                    const Syy = D[1][0]*Exx + D[1][1]*Eyy;
                    const Sxy = D[2][2]*Exy2;

                    const factor = detJ * wu * wv * this.solverFOM.thickness;
                    const axis = sensorIdx % 2;
                    const dRdx = B_param[cpIdx][0], dRdy = B_param[cpIdx][1];
                    
                    let val = 0;
                    if (axis === 0) {
                        val = ((1 + dudx)*dRdx * Sxx + (dudy)*dRdy * Syy + ((1 + dudx)*dRdy + dudy*dRdx) * Sxy) * factor;
                    } else {
                        val = ((dvdx)*dRdx * Sxx + (1 + dvdy)*dRdy * Syy + ((1 + dvdy)*dRdx + dvdx*dRdy) * Sxy) * factor;
                    }
                    if (!isNaN(val)) f_fom += val;
                }
            }
        }
    }

    // 2. DEIM Calculation for this same DOF (Optimized Assembly)
    const f_deim_vec = this.deimEngine.calculateSampledInternalForce(this.solverFOM, patch, u);
    const f_deim = f_deim_vec[0]; // indices[0] corresponds to f_deim_vec[0]

    console.log(`FOM  Value: ${f_fom.toExponential(6)}`);
    console.log(`DEIM Value: ${f_deim.toExponential(6)}`);
    console.log(`Ratio (DEIM/FOM): ${(f_deim/f_fom).toFixed(4)}`);
    
    if (Math.abs(f_deim - f_fom) > 1e-8 * Math.abs(f_fom)) {
        console.error("!!! DISCREPANCY DETECTED !!!");
        console.log("Check if activeElements in DEIM is missing some knot spans.");
    } else {
        console.log("%c✓ Physical Consistency Verified", "color:green");
    }
};
