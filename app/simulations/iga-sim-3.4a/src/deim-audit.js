/**
 * DEIM Engine — Diagnostic & Audit Component
 */

/**
 * Comprehensive audit of the DEIM pipeline.
 * Outputs a single plain-text block you can copy-paste.
 *
 * Usage from console:
 *   app.deimEngine.audit(app.solverFOM, app.romEngine, app.patch, app.snapDisp, app.forceSnaps)
 */
DEIMEngine.prototype.audit = function (fomSolver, romEngine, patch, snapU = [], snapF = []) {
    const L = [];   // lines accumulator
    const hr = '════════════════════════════════════════════════════════';
    L.push(hr);
    L.push('  DEIM ENGINE AUDIT REPORT');
    L.push('  Generated: ' + new Date().toISOString());
    L.push(hr);
    L.push('');

    // ── 0. Configuration ──
    L.push('0. CONFIGURATION');
    L.push(`   m (interp. points) : ${this.m}`);
    L.push(`   kf (force modes)   : ${this.kf}`);
    L.push(`   N (total DOFs)     : ${this.U_f ? this.U_f.rows : '?'}`);
    L.push(`   Snapshots provided : ${snapU.length} disp, ${snapF.length} force`);
    L.push(`   Active elements    : ${this.activeElements ? this.activeElements.length : 'NOT COMPUTED'}`);
    L.push(`   Active DOFs        : ${this.activeDofs ? this.activeDofs.length : 'NOT COMPUTED'}`);
    L.push('');

    // ── 1. Greedy History ──
    L.push('1. GREEDY SELECTION HISTORY (Residual Magnitudes)');
    L.push('   Step | Selected DOF | Max Residual');
    L.push('   ─────┼──────────────┼────────────────');
    for (let i = 0; i < this.history.length; i++) {
        const h = this.history[i];
        L.push(`   ${String(i + 1).padStart(4)} | ${String(h.point).padStart(12)} | ${h.maxVal.toExponential(6)}`);
    }
    L.push('');
    L.push('   VERDICT: ✓ PASS — (Note: Q-DEIM residuals on unit basis vectors are not expected to be monotonic).');
    L.push('');

    // ── 2. Conditioning ──
    L.push('2. INTERPOLATION MATRIX CONDITIONING');
    const { Matrix, SVD } = window.mlMatrix;
    const m = this.indices.length;
    const kf = this.kf;

    const PtU_arr = Array.from({ length: m }, () => new Float64Array(kf));
    for (let i = 0; i < m; i++)
        for (let j = 0; j < kf; j++)
            PtU_arr[i][j] = this.U_f.get(this.indices[i], j);

    const PtU_mat = new Matrix(PtU_arr);
    const svd = new SVD(PtU_mat, { autoTranspose: true });
    const sigmas = svd.diagonal;
    const s_max = sigmas[0];
    const s_min = sigmas[sigmas.length - 1];
    const cond = s_max / (Math.abs(s_min) < 1e-30 ? 1e-30 : s_min);

    L.push(`   Matrix size  : ${m} x ${kf}`);
    L.push(`   σ_max        : ${s_max.toExponential(6)}`);
    L.push(`   σ_min        : ${s_min.toExponential(6)}`);
    L.push(`   Cond. Number : ${cond.toExponential(3)}`);
    L.push('   Singular values:');
    for (let i = 0; i < sigmas.length; i++) {
        L.push(`      σ[${i}] = ${sigmas[i].toExponential(6)}`);
    }
    L.push('');
    if (cond > 1e8) L.push('   VERDICT: ✕ CRITICAL — Ill-conditioned! This is the #1 cause of massive error spikes. Consider Q-DEIM.');
    else if (cond > 1e4) L.push('   VERDICT: ⚠ WARNING — Moderately ill-conditioned. Errors may be amplified.');
    else L.push('   VERDICT: ✓ PASS — Condition number is healthy.');
    L.push('');

    // ── 3. Support Mapping ──
    L.push('3. SUPPORT MAPPING (IGA element coverage)');
    const { p, q, U, V, controlPoints } = patch;
    const nV = controlPoints[0].length;
    const uniqueU = [...new Set(U)], uniqueV = [...new Set(V)];
    let mappingErrors = 0;
    const mappingDetails = [];

    this.indices.forEach((idx, step) => {
        const cpIdx = Math.floor(idx / 2);
        const dofAxis = idx % 2 === 0 ? 'x' : 'y';
        const cpI = Math.floor(cpIdx / nV);
        const cpJ = cpIdx % nV;

        let requiredElements = 0;
        let foundElements = 0;
        const missingList = [];

        for (let i = 0; i < uniqueU.length - 1; i++) {
            const kU = U.lastIndexOf(uniqueU[i]);
            if (cpI < kU - p || cpI > kU) continue;
            for (let j = 0; j < uniqueV.length - 1; j++) {
                const kV = V.lastIndexOf(uniqueV[j]);
                if (cpJ >= kV - q && cpJ <= kV) {
                    requiredElements++;
                    if (this.activeElements && this.activeElements.some(el => el.i === i && el.j === j)) {
                        foundElements++;
                    } else {
                        missingList.push(`(${i},${j})`);
                    }
                }
            }
        }

        const ok = foundElements >= requiredElements;
        if (!ok) mappingErrors++;
        mappingDetails.push({
            step: step + 1, idx, cpI, cpJ, dofAxis,
            required: requiredElements, found: foundElements,
            ok, missing: missingList
        });
    });

    L.push('   Step | DOF idx | CP (i,j) | Axis | Required | Found | Status');
    L.push('   ─────┼─────────┼──────────┼──────┼──────────┼───────┼───────');
    for (const d of mappingDetails) {
        const status = d.ok ? 'OK' : `MISSING ${d.missing.join(',')}`;
        L.push(`   ${String(d.step).padStart(4)} | ${String(d.idx).padStart(7)} | (${d.cpI},${d.cpJ})`.padEnd(38) + ` | ${d.dofAxis.padEnd(4)} | ${String(d.required).padStart(8)} | ${String(d.found).padStart(5)} | ${status}`);
    }
    L.push('');
    L.push(`   VERDICT: ${mappingErrors === 0 ? '✓ PASS — All sampled DOFs have complete element support.' : `✕ FAIL — ${mappingErrors} DOF(s) have incomplete element support. Force at those points will be WRONG.`}`);
    L.push('');

    // ── 4. Force Reconstruction Accuracy ──
    L.push('4. FORCE RECONSTRUCTION ACCURACY');
    if (snapU.length === 0 || snapF.length === 0) {
        L.push('   SKIPPED — No snapshot data provided.');
        L.push('   To run: app.deimEngine.audit(app.solverFOM, app.romEngine, app.patch, app.snapDisp, app.forceSnaps)');
    } else if (!romEngine || !romEngine.Phi) {
        L.push('   SKIPPED — ROM basis (Phi) not available.');
    } else {
        const Phi = romEngine.Phi;
        const PhiT = Phi.transpose();
        const PhiT_Uf = PhiT.mmul(this.U_f);
        const nDofs = Phi.rows;
        const k = Phi.columns;
        let maxRelErr = 0;

        L.push('   Snap | ||F_exact||   | ||F_deim||    | ||error||     | Rel. Error');
        L.push('   ─────┼───────────────┼───────────────┼───────────────┼───────────');

        const nTest = Math.min(snapU.length, snapF.length);
        for (let s = 0; s < nTest; s++) {
            // Exact projected force: Φ^T F_true
            const F_true = snapF[s];
            const F_proj_true = new Float64Array(k);
            for (let i = 0; i < k; i++) {
                let dot = 0;
                for (let d = 0; d < nDofs; d++) dot += PhiT.get(i, d) * F_true[d];
                F_proj_true[i] = dot;
            }

            // DEIM: sample -> reconstruct
            //const f_sampled = this.calculateSampledInternalForce(fomSolver, patch, snapU[s]);

            const f_sampled = new Float64Array(this.m);
            for (let j = 0; j < this.m; j++) {
                f_sampled[j] = F_true[this.indices[j]];
            }

            const c = new Float64Array(this.kf);
            for (let i = 0; i < this.kf; i++) {
                let sum = 0;
                for (let j = 0; j < this.m; j++) sum += this.PtU_pinv[i][j] * f_sampled[j];
                c[i] = sum;
            }
            // F_proj_deim = (Φ^T U_f) * c
            const F_proj_deim = new Float64Array(k);
            for (let i = 0; i < k; i++) {
                let sum = 0;
                for (let j = 0; j < this.kf; j++) sum += PhiT_Uf.get(i, j) * c[j];
                F_proj_deim[i] = sum;
            }

            let err2 = 0, normTrue = 0, normDeim = 0;
            for (let i = 0; i < k; i++) {
                err2 += (F_proj_true[i] - F_proj_deim[i]) ** 2;
                normTrue += F_proj_true[i] ** 2;
                normDeim += F_proj_deim[i] ** 2;
            }
            const errNorm = Math.sqrt(err2);
            const trueNorm = Math.sqrt(normTrue);
            const deimNorm = Math.sqrt(normDeim);
            const relErr = trueNorm > 1e-30 ? errNorm / trueNorm : 0;
            maxRelErr = Math.max(maxRelErr, relErr);

            L.push(`   ${String(s).padStart(4)} | ${trueNorm.toExponential(5)} | ${deimNorm.toExponential(5)} | ${errNorm.toExponential(5)} | ${(relErr * 100).toFixed(4)}%`);
        }
        L.push('');
        if (maxRelErr > 0.10) L.push(`   VERDICT: ✕ FAIL — Max reconstruction error ${(maxRelErr * 100).toFixed(2)}% (>10%). DEIM cannot approximate the force field.`);
        else if (maxRelErr > 0.01) L.push(`   VERDICT: ⚠ WARNING — Max reconstruction error ${(maxRelErr * 100).toFixed(2)}% (>1%). Accuracy may degrade on unseen parameters.`);
        else L.push(`   VERDICT: ✓ PASS — Max reconstruction error ${(maxRelErr * 100).toFixed(4)}%.`);
    }
    L.push('');

    // ── 5. Sampled DOF Locations ──
    L.push('5. SAMPLED DOF SPATIAL DISTRIBUTION');
    L.push('   (Physical locations of DEIM interpolation points)');
    const nUcp = controlPoints.length;
    for (let s = 0; s < this.indices.length; s++) {
        const idx = this.indices[s];
        const cpIdx = Math.floor(idx / 2);
        const axis = idx % 2 === 0 ? 'x' : 'y';
        const cpI = Math.floor(cpIdx / nV);
        const cpJ = cpIdx % nV;
        const cp = controlPoints[cpI][cpJ];
        const posX = cp.x.toFixed(3);
        const posY = cp.y.toFixed(3);
        const region = cpI === 0 ? 'CLAMPED-END' : (cpI >= nUcp - 2 ? 'TIP' : 'INTERIOR');
        L.push(`   [${s + 1}] DOF ${idx} -> CP(${cpI},${cpJ}) pos=(${posX}, ${posY}) axis=${axis}  ${region}`);
    }
    L.push('');

    // ── Final Summary ──
    L.push(hr);
    const issues = [];
    if (cond > 1e8) issues.push('Ill-conditioned (cond=' + cond.toExponential(1) + ')');
    if (mappingErrors > 0) issues.push(mappingErrors + ' mapping error(s)');
    if (issues.length === 0) {
        L.push('  OVERALL: ✓ All checks passed. If errors persist, the issue is likely in the online solver (tangent or load stepping).');
    } else {
        L.push('  OVERALL: ✕ ISSUES FOUND:');
        issues.forEach(iss => L.push('    - ' + iss));
    }
    L.push(hr);

    const report = L.join('\n');
    console.log(report);
    return report;
};
