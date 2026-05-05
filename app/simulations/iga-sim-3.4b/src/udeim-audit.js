/**
 * U-DEIM Engine — Advanced Mathematical Audit (v2)
 * Phase 3.4b | Deep Pipeline Verification
 */

UDEIMEngine.prototype.auditMath = function (fomSolver, romEngine, patch, snapU = []) {
    const L = [];
    const hr = '════════════════════════════════════════════════════════════';
    const subhr = '────────────────────────────────────────────────────────────';
    const { Matrix, SVD } = window.mlMatrix;

    L.push(hr);
    L.push('  U-DEIM ADVANCED MATHEMATICAL AUDIT');
    L.push('  Generated: ' + new Date().toISOString());
    L.push(hr);

    // --- 1. CONFIGURATION ---
    const N_u = this.numElements * this.numLocalDofs;
    const N = romEngine.Phi.rows;
    const k = romEngine.Phi.columns;
    const kf = this.kf;
    const m = this.m;

    L.push('\n1. DIMENSIONALITY & BASIS');
    L.push(`   Global DOFs (N)     : ${N}`);
    L.push(`   Unassembled DOFs(Nu): ${N_u}`);
    L.push(`   POD Modes (k)       : ${k}`);
    L.push(`   Force Modes (kf)    : ${kf}`);
    L.push(`   Sensors (m)         : ${m}`);
    L.push(`   Compression Ratio   : ${(N / k).toFixed(1)}x (State) | ${(N_u / m).toFixed(1)}x (Force)`);

    // --- 2. INTERPOLATION CONDITIONING ---
    L.push('\n2. INTERPOLATION MATRIX CONDITIONING (Pᵀ U_f)');
    const PtU_arr = Array.from({ length: m }, () => new Float64Array(kf));
    for (let i = 0; i < m; i++)
        for (let j = 0; j < kf; j++)
            PtU_arr[i][j] = this.U_f.get(this.indices[i], j);

    const PtU_mat = new Matrix(PtU_arr);
    const svd = new SVD(PtU_mat);
    const sigmas = svd.diagonal;
    const cond = sigmas[0] / (sigmas[sigmas.length - 1] || 1e-30);

    L.push(`   Matrix Size  : ${m} x ${kf}`);
    L.push(`   Cond. Number : ${cond.toExponential(3)}`);
    L.push(`   Min Sigma    : ${sigmas[sigmas.length - 1].toExponential(4)}`);
    
    if (cond > 1e10) L.push('   VERDICT: ❌ CRITICAL — Interpolation is numerically unstable!');
    else if (cond > 1e6) L.push('   VERDICT: ⚠️ WARNING — High sensitivity to noise.');
    else L.push('   VERDICT: ✅ PASS — Matrix is well-conditioned.');

    // --- 3. ASSEMBLY CONSISTENCY ---
    L.push('\n3. ASSEMBLY CONSISTENCY (Φᵀ A U_f)');
    // We check if the precomputed M matrix matches a manual assembly of the first force mode
    if (this.M && this.U_f) {
        const Phi = romEngine.Phi;
        const PhiT = Phi.transpose();
        
        // Manual assembly of Mode 0
        const f_mode0_u = new Float64Array(N_u);
        for(let i=0; i<N_u; i++) f_mode0_u[i] = this.U_f.get(i, 0);
        
        // Global assembly
        const f_mode0_assembled = new Float64Array(N);
        for(let e=0; e<this.numElements; e++) {
            const map = this.elementDofMap[e];
            for(let ld=0; ld<this.numLocalDofs; ld++) {
                f_mode0_assembled[map[ld]] += f_mode0_u[e * this.numLocalDofs + ld];
            }
        }
        
        // Projection
        const proj_true = new Float64Array(k);
        for(let i=0; i<k; i++) {
            for(let d=0; d<N; d++) proj_true[i] += PhiT.get(i, d) * f_mode0_assembled[d];
        }
        
        // DEIM Projection (M * Pᵀ * f_mode0_u)
        // Since f_mode0_u is the 0-th column of U_f, its sample at indices is the 0-th column of PtU
        const samples = new Float64Array(m);
        for(let j=0; j<m; j++) samples[j] = f_mode0_u[this.indices[j]];
        
        const proj_deim = new Float64Array(k);
        for(let i=0; i<k; i++) {
            for(let j=0; j<m; j++) proj_deim[i] += this.M[i][j] * samples[j];
        }
        
        let diff2 = 0, norm2 = 0;
        for(let i=0; i<k; i++) {
            diff2 += (proj_true[i] - proj_deim[i])**2;
            norm2 += proj_true[i]**2;
        }
        const err = Math.sqrt(diff2) / Math.sqrt(norm2 || 1e-30);
        L.push(`   Basis Projection Error: ${(err * 100).toFixed(6)}%`);
        if (err > 1e-8) L.push('   VERDICT: ❌ FAIL — M matrix assembly logic is inconsistent with FOM!');
        else L.push('   VERDICT: ✅ PASS — Assembly math is perfect.');
    }

    // --- 4. SNAPSHOT RECONSTRUCTION ---
    L.push('\n4. SNAPSHOT RECONSTRUCTION (Offline Check)');
    if (snapU.length > 0) {
        const nTest = Math.min(snapU.length, 3);
        const gRule = window.GaussQuadrature2D.getPoints(Math.max(patch.p, patch.q) + 1);
        
        for (let s = 0; s < nTest; s++) {
            const u = snapU[s];
            // 1. True Internal Force
            const f_true = fomSolver.calculateInternalForce(patch, u);
            // 2. True Projection
            const proj_true = new Float64Array(k);
            const PhiT = romEngine.Phi.transpose();
            for(let i=0; i<k; i++) {
                for(let d=0; d<N; d++) proj_true[i] += PhiT.get(i, d) * f_true[d];
            }
            
            // 3. U-DEIM Sampling
            const f_sampled = new Float64Array(m);
            for (let i = 0; i < m; i++) {
                const { e, ld } = this.sampledDofs[i];
                // Manually re-compute local force for this element to be independent of solver state
                const el = {
                    uMin: patch.U[Math.floor(e / (patch.V.length - 1))], // Rough mapping
                    uMax: patch.U[Math.floor(e / (patch.V.length - 1)) + 1],
                    vMin: patch.V[e % (patch.V.length - 1)],
                    vMax: patch.V[e % (patch.V.length - 1) + 1]
                };
                // Simplified: use element index directly
                const f_e = this._computeLocalForce(fomSolver, patch, this.activeElements.find(ael => ael.index === e) || this.elementDofMap[e], u, gRule);
                f_sampled[i] = f_e[ld];
            }
            
            // 4. DEIM Projection
            const proj_deim = new Float64Array(k);
            for(let i=0; i<k; i++) {
                for(let j=0; j<m; j++) proj_deim[i] += this.M[i][j] * f_sampled[j];
            }
            
            let err2 = 0, n2 = 0;
            for(let i=0; i<k; i++) {
                err2 += (proj_true[i] - proj_deim[i])**2;
                n2 += proj_true[i]**2;
            }
            const relErr = Math.sqrt(err2) / Math.sqrt(n2 || 1e-30);
            L.push(`   Snap ${s} | Galerkin Projection Error: ${(relErr * 100).toFixed(4)}%`);
        }
    } else {
        L.push('   (Skipped: No training snapshots provided for offline check)');
    }

    // --- 5. BASIS BOUNDARY LEAKAGE ---
    L.push('\n5. BASIS BOUNDARY LEAKAGE');
    let maxLeak = 0;
    const Phi = romEngine.Phi;
    // Heuristic for cantilever beam root
    const nV = patch.controlPoints[0].length;
    for (let j = 0; j < nV; j++) {
        const dofX = (0 * nV + j) * 2;
        const dofY = (0 * nV + j) * 2 + 1;
        for (let mode = 0; mode < k; mode++) {
            maxLeak = Math.max(maxLeak, Math.abs(Phi.get(dofX, mode)));
            maxLeak = Math.max(maxLeak, Math.abs(Phi.get(dofY, mode)));
        }
    }
    L.push(`   Max Boundary Displacement: ${maxLeak.toExponential(4)}`);
    if (maxLeak > 1e-12) L.push('   VERDICT: ❌ FAIL — Basis violates constraints!');
    else L.push('   VERDICT: ✅ PASS — Basis is strictly constrained.');

    L.push('\n' + hr);
    const report = L.join('\n');
    console.log(report);
    return report;
};
