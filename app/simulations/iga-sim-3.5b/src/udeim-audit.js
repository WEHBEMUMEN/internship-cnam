/**
 * Phase 3.5b: Advanced Mathematical Audit for U-DEIM
 * Performs deep verification of the unassembled force projection logic.
 */
DEIMBenchmarkApp.prototype.auditMath = function() {
    const reportWidth = 60;
    const line = "═".repeat(reportWidth);
    const subLine = "─".repeat(reportWidth);

    console.log(`%c${line}`, "color:#3b82f6; font-weight:bold;");
    console.log("%c  ADVANCED U-DEIM MATHEMATICAL AUDIT", "color:#3b82f6; font-weight:bold; font-size:1.2em;");
    console.log(`  Generated: ${new Date().toISOString()}`);
    console.log(`%c${line}`, "color:#3b82f6; font-weight:bold;");

    if (!this.isTrained || !this.deimEngine) {
        console.error("❌ ERROR: Engine not trained.");
        return;
    }

    const { deimEngine, romEngine, solverFOM, patch } = this;
    const Phi = romEngine.Phi;
    const PhiT = Phi.transpose();
    const k = Phi.columns;
    const nDofs = Phi.rows;

    // --- SECTION 1: BASIS ORTHOGONALITY ---
    console.log("\n1. REDUCED BASIS (Phi) AUDIT");
    console.log(subLine);
    const IT = PhiT.mmul(Phi);
    let orthErr = 0;
    for(let i=0; i<k; i++) {
        for(let j=0; j<k; j++) {
            const target = (i === j) ? 1.0 : 0.0;
            orthErr = Math.max(orthErr, Math.abs(IT.get(i,j) - target));
        }
    }
    console.log(`   Orthogonality Error (||Phi^T Phi - I||_max): ${orthErr.toExponential(4)}`);
    console.log(`   VERDICT: ${orthErr < 1e-10 ? "✓ PASS" : "❌ FAIL"}`);

    // --- SECTION 2: INTERPOLATION CONDITIONING ---
    console.log("\n2. INTERPOLATION MATRIX (P^T U_f) AUDIT");
    console.log(subLine);
    // P^T U_f is the m x kf matrix of basis values at sampled DOFs
    const PtU_arr = Array.from({ length: deimEngine.m }, () => new Float64Array(deimEngine.kf));
    for (let i = 0; i < deimEngine.m; i++) {
        for (let j = 0; j < deimEngine.kf; j++) PtU_arr[i][j] = deimEngine.U_f.get(deimEngine.indices[i], j);
    }
    const PtU_mat = new window.mlMatrix.Matrix(PtU_arr.map(r => Array.from(r)));
    const svd = new window.mlMatrix.SVD(PtU_mat);
    const s = svd.diagonal;
    const cond = s[0] / Math.max(s[s.length-1], 1e-30);
    console.log(`   Condition Number: ${cond.toExponential(2)}`);
    console.log(`   Min Singular Val: ${s[s.length-1].toExponential(4)}`);
    console.log(`   VERDICT: ${cond < 1e6 ? "✓ HEALTHY" : "⚠ DANGEROUS (Poor sensor placement)"}`);

    // --- SECTION 3: ASSEMBLY CONSISTENCY (THE "GOLDEN" CHECK) ---
    console.log("\n3. ASSEMBLY CONSISTENCY AUDIT");
    console.log(subLine);
    console.log("   Checking if M * f_sampled == Phi^T * Assemble(U_f * PtU_pinv * f_sampled)");

    // Test with a random synthetic sampled force
    const f_test = new Float64Array(deimEngine.m).map(() => Math.random());
    
    // Path A: The Direct U-DEIM Projection (Offline path)
    const F_red_direct = new Float64Array(k);
    for(let i=0; i<k; i++) {
        for(let j=0; j<deimEngine.m; j++) F_red_direct[i] += deimEngine.M[i][j] * f_test[j];
    }

    // Path B: Manual Assembly (The expensive "True" path)
    const f_recon_full = deimEngine.reconstructFullForce(f_test);
    const F_red_manual = new Float64Array(k);
    for(let i=0; i<k; i++) {
        for(let d=0; d<nDofs; d++) F_red_manual[i] += PhiT.get(i, d) * f_recon_full[d];
    }

    let diffNorm = 0, manualNorm = 0;
    for(let i=0; i<k; i++) {
        const diff = F_red_direct[i] - F_red_manual[i];
        diffNorm += diff * diff;
        manualNorm += F_red_manual[i] * F_red_manual[i];
    }
    const assemblyErr = Math.sqrt(diffNorm) / Math.max(Math.sqrt(manualNorm), 1e-12);
    console.log(`   Relative Assembly Mismatch: ${(assemblyErr * 100).toFixed(6)}%`);
    console.log(`   VERDICT: ${assemblyErr < 1e-10 ? "✓ PASS (Math is perfect)" : "❌ FAIL (Projection logic bug!)"}`);

    // --- SECTION 4: SNAPSHOT RECONSTRUCTION ---
    console.log("\n4. SNAPSHOT RECONSTRUCTION AUDIT");
    console.log(subLine);
    console.log("   Verifying U-DEIM against first training snapshot...");
    
    const u_snap = this.snapDisp[0];
    const f_true_snap = solverFOM.calculateInternalForce(patch, u_snap);
    const f_sampled_snap = deimEngine.calculateSampledInternalForce(solverFOM, patch, u_snap);
    const f_recon_snap = deimEngine.reconstructFullForce(f_sampled_snap);

    let rNum = 0, rDen = 0;
    for(let i=0; i<nDofs; i++) {
        rNum += (f_true_snap[i] - f_recon_snap[i])**2;
        rDen += f_true_snap[i]**2;
    }
    const snapErr = Math.sqrt(rNum / Math.max(rDen, 1e-30));
    console.log(`   Snapshot Force Error: ${(snapErr * 100).toFixed(4)}%`);
    console.log(`   VERDICT: ${snapErr < 0.05 ? "✓ PASS" : "⚠ WARNING (Lossy basis)"}`);

    console.log(`\n%c${line}`, "color:#3b82f6; font-weight:bold;");
};
