/**
 * DEIM Online Audit Module (Advanced Reporter)
 * Generates a detailed mathematical breakdown of the online hyper-reduction performance.
 */
DEIMBenchmarkApp.prototype.runOnlineAudit = async function() {
    const reportWidth = 60;
    const line = "═".repeat(reportWidth);
    const subLine = "─".repeat(reportWidth);

    console.log(`%c${line}`, "color:#ec4899; font-weight:bold;");
    console.log("%c  DEIM ONLINE SOLVER AUDIT", "color:#ec4899; font-weight:bold; font-size:1.2em;");
    console.log(`  Generated: ${new Date().toISOString()}`);
    console.log(`%c${line}`, "color:#ec4899; font-weight:bold;");

    if (!this.isTrained) {
        console.error("❌ ERROR: DEIM not trained. Please run 'Train DEIM' first.");
        return;
    }

    const loadMag = this.loadMag;
    console.log("\n0. CONFIGURATION");
    console.log(`   Target Load    : ${loadMag} N (${this.loadType})`);
    console.log(`   Mesh Level (h) : ${this.meshLevel}`);
    console.log(`   Basis Modes (k): ${this.k}`);
    console.log(`   Sensors (m)    : ${this.deimM}`);
    console.log(`   N (total DOFs) : ${this.patch.controlPoints.length * this.patch.controlPoints[0].length * 2}`);

    // Helper: extract final residual norm from history
    const getFinalNorm = (hist) => {
        if (!hist || hist.length === 0) return NaN;
        const last = hist[hist.length - 1];
        return typeof last === 'number' ? last : (last.norm !== undefined ? last.norm : NaN);
    };

    // --- STEP 1: CONVERGENCE AUDIT ---
    console.log(`\n1. CONVERGENCE AUDIT`);
    console.log(subLine);
    
    const fomWrap = this.solve('fom', loadMag);
    const fom = fomWrap.result;
    const dt_fom = fomWrap.meta.time;

    const deimWrap = this.solve('deim', loadMag);
    const deim = deimWrap.result;
    const dt_deim = deimWrap.meta.time;

    const fomFinalRes = getFinalNorm(fom.residualHistory);
    const deimFinalRes = getFinalNorm(deim.residualHistory);

    console.log(`   Method   | Iterations | Final Residual | Time (ms)`);
    console.log(`   ─────────┼────────────┼────────────────┼──────────`);
    console.log(`   FOM      | ${fom.residualHistory.length.toString().padEnd(10)} | ${fomFinalRes.toExponential(3).padEnd(14)} | ${dt_fom.toFixed(1)}`);
    console.log(`   DEIM     | ${deim.residualHistory.length.toString().padEnd(10)} | ${deimFinalRes.toExponential(3).padEnd(14)} | ${dt_deim.toFixed(1)}`);

    // Show full convergence trace for DEIM
    console.log(`\n   DEIM Convergence Trace:`);
    console.log(`   Iter | Residual Norm`);
    console.log(`   ─────┼───────────────`);
    deim.residualHistory.forEach((entry, idx) => {
        const norm = typeof entry === 'number' ? entry : entry.norm;
        console.log(`   ${idx.toString().padEnd(4)} | ${norm.toExponential(6)}`);
    });

    const iterMatch = fom.residualHistory.length === deim.residualHistory.length;
    console.log(`\n   VERDICT: ${iterMatch ? "✓ PASS" : "⚠ NOTE"} — DEIM used ${deim.residualHistory.length} iterations vs FOM's ${fom.residualHistory.length}.`);

    // --- STEP 2: DISPLACEMENT ACCURACY ---
    console.log(`\n2. DISPLACEMENT ACCURACY (State Error)`);
    console.log(subLine);
    
    if (!deim.u || !fom.u) {
        console.error("   ❌ DEIM solver returned null displacement. Likely diverged.");
        return;
    }

    let num = 0, den = 0, maxErr = 0, maxErrDof = -1, leak = 0, leakDof = -1;
    const bcs = this.getBCs();
    const fixedDofs = bcs.filter(b => b.type === 'dirichlet').map(b => b.dof);

    for (let i = 0; i < fom.u.length; i++) {
        const err = Math.abs(fom.u[i] - deim.u[i]);
        num += err**2;
        den += fom.u[i]**2;
        if (err > maxErr) { maxErr = err; maxErrDof = i; }
        
        if (fixedDofs.includes(i)) {
            if (Math.abs(deim.u[i]) > leak) { leak = Math.abs(deim.u[i]); leakDof = i; }
        }
    }
    const relError = Math.sqrt(num / Math.max(den, 1e-30));

    // Tip displacement comparison
    const nU = this.patch.controlPoints.length, nV = this.patch.controlPoints[0].length;
    const tipIdx = ((nU - 1) * nV + Math.floor(nV / 2)) * 2 + 1;
    const tipFom = fom.u[tipIdx];
    const tipDeim = deim.u[tipIdx];

    console.log(`   L2 Relative Error : ${(relError * 100).toFixed(6)}%`);
    console.log(`   Max Point Error   : ${maxErr.toExponential(4)} (DOF ${maxErrDof})`);
    console.log(`   Tip Displacement  : FOM=${tipFom.toFixed(6)}, DEIM=${tipDeim.toFixed(6)}, Δ=${(tipFom - tipDeim).toExponential(4)}`);
    console.log(`   Boundary Leakage  : ${leak.toExponential(4)} (DOF ${leakDof})`);

    const leakPass = leak < 1e-10;
    const accPass = relError < 0.01;
    console.log(`\n   ACCURACY:  ${accPass ? "✓ PASS" : "❌ FAIL"} — ${accPass ? "Healthy (<1%)" : "Poor (>1%)"}.`);
    console.log(`   BOUNDARY:  ${leakPass ? "✓ PASS" : "❌ FAIL"} — ${leakPass ? "Strictly enforced" : `Leaking at DOF ${leakDof}`}.`);

    // --- STEP 3: PHYSICAL FORCE AUDIT ---
    console.log(`\n3. PHYSICAL FORCE AUDIT (Hyper-Reduction Quality)`);
    console.log(subLine);
    
    // 1. Get the exact internal force at the DEIM converged state
    const f_true = this.solverFOM.calculateInternalForce(this.patch, deim.u);
    // FIX: To be fair, zero out the constrained DOFs in the true force,
    // since the DEIM basis U_f was trained on 'clean' forces.
    fixedDofs.forEach(d => f_true[d] = 0);
    
    // 2. Get the DEIM-sampled force at the same state
    const f_sampled = this.deimEngine.calculateSampledInternalForce(this.solverFOM, this.patch, deim.u);
    // calculateSampledInternalForce already zeroes out its internal buffer before extracting indices.
    
    // 3. Reconstruct the full N-dimensional force from those M samples
    const f_deim_full = this.deimEngine.reconstructFullForce(f_sampled);

    if (f_deim_full.length === 0) {
        console.warn("   ⚠ reconstructFullForce returned empty. Skipping force audit.");
    } else {
        let fNum = 0, fDen = 0, fMaxErr = 0, fMaxDof = -1;
        for (let i = 0; i < f_true.length; i++) {
            if (fixedDofs.includes(i)) continue;
            const err = Math.abs(f_true[i] - f_deim_full[i]);
            fNum += err**2;
            fDen += f_true[i]**2;
            if (err > fMaxErr) { fMaxErr = err; fMaxDof = i; }
        }
        const forceRecError = Math.sqrt(fNum / Math.max(fDen, 1e-30));

        console.log(`   Force Reconstruction Error: ${(forceRecError * 100).toFixed(6)}%`);
        console.log(`   Max Force Error at DOF    : ${fMaxDof} (|Δf| = ${fMaxErr.toExponential(4)})`);
        console.log(`   ||f_true|| = ${Math.sqrt(fDen).toExponential(4)}, ||f_deim|| = ${Math.sqrt(f_deim_full.reduce((s,v,i) => fixedDofs.includes(i) ? s : s + v*v, 0)).toExponential(4)}`);

        const forcePass = forceRecError < 0.05;
        console.log(`\n   VERDICT: ${forcePass ? "✓ PASS" : "⚠ WARNING"} — ${forcePass ? "Physically consistent" : "High interpolation error"}.`);
        
        // 4. NEW: Projection Consistency Check (Is Φ^T F_true ≈ M_deim * f_sampled?)
        const Phi = this.romEngine.Phi;
        const PhiT = Phi.transpose();
        const k_safe = Math.min(this.k, Phi.columns);
        const f_proj_exact = new Float64Array(k_safe);
        for(let i=0; i<k_safe; i++) {
            for(let d=0; d<f_true.length; d++) f_proj_exact[i] += PhiT.get(i, d) * f_true[d];
        }

        // DEIM reduced force using precomputed M_deim matrix (exactly as solveReduced does)
        const U_f = this.deimEngine.U_f;
        const PtU_pinv = this.deimEngine.PtU_pinv;
        const m = this.deimEngine.m;
        const f_proj_deim = new Float64Array(this.k);
        
        // Reconstruct c first (coefficients in U_f space)
        const c = new Float64Array(U_f.columns);
        for(let i=0; i<U_f.columns; i++) {
            for(let j=0; j<m; j++) {
                const val = f_sampled[j];
                if (isFinite(val)) c[i] += PtU_pinv[i][j] * val;
            }
        }
        // Project U_f * c
        const PhiT_Uf = PhiT.mmul(U_f);
        for(let i=0; i<this.k; i++) {
            for(let j=0; j<U_f.columns; j++) f_proj_deim[i] += PhiT_Uf.get(i, j) * c[j];
        }

        let pNum = 0, pDen = 0;
        for(let i=0; i<this.k; i++) {
            const err = f_proj_exact[i] - f_proj_deim[i];
            pNum += err * err;
            pDen += f_proj_exact[i] * f_proj_exact[i];
        }
        const projError = Math.sqrt(pNum / Math.max(pDen, 1e-30));
        console.log(`\n   DEIM Projection Error : ${(projError * 100).toFixed(6)}%`);
        console.log(`   (Compares ΦᵀF_true vs M_deim * f_sampled)`);

        // 5. DIAGNOSTIC: Direct Sample Comparison
        console.log(`\n   DIAGNOSTIC: First 5 Samples`);
        console.log(`   Idx | DOF | FOM Value | Sampled Value | Match?`);
        console.log(`   ────┼─────┼───────────┼───────────────┼───────`);
        for(let i=0; i<Math.min(5, m); i++) {
            const dof = this.deimEngine.indices[i];
            const v1 = f_true[dof], v2 = f_sampled[i];
            const match = Math.abs(v1 - v2) < 1e-8 ? "YES" : "NO!!";
            console.log(`   ${i+1} | ${dof.toString().padStart(3)} | ${v1.toExponential(4)} | ${v2.toExponential(4)} | ${match}`);
        }
    }

    // --- STEP 4: REDUCED SPACE DIAGNOSTICS ---
    console.log(`\n4. REDUCED SPACE DIAGNOSTICS`);
    console.log(subLine);

    const Phi = this.romEngine.Phi || this.romEngine.basis;
    if (Phi) {
        // Projection error: how well does Phi * Phi^T * u_fom approximate u_fom?
        const k = Phi.columns;
        const ur_proj = new Float64Array(k);
        for (let j = 0; j < k; j++) {
            for (let i = 0; i < Phi.rows; i++) {
                ur_proj[j] += Phi.get(i, j) * fom.u[i];
            }
        }
        const u_proj = new Float64Array(Phi.rows);
        for (let i = 0; i < Phi.rows; i++) {
            for (let j = 0; j < k; j++) {
                u_proj[i] += Phi.get(i, j) * ur_proj[j];
            }
        }
        let projNum = 0, projDen = 0;
        for (let i = 0; i < fom.u.length; i++) {
            projNum += (fom.u[i] - u_proj[i])**2;
            projDen += fom.u[i]**2;
        }
        const projError = Math.sqrt(projNum / Math.max(projDen, 1e-30));
        console.log(`   POD Projection Error : ${(projError * 100).toFixed(6)}% (k=${k} modes)`);
        console.log(`   (Best achievable error with this basis)`)
        
        // Check if basis is zero at boundaries
        let basisLeak = 0;
        fixedDofs.forEach(d => {
            for (let j = 0; j < k; j++) basisLeak = Math.max(basisLeak, Math.abs(Phi.get(d, j)));
        });
        console.log(`   Basis Boundary Leak  : ${basisLeak.toExponential(4)}`);
        console.log(`   VERDICT: ${basisLeak < 1e-12 ? "✓ PASS" : "❌ FAIL"} — Basis ${basisLeak < 1e-12 ? "is clean at boundaries" : "has non-zero values at fixed DOFs"}.`);
    }

    // --- STEP 5: PERFORMANCE BENCHMARK ---
    console.log(`\n5. PERFORMANCE BENCHMARK`);
    console.log(subLine);
    const speedup = dt_fom / dt_deim;
    console.log(`   FOM Assembly Time  : ${dt_fom.toFixed(1)} ms`);
    console.log(`   DEIM Solve Time    : ${dt_deim.toFixed(1)} ms`);
    console.log(`   Speedup            : ${speedup.toFixed(2)}x`);
    console.log(`   Active Elements    : ${this.deimEngine.activeElements ? this.deimEngine.activeElements.length : 'N/A'} / ${this.patch.elements ? this.patch.elements.length : 'N/A'}`);
    console.log(`   Status: ${speedup > 1.0 ? "🚀 ACCELERATED" : "🐢 OVERHEAD-DOMINATED (expected at low mesh levels)"}`);

    // --- OVERALL ---
    console.log(`\n%c${line}`, "color:#ec4899; font-weight:bold;");
    const overall = accPass && leakPass;
    console.log(`%c  OVERALL: ${overall ? "✓ ALL CHECKS PASSED" : "❌ ISSUES DETECTED — See details above"}`, `color:${overall ? "#10b981" : "#ef4444"}; font-weight:bold;`);
    console.log(`%c${line}`, "color:#ec4899; font-weight:bold;");
};
