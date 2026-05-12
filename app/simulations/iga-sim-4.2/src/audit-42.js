/**
 * Phase 4.2 Audit & Diagnostic Tool (CSM Benchmark)
 */

class Audit42 {
    static run(app) {
        let report = [];
        const hr = '========================================';
        report.push(hr);
        report.push('    PHASE 4.2 CSM BENCHMARK AUDIT    ');
        report.push(hr);

        try {
            report.push(`1. STATE INFO`);
            report.push(`   Time        : ${app.currentTime.toFixed(4)} s`);
            report.push(`   Running     : ${app.isRunning}`);
            
            const patch = app.patch;
            report.push(`\n2. MESH & REFINEMENT`);
            report.push(`   Degree (p,q): ${patch.p}, ${patch.q}`);
            report.push(`   Control Pts : ${patch.controlPoints.length} x ${patch.controlPoints[0].length}`);
            report.push(`   Total DOFs  : ${app.dyn.dofs}`);

            report.push(`\n3. DYNAMICS ENGINE`);
            report.push(`   Alpha (M)   : ${app.dyn.alpha}`);
            report.push(`   Beta (K)    : ${app.dyn.betaR}`);
            
            // Check State Vectors
            let uNaN = false;
            let uMax = -Infinity;
            for (let i = 0; i < app.dyn.dofs; i++) {
                const ui = app.dyn.u[i];
                if (isNaN(ui)) uNaN = true;
                if (Math.abs(ui) > uMax) uMax = Math.abs(ui);
            }
            report.push(`   Disp Vector : ${uNaN ? '❌ NaN' : '✅ OK'}`);
            report.push(`   Max |U|      : ${(uMax * 1000).toFixed(4)} mm`);

            report.push(`\n4. BENCHMARK VERIFICATION`);
            const gy = parseFloat(document.getElementById('input-gravity').value) || 0;
            const currentTipY = Math.abs(app.dyn.u[app.dyn.dofs - 1]) * 1000; // in mm
            const targetY = Math.abs(app.targetY);

            report.push(`   Benchmark Target: ${targetY.toFixed(3)} mm`);
            report.push(`   Current Results : ${currentTipY.toFixed(3)} mm`);
            
            const error = Math.abs((currentTipY - targetY) / targetY) * 100;
            report.push(`   Relative Error  : ${error.toFixed(4)} %`);

            if (error < 0.5) report.push(`   ✅ EXCELLENT ACCURACY (<0.5%)`);
            else if (error < 2.0) report.push(`   ✅ GOOD ACCURACY (<2.0%)`);
            else if (error < 10.0) report.push(`   ⚠️ MODERATE ERROR: Consider p-refinement.`);
            else report.push(`   ❌ LARGE ERROR: Check mesh density and material props.`);

            // Skip Static Truth solve for simplicity in audit (too slow to run two simulations)
            report.push(`\n5. PHYSICAL CONSTANTS (Turek CSM)`);
            report.push(`   Length (L)   : 0.35 m`);
            report.push(`   Height (H)   : 0.02 m`);
            report.push(`   E-Modulus    : ${app.fom.E.toExponential(2)} Pa`);
            report.push(`   Poisson (v)  : ${app.fom.nu.toFixed(2)}`);
            report.push(`   Density (rho): ${app.rho} kg/m³`);

            // Dump to console
            const out = report.join('\n');
            console.log(out);
            
            const statusEl = document.getElementById('run-status');
            if (statusEl) {
                statusEl.classList.remove('hidden');
                statusEl.innerHTML = `<pre style="font-size:0.6rem; color:#f59e0b; margin:10px 0; max-height:200px; overflow-y:auto; background:rgba(0,0,0,0.3); padding:10px; border-radius:12px; border:1px solid rgba(245,158,11,0.2);">${out}</pre>`;
            } else {
                alert("Audit printed to console.");
            }

        } catch (e) {
            report.push(`\n❌ FATAL AUDIT ERROR: ${e.message}`);
            console.error("Audit Error:", e);
        }
    }
}
window.Audit42 = Audit42;
