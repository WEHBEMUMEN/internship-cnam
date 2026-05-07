/**
 * Phase 4.0 Audit & Diagnostic Tool
 */

class Audit40 {
    static run(app) {
        let report = [];
        const hr = '========================================';
        report.push(hr);
        report.push('    PHASE 4.0 DIAGNOSTIC AUDIT    ');
        report.push(hr);

        try {
            report.push(`1. STATE INFO`);
            report.push(`   Time        : ${app.currentTime.toFixed(4)} s`);
            report.push(`   Running     : ${app.isRunning}`);
            
            const patch = app.patch;
            const nV = patch.controlPoints[0].length;
            const q = patch.q;
            const V = patch.V;

            report.push(`\n2. MESH & REFINEMENT`);
            report.push(`   Degree (p,q): ${patch.p}, ${patch.q}`);
            report.push(`   Knots U len : ${patch.U.length}`);
            report.push(`   Knots V len : ${patch.V.length}`);
            report.push(`   Control Pts : ${patch.controlPoints.length} x ${patch.controlPoints[0].length}`);
            report.push(`   Total DOFs  : ${app.dyn.dofs}`);

            report.push(`\n3. DYNAMICS ENGINE`);
            report.push(`   Alpha (M)   : ${app.dyn.alpha}`);
            report.push(`   Beta (K)    : ${app.dyn.betaR}`);
            
            // Check Mass Matrix
            let mNaN = false, mInf = false, mZero = true;
            for (let i = 0; i < app.dyn.dofs; i++) {
                for (let j = 0; j < app.dyn.dofs; j++) {
                    const val = app.dyn.M[i][j];
                    if (isNaN(val)) mNaN = true;
                    if (!isFinite(val)) mInf = true;
                    if (Math.abs(val) > 1e-12) mZero = false;
                }
            }
            report.push(`   Mass Matrix : ${mNaN ? '❌ Contains NaN' : (mInf ? '❌ Contains Inf' : (mZero ? '❌ All Zero' : '✅ OK'))}`);

            // Check State Vectors
            let uNaN = false, uInf = false;
            let vNaN = false, vInf = false;
            let aNaN = false, aInf = false;
            let uMax = -Infinity, uMin = Infinity;
            
            for (let i = 0; i < app.dyn.dofs; i++) {
                const ui = app.dyn.u[i];
                if (isNaN(ui)) uNaN = true;
                if (!isFinite(ui)) uInf = true;
                if (Math.abs(ui) > uMax) uMax = Math.abs(ui);
                if (Math.abs(ui) < uMin && Math.abs(ui) > 1e-18) uMin = Math.abs(ui);
                
                if (isNaN(app.dyn.v[i])) vNaN = true;
                if (!isFinite(app.dyn.v[i])) vInf = true;
                if (isNaN(app.dyn.a[i])) aNaN = true;
                if (!isFinite(app.dyn.a[i])) aInf = true;
            }
            report.push(`   Disp Vector : ${uNaN ? '❌ NaN' : (uInf ? '❌ Inf' : '✅ OK')}`);
            report.push(`   Max |U|      : ${(uMax * 1000).toFixed(4)} mm`);
            report.push(`   Min |U|>0    : ${(uMin * 1000).toExponential(2)} mm`);
            report.push(`   Vel Vector  : ${vNaN ? '❌ NaN' : (vInf ? '❌ Inf' : '✅ OK')}`);
            report.push(`   Acc Vector  : ${aNaN ? '❌ NaN' : (aInf ? '❌ Inf' : '✅ OK')}`);

            // NEW: Boundary Leakage Check
            let maxLeak = 0;
            for (let j = 0; j < nV; j++) {
                const idx = j * 2;
                maxLeak = Math.max(maxLeak, Math.abs(app.dyn.u[idx]), Math.abs(app.dyn.u[idx+1]));
            }
            report.push(`   Boundary Leak: ${(maxLeak * 1000).toExponential(2)} mm`);
            if (maxLeak > 1e-5) report.push(`   ❌ WARNING: BC Leakage detected! Support is moving.`);
            else report.push(`   ✅ Boundary is strictly fixed.`);

            report.push(`\n4. PHYSICAL CONSISTENCY (U-OUTCOME)`);
            // 1. Nonlinear Static Ground Truth
            const F0 = parseFloat(document.getElementById('input-f0').value) || 0;
            const bcs = [];
            for (let j = 0; j < nV; j++) bcs.push({ i: 0, j: j, axis: 'both', value: 0 });
            
            const nU_f = patch.controlPoints.length;
            const loads = [];
            for (let j = 0; j < nV; j++) {
                const weight = (V[j + q + 1] - V[j]) / (q + 1);
                loads.push({ i: nU_f - 1, j: j, fx: 0, fy: -F0 * weight, type: 'nodal' });
            }

            const staticRes = app.fom.solveNonlinear(patch, bcs, loads, { iterations: 10, tolerance: 1e-6 });
            const uStatic = staticRes.u;
            const staticTipY = Math.abs(uStatic[app.dyn.dofs - 1]);
            const currentTipY = Math.abs(app.dyn.u[app.dyn.dofs - 1]);

            // 2. Linear Reference (Simple check for sanity)
            const L = 10.0, H = 2.0, t = app.fom.thickness || 1.0, E = app.fom.E || 100000;
            const I = (t * Math.pow(H, 3)) / 12;
            const linearTheory = (F0 * Math.pow(L, 3)) / (3 * E * I);

            report.push(`   NL Static Tip (Truth): ${(staticTipY * 1000).toFixed(4)} mm`);
            report.push(`   Dynamic Tip (Current): ${(currentTipY * 1000).toFixed(4)} mm`);
            
            const nlRatio = staticTipY > 0 ? (currentTipY / staticTipY) : 0;
            report.push(`   Convergence Ratio   : ${nlRatio.toFixed(4)} (Dynamic/Static NL)`);

            if (nlRatio > 0.98 && nlRatio < 1.02) {
                report.push(`   ✅ Steady-state reached (<2% diff from Static NL).`);
            } else {
                report.push(`   ℹ️ Transient state: Oscillating around equilibrium.`);
            }

            if (staticTipY < linearTheory * 0.5 && patch.p === 1) {
                report.push(`   ⚠️ DIAGNOSIS: Massive Shear Locking (Static NL is 2x stiffer than Linear).`);
            } else if (staticTipY > linearTheory * 1.5) {
                report.push(`   ℹ️ Geometric Nonlinearity detected (Large Displacements).`);
            }

            report.push(`\n5. L2 CONVERGENCE ERROR`);
            if (F0 > 0) {
                let diff2 = 0, norm2 = 0;
                for (let i = 0; i < app.dyn.dofs; i++) {
                    diff2 += Math.pow(app.dyn.u[i] - uStatic[i], 2);
                    norm2 += Math.pow(uStatic[i], 2);
                }
                const relError = Math.sqrt(diff2) / Math.sqrt(norm2 || 1e-18);
                report.push(`   Global L2 Error  : ${(relError * 100).toFixed(4)}% (Dynamic vs Static NL)`);
            } else {
                report.push(`   (Skipped: Load is zero)`);
            }

            report.push(`\n6. INTERNAL FORCE & STIFFNESS`);
            const Fext = new Float64Array(app.dyn.dofs);
            const nU = patch.controlPoints.length;
            let sumWeights = 0;
            report.push(`   Testing Edge Loading (Right edge)`);
            for (let j = 0; j < nV; j++) {
                const dofIdx = ((nU - 1) * nV + j) * 2 + 1; 
                const weight = (V[j + q + 1] - V[j]) / (q + 1);
                sumWeights += weight;
                if (isNaN(weight)) report.push(`   ❌ Weight [${j}] is NaN (V[${j+q+1}]=${V[j+q+1]}, V[${j}]=${V[j]})`);
                else report.push(`   - Node ${j} weight: ${weight.toFixed(4)}`);
            }
            report.push(`   Total Integration Sum: ${sumWeights.toFixed(6)} (Should be 1.0000)`);
            if (Math.abs(sumWeights - 1.0) > 1e-6) {
                report.push(`   ❌ ERROR: Force distribution does not sum to 1.0!`);
            } else {
                report.push(`   ✅ Force integration is mathematically consistent.`);
            }

            // Dump to console
            const out = report.join('\n');
            console.log(out);
            
            // Try updating UI status if exists
            const statusEl = document.getElementById('run-status');
            if (statusEl) {
                statusEl.classList.remove('hidden');
                statusEl.innerHTML = `<pre style="font-size:0.6rem; color:#10b981; margin:10px 0; max-height:150px; overflow-y:auto; background:rgba(0,0,0,0.3); padding:8px; border-radius:8px;">${out}</pre>`;
            } else {
                alert("Audit printed to console. Check F12 Developer Tools.");
            }

        } catch (e) {
            report.push(`\n❌ FATAL AUDIT ERROR: ${e.message}`);
            report.push(e.stack);
            console.error("Audit Error:", e);
            alert("Audit failed! See console.");
        }
    }
}
window.Audit40 = Audit40;
