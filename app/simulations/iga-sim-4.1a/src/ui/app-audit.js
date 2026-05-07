/**
 * Phase 4.1 Audit System
 * Centralized logging and benchmarking for ROM Offline Training.
 */

class OfflineAudit {
    constructor(app) {
        this.app = app;
        this.startTime = 0;
        this.phases = {
            1: { name: "Parametric Load Sweep", status: "Pending" },
            2: { name: "POD Basis Extraction", status: "Pending" },
            3: { name: "ECSW Sparse Sampling", status: "Pending" },
            4: { name: "Reduced Matrix Assembly", status: "Pending" },
            5: { name: "Verification & Export", status: "Pending" }
        };
    }

    start() {
        this.startTime = performance.now();
        console.clear();
        console.log("%c--- ROM OFFLINE AUDIT START ---", "color: #f43f5e; font-weight: 900; font-size: 14px;");
    }

    log(phaseNum, message, data = null) {
        const phase = this.phases[phaseNum];
        const time = ((performance.now() - this.startTime) / 1000).toFixed(2);
        
        console.log(
            `%c[PHASE ${phaseNum}] %c${phase.name} %c(${time}s)`,
            "color: #f43f5e; font-weight: bold;",
            "color: #f8fafc; font-weight: 600;",
            "color: #94a3b8; font-style: italic;",
            message
        );
        
        if (data) {
            console.table(data);
        }

        // Update UI status if possible
        const statusEl = document.getElementById('train-status');
        if (statusEl) {
            statusEl.innerText = `[${phaseNum}/5] ${phase.name}: ${message}`;
        }
    }

    reportSVD(energyTrace, k) {
        const actualK = Math.min(k, energyTrace.length);
        const top5 = energyTrace.slice(0, 5).map((e, i) => ({ Mode: i+1, "Cumul. Energy": (e*100).toFixed(4) + "%" }));
        this.log(2, `Extracted ${actualK} modes. Cumulative Energy: ${(energyTrace[actualK-1]*100).toFixed(4)}%`, top5);
    }

    reportECSW(total, sampled, tol) {
        const ratio = ((sampled / total) * 100).toFixed(2);
        this.log(3, `ECSW Training Complete (Tol: ${tol}).`, {
            "Total Elements": total,
            "Sampled Elements": sampled,
            "Sparsity Ratio": ratio + "%"
        });
    }

    reportReconstruction(avgErr, maxErr) {
        this.log(2, `Basis Reconstruction Quality Check:`, {
            "Average L2 Error": (avgErr * 100).toExponential(4) + "%",
            "Maximum L2 Error": (maxErr * 100).toExponential(4) + "%",
            "Health": maxErr < 0.01 ? "✅ EXCELLENT" : (maxErr < 0.05 ? "⚠️ ACCEPTABLE" : "❌ POOR")
        });
    }

    success() {
        const totalTime = ((performance.now() - this.startTime) / 1000).toFixed(2);
        console.log("%c--- AUDIT COMPLETE ---", "color: #10b981; font-weight: 900; font-size: 14px;");
        console.log(`Total Training Time: ${totalTime}s`);
    }
}

window.OfflineAudit = OfflineAudit;

/**
 * Diagnostic Audit for Phase 4.1a (Ported from 4.0)
 */
class Audit41a {
    static run(app) {
        let report = [];
        const hr = '========================================';
        report.push(hr);
        report.push('    PHASE 4.1a DIAGNOSTIC AUDIT    ');
        report.push(hr);

        try {
            const patch = app.patch;
            const nV = patch.controlPoints[0].length;
            const q = patch.q;
            const V = patch.V;

            report.push(`1. MESH & REFINEMENT`);
            report.push(`   Degree (p,q): ${patch.p}, ${patch.q}`);
            report.push(`   Knots U len : ${patch.U.length}`);
            report.push(`   Knots V len : ${patch.V.length}`);
            report.push(`   Control Pts : ${patch.controlPoints.length} x ${nV}`);
            report.push(`   Total DOFs  : ${app.dyn.dofs}`);

            report.push(`\n2. DYNAMICS STATE`);
            let uNaN = false, uMax = 0;
            for (let i = 0; i < app.dyn.dofs; i++) {
                if (isNaN(app.dyn.u[i])) uNaN = true;
                uMax = Math.max(uMax, Math.abs(app.dyn.u[i]));
            }
            report.push(`   Disp Vector : ${uNaN ? '❌ NaN' : '✅ OK'}`);
            report.push(`   Max |U|      : ${(uMax * 1000).toFixed(4)} mm`);

            report.push(`\n3. PHYSICAL CONSISTENCY`);
            const L = 10.0, H = 2.0, t = app.fom.thickness || 1.0, E = app.fom.E || 100000;
            const I = (t * Math.pow(H, 3)) / 12;
            const F0 = 400; // Hardcoded mag in app-core for 4.1a
            const staticTheory = (F0 * Math.pow(L, 3)) / (3 * E * I);
            const currentTipY = Math.abs(app.dyn.u[app.dyn.dofs - 1]);
            const ratio = staticTheory > 0 ? (currentTipY / staticTheory) : 0;

            report.push(`   Theory (Linear): ${(staticTheory * 1000).toFixed(4)} mm`);
            report.push(`   Current Tip Y  : ${(currentTipY * 1000).toFixed(4)} mm`);
            report.push(`   Accuracy Ratio : ${ratio.toFixed(4)}`);
            
            if (patch.p === 1 && ratio < 0.2) {
                report.push(`   ⚠️ DIAGNOSIS: Massive Shear Locking (p=1).`);
            }

            report.push(`\n4. FORCE INTEGRATION`);
            let sumWeights = 0;
            for (let j = 0; j < nV; j++) {
                const weight = (V[j + q + 1] - V[j]) / (q + 1);
                sumWeights += weight;
            }
            report.push(`   Total Integration Sum: ${sumWeights.toFixed(6)}`);
            report.push(sumWeights > 0.99 && sumWeights < 1.01 ? "   ✅ Force math is perfect." : "   ❌ Force math ERROR.");

            const out = report.join('\n');
            console.log(out);
            const statusEl = document.getElementById('run-status');
            if (statusEl) {
                statusEl.classList.remove('hidden');
                statusEl.innerHTML = `<pre style="font-size:0.6rem; color:#10b981; margin:10px 0; max-height:150px; overflow-y:auto; background:rgba(0,0,0,0.3); padding:8px; border-radius:8px;">${out}</pre>`;
            }
        } catch (e) {
            console.error("Audit Error:", e);
        }
    }
}
window.Audit41a = Audit41a;
