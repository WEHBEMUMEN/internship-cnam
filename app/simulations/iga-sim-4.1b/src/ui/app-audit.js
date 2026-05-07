/**
 * Phase 4.1 Audit System
 * Centralized logging and benchmarking for ROM Offline Training.
 */

class OfflineAudit {
    constructor(app) {
        this.app = app;
        this.startTime = 0;
        this.phases = {
            1: { name: "Transient FOM Training", status: "Pending" },
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

    success() {
        const totalTime = ((performance.now() - this.startTime) / 1000).toFixed(2);
        console.log("%c--- AUDIT COMPLETE ---", "color: #10b981; font-weight: 900; font-size: 14px;");
        console.log(`Total Training Time: ${totalTime}s`);
    }
}

window.OfflineAudit = OfflineAudit;
