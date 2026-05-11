/**
 * Phase 5.1a - Offline Trainer
 * Orchestrates the geometric sweep and POD extraction.
 */

class OfflineTrainer {
    constructor(app) {
        this.app = app;
    }

    async runGeometricSweep(R_range, L1_range, L2_range) {
        const snapshots = [];
        const total = R_range.length * L1_range.length * L2_range.length;
        let count = 0;

        for (const r of R_range) {
            for (const l1 of L1_range) {
                for (const l2 of L2_range) {
                    count++;
                    window.dom.log('system', `Solving shape: R=${r}, L1=${l1}, L2=${l2} (${count}/${total})`);
                    
                    const patch = window.GeometryFactory.generateNotchPlate(r, l1, l2);
                    const refined = window.GeometryFactory.refine(patch, 2, 0);
                    const res = window.app.bridge.solveLinearStatic(refined, 100.0);
                    
                    snapshots.push(res.u);
                    await new Promise(r => setTimeout(r, 20));
                }
            }
        }
        return snapshots;
    }
}

window.OfflineTrainer = OfflineTrainer;
