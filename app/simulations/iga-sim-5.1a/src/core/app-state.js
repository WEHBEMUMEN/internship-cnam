/**
 * Phase 5.1a - App State
 * The central orchestrator for the Geometric Offline Trainer.
 */

class AppState {
    constructor() {
        // 1. Initialize Core Engines
        this.nurbs = new window.NURBS2D();
        this.factory = window.GeometryFactory;
        this.solver = new window.IGA2DSolver(this.nurbs);
        this.bridge = new window.SolverBridge(this.solver);
        this.collector = new window.SnapshotCollector();
        this.reporter = new window.AuditReporter('audit-log');
        
        // 2. State Properties
        this.params = { R: 1.0, L1: 4.0, L2: 4.0 };
        this.k = 5;
        this.isTrained = false;
        
        // 3. Current Data
        this.basePatch = null;
        this.activePatch = null;
        this.lastResult = null;
        
        // 4. ROM Data
        this.snapshots = [];
        this.basis = null;

        this.init();
    }

    async init() {
        console.log("[Phase 5] Initializing App State...");
        this.updateGeometry();
    }

    /**
     * Updates the base geometry and solves the system
     */
    async updateGeometry() {
        // Generate base 2-element notch plate
        this.basePatch = this.factory.generateNotchPlate(this.params.R, this.params.L1, this.params.L2);
        
        // Refine for solver accuracy (h-refinement)
        // We'll use 2 subdivisions for a 16-element model by default
        this.activePatch = this.factory.refine(this.basePatch, 2, 0);
        
        // Solve linear static
        const result = this.bridge.solveLinearStatic(this.activePatch, 100.0);
        this.lastResult = result;

        // Update UI/Visuals if they exist
        if (window.viz) window.viz.update(this.activePatch, result.u);
        if (window.app.reporter) window.app.reporter.reportStats({
            nDofs: this.activePatch.controlPoints.length * this.activePatch.controlPoints[0].length * 2,
            nSnaps: this.collector.count
        });
    }

    /**
     * Runs the Geometric Sweep
     */
    async runSweep() {
        console.log("[Phase 5] Starting Geometric Sweep...");
        this.snapshots = [];
        
        const R_range = [0.8, 1.0, 1.2, 1.5];
        const L_range = [4.0, 5.0, 6.0];
        
        let count = 0;
        const total = R_range.length * L_range.length;

        for (const r of R_range) {
            for (const l of L_range) {
                count++;
                if (window.dom) window.dom.setLoading(`Sweeping: ${count}/${total}`);
                
                const p = this.factory.generateNotchPlate(r, l, l);
                const ap = this.factory.refine(p, 2, 0);
                const res = this.bridge.solveLinearStatic(ap, 100.0);
                
                this.snapshots.push(res.u);
                
                // Yield to UI
                await new Promise(r => setTimeout(r, 10));
            }
        }

        console.log(`[Phase 5] Sweep complete. Collected ${this.snapshots.length} snapshots.`);
        this.isTrained = true;
        if (window.dom) {
            window.dom.setLoading("Sweep Complete");
            window.dom.enableExport(true);
        }
    }
}

// Global instance
window.app = new AppState();
