/**
 * Phase 5.1a - App State
 * The central orchestrator for the Geometric Offline Trainer.
 */

class AppState {
    constructor() {
        window.app = this;
        // 1. Initialize Core Engines
        this.nurbs = new window.NURBS2D();
        this.factory = window.GeometryFactory;
        this.solver = new window.IGA2DSolver(this.nurbs);
        this.bridge = new window.SolverBridge(this.solver);
        this.collector = new window.SnapshotCollector();
        this.reporter = new window.AuditReporter('audit-log');
        this.charts = new window.AppCharts();
        
        // 1.5 Initialize Visuals
        window.viz = new window.VisualsEngine('canvas-container', this.nurbs, this.solver);
        
        // 2. State Properties
        this.params = { R: 1.0, L1: 4.0, L2: 4.0, h: 2, p: 0 };
        this.k = 5;
        this.isTrained = false;
        this.evaluationMode = 'fom'; // 'fom' or 'rom'
        
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
        // Patch visuals if it exists but lacks solver
        if (window.viz && !window.viz.solver) window.viz.solver = this.solver;
        this.updateGeometry();
    }

    /**
     * Updates the base geometry and solves the system
     */
    async updateGeometry() {
        if (window.dom) window.dom.setLoading("Solving...");
        
        // Yield to allow UI to show "Solving..."
        await new Promise(r => setTimeout(r, 20));

        // Generate base 2-element notch plate
        this.basePatch = this.factory.generateNotchPlate(this.params.R, this.params.L1, this.params.L2);
        
        // Refine for solver accuracy (h-refinement and p-elevation)
        this.activePatch = this.factory.refine(this.basePatch, this.params.h, this.params.p);
        
        // Solve linear static (Check mode)
        let result;
        if (this.evaluationMode === 'rom' && this.isTrained && this.basis) {
            console.log("[AppState] Solving in REDUCED space...");
            result = this.bridge.solveReduced(this.activePatch, this.basis, 100.0);
        } else {
            console.log("[AppState] Solving in FULL space...");
            result = this.bridge.solveLinearStatic(this.activePatch, 100.0);
        }
        
        this.lastResult = result;

        // Update UI/Visuals if they exist
        if (window.viz) window.viz.update(this.activePatch, result.u);
        if (this.charts) this.charts.updateCurrent(this.params.R, this.params.L1);
        if (this.reporter) this.reporter.reportStats({
            nDofs: this.activePatch.controlPoints.length * this.activePatch.controlPoints[0].length * 2,
            nSnaps: this.collector.count
        });
    }

    /**
     * Runs the Full Geometric Training Job
     */
    async runSweep() {
        this.reporter.log('system', 'Starting Full Geometric Training Job...');
        this.collector.clear();
        this.charts.clearSamples();
        
        // 1. Snapshot Harvesting (Sweep)
        const R_range = [0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0];
        const L_range = [4.0, 5.0, 6.0, 7.0];
        
        let count = 0;
        const total = R_range.length * L_range.length;

        for (const r of R_range) {
            for (const l of L_range) {
                count++;
                if (window.dom) window.dom.setLoading(`Sweeping: ${count}/${total}`);
                
                const p = this.factory.generateNotchPlate(r, l, l);
                const ap = this.factory.refine(p, this.params.h, this.params.p);
                const res = this.bridge.solveLinearStatic(ap, 100.0);
                
                this.collector.add({ R: r, L1: l, L2: l }, res.u);
                this.charts.addSample(r, l);

                // Live update visual for the first and last
                if (count === 1 || count === total) window.viz.update(ap, res.u);
                
                await new Promise(r => setTimeout(r, 5));
            }
        }

        this.reporter.log('success', `Harvested ${this.collector.count} snapshots.`);

        // 2. Basis Reduction (POD)
        console.log("[AppState] Starting POD calculation...");
        this.reporter.log('system', `Computing POD Basis (k=${this.k})...`);
        try {
            const result = window.PODEngine.computeBasis(this.collector.snapshots, this.k);
            console.log("[AppState] POD Result:", result);
            this.basis = result.Phi;

            this.charts.updateEnergy(result.singularValues);
            const finalEnergy = result.energy[result.energy.length - 1];
            this.reporter.log('success', `Reduced Basis Extracted. Energy: ${(finalEnergy * 100).toFixed(6)}%`);
        } catch (err) {
            console.error("[AppState] POD Error:", err);
            this.reporter.log('warning', `POD Error: ${err.message}`);
        }

        this.isTrained = true;
        if (window.dom) {
            window.dom.setLoading("Training Complete");
            window.dom.enableExport(true);
            window.dom.enableROM(true);
        }
    }
}

// Global instance
window.app = new AppState();
