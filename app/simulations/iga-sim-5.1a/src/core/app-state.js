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
        this.params = { R: 1.0, L1: 4.0, L2: 4.0, h: 2, p: 0, Tx: 100.0 };
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
            result = this.bridge.solveReduced(this.activePatch, this.basis, this.params.Tx);
        } else {
            console.log("[AppState] Solving in FULL space...");
            result = this.bridge.solveLinearStatic(this.activePatch, this.params.Tx);
        }
        
        this.lastResult = result;

        // Update UI/Visuals if they exist
        if (window.viz) window.viz.update(this.activePatch, result.u);
        if (this.charts) this.charts.updateCurrent(this.params.R, this.params.L1);
        
        // Calculate Error if in ROM mode
        let romError = 0;
        if (this.evaluationMode === 'rom' && this.isTrained) {
            const fomRes = this.bridge.solveLinearStatic(this.activePatch, this.params.Tx);
            romError = this.calculateError(fomRes.u, result.u);
            const errorStr = (romError * 100).toFixed(4);
            console.log(`%c[AppState] ROM Update | R: ${this.params.R.toFixed(2)}, L1: ${this.params.L1.toFixed(2)}, L2: ${this.params.L2.toFixed(2)} | Error: ${errorStr}%`, "color: #f59e0b; font-weight: bold; background: rgba(245, 158, 11, 0.1); padding: 2px 5px; border-radius: 3px;");
        } else {
            console.log(`%c[AppState] FOM Update | R: ${this.params.R.toFixed(2)}, L1: ${this.params.L1.toFixed(2)}, L2: ${this.params.L2.toFixed(2)}`, "color: #8b5cf6; font-weight: bold;");
        }

        const nDofs = this.activePatch.controlPoints.length * this.activePatch.controlPoints[0].length * 2;
        const k = (this.evaluationMode === 'rom' && this.basis) ? this.basis[0].length : nDofs;

        if (this.reporter) this.reporter.reportStats({
            nDofs: nDofs,
            nSnaps: this.collector.count,
            k: k,
            error: romError,
            mode: this.evaluationMode
        });
    }

    /**
     * Calculates Relative L2 Error between two vectors
     */
    calculateError(u_fom, u_rom) {
        let diffNorm = 0;
        let fomNorm = 0;
        for (let i = 0; i < u_fom.length; i++) {
            const diff = u_fom[i] - u_rom[i];
            diffNorm += diff * diff;
            fomNorm += u_fom[i] * u_fom[i];
        }
        if (fomNorm < 1e-18) return 0;
        return Math.sqrt(diffNorm) / Math.sqrt(fomNorm);
    }

    /**
     * Runs the Full Geometric Training Job (3D Parameter Space)
     */
    async runSweep() {
        this.reporter.log('system', 'Starting Full Geometric Training Job...');
        this.collector.clear();
        this.charts.clearSamples();
        
        // 1. Snapshot Harvesting (3D Sweep: R x L1 x L2)
        const R_range = [0.8, 1.2, 1.6, 2.0];
        const L1_range = [4.0, 6.0, 8.0, 10.0];
        const L2_range = [4.0, 6.0, 8.0, 10.0];
        
        let count = 0;
        const total = R_range.length * L1_range.length * L2_range.length;
        this.reporter.log('system', `Harvesting ${total} snapshots across 3D space...`);

        for (const r of R_range) {
            for (const l1 of L1_range) {
                for (const l2 of L2_range) {
                    count++;
                    if (window.dom) window.dom.setLoading(`Sweeping: ${count}/${total}`);
                    
                    const p = this.factory.generateNotchPlate(r, l1, l2);
                    const ap = this.factory.refine(p, this.params.h, this.params.p);
                    const res = this.bridge.solveLinearStatic(ap, this.params.Tx);
                    
                    this.collector.add({ R: r, L1: l1, L2: l2 }, res.u);
                    this.charts.addSample(r, l1); // Visualizing R vs L1 map

                    // Live update visual for the first, middle and last
                    if (count === 1 || count === Math.floor(total/2) || count === total) {
                        window.viz.update(ap, res.u);
                    }
                    
                    // Small yield for UI responsiveness
                    if (count % 4 === 0) await new Promise(r => setTimeout(r, 5));
                }
            }
        }

        this.reporter.log('success', `Harvested ${this.collector.count} snapshots.`);

        // 2. Basis Reduction (POD)
        const k_target = 8; // Higher k to capture independent L1/L2 variations
        console.log(`[AppState] Starting POD calculation (k=${k_target})...`);
        this.reporter.log('system', `Computing POD Basis (k=${k_target})...`);
        try {
            const result = window.PODEngine.computeBasis(this.collector.snapshots, k_target);
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
