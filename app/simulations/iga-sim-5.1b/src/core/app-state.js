/**
 * Phase 5.1b - App State
 * The central orchestrator for the Geometric ECSW Trainer.
 */

class AppState {
    constructor() {
        window.app = this;
        // 1. Initialize Core Engines
        this.nurbs = new window.NURBS2D();
        this.factory = window.GeometryFactory;
        this.solver = new window.ECSWSolver(this.nurbs);
        this.bridge = new window.SolverBridge(this.solver);
        this.collector = new window.SnapshotCollector();
        this.reporter = new window.AuditReporter('audit-log');
        this.charts = new window.AppCharts();
        
        // 1.5 Initialize Visuals
        window.viz = new window.VisualsEngine('canvas-container', this.nurbs, this.solver);
        
        // 2. State Properties
        this.params = { R: 1.0, L1: 4.0, L2: 4.0, h: 2, p: 0, Tx: 100.0 };
        this.k = 8;
        this.isTrained = false;
        this.evaluationMode = 'fom'; // 'fom', 'rom', or 'ecsw'
        
        // 3. Current Data
        this.basePatch = null;
        this.activePatch = null;
        this.lastResult = null;
        
        // 4. ROM Data
        this.snapshots = [];
        this.basis = null;
        this.ecswData = null;

        // 5. Update Loop State
        this.isDirty = false;
        this.isSolving = false;

        this.init();
        this.startUpdateLoop();
    }

    async init() {
        console.log("[Phase 5.1b] Initializing App State...");
        if (window.viz && !window.viz.solver) window.viz.solver = this.solver;
        this.requestUpdate();
    }

    requestUpdate() {
        this.isDirty = true;
    }

    startUpdateLoop() {
        const loop = async () => {
            if (this.isDirty && !this.isSolving) {
                this.isDirty = false;
                this.isSolving = true;
                try {
                    await this.updateGeometry();
                } catch (err) {
                    console.error("[UpdateLoop] Error:", err);
                } finally {
                    this.isSolving = false;
                }
            }
            requestAnimationFrame(loop);
        };
        requestAnimationFrame(loop);
    }

    async updateGeometry() {
        if (window.dom) window.dom.setLoading("Solving...");
        await new Promise(r => setTimeout(r, 20));

        this.basePatch = this.factory.generateNotchPlate(this.params.R, this.params.L1, this.params.L2);
        this.activePatch = this.factory.refine(this.basePatch, this.params.h, this.params.p);
        
        let result;
        let t0 = performance.now();
        if (this.evaluationMode === 'rom' && this.isTrained && this.basis) {
            result = this.bridge.solveReduced(this.activePatch, this.basis, this.params.Tx);
        } else if (this.evaluationMode === 'ecsw' && this.isTrained && this.ecswData) {
            result = this.bridge.solveECSW(this.activePatch, this.basis, this.ecswData, this.params.Tx);
        } else {
            result = this.bridge.solveLinearStatic(this.activePatch, this.params.Tx);
        }
        let solveTime = performance.now() - t0;
        
        this.lastResult = result;
        
        // --- FPS & Speedup Calculation ---
        const now = performance.now();
        if (this.lastFrameTime) {
            const dt = now - this.lastFrameTime;
            this.fps = 0.9 * (this.fps || 30) + 0.1 * (1000 / dt);
        }
        this.lastFrameTime = now;

        if (window.viz) window.viz.update(this.activePatch, result.u);
        if (this.charts) this.charts.updateCurrent(this.params.R, this.params.L1);
        
        let romError = 0;
        let speedup = 1.0;
        if (this.evaluationMode !== 'fom' && this.isTrained) {
            let tStartFom = performance.now();
            const fomRes = this.bridge.solveLinearStatic(this.activePatch, this.params.Tx);
            let fomTime = performance.now() - tStartFom;
            
            speedup = fomTime / Math.max(0.01, solveTime);
            romError = this.calculateError(fomRes.u, result.u);
            
            if (this.reporter) this.reporter.logSpeedup(fomTime, solveTime, speedup);
        }

        const nDofs = this.activePatch.controlPoints.length * this.activePatch.controlPoints[0].length * 2;
        const k_eff = (this.evaluationMode !== 'fom' && this.basis) ? this.basis[0].length : nDofs;

        if (this.reporter) this.reporter.reportStats({
            nDofs: nDofs,
            nSnaps: this.collector.count,
            k: k_eff,
            error: romError,
            speedup: speedup,
            mode: this.evaluationMode,
            fps: this.fps
        });
    }

    calculateError(u_fom, u_rom) {
        let diffNorm = 0, fomNorm = 0;
        for (let i = 0; i < u_fom.length; i++) {
            const diff = u_fom[i] - u_rom[i];
            diffNorm += diff * diff;
            fomNorm += u_fom[i] * u_fom[i];
        }
        return fomNorm < 1e-18 ? 0 : Math.sqrt(diffNorm) / Math.sqrt(fomNorm);
    }

    async runSweep() {
        this.reporter.log('system', 'Starting ECSW Training Job...');
        this.collector.clear();
        this.charts.clearSamples();
        
        const R_range = [0.8, 1.4, 2.0];
        const L1_range = [4.0, 7.0, 10.0];
        const L2_range = [4.0, 7.0, 10.0];
        
        let count = 0;
        const total = R_range.length * L1_range.length * L2_range.length;

        for (const r of R_range) {
            for (const l1 of L1_range) {
                for (const l2 of L2_range) {
                    count++;
                    if (window.dom) window.dom.setLoading(`Sweeping: ${count}/${total}`);
                    const p = this.factory.generateNotchPlate(r, l1, l2);
                    const ap = this.factory.refine(p, this.params.h, this.params.p);
                    const res = this.bridge.solveLinearStatic(ap, this.params.Tx);
                    this.collector.add({ R: r, L1: l1, L2: l2 }, res.u);
                    this.charts.addSample(r, l1);
                    if (count === 1 || count === total) window.viz.update(ap, res.u);
                    if (count % 5 === 0) await new Promise(r => setTimeout(r, 5));
                }
            }
        }

        this.reporter.log('system', `Computing POD Basis...`);
        const podRes = window.PODEngine.computeBasis(this.collector.snapshots, this.k);
        this.basis = podRes.Phi;
        this.charts.updateEnergy(podRes.singularValues);

        this.reporter.log('system', `Training ECSW Hyper-reduction (tol=${this.params.ecswTol || 1e-4})...`);
        const PhiMat = new (window.mlMatrix || window.ML).Matrix(this.basis);
        this.ecswData = window.ECSWEngine.train(this.collector.snapshots, this.activePatch, PhiMat, this.solver, this.params.ecswTol || 1e-4);
        
        this.reporter.log('success', `ECSW Complete. Active Elements: ${this.ecswData.indices.length}/${this.activePatch.elements.length}`);
        
        this.isTrained = true;
        if (window.dom) {
            window.dom.setLoading("Training Complete");
            window.dom.enableROM(true);
            window.dom.enableECSW(true);
        }
    }
}

window.app = new AppState();
