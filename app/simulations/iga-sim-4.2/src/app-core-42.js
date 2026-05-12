/**
 * Phase 4.0 Application Core — Nonlinear Transient FOM
 */

class TransientLab {
    constructor() {
        this.canvas = document.getElementById('canvas-container');
        this.patch = null;
        this.fom = null; // Nonlinear Static Solver (for K, Fint)
        this.dyn = null; // Dynamics Solver (for M, C, Time Loop)
        this.viz = null; // 3D Visuals
        
        this.isRunning = false;
        this.currentTime = 0;
        this.trace = []; // {t, disp}
        this.h = 0; // Refinement level
        
        // Material State (CSM1 Defaults)
        this.E = 1.4e6;
        this.nu = 0.4;
        this.rho = 1100;
        this.gy = -2.0;

        this.targetY = -7.187; // mm
        this.init();
    }

    async init() {
        const preset = window.NURBSPresets.generateRectangle(0.35, 0.02);
        preset.constraints = [{ type: 'clamp', side: 'left' }];
        this.engine = new window.NURBS2D();
        this.patch = preset;
        this.h = 2; // Default to level 2

        // Apply initial refinement
        window.RefineUtils.apply(this.engine, this.patch, { h: this.h, p: 2 });

        this.fom = new window.IGANonlinearSolver(this.engine);
        this.fom.E = this.E;
        this.fom.nu = this.nu;
        this.fom.analysisType = 'plane-strain';

        this.dyn = new DynamicsSolver(this.patch, this.fom);
        this.dyn.rho = this.rho;
        this.dyn.assembleMass();
        this.viz = new TransientVisuals(this);
        this.updateTrace(0, 0);
        this.initCharts();
        setTimeout(() => this.runSimulation(), 1000);
    }

    updateTrace(t, disp) {
        this.viz.updateMesh(this.dyn.u, t);
    }

    initCharts() {
        const ctx = document.getElementById('chart-full-trace').getContext('2d');
        this.mainChart = new Chart(ctx, {
            type: 'line',
            data: { labels: [], datasets: [{ label: 'Tip Displacement (mm)', data: [], borderColor: '#0ea5e9', backgroundColor: 'rgba(14,165,233,0.1)', fill: true, tension: 0.4, borderWidth: 2, pointRadius: 0 }] },
            options: { 
                responsive: true, 
                maintainAspectRatio: false, 
                onClick: (e) => this.handleChartClick(e),
                scales: { 
                    x: { display: true, title: { display: true, text: 'Time (s)', color: '#94a3b8' }, grid: { color: 'rgba(148,163,184,0.1)' } }, 
                    y: { display: true, title: { display: true, text: 'Disp (mm)', color: '#94a3b8' }, grid: { color: 'rgba(148,163,184,0.1)' } } 
                }, 
                plugins: { 
                    legend: { display: false },
                    zoom: {
                        zoom: { wheel: { enabled: true }, pinch: { enabled: true }, mode: 'x' },
                        pan: { enabled: true, mode: 'x', modifierKey: 'ctrl' }
                    }
                } 
            }
        });

        const ctxRes = document.getElementById('chart-response').getContext('2d');
        this.responseChart = new Chart(ctxRes, {
            type: 'line',
            data: { labels: [], datasets: [{ label: 'Force (N)', data: [], borderColor: '#10b981', backgroundColor: 'rgba(16,185,129,0.1)', fill: true, borderWidth: 2, pointRadius: 0, tension: 0.4 }] },
            options: { responsive: true, maintainAspectRatio: false, plugins: { legend: { display: false } }, scales: { x: { display: false }, y: { grid: { color: 'rgba(16,185,129,0.1)' }, ticks: { font: { size: 9 } } } }, animation: false }
        });
    }

    handleChartClick(event) {
        if (this.isRunning || event.ctrlKey) return; // Ignore if running or Ctrl is held (panning)
        const points = this.mainChart.getElementsAtEventForMode(event, 'nearest', { intersect: false }, true);
        if (points.length) {
            const index = points[0].index;
            const traceOffset = this.trace.length - this.mainChart.data.labels.length;
            const historyIndex = traceOffset + index;
            if (this.trace[historyIndex]) {
                const state = this.trace[historyIndex];
                this.viz.updateMesh(state.u);
                document.getElementById('curr-time').innerText = state.t.toFixed(3) + 's (History)';
                document.getElementById('tip-disp').innerText = (state.disp * 1000).toFixed(2) + 'mm';
            }
        }
    }

    async runSimulation() {
        if (this.isRunning) return;
        this.isRunning = true;
        
        if (this.currentTime === 0) {
            this.dyn.u.fill(0);
            this.dyn.v.fill(0);
            this.dyn.a.fill(0);
            this.trace = [];
            this.mainChart.data.labels = [];
            this.mainChart.data.datasets[0].data = [];
            this.responseChart.data.labels = [];
            this.responseChart.data.datasets[0].data = [];
        }
        
        const timeScaleEl = document.getElementById('input-timescale');
        const timeScale = timeScaleEl ? (parseFloat(timeScaleEl.value) || 1.0) : 1.0;
        const startTime = performance.now() - (this.currentTime / timeScale) * 1000;
        
        while (this.isRunning) {
            const gy = parseFloat(document.getElementById('input-gravity').value);
            this.dyn.alpha = parseFloat(document.getElementById('input-alpha').value);
            this.dyn.betaR = parseFloat(document.getElementById('input-beta').value);
            const dt = 0.005; // Force dt for benchmark stability

            // 1. Calculate Body Forces (Gravity)
            const Fext = this.fom.calculateBodyForces(this.patch, this.rho, 0, gy);

            // 2. Solve Step
            const result = await this.dyn.solveStep(dt, Fext);
            this.currentTime += dt;
            
            const tipDisp = this.dyn.u[this.dyn.dofs - 1];
            this.trace.push({ t: this.currentTime, disp: tipDisp, u: new Float64Array(this.dyn.u) });
            if (this.trace.length > 10000) this.trace.shift();
            
            const currentTimeScaleEl = document.getElementById('input-timescale');
            const currentTimeScale = currentTimeScaleEl ? (parseFloat(currentTimeScaleEl.value) || 1.0) : 1.0;
            const elapsedWallTime = (performance.now() - startTime) / 1000;
            const targetWallTime = this.currentTime / currentTimeScale;
            
            if (targetWallTime > elapsedWallTime) {
                const delay = (targetWallTime - elapsedWallTime) * 1000;
                if (delay > 1) await new Promise(r => setTimeout(r, Math.min(delay, 100))); 
            }

            // Dynamic UI Skip: Render less frequently as timescale increases to save CPU
            const skipRate = Math.max(1, Math.floor(timeScale * 2));
            if (Math.round(this.currentTime / dt) % skipRate === 0) {
                this.updateUI(result, tipDisp, gy);
                this.viz.updateMesh(this.dyn.u);
                await new Promise(requestAnimationFrame);
            }
        }
    }

    updateUI(result, tipDisp, gravity) {
        document.getElementById('curr-time').innerText = this.currentTime.toFixed(3) + 's';
        document.getElementById('iters-val').innerText = result.iters;
        document.getElementById('tip-disp').innerText = (tipDisp * 1000).toFixed(3) + 'mm';
        
        // Update Error
        const currentY = tipDisp * 1000;
        const error = Math.abs((currentY - this.targetY) / this.targetY) * 100;
        document.getElementById('error-val').innerText = error.toFixed(2) + '%';

        // Frequency Estimation for CSM3
        if (this.trace.length > 500) {
            const freq = this.calculateFrequency();
            const freqEl = document.getElementById('freq-val');
            if (freqEl) freqEl.innerText = freq.toFixed(3) + ' Hz';
        }

        const timeStr = this.currentTime.toFixed(2);
        this.mainChart.data.labels.push(timeStr);
        this.mainChart.data.datasets[0].data.push(tipDisp * 1000);
        
        this.responseChart.data.labels.push(timeStr);
        this.responseChart.data.datasets[0].data.push(gravity);

        if (this.mainChart.data.labels.length > 5000) {
            this.mainChart.data.labels.shift();
            this.mainChart.data.datasets[0].data.shift();
        }
        if (this.responseChart.data.labels.length > 40) {
            this.responseChart.data.labels.shift();
            this.responseChart.data.datasets[0].data.shift();
        }
        
        this.mainChart.update('none');
        this.responseChart.update('none');
    }

    async updateRefinement(h, p = null) {
        this.isRunning = false;
        this.h = h;
        if (p !== null) this.p = p;
        const preset = window.NURBSPresets.generateRectangle(0.35, 0.02);
        preset.constraints = [{ type: 'clamp', side: 'left' }];
        this.patch = preset;
        window.RefineUtils.apply(this.engine, this.patch, { h: this.h, p: this.p || this.patch.p });
        this.fom = new window.IGANonlinearSolver(this.engine);
        this.fom.E = this.E;
        this.fom.nu = this.nu;
        this.fom.analysisType = 'plane-strain';

        this.dyn = new DynamicsSolver(this.patch, this.fom);
        this.dyn.rho = this.rho;
        this.dyn.assembleMass();
        this.currentTime = 0;
        this.trace = [];
        if (this.viz.cpPoints) {
            this.viz.scene.remove(this.viz.cpPoints, this.viz.cpLattice);
            this.viz.cpPoints = null;
            this.viz.cpLattice = null;
        }
        this.viz.updateMesh(this.dyn.u);
        this.mainChart.data.labels = [];
        this.mainChart.data.datasets[0].data = [];
        this.mainChart.update();
    }

    calculateFrequency() {
        if (this.trace.length < 100) return 0;
        const data = this.trace.slice(-500); 
        let peaks = [];
        for (let i = 1; i < data.length - 1; i++) {
            if (data[i].disp > data[i-1].disp && data[i].disp > data[i+1].disp) {
                peaks.push(data[i].t);
            }
        }
        if (peaks.length < 2) return 0;
        const avgPeriod = (peaks[peaks.length - 1] - peaks[0]) / (peaks.length - 1);
        return 1.0 / avgPeriod;
    }

}

window.app = new TransientLab();
