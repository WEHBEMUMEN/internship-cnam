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
        
        this.init();
    }

    async init() {
        // Setup NURBS Patch (Cantilever by default)
        const preset = window.nurbsPresets.cantilever();
        this.patch = new window.NURBS2D(preset);
        
        // Setup Nonlinear Solver (FOM)
        this.fom = new window.IGANonlinearSolver(this.patch);
        
        // Setup Dynamics Solver
        this.dyn = new DynamicsSolver(this.patch, this.fom);
        
        // Assemble Initial Mass
        this.dyn.assembleMass();
        
        // Setup Visuals
        this.viz = new TransientVisuals(this);
        this.viz.updateMesh(this.dyn.u);
        
        this.initCharts();
        this.render();
    }

    initCharts() {
        const ctx = document.getElementById('chart-full-trace').getContext('2d');
        this.mainChart = new Chart(ctx, {
            type: 'line',
            data: { labels: [], datasets: [{ label: 'Tip Displacement (mm)', data: [], borderColor: '#0ea5e9', backgroundColor: 'rgba(14,165,233,0.1)', fill: true, tension: 0.4, borderWidth: 2, pointRadius: 0 }] },
            options: { responsive: true, maintainAspectRatio: false, scales: { x: { display: true, title: { display: true, text: 'Time (s)', color: '#94a3b8' }, grid: { color: 'rgba(148,163,184,0.1)' } }, y: { display: true, title: { display: true, text: 'Disp (mm)', color: '#94a3b8' }, grid: { color: 'rgba(148,163,184,0.1)' } } }, plugins: { legend: { display: false } } }
        });
    }

    async runSimulation() {
        if (this.isRunning) return;
        this.isRunning = true;
        
        // Reset state
        this.currentTime = 0;
        this.dyn.u.fill(0);
        this.dyn.v.fill(0);
        this.dyn.a.fill(0);
        this.trace = [];
        this.mainChart.data.labels = [];
        this.mainChart.data.datasets[0].data = [];
        
        const T = parseFloat(document.getElementById('input-time').value);
        const dt = parseFloat(document.getElementById('input-dt').value);
        const F0 = parseFloat(document.getElementById('input-f0').value);
        const forceType = document.getElementById('input-force-type').value;

        // Update Rayleigh params
        this.dyn.alpha = parseFloat(document.getElementById('input-alpha').value);
        this.dyn.betaR = parseFloat(document.getElementById('input-beta').value);

        while (this.currentTime < T && this.isRunning) {
            // Compute External Force at current time step
            const Fext = new Float64Array(this.dyn.dofs);
            let mag = 0;
            if (forceType === 'step') mag = F0;
            else if (forceType === 'impulse') mag = (this.currentTime < 0.05) ? F0 : 0;
            else if (forceType === 'harmonic') mag = F0 * Math.sin(2 * Math.PI * 10 * this.currentTime);

            // Apply to tip (last control point, Y DOF)
            Fext[this.dyn.dofs - 1] = -mag; 

            // Solve Step
            const result = await this.dyn.solveStep(dt, Fext);
            
            this.currentTime += dt;
            
            // Record Trace
            const tipDisp = this.dyn.u[this.dyn.dofs - 1];
            this.trace.push({ t: this.currentTime, disp: tipDisp });
            
            // UI Updates
            if (Math.round(this.currentTime / dt) % 2 === 0) {
                this.updateUI(result, tipDisp);
                this.viz.updateMesh(this.dyn.u);
            }
        }
        
        this.isRunning = false;
        document.getElementById('run-label').innerText = 'Complete';
        document.getElementById('run-dot').style.background = '#94a3b8';
    }

    updateUI(result, tipDisp) {
        document.getElementById('curr-time').innerText = this.currentTime.toFixed(3) + 's';
        document.getElementById('iters-val').innerText = result.iters;
        document.getElementById('tip-disp').innerText = (tipDisp * 1000).toFixed(2) + 'mm';
        
        // Update Chart (limit to last 100 pts for performance)
        this.mainChart.data.labels.push(this.currentTime.toFixed(3));
        this.mainChart.data.datasets[0].data.push(tipDisp * 1000);
        if (this.mainChart.data.labels.length > 200) {
            this.mainChart.data.labels.shift();
            this.mainChart.data.datasets[0].data.shift();
        }
        this.mainChart.update('none');
    }

}

window.app = new TransientLab();
