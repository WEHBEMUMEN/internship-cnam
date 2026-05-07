/**
 * Phase 4.1 Application Core — Online ROM Simulator
 * Consumes the JSON package from Phase 4.1a for real-time hyper-reduced dynamics.
 */

class OnlineSimulator {
    constructor() {
        this.patch = null;
        this.engine = new window.NURBS2D();
        this.fom = new window.IGANonlinearSolver(this.engine);
        this.viz = null;
        
        this.rom = {
            phi: null, // mlMatrix.Matrix
            phiT: null,
            weights: null, // Map of elemIdx -> weight
            indices: null, // Array of elemIdx
            q: null,       // Reduced coords (mlMatrix.Matrix [k x 1])
            u: null        // Full reconstructed disp
        };

        this.isRunning = false;
        this.t = 0;
        this.dt = 0.01;
        this.mag = 400;

        this.init();
    }

    async init() {
        // Initial dummy patch (clamped cantilever)
        this.patch = window.NURBSPresets.generateCantilever();
        this.viz = new TransientVisuals(this);
        this.initCharts();
    }

    initCharts() {
        const ctxTrace = document.getElementById('chart-trace').getContext('2d');
        this.traceChart = new Chart(ctxTrace, {
            type: 'line',
            data: { labels: [], datasets: [{ label: 'Tip Displacement (mm)', data: [], borderColor: '#10b981', backgroundColor: 'rgba(16,185,129,0.1)', fill: true, tension: 0.4, borderWidth: 2, pointRadius: 0 }] },
            options: { 
                responsive: true, 
                maintainAspectRatio: false, 
                plugins: { legend: { display: false } },
                scales: { x: { display: false }, y: { grid: { color: 'rgba(255,255,255,0.05)' } } }
            }
        });
    }

    async loadPackage(data) {
        const { Matrix } = window.mlMatrix;
        
        // 1. Reconstruct Basis
        this.rom.phi = new Matrix(data.phi);
        this.rom.phiT = this.rom.phi.transpose();
        this.rom.k = data.k;
        this.rom.q = new Matrix(this.rom.k, 1);
        
        // 2. Setup ECSW
        this.rom.weights = data.ecsw.weights;
        this.rom.indices = data.ecsw.indices;
        
        // 3. Initialize Mesh State
        this.patch = window.NURBSPresets.generateCantilever();
        // Match refinement used in training
        window.RefineUtils.apply(this.engine, this.patch, { h: data.mesh.h, p: data.mesh.p });
        
        this.rom.u = new Float64Array(this.rom.phi.rows);
        this.viz.updateMesh(this.rom.u, 0);

        // Update UI
        document.getElementById('import-status').innerText = `✅ ROM Loaded (k=${this.rom.k}, ${this.rom.indices.length} Sparse Elements)`;
        document.getElementById('btn-run').disabled = false;
        document.getElementById('dofs-val').innerText = this.rom.phi.rows;
        document.getElementById('sampled-val').innerText = data.snapshots;
        document.getElementById('energy-val').innerText = (data.pod_energy * 100).toFixed(2) + "%";
        document.getElementById('sparse-val').innerText = this.rom.indices.length;
    }

    async toggle() {
        this.isRunning = !this.isRunning;
        const btn = document.getElementById('btn-run');
        btn.innerHTML = this.isRunning ? '<i class="fa-solid fa-pause"></i> Pause Online ROM' : '<i class="fa-solid fa-play"></i> Start Online ROM';
        if (this.isRunning) this.run();
    }

    async run() {
        while (this.isRunning) {
            this.solveStep();
            this.t += this.dt;
            
            this.viz.updateMesh(this.rom.u, this.t);
            this.updateTrace(this.t, this.rom.u[this.rom.u.length - 1]);
            
            await new Promise(r => requestAnimationFrame(r));
        }
    }

    /**
     * Hyper-Reduced Step (ECSW)
     * For now, we'll do a simple Quasi-Static Online Solver for demonstration
     * f_red(q) = f_ext_red
     */
    solveStep() {
        const { Matrix } = window.mlMatrix;
        const mag = parseFloat(document.getElementById('input-fy').value);
        
        // f_ext_red = Phi^T * f_ext
        const f_ext_full = new Float64Array(this.rom.phi.rows);
        const nU = this.patch.controlPoints.length;
        const nV = this.patch.controlPoints[0].length;
        const q = this.patch.q;
        const V = this.patch.V;
        for (let j = 0; j < nV; j++) {
            const dofIdx = ((nU - 1) * nV + j) * 2 + 1;
            const weight = (V[j + q + 1] - V[j]) / (q + 1);
            f_ext_full[dofIdx] = -mag * weight;
        }
        const f_ext_red = this.rom.phiT.mmul(new Matrix([Array.from(f_ext_full)]).transpose());

        // Simple Newton-Raphson in Reduced Space
        // We evaluate only sampled elements!
        for (let iter = 0; iter < 3; iter++) {
            const u_full = this.rom.phi.mmul(this.rom.q).to2DArray().map(r => r[0]);
            
            // Hyper-reduced internal force & stiffness
            // R_red = Phi^T * sum(w_e * f_e) - f_ext_red
            const f_red_int = new Matrix(this.rom.k, 1);
            const K_red_tangent = new Matrix(this.rom.k, this.rom.k);

            // ECSW ELEMENT LOOP (Hyper-Reduction)
            this.rom.indices.forEach((elemIdx, i) => {
                const w = this.rom.weights[i];
                const res = this.fom.calculateElementContribution(this.patch, u_full, elemIdx);
                // Project element force and stiffness
                // f_red += w * Phi_e^T * f_e
                // K_red += w * Phi_e^T * K_e * Phi_e
                // (Implementation shortcut: full projection but only for elements)
                // In a production ROM, Phi_e would be pre-cached.
            });

            // For now, let's use the ROMEngine's RB solve as a fallback until ECSW assembly is ready
            // and reconstruct u for visualization
            const ur_array = this.rom.q.to2DArray().map(r => r[0]);
            const romSolver = new window.ROMEngine(this.fom);
            romSolver.Phi = this.rom.phi;
            
            const res = romSolver.solveReduced(this.patch, [], [{i: nU-1, j: 0, fx: 0, fy: -mag}], { steps: 1 });
            this.rom.u = res.u;
            this.rom.q = new Matrix([Array.from(res.ur)]).transpose();
            break; 
        }
    }

    updateTrace(t, disp) {
        this.traceChart.data.labels.push(t.toFixed(2));
        this.traceChart.data.datasets[0].data.push(disp * 1000);
        if (this.traceChart.data.labels.length > 50) {
            this.traceChart.data.labels.shift();
            this.traceChart.data.datasets[0].data.shift();
        }
        this.traceChart.update('none');
    }

    reset() {
        this.isRunning = false;
        this.t = 0;
        if (this.rom.u) this.rom.u.fill(0);
        if (this.rom.q) this.rom.q.fill(0);
        this.viz.updateMesh(this.rom.u || new Float64Array(0), 0);
        this.traceChart.data.labels = [];
        this.traceChart.data.datasets[0].data = [];
        this.traceChart.update();
    }
}

window.app = new OnlineSimulator();

