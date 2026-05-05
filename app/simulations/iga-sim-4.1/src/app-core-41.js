/**
 * Phase 4.1 Application Core — Offline ROM Trainer
 */

class OfflineLab {
    constructor() {
        this.patch = null;
        this.fom = null;
        this.dyn = null;
        this.viz = null;
        this.trainer = null;
        
        this.snapshots = []; // Displacement snapshots
        this.isTraining = false;
        
        this.init();
    }

    async init() {
        const preset = window.NURBSPresets.generateCantilever();
        this.engine = new window.NURBS2D();
        this.patch = preset;
        this.fom = new window.IGANonlinearSolver(this.engine);
        this.dyn = new DynamicsSolver(this.patch, this.fom);
        this.dyn.assembleMass();
        
        this.viz = new TransientVisuals(this);
        this.trainer = new OfflineTrainer(this);
        
        this.initCharts();
    }

    initCharts() {
        const ctxTrace = document.getElementById('chart-trace').getContext('2d');
        this.traceChart = new Chart(ctxTrace, {
            type: 'line',
            data: { labels: [], datasets: [{ label: 'Tip Displacement', data: [], borderColor: '#f43f5e', fill: false }] },
            options: { responsive: true, maintainAspectRatio: false, plugins: { legend: { display: false } } }
        });

        const ctxEnergy = document.getElementById('chart-energy').getContext('2d');
        this.energyChart = new Chart(ctxEnergy, {
            type: 'bar',
            data: { labels: [], datasets: [{ label: 'Cumulative Energy', data: [], backgroundColor: '#10b981' }] },
            options: { responsive: true, maintainAspectRatio: false, plugins: { legend: { display: false } }, scales: { y: { min: 0.9, max: 1.0 } } }
        });
    }

    async runTraining() {
        if (this.isTraining) return;
        this.isTraining = true;
        
        const T = parseFloat(document.getElementById('input-time').value);
        const dt = parseFloat(document.getElementById('input-dt')?.value || 0.01);
        const steps = parseInt(document.getElementById('input-steps').value);
        const dt_calc = T / steps;

        this.snapshots = [];
        this.dyn.u.fill(0);
        this.dyn.v.fill(0);
        this.dyn.a.fill(0);

        document.getElementById('train-status').classList.remove('hidden');
        document.getElementById('train-status').innerText = "Running Transient FOM...";

        for (let s = 0; s < steps; s++) {
            const Fext = new Float64Array(this.dyn.dofs);
            // Simple step force for training
            Fext[this.dyn.dofs - 1] = -400; 

            await this.dyn.solveStep(dt_calc, Fext);
            
            // Store snapshot
            this.snapshots.push(new Float64Array(this.dyn.u));
            
            if (s % 5 === 0) {
                this.viz.updateMesh(this.dyn.u);
                this.updateTrace(s * dt_calc, this.dyn.u[this.dyn.dofs - 1]);
                document.getElementById('train-status').innerText = `Snapshot ${s}/${steps}...`;
            }
        }

        document.getElementById('train-status').innerText = "Computing POD Basis...";
        const k = parseInt(document.getElementById('input-k').value);
        const { Phi, energyTrace } = this.trainer.computePOD(this.snapshots, k);
        
        this.updateEnergyChart(energyTrace);
        this.trainer.computeReducedMatrices(Phi);
        
        document.getElementById('train-status').innerText = "Training Complete.";
        document.getElementById('btn-export').disabled = false;
        document.getElementById('btn-verify').disabled = false;
        
        this.isTraining = false;
        this._updateStats(Phi.columns, energyTrace[k-1]);
    }

    updateTrace(t, disp) {
        this.traceChart.data.labels.push(t.toFixed(2));
        this.traceChart.data.datasets[0].data.push(disp * 1000);
        this.traceChart.update('none');
    }

    updateEnergyChart(trace) {
        this.energyChart.data.labels = trace.map((_, i) => i + 1);
        this.energyChart.data.datasets[0].data = trace;
        this.energyChart.update();
    }

    _updateStats(k, energy) {
        document.getElementById('dofs-val').innerText = this.dyn.dofs;
        document.getElementById('sampled-val').innerText = this.snapshots.length;
        document.getElementById('energy-val').innerText = (energy * 100).toFixed(4) + "%";
        document.getElementById('sparse-val').innerText = "Pending";
    }
}

window.app = new OfflineLab();
