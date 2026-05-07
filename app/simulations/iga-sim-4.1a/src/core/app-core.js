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
        this.h = 0;
        this.p = 1;
        
        this.snapshots = []; 
        this.isTraining = false;
        
        this.init();
    }

    async updateRefinement(h, p = null) {
        if (this.isTraining) return;
        this.h = h;
        if (p !== null) this.p = p;
        
        const preset = window.NURBSPresets.generateCantilever();
        this.patch = preset;
        window.RefineUtils.apply(this.engine, this.patch, { h: this.h, p: this.p || this.patch.p });
        
        this.fom = new window.IGANonlinearSolver(this.engine);
        this.dyn = new DynamicsSolver(this.patch, this.fom);
        this.dyn.assembleMass();
        
        if (this.viz.cpPoints) {
            this.viz.scene.remove(this.viz.cpPoints, this.viz.cpLattice);
            this.viz.cpPoints = null;
            this.viz.cpLattice = null;
        }
        this.viz.updateMesh(this.dyn.u);
        document.getElementById('dofs-val').innerText = this.dyn.dofs;
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
            data: { labels: [], datasets: [{ label: 'Tip Displacement (mm)', data: [], borderColor: '#f43f5e', backgroundColor: 'rgba(244,63,94,0.1)', fill: true, tension: 0.4, borderWidth: 2, pointRadius: 0 }] },
            options: { 
                responsive: true, 
                maintainAspectRatio: false, 
                plugins: { 
                    legend: { display: false },
                    zoom: {
                        zoom: { wheel: { enabled: true }, pinch: { enabled: true }, mode: 'x' },
                        pan: { enabled: true, mode: 'x' }
                    }
                },
                scales: {
                    x: { grid: { color: 'rgba(255,255,255,0.05)' } },
                    y: { grid: { color: 'rgba(255,255,255,0.05)' } }
                }
            }
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
        this.audit.start();
        this.audit.log(1, "Starting high-fidelity FOM simulation...");
        
        const T = parseFloat(document.getElementById('input-time').value);
        const steps = parseInt(document.getElementById('input-steps').value);
        const dt_calc = T / steps;

        this.snapshots = [];
        this.dyn.u.fill(0);
        this.dyn.v.fill(0);
        this.dyn.a.fill(0);

        for (let s = 0; s < steps; s++) {
            const Fext = new Float64Array(this.dyn.dofs);
            const mag = 400; 

            const nU = this.patch.controlPoints.length;
            const nV = this.patch.controlPoints[0].length;
            const q = this.patch.q;
            const V = this.patch.V;

            for (let j = 0; j < nV; j++) {
                const dofIdx = ((nU - 1) * nV + j) * 2 + 1;
                // EXACT IGA CONSISTENT FORCE DISTRIBUTION:
                // \int N_{j,q} d\eta = (V_{j+q+1} - V_j) / (q + 1)
                const weight = (V[j + q + 1] - V[j]) / (q + 1);
                Fext[dofIdx] = -mag * weight;
            }

            await this.dyn.solveStep(dt_calc, Fext);
            this.snapshots.push(new Float64Array(this.dyn.u));
            
            // UI Performance Optimization: Render every 5 steps or skip more if many steps
            const skipRate = Math.max(1, Math.floor(steps / 20)); 
            if (s % skipRate === 0) {
                this.viz.updateMesh(this.dyn.u, s * dt_calc);
                this.updateTrace(s * dt_calc, this.dyn.u[this.dyn.dofs - 1]);
                this.audit.log(1, `Snapshot collection: ${s}/${steps}...`);
            }
        }
        this.audit.log(1, `Phase 1 Complete. Harvested ${this.snapshots.length} snapshots.`);

        // Phase 2: POD
        this.audit.log(2, "Extracting POD Basis via SVD...");
        const k = parseInt(document.getElementById('input-k').value);
        const { Phi, energyTrace } = this.trainer.computePOD(this.snapshots, k);
        this.updateEnergyChart(energyTrace);
        this.audit.reportSVD(energyTrace, Phi.columns);

        // Phase 3: ECSW
        this.audit.log(3, "Starting ECSW Sparse Training...");
        const mScale = parseInt(document.getElementById('input-m').value);
        const ecswTol = Math.pow(10, -mScale);
        const ecswResult = await this.trainer.trainECSW(this.snapshots, Phi, ecswTol);
        
        const totalElements = ([...new Set(this.patch.U)].length - 1) * ([...new Set(this.patch.V)].length - 1);
        this.audit.reportECSW(totalElements, ecswResult.indices.length, ecswTol);

        // Phase 4: Projection
        this.audit.log(4, "Assembling Reduced Matrices...");
        this.trainer.computeReducedMatrices(Phi);
        
        // Phase 5: Verification
        this.audit.log(5, "Verifying Package Integrity...");
        this.trainer.verifyPackage();
        
        document.getElementById('btn-export').disabled = false;
        document.getElementById('btn-verify').disabled = false;
        
        this.isTraining = false;
        this.audit.success();
        
        // Update Final UI Stats
        const actualK = Phi.columns;
        const finalEnergy = energyTrace[actualK - 1];
        this._updateStats(actualK, finalEnergy);
    }

    updateTrace(t, disp) {
        this.viz.updateMesh(this.dyn.u, t);
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
