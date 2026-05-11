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
        this.snapshotMetadata = []; // Tracks { t, mu } for each snapshot
        this.isTraining = false;
        this.audit = new window.OfflineAudit(this);
        
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
        const uniqueU = [...new Set(this.patch.U)].length - 1;
        const uniqueV = [...new Set(this.patch.V)].length - 1;
        this.elementCount = uniqueU * uniqueV;

        document.getElementById('dofs-val').innerText = this.dyn.dofs;
        document.getElementById('elems-val').innerText = this.elementCount;
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

        const ctxSampling = document.getElementById('chart-sampling').getContext('2d');
        this.samplingChart = new Chart(ctxSampling, {
            type: 'scatter',
            data: { 
                datasets: [
                    { label: 'Path Snapshots', data: [], backgroundColor: 'rgba(244,63,94,0.3)', pointRadius: 3 },
                    { label: 'Convergence Point', data: [], backgroundColor: '#10b981', pointRadius: 6, pointStyle: 'rectRot' }
                ] 
            },
            options: { 
                responsive: true, 
                maintainAspectRatio: false, 
                plugins: { legend: { display: true, labels: { color: '#94a3b8', boxWidth: 10, font: { size: 9 } } } },
                scales: { 
                    x: { title: { display: true, text: 'Time (t)', color: '#94a3b8', font: { size: 10 } }, grid: { color: 'rgba(255,255,255,0.05)' } },
                    y: { title: { display: true, text: 'Load (\u03bc)', color: '#94a3b8', font: { size: 10 } }, grid: { color: 'rgba(255,255,255,0.05)' } }
                }
            }
        });
    }

    async runTraining() {
        if (this.isTraining) return;
        this.isTraining = true;
        this.audit.start();
        this.audit.log(1, "Starting parametric Load Sweep (High-Fidelity FOM)...");
        
        const progContainer = document.getElementById('progress-container');
        const progStatus = document.getElementById('progress-status');
        const progFill = document.getElementById('progress-fill');
        const progPercent = document.getElementById('progress-percent');
        
        if (progContainer) progContainer.style.display = 'block';

        const fMin = parseFloat(document.getElementById('input-fmin').value);
        const fMax = parseFloat(document.getElementById('input-fmax').value);
        const samples = parseInt(document.getElementById('input-samples').value);

        this.snapshots = [];
        this.snapshotMetadata = [];
        this.samplingChart.data.datasets[0].data = [];
        this.samplingChart.data.datasets[1].data = [];
        this.samplingChart.update();

        this.dyn.u.fill(0);
        this.dyn.v.fill(0);
        this.dyn.a.fill(0);

        const convTol = 1e-7;

        for (let s = 0; s < samples; s++) {
            const mag = fMin + (fMax - fMin) * (s / (samples - 1));
            const Fext = new Float64Array(this.dyn.dofs);

            const nU = this.patch.controlPoints.length;
            const nV = this.patch.controlPoints[0].length;
            const q = this.patch.q;
            const V = this.patch.V;

            for (let j = 0; j < nV; j++) {
                const dofIdx = ((nU - 1) * nV + j) * 2 + 1;
                const weight = (V[j + q + 1] - V[j]) / (q + 1);
                Fext[dofIdx] = -mag * weight;
            }

            // Adaptive Sampling: Run until steady state
            let converged = false;
            let uPrev = new Float64Array(this.dyn.u);
            
            for(let step=1; step<=150; step++) { // Allow more room for high loads
                const dt = 0.1;
                
                // Solve with iteration callback for real-time visualization
                await this.dyn.solveStep(dt, Fext, {
                    onIter: async (u_next, iter, norm) => {
                        this.viz.updateMesh(u_next, s);
                        if (progStatus) progStatus.innerText = `Solving ${mag.toFixed(0)}N (Iter ${iter})...`;
                        await new Promise(r => requestAnimationFrame(r));
                    }
                });
                
                // Convergence check (Steady state)
                let diffNorm = 0;
                for(let i=0; i<this.dyn.u.length; i++) {
                    const d = this.dyn.u[i] - uPrev[i];
                    diffNorm += d*d;
                    uPrev[i] = this.dyn.u[i];
                }
                const err = Math.sqrt(diffNorm);
                const time = step * dt;

                // 1. Path snapshots (every 10 steps)
                if (step % 10 === 0 && !converged) {
                    this.snapshots.push(new Float64Array(this.dyn.u));
                    this.snapshotMetadata.push({ t: time, mu: mag, type: 'path' });
                    this.updateSamplingChart(time, mag, false);
                }

                // 2. Convergence point detection
                if (err < convTol && step > 5) {
                    converged = true;
                    this.snapshots.push(new Float64Array(this.dyn.u));
                    this.snapshotMetadata.push({ t: time, mu: mag, type: 'converged' });
                    this.updateSamplingChart(time, mag, true);
                    this.audit.log(1, `Load Level ${mag.toFixed(0)}N Converged at t=${time.toFixed(1)}s.`);
                    break; 
                }

                if (step === 150) {
                    this.audit.log(1, `Warning: Load Level ${mag.toFixed(0)}N failed to converge in 150 steps.`);
                }
            }
            
            this.viz.updateMesh(this.dyn.u, s);
            this.updateTrace(s, this.dyn.u[this.dyn.dofs - 1]);
            
            // Update Global Progress Bar
            const totalProgress = ((s + 1) / samples) * 100;
            if (progFill) progFill.style.width = totalProgress + '%';
            if (progPercent) progPercent.innerText = Math.round(totalProgress) + '%';
        }
        
        if (progContainer) progContainer.style.display = 'none';
        this.audit.log(1, `Phase 1 Complete. Harvested ${this.snapshots.length} snapshots with adaptive convergence tracking.`);

        // Phase 2: POD
        this.audit.log(2, "Extracting POD Basis via SVD...");
        const k = parseInt(document.getElementById('input-k').value);
        const { Phi, energyTrace } = this.trainer.computePOD(this.snapshots, k);
        this.updateEnergyChart(energyTrace);
        this.audit.reportSVD(energyTrace, Phi.columns);

        // Verification: Compare Basis against ALL snapshots (Audit requirement)
        let totalErr = 0, maxErr = 0;
        const { Matrix } = window.mlMatrix;
        const PhiT = Phi.transpose();

        this.snapshots.forEach(u_fom => {
            const u_mat = new Matrix([Array.from(u_fom)]).transpose();
            // Project & Reconstruct: u_rom = Phi * Phi^T * u_fom
            const q = PhiT.mmul(u_mat);
            const u_rec_mat = Phi.mmul(q);
            const u_rec = u_rec_mat.to2DArray().map(r => r[0]);

            let diffNorm = 0, fomNorm = 0;
            for(let i=0; i<u_fom.length; i++) {
                const d = u_fom[i] - u_rec[i];
                diffNorm += d*d;
                fomNorm += u_fom[i] * u_fom[i];
            }
            const relErr = Math.sqrt(diffNorm) / (Math.sqrt(fomNorm) || 1e-12);
            totalErr += relErr;
            maxErr = Math.max(maxErr, relErr);
        });
        this.audit.reportReconstruction(totalErr / this.snapshots.length, maxErr);

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

    async recomputeECSW(mScale) {
        if (!this.snapshots.length || !this.trainer.package.phi) return;
        
        this.audit.log(3, `Retraining ECSW with tol=1e-${mScale}...`);
        const ecswTol = Math.pow(10, -mScale);
        const Phi = new window.mlMatrix.Matrix(this.trainer.package.phi);
        
        const ecswResult = await this.trainer.trainECSW(this.snapshots, Phi, ecswTol);
        
        const totalElements = ([...new Set(this.patch.U)].length - 1) * ([...new Set(this.patch.V)].length - 1);
        this.audit.reportECSW(totalElements, ecswResult.indices.length, ecswTol);
        
        // Re-project reduced matrices
        this.trainer.computeReducedMatrices(Phi);
        this._updateStats(Phi.columns, this.trainer.package.pod_energy);
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

    updateSamplingChart(t, mu, isConverged) {
        const datasetIdx = isConverged ? 1 : 0;
        this.samplingChart.data.datasets[datasetIdx].data.push({ x: t, y: mu });
        this.samplingChart.update('none');
    }

    _updateStats(k, energy) {
        document.getElementById('dofs-val').innerText = this.dyn.dofs;
        document.getElementById('elems-val').innerText = this.elementCount || "—";
        document.getElementById('sampled-val').innerText = this.snapshots.length;
        document.getElementById('energy-val').innerText = (energy * 100).toFixed(4) + "%";
        document.getElementById('sparse-val').innerText = "Pending";
    }
}

window.app = new OfflineLab();
