/**
 * Phase 3.4b: U-DEIM Benchmark — Training Logic
 */

UDEIMBenchmarkApp.prototype.trainAll = async function() {
    const btn = document.getElementById('btn-train');
    const status = document.getElementById('train-status');
    btn.disabled = true;
    status.classList.remove('hidden');
    this.romEngine.clearSnapshots();

    // --- DYNAMIC SCALING LOGIC (DISABLED TO ALLOW MANUAL TUNING) ---
    const nTotalDOFs = this.patch.controlPoints.length * this.patch.controlPoints[0].length * 2;
    const nConstrained = this._getConstrainedDofs().length;
    const nActive = nTotalDOFs - nConstrained;

    const uiM = document.getElementById('input-m');
    const uiK = document.getElementById('input-k');
    if (uiM && uiK) {
        uiM.max = nActive;
        // Respect current UI values instead of overriding them
        this.deimM = parseInt(uiM.value);
        this.k = parseInt(uiK.value);
        uiK.max = this.deimM; 
    }

    // ----------------------------------------------------------------


    const strategy = document.getElementById('input-sampling-strategy').value;
    const bcs = this.getBCs();
    const snapDisp = [];

    const E_vals = [50000, 100000, 150000];
    const nu_vals = [0.3, 0.4];
    const load_fracs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];

    const allCandidates = [];
    for (const E of E_vals) for (const nu of nu_vals) for (const frac of load_fracs) {
        allCandidates.push({ E, nu, f: this.loadMag * frac });
    }

    this.testSet8020 = null;
    this.samplingStrategy = strategy;

    if (strategy === '8020') {
        for (let i = allCandidates.length - 1; i > 0; i--) {
            const j = Math.floor(Math.random() * (i + 1));
            [allCandidates[i], allCandidates[j]] = [allCandidates[j], allCandidates[i]];
        }
        const trainCount = Math.floor(allCandidates.length * 0.8);
        this.testSet8020 = allCandidates.slice(trainCount);
        const trainSet = allCandidates.slice(0, trainCount);

        let snapIdx = 1;
        for (const state of trainSet) {
            this.solverFOM.E = state.E;
            this.solverFOM.nu = state.nu;
            status.textContent = `FOM snapshot ${snapIdx}/${trainSet.length} (E=${state.E/1000}k, nu=${state.nu}, F=${state.f.toFixed(0)})...`;
            
            const res = this.solverFOM.solveNonlinear(this.patch, bcs, this.getLoads(state.f), {
                steps: 3, 
                iterations: 12,
                onProgress: (info) => {
                    if (info.u) {
                        this.romEngine.addSnapshot(info.u);
                        snapDisp.push(info.u);
                    }
                }
            });
            this.romEngine.addSnapshot(res.u);
            snapDisp.push(new Float64Array(res.u));
            snapIdx++;
            await new Promise(r => setTimeout(r, 5));
        }
    } else {
        const trainSet = [allCandidates[0], allCandidates[allCandidates.length-1]];
        const candidatePool = allCandidates.slice(1, allCandidates.length-1);
        const targetSnapshots = 25;


        for(let i=0; i<trainSet.length; i++) {
            const state = trainSet[i];
            this.solverFOM.E = state.E; this.solverFOM.nu = state.nu;
            status.textContent = `Greedy Init ${i+1}/2...`;
            const res = this.solverFOM.solveNonlinear(this.patch, bcs, this.getLoads(state.f), {
                steps: 3, 
                iterations: 12,
                onProgress: (info) => {
                    if (info.u) {
                        this.romEngine.addSnapshot(info.u);
                        snapDisp.push(info.u);
                    }
                }
            });
            this.romEngine.addSnapshot(res.u);
            snapDisp.push(new Float64Array(res.u));
            await new Promise(r => setTimeout(r, 5));
        }

        for(let step = trainSet.length; step < targetSnapshots; step++) {
            status.textContent = `Greedy Selection Step ${step+1}/${targetSnapshots}...`;
            await new Promise(r => setTimeout(r, 10));

            this.romEngine.computePOD(this.k, this.patch, bcs);

            const tempDeim = new UDEIMEngine();
            await tempDeim.train(this.solverFOM, this.romEngine, this.patch, snapDisp, this.deimM, this.k, bcs);


            let maxError = -1, worstIdx = -1;
            for (let i = 0; i < candidatePool.length; i++) {
                const cand = candidatePool[i];
                this.solverFOM.E = cand.E; this.solverFOM.nu = cand.nu;
                try {
                    const loads = this.getLoads(cand.f);
                    const res = tempDeim.solveReduced(this.solverFOM, this.romEngine, this.patch, bcs, loads, {iterations:10, steps:1});
                    const Fint_exact = this.solverFOM.calculateInternalForce(this.patch, res.u);
                    this.solverFOM.applyPenaltyConstraints(null, Fint_exact, res.u, this.patch, bcs);

                    const Fext = new Float64Array(Fint_exact.length);
                    loads.forEach(l => {
                        const idx = (l.i * this.patch.controlPoints[0].length + l.j) * 2;
                        Fext[idx] += l.fx; Fext[idx+1] += l.fy;
                    });
                    let err2 = 0;
                    for(let j=0; j<Fext.length; j++) err2 += Math.pow(Fext[j] - Fint_exact[j], 2);
                    const err = Math.sqrt(err2);
                    if (err > maxError) { maxError = err; worstIdx = i; }
                } catch(e) { maxError = Infinity; worstIdx = i; break; }
            }

            if (worstIdx !== -1) {
                const worst = candidatePool[worstIdx];
                candidatePool.splice(worstIdx, 1);
                this.solverFOM.E = worst.E; this.solverFOM.nu = worst.nu;
                const res = this.solverFOM.solveNonlinear(this.patch, bcs, this.getLoads(worst.f), {steps:3, iterations:12});
                this.romEngine.addSnapshot(res.u);
                snapDisp.push(new Float64Array(res.u));
            }
        }
    }
    
    this.solverFOM.E = 100000; this.solverFOM.nu = 0.3;
    status.textContent = 'Computing POD basis...';
    await new Promise(r => setTimeout(r, 10));
    const podInfo = this.romEngine.computePOD(this.k, this.patch, bcs);

    document.getElementById('energy-val').textContent = `Energy: ${(podInfo.energy*100).toFixed(2)}%`;
    document.getElementById('input-k').disabled = false;

    status.textContent = `Training U-DEIM (m=${this.deimM}, kf=${this.k})...`;
    await new Promise(r => setTimeout(r, 10));
    const deimInfo = await this.deimEngine.train(this.solverFOM, this.romEngine, this.patch, snapDisp, this.deimM, this.k, bcs);

    document.getElementById('input-m').disabled = false;
    document.getElementById('deim-info').textContent = `${deimInfo.m} interpolation points selected`;

    this.isTrained = true;
    this.snapDisp = snapDisp;
    btn.disabled = false;
    document.getElementById('btn-compare').disabled = false;
    document.getElementById('btn-explorer').disabled = false;
    status.textContent = `Training complete ✓`;
    setTimeout(() => status.classList.add('hidden'), 5000);
    this.updatePhysics();
    this.deimEngine.audit(this.solverFOM, this.romEngine, this.patch, this.snapDisp);
};

UDEIMBenchmarkApp.prototype.runComparison = async function() {
    if (!this.isTrained) return;
    document.getElementById('btn-compare').disabled = true;
    const data = {};
    for (const m of Object.keys(METHODS)) {
        try {
            const { result, meta } = this.solve(m, this.loadMag);
            data[m] = meta;
        } catch(e) { data[m] = { error: e.message }; }
        await new Promise(r => setTimeout(r, 5));
    }
    this.updateSpeedupChart(data);
    document.getElementById('btn-compare').disabled = false;
};

UDEIMBenchmarkApp.prototype.runFDCurves = async function() {
    if (!this.isTrained) return;
    const loads = [20, 50, 100, 150, 200, 300, 400, 500];
    const datasets_fd = {};
    const datasets_err = {};
    const fom_results = [];
    for (const f of loads) {
        try { const { result } = this.solve('fom', f); fom_results.push(result.u); }
        catch(e) { fom_results.push(null); }
    }
    for (const m of Object.keys(METHODS)) {
        datasets_fd[m] = []; datasets_err[m] = [];
        for (let i = 0; i < loads.length; i++) {
            const f = loads[i];
            try { 
                const { result, meta } = this.solve(m, f); 
                datasets_fd[m].push(Math.abs(meta.tipDisp || 0));
                if (m === 'fom') datasets_err[m].push(1e-12);
                else if (fom_results[i]) {
                    let err2 = 0, norm2 = 0;
                    for(let j=0; j<result.u.length; j++) {
                        err2 += Math.pow(result.u[j] - fom_results[i][j], 2);
                        norm2 += Math.pow(fom_results[i][j], 2);
                    }
                    datasets_err[m].push((Math.sqrt(err2) / Math.max(Math.sqrt(norm2), 1e-12)) * 100);
                } else datasets_err[m].push(null);
            } catch(e) { datasets_fd[m].push(0); datasets_err[m].push(null); }
        }
    }
    this.fdChart.data.labels = loads;
    this.fdChart.data.datasets = Object.keys(METHODS).map(m => ({
        label: METHODS[m].label, data: datasets_fd[m],
        borderColor: METHODS[m].color, backgroundColor: METHODS[m].color + '33',
        borderWidth: 2, pointRadius: 3, tension: 0.3
    }));
    this.fdChart.update();
    this.errorChart.data.labels = loads;
    this.errorChart.data.datasets = Object.keys(METHODS).filter(m => m !== 'fom').map(m => ({
        label: METHODS[m].label, data: datasets_err[m],
        borderColor: METHODS[m].color, backgroundColor: METHODS[m].color + '33',
        borderWidth: 2, pointRadius: 3, tension: 0.3
    }));
    this.errorChart.update();
};
