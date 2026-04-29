/**
 * DEIM Benchmark — Offline Training & Sampling
 */

DEIMBenchmarkApp.prototype.getBCs = function() {
    const nV = this.patch.controlPoints[0].length;
    const bcs = [];
    for (let j = 0; j < nV; j++) bcs.push({i:0, j, axis:'both', value:0});
    return bcs;
};

DEIMBenchmarkApp.prototype._getConstrainedDofs = function() {
    const nV = this.patch.controlPoints[0].length;
    const bcs = this.getBCs();
    const dofs = [];
    bcs.forEach(bc => {
        const base = (bc.i * nV + bc.j) * 2;
        if (bc.axis === 'x' || bc.axis === 'both') dofs.push(base);
        if (bc.axis === 'y' || bc.axis === 'both') dofs.push(base + 1);
    });
    return dofs;
};

DEIMBenchmarkApp.prototype.getLoads = function(mag) {
    const nU = this.patch.controlPoints.length, nV = this.patch.controlPoints[0].length;
    const loads = [];
    if (this.loadType === 'tip') {
        for (let j = 0; j < nV; j++) loads.push({type:'nodal', i:nU-1, j, fx:0, fy:-mag/nV});
    } else {
        for (let i = 0; i < nU; i++) loads.push({type:'nodal', i, j:nV-1, fx:0, fy:-mag/nU});
    }
    return loads;
};

DEIMBenchmarkApp.prototype.trainAll = async function() {
    const btn = document.getElementById('btn-train');
    const status = document.getElementById('train-status');
    btn.disabled = true;
    status.classList.remove('hidden');
    this.romEngine.clearSnapshots();

    const strategy = document.getElementById('input-sampling-strategy').value;
    const bcs = this.getBCs();
    const forceSnaps = [];
    const snapDisp = [];

    const E_vals = [50000, 100000, 150000];
    const nu_vals = [0.25, 0.4];
    const load_fracs = [0.25, 0.50, 0.75, 1.00];
    const maxTrainLoad = Math.max(600, this.loadMag);
    const allCandidates = [];
    for (const E of E_vals) for (const nu of nu_vals) for (const frac of load_fracs) {
        allCandidates.push({ E, nu, f: maxTrainLoad * frac });
    }

    this.testSet8020 = null;
    this.samplingStrategy = strategy;

    if (strategy === '8020') {
        console.group("DEIM Training: 80/20 Train/Test Split");
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
            
            const res = this.solverFOM.solveNonlinear(this.patch, bcs, this.getLoads(state.f), {steps:3, iterations:12});
            this.romEngine.addSnapshot(res.u);
            snapDisp.push(new Float64Array(res.u));

            const Fint = this.solverFOM.calculateInternalForce(this.patch, res.u);
            this._getConstrainedDofs().forEach(dof => Fint[dof] = 0);
            forceSnaps.push(Fint);
            
            snapIdx++;
            await new Promise(r => setTimeout(r, 5));
        }
        console.groupEnd();
    } else {
        // Greedy Sampling Logic (omitted for brevity in this split or moved here)
        // [Simplified Greedy logic from main.js]
        console.group("DEIM Training: Greedy Parametric Sampling");
        const trainSet = [allCandidates[0], allCandidates[allCandidates.length-1]];
        const candidatePool = allCandidates.slice(1, allCandidates.length-1);
        const targetSnapshots = 12;

        for(let i=0; i<trainSet.length; i++) {
            const state = trainSet[i];
            this.solverFOM.E = state.E; this.solverFOM.nu = state.nu;
            status.textContent = `Greedy Init ${i+1}/2...`;
            const res = this.solverFOM.solveNonlinear(this.patch, bcs, this.getLoads(state.f), {steps:3, iterations:12});
            this.romEngine.addSnapshot(res.u);
            snapDisp.push(new Float64Array(res.u));
            const Fint = this.solverFOM.calculateInternalForce(this.patch, res.u);
            forceSnaps.push(Fint);
            await new Promise(r => setTimeout(r, 5));
        }

        for(let step = trainSet.length; step < targetSnapshots; step++) {
            status.textContent = `Greedy Selection Step ${step+1}/${targetSnapshots}...`;
            await new Promise(r => setTimeout(r, 10));

            this.romEngine.computePOD(this.k);
            const tempDeim = new DEIMEngine();
            tempDeim.train(forceSnaps, this.deimM, this.k);
            tempDeim.computeActiveElements(this.patch);
            tempDeim.precomputeReducedTangent(this.solverFOM, this.romEngine, this.patch, snapDisp);
            tempDeim.precomputeReducedPenalty(this.solverFOM, this.romEngine, this.patch, bcs);

            let maxError = -1, worstIdx = -1;
            for (let i = 0; i < candidatePool.length; i++) {
                const cand = candidatePool[i];
                this.solverFOM.E = cand.E; this.solverFOM.nu = cand.nu;
                try {
                    const loads = this.getLoads(cand.f);
                    const res = tempDeim.solveReduced(this.solverFOM, this.romEngine, this.patch, bcs, loads, {iterations:10, steps:1});
                    const Fint_exact = this.solverFOM.calculateInternalForce(this.patch, res.u);
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
                const Fint = this.solverFOM.calculateInternalForce(this.patch, res.u);
                forceSnaps.push(Fint);
            }
        }
        console.groupEnd();
    }
    
    this.solverFOM.E = 100000; this.solverFOM.nu = 0.3;
    status.textContent = 'Computing POD basis...';
    await new Promise(r => setTimeout(r, 10));
    const podInfo = this.romEngine.computePOD(this.k);
    document.getElementById('energy-val').textContent = `Energy: ${(podInfo.energy*100).toFixed(2)}%`;
    document.getElementById('input-k').disabled = false;
    document.getElementById('input-k').max = forceSnaps.length;

    status.textContent = `Training DEIM (m=${this.deimM}, kf=${this.k})...`;
    await new Promise(r => setTimeout(r, 10));
    const constrainedDofs = this._getConstrainedDofs();
    const deimInfo = this.deimEngine.train(forceSnaps, this.deimM, this.k, constrainedDofs);
    document.getElementById('input-m').disabled = false;
    document.getElementById('deim-info').textContent = `${deimInfo.m} interpolation points selected`;

    status.textContent = 'Mapping active elements...';
    this.deimEngine.computeActiveElements(this.patch);
    status.textContent = 'Pre-computing reduced tangent & penalty...';
    await new Promise(r => setTimeout(r, 10));
    this.deimEngine.precomputeReducedTangent(this.solverFOM, this.romEngine, this.patch, snapDisp);
    this.deimEngine.precomputeReducedPenalty(this.solverFOM, this.romEngine, this.patch, bcs);

    this.isTrained = true;
    this.forceSnaps = forceSnaps;
    this.snapDisp = snapDisp;
    btn.disabled = false;
    document.getElementById('btn-compare').disabled = false;
    document.getElementById('btn-explorer').disabled = false;
    status.textContent = `Training complete ✓`;
    setTimeout(() => status.classList.add('hidden'), 5000);

    this.deimEngine.audit(this.solverFOM, this.romEngine, this.patch, this.snapDisp, this.forceSnaps, this._getConstrainedDofs());
    this.updatePhysics();
    this.runOnlineAudit();
};
