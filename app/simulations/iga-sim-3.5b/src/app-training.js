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
    const nSamplesInput = document.getElementById('input-samples');
    const maxForceInput = document.getElementById('input-max-force');
    
    if (!nSamplesInput || !maxForceInput) {
        console.error("3.5 UI elements missing");
        return;
    }

    const nSamples = parseInt(nSamplesInput.value);
    const maxForce = parseFloat(maxForceInput.value);
    
    btn.disabled = true;
    status.classList.remove('hidden');
    this.romEngine.clearSnapshots();

    // System size logic
    const nTotalDOFs = this.patch.controlPoints.length * this.patch.controlPoints[0].length * 2;
    const nConstrained = this._getConstrainedDofs().length;
    const nActive = nTotalDOFs - nConstrained;

    // Scale defaults based on active system size
    this.deimM = Math.min(Math.max(10, Math.floor(nActive * 0.15)), 80); 
    this.k = Math.min(Math.max(5, Math.floor(this.deimM * 0.7)), 50); 

    // UI Updates
    const uiM = document.getElementById('input-m');
    const uiK = document.getElementById('input-k');
    if (uiM && uiK) {
        uiM.max = nActive;
        uiM.value = this.deimM;
        document.getElementById('m-val').textContent = this.deimM;
        uiK.max = this.deimM;
        uiK.value = this.k;
        document.getElementById('k-val').textContent = this.k;
    }

    const bcs = this.getBCs();
    const forceSnaps = [];
    const snapDisp = [];

    // --- EQUIDISTANT SAMPLING ---
    const forceVals = [];
    if (nSamples === 1) {
        forceVals.push(maxForce);
    } else {
        for (let i = 0; i < nSamples; i++) {
            const f = (maxForce * 0.1) + (maxForce * 0.9 * i) / (nSamples - 1);
            forceVals.push(f);
        }
    }

    // Material
    this.solverFOM.E = 100000;
    this.solverFOM.nu = 0.3;

    let snapIdx = 1;
    for (const f of forceVals) {
        status.textContent = `Snapshot ${snapIdx}/${nSamples} (F=${f.toFixed(0)})...`;
        
        const res = this.solverFOM.solveNonlinear(this.patch, bcs, this.getLoads(f), {steps:4, iterations:15});
        this.romEngine.addSnapshot(res.u);
        snapDisp.push(new Float64Array(res.u));

        const Fint_eq = this.solverFOM.calculateInternalForce(this.patch, res.u);
        const constrained = this._getConstrainedDofs();
        constrained.forEach(dof => Fint_eq[dof] = 0);
        forceSnaps.push(Fint_eq);

        [0.9, 1.1].forEach(scale => {
            const u_pert = res.u.map(v => v * scale);
            const Fint_pert = this.solverFOM.calculateInternalForce(this.patch, u_pert);
            constrained.forEach(dof => Fint_pert[dof] = 0);
            forceSnaps.push(Fint_pert);
        });
        
        snapIdx++;
        await new Promise(r => setTimeout(r, 5));
    }

    status.textContent = 'Computing POD basis...';
    await new Promise(r => setTimeout(r, 10));
    const podInfo = this.romEngine.computePOD(this.k);
    
    const energyEl = document.getElementById('energy-val');
    if (energyEl) energyEl.textContent = `Energy: ${(podInfo.energy * 100).toFixed(4)}%`;
    
    document.getElementById('input-k').disabled = false;
    document.getElementById('input-k').max = Math.min(forceSnaps.length, 100);

    // Force Enrichment
    const Phi = this.romEngine.Phi;
    const bcs_dofs = this._getConstrainedDofs();
    for (let j = 0; j < Math.min(Phi.columns, 20); j++) {
        const u_mode = new Float64Array(Phi.rows);
        for (let d = 0; d < Phi.rows; d++) u_mode[d] = Phi.get(d, j);
        [10.0, -10.0].forEach(scale => {
            const u_pert = u_mode.map(v => v * scale);
            const f_pert = this.solverFOM.calculateInternalForce(this.patch, u_pert);
            bcs_dofs.forEach(d => f_pert[d] = 0);
            forceSnaps.push(f_pert);
        });
    }

    status.textContent = `Training U-DEIM (m=${this.deimM}, kf=${this.k})...`;
    await new Promise(r => setTimeout(r, 10));
    const deimInfo = await this.deimEngine.train(this.solverFOM, this.romEngine, this.patch, snapDisp, this.deimM, this.k, bcs);
    document.getElementById('input-m').disabled = false;
    document.getElementById('deim-info').textContent = `${deimInfo.m} interpolation points | ${deimInfo.elementCount} active elements`;

    this.isTrained = true;
    this.forceSnaps = forceSnaps;
    this.snapDisp = snapDisp;
    btn.disabled = false;
    document.getElementById('btn-compare').disabled = false;
    document.getElementById('btn-explorer').disabled = false;
    status.textContent = `Training complete ✓`;
    
    this.updatePhysics();
    this.runOnlineAudit();
};
