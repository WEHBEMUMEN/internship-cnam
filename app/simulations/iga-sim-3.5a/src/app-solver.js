/**
 * DEIM Benchmark — Solver & Comparison logic
 */

DEIMBenchmarkApp.prototype.solve = function(method, mag) {
    const bcs = this.getBCs(), loads = this.getLoads(mag);
    const t0 = performance.now();
    let result, meta = {};

    if (method === 'fom') {
        result = this.solverFOM.solveNonlinear(this.patch, bcs, loads, {iterations:15, steps:3});
    } else if (method === 'galerkin') {
        result = this.romEngine.solveReduced(this.patch, bcs, loads, {iterations:15});
    } else if (method === 'deim') {
        result = this.deimEngine.solveReduced(this.solverFOM, this.romEngine, this.patch, bcs, loads, {iterations:15, steps:1});
        meta.sampled = `${result.sampledCount} / ${result.totalDofs} dofs`;
    }

    const dt = performance.now() - t0;
    meta.time = dt;

    const nU = this.patch.controlPoints.length, nV = this.patch.controlPoints[0].length;
    const tipIdx = ((nU-1)*nV + Math.floor(nV/2)) * 2 + 1;
    
    if (result && result.u) {
        meta.tipDisp = result.u[tipIdx];
        if (this.lastFomResult && method !== 'fom') {
            let num = 0, den = 0;
            for (let i = 0; i < result.u.length; i++) {
                num += (result.u[i] - this.lastFomResult.u[i])**2;
                den += this.lastFomResult.u[i]**2;
            }
            meta.error = den > 0 ? Math.sqrt(num/den) : 0;
        }
    } else {
        meta.tipDisp = null;
        meta.error = null;
    }
    
    return { result, meta };
};

DEIMBenchmarkApp.prototype.updatePhysics = function() {
    if (!this.patch) return;
    if (this.method !== 'fom' && !this.isTrained) { alert('Train first!'); return; }

    if (this.method === 'fom' || !this.lastFomResult) {
        const fom = this.solve('fom', this.loadMag);
        this.lastFomTime = fom.meta.time;
        this.lastFomResult = fom.result;
        if (this.method === 'fom') {
            this.lastResult = fom.result;
            this._updateStats('fom', fom.meta);
            this.sparkline.update(fom.result.residualHistory);
            this.updateMesh(fom.result.u);
            this._render();
            return;
        }
    }

    const { result, meta } = this.solve(this.method, this.loadMag);
    this.lastResult = result;
    this._updateStats(this.method, meta);
    this.sparkline.update(result.residualHistory);
    if (!result.u.some(v => !isFinite(v))) this.updateMesh(result.u);
    this._render();
};

DEIMBenchmarkApp.prototype.runComparison = async function() {
    if (!this.isTrained) return;
    document.getElementById('btn-compare').disabled = true;
    const data = {};
    for (const m of Object.keys(METHODS)) {
        try {
            const { result, meta } = this.solve(m, this.loadMag);
            data[m] = meta;
            data[m].converged = (result.residualHistory && result.residualHistory.length > 0) ? result.residualHistory[result.residualHistory.length-1].norm < 1e-4 : false;
        } catch(e) { 
            console.warn(`${m} failed:`, e); 
            data[m] = { error: e.message };
        }
        await new Promise(r => setTimeout(r, 5));
    }
    this.updateSpeedupChart(data);
    document.getElementById('btn-compare').disabled = false;
};

DEIMBenchmarkApp.prototype.runFDCurves = async function() {
    if (!this.isTrained) return;
    document.getElementById('btn-compare').disabled = true;
    const loads = [];
    for (let f = 10; f <= 600; f += 20) loads.push(f);
    
    const results = { fom: { d: [], e: [] }, galerkin: { d: [], e: [] }, deim: { d: [], e: [] } };
    const originalFom = this.lastFomResult; // Save to prevent UI state pollution
    
    for (const f of loads) {
        // Run FOM
        try {
            const fom = this.solve('fom', f);
            this.lastFomResult = fom.result; // Set temporarily for error calc
            results.fom.d.push(fom.meta.tipDisp);
            results.fom.e.push(0);
        } catch (e) { results.fom.d.push(NaN); results.fom.e.push(NaN); }

        // Run Galerkin
        try {
            const gal = this.solve('galerkin', f);
            results.galerkin.d.push(gal.meta.tipDisp);
            results.galerkin.e.push(gal.meta.error);
        } catch (e) { results.galerkin.d.push(NaN); results.galerkin.e.push(NaN); }

        // Run DEIM
        try {
            const deim = this.solve('deim', f);
            results.deim.d.push(deim.meta.tipDisp);
            results.deim.e.push(deim.meta.error);
        } catch (e) { results.deim.d.push(NaN); results.deim.e.push(NaN); }
        
        await new Promise(r => setTimeout(r, 1)); // Yield UI
    }
    
    this.lastFomResult = originalFom; // Restore original FOM result for the current UI load

    // Update FD Chart
    this.fdChart.data.labels = loads;
    this.fdChart.data.datasets = [
        { label: 'FOM', data: results.fom.d, borderColor: '#64748b', tension: 0.1 },
        { label: 'Galerkin', data: results.galerkin.d, borderColor: '#0ea5e9', borderDash: [5,5], tension: 0.1 },
        { label: 'DEIM', data: results.deim.d, borderColor: '#8b5cf6', tension: 0.1 }
    ];
    this.fdChart.update();

    // Update Error Chart
    this.errorChart.data.labels = loads;
    this.errorChart.data.datasets = [
        { label: 'Galerkin Error', data: results.galerkin.e, borderColor: '#0ea5e9', tension: 0.1 },
        { label: 'DEIM Error', data: results.deim.e, borderColor: '#8b5cf6', tension: 0.1 }
    ];
    this.errorChart.update();

    document.getElementById('btn-compare').disabled = false;
    this.updatePhysics(); // Refresh the UI stats to clear the 64% artifact
};
