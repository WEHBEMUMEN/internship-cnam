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
    meta.tipDisp = result.u[tipIdx];

    if (this.lastFomResult && method !== 'fom') {
        let num = 0, den = 0;
        for (let i = 0; i < result.u.length; i++) {
            num += (result.u[i] - this.lastFomResult.u[i])**2;
            den += this.lastFomResult.u[i]**2;
        }
        meta.error = den > 0 ? Math.sqrt(num/den) : 0;
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
