/**
 * Phase 3.4b: U-DEIM Benchmark — Solver Logic
 */

UDEIMBenchmarkApp.prototype.loadBenchmark = function() {
    this.clearMesh();
    this.patch = NURBSPresets.generateCantilever(10, 2);
    RefineUtils.apply(this.engine, this.patch, { p: 3, h: this.meshLevel });
    this.solverFOM.E = 100000;
    this.solverFOM.nu = 0.3;
    this.solverFOM.thickness = 1.0;

    this.isTrained = false;
    this.method = 'fom';
    this.lastFomResult = null;
    this.lastFomTime = null;
    document.querySelectorAll('[data-method]').forEach(b => {
        b.classList.remove('active');
        if(b.dataset.method === 'fom') b.classList.add('active');
    });
    document.getElementById('btn-compare').disabled = true;
    document.getElementById('input-k').disabled = true;
    document.getElementById('input-m').disabled = true;
    this.updateMesh(null);
    this._render();
};

UDEIMBenchmarkApp.prototype.clearMesh = function() {
    [this.deformedMesh, this.ghostMesh].forEach(m => { if (m) { this.scene.remove(m); m.geometry.dispose(); } });
    this.deformedMesh = this.ghostMesh = null;
};

UDEIMBenchmarkApp.prototype.getBCs = function() {
    const nV = this.patch.controlPoints[0].length;
    const bcs = [];
    for (let j = 0; j < nV; j++) bcs.push({i:0, j, axis:'both', value:0});
    return bcs;
};

UDEIMBenchmarkApp.prototype._getConstrainedDofs = function() {
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

UDEIMBenchmarkApp.prototype.getLoads = function(mag) {
    const nU = this.patch.controlPoints.length, nV = this.patch.controlPoints[0].length;
    const loads = [];
    if (this.loadType === 'tip') {
        for (let j = 0; j < nV; j++) loads.push({type:'nodal', i:nU-1, j, fx:0, fy:-mag/nV});
    } else {
        for (let i = 0; i < nU; i++) loads.push({type:'nodal', i, j:nV-1, fx:0, fy:-mag/nU});
    }
    return loads;
};

UDEIMBenchmarkApp.prototype.solve = function(method, mag) {
    const bcs = this.getBCs(), loads = this.getLoads(mag);
    const t0 = performance.now();
    let result, meta = {};

    if (method === 'fom') {
        result = this.solverFOM.solveNonlinear(this.patch, bcs, loads, {iterations:15, steps:3});
    } else if (method === 'galerkin') {
        result = this.romEngine.solveReduced(this.patch, bcs, loads, {iterations:15});
    } else if (method === 'deim') {
        result = this.deimEngine.solveReduced(this.solverFOM, this.romEngine, this.patch, bcs, loads, {iterations:15, steps:1});
        meta.sampled = `${result.sampledCount} elements`;
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

UDEIMBenchmarkApp.prototype.updatePhysics = function() {
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
