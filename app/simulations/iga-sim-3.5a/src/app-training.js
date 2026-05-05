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

    status.textContent = `Training DEIM (m=${this.deimM}, kf=${this.k})...`;
    await new Promise(r => setTimeout(r, 10));
    const deimInfo = this.deimEngine.train(forceSnaps, this.deimM, this.k, bcs_dofs);
    document.getElementById('input-m').disabled = false;
    document.getElementById('deim-info').textContent = `${deimInfo.m} interpolation points`;

    status.textContent = 'Mapping active elements...';
    this.deimEngine.computeActiveElements(this.patch);
    status.textContent = 'Pre-computing reduced tangent...';
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
    
    // Update Points vs Error chart with TRUE Reconstruction Error
    if (this.pointsChart && this.deimEngine && this.deimEngine.history) {
        const labels = [];
        const maxErrors = [];
        
        const PhiT = this.romEngine.Phi.transpose();
        const PhiT_Uf = PhiT.mmul(this.deimEngine.U_f);
        const { Matrix, SVD } = window.mlMatrix;

        // Precalculate true force projections to save time
        const F_true_proj = [];
        const trueNorms = [];
        this.forceSnaps.forEach(snap => {
            const F_proj = new Float64Array(this.romEngine.Phi.columns);
            let normSq = 0;
            for (let i = 0; i < F_proj.length; i++) {
                let dot = 0;
                for (let d = 0; d < snap.length; d++) dot += PhiT.get(i, d) * snap[d];
                F_proj[i] = dot;
                normSq += dot * dot;
            }
            F_true_proj.push(F_proj);
            trueNorms.push(Math.sqrt(normSq));
        });

        // Compute error for each m
        for (let m = 1; m <= this.deimEngine.m; m++) {
            const kf = Math.min(this.k, m);
            const indices = this.deimEngine.indices.slice(0, m);
            
            const PtU_arr = Array.from({ length: m }, () => new Float64Array(kf));
            for (let i = 0; i < m; i++)
                for (let j = 0; j < kf; j++)
                    PtU_arr[i][j] = this.deimEngine.U_f.get(indices[i], j);
            
            const PtU_mat = new Matrix(PtU_arr);
            const svd = new SVD(PtU_mat, { computeLeftSingularVectors: true, computeRightSingularVectors: true });
            const U = svd.leftSingularVectors;
            const V = svd.rightSingularVectors;
            const sigmas = svd.diagonal;
            
            const Sinv = new Matrix(V.columns, U.columns);
            const tol = sigmas[0] * 1e-12;
            for (let i = 0; i < sigmas.length; i++) {
                if (sigmas[i] > tol) Sinv.set(i, i, 1.0 / sigmas[i]);
            }
            const PtU_pinv_mat = V.mmul(Sinv).mmul(U.transpose());
            const PtU_pinv = PtU_pinv_mat.to2DArray();

            let maxRelErr = 0;
            for (let s = 0; s < this.forceSnaps.length; s++) {
                const f_sampled = new Float64Array(m);
                for (let j = 0; j < m; j++) f_sampled[j] = this.forceSnaps[s][indices[j]];

                const c = new Float64Array(kf);
                for (let i = 0; i < kf; i++) {
                    let sum = 0;
                    for (let j = 0; j < m; j++) sum += PtU_pinv[i][j] * f_sampled[j];
                    c[i] = sum;
                }
                
                const F_proj_deim = new Float64Array(this.romEngine.Phi.columns);
                for (let i = 0; i < F_proj_deim.length; i++) {
                    let sum = 0;
                    for (let j = 0; j < kf; j++) sum += PhiT_Uf.get(i, j) * c[j];
                    F_proj_deim[i] = sum;
                }

                let err2 = 0;
                for (let i = 0; i < F_proj_deim.length; i++) {
                    err2 += (F_true_proj[s][i] - F_proj_deim[i]) ** 2;
                }
                const relErr = trueNorms[s] > 1e-30 ? Math.sqrt(err2) / trueNorms[s] : 0;
                maxRelErr = Math.max(maxRelErr, relErr);
            }
            labels.push(m);
            // Cap at 100% so initial terrible errors don't squash the chart
            maxErrors.push(Math.min(100, maxRelErr * 100));
        }

        this.pointsChart.data.labels = labels;
        this.pointsChart.data.datasets[0].data = maxErrors;
        this.pointsChart.data.datasets[0].label = 'Force Reconstruction Error (%)';
        this.pointsChart.update();
    }
    
    this.updatePhysics();
    this.runOnlineAudit();
};
