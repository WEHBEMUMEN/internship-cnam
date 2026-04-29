/**
 * DEIM Benchmark — DEIM Explorer Component
 */

DEIMBenchmarkApp.prototype.runDEIMExplorer = function(step) {
    if (!this.deimEngine.history || this.deimEngine.history.length === 0) return;
    const h = this.deimEngine.history[step - 1];
    
    document.getElementById('exp-step-label').textContent = `${step} / ${this.deimEngine.history.length}`;
    document.getElementById('btn-exp-prev').className = step === 1 ? 'btn-pill' : 'btn-pill active';
    document.getElementById('btn-exp-next').className = step === this.deimEngine.history.length ? 'btn-pill' : 'btn-pill active';

    const labels = Array.from({length: h.residual.length}, (_, i) => i);
    const selData = Array(h.residual.length).fill(null);
    selData[h.point] = h.residual[h.point];
    
    this.residualChart.data.labels = labels;
    this.residualChart.data.datasets[0].data = h.residual;
    this.residualChart.data.datasets[1].data = selData;
    this.residualChart.update();

    this.sensorSpheres.clear();
    this.scene.add(this.sensorSpheres);
    
    const nV = this.patch.controlPoints[0].length;
    const currentIndices = this.deimEngine.history.slice(0, step).map(s => s.point);
    
    currentIndices.forEach((dofIdx, idx) => {
        const cpIdx = Math.floor(dofIdx / 2);
        const cpI = Math.floor(cpIdx / nV);
        const cpJ = cpIdx % nV;
        const u = this.patch.U[cpI + this.patch.p];
        const v = this.patch.V[cpJ + this.patch.q];
        const st = this.engine.getSurfaceState(this.patch, u, v);
        const isNewest = (idx === currentIndices.length - 1);
        
        const geo = new THREE.SphereGeometry(isNewest ? 0.3 : 0.2, 16, 16);
        const mat = new THREE.MeshBasicMaterial({ color: isNewest ? 0xef4444 : 0xf59e0b });
        const sphere = new THREE.Mesh(geo, mat);
        sphere.position.set(st.position.x, st.position.y, 0.1);
        this.sensorSpheres.add(sphere);
    });

    const uniqueU = [...new Set(this.patch.U)], uniqueV = [...new Set(this.patch.V)];
    const activeSpans = [];
    currentIndices.forEach((dofIdx) => {
        const cpIdx = Math.floor(dofIdx / 2);
        const cpI = Math.floor(cpIdx / nV);
        const cpJ = cpIdx % nV;
        for (let eI = cpI; eI <= cpI + this.patch.p; eI++) {
            for (let eJ = cpJ; eJ <= cpJ + this.patch.q; eJ++) {
                if (eI >= 0 && eI < uniqueU.length - 1 && eJ >= 0 && eJ < uniqueV.length - 1) {
                    activeSpans.push({ uMin: uniqueU[eI], uMax: uniqueU[eI+1], vMin: uniqueV[eJ], vMax: uniqueV[eJ+1] });
                }
            }
        }
    });

    if (this.deformedMesh) {
        const colors = [];
        const res = 32;
        for (let i = 0; i <= res; i++) {
            const u = Math.min(i/res, 0.9999);
            for (let j = 0; j <= res; j++) {
                const v = Math.min(j/res, 0.9999);
                let isActive = false;
                for (const span of activeSpans) {
                    if (u >= span.uMin && u < span.uMax + 1e-6 && v >= span.vMin && v < span.vMax + 1e-6) { isActive = true; break; }
                }
                colors.push(isActive ? 0.98 : 0.9, isActive ? 0.8 : 0.9, isActive ? 0.1 : 0.9);
            }
        }
        this.deformedMesh.geometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
        this.deformedMesh.geometry.attributes.color.needsUpdate = true;
    }
    this._render();
};
