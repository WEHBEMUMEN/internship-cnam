/**
 * Phase 3.4b: U-DEIM Benchmark — Visuals & Rendering
 */

UDEIMBenchmarkApp.prototype._render = function() { this.renderer.render(this.scene, this.camera); };

UDEIMBenchmarkApp.prototype.updateMesh = function(uDisp) {
    const res = 32;
    const posDef = [], posUnd = [], colors = [];
    let maxF = 0;

    for (let i = 0; i <= res; i++) {
        const u = Math.min(i/res, 0.9999);
        for (let j = 0; j <= res; j++) {
            const v = Math.min(j/res, 0.9999);
            const st = this.engine.getSurfaceState(this.patch, u, v);
            let px = st.position.x, py = st.position.y;
            posUnd.push(px, py, 0);
            let field = 0;
            if (uDisp) {
                const d = this._interpDisp(u, v, st.denominator, uDisp);
                if (isFinite(d.x) && isFinite(d.y)) { px += d.x; py += d.y; }
                field = Math.sqrt(d.x*d.x + d.y*d.y);
                maxF = Math.max(maxF, field);
            }
            posDef.push(px, py, 0);
        }
    }

    for (let i = 0; i <= res; i++) {
        const u = Math.min(i/res, 0.9999);
        for (let j = 0; j <= res; j++) {
            const v = Math.min(j/res, 0.9999);
            if (uDisp && maxF > 0) {
                const st = this.engine.getSurfaceState(this.patch, u, v);
                const d = this._interpDisp(u, v, st.denominator, uDisp);
                const f = Math.sqrt(d.x*d.x + d.y*d.y) / maxF;
                const [r,g,b] = jet(f);
                colors.push(r, g, b);
            } else colors.push(0.23, 0.51, 0.96);
        }
    }

    if (!this.deformedMesh) {
        const idx = [];
        for (let i = 0; i < res; i++) for (let j = 0; j < res; j++) {
            const a = i*(res+1)+j, b = (i+1)*(res+1)+j, c = (i+1)*(res+1)+j+1, d = i*(res+1)+j+1;
            idx.push(a,b,d,b,c,d);
        }
        const gD = new THREE.BufferGeometry(); gD.setIndex(idx);
        gD.setAttribute('position', new THREE.Float32BufferAttribute(posDef, 3));
        gD.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
        this.deformedMesh = new THREE.Mesh(gD, new THREE.MeshStandardMaterial({vertexColors:true, side:THREE.DoubleSide, roughness:0.5, metalness:0.1}));
        this.deformedMesh.geometry.computeVertexNormals();

        const gG = new THREE.BufferGeometry(); gG.setIndex(idx);
        gG.setAttribute('position', new THREE.Float32BufferAttribute(posUnd, 3));
        this.ghostMesh = new THREE.Mesh(gG, new THREE.MeshPhongMaterial({color:0xcbd5e1, transparent:true, opacity:0.3, side:THREE.DoubleSide}));
        this.ghostMesh.position.z = -0.05;
        this.scene.add(this.ghostMesh);
        this.scene.add(this.deformedMesh);
    } else {
        const pa = this.deformedMesh.geometry.getAttribute('position');
        const ca = this.deformedMesh.geometry.getAttribute('color');
        for (let i = 0; i < posDef.length; i++) pa.array[i] = posDef[i];
        for (let i = 0; i < colors.length; i++) ca.array[i] = colors[i];
        pa.needsUpdate = true; ca.needsUpdate = true;
        this.deformedMesh.geometry.computeVertexNormals();
    }

    while(this.loadArrows.children.length > 0) this.loadArrows.remove(this.loadArrows.children[0]);
    const loads = this.getLoads(this.loadMag);
    const dir = new THREE.Vector3(0, -1, 0);
    loads.forEach(load => {
        if(load.fy === 0) return;
        let px = this.patch.controlPoints[load.i][load.j].x;
        let py = this.patch.controlPoints[load.i][load.j].y;
        if (uDisp) {
            const idx = (load.i * this.patch.controlPoints[0].length + load.j) * 2;
            px += uDisp[idx]; py += uDisp[idx+1];
        }
        const length = Math.min(2.0, Math.max(0.5, Math.abs(load.fy) * 0.05));
        const origin = new THREE.Vector3(px, py + length, 0.1);
        this.loadArrows.add(new THREE.ArrowHelper(dir, origin, length, 0xef4444, 0.3, 0.2));
    });

    if (this.dofPoints) { this.scene.remove(this.dofPoints); this.dofPoints.geometry.dispose(); this.dofPoints = null; }
    if (this.showDOFs) {
        const dofGeo = new THREE.BufferGeometry();
        const dofPos = [];
        for (let i = 0; i < this.patch.controlPoints.length; i++)
            for (let j = 0; j < this.patch.controlPoints[0].length; j++) {
                let px = this.patch.controlPoints[i][j].x, py = this.patch.controlPoints[i][j].y;
                if (uDisp) { const idx = (i * this.patch.controlPoints[0].length + j) * 2; px += uDisp[idx]; py += uDisp[idx+1]; }
                dofPos.push(px, py, 0.1);
            }
        dofGeo.setAttribute('position', new THREE.Float32BufferAttribute(dofPos, 3));
        this.dofPoints = new THREE.Points(dofGeo, new THREE.PointsMaterial({ color: 0x10b981, size: 0.15 }));
        this.scene.add(this.dofPoints);
    }
};

UDEIMBenchmarkApp.prototype._interpDisp = function(u, v, denom, uDisp) {
    const nU = this.patch.controlPoints.length, nV = this.patch.controlPoints[0].length;
    let dx = 0, dy = 0;
    for (let i = 0; i < nU; i++) {
        const Ni = this.engine.basis1D(i, this.patch.p, this.patch.U, u);
        if (!Ni) continue;
        for (let j = 0; j < nV; j++) {
            const Mj = this.engine.basis1D(j, this.patch.q, this.patch.V, v);
            const w = this.patch.weights[i][j] || 1;
            const R = Math.abs(denom) > 1e-12 ? (Ni*Mj*w)/denom : 0;
            dx += R * uDisp[(i*nV+j)*2];
            dy += R * uDisp[(i*nV+j)*2+1];
        }
    }
    return {x: dx, y: dy};
};

UDEIMBenchmarkApp.prototype.runDEIMExplorer = function(step) {
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
        const eIdx = Math.floor(dofIdx / this.deimEngine.numLocalDofs);
        const ld = dofIdx % this.deimEngine.numLocalDofs;
        const el = this.deimEngine.activeElements.find(e => e.index === eIdx);
        if (!el) return;
        const map = this.deimEngine.elementDofMap[eIdx];
        const globalDof = map[ld];
        const cpIdx = Math.floor(globalDof / 2);
        const cpI = Math.floor(cpIdx / nV), cpJ = cpIdx % nV;
        const u = this.patch.U[cpI + this.patch.p], v = this.patch.V[cpJ + this.patch.q];
        const st = this.engine.getSurfaceState(this.patch, u, v);
        const isNewest = (idx === currentIndices.length - 1);
        const geo = new THREE.SphereGeometry(isNewest ? 0.3 : 0.2, 16, 16);
        const mat = new THREE.MeshBasicMaterial({ color: isNewest ? 0xef4444 : 0xf59e0b });
        const sphere = new THREE.Mesh(geo, mat);
        sphere.position.set(st.position.x, st.position.y, 0.1);
        this.sensorSpheres.add(sphere);
    });

    const activeSpans = [];
    currentIndices.forEach((dofIdx) => {
        const eIdx = Math.floor(dofIdx / this.deimEngine.numLocalDofs);
        const el = this.deimEngine.activeElements.find(e => e.index === eIdx);
        if (el) activeSpans.push(el);
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
                if (isActive) colors.push(0.98, 0.8, 0.1); else colors.push(0.9, 0.9, 0.9);
            }
        }
        this.deformedMesh.geometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
        this.deformedMesh.geometry.attributes.color.needsUpdate = true;
    }
    this._render();
};
