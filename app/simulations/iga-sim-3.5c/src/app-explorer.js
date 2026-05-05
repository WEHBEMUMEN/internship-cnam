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
    const previousSpans = [], newestSpans = [];
    
    currentIndices.forEach((dofIdx, idx) => {
        const cpIdx = Math.floor(dofIdx / 2);
        const cpI = Math.floor(cpIdx / nV);
        const cpJ = cpIdx % nV;
        const isNewest = (idx === currentIndices.length - 1);
        
        for (let i = 0; i < uniqueU.length - 1; i++) {
            const uMid = (uniqueU[i] + uniqueU[i+1])/2;
            const spanU = this.engine.findSpan(this.patch.controlPoints.length - 1, this.patch.p, uMid, this.patch.U);
            if (cpI >= spanU - this.patch.p && cpI <= spanU) {
                for (let j = 0; j < uniqueV.length - 1; j++) {
                    const vMid = (uniqueV[j] + uniqueV[j+1])/2;
                    const spanV = this.engine.findSpan(this.patch.controlPoints[0].length - 1, this.patch.q, vMid, this.patch.V);
                    if (cpJ >= spanV - this.patch.q && cpJ <= spanV) {
                        const span = { uMin: uniqueU[i], uMax: uniqueU[i+1], vMin: uniqueV[j], vMax: uniqueV[j+1] };
                        if (isNewest) newestSpans.push(span);
                        else previousSpans.push(span);
                    }
                }
            }
        }
    });

    if (this.deformedMesh) {
        const colors = [];
        const res = 32;
        this.updateMesh(null);
        
        for (let i = 0; i <= res; i++) {
            const u = Math.min(i/res, 0.9999);
            for (let j = 0; j <= res; j++) {
                const v = Math.min(j/res, 0.9999);
                
                let isNew = false, isPrev = false;
                for (const span of newestSpans) {
                    if (u >= span.uMin && u < span.uMax + 1e-6 && v >= span.vMin && v < span.vMax + 1e-6) { isNew = true; break; }
                }
                if (!isNew) {
                    for (const span of previousSpans) {
                        if (u >= span.uMin && u < span.uMax + 1e-6 && v >= span.vMin && v < span.vMax + 1e-6) { isPrev = true; break; }
                    }
                }

                if (isNew) colors.push(1.0, 0.84, 0.0);       // Bright Gold
                else if (isPrev) colors.push(0.96, 0.5, 0.15); // Deep Orange
                else colors.push(0.85, 0.88, 0.9);           // Slate grey
            }
        }
        this.deformedMesh.geometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
        this.deformedMesh.geometry.attributes.color.needsUpdate = true;
        this.deformedMesh.visible = true;
        
        // Phase 3.5 Polish: Hide original wireframe and draw span borders
        if (this.wireframe) this.wireframe.visible = false;
        this._updateSpanBorders();
    }
    this._render();
};

DEIMBenchmarkApp.prototype.runECSWExplorer = function(step) {
    if (!this.ecswEngine.sampleElements) return;
    
    // For ECSW, we don't always have a step-by-step history like DEIM,
    // so we just show the active elements.
    const activeElements = this.ecswEngine.sampleElements;
    document.getElementById('exp-step-label').textContent = `All Active Elements (${activeElements.length})`;
    
    const uniqueU = [...new Set(this.patch.U)], uniqueV = [...new Set(this.patch.V)];
    const activeSpans = activeElements.map(el => ({
        uMin: uniqueU[el.i], uMax: uniqueU[el.i+1],
        vMin: uniqueV[el.j], vMax: uniqueV[el.j+1]
    }));

    if (this.deformedMesh) {
        const colors = [];
        const res = 32;
        this.updateMesh(null);
        
        for (let i = 0; i <= res; i++) {
            const u = Math.min(i/res, 0.9999);
            for (let j = 0; j <= res; j++) {
                const v = Math.min(j/res, 0.9999);
                
                let isActive = false;
                for (const span of activeSpans) {
                    if (u >= span.uMin && u < span.uMax + 1e-6 && v >= span.vMin && v < span.vMax + 1e-6) { isActive = true; break; }
                }

                if (isActive) colors.push(0.1, 0.8, 0.4); // ECSW Green
                else colors.push(0.85, 0.88, 0.9);         // Slate grey
            }
        }
        this.deformedMesh.geometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
        this.deformedMesh.geometry.attributes.color.needsUpdate = true;
        this.deformedMesh.visible = true;
        
        if (this.wireframe) this.wireframe.visible = false;
        this._updateSpanBorders();
    }
    this._render();
};

DEIMBenchmarkApp.prototype._updateSpanBorders = function() {
    if (this.spanBorders) {
        this.scene.remove(this.spanBorders);
        this.spanBorders = null;
    }
    if (!this.isExplorerActive) return;

    const uniqueU = [...new Set(this.patch.U)], uniqueV = [...new Set(this.patch.V)];
    
    // Support both DEIM and ECSW in border rendering
    let activeSpansSet = new Set();
    if (this.method === 'ecsw' && this.ecswEngine) {
        this.ecswEngine.sampleElements.forEach(el => activeSpansSet.add(`${el.i},${el.j}`));
    } else if (this.deimEngine && this.deimEngine.history) {
        const activeIndices = this.deimEngine.history.map(h => h.point);
        const nV = this.patch.controlPoints[0].length;
        activeIndices.forEach(dofIdx => {
            const cpIdx = Math.floor(dofIdx / 2);
            const cpI = Math.floor(cpIdx / nV);
            const cpJ = cpIdx % nV;
            for (let i = 0; i < uniqueU.length - 1; i++) {
                const uMid = (uniqueU[i] + uniqueU[i+1])/2;
                const spanU = this.engine.findSpan(this.patch.controlPoints.length - 1, this.patch.p, uMid, this.patch.U);
                if (cpI >= spanU - this.patch.p && cpI <= spanU) {
                    for (let j = 0; j < uniqueV.length - 1; j++) {
                        const vMid = (uniqueV[j] + uniqueV[j+1])/2;
                        const spanV = this.engine.findSpan(this.patch.controlPoints[0].length - 1, this.patch.q, vMid, this.patch.V);
                        if (cpJ >= spanV - this.patch.q && cpJ <= spanV) activeSpansSet.add(`${i},${j}`);
                    }
                }
            }
        });
    }

    const activeLines = [], inactiveLines = [];
    for (let i = 0; i < uniqueU.length - 1; i++) {
        for (let j = 0; j < uniqueV.length - 1; j++) {
            const corners = [
                {u: uniqueU[i], v: uniqueV[j]}, {u: uniqueU[i+1], v: uniqueV[j]},
                {u: uniqueU[i+1], v: uniqueV[j+1]}, {u: uniqueU[i], v: uniqueV[j+1]}
            ].map(c => {
                const st = this.engine.getSurfaceState(this.patch, Math.min(c.u, 0.999), Math.min(c.v, 0.999));
                return new THREE.Vector3(st.position.x, st.position.y, 0.11); 
            });
            const isActive = activeSpansSet.has(`${i},${j}`);
            const list = isActive ? activeLines : inactiveLines;
            for (let k = 0; k < 4; k++) {
                list.push(corners[k].x, corners[k].y, corners[k].z);
                list.push(corners[(k+1)%4].x, corners[(k+1)%4].y, corners[(k+1)%4].z);
            }
        }
    }

    const group = new THREE.Group();
    if (activeLines.length) {
        const geo = new THREE.BufferGeometry();
        geo.setAttribute('position', new THREE.Float32BufferAttribute(activeLines, 3));
        const mat = new THREE.LineBasicMaterial({ color: 0x000000 });
        group.add(new THREE.LineSegments(geo, mat));
        const geo2 = geo.clone(); geo2.translate(0.01, 0.01, 0);
        group.add(new THREE.LineSegments(geo2, mat));
    }
    if (inactiveLines.length) {
        const geo = new THREE.BufferGeometry();
        geo.setAttribute('position', new THREE.Float32BufferAttribute(inactiveLines, 3));
        const mat = new THREE.LineBasicMaterial({ color: 0x94a3b8, transparent: true, opacity: 0.3 });
        group.add(new THREE.LineSegments(geo, mat));
    }
    this.spanBorders = group;
    this.scene.add(this.spanBorders);
};
