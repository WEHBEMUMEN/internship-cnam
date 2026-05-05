/**
 * DEIM Benchmark — Visualizations & Mesh Rendering
 */

DEIMBenchmarkApp.prototype.initCharts = function() {
    this.sparkline = new Sparkline(document.getElementById('sparkline-canvas'));
    this.speedupChart = new Chart(document.getElementById('chart-speedup'), {
        type: 'bar',
        data: { labels: [], datasets: [{ label: 'Speedup (×)', data: [], backgroundColor: [
            'rgba(148, 163, 184, 0.5)', // FOM
            'rgba(14, 165, 233, 0.5)',  // Galerkin
            'rgba(139, 92, 246, 0.5)'   // DEIM
        ], borderColor: ['#64748b', '#0ea5e9', '#8b5cf6'], borderWidth: 1 }] },
        options: { responsive:true, maintainAspectRatio:false, plugins:{legend:{display:false}},
            scales:{y:{beginAtZero:true, title:{display:true, text:'Speedup (×)', font:{size:10}}}, x:{ticks:{font:{size:9}}}}}
    });
    this.fdChart = new Chart(document.getElementById('chart-fd'), {
        type: 'line', data: { labels: [], datasets: [] },
        options: { responsive:true, maintainAspectRatio:false,
            plugins:{legend:{position:'top', labels:{font:{size:9}}}},
            scales:{x:{title:{display:true, text:'Load F', font:{size:10}}}, y:{title:{display:true, text:'Tip Displacement', font:{size:10}}}}}
    });

    document.querySelectorAll('.chart-tab').forEach(t => {
        t.onclick = () => {
            document.querySelectorAll('.chart-tab').forEach(b => b.classList.remove('active'));
            t.classList.add('active');
            document.getElementById('chart-speedup-wrap').classList.toggle('hidden', t.dataset.chart !== 'speedup');
            document.getElementById('chart-fd-wrap').classList.toggle('hidden', t.dataset.chart !== 'fd');
            document.getElementById('chart-error-wrap').classList.toggle('hidden', t.dataset.chart !== 'error');
            document.getElementById('chart-points-wrap').classList.toggle('hidden', t.dataset.chart !== 'points');
            document.getElementById('explorer-wrap').classList.toggle('hidden', t.dataset.chart !== 'explorer');
            
            if ((t.dataset.chart === 'fd' || t.dataset.chart === 'error') && this.isTrained) this.runFDCurves();
            if (t.dataset.chart === 'explorer') {
                this.isExplorerActive = true;
                this.runDEIMExplorer(this.explorerStep);
            } else {
                this.isExplorerActive = false;
                this.sensorSpheres.clear();
                if (this.spanBorders) {
                    this.scene.remove(this.spanBorders);
                    this.spanBorders = null;
                }
                if (this.activeElementMesh) this.activeElementMesh.visible = false;
                if (this.lastResult) this.updateMesh(this.lastResult.u);
                this._render();
            }
        };
    });

    this.errorChart = new Chart(document.getElementById('chart-error'), {
        type: 'line', data: { labels: [], datasets: [] },
        options: { responsive:true, maintainAspectRatio:false,
            plugins:{legend:{position:'top', labels:{font:{size:9}}}},
            scales:{x:{title:{display:true, text:'Load F', font:{size:10}}}, y:{type: 'logarithmic', title:{display:true, text:'L2 Rel Error', font:{size:10}}}}}
    });

    this.pointsChart = new Chart(document.getElementById('chart-points'), {
        type: 'line',
        data: { labels: [], datasets: [
            { label: 'Max Interpolation Error', data: [], borderColor: '#f59e0b', backgroundColor: 'rgba(245,158,11,0.2)', borderWidth: 2, tension: 0.1, fill: true, pointRadius: 3 }
        ]},
        options: { responsive:true, maintainAspectRatio:false, plugins:{legend:{position:'top', labels:{font:{size:9}}}},
            scales:{x:{title:{display:true, text:'Number of DEIM Points (m)', font:{size:10}}}, y:{type:'logarithmic', title:{display:true, text:'Max Error', font:{size:10}}}}
        }
    });

    this.residualChart = new Chart(document.getElementById('chart-residual'), {
        type: 'line',
        data: { labels: [], datasets: [
            { label: 'Residual', data: [], borderColor: '#f59e0b', backgroundColor: 'rgba(245,158,11,0.2)', borderWidth: 2, tension: 0.1, fill: true, pointRadius: 0 },
            { label: 'Selected Node', data: [], pointBackgroundColor: '#ef4444', pointBorderColor: '#fff', pointRadius: 5, pointHoverRadius: 7, showLine: false }
        ]},
        options: { responsive:true, maintainAspectRatio:false, plugins:{legend:{display:false}}, animation: false,
            scales:{x:{display:false}, y:{type:'logarithmic', title:{display:true, text:'Magnitude', font:{size:9}}}}
        }
    });
};

DEIMBenchmarkApp.prototype.updateSpeedupChart = function(data) {
    const labels = ['FOM', 'Galerkin', 'DEIM'];
    const values = [1.0];
    
    // Galerkin
    if (data.galerkin && this.lastFomTime) {
        values.push(this.lastFomTime / data.galerkin.time);
    } else values.push(0);

    // DEIM
    if (data.deim && this.lastFomTime) {
        values.push(this.lastFomTime / data.deim.time);
    } else values.push(0);

    this.speedupChart.data.labels = labels;
    this.speedupChart.data.datasets[0].data = values;
    this.speedupChart.update();
};

DEIMBenchmarkApp.prototype._updateStats = function(method, meta) {
    document.getElementById('method-name').textContent = METHODS[method].label;
    document.getElementById('time-val').textContent = meta.time.toFixed(1);
    const speedup = this.lastFomTime ? (this.lastFomTime / meta.time).toFixed(1) : '—';
    document.getElementById('speedup-val').textContent = method === 'fom' ? '1.0' : speedup;
    document.getElementById('sampled-val').textContent = meta.sampled || 'All';
    document.getElementById('tip-val').textContent = meta.tipDisp ? Math.abs(meta.tipDisp).toFixed(3) : '—';
    document.getElementById('error-val').textContent = meta.error !== undefined ? (meta.error*100).toFixed(2)+'%' : '—';
};

DEIMBenchmarkApp.prototype.updateMesh = function(uDisp) {
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
            } else {
                colors.push(0.23, 0.51, 0.96);
            }
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
        this.ghostMesh = new THREE.Mesh(gD, new THREE.MeshPhongMaterial({color:0xcbd5e1, transparent:true, opacity:0.3, side:THREE.DoubleSide}));
        
        // Add Wireframe for element visibility
        this.wireframe = new THREE.LineSegments(
            new THREE.WireframeGeometry(gD),
            new THREE.LineBasicMaterial({ color: 0x475569, transparent: true, opacity: 0.2 })
        );
        
        this.scene.add(this.ghostMesh);
        this.scene.add(this.deformedMesh);
        this.scene.add(this.wireframe);
    } else {
        const pa = this.deformedMesh.geometry.getAttribute('position');
        const ca = this.deformedMesh.geometry.getAttribute('color');
        for (let i = 0; i < posDef.length; i++) pa.array[i] = posDef[i];
        for (let i = 0; i < colors.length; i++) ca.array[i] = colors[i];
        pa.needsUpdate = true; ca.needsUpdate = true;
        this.deformedMesh.geometry.computeVertexNormals();
        
        // Update wireframe geometry too
        this.wireframe.geometry.dispose();
        this.wireframe.geometry = new THREE.WireframeGeometry(this.deformedMesh.geometry);
    }
};

DEIMBenchmarkApp.prototype._interpDisp = function(u, v, denom, uDisp) {
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
