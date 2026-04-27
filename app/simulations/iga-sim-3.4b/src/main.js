// Phase 3.4b: ECSW Benchmark — Main Application

const JET = [{t:0,r:.23,g:.51,b:.96},{t:.2,r:.02,g:.71,b:.83},{t:.4,r:.06,g:.73,b:.51},{t:.6,r:.98,g:.8,b:.08},{t:.8,r:.98,g:.45,b:.09},{t:1,r:.94,g:.27,b:.27}];
function jet(t){t=Math.max(0,Math.min(1,t));for(let i=0;i<JET.length-1;i++){const a=JET[i],b=JET[i+1];if(t>=a.t&&t<=b.t){const f=(t-a.t)/(b.t-a.t);return[a.r+f*(b.r-a.r),a.g+f*(b.g-a.g),a.b+f*(b.b-a.b)];}}return[.94,.27,.27];}

class Sparkline{constructor(c){this.c=c;this.ctx=c.getContext('2d');this.data=[];}
update(d){this.data=d;this.draw();}
draw(){const{c,ctx,data}=this;const W=c.width,H=c.height;ctx.clearRect(0,0,W,H);ctx.fillStyle='rgba(241,245,249,0.8)';ctx.fillRect(0,0,W,H);if(!data.length)return;const vals=data.map(d=>Math.max(d.norm,1e-14));const mn=Math.log10(Math.min(...vals)),mx=Math.log10(Math.max(...vals))+.5;const toY=v=>H-((Math.log10(Math.max(v,1e-14))-mn)/(mx-mn))*H;const step=W/Math.max(data.length-1,1);ctx.strokeStyle='rgba(0,0,0,0.05)';ctx.lineWidth=1;for(let e=Math.floor(mn);e<=Math.ceil(mx);e++){const y=toY(Math.pow(10,e));if(y<0||y>H)continue;ctx.beginPath();ctx.moveTo(0,y);ctx.lineTo(W,y);ctx.stroke();ctx.fillStyle='#64748b';ctx.font='8px monospace';ctx.fillText(`10^${e}`,2,y-2);}ctx.strokeStyle='#10b981';ctx.lineWidth=1.5;ctx.lineJoin='round';ctx.beginPath();data.forEach((d,i)=>{const x=i*step,y=toY(d.norm);i===0?ctx.moveTo(x,y):ctx.lineTo(x,y);});ctx.stroke();data.forEach((d,i)=>{ctx.fillStyle='#10b981';ctx.beginPath();ctx.arc(i*step,toY(d.norm),2,0,Math.PI*2);ctx.fill();});if(data.length){const last=data[data.length-1].norm;ctx.fillStyle=last<1e-4?'#059669':'#f59e0b';ctx.font='bold 9px monospace';ctx.fillText(`Res: ${last.toExponential(2)}`,4,H-6);}}}

const METHODS = {
    fom:      {label:'FOM',      color:'#334155'},
    galerkin: {label:'Galerkin', color:'#3b82f6'},
    ecsw:     {label:'ECSW',     color:'#10b981'}
};

class ECSWBenchmarkApp {
    constructor() {
        this.engine = new NURBS2D();
        this.solverFOM = new IGANonlinearSolver(this.engine);
        this.romEngine = new ROMEngine(this.solverFOM);
        this.ecswEngine = new ECSWEngine();

        this.patch = null;
        this.method = 'fom';
        this.loadMag = 200;
        this.k = 5;
        this.ecswM = 20;
        this.meshLevel = 1;
        this.loadType = 'tip';
        this.showDOFs = true;
        this.isTrained = false;
        this.lastResult = null;
        this.lastFomTime = null;
        this.lastFomResult = null;
        this._timer = null;

        this.scene = new THREE.Scene();
        this.loadArrows = new THREE.Group();
        this.scene.add(this.loadArrows);
        this.dofPoints = null;
        this.camera = new THREE.PerspectiveCamera(45, window.innerWidth/window.innerHeight, 0.1, 1000);
        this.renderer = new THREE.WebGLRenderer({antialias:true});
        this.deformedMesh = null;
        this.ghostMesh = null;

        this.initThree();
        this.initUI();
        this.initCharts();
        this.loadBenchmark();
    }

    initThree() {
        this.scene.background = new THREE.Color(0xffffff);
        this.camera.position.set(5, 1, 14);
        this.camera.lookAt(5, 1, 0);
        this.renderer.setSize(window.innerWidth, window.innerHeight);
        this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
        document.getElementById('canvas-container').appendChild(this.renderer.domElement);
        this.controls = new THREE.OrbitControls(this.camera, this.renderer.domElement);
        this.controls.enableRotate = false;
        this.controls.target.set(5, 1, 0);
        this.controls.update();
        this.controls.addEventListener('change', () => this._render());
        this.scene.add(new THREE.AmbientLight(0xffffff, 1.2));
        const dl = new THREE.DirectionalLight(0xffffff, 0.6);
        dl.position.set(5, 5, 10);
        this.scene.add(dl);
        window.addEventListener('resize', () => {
            this.camera.aspect = window.innerWidth/window.innerHeight;
            this.camera.updateProjectionMatrix();
            this.renderer.setSize(window.innerWidth, window.innerHeight);
            this._render();
        });
    }

    _render() { this.renderer.render(this.scene, this.camera); }

    loadBenchmark() {
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
    }

    clearMesh() {
        [this.deformedMesh, this.ghostMesh].forEach(m => { if (m) { this.scene.remove(m); m.geometry.dispose(); } });
        this.deformedMesh = this.ghostMesh = null;
    }

    getBCs() {
        const nV = this.patch.controlPoints[0].length;
        const bcs = [];
        for (let j = 0; j < nV; j++) bcs.push({i:0, j, axis:'both', value:0});
        return bcs;
    }

    getLoads(mag) {
        const nU = this.patch.controlPoints.length, nV = this.patch.controlPoints[0].length;
        const loads = [];
        if (this.loadType === 'tip') {
            for (let j = 0; j < nV; j++) loads.push({type:'nodal', i:nU-1, j, fx:0, fy:-mag/nV});
        } else {
            for (let i = 0; i < nU; i++) loads.push({type:'nodal', i, j:nV-1, fx:0, fy:-mag/nU});
        }
        return loads;
    }

    solve(method, mag) {
        const bcs = this.getBCs(), loads = this.getLoads(mag);
        const t0 = performance.now();
        let result, meta = {};

        if (method === 'fom') {
            result = this.solverFOM.solveNonlinear(this.patch, bcs, loads, {iterations:15, steps:3});
        } else if (method === 'galerkin') {
            result = this.romEngine.solveReduced(this.patch, bcs, loads, {iterations:15});
        } else if (method === 'ecsw') {
            result = this.ecswEngine.solveReduced(this.solverFOM, this.patch, loads, {iterations:15, steps:3});
            meta.sampled = `${this.ecswEngine.sampleElements.length} / ${this.ecswEngine.weights.length} elems`;
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
    }

    updatePhysics() {
        if (!this.patch) return;
        const m = METHODS[this.method];
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
    }

    _updateStats(method, meta) {
        document.getElementById('method-name').textContent = METHODS[method].label;
        document.getElementById('time-val').textContent = meta.time.toFixed(1);
        const speedup = this.lastFomTime ? (this.lastFomTime / meta.time).toFixed(1) : '—';
        document.getElementById('speedup-val').textContent = method === 'fom' ? '1.0' : speedup;
        document.getElementById('sampled-val').textContent = meta.sampled || 'All';
        document.getElementById('tip-val').textContent = meta.tipDisp ? Math.abs(meta.tipDisp).toFixed(3) : '—';
        document.getElementById('error-val').textContent = meta.error !== undefined ? (meta.error*100).toFixed(2)+'%' : '—';
    }

    updateMesh(uDisp) {
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
                    colors.push(0.06, 0.73, 0.51); // ECSW Greenish
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
            this.deformedMesh.geometry.computeVertexNormals();

            const gG = new THREE.BufferGeometry(); gG.setIndex(idx);
            gG.setAttribute('position', new THREE.BufferAttribute(new Float32Array(posUnd), 3));
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
    }

    _interpDisp(u, v, denom, uDisp) {
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
    }

    async trainAll() {
        const btn = document.getElementById('btn-train');
        const status = document.getElementById('train-status');
        btn.disabled = true;
        status.classList.remove('hidden');
        this.romEngine.clearSnapshots();

        const bcs = this.getBCs();
        const nSnaps = 15;
        const snapDisp = [];

        for (let i = 1; i <= nSnaps; i++) {
            const f = (i/nSnaps) * this.loadMag * 1.5;
            status.textContent = `FOM Snapshot ${i}/${nSnaps} (F=${f.toFixed(0)})...`;
            const res = this.solverFOM.solveNonlinear(this.patch, bcs, this.getLoads(f), {steps:3, iterations:12});
            this.romEngine.addSnapshot(res.u);
            snapDisp.push(res.u);
            await new Promise(r => setTimeout(r, 5));
        }

        status.textContent = 'Computing POD basis...';
        await new Promise(r => setTimeout(r, 10));
        const podInfo = this.romEngine.computePOD(this.k);
        document.getElementById('energy-val').textContent = `Energy: ${(podInfo.energy*100).toFixed(2)}%`;
        document.getElementById('input-k').disabled = false;
        document.getElementById('input-k').max = nSnaps;

        status.textContent = `Training ECSW Weights...`;
        await new Promise(r => setTimeout(r, 10));
        const ecswInfo = await this.ecswEngine.train(this.solverFOM, this.romEngine, this.patch, snapDisp);
        document.getElementById('input-m').disabled = false;
        document.getElementById('ecsw-info').textContent = `${ecswInfo.elementCount} elements active`;

        this.isTrained = true;
        this.snapDisp = snapDisp;
        btn.disabled = false;
        document.getElementById('btn-compare').disabled = false;
        status.textContent = `Training complete ✓ — ${ecswInfo.elementCount} elements, POD: ${(podInfo.energy*100).toFixed(1)}% energy`;
        setTimeout(() => status.classList.add('hidden'), 5000);
        this.updatePhysics();
    }

    async runComparison() {
        if (!this.isTrained) return;
        document.getElementById('btn-compare').disabled = true;
        const data = {};
        for (const m of Object.keys(METHODS)) {
            try {
                const { result, meta } = this.solve(m, this.loadMag);
                data[m] = meta;
            } catch(e) { console.warn(`${m} failed:`, e); }
            await new Promise(r => setTimeout(r, 5));
        }
        this.updateSpeedupChart(data);
        document.getElementById('btn-compare').disabled = false;
    }

    initCharts() {
        this.sparkline = new Sparkline(document.getElementById('sparkline-canvas'));
        this.speedupChart = new Chart(document.getElementById('chart-speedup'), {
            type: 'bar',
            data: { labels: [], datasets: [{ label: 'Speedup (×)', data: [], backgroundColor: [] }] },
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
                if (t.dataset.chart === 'fd' && this.isTrained) this.runFDCurves();
            };
        });
    }

    updateSpeedupChart(data) {
        const labels = [], values = [], bgColors = [];
        const fomTime = data.fom?.time || 1;
        for (const [m, meta] of Object.entries(data)) {
            labels.push(METHODS[m].label);
            values.push(parseFloat((fomTime / meta.time).toFixed(1)));
            bgColors.push(METHODS[m].color + '99');
        }
        this.speedupChart.data.labels = labels;
        this.speedupChart.data.datasets[0].data = values;
        this.speedupChart.data.datasets[0].backgroundColor = bgColors;
        this.speedupChart.update();
    }

    async runFDCurves() {
        if (!this.isTrained) return;
        const loads = [20, 50, 100, 150, 200, 300, 400, 500];
        const datasets = {};
        for (const m of Object.keys(METHODS)) {
            datasets[m] = [];
            for (const f of loads) {
                try { const { meta } = this.solve(m, f); datasets[m].push(Math.abs(meta.tipDisp || 0)); }
                catch(e) { datasets[m].push(0); }
            }
        }
        this.fdChart.data.labels = loads;
        this.fdChart.data.datasets = Object.keys(METHODS).map(m => ({
            label: METHODS[m].label, data: datasets[m],
            borderColor: METHODS[m].color, backgroundColor: METHODS[m].color + '33',
            borderWidth: 2, pointRadius: 3, tension: 0.3
        }));
        this.fdChart.update();
    }

    initUI() {
        document.getElementById('btn-train').onclick = () => this.trainAll();
        document.getElementById('btn-compare').onclick = () => this.runComparison();

        document.querySelectorAll('[data-method]').forEach(btn => {
            btn.onclick = () => {
                document.querySelectorAll('[data-method]').forEach(b => b.classList.remove('active'));
                btn.classList.add('active');
                this.method = btn.dataset.method;
                if (this.method === 'fom' || this.isTrained) this.updatePhysics();
            };
        });

        document.getElementById('input-k').oninput = e => {
            this.k = parseInt(e.target.value);
            document.getElementById('k-val').textContent = this.k;
            if (this.isTrained) {
                this.romEngine.computePOD(this.k);
                this._scheduleUpdate();
            }
        };

        document.getElementById('input-mesh').oninput = e => {
            this.meshLevel = parseInt(e.target.value);
            document.getElementById('mesh-val').textContent = this.meshLevel;
            this.loadBenchmark();
        };

        document.getElementById('input-load-type').onchange = e => {
            this.loadType = e.target.value;
            this.lastFomResult = null;
            this._scheduleUpdate();
        };

        document.getElementById('input-show-dofs').onchange = e => {
            this.showDOFs = e.target.checked;
            if (this.lastResult) this.updateMesh(this.lastResult.u);
            else this.updateMesh(null);
            this._render();
        };

        document.getElementById('input-load').oninput = e => {
            this.loadMag = parseFloat(e.target.value);
            document.getElementById('load-val').textContent = this.loadMag;
            this.lastFomResult = null;
            this._scheduleUpdate();
        };
    }

    _scheduleUpdate() {
        clearTimeout(this._timer);
        this._timer = setTimeout(() => this.updatePhysics(), 200);
    }
}

new ECSWBenchmarkApp();
