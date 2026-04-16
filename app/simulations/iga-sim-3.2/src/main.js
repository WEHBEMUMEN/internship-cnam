import { NURBS2D } from '../../phase-2-core/src/nurbs-2d.js';
import { NURBSPresets } from '../../phase-2-core/src/nurbs-presets.js';
import { RefineUtils } from '../../phase-2-core/src/refine-utils.js';
import { IGANonlinearSolver } from '../../phase-3-core/src/iga-nonlinear.js';
import { ROMEngine } from '../../phase-3-core/src/rom-engine.js';

class ROMApp32 {
    constructor() {
        this.engine = new NURBS2D();
        this.solverFOM = new IGANonlinearSolver(this.engine);
        this.romEngine = new ROMEngine(this.solverFOM);
        
        // State
        this.currentBenchmark = 'beam';
        this.viewMode = 'fom';
        this.loadMag = 100;
        this.k = 5;
        this.patch = null;
        this.isTrained = false;
        
        // Three.js
        this.scene = new THREE.Scene();
        this.camera = new THREE.PerspectiveCamera(45, window.innerWidth / window.innerHeight, 0.1, 1000);
        this.renderer = new THREE.WebGLRenderer({ antialias: true });
        this.initThree();
        
        // UI
        this.initCharts();
        this.setupEventListeners();
        
        this.loadBenchmark('beam');
        this.animate();
    }

    initThree() {
        this.scene.background = new THREE.Color(0x020617);
        // Default for Cantilever (L=10)
        this.camera.position.set(5, 1, 12);
        this.camera.lookAt(5, 1, 0);
        this.renderer.setSize(window.innerWidth, window.innerHeight);
        this.renderer.setPixelRatio(window.devicePixelRatio);
        document.getElementById('canvas-container').appendChild(this.renderer.domElement);
        
        this.controls = new THREE.OrbitControls(this.camera, this.renderer.domElement);
        this.controls.enableRotate = false; // Lock to 2D
        this.controls.target.set(5, 1, 0);
        
        const ambient = new THREE.AmbientLight(0xffffff, 1.0);
        this.scene.add(ambient);
    }

    initCharts() {
        const ctx = document.getElementById('metrics-chart').getContext('2d');
        this.chart = new Chart(ctx, {
            type: 'line',
            data: { labels: [], datasets: [{ label: 'Convergence Norm', data: [], borderColor: '#f43f5e', tension: 0.1 }] },
            options: {
                scales: { 
                    y: { type: 'logarithmic', grid: { color: 'rgba(255,255,255,0.05)' } },
                    x: { grid: { display: false } }
                },
                plugins: { legend: { display: false } },
                maintainAspectRatio: false
            }
        });
    }

    loadBenchmark(type) {
        this.currentBenchmark = type;
        this.isTrained = false;
        this.romEngine.clearSnapshots();
        document.getElementById('input-k').disabled = true;
        document.getElementById('train-status').classList.add('hidden');
        
        if (type === 'beam') {
            this.patch = NURBSPresets.generateCantilever(10, 2);
            RefineUtils.apply(this.engine, this.patch, { p: 3, h: 1 });
            this.camera.position.set(5, 1, 12);
            this.controls.target.set(5, 1, 0);
        } else {
            this.patch = NURBSPresets.generatePlateWithHole(1.0, 4.0);
            RefineUtils.apply(this.engine, this.patch, { p: 3, h: 0 });
            this.patch = this.engine.elevateDegree(this.patch);
            this.patch = this.engine.subdivideGlobal(this.patch);
            this.camera.position.set(2, 2, 8);
            this.controls.target.set(2, 2, 0);
        }
        this.camera.updateProjectionMatrix();
        this.controls.update();
        
        this.updateMesh();
    }

    async trainROM() {
        const btn = document.getElementById('btn-train');
        const status = document.getElementById('train-status');
        btn.disabled = true;
        status.classList.remove('hidden');
        status.textContent = "Sampling Parameter Space...";
        
        this.romEngine.clearSnapshots();
        
        // BCs & Loads Setup
        const bcs = this.getBCs();
        const maxLoad = this.loadMag * 2;
        const nSnapshots = 10;
        
        for (let i = 1; i <= nSnapshots; i++) {
            const stepLoad = (i / nSnapshots) * maxLoad;
            status.textContent = `Solving FOM Step ${i}/${nSnapshots} (F=${stepLoad.toFixed(0)})...`;
            
            const results = this.solverFOM.solveNonlinear(this.patch, bcs, this.getLoads(stepLoad), { steps: 2 });
            this.romEngine.addSnapshot(results.u);
            await new Promise(r => setTimeout(r, 10));
        }
        
        status.textContent = "Computing POD Basis (SVD)...";
        const podInfo = this.romEngine.computePOD(this.k);
        
        document.getElementById('energy-val').textContent = `Energy: ${(podInfo.energy * 100).toFixed(4)}%`;
        document.getElementById('input-k').disabled = false;
        document.getElementById('input-k').max = nSnapshots;
        
        this.isTrained = true;
        btn.disabled = false;
        status.textContent = "ROM Trained Successfully!";
        setTimeout(() => status.classList.add('hidden'), 3000);
    }

    getBCs() {
        const bcs = [];
        const nBasisU = this.patch.controlPoints.length;
        const nBasisV = this.patch.controlPoints[0].length;
        
        if (this.currentBenchmark === 'beam') {
            for (let j = 0; j < nBasisV; j++) bcs.push({ i: 0, j, axis: 'both', value: 0 });
        } else {
            // Half symmetry
            for (let j = 0; j < nBasisV; j++) bcs.push({ i: nBasisU - 1, j, axis: 'x', value: 0 });
            for (let i = 0; i < nBasisU; i++) bcs.push({ i, j: 0, axis: 'y', value: 0 });
        }
        return bcs;
    }

    getLoads(mag) {
        const loads = [];
        const nBasisU = this.patch.controlPoints.length;
        const nBasisV = this.patch.controlPoints[0].length;
        
        if (this.currentBenchmark === 'beam') {
            for (let j = 0; j < nBasisV; j++) loads.push({ i: nBasisU - 1, j, fx: 0, fy: -mag / nBasisV });
        } else {
            // Right traction
            for (let j = 0; j < nBasisV; j++) loads.push({ i: nBasisU - 1, j, fx: mag / nBasisV, fy: 0 });
        }
        return loads;
    }

    updatePhysics() {
        const t0 = performance.now();
        let result;
        
        const bcs = this.getBCs();
        const loads = this.getLoads(this.loadMag);

        const onProgress = (p) => {
            this.chart.data.labels.push(p.iter);
            this.chart.data.datasets[0].data.push(p.norm);
            this.chart.update();
        };

        this.chart.data.labels = [];
        this.chart.data.datasets[0].data = [];

        if (this.viewMode === 'rom' && this.isTrained) {
            result = this.romEngine.solveReduced(this.patch, bcs, loads, { onProgress });
        } else {
            result = this.solverFOM.solveNonlinear(this.patch, bcs, loads, { onProgress });
        }
        
        const t1 = performance.now();
        const dt = t1 - t0;
        
        document.getElementById('time-val').textContent = dt.toFixed(2);
        
        // Update speedup display (rough estimate if we have a baseline)
        if (this.viewMode === 'fom') this.lastFomTime = dt;
        if (this.viewMode === 'rom' && this.lastFomTime) {
            const speedup = (this.lastFomTime / dt).toFixed(1);
            document.getElementById('speedup-val').textContent = speedup;
        }

        this.updateMesh(result.u);
    }

    updateMesh(uDisp = null) {
        if (this.surfaceMesh) this.scene.remove(this.surfaceMesh);
        if (this.wireMesh) this.scene.remove(this.wireMesh);

        const res = 40;
        const geometry = new THREE.BufferGeometry();
        const positions = [];
        const colors = [];

        for (let i = 0; i <= res; i++) {
            const u = i / res;
            for (let j = 0; j <= res; j++) {
                const v = j / res;
                const state = this.engine.getSurfaceState(this.patch, u, v);
                let p = { ...state.position };

                if (uDisp) {
                    const disp = this.interpolateDisplacement(u, v, state.denominator, uDisp);
                    p.x += disp.x; p.y += disp.y;
                    
                    const mag = Math.sqrt(disp.x**2 + disp.y**2) * 5;
                    colors.push(0.2, 0.4 + mag, 0.8);
                } else {
                    colors.push(0.1, 0.2, 0.3);
                }
                positions.push(p.x, p.y, p.z);
            }
        }

        const indices = [];
        for (let i = 0; i < res; i++) {
            for (let j = 0; j < res; j++) {
                const a = i * (res + 1) + j;
                const b = (i + 1) * (res + 1) + j;
                const c = (i + 1) * (res + 1) + (j + 1);
                const d = i * (res + 1) + (j + 1);
                indices.push(a, b, d, b, c, d);
            }
        }

        geometry.setIndex(indices);
        geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
        geometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
        
        this.surfaceMesh = new THREE.Mesh(geometry, new THREE.MeshPhongMaterial({ vertexColors: true, side: THREE.DoubleSide, shininess: 30 }));
        this.wireMesh = new THREE.LineSegments(new THREE.WireframeGeometry(geometry), new THREE.LineBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0.1 }));
        
        this.scene.add(this.surfaceMesh);
        this.scene.add(this.wireMesh);
    }

    interpolateDisplacement(u, v, denom, u_disp) {
        const nU = this.patch.controlPoints.length;
        const nV = this.patch.controlPoints[0].length;
        let dx = 0, dy = 0;
        for (let i = 0; i < nU; i++) {
            const Ni = this.engine.basis1D(i, this.patch.p, this.patch.U, u);
            if (Ni === 0) continue;
            for (let j = 0; j < nV; j++) {
                const Mj = this.engine.basis1D(j, this.patch.q, this.patch.V, v);
                const R = (Ni * Mj * this.patch.weights[i][j]) / denom;
                dx += R * u_disp[(i * nV + j) * 2];
                dy += R * u_disp[(i * nV + j) * 2 + 1];
            }
        }
        return { x: dx, y: dy };
    }

    setupEventListeners() {
        document.getElementById('btn-beam').onclick = (e) => {
            document.getElementById('btn-beam').classList.add('active');
            document.getElementById('btn-plate').classList.remove('active');
            this.loadBenchmark('beam');
        };
        document.getElementById('btn-plate').onclick = (e) => {
            document.getElementById('btn-plate').classList.add('active');
            document.getElementById('btn-beam').classList.remove('active');
            this.loadBenchmark('plate');
        };
        
        document.getElementById('btn-train').onclick = () => this.trainROM();
        
        document.getElementById('input-k').oninput = (e) => {
            this.k = parseInt(e.target.value);
            document.getElementById('k-val').textContent = this.k;
            if (this.isTrained) this.romEngine.computePOD(this.k);
            if (this.viewMode === 'rom') this.updatePhysics();
        };

        document.getElementById('input-load').oninput = (e) => {
            this.loadMag = parseFloat(e.target.value);
            document.getElementById('load-val').textContent = this.loadMag;
            this.updatePhysics();
        };

        document.getElementById('btn-view-fom').onclick = () => {
            this.viewMode = 'fom';
            document.getElementById('btn-view-fom').classList.add('active');
            document.getElementById('btn-view-rom').classList.remove('active');
            this.updatePhysics();
        };

        document.getElementById('btn-view-rom').onclick = () => {
            this.viewMode = 'rom';
            document.getElementById('btn-view-rom').classList.add('active');
            document.getElementById('btn-view-fom').classList.remove('active');
            if (this.isTrained) this.updatePhysics();
            else alert("Please generate snapshots first!");
        };

        window.addEventListener('resize', () => {
            this.camera.aspect = window.innerWidth / window.innerHeight;
            this.camera.updateProjectionMatrix();
            this.renderer.setSize(window.innerWidth, window.innerHeight);
        });
    }

    animate() {
        requestAnimationFrame(() => this.animate());
        this.renderer.render(this.scene, this.camera);
    }
}

new ROMApp32();
