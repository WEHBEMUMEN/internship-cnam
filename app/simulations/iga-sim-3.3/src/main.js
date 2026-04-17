// IGANonlinearSolver and ROMEngine are loaded as global scripts via index.html


// ─── Colour maps ────────────────────────────────────────────────────────────
// Jet (Blue→Cyan→Green→Yellow→Orange→Red) — same as Phase 3.1
const JET_STOPS = [
    { t: 0.0, r: 0.23, g: 0.51, b: 0.96 },
    { t: 0.2, r: 0.02, g: 0.71, b: 0.83 },
    { t: 0.4, r: 0.06, g: 0.73, b: 0.51 },
    { t: 0.6, r: 0.98, g: 0.80, b: 0.08 },
    { t: 0.8, r: 0.98, g: 0.45, b: 0.09 },
    { t: 1.0, r: 0.94, g: 0.27, b: 0.27 }
];
function jet(t) {
    t = Math.max(0, Math.min(1, t));
    for (let i = 0; i < JET_STOPS.length - 1; i++) {
        const s1 = JET_STOPS[i], s2 = JET_STOPS[i + 1];
        if (t >= s1.t && t <= s2.t) {
            const f = (t - s1.t) / (s2.t - s1.t);
            return [s1.r + f*(s2.r-s1.r), s1.g + f*(s2.g-s1.g), s1.b + f*(s2.b-s1.b)];
        }
    }
    return [0.94, 0.27, 0.27];
}
// Stress colormap: Blue → Cyan → Green → Yellow → Red  (no white)
const STRESS_STOPS = [
    { t: 0.00, r: 0.09, g: 0.26, b: 0.85 }, // deep blue
    { t: 0.25, r: 0.02, g: 0.71, b: 0.83 }, // cyan
    { t: 0.50, r: 0.06, g: 0.73, b: 0.31 }, // green
    { t: 0.75, r: 0.98, g: 0.80, b: 0.08 }, // yellow
    { t: 1.00, r: 0.94, g: 0.10, b: 0.10 }, // red
];
function coolwarm(t) {
    t = Math.max(0, Math.min(1, t));
    for (let i = 0; i < STRESS_STOPS.length - 1; i++) {
        const s1 = STRESS_STOPS[i], s2 = STRESS_STOPS[i + 1];
        if (t >= s1.t && t <= s2.t) {
            const f = (t - s1.t) / (s2.t - s1.t);
            return [s1.r + f*(s2.r-s1.r), s1.g + f*(s2.g-s1.g), s1.b + f*(s2.b-s1.b)];
        }
    }
    return [0.94, 0.10, 0.10];
}

// ─── Sparkline (no library — no animation — no flicker) ─────────────────────
class Sparkline {
    constructor(canvas) {
        this.canvas = canvas;
        this.ctx = canvas.getContext('2d');
        this.data = [];
    }
    update(data) {
        this.data = data;
        this.draw();
    }
    draw() {
        const { canvas, ctx, data } = this;
        const W = canvas.width, H = canvas.height;
        ctx.clearRect(0, 0, W, H);
        ctx.fillStyle = 'rgba(241, 245, 249, 0.8)';
        ctx.fillRect(0, 0, W, H);
        if (!data.length) return;

        // Log scale
        const vals = data.map(d => Math.max(d.norm, 1e-14));
        const minLog = Math.log10(Math.min(...vals));
        const maxLog = Math.log10(Math.max(...vals)) + 0.5;

        const toY = v => H - ((Math.log10(Math.max(v, 1e-14)) - minLog) / (maxLog - minLog)) * H;
        const step = W / Math.max(data.length - 1, 1);

        // Grid lines
        ctx.strokeStyle = 'rgba(0,0,0,0.05)';
        ctx.lineWidth = 1;
        for (let e = Math.floor(minLog); e <= Math.ceil(maxLog); e++) {
            const y = toY(Math.pow(10, e));
            if (y < 0 || y > H) continue;
            ctx.beginPath(); ctx.moveTo(0, y); ctx.lineTo(W, y); ctx.stroke();
            ctx.fillStyle = '#64748b';
            ctx.font = '8px monospace';
            ctx.fillText(`10^${e}`, 2, y - 2);
        }

        // Line
        ctx.strokeStyle = '#f43f5e';
        ctx.lineWidth = 1.5;
        ctx.lineJoin = 'round';
        ctx.beginPath();
        data.forEach((d, i) => {
            const x = i * step, y = toY(d.norm);
            i === 0 ? ctx.moveTo(x, y) : ctx.lineTo(x, y);
        });
        ctx.stroke();

        // Dots
        data.forEach((d, i) => {
            const x = i * step, y = toY(d.norm);
            ctx.fillStyle = '#f43f5e';
            ctx.beginPath(); ctx.arc(x, y, 2, 0, Math.PI * 2); ctx.fill();
        });

        // Convergence label
        if (data.length > 0) {
            const last = data[data.length - 1].norm;
            ctx.fillStyle = last < 1e-4 ? '#10b981' : '#f59e0b';
            ctx.font = 'bold 9px monospace';
            ctx.fillText(`Res: ${last.toExponential(2)}`, 4, H - 6);
        }
    }
}

// ─── Colorbar ────────────────────────────────────────────────────────────────
class Colorbar {
    constructor(canvas, labelEl, colorFn) {
        this.canvas = canvas;
        this.label = labelEl;
        this.colorFn = colorFn;
        this._draw();
    }
    _draw() {
        const ctx = this.canvas.getContext('2d');
        const W = this.canvas.width, H = this.canvas.height;
        for (let y = 0; y < H; y++) {
            const t = 1 - y / H;
            const [r, g, b] = this.colorFn(t);
            ctx.fillStyle = `rgb(${(r*255)|0},${(g*255)|0},${(b*255)|0})`;
            ctx.fillRect(0, y, W, 1);
        }
    }
    update(min, max, unit) {
        const mid = (min + max) / 2;
        this.label.innerHTML =
            `<div class="cb-tick">${max.toExponential(1)}</div>` +
            `<div class="cb-tick">${mid.toExponential(1)}</div>` +
            `<div class="cb-tick">${min.toExponential(1)}</div>` +
            `<div class="cb-unit">${unit}</div>`;
    }
}

// ─── Main App ────────────────────────────────────────────────────────────────
class App33 {
    constructor() {
        this.engine = new NURBS2D();
        this.solverFOM = new IGANonlinearSolver(this.engine);

        // State
        this.currentBenchmark = 'twist'; // 'twist' | 'torque'
        this.viewMode = 'disp';    // 'disp' | 'stress'
        this.showMap = false;
        this.loadMag = 45; // degrees or kN
        this.steps = 5;
        this.patch = null;
        this.lastResult = null;
        this._debounceTimer = null;

        // Three.js
        this.scene = new THREE.Scene();
        this.camera = new THREE.PerspectiveCamera(45, window.innerWidth / window.innerHeight, 0.1, 1000);
        this.renderer = new THREE.WebGLRenderer({ antialias: true });
        this.deformedMesh = null;
        this.ghostMesh = null;
        this.wireMesh = null;

        this.initThree();
        this.initUI();
        this.loadBenchmark('twist');
    }

    initThree() {
        this.scene.background = new THREE.Color(0xffffff);
        this.camera.position.set(0, 0, 16);
        this.camera.lookAt(0, 0, 0);
        this.renderer.setSize(window.innerWidth, window.innerHeight);
        this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
        this.renderer.setClearColor(0xffffff, 1);
        document.getElementById('canvas-container').appendChild(this.renderer.domElement);

        this.controls = new THREE.OrbitControls(this.camera, this.renderer.domElement);
        this.controls.enableRotate = false;
        this.controls.target.set(0, 0, 0);
        this.controls.update();

        this.controls.addEventListener('change', () => this._render());

        this.scene.add(new THREE.AmbientLight(0xffffff, 1.2));
        const dirLight = new THREE.DirectionalLight(0xffffff, 0.6);
        dirLight.position.set(5, 5, 10);
        this.scene.add(dirLight);

        window.addEventListener('resize', () => {
            this.camera.aspect = window.innerWidth / window.innerHeight;
            this.camera.updateProjectionMatrix();
            this.renderer.setSize(window.innerWidth, window.innerHeight);
            this._render();
        });
    }

    _render() {
        this.renderer.render(this.scene, this.camera);
    }

    loadBenchmark(type) {
        this.clearMesh();
        this.currentBenchmark = type;
        this.lastResult = null;
        
        // Annulus: Ri=2, Ro=4
        this.patch = NURBSPresets.generateAnnulus(2.0, 4.0);
        RefineUtils.apply(this.engine, this.patch, { p: 3, h: 1 }); // refine slightly

        const loadSlider = document.getElementById('input-load');
        const loadLabelRow = document.getElementById('load-label-text');
        
        if (type === 'twist') {
            loadSlider.max = 180;
            loadSlider.value = 45;
            this.loadMag = 45;
            loadLabelRow.textContent = "Twist Angle θ (deg)";
        } else {
            loadSlider.max = 500000;
            loadSlider.value = 100000;
            this.loadMag = 100000;
            loadLabelRow.textContent = "Applied Torque M (kN·mm)";
        }
        document.getElementById('load-val').textContent = this.loadMag;

        this.updateMesh(null);
        this.updateOverlays();
        this.updatePhysics();
    }

    clearMesh() {
        if (this.deformedMesh) {
            this.scene.remove(this.deformedMesh);
            this.deformedMesh.geometry.dispose();
            this.deformedMesh = null;
        }
        if (this.ghostMesh) {
            this.scene.remove(this.ghostMesh);
            this.ghostMesh.geometry.dispose();
            this.ghostMesh = null;
        }
        if (this.wireMesh) {
            this.scene.remove(this.wireMesh);
            this.wireMesh.geometry.dispose();
            this.wireMesh = null;
        }
        this.clearOverlays();
    }

    clearOverlays() {
        (this._overlayObjects || []).forEach(o => this.scene.remove(o));
        this._overlayObjects = [];
    }

    updateOverlays() {
        this.clearOverlays();
        const cps = this.patch.controlPoints;
        const nU = cps.length;

        // Force arrows on outer boundary
        if (this.currentBenchmark === 'torque') {
            for (let i = 0; i < nU; i += 2) {
                const cp = cps[i][1];
                const r = Math.sqrt(cp.x*cp.x + cp.y*cp.y);
                const tx = -cp.y / r, ty = cp.x / r;
                const arrowDir = new THREE.Vector3(tx, ty, 0);
                const a = new THREE.ArrowHelper(arrowDir,
                    new THREE.Vector3(cp.x, cp.y, 0.01), 1.2, 0xef4444, 0.3, 0.2);
                this.scene.add(a);
                this._overlayObjects.push(a);
            }
        }
    }

    getBCs(magAngleDeg = 0) {
        const bcs = [];
        const nU = this.patch.controlPoints.length;
        
        // Fixed inner boundary
        for (let i = 0; i < nU; i++) {
            bcs.push({ i, j: 0, axis: 'both', value: 0 });
        }

        if (this.currentBenchmark === 'twist') {
            const theta = magAngleDeg * Math.PI / 180;
            for (let i = 0; i < nU; i++) {
                const cp = this.patch.controlPoints[i][1];
                const nx = cp.x * Math.cos(theta) - cp.y * Math.sin(theta);
                const ny = cp.x * Math.sin(theta) + cp.y * Math.cos(theta);
                bcs.push({ i, j: 1, axis: 'x', value: nx - cp.x });
                bcs.push({ i, j: 1, axis: 'y', value: ny - cp.y });
            }
        }
        return bcs;
    }

    getFollowerLoads(totalMag, uDisp) {
        const loads = [];
        const nU = this.patch.controlPoints.length;
        const nV = this.patch.controlPoints[0].length;
        
        for (let i = 0; i < nU; i++) {
            const cp = this.patch.controlPoints[i][1];
            let px = cp.x, py = cp.y;
            if (uDisp) {
                const idx = (i * nV + 1) * 2;
                px += uDisp[idx];
                py += uDisp[idx+1];
            }
            const r = Math.sqrt(px*px + py*py);
            if (r < 1e-12) continue;
            
            const tx = -py / r;
            const ty = px / r;
            
            let nodeForce = totalMag / (nU - 1);
            if (i === 0 || i === nU - 1) nodeForce /= 2;

            loads.push({ type: 'nodal', i, j: 1, fx: nodeForce * tx, fy: nodeForce * ty });
        }
        return loads;
    }

    updatePhysics() {
        if (!this.patch) return;
        const t0 = performance.now();
        const steps = this.steps;
        
        let result = null;
        let u = null;
        let residualHist = [];

        if (this.currentBenchmark === 'twist') {
            // Built-in solver steps handle prescribed increments
            const bcs = this.getBCs(this.loadMag);
            result = this.solverFOM.solveNonlinear(this.patch, bcs, [], { iterations: 15, steps: steps });
            residualHist = result.residualHistory;
        } else {
            // Torque: External increment loop to update follower forces based on current deformation
            const bcs = this.getBCs(0);
            u = new Float64Array(this.patch.controlPoints.length * this.patch.controlPoints[0].length * 2);
            for (let s = 1; s <= steps; s++) {
                const fraction = s / steps;
                const currentLoads = this.getFollowerLoads(this.loadMag * fraction, u);
                
                // One step in the solver since we're stepping externally
                const stepRes = this.solverFOM.solveNonlinear(this.patch, bcs, currentLoads, { iterations: 15, steps: 1, initialU: u });
                
                u = stepRes.u;
                // Accumulate residuals for sparkline
                stepRes.residualHistory.forEach(h => residualHist.push({step: s, iter: h.iter, norm: h.norm}));
                
                if (stepRes.residualHistory.length > 0) {
                     const lastNorm = stepRes.residualHistory[stepRes.residualHistory.length - 1].norm;
                     if (lastNorm > 100 || isNaN(lastNorm)) break; // Divergence stop
                }
            }
            result = { u, residualHistory: residualHist };
        }

        const dt = performance.now() - t0;
        this.lastResult = result;

        document.getElementById('time-val').textContent = dt.toFixed(1);
        
        // Calculate max reaction force for twist mode
        let maxR = 0;
        if (this.currentBenchmark === 'twist') {
            const F_int = this.solverFOM.calculateInternalForce(this.patch, result.u);
            const nU = this.patch.controlPoints.length;
            const nV = this.patch.controlPoints[0].length;
            let sumRx = 0, sumRy = 0;
            for (let i=0; i<nU; i++) {
                sumRx += Math.abs(F_int[(i*nV+1)*2]);
                sumRy += Math.abs(F_int[(i*nV+1)*2+1]);
            }
            maxR = Math.sqrt(sumRx*sumRx + sumRy*sumRy) / 1000; // to kN vaguely
        }
        document.getElementById('reaction-val').textContent = this.currentBenchmark === 'twist' ? maxR.toFixed(1) : (this.loadMag/1000).toFixed(1);

        this.sparkline.update(result.residualHistory);

        const hasBadValues = result.u.some(v => !isFinite(v));
        if (hasBadValues) {
            console.warn('Solver diverged — mesh frozen at last stable state');
            this._render();
            return;
        }

        this.updateMesh(result.u);
        this._render();
    }

    computeFieldValues(res, uDisp) {
        const positions = [], values = [];
        for (let i = 0; i <= res; i++) {
            const u = Math.min(i / res, 0.9999);
            for (let j = 0; j <= res; j++) {
                const v = Math.min(j / res, 0.9999);
                const state = this.engine.getSurfaceState(this.patch, u, v);
                let px = isFinite(state.position.x) ? state.position.x : 0;
                let py = isFinite(state.position.y) ? state.position.y : 0;

                let field = 0;
                if (uDisp) {
                    const disp = this.interpolateDisp(u, v, state.denominator, uDisp);
                    if (isFinite(disp.x) && isFinite(disp.y)) {
                        px += disp.x; py += disp.y;
                    }
                    if (this.viewMode === 'disp') {
                        field = Math.sqrt(disp.x * disp.x + disp.y * disp.y);
                    } else {
                        try {
                            const s = this.solverFOM.getNumericalStress(
                                this.patch, uDisp, u, v,
                                this.solverFOM.E, this.solverFOM.nu);
                            const vm = s.vonMises;
                            field = isFinite(vm) ? vm : 0;
                        } catch (_) { field = 0; }
                    }
                }
                if (!isFinite(px)) px = state.position.x;
                if (!isFinite(py)) py = state.position.y;
                if (!isFinite(field)) field = 0;

                positions.push(px, py, 0);
                values.push(field);
            }
        }
        return { positions, values };
    }

    updateMesh(uDisp = null) {
        const res = 48; // High res for circular object
        const { positions: posDef, values } = this.computeFieldValues(res, uDisp);
        const { positions: posUndef } = this.computeFieldValues(res, null);

        const finite = values.filter(v => isFinite(v) && v >= 0);
        const minV = finite.length ? Math.min(...finite) : 0;
        const maxV = finite.length ? Math.max(...finite) : 1;
        const range = maxV - minV || 1;

        const colorFn = this.viewMode === 'stress' ? coolwarm : jet;
        const colors = [];
        values.forEach(v => {
            if (this.showMap && uDisp) {
                const param = isFinite(v) ? (v - minV) / range : 0;
                const [r, g, b] = colorFn(Math.max(0, Math.min(1, param)));
                colors.push(r, g, b);
            } else {
                colors.push(0.23, 0.51, 0.96);
            }
        });

        const unit = this.viewMode === 'disp' ? 'mm' : 'MPa';
        const label = this.viewMode === 'disp' ? 'Displacement' : 'Von Mises';
        document.getElementById('colorbar-label').textContent = label;
        const colorFnCb = this.viewMode === 'stress' ? coolwarm : jet;
        this.colorbar = new Colorbar(
            document.getElementById('colorbar-canvas'),
            document.getElementById('colorbar-ticks'),
            colorFnCb
        );
        this.colorbar.update(minV, maxV, unit);

        if (!this.deformedMesh) {
            const geoMain = new THREE.BufferGeometry();
            const geoGhost = new THREE.BufferGeometry();
            const idx = [];
            for (let i = 0; i < res; i++) for (let j = 0; j < res; j++) {
                const a = i * (res + 1) + j, b = (i + 1) * (res + 1) + j;
                const c = (i + 1) * (res + 1) + (j + 1), d = i * (res + 1) + (j + 1);
                idx.push(a, b, d, b, c, d);
            }

            geoMain.setIndex(idx);
            geoMain.setAttribute('position', new THREE.Float32BufferAttribute(posDef, 3));
            geoMain.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
            this.deformedMesh = new THREE.Mesh(geoMain, new THREE.MeshStandardMaterial({
                vertexColors: true, side: THREE.DoubleSide, roughness: 0.5, metalness: 0.1
            }));
            this.deformedMesh.geometry.computeVertexNormals();

            geoGhost.setIndex(idx);
            geoGhost.setAttribute('position', new THREE.Float32BufferAttribute(posUndef, 3));
            this.ghostMesh = new THREE.Mesh(geoGhost, new THREE.MeshPhongMaterial({
                color: 0xcbd5e1, transparent: true, opacity: 0.35, side: THREE.DoubleSide
            }));
            this.ghostMesh.position.z = -0.05;
            this.ghostMesh.renderOrder = -1;

            this.wireMesh = new THREE.LineSegments(
                new THREE.WireframeGeometry(geoMain),
                new THREE.LineBasicMaterial({ color: 0x000000, transparent: true, opacity: 0.06 })
            );

            this.scene.add(this.ghostMesh);
            this.scene.add(this.deformedMesh);
            this.scene.add(this.wireMesh);
        } else {
            const posAttr = this.deformedMesh.geometry.getAttribute('position');
            const colAttr = this.deformedMesh.geometry.getAttribute('color');
            for (let i = 0; i < posDef.length; i++) posAttr.array[i] = posDef[i];
            for (let i = 0; i < colors.length; i++) colAttr.array[i] = colors[i];
            posAttr.needsUpdate = true;
            colAttr.needsUpdate = true;
            this.deformedMesh.geometry.computeVertexNormals();
            this.deformedMesh.geometry.computeBoundingSphere();

            this.scene.remove(this.wireMesh);
            this.wireMesh = new THREE.LineSegments(
                new THREE.WireframeGeometry(this.deformedMesh.geometry),
                new THREE.LineBasicMaterial({ color: 0x000000, transparent: true, opacity: 0.06 })
            );
            this.scene.add(this.wireMesh);
        }
    }

    interpolateDisp(u, v, denom, uDisp) {
        const nU = this.patch.controlPoints.length, nV = this.patch.controlPoints[0].length;
        let dx = 0, dy = 0;
        for (let i = 0; i < nU; i++) {
            const Ni = this.engine.basis1D(i, this.patch.p, this.patch.U, u);
            if (!Ni) continue;
            for (let j = 0; j < nV; j++) {
                const Mj = this.engine.basis1D(j, this.patch.q, this.patch.V, v);
                const w = (this.patch.weights && this.patch.weights[i]) ? this.patch.weights[i][j] : 1.0;
                const R = (Math.abs(denom) > 1e-12) ? (Ni * Mj * w) / denom : 0;
                dx += R * uDisp[(i * nV + j) * 2];
                dy += R * uDisp[(i * nV + j) * 2 + 1];
            }
        }
        return { x: dx, y: dy };
    }

    initUI() {
        const spCanvas = document.getElementById('sparkline-canvas');
        this.sparkline = new Sparkline(spCanvas);

        const cbCanvas = document.getElementById('colorbar-canvas');
        const cbLabel = document.getElementById('colorbar-ticks');
        this.colorbar = new Colorbar(cbCanvas, cbLabel, jet);

        document.getElementById('btn-twist').onclick = () => {
            this._setActive('btn-twist', 'btn-torque');
            this.loadBenchmark('twist');
        };
        document.getElementById('btn-torque').onclick = () => {
            this._setActive('btn-torque', 'btn-twist');
            this.loadBenchmark('torque');
        };

        const loadSlider = document.getElementById('input-load');
        loadSlider.oninput = (e) => {
            this.loadMag = parseFloat(e.target.value);
            document.getElementById('load-val').textContent = this.loadMag;
            this._scheduleUpdate();
        };

        const stepsSlider = document.getElementById('input-steps');
        stepsSlider.oninput = (e) => {
            this.steps = parseInt(e.target.value);
            document.getElementById('steps-val').textContent = this.steps;
            this._scheduleUpdate();
        };

        document.getElementById('btn-view-disp').onclick = () => {
            this.viewMode = 'disp';
            this._setActive('btn-view-disp', 'btn-view-stress');
            this.colorbar = new Colorbar(
                document.getElementById('colorbar-canvas'),
                document.getElementById('colorbar-ticks'), jet);
            if (this.lastResult) { this.updateMesh(this.lastResult.u); this._render(); }
        };
        document.getElementById('btn-view-stress').onclick = () => {
            this.viewMode = 'stress';
            this._setActive('btn-view-stress', 'btn-view-disp');
            this.colorbar = new Colorbar(
                document.getElementById('colorbar-canvas'),
                document.getElementById('colorbar-ticks'), coolwarm);
            if (this.lastResult) { this.updateMesh(this.lastResult.u); this._render(); }
        };

        document.getElementById('btn-map-on').onclick = () => {
            this.showMap = true;
            this._setActive('btn-map-on', 'btn-map-off');
            document.getElementById('viz-selectors').classList.replace('hidden-fade', 'visible-fade');
            document.getElementById('colorbar-panel').classList.replace('hidden-fade', 'visible-fade');
            this.updateMesh(this.lastResult ? this.lastResult.u : null);
            this._render();
        };
        document.getElementById('btn-map-off').onclick = () => {
            this.showMap = false;
            this._setActive('btn-map-off', 'btn-map-on');
            document.getElementById('viz-selectors').classList.replace('visible-fade', 'hidden-fade');
            document.getElementById('colorbar-panel').classList.replace('visible-fade', 'hidden-fade');
            this.updateMesh(this.lastResult ? this.lastResult.u : null);
            this._render();
        };
    }

    _scheduleUpdate() {
        clearTimeout(this._debounceTimer);
        this._debounceTimer = setTimeout(() => this.updatePhysics(), 160);
    }

    _setActive(onId, offId) {
        document.getElementById(onId).classList.add('active');
        document.getElementById(offId).classList.remove('active');
    }
}

new App33();
