import { NURBS2D } from '../../phase-2-core/src/nurbs-2d.js';
import { NURBSPresets } from '../../phase-2-core/src/nurbs-presets.js';
import { RefineUtils } from '../../phase-2-core/src/refine-utils.js';
import { IGANonlinearSolver } from '../../phase-3-core/src/iga-nonlinear.js';
import { ROMEngine } from '../../phase-3-core/src/rom-engine.js';

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
function coolwarm(t) {
    const r = Math.min(1, 2 * t);
    const b = Math.min(1, 2 * (1 - t));
    const g = 1 - 0.6 * Math.abs(t - 0.5) * 2;
    return [r, g, b];
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
        ctx.fillStyle = 'rgba(0,0,0,0.3)';
        ctx.fillRect(0, 0, W, H);
        if (!data.length) return;

        // Log scale
        const vals = data.map(d => Math.max(d.norm, 1e-14));
        const minLog = Math.log10(Math.min(...vals));
        const maxLog = Math.log10(Math.max(...vals)) + 0.5;

        const toY = v => H - ((Math.log10(Math.max(v, 1e-14)) - minLog) / (maxLog - minLog)) * H;
        const step = W / Math.max(data.length - 1, 1);

        // Grid lines
        ctx.strokeStyle = 'rgba(255,255,255,0.05)';
        ctx.lineWidth = 1;
        for (let e = Math.floor(minLog); e <= Math.ceil(maxLog); e++) {
            const y = toY(Math.pow(10, e));
            if (y < 0 || y > H) continue;
            ctx.beginPath(); ctx.moveTo(0, y); ctx.lineTo(W, y); ctx.stroke();
            ctx.fillStyle = '#475569';
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
class ROMApp32 {
    constructor() {
        this.engine = new NURBS2D();
        this.solverFOM = new IGANonlinearSolver(this.engine);
        this.romEngine = new ROMEngine(this.solverFOM);

        // State
        this.currentBenchmark = 'beam';
        this.solverMode = 'fom';   // 'fom' | 'rom'
        this.viewMode = 'disp';    // 'disp' | 'stress'
        this.loadMag = 100;
        this.k = 5;
        this.patch = null;
        this.isTrained = false;
        this.lastResult = null;
        this.lastFomTime = null;
        this._debounceTimer = null;

        // Three.js
        this.scene = new THREE.Scene();
        this.camera = new THREE.PerspectiveCamera(45, window.innerWidth / window.innerHeight, 0.1, 1000);
        this.renderer = new THREE.WebGLRenderer({ antialias: true });
        this.surfaceMesh = null;
        this.wireMesh = null;

        this.initThree();
        this.initUI();
        this.loadBenchmark('beam');
        // No animate() loop — rendering is on-demand only
    }

    // ── Three.js ─────────────────────────────────────────────────────────────
    initThree() {
        this.scene.background = new THREE.Color(0x020617);
        this.camera.position.set(5, 1, 12);
        this.camera.lookAt(5, 1, 0);
        this.renderer.setSize(window.innerWidth, window.innerHeight);
        this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
        document.getElementById('canvas-container').appendChild(this.renderer.domElement);

        this.controls = new THREE.OrbitControls(this.camera, this.renderer.domElement);
        this.controls.enableRotate = false;
        this.controls.target.set(5, 1, 0);
        this.controls.update();

        // On-demand rendering: only render when camera moves
        this.controls.addEventListener('change', () => this._render());

        this.scene.add(new THREE.AmbientLight(0xffffff, 1.0));

        window.addEventListener('resize', () => {
            this.camera.aspect = window.innerWidth / window.innerHeight;
            this.camera.updateProjectionMatrix();
            this.renderer.setSize(window.innerWidth, window.innerHeight);
            this._render(); // Re-render on resize
        });
    }

    // Single render call — only fired on-demand
    _render() {
        this.renderer.render(this.scene, this.camera);
    }

    // ── Benchmark ────────────────────────────────────────────────────────────
    loadBenchmark(type) {
        this.clearMesh();
        this.currentBenchmark = type;
        this.isTrained = false;
        this.lastResult = null;
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
        this.updateMesh(null);
        this.updateOverlays();
        this._render();
    }

    clearMesh() {
        if (this.surfaceMesh) {
            this.scene.remove(this.surfaceMesh);
            this.surfaceMesh.geometry.dispose();
            this.surfaceMesh = null;
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

    // ── BC + Force Arrows  ───────────────────────────────────────────────────
    _makeBCMarker(type) {
        const group = new THREE.Group();
        if (type === 'both') {
            // Hatched wall (vertical line + hatches)
            const pts = [];
            for (let k = 0; k < 5; k++) {
                const y = -0.25 + k * 0.12;
                pts.push(new THREE.Vector3(-0.18, y, 0), new THREE.Vector3(-0.06, y + 0.1, 0));
            }
            pts.push(new THREE.Vector3(-0.12, -0.3, 0), new THREE.Vector3(-0.12, 0.3, 0));
            const geo = new THREE.BufferGeometry().setFromPoints(pts);
            group.add(new THREE.LineSegments(geo, new THREE.LineBasicMaterial({ color: 0x334155 })));

            const coneGeo = new THREE.ConeGeometry(0.09, 0.18, 3);
            const cone = new THREE.Mesh(coneGeo, new THREE.MeshBasicMaterial({ color: 0x334155 }));
            cone.rotation.z = Math.PI / 2;
            cone.position.x = 0.09;
            group.add(cone);
        } else {
            const coneGeo = new THREE.ConeGeometry(0.09, 0.18, 3);
            const cone = new THREE.Mesh(coneGeo, new THREE.MeshBasicMaterial({ color: 0x475569 }));
            if (type === 'y') {
                cone.position.y = -0.09;
                // roller line
                const lg = new THREE.BufferGeometry().setFromPoints([
                    new THREE.Vector3(-0.15, -0.2, 0), new THREE.Vector3(0.15, -0.2, 0)
                ]);
                group.add(new THREE.Line(lg, new THREE.LineBasicMaterial({ color: 0x64748b })));
            } else {
                cone.rotation.z = Math.PI / 2;
                cone.position.x = -0.09;
                const lg = new THREE.BufferGeometry().setFromPoints([
                    new THREE.Vector3(-0.2, -0.15, 0), new THREE.Vector3(-0.2, 0.15, 0)
                ]);
                group.add(new THREE.Line(lg, new THREE.LineBasicMaterial({ color: 0x64748b })));
            }
            group.add(cone);
        }
        return group;
    }

    updateOverlays() {
        this.clearOverlays();
        const cps = this.patch.controlPoints;
        const nU = cps.length, nV = cps[0].length;

        const addMarker = (type, cp) => {
            const m = this._makeBCMarker(type);
            m.position.set(cp.x, cp.y, 0.01);
            this.scene.add(m);
            this._overlayObjects.push(m);
        };

        if (this.currentBenchmark === 'beam') {
            for (let j = 0; j < nV; j++) addMarker('both', cps[0][j]);
            // Force arrows on right edge
            const arrowDir = new THREE.Vector3(0, -1, 0);
            [cps[nU-1][0], cps[nU-1][nV-1]].forEach(cp => {
                const a = new THREE.ArrowHelper(arrowDir,
                    new THREE.Vector3(cp.x, cp.y + 0.8, 0.01), 0.8, 0xef4444, 0.25, 0.12);
                this.scene.add(a);
                this._overlayObjects.push(a);
            });
        } else {
            // Symmetry bottom edge — fix Y
            for (let j = 0; j < nV; j++) addMarker('y', cps[0][j]);
            // Symmetry left edge (last row in U) — fix X
            for (let j = 0; j < nV; j++) addMarker('x', cps[nU-1][j]);
            // Force arrows on right edge (i=0, outer radial)
            const arrowDir = new THREE.Vector3(1, 0, 0);
            [cps[0][nV-1], cps[1][nV-1]].forEach(cp => {
                const a = new THREE.ArrowHelper(arrowDir,
                    new THREE.Vector3(cp.x - 0.8, cp.y, 0.01), 0.8, 0xef4444, 0.25, 0.12);
                this.scene.add(a);
                this._overlayObjects.push(a);
            });
        }
    }

    // ── Physics ───────────────────────────────────────────────────────────────
    getBCs() {
        const bcs = [];
        const nU = this.patch.controlPoints.length;
        const nV = this.patch.controlPoints[0].length;
        if (this.currentBenchmark === 'beam') {
            for (let j = 0; j < nV; j++) bcs.push({ i: 0, j, axis: 'both', value: 0 });
        } else {
            for (let j = 0; j < nV; j++) bcs.push({ i: nU - 1, j, axis: 'x', value: 0 });
            for (let i = 0; i < nU; i++) bcs.push({ i, j: 0, axis: 'y', value: 0 });
        }
        return bcs;
    }

    getLoads(mag) {
        const loads = [];
        const nU = this.patch.controlPoints.length;
        const nV = this.patch.controlPoints[0].length;
        if (this.currentBenchmark === 'beam') {
            for (let j = 0; j < nV; j++) loads.push({ type: 'nodal', i: nU - 1, j, fx: 0, fy: -mag / nV });
        } else {
            for (let j = 0; j < nV; j++) loads.push({ type: 'nodal', i: nU - 1, j, fx: mag / nV, fy: 0 });
        }
        return loads;
    }

    updatePhysics() {
        if (!this.patch) return;
        const bcs = this.getBCs();
        const loads = this.getLoads(this.loadMag);
        const t0 = performance.now();

        let result;
        if (this.solverMode === 'rom' && this.isTrained) {
            result = this.romEngine.solveReduced(this.patch, bcs, loads, { iterations: 15 });
        } else {
            result = this.solverFOM.solveNonlinear(this.patch, bcs, loads, { iterations: 15 });
        }

        const dt = performance.now() - t0;
        this.lastResult = result;

        // Speedup tracking
        if (this.solverMode === 'fom') this.lastFomTime = dt;
        const speedup = (this.solverMode === 'rom' && this.lastFomTime)
            ? (this.lastFomTime / dt).toFixed(1) : '—';

        document.getElementById('time-val').textContent = dt.toFixed(1);
        document.getElementById('speedup-val').textContent = speedup;

        // Sparkline — no animation, instant paint
        this.sparkline.update(result.residualHistory);

        this.updateMesh(result.u);
        this._render(); // Fire exactly one render after physics
    }

    // ── Mesh / Colours ────────────────────────────────────────────────────────
    computeFieldValues(res, uDisp) {
        const positions = [], values = [];
        for (let i = 0; i <= res; i++) {
            const u = Math.min(i / res, 0.9999); // avoid degenerate endpoints
            for (let j = 0; j <= res; j++) {
                const v = Math.min(j / res, 0.9999);
                const state = this.engine.getSurfaceState(this.patch, u, v);
                let px = state.position.x, py = state.position.y;

                let field = 0;
                if (uDisp) {
                    const disp = this.interpolateDisp(u, v, state.denominator, uDisp);
                    px += disp.x; py += disp.y;
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
                positions.push(px, py, 0);
                values.push(field);
            }
        }
        return { positions, values };
    }

    updateMesh(uDisp = null) {
        const res = 36;
        const { positions, values } = this.computeFieldValues(res, uDisp);

        // NaN-safe min/max (a single bad point won't corrupt the whole colormap)
        const finite = values.filter(v => isFinite(v) && v >= 0);
        const minV = finite.length ? Math.min(...finite) : 0;
        const maxV = finite.length ? Math.max(...finite) : 1;
        const range = maxV - minV || 1;
        // Use jet (displacement) or coolwarm (stress) — same as Phase 3.1 palette
        const colorFn = this.viewMode === 'stress' ? coolwarm : jet;
        const colors = [];
        values.forEach(v => {
            const [r, g, b] = uDisp ? colorFn((v - minV) / range) : [0.05, 0.12, 0.22];
            colors.push(r, g, b);
        });

        // Update colorbar
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

        if (!this.surfaceMesh) {
            // First build — allocate geometry
            const geo = new THREE.BufferGeometry();
            const idx = [];
            for (let i = 0; i < res; i++) for (let j = 0; j < res; j++) {
                const a = i * (res + 1) + j, b = (i + 1) * (res + 1) + j;
                const c = (i + 1) * (res + 1) + (j + 1), d = i * (res + 1) + (j + 1);
                idx.push(a, b, d, b, c, d);
            }
            geo.setIndex(idx);
            geo.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
            geo.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));

            this.surfaceMesh = new THREE.Mesh(geo, new THREE.MeshPhongMaterial({
                vertexColors: true, side: THREE.DoubleSide, shininess: 30
            }));
            this.wireMesh = new THREE.LineSegments(
                new THREE.WireframeGeometry(geo),
                new THREE.LineBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0.08 })
            );
            this.scene.add(this.surfaceMesh);
            this.scene.add(this.wireMesh);
        } else {
            // In-place GPU update — no scene add/remove, no flicker
            const posAttr = this.surfaceMesh.geometry.getAttribute('position');
            const colAttr = this.surfaceMesh.geometry.getAttribute('color');
            for (let i = 0; i < positions.length; i++) posAttr.array[i] = positions[i];
            for (let i = 0; i < colors.length; i++) colAttr.array[i] = colors[i];
            posAttr.needsUpdate = true;
            colAttr.needsUpdate = true;
            this.surfaceMesh.geometry.computeBoundingSphere();
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
                const R = (Ni * Mj * this.patch.weights[i][j]) / denom;
                dx += R * uDisp[(i * nV + j) * 2];
                dy += R * uDisp[(i * nV + j) * 2 + 1];
            }
        }
        return { x: dx, y: dy };
    }

    // ── Training ─────────────────────────────────────────────────────────────
    async trainROM() {
        const btn = document.getElementById('btn-train');
        const status = document.getElementById('train-status');
        btn.disabled = true;
        status.classList.remove('hidden');
        this.romEngine.clearSnapshots();

        const bcs = this.getBCs();
        const n = 10;
        for (let i = 1; i <= n; i++) {
            const f = (i / n) * this.loadMag * 2;
            status.textContent = `FOM Step ${i}/${n}  F=${f.toFixed(0)}...`;
            const res = this.solverFOM.solveNonlinear(this.patch, bcs, this.getLoads(f), { steps: 2 });
            this.romEngine.addSnapshot(res.u);
            await new Promise(r => setTimeout(r, 10));
        }
        status.textContent = 'Computing POD (SVD)...';
        const podInfo = this.romEngine.computePOD(this.k);

        document.getElementById('energy-val').textContent = `Energy: ${(podInfo.energy * 100).toFixed(2)}%`;
        document.getElementById('input-k').disabled = false;
        document.getElementById('input-k').max = n;

        this.isTrained = true;
        btn.disabled = false;
        status.textContent = 'ROM Trained ✓';
        setTimeout(() => status.classList.add('hidden'), 3000);
    }

    // ── UI ────────────────────────────────────────────────────────────────────
    initUI() {
        // Sparkline — raw canvas, no lib
        const spCanvas = document.getElementById('sparkline-canvas');
        this.sparkline = new Sparkline(spCanvas);

        // Colorbar
        const cbCanvas = document.getElementById('colorbar-canvas');
        const cbLabel = document.getElementById('colorbar-ticks');
        this.colorbar = new Colorbar(cbCanvas, cbLabel, jet);

        // Benchmark buttons
        document.getElementById('btn-beam').onclick = () => {
            this._setActive('btn-beam', 'btn-plate');
            this.loadBenchmark('beam');
        };
        document.getElementById('btn-plate').onclick = () => {
            this._setActive('btn-plate', 'btn-beam');
            this.loadBenchmark('plate');
        };

        // Train
        document.getElementById('btn-train').onclick = () => this.trainROM();

        // POD k slider
        document.getElementById('input-k').oninput = (e) => {
            this.k = parseInt(e.target.value);
            document.getElementById('k-val').textContent = this.k;
            if (this.isTrained) this.romEngine.computePOD(this.k);
            if (this.solverMode === 'rom') this._scheduleUpdate();
        };

        // Load slider — label instantly, solve debounced
        document.getElementById('input-load').oninput = (e) => {
            this.loadMag = parseFloat(e.target.value);
            document.getElementById('load-val').textContent = this.loadMag;
            this._scheduleUpdate();
        };

        // Solver mode (FOM / ROM)
        document.getElementById('btn-view-fom').onclick = () => {
            this.solverMode = 'fom';
            this._setActive('btn-view-fom', 'btn-view-rom');
            this.updatePhysics();
        };
        document.getElementById('btn-view-rom').onclick = () => {
            if (!this.isTrained) { alert('Generate snapshots first!'); return; }
            this.solverMode = 'rom';
            this._setActive('btn-view-rom', 'btn-view-fom');
            this.updatePhysics();
        };

        // View mode (Displacement / Stress)
        document.getElementById('btn-view-disp').onclick = () => {
            this.viewMode = 'disp';
            this._setActive('btn-view-disp', 'btn-view-stress');
            // Update colorbar gradient
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

new ROMApp32();
