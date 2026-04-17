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
class ROMApp32 {
    constructor() {
        this.engine = new NURBS2D();
        this.solverFOM = new IGANonlinearSolver(this.engine);
        this.romEngine = new ROMEngine(this.solverFOM);

        // State
        this.currentBenchmark = 'beam';
        this.solverMode = 'fom';   // 'fom' | 'rom'
        this.viewMode = 'disp';    // 'disp' | 'stress'
        this.showMap = false;     // Toggle for colormapping
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
        this.deformedMesh = null;
        this.ghostMesh = null;
        this.wireMesh = null;

        this.initThree();
        this.initUI();
        this.loadBenchmark('beam');
        // No animate() loop — rendering is on-demand only
    }

    // ── Three.js ─────────────────────────────────────────────────────────────
    initThree() {
        this.scene.background = new THREE.Color(0xffffff);
        this.camera.position.set(5, 1, 12);
        this.camera.lookAt(5, 1, 0);
        this.renderer.setSize(window.innerWidth, window.innerHeight);
        this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
        this.renderer.setClearColor(0xffffff, 1);
        document.getElementById('canvas-container').appendChild(this.renderer.domElement);

        this.controls = new THREE.OrbitControls(this.camera, this.renderer.domElement);
        this.controls.enableRotate = false;
        this.controls.target.set(5, 1, 0);
        this.controls.update();

        // On-demand rendering: only render when camera moves
        this.controls.addEventListener('change', () => this._render());

        this.scene.add(new THREE.AmbientLight(0xffffff, 1.2));
        const dirLight = new THREE.DirectionalLight(0xffffff, 0.6);
        dirLight.position.set(5, 5, 10);
        this.scene.add(dirLight);

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
            RefineUtils.apply(this.engine, this.patch, { p: 2, h: 0 }); // Match Phase 2
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
            // Symmetry bottom edge (i=0) — fix Y
            for (let j = 0; j < nV; j++) addMarker('y', cps[0][j]);
            // Symmetry left edge (i=nU-1) — fix X
            for (let j = 0; j < nV; j++) addMarker('x', cps[nU-1][j]);
            
            // Force arrows on outer boundary pulling right (j=nV-1)
            const arrowDir = new THREE.Vector3(1, 0, 0);
            // Targeted arrows along the vertical part of the outer edge
            for (let i = 0; i < nU; i++) {
                const cp = cps[i][nV - 1];
                if (cp.x > 0.1) { // Only on vertical or transition parts, avoid left axis
                    const a = new THREE.ArrowHelper(arrowDir,
                        new THREE.Vector3(cp.x - 0.8, cp.y, 0.01), 0.8, 0xef4444, 0.25, 0.12);
                    this.scene.add(a);
                    this._overlayObjects.push(a);
                }
            }
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
            // Symmetry bottom edge (i=0) - Fix Y
            for (let j = 0; j < nV; j++) bcs.push({ i: 0, j, axis: 'y', value: 0 });
            // Symmetry left edge (i=nU-1) - Fix X
            for (let j = 0; j < nV; j++) bcs.push({ i: nU - 1, j, axis: 'x', value: 0 });
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
            // Apply high-fidelity Kirsch analytical traction to truncated boundary
            const nodalForces = this._calculateNodalTraction(this.patch, mag);
            for (let i = 0; i < nU; i++) {
                const idx = (i * nV + (nV - 1)) * 2;
                loads.push({ type: 'nodal', i: i, j: nV - 1, fx: nodalForces[idx], fy: nodalForces[idx + 1] });
            }
        }
        return loads;
    }

    _getExactStress(x, y, R, Tx) {
        const r = Math.sqrt(x * x + y * y);
        const theta = Math.atan2(y, x);
        const R2 = R * R, R4 = R2 * R2, r2 = r * r, r4 = r2 * r2;
        const cos2t = Math.cos(2 * theta), sin2t = Math.sin(2 * theta);

        const s_rr = (Tx / 2) * (1 - R2/r2) + (Tx / 2) * (1 - 4*R2/r2 + 3*R4/r4) * cos2t;
        const s_tt = (Tx / 2) * (1 + R2/r2) - (Tx / 2) * (1 + 3*R4/r4) * cos2t;
        const s_rt = -(Tx / 2) * (1 + 2*R2/r2 - 3*R4/r4) * sin2t;

        const c = Math.cos(theta), s = Math.sin(theta);
        const sxx = s_rr * c*c + s_tt * s*s - 2 * s_rt * s*c;
        const syy = s_rr * s*s + s_tt * c*c + 2 * s_rt * s*c;
        const sxy = (s_rr - s_tt) * s * c + s_rt * (c*c - s*s);
        return { sxx, syy, sxy };
    }

    _calculateNodalTraction(patch, loadValue) {
        const { p, q, U, V, controlPoints, weights } = patch;
        const nU = controlPoints.length;
        const nV = controlPoints[0].length;
        const forces = new Float64Array(nU * nV * 2);
        
        const gPoints = [-0.7745966692414833, 0.0, 0.7745966692414833];
        const gWeights = [0.5555555555555556, 0.8888888888888888, 0.5555555555555556];
        const uniqueU = [...new Set(U)];
        const vParam = 1.0 - 1e-6; // Right boundary (j=nV-1)

        for (let e = 0; e < uniqueU.length - 1; e++) {
            const uMin = uniqueU[e], uMax = uniqueU[e+1];
            if (uMax - uMin < 1e-9) continue;
            for (let qu = 0; qu < gPoints.length; qu++) {
                const u = ((uMax - uMin) * gPoints[qu] + (uMax + uMin)) / 2;
                const gw = gWeights[qu] * (uMax - uMin) / 2;
                const deriv = this.engine.getSurfaceDerivatives(patch, u, vParam);
                const tangent = deriv.dU;
                const detJ_1D = Math.sqrt(tangent.x**2 + tangent.y**2);
                if (detJ_1D < 1e-12) continue;

                const pos = deriv.pos;
                const sigma = this._getExactStress(pos.x, pos.y, 1.0, loadValue);
                const nx = tangent.y / detJ_1D;
                const ny = -tangent.x / detJ_1D;
                const tx = sigma.sxx * nx + sigma.sxy * ny;
                const ty = sigma.sxy * nx + sigma.syy * ny;

                let W_1D = 0;
                for (let a = 0; a < nU; a++) {
                    W_1D += this.engine.basis1D(a, p, U, u) * weights[a][nV-1];
                }
                if (W_1D < 1e-12) W_1D = 1e-12;

                for (let a = 0; a < nU; a++) {
                    const Na = this.engine.basis1D(a, p, U, u);
                    if (Na === 0) continue;
                    const Ra = (Na * weights[a][nV-1]) / W_1D;
                    const idx = (a * nV + (nV - 1)) * 2;
                    forces[idx]     += Ra * tx * detJ_1D * gw;
                    forces[idx + 1] += Ra * ty * detJ_1D * gw;
                }
            }
        }
        return forces;
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
            // Increase steps to 3 for better nonlinear stability
            result = this.solverFOM.solveNonlinear(this.patch, bcs, loads, { iterations: 15, steps: 3 });
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

        // NaN safety: if solver diverged, freeze at last good state
        const hasBadValues = result.u.some(v => !isFinite(v));
        if (hasBadValues) {
            console.warn('[IGA-3.2] Solver diverged — mesh frozen at last stable state');
            this._render();
            return;
        }

        this.updateMesh(result.u);
        const maxU = Math.max(...Array.from(result.u).map(Math.abs));
        if (isFinite(maxU)) console.log(`Physics Update: load=${this.loadMag}, maxDisp=${maxU.toExponential(2)}, time=${dt.toFixed(1)}ms`);
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
                // Guard against NaN from NURBS engine state
                let px = isFinite(state.position.x) ? state.position.x : 0;
                let py = isFinite(state.position.y) ? state.position.y : 0;

                let field = 0;
                if (uDisp) {
                    const disp = this.interpolateDisp(u, v, state.denominator, uDisp);
                    // NaN-safety: only apply displacement if it's finite
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
                // Final safety for position and field
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
        const res = 36;
        const { positions: posDef, values } = this.computeFieldValues(res, uDisp);
        const { positions: posUndef } = this.computeFieldValues(res, null);

        // NaN-safe min/max
        const finite = values.filter(v => isFinite(v) && v >= 0);
        const minV = finite.length ? Math.min(...finite) : 0;
        const maxV = finite.length ? Math.max(...finite) : 1;
        const range = maxV - minV || 1;

        // Compute Main Mesh Colors
        const colorFn = this.viewMode === 'stress' ? coolwarm : jet;
        const colors = [];
        values.forEach(v => {
            if (this.showMap && uDisp) {
                const param = isFinite(v) ? (v - minV) / range : 0;
                const [r, g, b] = colorFn(Math.max(0, Math.min(1, param)));
                colors.push(r, g, b);
            } else {
                // Flat Light Blue for deformed state when map is OFF
                colors.push(0.23, 0.51, 0.96);
            }
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

        // ── Dual Mesh Initialization/Update ──
        if (!this.deformedMesh) {
            const geoMain = new THREE.BufferGeometry();
            const geoGhost = new THREE.BufferGeometry();
            const idx = [];
            for (let i = 0; i < res; i++) for (let j = 0; j < res; j++) {
                const a = i * (res + 1) + j, b = (i + 1) * (res + 1) + j;
                const c = (i + 1) * (res + 1) + (j + 1), d = i * (res + 1) + (j + 1);
                idx.push(a, b, d, b, c, d);
            }

            // Deformed Mesh
            geoMain.setIndex(idx);
            geoMain.setAttribute('position', new THREE.Float32BufferAttribute(posDef, 3));
            geoMain.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
            // Use MeshStandardMaterial or MeshBasicMaterial for better colors
            this.deformedMesh = new THREE.Mesh(geoMain, new THREE.MeshStandardMaterial({
                vertexColors: true, side: THREE.DoubleSide, roughness: 0.5, metalness: 0.1
            }));
            this.deformedMesh.geometry.computeVertexNormals();

            // Ghost Mesh (Undeformed)
            geoGhost.setIndex(idx);
            geoGhost.setAttribute('position', new THREE.Float32BufferAttribute(posUndef, 3));
            this.ghostMesh = new THREE.Mesh(geoGhost, new THREE.MeshPhongMaterial({
                color: 0xcbd5e1, transparent: true, opacity: 0.35, side: THREE.DoubleSide
            }));
            this.ghostMesh.position.z = -0.05; // Move behind deformed shape
            this.ghostMesh.renderOrder = -1;   // Ensure it's rendered behind

            this.wireMesh = new THREE.LineSegments(
                new THREE.WireframeGeometry(geoMain),
                new THREE.LineBasicMaterial({ color: 0x000000, transparent: true, opacity: 0.06 })
            );

            this.scene.add(this.ghostMesh);
            this.scene.add(this.deformedMesh);
            this.scene.add(this.wireMesh);
        } else {
            // Update Deformed
            const posAttr = this.deformedMesh.geometry.getAttribute('position');
            const colAttr = this.deformedMesh.geometry.getAttribute('color');
            for (let i = 0; i < posDef.length; i++) posAttr.array[i] = posDef[i];
            for (let i = 0; i < colors.length; i++) colAttr.array[i] = colors[i];
            posAttr.needsUpdate = true;
            colAttr.needsUpdate = true;
            this.deformedMesh.geometry.computeVertexNormals();
            this.deformedMesh.geometry.computeBoundingSphere();

            // Wireframe follows deformed
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
            // More steps during training for better snapshots
            const res = this.solverFOM.solveNonlinear(this.patch, bcs, this.getLoads(f), { steps: 4, iterations: 10 });
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

        // Map Toggle
        document.getElementById('btn-map-on').onclick = () => {
            this.showMap = true;
            this._setActive('btn-map-on', 'btn-map-off');
            document.getElementById('viz-selectors').classList.replace('hidden-fade', 'visible-fade');
            document.getElementById('colorbar-panel').classList.replace('hidden-fade', 'visible-fade');
            // Always update mesh on toggle
            this.updateMesh(this.lastResult ? this.lastResult.u : null);
            this._render();
        };
        document.getElementById('btn-map-off').onclick = () => {
            this.showMap = false;
            this._setActive('btn-map-off', 'btn-map-on');
            document.getElementById('viz-selectors').classList.replace('visible-fade', 'hidden-fade');
            document.getElementById('colorbar-panel').classList.replace('visible-fade', 'hidden-fade');
            // Always update mesh on toggle
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

new ROMApp32();
