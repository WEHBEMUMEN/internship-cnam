/**
 * Phase 2C.1 | Validated Snapshot Generator
 * LITERAL 1:1 Port from Phase 2B.2 (Gold Standard)
 */

const engine = new NURBS2D();
const solver = new IGA2DSolver(engine);

// Benchmark Constants
const R = 1.0;
const L = 4.0;
let basePatch = NURBSPresets.generatePlateWithHole(R, L);
let patch = { ...basePatch };

let analysisData = { u: null, load: 100 };
let targetState = { load: 100, E: 100000, nu: 0.3, h: 0, p: 2, k: 0 };
let activeState = { load: -1, E: -1, nu: -1, h: -1, p: -1 };
let visibilityState = { cp: false, bc: true, force: true };
let viewMode = 'displacement';
let activeSnapshotType = 'displacement';
let activity = { isGenerating: false };

// Snapshot Storage
let snapshots = [];
let parameters = [];
let l2_errors = [];

// Solver Cache (from 2B.2 logic)
let solverCache = {
    meshKey: "",
    referenceData: null,
    useScalarScaling: true 
};

// --- Three.js Setup (Literal from 2B.2) ---
const scene = new THREE.Scene();
scene.background = new THREE.Color(0x020617);
const camera = new THREE.PerspectiveCamera(45, window.innerWidth / window.innerHeight, 0.1, 1000);
camera.position.set(2.5, 2.5, 12); 
camera.lookAt(2.5, 2.5, 0);
const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setSize(window.innerWidth, window.innerHeight);
document.getElementById('canvas-container').appendChild(renderer.domElement);
const controls = new THREE.OrbitControls(camera, renderer.domElement);
controls.enableRotate = false;
controls.target.set(2.5, 2.5, 0);
scene.add(new THREE.AmbientLight(0xffffff, 1.0));

let surfaceMesh, wireframeOverlay;
let boundaryVisuals = [], forceVisuals = [], controlPointVisuals = [];
let chart;

function init() {
    initChart();
    setupUI();
    applyRefinements(); // Triggers initial render
    animate();
}

function initChart() {
    const ctx = document.getElementById('lhs-chart').getContext('2d');
    chart = new Chart(ctx, {
        type: 'scatter',
        data: { datasets: [{ label: 'Samples', data: [], backgroundColor: '#6366f1', pointRadius: 4 }] },
        options: {
            scales: {
                x: { title: { display: true, text: 'E (GPa)', color: '#64748b', font: { size: 10 } }, grid: { color: '#1e293b' } },
                y: { title: { display: true, text: 'Load (N)', color: '#64748b', font: { size: 10 } }, grid: { color: '#1e293b' } }
            },
            plugins: { legend: { display: false } },
            maintainAspectRatio: false
        }
    });
}

function updateStatusLabels() {
    const nU = patch.controlPoints.length, nV = patch.controlPoints[0].length;
    if (document.getElementById('status-degree')) document.getElementById('status-degree').textContent = `p=${patch.p}, q=${patch.q}`;
    if (document.getElementById('status-spans')) document.getElementById('status-spans').textContent = `${patch.U.length - 2 * patch.p - 1} x ${patch.V.length - 2 * patch.q - 1}`;
    if (document.getElementById('status-points')) document.getElementById('status-points').textContent = `${nU} x ${nV}`;
}

function applyRefinements() {
    patch = JSON.parse(JSON.stringify(basePatch));
    RefineUtils.apply(engine, patch, { h: targetState.h, p: targetState.p });
    // Handle K-Refinement (accumulated)
    for(let i=0; i<targetState.k; i++) {
        patch = engine.elevateDegree(patch);
        patch = engine.subdivideGlobal(patch);
    }
    updateStatusLabels();
    updateSurface();
    updateBoundaryVisuals();
    updateControlPoints();
    updateForceArrows();
}

function updateSurface(u_disp = null) {
    if (surfaceMesh) scene.remove(surfaceMesh);
    if (wireframeOverlay) scene.remove(wireframeOverlay);

    const res = 50;
    const geometry = new THREE.BufferGeometry();
    const positions = [], colors = [];
    const u_active = u_disp || (solverCache.referenceData ? solverCache.referenceData.u_ref : null);
    const scale = (u_disp) ? 1.0 : (solverCache.referenceData ? (targetState.load / solverCache.referenceData.load_ref) : 1);

    let maxVal = 1;
    if (u_active) maxVal = (viewMode === 'displacement') ? calculateMaxDisp(u_active) * scale : calculateMaxStress(u_active) * scale;
    if (document.getElementById('legend-val-max')) document.getElementById('legend-val-max').innerText = u_active ? `${maxVal.toFixed(viewMode === 'displacement' ? 4 : 2)}` : "---";
    
    const titleEl = document.getElementById('legend-mode-title');
    const subEl = document.getElementById('legend-mode-sub');
    if (titleEl) titleEl.innerText = (viewMode === 'displacement') ? "Displacement" : "Von-Mises Stress";
    if (subEl) subEl.innerText = (viewMode === 'displacement') ? "Magnitude (mm)" : "Stress (MPa)";

    for (let i = 0; i <= res; i++) {
        const u = i / res;
        for (let j = 0; j <= res; j++) {
            const v = j / res;
            const state = engine.getSurfaceState(patch, u, v);
            let p = { ...state.position };
            let val = 0;

            if (u_active && state.denominator > 0) {
                const interp = interpolateDisplacement(u, v, state.denominator, u_active);
                p.x += interp.x * scale; p.y += interp.y * scale;
                val = (viewMode === 'displacement') ? Math.abs(interp.x * scale) : solver.getNumericalStress(patch, u_active, u, v, targetState.E, targetState.nu).vonMises * scale;
                const t = Math.min(Math.max(val / (maxVal || 1e-6), 0), 1.0);
                const c = getHeatmapColor(t);
                colors.push(c.r, c.g, c.b);
            } else {
                colors.push(0.1, 0.2, 0.4); // Dark blue initial
            }
            positions.push(p.x, p.y, p.z);
        }
    }

    const indices = [];
    for (let i = 0; i < res; i++) for (let j = 0; j < res; j++) {
        const a = i * (res+1) + j, b = (i+1) * (res+1) + j, c = (i+1) * (res+1) + (j+1), d = i * (res+1) + (j+1);
        indices.push(a, b, d, b, c, d);
    }

    geometry.setIndex(indices);
    geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
    geometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
    surfaceMesh = new THREE.Mesh(geometry, new THREE.MeshBasicMaterial({ vertexColors: true, side: THREE.DoubleSide }));
    scene.add(surfaceMesh);
    wireframeOverlay = new THREE.LineSegments(new THREE.WireframeGeometry(geometry), new THREE.LineBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0.15 }));
    scene.add(wireframeOverlay);
}

function updateBoundaryVisuals() {
    boundaryVisuals.forEach(obj => scene.remove(obj));
    boundaryVisuals = [];
    const nU = patch.controlPoints.length, nV = patch.controlPoints[0].length; 
    for(let j=0; j<nV; j++) {
        const cp = patch.controlPoints[0][j]; const m = createBCMarker('y'); m.position.set(cp.x, cp.y, cp.z);
        scene.add(m); boundaryVisuals.push(m);
    }
    for(let j=0; j<nV; j++) {
        const cp = patch.controlPoints[nU-1][j]; const m = createBCMarker('x'); m.position.set(cp.x, cp.y, cp.z);
        scene.add(m); boundaryVisuals.push(m);
    }
}

function updateForceArrows() {
    forceVisuals.forEach(obj => scene.remove(obj));
    forceVisuals = [];
    if (targetState.load === 0) return;
    const posStart = engine.evaluateSurface(patch, 0.0, 1.0), posEnd = engine.evaluateSurface(patch, 0.5, 1.0);
    const yMin = Math.min(posStart.y, posEnd.y), yMax = Math.max(posStart.y, posEnd.y);
    for (let i = 0; i <= 8; i++) {
        const yTarget = yMin + (i / 8) * (yMax - yMin);
        let bestU = 0, minDist = Infinity;
        for (let j = 0; j <= 200; j++) {
            const uTry = j / 200, pos = engine.evaluateSurface(patch, uTry, 1.0), dist = Math.abs(pos.y - yTarget);
            if (Math.abs(pos.x - L) < 0.2 && dist < minDist) { minDist = dist; bestU = uTry; }
        }
        const state = engine.getSurfaceState(patch, bestU, 1.0);
        const arrow = new THREE.ArrowHelper(new THREE.Vector3(1,0,0), new THREE.Vector3(state.position.x, state.position.y, state.position.z), targetState.load/5000 + 0.3, 0xef4444);
        scene.add(arrow); forceVisuals.push(arrow);
    }
}

function updateControlPoints() {
    controlPointVisuals.forEach(obj => scene.remove(obj));
    controlPointVisuals = [];
    if (!visibilityState.cp) return;
    const geo = new THREE.SphereGeometry(0.06), mat = new THREE.MeshBasicMaterial({ color: 0x10b981 });
    patch.controlPoints.forEach(row => row.forEach(cp => {
        const s = new THREE.Mesh(geo, mat); s.position.set(cp.x, cp.y, cp.z);
        scene.add(s); controlPointVisuals.push(s);
    }));
}

function getHeatmapColor(t) {
    const stops = [{ t: 0.0, r: 0.23, g: 0.51, b: 0.96 }, { t: 0.2, r: 0.02, g: 0.71, b: 0.83 }, { t: 0.4, r: 0.06, g: 0.73, b: 0.51 }, { t: 0.6, r: 0.98, g: 0.80, b: 0.08 }, { t: 0.8, r: 0.98, g: 0.45, b: 0.09 }, { t: 1.0, r: 0.94, g: 0.27, b: 0.27 }];
    for (let i = 0; i < stops.length - 1; i++) {
        const s1 = stops[i], s2 = stops[i+1];
        if (t >= s1.t && t <= s2.t) {
            const factor = (t - s1.t) / (s2.t - s1.t);
            return { r: s1.r + factor * (s2.r - s1.r), g: s1.g + factor * (s2.g - s1.g), b: s1.b + factor * (s2.b - s1.b) };
        }
    }
    return stops[stops.length-1];
}

function createBCMarker(type = 'x') {
    const group = new THREE.Group();
    const triGeo = new THREE.ConeGeometry(0.12, 0.2, 3), triMat = new THREE.MeshBasicMaterial({ color: 0x334155 });
    const triangle = new THREE.Mesh(triGeo, triMat), rollGeo = new THREE.CircleGeometry(0.04, 16), rollMat = new THREE.MeshBasicMaterial({ color: 0x64748b });
    const roller = new THREE.Mesh(rollGeo, rollMat);
    if (type === 'y') { triangle.position.y = -0.15; roller.position.y = -0.28; }
    else { triangle.rotation.z = -Math.PI / 2; triangle.position.x = -0.15; roller.position.x = -0.28; }
    group.add(triangle); group.add(roller);
    return group;
}

function interpolateDisplacement(u, v, denom, u_disp) {
    const nU = patch.controlPoints.length, nV = patch.controlPoints[0].length;
    let dx = 0, dy = 0;
    for (let i = 0; i < nU; i++) {
        const Ni = engine.basis1D(i, patch.p, patch.U, u);
        if (Ni === 0) continue;
        for (let j = 0; j < nV; j++) {
            const Mj = engine.basis1D(j, patch.q, patch.V, v);
            const R = (Ni * Mj * patch.weights[i][j]) / denom;
            dx += R * u_disp[(i * nV + j) * 2]; dy += R * u_disp[(i * nV + j) * 2 + 1];
        }
    }
    return { x: dx, y: dy };
}

function calculateMaxDisp(u) { if(!u) return 1; let max = 0; for (let i = 0; i < u.length; i+=2) { const m = Math.abs(u[i]); if (m > max) max = m; } return max || 1; }
function calculateMaxStress(u_disp) { if(!u_disp) return 1; let max = 0; const res = 15; for (let i = 0; i <= res; i++) for (let j = 0; j <= res; j++) { const s = solver.getNumericalStress(patch, u_disp, i/res, (j/res)*0.95, targetState.E, targetState.nu); if (s.vonMises > max) max = s.vonMises; } return max || 1; }

// --- Generator & Probe Logic ---

async function solveVerified() {
    const nU = patch.controlPoints.length, nV = patch.controlPoints[0].length;
    const bcs = [];
    for(let j=0; j<nV; j++) bcs.push({ i: 0, j: j, axis: 'y', value: 0 });
    for(let j=0; j<nV; j++) bcs.push({ i: nU-1, j: j, axis: 'x', value: 0 });
    const integratedForces = solver.calculateNodalTraction(patch, targetState.load, 'right');
    const loads = [];
    for (let i = 0; i < integratedForces.length; i++) {
        if (integratedForces[i] !== 0) {
            const nodeIdx = Math.floor(i / 2), a = Math.floor(nodeIdx / nV), b = nodeIdx % nV;
            loads.push({ i: a, j: b, fx: (i%2==0 ? integratedForces[i] : 0), fy: (i%2==0 ? 0 : integratedForces[i]) });
        }
    }
    solver.E = targetState.E;
    solver.nu = targetState.nu;
    const u_sol = solver.solve(patch, bcs, loads);
    const l2u = solver.calculateRelativeL2DisplacementError(patch, u_sol, targetState.E, targetState.nu, targetState.load, R);
    const l2s = solver.calculateRelativeL2Error(patch, u_sol, targetState.E, targetState.nu, targetState.load, R);
    return { u: u_sol, l2u: l2u, l2s: l2s };
}

async function runProbe() {
    if (activity.isGenerating) return;
    const btn = document.getElementById('probe-btn');
    btn.disabled = true; btn.innerHTML = '<i class="fa-solid fa-spinner fa-spin mr-2"></i> SOLVING...';
    await new Promise(r => setTimeout(r, 10));
    const res = await solveVerified();
    solverCache.referenceData = { u_ref: res.u, load_ref: targetState.load };
    
    // Automatically switch view if stress type is selected
    if (activeSnapshotType === 'stress') viewMode = 'stress';
    updateSurface(res.u);
    
    document.getElementById('probe-results').classList.remove('hidden');
    document.getElementById('max-disp').textContent = `${calculateMaxDisp(res.u).toFixed(4)} mm`;
    document.getElementById('max-stress').textContent = `${calculateMaxStress(res.u).toFixed(2)} MPa`;
    
    const activeL2 = (viewMode === 'displacement') ? res.l2u : res.l2s;
    const l2Label = (viewMode === 'displacement') ? "L2 Error (U)" : "L2 Error (S)";
    const l2El = document.getElementById('l2-error');
    if (l2El) {
        l2El.textContent = `${(activeL2 * 100).toFixed(4)}%`;
        l2El.previousElementSibling.textContent = l2Label + ":";
    }
    btn.disabled = false; btn.innerHTML = '<i class="fa-solid fa-microscope mr-2"></i> RUN SINGLE VALIDATION';
}

async function startBatch() {
    if (activity.isGenerating) return;
    activity.isGenerating = true;
    snapshots = []; parameters = []; l2_errors = [];
    const nSamples = parseInt(document.getElementById('n-slider').value);
    document.getElementById('start-btn').disabled = true;
    document.getElementById('process-ui').classList.remove('hidden');
    document.getElementById('log-console').innerHTML = "";
    
    const dims = [ { name: 'E', min: 100, max: 500 }, { name: 'Load', min: 1000, max: 8000 } ];
    const samples = SamplingUtils.generateLHS(dims, nSamples);
    chart.data.datasets[0].data = samples.map(s => ({ x: s.E, y: s.Load }));
    chart.update();
    
    log(`Batch initialized (${nSamples} samples). Type: ${activeSnapshotType}`);

    for (let i = 0; i < nSamples; i++) {
        const mu = samples[i];
        targetState.E = mu.E * 1000; targetState.load = mu.Load;
        const res = await solveVerified();
        
        let snapshotData = [];
        if (activeSnapshotType === 'displacement') {
            snapshotData = Array.from(res.u);
        } else {
            // Stress Snapshot: Sample von Mises on a 30x30 grid
            const resGrid = 30;
            for (let gu = 0; gu <= resGrid; gu++) {
                for (let gv = 0; gv <= resGrid; gv++) {
                    const s = solver.getNumericalStress(patch, res.u, gu/resGrid, gv/resGrid, targetState.E, targetState.nu);
                    snapshotData.push(s.vonMises);
                }
            }
        }

        snapshots.push(snapshotData); 
        parameters.push({ ...mu }); 
        l2_errors.push(activeSnapshotType === 'displacement' ? res.l2u : res.l2s);
        
        updateProgress(i + 1, nSamples);
        updateSurface(res.u);
        log(`[${i+1}/${nSamples}] E=${mu.E.toFixed(0)}k F=${mu.Load.toFixed(0)} | Err: ${(res.l2*100).toFixed(3)}%`);
    }

    log(`Database generation complete.`);
    activity.isGenerating = false;
    document.getElementById('start-btn').disabled = false;
    document.getElementById('download-btn').classList.remove('hidden');
}

function updateProgress(curr, total) {
    const p = Math.round((curr / total) * 100);
    document.getElementById('progress-bar').style.width = `${p}%`;
    document.getElementById('percent-text').textContent = `${p}%`;
    document.getElementById('status-text').textContent = `Status: ${curr}/${total}`;
}

function log(msg) {
    const box = document.getElementById('log-console');
    box.innerHTML += `> ${msg}<br>`;
    box.scrollTop = box.scrollHeight;
}

function setupUI() {
    document.getElementById('start-btn').onclick = startBatch;
    document.getElementById('probe-btn').onclick = runProbe;
    document.getElementById('n-slider').oninput = (e) => document.getElementById('n-val').innerText = e.target.value;
    document.getElementById('h-slider').oninput = (e) => { targetState.h = parseInt(e.target.value); applyRefinements(); };
    document.getElementById('p-slider').oninput = (e) => { targetState.p = parseInt(e.target.value); applyRefinements(); };
    document.getElementById('apply-k').onclick = () => { targetState.k++; applyRefinements(); };
    
    document.getElementById('e-slider').oninput = (e) => { targetState.E = parseInt(e.target.value); document.getElementById('e-val').innerText = `${(targetState.E/1000).toFixed(0)}k`; };
    document.getElementById('load-slider').oninput = (e) => { targetState.load = parseInt(e.target.value); document.getElementById('load-val').innerText = `${targetState.load} N`; };
    document.getElementById('download-btn').onclick = () => {
        const data = { 
            meta: { 
                name: "IGA Snapshots", 
                h: targetState.h, 
                p: targetState.p,
                k: targetState.k,
                type: activeSnapshotType,
                grid: (activeSnapshotType === 'stress') ? 30 : null
            }, 
            parameters, 
            snapshots, 
            l2_errors 
        };
        const blob = new Blob([JSON.stringify(data)], { type: 'application/json' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a'); a.href = url; a.download = `snapshots_${activeSnapshotType}.json`; a.click();
    };
    window.setViewMode = (m) => { viewMode = m; updateSurface(); };
    window.setSnapshotMeta = (t) => { activeSnapshotType = t; };
    updateStatusLabels();
}

function animate() { requestAnimationFrame(animate); controls.update(); renderer.render(scene, camera); }
init();
