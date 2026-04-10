/**
 * Phase 2D.1 | SVD Reduction & Operator Projection
 */

const engine = new NURBS2D();
const solver = new IGA2DSolver(engine);

let snapshotData = null;
let patch = null;

// SVD Results
let svdResult = null;
let singularValues = [];
let cumulativeEnergy = [];
let PODModes = [];

// UI state
let targetState = { k: 1 };
let viewMode = 'displacement';
let currentModeIndex = 0;

// Three.js Setup
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
let chart;

function init() {
    initChart();
    setupUI();
    animate();
}

function initChart() {
    const ctx = document.getElementById('analysis-chart').getContext('2d');
    chart = new Chart(ctx, {
        type: 'line',
        data: { labels: [], datasets: [{ label: 'Cumulative Energy %', data: [], borderColor: '#38bdf8', backgroundColor: 'rgba(56, 189, 248, 0.2)', fill: true, tension: 0.2 }] },
        options: {
            scales: {
                x: { title: { display: true, text: 'Mode Count (k)', color: '#64748b', font: { size: 10 } }, grid: { color: '#1e293b' } },
                y: { title: { display: true, text: 'Energy (%)', color: '#64748b', font: { size: 10 } }, grid: { color: '#1e293b' }, min: 90, max: 100 }
            },
            plugins: { legend: { display: false } },
            maintainAspectRatio: false
        }
    });
}

function switchChart(type) {
    document.getElementById('tab-scree').className = `flex-1 py-1.5 text-[9px] font-bold rounded-md transition-all ${type === 'scree' ? 'bg-sky-600 text-white' : 'text-slate-400 hover:text-slate-200'}`;
    document.getElementById('tab-error').className = `flex-1 py-1.5 text-[9px] font-bold rounded-md transition-all ${type === 'error' ? 'bg-rose-600 text-white' : 'text-slate-400 hover:text-slate-200'}`;
    document.getElementById('tab-time').className = `flex-1 py-1.5 text-[9px] font-bold rounded-md transition-all ${type === 'time' ? 'bg-emerald-600 text-white' : 'text-slate-400 hover:text-slate-200'}`;
    
    if (!svdResult) return;

    if (type === 'scree') {
        chart.data.labels = singularValues.map((_, i) => i + 1);
        chart.data.datasets[0] = { label: 'Cumulative Energy %', data: cumulativeEnergy.map(e => e * 100), borderColor: '#38bdf8', backgroundColor: 'rgba(56, 189, 248, 0.2)', fill: true, tension: 0.2 };
        chart.options.scales.y.min = cumulativeEnergy[0] * 100 * 0.99;
        chart.options.scales.y.max = 100;
        chart.options.scales.y.title.text = 'Energy (%)';
    } else if (type === 'error') {
        // Simplified error curve: 1 - energy
        const err = cumulativeEnergy.map(e => (1 - e) * 100);
        chart.data.labels = singularValues.map((_, i) => i + 1);
        chart.data.datasets[0] = { label: 'Truncation Error %', data: err, borderColor: '#fb7185', backgroundColor: 'rgba(251, 113, 133, 0.2)', fill: true, tension: 0.2 };
        chart.options.scales.y.min = 0;
        chart.options.scales.y.max = err[0] * 1.1;
        chart.options.scales.y.title.text = 'Error (%)';
    } else {
        // Speed up illustration
        const modes = singularValues.map((_, i) => i + 1);
        const fomTime = modes.map(() => 100); // 100% time
        const romTime = modes.map(k => (k/100) * 100 + 1); // rough scaling
        chart.data.labels = modes;
        chart.data.datasets = [
            { label: 'FOM Time', data: fomTime, borderColor: '#ef4444', borderDash: [5, 5], fill: false, tension: 0 },
            { label: 'ROM Time', data: romTime, borderColor: '#10b981', fill: true, backgroundColor: 'rgba(16, 185, 129, 0.2)', tension: 0.2 }
        ];
        chart.options.scales.y.min = 0;
        chart.options.scales.y.max = 110;
        chart.options.scales.y.title.text = 'Relative Time (%)';
    }
    
    chart.update();
}

function processSnapshots(data) {
    if (!data.snapshots || data.snapshots.length === 0) return;
    document.getElementById('upload-status').innerText = "Processing SVD...";
    
    // Build S matrix (each column is a snapshot)
    const { Matrix, SVD } = mlMatrix;
    let Sarray = [];
    
    if (data.meta.type === 'displacement') {
        for (let i = 0; i < data.snapshots.length; i++) {
            Sarray.push(data.snapshots[i]);
        }
    } else {
        // mixed or displacement only?
        for (let i = 0; i < data.snapshots.length; i++) {
            let vec = data.snapshots[i].disp || data.snapshots[i];
            Sarray.push(vec);
        }
    }
    
    // S is nDofs x nSamples
    const S = new Matrix(Sarray).transpose();
    
    // Perform SVD
    const t0 = performance.now();
    const svd = new SVD(S, { computeLeftSingularVectors: true, computeRightSingularVectors: false, autoTranspose: true });
    const t1 = performance.now();
    
    svdResult = svd;
    singularValues = svd.diagonal;
    
    // Extract POD Modes (left singular vectors)
    PODModes = svd.leftSingularVectors.to2DArray();
    // Transpose back so PODModes[i] is a vector
    PODModes = new Matrix(PODModes).transpose().to2DArray();
    
    // Compute Energy
    const totalEnergy = singularValues.reduce((sum, v) => sum + v*v, 0);
    cumulativeEnergy = [];
    let currentEnergy = 0;
    for (let i = 0; i < singularValues.length; i++) {
        currentEnergy += singularValues[i]*singularValues[i];
        cumulativeEnergy.push(currentEnergy / totalEnergy);
    }
    
    rebuildMeshFromMeta(data.meta);
    document.getElementById('upload-status').innerText = `SVD complete in ${(t1-t0).toFixed(0)} ms. Modes: ${singularValues.length}`;
    
    document.getElementById('k-slider').disabled = false;
    document.getElementById('k-slider').max = singularValues.length;
    document.getElementById('visualize-mode-btn').disabled = false;
    document.getElementById('project-btn').disabled = false;
    document.getElementById('mode-info').classList.remove('hidden');
    
    updateKSlider();
    switchChart('scree');
}

function rebuildMeshFromMeta(meta) {
    if (!meta.mesh) {
        console.warn("No mesh metadata found. Using fallback.");
        patch = { ...NURBSPresets.generatePlateWithHole(meta.geometry?.R || 1.0, meta.geometry?.L || 4.0) };
        RefineUtils.apply(engine, patch, { h: meta.analysis?.h || 0, p: meta.analysis?.p || 2 });
        return;
    }
    patch = {
        U: meta.mesh.knotU,
        V: meta.mesh.knotV,
        p: meta.mesh.p,
        q: meta.mesh.q,
        weights: meta.mesh.weights,
        controlPoints: meta.mesh.controlPoints_base
    };
    updateSurface(new Float64Array(patch.controlPoints.length * patch.controlPoints[0].length * 2));
}

function updateKSlider() {
    targetState.k = parseInt(document.getElementById('k-slider').value);
    document.getElementById('k-val').innerText = targetState.k;
    const energy = cumulativeEnergy[targetState.k - 1];
    document.getElementById('energy-val').innerText = (energy * 100).toFixed(5) + "%";
    document.getElementById('error-val').innerText = ((1 - energy) * 100).toFixed(5) + "%";
}

function visualizeCurrentMode() {
    if (!PODModes || PODModes.length === 0) return;
    const modeIdx = currentModeIndex % targetState.k;
    document.getElementById('btn-mode-num').innerText = modeIdx + 1;
    
    // Scale mode for visualization (modes are normalized to 1)
    const modeVec = PODModes[modeIdx];
    const maxVal = Math.max(...modeVec.map(Math.abs));
    const scale = 0.5 / (maxVal || 1); // Target 0.5 visual size
    
    const scaledMode = modeVec.map(v => v * scale);
    updateSurface(scaledMode, true);
    currentModeIndex++;
}

async function projectOperators() {
    const btn = document.getElementById('project-btn');
    btn.disabled = true;
    btn.innerHTML = '<i class="fa-solid fa-spinner fa-spin mr-2"></i> PROJECTING...';
    
    await new Promise(r => setTimeout(r, 10)); // UI refresh
    
    const { Matrix } = mlMatrix;
    // Basis matrix slice
    const Phi = new Matrix(PODModes.slice(0, targetState.k)).transpose(); // nDofs x k
    
    // We compute K1 and K2. K1 = K evaluated with E=1, nu=0. K2 = K evaluated with E=1, nu=delta
    // Affine: K(E, nu) = E / (1-nu^2) * [ K_main + nu * K_cross ]
    // We'll extract K_main and K_cross from the solver.
    const K_main = solver.assembleStiffness(patch, 1.0, 0.0);
    const K_nu_eval = solver.assembleStiffness(patch, 1.0, 0.3); // evaluating at nu=0.3 to extract K_cross
    
    const kDim = K_main.length;
    const K1mat = new Matrix(kDim, kDim);
    const K2mat = new Matrix(kDim, kDim);
    
    for (let i = 0; i < kDim; i++) {
        for (let j = 0; j < kDim; j++) {
            K1mat.set(i, j, K_main[i][j]);
            const k2_val = (K_nu_eval[i][j] * (1 - 0.3*0.3) - K_main[i][j]) / 0.3;
            K2mat.set(i, j, k2_val);
        }
    }
    
    // Project K1_red = Phi^T * K1 * Phi
    const PhiT = Phi.transpose();
    const K1_red = PhiT.mmul(K1mat).mmul(Phi).to2DArray();
    const K2_red = PhiT.mmul(K2mat).mmul(Phi).to2DArray();
    
    // Compute Forces (Default Traction)
    const integratedForces = solver.calculateNodalTraction(patch, 100, 'right'); // 100N base force
    const F_full = new Float64Array(kDim).fill(0);
    const nV = patch.controlPoints[0].length;
    let idx = 0;
    for (let i = 0; i < integratedForces.length; i++) {
        if (integratedForces[i] !== 0) {
            const nodeIdx = Math.floor(i / 2);
            if (i%2==0) F_full[nodeIdx*2] = integratedForces[i];
            else F_full[nodeIdx*2 + 1] = integratedForces[i];
        }
    }
    const Fmat = new Matrix([Array.from(F_full)]).transpose();
    const F_red = PhiT.mmul(Fmat).to2DArray().map(r => r[0]);
    
    // Apply boundary conditions directly to reduced system via Penalty on full system before projecting?
    // Actually, penalty on K1:
    const K_penalty = solver.assembleStiffness(patch, 0, 0); // zero stiffness
    solver.applyPenaltyConstraints(K_penalty, patch);
    const KpMat = new Matrix(K_penalty);
    const Kp_red = PhiT.mmul(KpMat).mmul(Phi).to2DArray();
    for (let i=0; i<targetState.k; i++) {
        for(let j=0; j<targetState.k; j++) {
             K1_red[i][j] += Kp_red[i][j];
        }
    }

    document.getElementById('projection-status').classList.remove('hidden');
    btn.innerHTML = '<i class="fa-solid fa-check mr-2"></i> PROJECT OPERATORS';
    btn.disabled = false;
    document.getElementById('export-basis-btn').disabled = false;
    
    const vBtn = document.getElementById('validate-btn');
    vBtn.disabled = false;
    vBtn.classList.remove('hidden');
    
    // Save for export
    window.ROMData = { Phi: Phi.to2DArray(), K1_red, K2_red, Fr: F_red };
}

async function validateReducedModel() {
    const btn = document.getElementById('validate-btn');
    btn.disabled = true;
    btn.innerHTML = '<i class="fa-solid fa-spinner fa-spin mr-2"></i> VALIDATING...';
    
    await new Promise(r => setTimeout(r, 10)); // UI refresh
    
    const { Matrix } = mlMatrix;
    const Phi = new Matrix(window.ROMData.Phi);
    const PhiT = Phi.transpose();
    
    // 1. Orthogonality Check: || Phi^T * Phi - I ||
    const identityTest = PhiT.mmul(Phi);
    let orthoErr = 0;
    for(let i=0; i<identityTest.rows; i++){
        for(let j=0; j<identityTest.columns; j++){
            let target = (i===j) ? 1.0 : 0.0;
            orthoErr += Math.pow(identityTest.get(i,j) - target, 2);
        }
    }
    orthoErr = Math.sqrt(orthoErr);
    document.getElementById('ortho-val').innerText = orthoErr.toExponential(3);
    
    // 2. Snapshot Reconstruction Error
    // Pick the middle snapshot to test
    const testIdx = Math.floor(snapshotData.snapshots.length / 2);
    let u_snap = snapshotData.snapshots[testIdx];
    if (u_snap.disp) u_snap = u_snap.disp;
    
    const u_mat = new Matrix([Array.from(u_snap)]).transpose(); // Nx1
    const q = PhiT.mmul(u_mat); // kx1
    const u_tilde = Phi.mmul(q); // Nx1
    
    let l2_err_num = 0;
    let l2_err_den = 0;
    for(let i=0; i<u_snap.length; i++) {
        let diff = u_snap[i] - u_tilde.get(i,0);
        l2_err_num += diff*diff;
        l2_err_den += u_snap[i]*u_snap[i];
    }
    const rel_err = Math.sqrt(l2_err_num / (l2_err_den || 1));
    document.getElementById('rom-l2-val').innerText = (rel_err * 100).toFixed(6) + "%";
    
    document.getElementById('validation-info').classList.remove('hidden');
    btn.innerHTML = '<i class="fa-solid fa-check-double mr-2"></i> VALIDATION COMPLETE';
}

function updateSurface(u_disp = null, skipHeatmap = false) {
    if (surfaceMesh) scene.remove(surfaceMesh);
    if (wireframeOverlay) scene.remove(wireframeOverlay);

    if (!patch) return;

    const res = 50;
    const geometry = new THREE.BufferGeometry();
    const positions = [], colors = [];
    
    let maxVal = 1;

    for (let i = 0; i <= res; i++) {
        const u = i / res;
        for (let j = 0; j <= res; j++) {
            const v = j / res;
            const state = engine.getSurfaceState(patch, u, v);
            let p = { ...state.position };
            
            if (u_disp && state.denominator > 0) {
                const interp = interpolateDisplacement(u, v, state.denominator, u_disp);
                p.x += interp.x; p.y += interp.y;
                if (!skipHeatmap) {
                    const c = getHeatmapColor(Math.abs(interp.x) / maxVal);
                    colors.push(c.r, c.g, c.b);
                } else {
                    // For mode shapes, use a cool mixed color based on magnitude
                    const mag = Math.sqrt(interp.x*interp.x + interp.y*interp.y);
                    colors.push(0.5 + mag, 0.2 + mag*0.5, 0.8);
                }
            } else {
                colors.push(0.1, 0.2, 0.4);
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

function setupUI() {
    document.getElementById('file-upload').addEventListener('change', function(e) {
        const file = e.target.files[0];
        if (!file) return;
        const reader = new FileReader();
        reader.onload = function(evt) {
            snapshotData = JSON.parse(evt.target.result);
            processSnapshots(snapshotData);
        };
        reader.readAsText(file);
    });

    document.getElementById('k-slider').oninput = updateKSlider;
    document.getElementById('visualize-mode-btn').onclick = visualizeCurrentMode;
    document.getElementById('project-btn').onclick = projectOperators;
    document.getElementById('validate-btn').onclick = validateReducedModel;
    document.getElementById('export-basis-btn').onclick = () => {
        if (!window.ROMData) return;
        ROMExporter.exportSmartBasis(snapshotData.meta, window.ROMData.Phi, window.ROMData.K1_red, window.ROMData.K2_red, window.ROMData.Fr);
    };
}

function animate() { requestAnimationFrame(animate); controls.update(); renderer.render(scene, camera); }
init();
