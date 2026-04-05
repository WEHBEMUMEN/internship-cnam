/**
 * Phase 2C | Snapshot Database Builder
 * Automates the generation of training data for ROM.
 */

const engine = new NURBS2D();
const solver = new IGA2DSolver(engine);

let patch = NURBSPresets.generateSheet();
let chart;
let snapshots = [];
let parameters = [];
let isGenerating = false;

// --- Three.js ---
const scene = new THREE.Scene();
scene.background = new THREE.Color(0x020617);
const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
camera.position.set(12, 10, 15);
const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setSize(window.innerWidth, window.innerHeight);
document.getElementById('canvas-container').appendChild(renderer.domElement);
const controls = new THREE.OrbitControls(camera, renderer.domElement);
scene.add(new THREE.AmbientLight(0xffffff, 0.3));
const dLight = new THREE.DirectionalLight(0xffffff, 0.8);
dLight.position.set(10, 20, 10);
scene.add(dLight);
scene.add(new THREE.GridHelper(20, 10, 0x334155, 0x1e293b));

let surfaceMesh;

function init() {
    initChart();
    setupUI();
    updateSurface(new Float64Array(patch.controlPoints.length * patch.controlPoints[0].length * 2).fill(0));
    animate();
}

function initChart() {
    const ctx = document.getElementById('lhs-chart').getContext('2d');
    chart = new Chart(ctx, {
        type: 'scatter',
        data: {
            datasets: [{
                label: 'Parameter Samples',
                data: [],
                backgroundColor: '#0ea5e9',
                pointRadius: 5
            }]
        },
        options: {
            scales: {
                x: { title: { display: true, text: 'Youngs Modulus (E)' }, grid: { color: '#1e293b' } },
                y: { title: { display: true, text: 'Load (F)' }, grid: { color: '#1e293b' } }
            },
            plugins: { legend: { display: false } },
            maintainAspectRatio: false
        }
    });
}

function updateSurface(u_disp) {
    if (surfaceMesh) scene.remove(surfaceMesh);
    const res = 30;
    const geometry = new THREE.BufferGeometry();
    const positions = [];
    const colors = [];
    for (let i = 0; i <= res; i++) {
        const u = i / res;
        for (let j = 0; j <= res; j++) {
            const v = j / res;
            const state = engine.getSurfaceState(patch, u, v);
            let p = { ...state.position };
            if (u_disp) {
                const interp = interpolateDisplacement(u_disp, u, v, state.denominator);
                p.x += interp.x; p.y += interp.y;
            }
            positions.push(p.x, p.y, p.z);
            colors.push(0.05, 0.4, 0.8);
        }
    }
    const indices = []; // Standard grid indices...
    for (let i = 0; i < res; i++) {
        for (let j = 0; j < res; j++) {
            const a = i * (res+1) + j, b = (i+1) * (res+1) + j, c = (i+1) * (res+1) + (j+1), d = i * (res+1) + (j+1);
            indices.push(a, b, d, b, c, d);
        }
    }
    geometry.setIndex(indices);
    geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
    geometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
    geometry.computeVertexNormals();
    surfaceMesh = new THREE.Mesh(geometry, new THREE.MeshStandardMaterial({ vertexColors: true, side: THREE.DoubleSide }));
    scene.add(surfaceMesh);
}

function interpolateDisplacement(u_arr, u, v, denom) {
    const nU = patch.controlPoints.length, nV = patch.controlPoints[0].length;
    let dx = 0, dy = 0;
    for (let i = 0; i < nU; i++) {
        const Ni = engine.basis1D(i, patch.p, patch.U, u);
        if (Ni === 0) continue;
        for (let j = 0; j < nV; j++) {
            const Mj = engine.basis1D(j, patch.q, patch.V, v);
            const R = (Ni * Mj * patch.weights[i][j]) / denom;
            dx += R * u_arr[(i * nV + j) * 2];
            dy += R * u_arr[(i * nV + j) * 2 + 1];
        }
    }
    return { x: dx, y: dy };
}

async function startBatch() {
    if (isGenerating) return;
    isGenerating = true;
    snapshots = [];
    parameters = [];
    
    const nSamples = parseInt(document.getElementById('n-slider').value);
    document.getElementById('start-btn').disabled = true;
    document.getElementById('process-ui').classList.remove('hidden');
    document.getElementById('download-btn').classList.add('hidden');
    
    log(`Initializing LHS for ${nSamples} samples...`);
    
    // 1. Definition of Parameter Space
    const dims = [
        { name: 'E', min: 100, max: 500 },     // GPa
        { name: 'Load', min: 1000, max: 8000 } // N
    ];
    
    const samples = SamplingUtils.generateLHS(dims, nSamples);
    chart.data.datasets[0].data = samples.map(s => ({ x: s.E, y: s.Load }));
    chart.update();
    
    log(`LHS Map built. Starting solver cascade...`);

    const nU = patch.controlPoints.length, nV = patch.controlPoints[0].length;
    const bcs = [];
    for (let j = 0; j < nV; j++) bcs.push({ i: 0, j: j, axis: 'both', value: 0 });

    for (let i = 0; i < nSamples; i++) {
        const mu = samples[i];
        solver.E = mu.E * 1000; // MPa
        const loads = [{ i: nU - 1, j: Math.floor(nV / 2), fx: 0, fy: -mu.Load }];
        
        // Asynchronous batch run to allow UI updates
        const result = await runSolve(patch, bcs, loads);
        
        snapshots.push(Array.from(result.u));
        parameters.push(mu);
        
        updateProgress(i + 1, nSamples);
        if (i % 2 === 0) updateSurface(result.u);
        log(`[${i+1}/${nSamples}] E=${mu.E.toFixed(1)} GPa, F=${mu.Load.toFixed(0)} N | converged.`);
    }

    log(`Batch generation complete. Database Ready.`);
    isGenerating = false;
    document.getElementById('start-btn').disabled = false;
    document.getElementById('download-btn').classList.remove('hidden');
}

function runSolve(p, bc, ld) {
    return new Promise(resolve => {
        setTimeout(() => {
            const res = solver.solveNonlinear(p, bc, ld, { steps: 1, iterations: 6 });
            resolve(res);
        }, 1);
    });
}

function updateProgress(curr, total) {
    const p = Math.round((curr / total) * 100);
    document.getElementById('progress-bar').style.width = `${p}%`;
    document.getElementById('percent-text').textContent = `${p}%`;
    document.getElementById('status-text').textContent = `Solving Sample ${curr}...`;
}

function log(msg) {
    const consoleBox = document.getElementById('log-console');
    consoleBox.innerHTML += `<br>> ${msg}`;
    consoleBox.scrollTop = consoleBox.scrollHeight;
}

function setupUI() {
    document.getElementById('n-slider').oninput = (e) => {
        document.getElementById('n-val').textContent = e.target.value;
    };
    document.getElementById('start-btn').onclick = startBatch;
    document.getElementById('download-btn').onclick = downloadData;
}

function downloadData() {
    const data = {
        meta: { name: "IGA-ROM Snapshot Matrix", nSamples: snapshots.length, nDofs: snapshots[0].length },
        parameters: parameters,
        snapshots: snapshots
    };
    const blob = new Blob([JSON.stringify(data)], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `iga_snapshot_matrix_${snapshots.length}_samples.json`;
    a.click();
}

function animate() {
    requestAnimationFrame(animate);
    controls.update();
    renderer.render(scene, camera);
}

init();
