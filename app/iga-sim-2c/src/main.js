/**
 * Phase 2C | Snapshot Database Builder
 * Automates the generation of training data for ROM.
 */

const engine = new NURBS2D();
const solver = new IGA2DSolver(engine);

let patch = NURBSPresets.generateSheet();
// Ensure exact 2D Flat Plate
for(let i=0; i<patch.controlPoints.length; i++) {
    for(let j=0; j<patch.controlPoints[0].length; j++) {
        patch.controlPoints[i][j].z = 0;
    }
}

let chart;
let snapshots = [];
let parameters = [];
let forceArrow;
let boundaryVisuals = [];
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

    // Add explicit wireframe to ensure visibility
    const wireframeOverlay = new THREE.LineSegments(
        new THREE.WireframeGeometry(geometry),
        new THREE.LineBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0.3 })
    );
    scene.add(wireframeOverlay);
    
    refreshVisuals(); // Initialize indicators
}

function refreshVisuals() {
    if (!patch || !patch.controlPoints) return;
    const boundEdge = document.getElementById('bound-edge').value;
    const forceU = parseInt(document.getElementById('force-u').value) || 0;
    const forceV = parseInt(document.getElementById('force-v').value) || 0;
    const forceDir = document.getElementById('force-dir').value;
    
    const nU = patch.controlPoints.length;
    const nV = patch.controlPoints[0].length;
    
    boundaryVisuals.forEach(mesh => scene.remove(mesh));
    boundaryVisuals = [];
    
    let bcs = [];
    if (boundEdge === 'left') {
        for (let j = 0; j < nV; j++) bcs.push({ i: 0, j: j });
    } else if (boundEdge === 'right') {
        for (let j = 0; j < nV; j++) bcs.push({ i: nU - 1, j: j });
    } else if (boundEdge === 'bottom') {
        for (let i = 0; i < nU; i++) bcs.push({ i: i, j: 0 });
    } else if (boundEdge === 'top') {
        for (let i = 0; i < nU; i++) bcs.push({ i: i, j: nV - 1 });
    }

    const geo = new THREE.BoxGeometry(0.8, 0.8, 0.8);
    const mat = new THREE.MeshBasicMaterial({ color: 0x0ea5e9 }); 
    
    bcs.forEach(bc => {
        const cp = patch.controlPoints[bc.i][bc.j];
        const mesh = new THREE.Mesh(geo, mat);
        mesh.position.set(cp.x, cp.y, cp.z);
        scene.add(mesh);
        boundaryVisuals.push(mesh);
    });
    
    const isHoriz = (forceDir === 'horizontal');
    const dirVector = isHoriz ? new THREE.Vector3(1, 0, 0) : new THREE.Vector3(0, -1, 0);

    if (!forceArrow) {
        forceArrow = new THREE.ArrowHelper(dirVector, new THREE.Vector3(0,0,0), 5, 0x0ea5e9, 1.5, 1);
        scene.add(forceArrow);
    } else {
        forceArrow.setDirection(dirVector);
    }
    
    const cl_U = Math.min(Math.max(forceU, 0), nU - 1);
    const cl_V = Math.min(Math.max(forceV, 0), nV - 1);
    const cp = patch.controlPoints[cl_U][cl_V];
    
    const length = 4;
    forceArrow.setLength(length, 1.0, 0.8);
    
    if (isHoriz) {
        forceArrow.position.set(cp.x - length, cp.y, cp.z);
    } else {
        forceArrow.position.set(cp.x, cp.y + length, cp.z);
    }
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

    log(`LHS Map built. Starting solver cascade...`);

    const nU = patch.controlPoints.length, nV = patch.controlPoints[0].length;
    const bcs = [];
    const boundEdge = document.getElementById('bound-edge').value;
    if (boundEdge === 'left') {
        for (let j = 0; j < nV; j++) bcs.push({ i: 0, j: j, axis: 'both', value: 0 });
    } else if (boundEdge === 'right') {
        for (let j = 0; j < nV; j++) bcs.push({ i: nU - 1, j: j, axis: 'both', value: 0 });
    } else if (boundEdge === 'bottom') {
        for (let i = 0; i < nU; i++) bcs.push({ i: i, j: 0, axis: 'both', value: 0 });
    } else if (boundEdge === 'top') {
        for (let i = 0; i < nU; i++) bcs.push({ i: i, j: nV - 1, axis: 'both', value: 0 });
    }

    const forceU = Math.min(Math.max(parseInt(document.getElementById('force-u').value) || 0, 0), nU - 1);
    const forceV = Math.min(Math.max(parseInt(document.getElementById('force-v').value) || 0, 0), nV - 1);
    const forceDir = document.getElementById('force-dir').value;

    for (let i = 0; i < nSamples; i++) {
        const mu = samples[i];
        solver.E = mu.E * 1000; // MPa
        
        let fx = 0, fy = 0;
        if (forceDir === 'horizontal') fx = mu.Load;
        else fy = -mu.Load;
        
        const loads = [{ i: forceU, j: forceV, fx: fx, fy: fy }];
        
        // Asynchronous batch run to allow UI updates
        const result = await runSolve(patch, bcs, loads);
        
        let hasNaN = false;
        if(result && result.u) {
            for(let k=0; k<result.u.length; k++) if(isNaN(result.u[k])) hasNaN = true;
        }

        if(hasNaN) {
            log(`[ERROR] Sample ${i+1} diverged! Skipping...`);
            continue;
        }

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
    document.getElementById('start-btn').onclick = startBatch;
    document.getElementById('download-btn').onclick = downloadData;
    
    document.getElementById('n-slider').oninput = (e) => {
        document.getElementById('n-val').textContent = e.target.value;
    };

    const uiIds = ['bound-edge', 'force-u', 'force-v', 'force-dir'];
    uiIds.forEach(id => {
        document.getElementById(id).addEventListener('change', refreshVisuals);
        document.getElementById(id).addEventListener('input', refreshVisuals);
    });

    const raycaster = new THREE.Raycaster();
    const mouse = new THREE.Vector2();
    renderer.domElement.addEventListener('pointerdown', (e) => {
        if (isGenerating) return; // Lock UI during batch solving
        const rect = renderer.domElement.getBoundingClientRect();
        mouse.x = ((e.clientX - rect.left) / rect.width) * 2 - 1;
        mouse.y = -((e.clientY - rect.top) / rect.height) * 2 + 1;
        
        raycaster.setFromCamera(mouse, camera);
        if (surfaceMesh) {
            const intersects = raycaster.intersectObject(surfaceMesh);
            if (intersects.length > 0) {
                const a = intersects[0].face.a;
                const res = 40;
                const i_grid = Math.floor(a / (res+1));
                const j_grid = a % (res+1);
                const nU = patch.controlPoints.length, nV = patch.controlPoints[0].length;
                const bestU = Math.min(Math.max(Math.round((i_grid / res) * (nU - 1)), 0), nU - 1);
                const bestV = Math.min(Math.max(Math.round((j_grid / res) * (nV - 1)), 0), nV - 1);
                
                document.getElementById('force-u').value = bestU;
                document.getElementById('force-v').value = bestV;
                refreshVisuals();
            }
        }
    });
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
