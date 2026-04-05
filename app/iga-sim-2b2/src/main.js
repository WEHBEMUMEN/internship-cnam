/**
 * Phase 2B.2 | Nonlinear Newton-Raphson Lab
 * 2D IGA Large Deflection Benchmark
 */

const engine = new NURBS2D();
const solver = new IGA2DSolver(engine);

let patch = NURBSPresets.generateSheet(); // 10x10 Quadratic Plate
let analysisData = {
    u: null,
    residuals: [],
    mode: 'nonlinear',
    load: 5000
};

let chart;

// --- Three.js ---
const scene = new THREE.Scene();
scene.background = new THREE.Color(0x020617);
const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
camera.position.set(15, 12, 18);

const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setSize(window.innerWidth, window.innerHeight);
document.getElementById('canvas-container').appendChild(renderer.domElement);

const controls = new THREE.OrbitControls(camera, renderer.domElement);
scene.add(new THREE.AmbientLight(0xffffff, 0.4));
const dLight = new THREE.DirectionalLight(0xffffff, 0.8);
dLight.position.set(10, 20, 10);
scene.add(dLight);

scene.add(new THREE.GridHelper(20, 10, 0x334155, 0x1e293b));

let surfaceMesh, wireframeOverlay;

function init() {
    initChart();
    setupUI();
    updateSurface();
    animate();
}

function initChart() {
    const ctx = document.getElementById('residual-chart').getContext('2d');
    chart = new Chart(ctx, {
        type: 'line',
        data: {
            labels: [],
            datasets: [{
                label: 'Residual Norm ||R||',
                data: [],
                borderColor: '#f43f5e',
                borderWidth: 2,
                tension: 0.2,
                fill: false,
                pointRadius: 4,
                pointBackgroundColor: '#fb7185'
            }]
        },
        options: {
            scales: {
                y: { type: 'logarithmic', grid: { color: 'rgba(255,255,255,0.05)' }, border: { dash: [4, 4] } },
                x: { grid: { display: false } }
            },
            plugins: { legend: { display: false } },
            maintainAspectRatio: false
        }
    });
}

function updateSurface() {
    if (surfaceMesh) scene.remove(surfaceMesh);
    if (wireframeOverlay) scene.remove(wireframeOverlay);

    const res = 40;
    const geometry = new THREE.BufferGeometry();
    const positions = [];
    const colors = [];
    
    // We render in REAL scale for Nonlinear, 
    // unless linear model is selected (then we might see distortion)
    const visualScale = 1.0; 

    const maxDisp = analysisData.u ? calculateMaxDisp(analysisData.u) : 1;

    for (let i = 0; i <= res; i++) {
        const u = i / res;
        for (let j = 0; j <= res; j++) {
            const v = j / res;
            const state = engine.getSurfaceState(patch, u, v);
            let p = { ...state.position };

            if (analysisData.u) {
                const interp = interpolateDisplacement(u, v, state.denominator);
                p.x += interp.x * visualScale;
                p.y += interp.y * visualScale;
                
                const mag = Math.sqrt(interp.x**2 + interp.y**2);
                const t = Math.min(mag / maxDisp, 1.0);
                const color = getHeatmapColor(t);
                colors.push(color.r, color.g, color.b);
            } else {
                colors.push(0.1, 0.2, 0.4);
            }
            positions.push(p.x, p.y, p.z);
        }
    }

    const indices = [];
    for (let i = 0; i < res; i++) {
        for (let j = 0; j < res; j++) {
            const a = i * (res+1) + j;
            const b = (i+1) * (res+1) + j;
            const c = (i+1) * (res+1) + (j+1);
            const d = i * (res+1) + (j+1);
            indices.push(a, b, d, b, c, d);
        }
    }

    geometry.setIndex(indices);
    geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
    geometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
    geometry.computeVertexNormals();

    const material = new THREE.MeshStandardMaterial({ vertexColors: true, side: THREE.DoubleSide });
    surfaceMesh = new THREE.Mesh(geometry, material);
    scene.add(surfaceMesh);
}

function interpolateDisplacement(u, v, denom) {
    const nU = patch.controlPoints.length;
    const nV = patch.controlPoints[0].length;
    let dx = 0, dy = 0;
    for (let i = 0; i < nU; i++) {
        const Ni = engine.basis1D(i, patch.p, patch.U, u);
        if (Ni === 0) continue;
        for (let j = 0; j < nV; j++) {
            const Mj = engine.basis1D(j, patch.q, patch.V, v);
            const R = (Ni * Mj * patch.weights[i][j]) / denom;
            const idx = (i * nV + j) * 2;
            dx += R * analysisData.u[idx];
            dy += R * analysisData.u[idx + 1];
        }
    }
    return { x: dx, y: dy };
}

function calculateMaxDisp(u) {
    let max = 0;
    for (let i = 0; i < u.length; i += 2) {
        const m = Math.sqrt(u[i]**2 + u[i+1]**2);
        if (m > max) max = m;
    }
    return max || 1;
}

function getHeatmapColor(t) {
    // Rose/Salmon heatmap
    return { r: 0.1 + t * 0.9, g: 0.1 + t * 0.2, b: 0.4 - t * 0.3 };
}

async function runAnalysis() {
    const btn = document.getElementById('solve-btn');
    btn.disabled = true;
    btn.innerHTML = '<i class="fa-solid fa-spinner fa-spin"></i> Converging...';
    
    document.getElementById('conv-stat').textContent = "Status: Executing Assembly...";

    // Setup BCs & Loads
    const nU = patch.controlPoints.length;
    const nV = patch.controlPoints[0].length;
    const bcs = [];
    for (let j = 0; j < nV; j++) bcs.push({ i: 0, j: j, axis: 'both', value: 0 }); // Cantilever
    
    const loadVal = parseFloat(document.getElementById('load-slider').value);
    const loads = [{ i: nU - 1, j: Math.floor(nV / 2), fx: 0, fy: -loadVal }];

    setTimeout(() => {
        const t0 = performance.now();
        if (analysisData.mode === 'linear') {
            analysisData.u = solver.solve(patch, bcs, loads);
            analysisData.residuals = [1.0, 1e-15];
        } else {
            const result = solver.solveNonlinear(patch, bcs, loads, { steps: 2, iterations: 10 });
            analysisData.u = result.u;
            analysisData.residuals = result.residualHistory;
        }
        const t1 = performance.now();

        updateChart();
        updateSurface();
        
        btn.disabled = false;
        btn.innerHTML = '<i class="fa-solid fa-bolt"></i> Run Analysis';
        document.getElementById('conv-stat').textContent = `Success: Solved in ${ (t1 - t0).toFixed(1) }ms`;
    }, 50);
}

function updateChart() {
    chart.data.labels = analysisData.residuals.map((_, i) => i + 1);
    chart.data.datasets[0].data = analysisData.residuals;
    chart.update();
}

function setupUI() {
    document.getElementById('solve-btn').onclick = runAnalysis;
    
    document.getElementById('mode-linear').onclick = () => {
        analysisData.mode = 'linear';
        document.getElementById('mode-linear').className = "py-2 px-3 rounded-lg border border-blue-500 bg-blue-500/10 text-blue-400 text-xs font-bold transition-colors";
        document.getElementById('mode-nonlinear').className = "py-2 px-3 rounded-lg border border-slate-700 text-xs font-bold hover:bg-slate-800 transition-colors";
    };
    
    document.getElementById('mode-nonlinear').onclick = () => {
        analysisData.mode = 'nonlinear';
        document.getElementById('mode-nonlinear').className = "py-2 px-3 rounded-lg border border-rose-500 bg-rose-500/10 text-rose-400 text-xs font-bold transition-colors";
        document.getElementById('mode-linear').className = "py-2 px-3 rounded-lg border border-slate-700 text-xs font-bold hover:bg-slate-800 transition-colors";
    };

    document.getElementById('load-slider').oninput = (e) => {
        document.getElementById('load-val').textContent = `${e.target.value} N`;
    };
}

function animate() {
    requestAnimationFrame(animate);
    controls.update();
    renderer.render(scene, camera);
}

window.addEventListener('resize', () => {
    renderer.setSize(window.innerWidth, window.innerHeight);
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
});

init();
