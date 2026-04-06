/**
 * Phase 2B.2 | Nonlinear Newton-Raphson Lab
 * 2D IGA Large Deflection Benchmark
 */

const engine = new NURBS2D();
const solver = new IGA2DSolver(engine);

let patch = NURBSPresets.generateSheet(); // 10x10 Quadratic Plate
// Ensure exact 2D Flat Plate
for(let i=0; i<patch.controlPoints.length; i++) {
    for(let j=0; j<patch.controlPoints[0].length; j++) {
        patch.controlPoints[i][j].z = 0;
    }
}

let analysisData = {
    u: null,
    residuals: [],
    mode: 'nonlinear',
    load: 5000
};

let targetState = {
    mode: 'nonlinear',
    load: 5000,
    forceU: 3,
    forceV: 1,
    forceDir: 'vertical',
    boundEdge: 'left'
};
let activeState = { mode: '', load: -1, forceU: -1, forceV: -1, forceDir: '', boundEdge: '' };
let isSolving = false;
let forceArrow;
let boundaryVisuals = [];
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

    // Add explicit wireframe to ensure visibility
    wireframeOverlay = new THREE.LineSegments(
        new THREE.WireframeGeometry(geometry),
        new THREE.LineBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0.3 })
    );
    scene.add(wireframeOverlay);
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

function updateForceArrow(uIndex, vIndex, load, dirString) {
    const isHoriz = (dirString === 'horizontal');
    const dirVector = isHoriz ? new THREE.Vector3(1, 0, 0) : new THREE.Vector3(0, -1, 0);

    if (!forceArrow) {
        forceArrow = new THREE.ArrowHelper(dirVector, new THREE.Vector3(0,0,0), 5, 0xef4444, 1.5, 1);
        scene.add(forceArrow);
    } else {
        forceArrow.setDirection(dirVector);
    }
    
    const cp = patch.controlPoints[Math.min(uIndex, patch.controlPoints.length - 1)][Math.min(vIndex, patch.controlPoints[0].length - 1)];
    let px = cp.x, py = cp.y, pz = cp.z;
    if (analysisData.u) {
        const idx = (uIndex * patch.controlPoints[0].length + vIndex) * 2;
        px += analysisData.u[idx];
        py += analysisData.u[idx + 1];
    }
    
    const length = 2 + (load / 2000);
    forceArrow.setLength(length, 1.0, 0.8);
    
    if (isHoriz) {
        // Arrow points to the RIGHT (1, 0, 0). So origin is to the LEFT.
        forceArrow.position.set(px - length, py, pz);
    } else {
        // Arrow points DOWN (0, -1, 0). So origin is ABOVE.
        forceArrow.position.set(px, py + length, pz);
    }
}

function updateBoundaryVisuals(boundEdge) {
    boundaryVisuals.forEach(mesh => scene.remove(mesh));
    boundaryVisuals = [];
    
    const nU = patch.controlPoints.length;
    const nV = patch.controlPoints[0].length;
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
    const mat = new THREE.MeshBasicMaterial({ color: 0x3b82f6 }); // Blue-500
    
    bcs.forEach(bc => {
        const cp = patch.controlPoints[bc.i][bc.j];
        const mesh = new THREE.Mesh(geo, mat);
        let px = cp.x, py = cp.y, pz = cp.z;
        if (analysisData.u) {
            const idx = (bc.i * nV + bc.j) * 2;
            px += analysisData.u[idx];
            py += analysisData.u[idx + 1];
        }
        mesh.position.set(px, py, pz);
        scene.add(mesh);
        boundaryVisuals.push(mesh);
    });
}

async function solverLoop() {
    if (!isSolving) {
        if (targetState.load !== activeState.load || targetState.forceU !== activeState.forceU || targetState.forceV !== activeState.forceV || targetState.mode !== activeState.mode || targetState.forceDir !== activeState.forceDir || targetState.boundEdge !== activeState.boundEdge) {
            isSolving = true;
            const stateToSolve = { ...targetState };
            
            // Yield to browser UI thread
            await new Promise(resolve => setTimeout(resolve, 5));
            
            const t0 = performance.now();
            document.getElementById('conv-stat').textContent = "Status: Executing Assembly...";
            
            const nU = patch.controlPoints.length;
            const nV = patch.controlPoints[0].length;
            const bcs = [];
            if (stateToSolve.boundEdge === 'left') {
                for (let j = 0; j < nV; j++) bcs.push({ i: 0, j: j, axis: 'both', value: 0 });
            } else if (stateToSolve.boundEdge === 'right') {
                for (let j = 0; j < nV; j++) bcs.push({ i: nU - 1, j: j, axis: 'both', value: 0 });
            } else if (stateToSolve.boundEdge === 'bottom') {
                for (let i = 0; i < nU; i++) bcs.push({ i: i, j: 0, axis: 'both', value: 0 });
            } else if (stateToSolve.boundEdge === 'top') {
                for (let i = 0; i < nU; i++) bcs.push({ i: i, j: nV - 1, axis: 'both', value: 0 });
            }
            
            let fx = 0, fy = 0;
            if (stateToSolve.forceDir === 'horizontal') fx = stateToSolve.load;
            else fy = -stateToSolve.load;

            const loads = [{ i: stateToSolve.forceU, j: stateToSolve.forceV, fx: fx, fy: fy }];

            try {
                if (stateToSolve.mode === 'linear') {
                    analysisData.u = solver.solve(patch, bcs, loads);
                    analysisData.residuals = [1.0, 1e-15];
                } else {
                    const result = solver.solveNonlinear(patch, bcs, loads, { steps: 2, iterations: 10 });
                    analysisData.u = result.u;
                    analysisData.residuals = result.residualHistory;
                }
                
                let hasNaN = false;
                if(analysisData.u) {
                    for(let i=0; i<analysisData.u.length; i++) if(isNaN(analysisData.u[i])) hasNaN = true;
                }
                if(hasNaN) console.error("Solver produced NaN displacements!");
                
                const t1 = performance.now();
                updateChart();
                updateSurface();
                updateForceArrow(stateToSolve.forceU, stateToSolve.forceV, stateToSolve.load, stateToSolve.forceDir);
                updateBoundaryVisuals(stateToSolve.boundEdge);
                
                document.getElementById('conv-stat').textContent = `Success: Solved in ${ (t1 - t0).toFixed(1) }ms`;
                activeState = { ...stateToSolve };
            } catch (err) {
                console.error("Solver exception:", err);
                document.getElementById('conv-stat').textContent = `Error: ${err.message}`;
            }
            isSolving = false;
        }
    }
    requestAnimationFrame(solverLoop);
}

function updateChart() {
    chart.data.labels = analysisData.residuals.map((_, i) => i + 1);
    chart.data.datasets[0].data = analysisData.residuals;
    chart.update();
}

function setupUI() {
    const forceUInput = document.getElementById('force-u');
    const forceVInput = document.getElementById('force-v');
    const forceDirInput = document.getElementById('force-dir');
    const boundEdgeInput = document.getElementById('bound-edge');

    forceUInput.oninput = (e) => targetState.forceU = parseInt(e.target.value) || 0;
    forceVInput.oninput = (e) => targetState.forceV = parseInt(e.target.value) || 0;
    forceDirInput.onchange = (e) => targetState.forceDir = e.target.value;
    boundEdgeInput.onchange = (e) => targetState.boundEdge = e.target.value;

    document.getElementById('solve-btn').onclick = () => { /* Now Realtime, button config optional */ };
    document.getElementById('solve-btn').style.display = 'none'; // Hide since it's realtime now
    
    document.getElementById('mode-linear').onclick = () => {
        targetState.mode = 'linear';
        document.getElementById('mode-linear').className = "py-2 px-3 rounded-lg border border-blue-500 bg-blue-500/10 text-blue-400 text-xs font-bold transition-colors";
        document.getElementById('mode-nonlinear').className = "py-2 px-3 rounded-lg border border-slate-700 text-xs font-bold hover:bg-slate-800 transition-colors";
    };
    
    document.getElementById('mode-nonlinear').onclick = () => {
        targetState.mode = 'nonlinear';
        document.getElementById('mode-nonlinear').className = "py-2 px-3 rounded-lg border border-rose-500 bg-rose-500/10 text-rose-400 text-xs font-bold transition-colors";
        document.getElementById('mode-linear').className = "py-2 px-3 rounded-lg border border-slate-700 text-xs font-bold hover:bg-slate-800 transition-colors";
    };

    document.getElementById('load-slider').oninput = (e) => {
        document.getElementById('load-val').textContent = `${e.target.value} N`;
        targetState.load = parseFloat(e.target.value);
    };

    // --- Raycasting Interaction ---
    const raycaster = new THREE.Raycaster();
    const mouse = new THREE.Vector2();
    
    renderer.domElement.addEventListener('pointerdown', (e) => {
        const rect = renderer.domElement.getBoundingClientRect();
        mouse.x = ((e.clientX - rect.left) / rect.width) * 2 - 1;
        mouse.y = -((e.clientY - rect.top) / rect.height) * 2 + 1;
        
        raycaster.setFromCamera(mouse, camera);
        if (surfaceMesh) {
            const intersects = raycaster.intersectObject(surfaceMesh);
            if (intersects.length > 0) {
                const a = intersects[0].face.a;
                const res = 40;
                
                // Map vertex index to parameter geometry
                const i_grid = Math.floor(a / (res+1));
                const j_grid = a % (res+1);
                
                const u_param = i_grid / res;
                const v_param = j_grid / res;
                
                const nU = patch.controlPoints.length;
                const nV = patch.controlPoints[0].length;
                
                const bestU = Math.min(Math.max(Math.round(u_param * (nU - 1)), 0), nU - 1);
                const bestV = Math.min(Math.max(Math.round(v_param * (nV - 1)), 0), nV - 1);
                
                forceUInput.value = bestU;
                forceVInput.value = bestV;
                targetState.forceU = bestU;
                targetState.forceV = bestV;
            }
        }
    });

    solverLoop(); // Start Async Polling
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
