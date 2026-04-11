/**
 * Phase 2D.3 | Online Real-Time Solver Digital Twin
 */

const engine = new NURBS2D();
const solver = new IGA2DSolver(engine);
let patch = null;
let last_u_tilde = null;

// Extracted ROM data
let ROM = null; 
let nDofs = 0;

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

// Render loop
function animate() {
    requestAnimationFrame(animate);
    controls.update();
    renderer.render(scene, camera);
}
animate();

// --- UI AND FILE PARSER ---

const uploadInput = document.getElementById('upload-input');
const uploadZone = document.getElementById('upload-zone');
const E_slider = document.getElementById('E-slider');
const nu_slider = document.getElementById('nu-slider');
const F_slider = document.getElementById('F-slider');

uploadZone.addEventListener('click', () => uploadInput.click());

uploadInput.addEventListener('change', (e) => {
    const file = e.target.files[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = (ev) => {
        try {
            const data = JSON.parse(ev.target.result);
            if (!data.ROM || !data.ROM.Phi) throw new Error("Invalid ROM format");
            loadROM(data);
        } catch(err) {
            alert("Error parsing JSON: " + err.message);
        }
    };
    reader.readAsText(file);
});

function loadROM(data) {
    // 1. Rebuild mesh patch from metadata so we can visualize
    const sm = data.meta.mesh;
    patch = {
        U: sm.knotU,
        V: sm.knotV,
        p: sm.p,
        q: sm.q,
        controlPoints: sm.controlPoints_base,
        weights: sm.weights
    };
    nDofs = patch.controlPoints.length * patch.controlPoints[0].length * 2;

    // 2. Load Operators
    const { Matrix } = mlMatrix;
    ROM = {
        k: data.meta.modesRetained,
        Phi: new Matrix(data.ROM.Phi),
        K1: new Matrix(data.ROM.K1_red),
        K2: new Matrix(data.ROM.K2_red),
        Fr: new Float64Array(data.ROM.Fr)
    };

    document.getElementById('upload-zone').classList.add('hidden');
    document.getElementById('load-status').classList.remove('hidden');
    document.getElementById('k-info').innerText = ROM.k;
    
    // Enable online controls
    const solverCtrl = document.getElementById('solver-controls');
    solverCtrl.classList.remove('opacity-50', 'pointer-events-none');

    // Initial render
    updateSurface(new Float64Array(nDofs));
    
    // Initial evaluation
    evaluateOnline();
}

[E_slider, nu_slider, F_slider].forEach(el => {
    el.addEventListener('input', () => {
        document.getElementById('E-val').innerText = parseInt(E_slider.value).toLocaleString() + " MPa";
        document.getElementById('nu-val').innerText = parseFloat(nu_slider.value).toFixed(2);
        document.getElementById('F-val').innerText = parseInt(F_slider.value) + " N";
        if (ROM) evaluateOnline();
    });
});

document.getElementById('view-mode').addEventListener('change', () => {
    if (ROM && last_u_tilde) {
        updateSurface(last_u_tilde, parseFloat(E_slider.value), parseFloat(nu_slider.value));
    }
});

// --- CORE ONLINE SOLVER MATH ---

function gaussianElimination(A, b) {
    const n = b.length;
    for (let i = 0; i < n; i++) {
        let max = i;
        for (let j = i + 1; j < n; j++) if (Math.abs(A[j][i]) > Math.abs(A[max][i])) max = j;
        [A[i], A[max]] = [A[max], A[i]]; [b[i], b[max]] = [b[max], b[i]];
        if (Math.abs(A[i][i]) < 1e-15) { A[i][i] = 1.0; b[i] = 0; for (let k = i+1; k < n; k++) A[i][k] = 0; }
        for (let j = i + 1; j < n; j++) {
            const f = A[j][i] / A[i][i];
            b[j] -= f * b[i];
            for (let k = i; k < n; k++) A[j][k] -= f * A[i][k];
        }
    }
    const x = new Float64Array(n);
    for (let i = n - 1; i >= 0; i--) {
        let s = 0;
        for (let j = i + 1; j < n; j++) s += A[i][j] * x[j];
        x[i] = (b[i] - s) / A[i][i];
    }
    return x;
}

function evaluateOnline() {
    if (!ROM) return;
    const t0 = performance.now();
    
    // Parse user parameters
    const E = parseFloat(E_slider.value);
    const nu = parseFloat(nu_slider.value);
    
    // Parameterized stiffness factor
    const factor = E / (1 - nu * nu);
    
    // Kr(E, nu) = factor * (K1 + nu * K2)
    // To do this fast, we extract data to arrays, avoiding heavy allocations
    const K1_arr = ROM.K1.to2DArray();
    const K2_arr = ROM.K2.to2DArray();
    const Kr_arr = [];
    for(let i=0; i<ROM.k; i++){
        Kr_arr[i] = [];
        for(let j=0; j<ROM.k; j++){
            Kr_arr[i][j] = factor * (K1_arr[i][j] + nu * K2_arr[i][j]);
        }
    }

    const forceMultiplier = parseFloat(F_slider.value) / 100.0;
    const Fr_arr = new Float64Array(ROM.Fr.length);
    for(let i=0; i<ROM.Fr.length; i++) Fr_arr[i] = ROM.Fr[i] * forceMultiplier;
    
    // Online Solve (k x k dense system)
    const q_arr = gaussianElimination(Kr_arr, Fr_arr);
    const t1 = performance.now();

    // Expansion (u = Phi * q)
    const { Matrix } = mlMatrix;
    const q_mat = new Matrix([Array.from(q_arr)]).transpose();
    const u_tilde = ROM.Phi.mmul(q_mat).to1DArray();
    const t2 = performance.now();

    // Update charts
    document.getElementById('solve-time').innerText = (t1 - t0).toFixed(2) + " ms";
    document.getElementById('recon-time').innerText = (t2 - t1).toFixed(2) + " ms";
    
    const approxFPS = Math.min(1000 / ((t2 - t0) || 1), 144);
    document.getElementById('fps-val').innerText = approxFPS.toFixed(0);

    // Warp visualization
    last_u_tilde = u_tilde;
    updateSurface(u_tilde, E, nu);
}

// --- VISUALIZATION HELPERS (No Physics!) ---

function updateSurface(u_disp = null, E = 100000, nu = 0.3) {
    if (surfaceMesh) scene.remove(surfaceMesh);
    if (wireframeOverlay) scene.remove(wireframeOverlay);
    if (!patch) return;

    const res = 40; 
    const geometry = new THREE.BufferGeometry();
    const positions = [], colors = [];
    
    const viewMode = document.getElementById('view-mode').value;
    const isStress = viewMode === 'stress';
    
    // First Pass: Calculate values and bounds
    let maxVal = 1e-6; // prevent division by zero
    const vertexData = [];

    for (let i = 0; i <= res; i++) {
        const u = i / res;
        for (let j = 0; j <= res; j++) {
            const v = j / res;
            const state = engine.getSurfaceState(patch, u, v);
            let p = { ...state.position };
            
            let val = 0;
            let interp = { x: 0, y: 0 };
            
            if (u_disp && state.denominator > 0) {
                interp = interpolateDisplacement(u, v, state.denominator, u_disp);
                
                if (isStress) {
                    try {
                        const stressState = solver.getNumericalStress(patch, u_disp, u, v, E, nu);
                        val = stressState.vonMises;
                        if (isNaN(val)) val = 0;
                    } catch (e) {
                        document.getElementById('legend-title').innerText = "ERROR: " + e.message;
                        val = 0;
                    }
                } else {
                    val = Math.abs(interp.x);
                }
                
                if (val > maxVal) maxVal = val;
            }
            
            vertexData.push({ p, interp, val, hasDisp: u_disp && state.denominator > 0 });
        }
    }

    // Legend Update
    if (u_disp) {
        document.getElementById('legend-container').classList.remove('hidden');
        if (!document.getElementById('legend-title').innerText.startsWith('ERROR')) {
            document.getElementById('legend-title').innerText = isStress ? "VON MISES STRESS (MPa)" : "X-DISPLACEMENT (mm)";
        }
        document.getElementById('legend-max').innerText = isStress ? maxVal.toFixed(1) : maxVal.toExponential(2);
    } else {
        document.getElementById('legend-container').classList.add('hidden');
    }

    // Second Pass: Color mapping and geometry building
    for (let idx = 0; idx < vertexData.length; idx++) {
        const data = vertexData[idx];
        let p = data.p;
        
        if (data.hasDisp) {
            // Simplify scaling purely for visuals
            p.x += data.interp.x * 50; 
            p.y += data.interp.y * 50;
            
            const c = getHeatmapColor(data.val / maxVal);
            colors.push(c.r, c.g, c.b);
        } else {
            colors.push(0.1, 0.2, 0.4);
        }
        positions.push(p.x, p.y, p.z);
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
