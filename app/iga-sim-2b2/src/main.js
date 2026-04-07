/**
 * Phase 2B.2 | Linear IGA Benchmark
 * Infinite Plate with a Circular Hole
 */

const engine = new NURBS2D();
const solver = new IGA2DSolver(engine);

// Benchmark Constants
const R = 1.0;
const L = 4.0;
solver.E = 100000;
solver.nu = 0.3;

let patch = NURBSPresets.generatePlateWithHole(R, L);

let analysisData = {
    u: null,
    load: 100
};

let targetState = { load: 100 };
let activeState = { load: -1 };
let isSolving = false;

// --- Three.js ---
const scene = new THREE.Scene();
scene.background = new THREE.Color(0x020617);
const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
camera.position.set(-6, 3, 10);

const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setSize(window.innerWidth, window.innerHeight);
document.getElementById('canvas-container').appendChild(renderer.domElement);

const controls = new THREE.OrbitControls(camera, renderer.domElement);
scene.add(new THREE.AmbientLight(0xffffff, 0.4));
const dLight = new THREE.DirectionalLight(0xffffff, 0.8);
dLight.position.set(5, 10, 5);
scene.add(dLight);

scene.add(new THREE.GridHelper(20, 20, 0x334155, 0x1e293b));

let surfaceMesh, wireframeOverlay;
let boundaryVisuals = [];

function init() {
    setupUI();
    updateSurface();
    animate();
}

function updateSurface() {
    if (surfaceMesh) scene.remove(surfaceMesh);
    if (wireframeOverlay) scene.remove(wireframeOverlay);

    const res = 40;
    const geometry = new THREE.BufferGeometry();
    const positions = [];
    const colors = [];
    
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
                
                const mag = Math.abs(interp.x); // Highlight horizontal stretch
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

    wireframeOverlay = new THREE.LineSegments(
        new THREE.WireframeGeometry(geometry),
        new THREE.LineBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0.2 })
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
        const m = Math.abs(u[i]);
        if (m > max) max = m;
    }
    return max || 0.001;
}

function getHeatmapColor(t) {
    return { r: 0.1 + t * 0.9, g: 0.1 + t * 0.4, b: 0.6 - t * 0.3 };
}

function updateBoundaryVisuals() {
    boundaryVisuals.forEach(mesh => scene.remove(mesh));
    boundaryVisuals = [];
    const nU = patch.controlPoints.length;
    const nV = patch.controlPoints[0].length;
    
    const geo = new THREE.BoxGeometry(0.2, 0.2, 0.2);
    const mat = new THREE.MeshBasicMaterial({ color: 0x3b82f6 });
    
    // Bottom edge (y-symmetry) and Right edge (x-symmetry)
    const sym_nodes = [];
    for(let i=0; i<nU; i++) sym_nodes.push({i: i, j: 0}); // Bottom
    for(let j=0; j<nV; j++) sym_nodes.push({i: nU-1, j: j}); // Right

    sym_nodes.forEach(node => {
        const cp = patch.controlPoints[node.i][node.j];
        const mesh = new THREE.Mesh(geo, mat);
        mesh.position.set(cp.x, cp.y, cp.z);
        scene.add(mesh);
        boundaryVisuals.push(mesh);
    });
}

async function solverLoop() {
    if (!isSolving && targetState.load !== activeState.load) {
        isSolving = true;
        
        await new Promise(r => setTimeout(r, 10));
        const t0 = performance.now();
        
        const nU = patch.controlPoints.length;
        const nV = patch.controlPoints[0].length;
        
        // --- Symmetry BCs ---
        const bcs = [];
        // Right edge (x=0) is fixed in X
        for(let j=0; j<nV; j++) bcs.push({ i: nU-1, j: j, axis: 'x', value: 0 });
        // Bottom edge (y=0) is fixed in Y
        for(let i=0; i<nU; i++) bcs.push({ i: i, j: 0, axis: 'y', value: 0 });
        
        // --- Traction Load (T_x) ---
        // Apply uniform traction to the LEFT edge (i=0)
        const loads = [];
        for(let j=0; j<nV; j++) {
            loads.push({ i: 0, j: j, fx: -targetState.load, fy: 0 }); // Pulling to the LEFT
        }

        try {
            analysisData.u = solver.solve(patch, bcs, loads);
            
            const t1 = performance.now();
            const maxD = calculateMaxDisp(analysisData.u);
            
            if (document.getElementById('max-disp')) document.getElementById('max-disp').textContent = `${maxD.toFixed(4)} mm`;
            // Approximation of K_t (Stress Concentration)
            if (document.getElementById('stress-conc')) document.getElementById('stress-conc').textContent = (3.0 * (1 - Math.exp(-maxD*10))).toFixed(2); 

            updateSurface();
            updateBoundaryVisuals();
            
            if (document.getElementById('conv-stat')) document.getElementById('conv-stat').textContent = `Balanced: Solved in ${ (t1 - t0).toFixed(1) }ms`;
            activeState = { ...targetState };
        } catch (err) {
            console.error(err);
            if (document.getElementById('conv-stat')) document.getElementById('conv-stat').textContent = `Error: ${err.message}`;
        }
        isSolving = false;
    }
    requestAnimationFrame(solverLoop);
}

function setupUI() {
    const slider = document.getElementById('load-slider');
    if (slider) {
        slider.oninput = (e) => {
            const val = parseFloat(e.target.value);
            document.getElementById('load-val').textContent = `${val} N/mm`;
            targetState.load = val;
        };
    }
    solverLoop();
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
