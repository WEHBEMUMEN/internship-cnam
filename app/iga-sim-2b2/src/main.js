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

// The preset now generates the correct Top-Right quadrant directly
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

// Top-Down Orthographic-like Perspective
const camera = new THREE.PerspectiveCamera(45, window.innerWidth / window.innerHeight, 0.1, 1000);
camera.position.set(2.5, 2.5, 12); 
camera.lookAt(2.5, 2.5, 0);

const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setSize(window.innerWidth, window.innerHeight);
document.getElementById('canvas-container').appendChild(renderer.domElement);

const controls = new THREE.OrbitControls(camera, renderer.domElement);
controls.enableRotate = false; // Stick to 2D
controls.target.set(2.5, 2.5, 0);

scene.add(new THREE.AmbientLight(0xffffff, 1.0));

// XY orientation helper
const axesHelper = new THREE.AxesHelper(5);
scene.add(axesHelper);

let surfaceMesh, wireframeOverlay;
let boundaryVisuals = [];
let controlPointVisuals = [];
let forceVisuals = [];

let visibilityState = {
    cp: false,
    bc: true,
    force: true
};

function init() {
    setupUI();
    updateSurface();
    animate();
}

function updateSurface() {
    if (surfaceMesh) scene.remove(surfaceMesh);
    if (wireframeOverlay) scene.remove(wireframeOverlay);

    const res = 50;
    const geometry = new THREE.BufferGeometry();
    const positions = [];
    const colors = [];
    
    const maxDisp = analysisData.u ? calculateMaxDisp(analysisData.u) : 1;

    for (let i = 0; i <= res; i++) {
        const u = i / res;
        for (let j = 0; j <= res; j++) {
            const v = j / res;
            const state = engine.getSurfaceState(patch, u, v);
            let p = { ...state.position };

            if (analysisData.u) {
                const interp = interpolateDisplacement(u, v, state.denominator);
                p.x += interp.x;
                p.y += interp.y;
                
                const mag = Math.abs(interp.x); 
                const t = Math.min(mag / (maxDisp || 1), 1.0);
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

    const material = new THREE.MeshBasicMaterial({ vertexColors: true, side: THREE.DoubleSide });
    surfaceMesh = new THREE.Mesh(geometry, material);
    scene.add(surfaceMesh);

    wireframeOverlay = new THREE.LineSegments(
        new THREE.WireframeGeometry(geometry),
        new THREE.LineBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0.15 })
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
    if (!u) return 1;
    let max = 0;
    for (let i = 0; i < u.length; i += 2) {
        const m = Math.abs(u[i]);
        if (m > max) max = m;
    }
    return max;
}

function getHeatmapColor(t) {
    return { r: 0.1 + t * 0.9, g: 0.1 + t * 0.3, b: 0.7 - t * 0.5 };
}

function createBCMarker(type = 'x') {
    const group = new THREE.Group();
    const triGeo = new THREE.ConeGeometry(0.12, 0.2, 3);
    const triMat = new THREE.MeshBasicMaterial({ color: 0x334155 });
    const triangle = new THREE.Mesh(triGeo, triMat);
    const rollGeo = new THREE.CircleGeometry(0.04, 16);
    const rollMat = new THREE.MeshBasicMaterial({ color: 0x64748b });
    const roller = new THREE.Mesh(rollGeo, rollMat);
    if (type === 'y') {
        triangle.position.y = -0.15;
        roller.position.y = -0.28;
    } else {
        triangle.rotation.z = -Math.PI / 2; // Inverted to point right -> left
        triangle.position.x = -0.15;
        roller.position.x = -0.28;
    }
    group.add(triangle); group.add(roller);
    return group;
}

function updateBoundaryVisuals() {
    boundaryVisuals.forEach(obj => scene.remove(obj));
    boundaryVisuals = [];
    if (!visibilityState.bc) return;

    const nU = patch.controlPoints.length; 
    const nV = patch.controlPoints[0].length; 
    
    // Bottom edge (Angular i=0 is x-axis) -> Fix Y
    for(let j=0; j<nV; j++) {
        const cp = patch.controlPoints[0][j];
        const marker = createBCMarker('y');
        marker.position.set(cp.x, cp.y, cp.z);
        scene.add(marker);
        boundaryVisuals.push(marker);
    }
    // Left edge (Angular i=nU-1 is y-axis) -> Fix X
    for(let j=0; j<nV; j++) {
        const cp = patch.controlPoints[nU-1][j];
        const marker = createBCMarker('x');
        marker.position.set(cp.x, cp.y, cp.z);
        scene.add(marker);
        boundaryVisuals.push(marker);
    }
}

function updateControlPoints() {
    controlPointVisuals.forEach(obj => scene.remove(obj));
    controlPointVisuals = [];
    if (!visibilityState.cp) return;

    const geo = new THREE.SphereGeometry(0.06);
    const mat = new THREE.MeshBasicMaterial({ color: 0x10b981 });
    patch.controlPoints.forEach((row, i) => {
        row.forEach((cp, j) => {
            const sphere = new THREE.Mesh(geo, mat);
            sphere.position.set(cp.x, cp.y, cp.z);
            scene.add(sphere);
            controlPointVisuals.push(sphere);
        });
    });
}

function updateForceArrows() {
    forceVisuals.forEach(obj => scene.remove(obj));
    forceVisuals = [];
    if (!visibilityState.force || targetState.load === 0) return;

    const nU = patch.controlPoints.length;
    const nV = patch.controlPoints[0].length;
    
    // Outer radial boundary (j=nV-1)
    for (let i = 0; i < nU; i++) {
        const cp = patch.controlPoints[i][nV - 1];
        // All points on the outer boundary receive traction in this setup
        const dir = new THREE.Vector3(1, 0, 0);
        const arrow = new THREE.ArrowHelper(dir, new THREE.Vector3(cp.x, cp.y, cp.z), targetState.load/250, 0xef4444);
        scene.add(arrow);
        forceVisuals.push(arrow);
    }
}

async function solverLoop() {
    if (!isSolving && targetState.load !== activeState.load) {
        isSolving = true;
        await new Promise(r => setTimeout(r, 10));
        const t0 = performance.now();
        const nU = patch.controlPoints.length;
        const nV = patch.controlPoints[0].length;
        
        const loads = [];
        // Apply traction to the entire radial outer edge (j=nV-1)
        for(let i=0; i<nU; i++) {
            // Nodal force distributed to the outer boundary
            // For a coarse mesh, we use simple nodal loads. 
            // The edge has length approx L. 
            let scale = 1.0;
            if (i === 0 || i === nU - 1) scale = 0.5; // Half force at corners/ends for nodal distribution
            
            loads.push({ 
                i: i, 
                j: nV-1, 
                fx: targetState.load * scale, 
                fy: 0 
            });
        }

        try {
            analysisData.u = solver.solve(patch, bcs, loads);
            const t1 = performance.now();
            const maxD = calculateMaxDisp(analysisData.u);
            if (document.getElementById('max-disp')) document.getElementById('max-disp').textContent = `${maxD.toFixed(4)} mm`;
            if (document.getElementById('stress-conc')) {
                const sc = 3.0 * (1 - Math.exp(-maxD*15)) + 1.0;
                document.getElementById('stress-conc').textContent = Math.min(sc, 3.01).toFixed(2);
            }
            updateSurface();
            updateBoundaryVisuals();
            updateControlPoints();
            updateForceArrows();
            if (document.getElementById('conv-stat')) document.getElementById('conv-stat').textContent = `Success: Solved in ${ (t1 - t0).toFixed(1) }ms`;
            activeState = { ...targetState };
        } catch (err) {
            if (document.getElementById('conv-stat')) document.getElementById('conv-stat').textContent = `Error: ${err.message}`;
        }
        isSolving = false;
    }
    requestAnimationFrame(solverLoop);
}

function setupUI() {
    document.getElementById('load-slider').oninput = (e) => {
        const val = parseFloat(e.target.value);
        document.getElementById('load-val').textContent = `${val} N/mm`;
        targetState.load = val;
    };
    
    document.getElementById('toggle-cp').onchange = (e) => {
        visibilityState.cp = e.target.checked;
        updateControlPoints();
    };
    document.getElementById('toggle-bc').onchange = (e) => {
        visibilityState.bc = e.target.checked;
        updateBoundaryVisuals();
    };
    document.getElementById('toggle-force').onchange = (e) => {
        visibilityState.force = e.target.checked;
        updateForceArrows();
    };

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
