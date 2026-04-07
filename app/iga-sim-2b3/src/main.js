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

// The preset now generates a refined 12-element quadrant
let basePatch = NURBSPresets.generatePlateWithHole12(R, L);
let patch = { ...basePatch };

let analysisData = {
    u: null,
    load: 100
};

let targetState = { 
    load: 100,
    h: 0,
    p: 2,
    k: false,
    E: 100000,
    nu: 0.30
};
let activeState = { load: -1, h: -1, p: -1, k: false, E: -1, nu: -1 };
let isSolving = false;

function updateStatusLabels() {
    const nU = patch.controlPoints.length;
    const nV = patch.controlPoints[0].length;
    if (document.getElementById('status-degree')) document.getElementById('status-degree').textContent = `p=${patch.p}, q=${patch.q}`;
    if (document.getElementById('status-spans')) document.getElementById('status-spans').textContent = `${patch.U.length - 2 * patch.p - 1} x ${patch.V.length - 2 * patch.q - 1}`;
    if (document.getElementById('status-points')) document.getElementById('status-points').textContent = `${nU} x ${nV}`;
    
    // Show warning if complexity is high
    const warning = document.getElementById('refinement-warning');
    if (warning) {
        warning.classList.toggle('hidden', targetState.h < 3 && targetState.p < 4);
    }
}

function applyRefinements() {
    // Reset to Base
    patch = JSON.parse(JSON.stringify(basePatch));
    
    // Apply Refinement Stack (p elevation then h subdivision)
    RefineUtils.apply(engine, patch, { 
        h: targetState.h, 
        p: targetState.p 
    });
    
    updateStatusLabels();
    updateSurface();
    updateBoundaryVisuals();
    updateControlPoints();
    updateForceArrows();
}

function init() {
    setupUI();
    updateSurface();
    updateBoundaryVisuals();
    updateForceArrows();
    animate();
}

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
    updateBoundaryVisuals();
    updateForceArrows();
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

    const samples = 15;
    
    // Sample along the outer boundary (v=1)
    for (let i = 0; i <= samples; i++) {
        const u = i / samples;
        const v = 1.0;
        const state = engine.getSurfaceState(patch, u, v);
        const cp = state.position;

        // Only show arrows on the RIGHT boundary (x = L)
        if (Math.abs(cp.x - L) < 1e-3) {
            const dir = new THREE.Vector3(1, 0, 0);
            const arrow = new THREE.ArrowHelper(dir, new THREE.Vector3(cp.x, cp.y, cp.z), targetState.load/250, 0xef4444);
            scene.add(arrow);
            forceVisuals.push(arrow);
        }
    }
}

async function solverLoop() {
    const stateChanged = targetState.load !== activeState.load || 
                         targetState.h !== activeState.h || 
                         targetState.p !== activeState.p ||
                         targetState.k !== activeState.k ||
                         targetState.E !== activeState.E ||
                         targetState.nu !== activeState.nu;

    if (!isSolving && stateChanged) {
        isSolving = true;

        // Sync Material Properties to Solver
        solver.E = targetState.E;
        solver.nu = targetState.nu;
        
        // If topology changed, re-apply refinements
        if (targetState.h !== activeState.h || targetState.p !== activeState.p || targetState.k !== activeState.k) {
            applyRefinements();
        }

        await new Promise(r => setTimeout(r, 10));
        const t0 = performance.now();
        const nU = patch.controlPoints.length;
        const nV = patch.controlPoints[0].length;

        const bcs = [];
        // Angular edge i=0 is bottom (y=0) -> Fix Y
        if (patch.controlPoints[0]) {
            for(let j=0; j<nV; j++) bcs.push({ i: 0, j: j, axis: 'y', value: 0 });
        }
        // Angular edge i=nU-1 is left (x=0) -> Fix X
        if (patch.controlPoints[nU-1]) {
            for(let j=0; j<nV; j++) bcs.push({ i: nU-1, j: j, axis: 'x', value: 0 });
        }
        
        const loads = [];
        // Apply traction only to the Right edge (x = L)
        for(let i=0; i<nU; i++) {
            const cp = patch.controlPoints[i][nV-1];
            if (cp && Math.abs(cp.x - L) < 1e-3) {
                // For uniform traction, internal nodes get 1.0 weight and ends get 0.5
                const weight = (i === 0 || i === nU-1) ? 0.5 : 1.0;
                loads.push({ i: i, j: nV-1, fx: targetState.load * weight, fy: 0 });
            }
        }

        try {
            analysisData.u = solver.solve(patch, bcs, loads);
            const t1 = performance.now();
            const maxD = calculateMaxDisp(analysisData.u);
            if (document.getElementById('max-disp')) document.getElementById('max-disp').textContent = `${maxD.toFixed(4)} mm`;
            
            updateSurface();
            updateBoundaryVisuals();
            updateControlPoints();
            updateForceArrows();
            if (document.getElementById('conv-stat')) document.getElementById('conv-stat').textContent = `Success: Solved in ${ (t1 - t0).toFixed(1) }ms`;
            activeState = JSON.parse(JSON.stringify(targetState));
        } catch (err) {
            if (document.getElementById('conv-stat')) document.getElementById('conv-stat').textContent = `Error: ${err.message}`;
        }
        isSolving = false;
    }
    requestAnimationFrame(solverLoop);
}

function setupUI() {
    // Load Slider
    document.getElementById('load-slider').oninput = (e) => {
        const val = parseFloat(e.target.value);
        if (document.getElementById('load-val')) document.getElementById('load-val').textContent = `${val} N/mm`;
        targetState.load = val;
    };

    // Material Sliders
    if (document.getElementById('e-slider')) {
        document.getElementById('e-slider').oninput = (e) => {
            const val = parseFloat(e.target.value);
            if (document.getElementById('e-val')) document.getElementById('e-val').textContent = val >= 10000 ? `${(val / 1000).toFixed(0)}k` : val;
            targetState.E = val;
        };
    }
    if (document.getElementById('nu-slider')) {
        document.getElementById('nu-slider').oninput = (e) => {
            const val = parseFloat(e.target.value);
            if (document.getElementById('nu-val')) document.getElementById('nu-val').textContent = val.toFixed(2);
            targetState.nu = val;
        };
    }

    // Refinement Sliders
    if (document.getElementById('h-slider')) {
        document.getElementById('h-slider').oninput = (e) => {
            targetState.h = parseInt(e.target.value);
        };
    }
    if (document.getElementById('p-slider')) {
        document.getElementById('p-slider').oninput = (e) => {
            targetState.p = parseInt(e.target.value);
        };
    }

    // Reset Button
    if (document.getElementById('reset-ref')) {
        document.getElementById('reset-ref').onclick = () => {
            targetState.h = 0;
            targetState.p = 2;
            if (document.getElementById('h-slider')) document.getElementById('h-slider').value = 0;
            if (document.getElementById('p-slider')) document.getElementById('p-slider').value = 2;
        };
    }

    // K-Refinement Shortcut
    if (document.getElementById('apply-k')) {
        document.getElementById('apply-k').onclick = () => {
            targetState.p = Math.min(targetState.p + 1, 5);
            targetState.h = Math.min(targetState.h + 1, 4);
            if (document.getElementById('p-slider')) document.getElementById('p-slider').value = targetState.p;
            if (document.getElementById('h-slider')) document.getElementById('h-slider').value = targetState.h;
        };
    }
    
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

    updateStatusLabels();
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
