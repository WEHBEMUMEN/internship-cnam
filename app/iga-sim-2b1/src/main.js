/**
 * Phase 2B.1 | Stiffness Matrix Assembly Lab
 * Implements 2D NURBS Linear Elastic Solver
 */

const engine = new NURBS2D();
const solver = new IGA2DSolver(engine);

let patch = NURBSPresets.generateSheet(); // 10x10 Square Plate, Quadratic
let displacements = null;
let maxDisplacement = 0;

// --- Three.js Setup ---
const scene = new THREE.Scene();
scene.background = new THREE.Color(0x020617);
const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
camera.position.set(15, 15, 20);

const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setSize(window.innerWidth, window.innerHeight);
document.getElementById('canvas-container').appendChild(renderer.domElement);

const controls = new THREE.OrbitControls(camera, renderer.domElement);
scene.add(new THREE.AmbientLight(0xffffff, 0.5));
const dLight = new THREE.DirectionalLight(0xffffff, 0.8);
dLight.position.set(10, 20, 10);
scene.add(dLight);
scene.add(new THREE.GridHelper(20, 10, 0x334155, 0x1e293b));

// --- Objects ---
let surfaceMesh, wireframeOverlay;
let latticeGroup = new THREE.Group();
scene.add(latticeGroup);

function fullRebuild() {
    createMesh();
    updateUI();
}

function updateUI() {
    const nDofs = patch.controlPoints.length * patch.controlPoints[0].length * 2;
    document.getElementById('num-dofs').textContent = nDofs;
}

function createMesh() {
    if (surfaceMesh) scene.remove(surfaceMesh);
    if (wireframeOverlay) scene.remove(wireframeOverlay);
    
    const res = 40; // Rendering resolution
    const geometry = new THREE.BufferGeometry();
    const positions = [];
    const colors = [];
    
    // Displacement Visual Scaling
    const visualScale = 10.0;

    for (let i = 0; i <= res; i++) {
        const u = i / res;
        for (let j = 0; j <= res; j++) {
            const v = j / res;
            
            // 1. Calculate base physical coordinates
            const state = engine.getSurfaceState(patch, u, v);
            let pos = { ...state.position };

            // 2. Interpolate Displacement Field if available
            if (displacements) {
                const nU = patch.controlPoints.length;
                const nV = patch.controlPoints[0].length;
                
                let du = 0, dv = 0;
                for (let k = 0; k < nU; k++) {
                    const Nk = engine.basis1D(k, patch.p, patch.U, u);
                    if (Nk === 0) continue;
                    for (let l = 0; l < nV; l++) {
                        const Ml = engine.basis1D(l, patch.q, patch.V, v);
                        const R = (Nk * Ml * patch.weights[k][l]) / state.denominator;
                        
                        const baseIdx = (k * nV + l) * 2;
                        du += R * displacements[baseIdx];
                        dv += R * displacements[baseIdx + 1];
                    }
                }
                
                pos.x += du * visualScale;
                pos.y += dv * visualScale;
                
                // Color based on displacement magnitude
                const mag = Math.sqrt(du * du + dv * dv);
                const t = Math.min(mag / (maxDisplacement || 1.0), 1.0);
                const color = getHeatmapColor(t);
                colors.push(color.r, color.g, color.b);
            } else {
                colors.push(0.2, 0.4, 0.8);
            }

            positions.push(pos.x, pos.y, pos.z);
        }
    }

    const indices = [];
    for (let i = 0; i < res; i++) {
        for (let j = 0; j < res; j++) {
            const a = i * (res + 1) + j;
            const b = (i + 1) * (res + 1) + j;
            const c = (i + 1) * (res + 1) + (j + 1);
            const d = i * (res + 1) + (j + 1);
            indices.push(a, b, d, b, c, d);
        }
    }

    geometry.setIndex(indices);
    geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
    geometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
    geometry.computeVertexNormals();

    const material = new THREE.MeshStandardMaterial({ 
        vertexColors: true, side: THREE.DoubleSide, transparent: true, opacity: 0.9 
    });
    surfaceMesh = new THREE.Mesh(geometry, material);
    scene.add(surfaceMesh);

    wireframeOverlay = new THREE.LineSegments(
        new THREE.WireframeGeometry(geometry),
        new THREE.LineBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0.15 })
    );
    scene.add(wireframeOverlay);
}

function getHeatmapColor(t) {
    // Blue to Red
    return {
        r: t,
        g: 0.2,
        b: 1 - t
    };
}

async function runAnalysis() {
    if (window.perfMonitor) window.perfMonitor.startMeasure();
    
    // Set Solver state from UI
    solver.E = parseFloat(document.getElementById('e-slider').value) * 1000; // to MPa
    solver.nu = parseFloat(document.getElementById('nu-slider').value) / 100.0;
    const loadValue = parseFloat(document.getElementById('load-slider').value);

    const nU = patch.controlPoints.length;
    const nV = patch.controlPoints[0].length;
    
    // 1. Setup BCs (Clamped Left Edge)
    const bcs = [];
    for (let j = 0; j < nV; j++) {
        bcs.push({ i: 0, j: j, axis: 'both', value: 0 });
    }

    // 2. Setup Loading (Point load at center-right)
    const loads = [{ i: nU - 1, j: Math.floor(nV / 2), fx: 0, fy: -loadValue }];

    // 3. Solve
    console.log("Starting 2D Stiffness Assembly...");
    displacements = solver.solve(patch, bcs, loads);
    console.log("Solution complete.");

    // 4. Calculate Max Displacement for Scaling
    maxDisplacement = 0;
    for (let i = 0; i < displacements.length; i += 2) {
        const mag = Math.sqrt(displacements[i]**2 + displacements[i+1]**2);
        if (mag > maxDisplacement) maxDisplacement = mag;
    }

    document.getElementById('max-disp').textContent = `${maxDisplacement.toFixed(4)} mm`;

    // 5. Compare with Beam Theory (Simple Cantilever Beam L=10, H=10, Elasticity)
    // d_max = (P * L^3) / (3 * E * I)
    const L = 10.0;
    const H = 10.0;
    const I = (H * H * H) / 12; // Per unit width, but our plate is H wide...
    // Actually, I_plate = (width * thickness^3)/12. 
    // In our 2D Plane stress, 'thickness' is solver.thickness (1.0).
    const I_beam = (H * 1.0 * 1.0 * 1.0) / 12; 
    const analytical = (loadValue * Math.pow(L, 3)) / (3 * solver.E * I_beam);
    
    const error = Math.abs((maxDisplacement - analytical) / analytical) * 100;
    document.getElementById('error-theory').textContent = `${error.toFixed(2)}%`;

    createMesh();
    if (window.perfMonitor) window.perfMonitor.endMeasure();
}

// UI Handlers
document.getElementById('solve-btn').onclick = runAnalysis;

document.getElementById('e-slider').oninput = (e) => {
    document.getElementById('e-val').textContent = `${e.target.value} GPa`;
};
document.getElementById('nu-slider').oninput = (e) => {
    document.getElementById('nu-val').textContent = (e.target.value / 100).toFixed(2);
};
document.getElementById('load-slider').oninput = (e) => {
    document.getElementById('load-val').textContent = e.target.value;
};

window.addEventListener('resize', () => {
    renderer.setSize(window.innerWidth, window.innerHeight);
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
});

function animate() {
    requestAnimationFrame(animate);
    controls.update();
    renderer.render(scene, camera);
}

fullRebuild();
animate();
