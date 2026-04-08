
// Verification Lab 2.5 | Differential Area & Jacobian
// Phase 2.5 | Computational Core

const engine = new NURBS2D();
let currentMode = 'sphere';
let patch = NURBSPresets.generateSphere();

// --- Three.js Setup ---
const scene = new THREE.Scene();
scene.background = new THREE.Color(0x020617);
const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 2000);
camera.position.set(50, 50, 50);

const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setSize(window.innerWidth, window.innerHeight);
document.getElementById('canvas-container').appendChild(renderer.domElement);

const orbitControls = new THREE.OrbitControls(camera, renderer.domElement);
const transformControls = new THREE.TransformControls(camera, renderer.domElement);
scene.add(transformControls);

scene.add(new THREE.AmbientLight(0xffffff, 0.4));
const dLight = new THREE.DirectionalLight(0xffffff, 0.8);
dLight.position.set(50, 100, 50);
scene.add(dLight);
scene.add(new THREE.GridHelper(100, 20, 0x334155, 0x1e293b));

// --- Objects ---
let surfaceMesh, wireframeOverlay, normalHelpers;
let gridGroup = new THREE.Group();
let pointMeshes = []; 
scene.add(gridGroup);

const raycaster = new THREE.Raycaster();
const mouse = new THREE.Vector2();

// --- Colormap Helper ---
function getTurboColor(t) {
    // Simplified Turbo-like colormap
    const r = Math.min(Math.max(1.5 * t - 0.5, 0), 1);
    const g = 1 - Math.abs(2 * t - 1);
    const b = Math.min(Math.max(1.5 * (1 - t) - 0.5, 0), 1);
    return { r, g, b };
}

function fullRebuild() {
    transformControls.detach();
    pointMeshes = [];
    gridGroup.clear();
    createSurface();
    updateGrid();
    syncUI();
}

function syncUI() {
    document.getElementById('state-points').textContent = `Lattice: ${patch.controlPoints.length} x ${patch.controlPoints[0].length}`;
}

function createSurface() {
    if (window.perfMonitor) window.perfMonitor.startMeasure();
    if (surfaceMesh) scene.remove(surfaceMesh);
    if (wireframeOverlay) scene.remove(wireframeOverlay);
    if (normalHelpers) scene.remove(normalHelpers);
    
    const res = parseInt(document.getElementById('resolution-slider').value);
    const geometry = new THREE.BufferGeometry();
    const positions = [];
    const colors = [];
    
    let detMin = Infinity, detMax = -Infinity;
    const dets = [];

    // First pass to find range
    for (let i = 0; i <= res; i++) {
        const u = i / res;
        for (let j = 0; j <= res; j++) {
            const v = j / res;
            const det = engine.getJacobianDeterminant(patch, u, v);
            dets.push(det);
            detMin = Math.min(detMin, det);
            detMax = Math.max(detMax, det);
        }
    }

    document.getElementById('state-det-range').textContent = `|J|: ${detMin.toFixed(2)} to ${detMax.toFixed(2)}`;

    let idx = 0;
    for (let i = 0; i <= res; i++) {
        const u = i / res;
        for (let j = 0; j <= res; j++) {
            const v = j / res;
            const state = engine.getSurfaceState(patch, u, v);
            positions.push(state.position.x, state.position.y, state.position.z);
            
            if (document.getElementById('show-jacobian').checked) {
                const det = dets[idx++];
                const range = detMax - detMin || 1.0;
                const t = (det - detMin) / range;
                const c = getTurboColor(t);
                colors.push(c.r, c.g, c.b);
            } else {
                colors.push(0.3, 0.4, 0.9); // Default blue
                idx++;
            }
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
        new THREE.LineBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0.1 })
    );
    scene.add(wireframeOverlay);

    if (document.getElementById('show-normals').checked) {
        normalHelpers = new THREE.Group();
        const skip = Math.ceil(res / 10);
        for (let i = 0; i <= res; i += skip) {
            for (let j = 0; j <= res; j += skip) {
                const u = i / res;
                const v = j / res;
                const pos = engine.evaluateSurface(patch, u, v);
                const norm = engine.getSurfaceNormal(patch, u, v);
                const arrow = new THREE.ArrowHelper(
                    new THREE.Vector3(norm.x, norm.y, norm.z),
                    new THREE.Vector3(pos.x, pos.y, pos.z),
                    4, 0xfacc15, 1, 0.5
                );
                normalHelpers.add(arrow);
            }
        }
        scene.add(normalHelpers);
    }

    if (window.perfMonitor) window.perfMonitor.endMeasure();
}

function updateGrid() {
    gridGroup.clear();
    const lineMat = new THREE.LineBasicMaterial({ color: 0x475569, transparent: true, opacity: 0.4 });
    const n = patch.controlPoints.length;
    const m = patch.controlPoints[0].length;

    if (pointMeshes.length === 0) {
        const sphereGeom = new THREE.SphereGeometry(0.5, 12, 12);
        const sphereMat = new THREE.MeshBasicMaterial({ color: 0xfbbf24 });
        for (let i = 0; i < n; i++) {
            for (let j = 0; j < m; j++) {
                const sp = new THREE.Mesh(sphereGeom, sphereMat.clone());
                sp.userData = { i, j };
                pointMeshes.push(sp);
            }
        }
    }

    if (document.getElementById('show-grid').checked) {
        pointMeshes.forEach(sp => {
            const { i, j } = sp.userData;
            const cp = patch.controlPoints[i][j];
            sp.position.set(cp.x, cp.y, cp.z);
            gridGroup.add(sp);
        });

        for (let i = 0; i < n; i++) {
            const pts = [];
            for (let j = 0; j < m; j++) {
                const cp = patch.controlPoints[i][j];
                pts.push(cp.x, cp.y, cp.z);
            }
            gridGroup.add(new THREE.Line(new THREE.BufferGeometry().setAttribute('position', new THREE.Float32BufferAttribute(pts, 3)), lineMat));
        }
        for (let j = 0; j < m; j++) {
            const pts = [];
            for (let i = 0; i < n; i++) {
                const cp = patch.controlPoints[i][j];
                pts.push(cp.x, cp.y, cp.z);
            }
            gridGroup.add(new THREE.Line(new THREE.BufferGeometry().setAttribute('position', new THREE.Float32BufferAttribute(pts, 3)), lineMat));
        }
    }
}

// Interaction
function onPointerDown(event) {
    const rect = renderer.domElement.getBoundingClientRect();
    mouse.x = ((event.clientX - rect.left) / rect.width) * 2 - 1;
    mouse.y = -((event.clientY - rect.top) / rect.height) * 2 + 1;
    raycaster.setFromCamera(mouse, camera);
    const intersects = raycaster.intersectObjects(pointMeshes);
    if (intersects.length > 0) {
        const s = intersects[0].object;
        transformControls.attach(s);
    } else if (!transformControls.dragging) {
        transformControls.detach();
    }
}

function onPointerMove(event) {
    if (transformControls.dragging) return;
    
    const rect = renderer.domElement.getBoundingClientRect();
    mouse.x = ((event.clientX - rect.left) / rect.width) * 2 - 1;
    mouse.y = -((event.clientY - rect.top) / rect.height) * 2 + 1;
    raycaster.setFromCamera(mouse, camera);
    
    if (!surfaceMesh) return;
    const intersects = raycaster.intersectObject(surfaceMesh);
    const overlay = document.getElementById('probe-overlay');
    
    if (intersects.length > 0) {
        const intersection = intersects[0];
        // Note: For BufferGeometry, we don't directly get U,V. 
        // We estimate it from the face and interpolate for this lab.
        // Simplified: use raycaster UV if available (requires proper UV mapping in geometry)
        if (intersection.uv) {
            const u = intersection.uv.x;
            const v = intersection.uv.y;
            updateProbe(u, v);
            overlay.style.opacity = "1";
        }
    } else {
        overlay.style.opacity = "0";
    }
}

function updateProbe(u, v) {
    const det = engine.getJacobianDeterminant(patch, u, v);
    const norm = engine.getSurfaceNormal(patch, u, v);
    
    document.getElementById('p-coords').textContent = `(${u.toFixed(2)}, ${v.toFixed(2)})`;
    document.getElementById('p-det').textContent = det.toFixed(4);
    document.getElementById('p-normal').textContent = `(${norm.x.toFixed(2)}, ${norm.y.toFixed(2)}, ${norm.z.toFixed(2)})`;
}

transformControls.addEventListener('dragging-changed', (e) => orbitControls.enabled = !e.value);
transformControls.addEventListener('objectChange', () => {
    const obj = transformControls.object;
    if (obj) {
        const { i, j } = obj.userData;
        patch.controlPoints[i][j] = { x: obj.position.x, y: obj.position.y, z: obj.position.z };
        createSurface();
        updateGrid();
    }
});

// UI Handlers
document.getElementById('load-plane').onclick = () => {
    currentMode = 'plane';
    patch = NURBSPresets.generateSheet();
    fullRebuild();
};
document.getElementById('load-sphere').onclick = () => {
    currentMode = 'sphere';
    patch = NURBSPresets.generateSphere();
    fullRebuild();
};

document.getElementById('resolution-slider').oninput = createSurface;
document.getElementById('show-jacobian').onchange = createSurface;
document.getElementById('show-normals').onchange = createSurface;
document.getElementById('show-grid').onchange = () => {
    updateGrid();
    createSurface();
};

window.addEventListener('resize', () => {
    renderer.setSize(window.innerWidth, window.innerHeight);
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
});

renderer.domElement.addEventListener('pointerdown', onPointerDown);
// renderer.domElement.addEventListener('pointermove', onPointerMove); // Disable for now to focus on static heatmap

function animate() {
    requestAnimationFrame(animate);
    orbitControls.update();
    renderer.render(scene, camera);
}

fullRebuild();
animate();
